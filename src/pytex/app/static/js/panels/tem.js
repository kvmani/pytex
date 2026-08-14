/**
 * The TEM solver panel: open a pattern, calibrate, pick, index, choose, tilt.
 *
 * The picking surface is the whole point. A diffraction pattern is indexed by
 * clicking on it, and everything a user needs while clicking — which pick is
 * the beam, how far each spot is from it, what d-spacing that corresponds to —
 * has to be visible without leaving the image. So the pattern is the plot, the
 * cursor readout reports the distance from the marked beam in both millimetres
 * and inverse angstroms, and every pick is a hoverable entity like any other.
 *
 * Two kinds of pattern arrive on that surface and are treated identically from
 * the first click onward:
 *
 * 1. **An uploaded micrograph.** The image never leaves the browser. Only the
 *    picked coordinates are sent, which keeps an unpublished micrograph on the
 *    machine it was opened on and makes the request small enough to be instant
 *    on a slow intranet.
 * 2. **A practice pattern from the gallery.** Simulated in Python and sent as
 *    coordinates and brightnesses rather than as a raster — a few kilobytes
 *    instead of a megabyte, drawn here as SVG, and crisp at any zoom. Because it
 *    was built from a known zone axis, the answer travels with it, and indexing
 *    can be checked rather than merely performed.
 *
 * The panel is deliberately ordered as the microscope session is: what am I
 * looking at, what is it, where should I go next, how do I get there.
 */

import { el, formatNumber, markdown, svg } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';
import { call } from '../core/api.js';

export const panel = {
  id: 'tem',
  title: 'TEM Solver',
  tagline: 'Open a pattern, index it, and plan the tilt to the next zone axis.',
};

/** Colour of a simulated reflection. A plate is monochrome; so is this. */
const SPOT_COLOUR = '#eaf2ff';

/*
 * Overlay colours are fixed, not theme tokens, and each is drawn over a dark
 * halo.
 *
 * The plate is always dark whatever the interface theme, so `var(--teal)` and
 * `var(--violet)` — which resolve to deep, saturated values in light mode —
 * would put a dark line on a near-black ground for half the users. And an
 * uploaded micrograph may be a light-ground print, where the reverse happens.
 * A bright stroke over a dark one is legible on both, which is the only way to
 * be right about a background this code does not control.
 */
const LATTICE_COLOUR = '#4fd3d3';
const CALCULATED_COLOUR = '#c9b0ff';
const HALO_COLOUR = '#05070d';

export function mount(context) {
  const operations = context.manifest.operations.filter((entry) => entry.panel === panel.id);
  const galleryOperation = operations.find((entry) => entry.id === 'tem.gallery_pattern');
  const solveOperation = operations.find((entry) => entry.id === 'tem.solve_pattern');
  const atlasOperation = operations.find((entry) => entry.id === 'tem.zone_axis_atlas');
  const tiltOperation = operations.find((entry) => entry.id === 'tem.plan_tilt');
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  // The gallery is declared once, in Python, as the options of the operation's
  // `pattern` parameter. Reading it from the manifest means a fourth practice
  // plate appears here the moment it is added there, with no edit in this file.
  const galleryEntries =
    galleryOperation?.parameters.find((parameter) => parameter.name === 'pattern')?.options ?? [];

  const hiddenFields = new Set();

  const state = {
    mode: 'tilt', // 'pick' once a pattern is on the canvas
    source: null, // {kind: 'image'} or {kind: 'gallery'}
    image: null,
    gallery: null,
    picks: { centre: null, spots: [] },
    pickedCentre: null, // where the beam was clicked, before any refinement
    fit: null, // the live 2D lattice fitted to the picks
    fitPending: null,
    showLattice: true,
    showCalculated: true,
    solutions: [],
    selected: 0,
    accepted: null,
    nudgeStep: 1,
    showAnswer: false,
    solveForm: null,
    atlasForm: null,
    tiltForm: null,
    teaches: null,
    solution: null,
  };

  const answerButton = el('button.button', {
    type: 'button',
    text: 'Show answer',
    hidden: true,
    'aria-pressed': 'false',
    title: 'Label every simulated spot with its indices',
    onclick: () => {
      state.showAnswer = !state.showAnswer;
      answerButton.setAttribute('aria-pressed', String(state.showAnswer));
      answerButton.textContent = state.showAnswer ? 'Hide answer' : 'Show answer';
      drawPattern();
    },
  });

  const autoPickButton = el('button.button', {
    type: 'button',
    text: 'Auto-pick',
    hidden: true,
    title: 'Place the beam and six well-separated strong spots',
    onclick: () => autoPick(),
  });

  const latticeButton = el('button.button', {
    type: 'button',
    text: 'Lattice',
    'aria-pressed': 'true',
    title: 'Draw the lattice fitted to the picked spots',
    onclick: () => {
      state.showLattice = !state.showLattice;
      latticeButton.setAttribute('aria-pressed', String(state.showLattice));
      drawPattern();
    },
  });

  const calculatedButton = el('button.button', {
    type: 'button',
    text: 'Calculated',
    hidden: true,
    'aria-pressed': 'true',
    title: 'Superimpose the pattern the selected solution predicts',
    onclick: () => {
      state.showCalculated = !state.showCalculated;
      calculatedButton.setAttribute('aria-pressed', String(state.showCalculated));
      drawPattern();
    },
  });

  const frame = plotFrame({
    title: 'Pattern',
    toolbar: [
      autoPickButton,
      latticeButton,
      calculatedButton,
      answerButton,
      el('button.button', {
        type: 'button',
        text: 'Undo pick',
        onclick: () => {
          if (state.picks.spots.length) state.picks.spots.pop();
          else state.picks.centre = null;
          scheduleFit();
          drawPattern();
        },
      }),
      el('button.button', {
        type: 'button',
        text: 'Clear picks',
        onclick: () => {
          state.picks = { centre: null, spots: [] };
          state.pickedCentre = null;
          state.fit = null;
          state.solutions = [];
          calculatedButton.hidden = true;
          renderCentreTool();
          drawPattern();
        },
      }),
    ],
  });

  const details = el('div');
  context.stage.append(frame.element, details);

  /* ------------------------------------------------------------ controls */

  const fileInput = el('input', {
    type: 'file',
    accept: 'image/*',
    onchange: (event) => loadImage(event.target.files?.[0]),
  });

  const galleryHost = el('div');
  const centreHost = el('div.centre-tool');
  const solveHost = el('div');
  const atlasHost = el('div');
  const tiltHost = el('div');

  const galleryButton = el('button.button.button--block', {
    type: 'button',
    text: 'Reload with these settings',
    onclick: () => loadGallery(state.gallery?.inputs.pattern ?? galleryEntries[0]?.value),
  });
  const solveButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Index the pattern',
    onclick: () => solve(),
  });
  const atlasButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'List the zone axes',
    onclick: () => runAtlas(),
  });
  const tiltButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Plan the tilt',
    onclick: () => planTilt(),
  });

  context.rail.append(
    el('details.group', { open: true }, [
      el('summary', { text: '1 · Open a pattern' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'Pick a practice plate to try the workflow on — each one is a real calculation with ' +
            'the answer attached — or open a micrograph of your own. Your image stays on this ' +
            'machine; only the coordinates you click are sent.',
        }),
        el(
          'div.examples',
          {},
          galleryEntries.map((entry) =>
            el('button.example', { type: 'button', onclick: () => loadGallery(entry.value) }, [
              el('strong', { text: entry.label }),
              el('span', { text: entry.help }),
            ]),
          ),
        ),
        el('p.field__help', { text: 'Or open your own pattern image:' }),
        fileInput,
        galleryHost,
        galleryButton,
      ]),
    ]),
    el('details.group', { open: true }, [
      el('summary', { text: '2 · Calibrate and index' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'Click the transmitted beam first, then the reflections. The beam is not a ' +
            'reflection: it is the origin every spot is measured from, so an error there biases ' +
            'every d-spacing in the pattern.',
        }),
        centreHost,
        solveHost,
        solveButton,
      ]),
    ]),
    el('details.group', { open: true }, [
      el('summary', { text: '3 · Where to go next' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'Every zone axis of this phase near the one on the beam, with how far it is, how ' +
            'much its pattern shows, and whether the holder can reach it. Choose a destination ' +
            'from the list and it becomes the tilt target below.',
        }),
        atlasHost,
        atlasButton,
      ]),
    ]),
    el('details.group', { open: true }, [
      el('summary', { text: '4 · Tilt to the target' }),
      el('div.group__body', {}, [tiltHost, tiltButton]),
    ]),
    el('details.group', {}, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'Complete scenarios you can run without a pattern of your own.',
        }),
        el(
          'div.examples',
          {},
          examples.map((example) =>
            el('button.example', { type: 'button', onclick: () => loadExample(example) }, [
              el('strong', { text: example.title }),
              el('span', { text: example.summary }),
            ]),
          ),
        ),
      ]),
    ]),
  );

  function renderForms({ gallery = {}, index = {}, atlas = {}, tilt = {} } = {}) {
    if (galleryOperation) {
      state.galleryForm = buildForm(galleryOperation, { initial: gallery });
      hideField(state.galleryForm.element, 'pattern');
      galleryHost.replaceChildren(state.galleryForm.element);
    }

    state.solveForm = buildForm(solveOperation, { initial: index });
    // The picks live on the canvas, not in a text box, and the expected axis
    // comes from the gallery entry; both are hidden from the generated form and
    // supplied at call time.
    hideField(state.solveForm.element, 'picks');
    hideField(state.solveForm.element, 'expected_zone_axis');
    solveHost.replaceChildren(state.solveForm.element);

    if (atlasOperation) {
      state.atlasForm = buildForm(atlasOperation, { initial: atlas });
      atlasHost.replaceChildren(state.atlasForm.element);
    }

    state.tiltForm = buildForm(tiltOperation, { initial: tilt });
    tiltHost.replaceChildren(state.tiltForm.element);
  }

  /**
   * Hide a generated field whose value comes from somewhere other than the form.
   *
   * The name is remembered, because a validation error on a hidden field has
   * nowhere to appear. "The transmitted beam has not been marked" is the most
   * likely error this panel produces and its field is the hidden picker: routed
   * to the form it landed on an invisible row, so pressing Index did nothing at
   * all. Errors on these fields go to the toast and the plot status instead.
   */
  function hideField(root, name) {
    for (const field of root.querySelectorAll('.field')) {
      const control = field.querySelector(`[id^="ctl-${name}-"]`);
      if (control) field.hidden = true;
    }
    hiddenFields.add(name);
  }

  /** Report a failure on the form, or visibly if its field is not on screen. */
  function reportError(form, error) {
    if (error?.field && hiddenFields.has(error.field)) {
      form.clearErrors();
      context.showError(error);
      frame.setStatus(error.message);
      return;
    }
    if (!form.showError(error)) context.showError(error);
    else context.showError(error, { quiet: true });
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    if (example.operation === 'tem.gallery_pattern') {
      loadGallery(example.request.pattern, example.request);
      return;
    }
    if (example.operation === 'tem.zone_axis_atlas') {
      renderForms({ atlas: example.request });
      runAtlas();
      return;
    }
    renderForms({ tilt: example.request });
    planTilt();
  }

  /* --------------------------------------------------------- the gallery */

  async function loadGallery(entryId, overrides = null) {
    if (!entryId || !galleryOperation) return;
    galleryButton.disabled = true;
    try {
      const request = { ...(state.galleryForm?.values() ?? {}), ...(overrides ?? {}) };
      request.pattern = entryId;
      const result = await call('tem.gallery_pattern', request);
      state.gallery = result;
      state.source = { kind: 'gallery' };
      state.mode = 'pick';
      state.picks = { centre: null, spots: [] };
      state.pickedCentre = null;
      state.fit = null;
      state.solutions = [];
      state.accepted = null;
      calculatedButton.hidden = true;
      state.image = null;
      state.showAnswer = false;
      answerButton.hidden = false;
      answerButton.textContent = 'Show answer';
      answerButton.setAttribute('aria-pressed', 'false');
      autoPickButton.hidden = false;

      const calibration = result.data.calibration;
      // Carrying the calibration across is the whole point of sending it: a
      // camera constant retyped by hand is a camera constant that can be
      // retyped wrongly, and a pattern indexed against the wrong one indexes
      // to a plausible, self-consistent, wrong material.
      state.solveForm.setValues({
        phase: calibration.phase,
        units: calibration.units,
        camera_constant_mm_angstrom: calibration.camera_constant_mm_angstrom,
        pixel_size_mm: calibration.pixel_size_mm,
      });
      state.atlasForm?.setValues({
        phase: calibration.phase,
        current_zone_axis: result.data.pattern.zone_axis,
      });
      state.tiltForm.setValues({
        phase: calibration.phase,
        current_zone_axis: result.data.pattern.zone_axis,
      });
      if (state.galleryForm) state.galleryForm.setValues({ pattern: entryId });

      drawPattern();
      renderResult(details, result, {
        teaches: state.teaches,
        extra: [galleryGuidance(result)],
      });
      state.teaches = null;
    } catch (error) {
      context.showError(error);
    } finally {
      galleryButton.disabled = false;
    }
  }

  /** What this plate teaches, and where it is worth going from here. */
  function galleryGuidance(result) {
    const targets = result.data.targets ?? [];
    return el('section.card', {}, [
      el('div.card__header', {}, [
        el('h2.card__title', { text: 'What this pattern teaches' }),
        el('p.card__subtitle', {
          text: 'Read this before picking; it is the reasoning the indexing will confirm.',
        }),
      ]),
      el('div.card__body', {}, [
        ...markdown(result.data.entry.teaches),
        targets.length
          ? el('h3.destinations__heading', { text: 'Worth going to from here' })
          : null,
        targets.length
          ? el(
              'div.examples',
              {},
              targets.map((target) =>
                el(
                  'button.example',
                  {
                    type: 'button',
                    onclick: () => chooseTarget(target.indices, target.label),
                  },
                  [el('strong', { text: target.label }), el('span', { text: target.reason })],
                ),
              ),
            )
          : null,
      ]),
    ]);
  }

  /* --------------------------------------------------------------- image */

  function loadImage(file) {
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const image = new Image();
      image.onload = () => {
        const hadGallery = state.gallery !== null;
        state.image = { source: reader.result, width: image.width, height: image.height };
        state.source = { kind: 'image' };
        state.gallery = null;
        state.mode = 'pick';
        state.picks = { centre: null, spots: [] };
        state.pickedCentre = null;
        state.fit = null;
        state.solutions = [];
        state.accepted = null;
        calculatedButton.hidden = true;
        answerButton.hidden = true;
        autoPickButton.hidden = true;
        drawPattern();
        // The calibration fields still hold whatever the last practice plate
        // put there, and a camera constant belonging to a different exposure is
        // the one error this panel cannot detect: it produces a
        // self-consistent, plausible, wrong answer rather than a failure.
        if (hadGallery) {
          frame.setStatus(
            'Your own pattern. Set the camera constant and pixel size for this image before ' +
              'indexing — the fields still hold the practice plate’s calibration, and a ' +
              'camera constant from another exposure indexes a pattern to the wrong material ' +
              'without ever looking wrong.',
          );
        }
      };
      image.src = reader.result;
    };
    reader.readAsDataURL(file);
  }

  /* -------------------------------------------------------------- picking */

  function autoPick() {
    const suggested = state.gallery?.data.suggested_picks;
    if (!suggested) return;
    state.picks = {
      centre: [...suggested.centre],
      spots: suggested.spots.map((spot) => ({ x: spot.x, y: spot.y })),
    };
    state.pickedCentre = null;
    scheduleFit();
    drawPattern();
    frame.setStatus(
      `Beam and ${state.picks.spots.length} spots placed for you — they are the strongest ` +
        'reflections whose directions are all different, which is what indexing needs. Index ' +
        'them, or clear and pick your own.',
    );
  }

  /**
   * Re-fit the 2D lattice to the current picks, and redraw.
   *
   * Debounced rather than immediate: picking is a burst of clicks, and a request
   * per click would put a dozen in flight to answer a question only the last one
   * asks. 200 ms is below the threshold where a user reads the overlay as
   * lagging and above the interval between two deliberate clicks.
   *
   * A failure here is not an error to shout about. Two spots on one row cannot
   * define a lattice, and that is a normal state halfway through picking, so the
   * overlay simply goes away until the picks can support it.
   */
  function scheduleFit() {
    if (state.fitPending) clearTimeout(state.fitPending);
    state.fitPending = setTimeout(() => {
      state.fitPending = null;
      refreshFit();
    }, 200);
  }

  async function refreshFit() {
    const size = frameSize();
    if (!state.picks.centre || state.picks.spots.length < 2 || !size) {
      state.fit = null;
      renderCentreTool();
      drawPattern();
      return;
    }
    try {
      const result = await call('tem.fit_lattice', {
        picks: state.picks,
        frame_width: size.width,
        frame_height: size.height,
      });
      state.fit = result;
    } catch {
      // Not yet a lattice. Say nothing; the picks are mid-flight.
      state.fit = null;
    }
    renderCentreTool();
    drawPattern();
  }

  /** Move the transmitted beam by a whole number of picking units. */
  function nudgeCentre(dx, dy) {
    if (!state.picks.centre) return;
    if (!state.pickedCentre) state.pickedCentre = [...state.picks.centre];
    state.picks.centre = [
      state.picks.centre[0] + dx * state.nudgeStep,
      state.picks.centre[1] + dy * state.nudgeStep,
    ];
    scheduleFit();
    drawPattern();
  }

  function adoptRefinedCentre() {
    if (!state.fit) return;
    if (!state.pickedCentre) state.pickedCentre = [...state.picks.centre];
    state.picks.centre = [...state.fit.data.centre];
    scheduleFit();
    drawPattern();
  }

  function restorePickedCentre() {
    if (!state.pickedCentre) return;
    state.picks.centre = [...state.pickedCentre];
    state.pickedCentre = null;
    scheduleFit();
    drawPattern();
  }

  /**
   * The beam-centre panel: where it is, where the spots say it should be, and
   * the controls to close the gap.
   *
   * The nudge buttons exist because the refinement cannot always be trusted
   * blindly — a centre wrong by an exact lattice vector fits perfectly, and only
   * a person looking at which spot is brightest can settle it. So the tool
   * offers both: solve for it, or move it by hand and watch the overlay.
   */
  function renderCentreTool() {
    if (!state.picks.centre) {
      centreHost.replaceChildren(
        el('p.field__help', {
          text: 'Click the transmitted beam to enable centre refinement.',
        }),
      );
      return;
    }
    const fit = state.fit?.data;
    const shift = fit ? fit.centre_shift : null;
    const arrow = (label, dx, dy, description) =>
      el('button.button.button--small', {
        type: 'button',
        text: label,
        title: description,
        'aria-label': description,
        onclick: () => nudgeCentre(dx, dy),
      });

    centreHost.replaceChildren(
      el('div.field__label', { text: 'Transmitted beam' }),
      el('p.field__help', {
        text:
          'The largest avoidable error in the whole workflow: it biases every d-spacing at once, ' +
          'and leaves the pattern self-consistent while doing it. Nudge it and watch the fitted ' +
          'lattice, or let the spots solve for it.',
      }),
      el('div.centre-tool__readout', {
        text:
          `Now at ${formatNumber(state.picks.centre[0], 1)}, ` +
          `${formatNumber(state.picks.centre[1], 1)}` +
          (shift === null
            ? ' · pick two more spots to fit a lattice'
            : ` · the spots put it ${formatNumber(shift, 1)} away`),
      }),
      el('div.centre-tool__pad', {}, [
        arrow('◀', -1, 0, 'Move the beam left'),
        el('div.centre-tool__column', {}, [
          arrow('▲', 0, -1, 'Move the beam up'),
          arrow('▼', 0, 1, 'Move the beam down'),
        ]),
        arrow('▶', 1, 0, 'Move the beam right'),
        el(
          'select.centre-tool__step',
          {
            'aria-label': 'Nudge step',
            onchange: (event) => {
              state.nudgeStep = Number(event.target.value);
            },
          },
          [0.5, 1, 2, 5, 10].map((step) =>
            el('option', {
              value: String(step),
              text: `${step} px`,
              selected: step === state.nudgeStep,
            }),
          ),
        ),
      ]),
      el('div.button-row', {}, [
        el('button.button.button--small', {
          type: 'button',
          text: 'Refine from the spots',
          disabled: !fit,
          title: 'Move the beam to where the fitted lattice puts it',
          onclick: () => adoptRefinedCentre(),
        }),
        el('button.button.button--small', {
          type: 'button',
          text: 'Undo refinement',
          disabled: !state.pickedCentre,
          title: 'Put the beam back where it was clicked',
          onclick: () => restorePickedCentre(),
        }),
      ]),
      fit && fit.basis_vectors?.length
        ? el('p.field__hint', {
            text:
              'The arrows on the pattern are the two lattice vectors ' +
              fit.basis_vectors
                .map((vector) =>
                  vector.on_a_pick
                    ? `${vector.label} (spot ${vector.spot}, ` +
                      `${formatNumber(vector.length, 1)} px)`
                    : `${vector.label} (no pick on this node)`,
                )
                .join(' and ') +
              '. Every line in the grid turns with them, so those are the two picks worth ' +
              'adjusting first.',
          })
        : null,
      fit && fit.fit.notes.length
        ? el('p.field__hint', { text: fit.fit.notes.join(' ') })
        : null,
      fit && fit.outliers.length
        ? el('p.field__hint', {
            text: `Spot(s) ${fit.outliers.join(', ')} do not sit on the fitted lattice.`,
          })
        : null,
    );
  }

  function calibrationValues() {
    const values = state.solveForm ? state.solveForm.values() : {};
    return {
      units: values.units ?? 'px',
      camera: Number(values.camera_constant_mm_angstrom ?? 180),
      pixel: Number(values.pixel_size_mm ?? 0.05),
    };
  }

  /** Convert a picked coordinate offset into the reciprocal-space radius. */
  function reciprocalRadius(dx, dy) {
    const { units, camera, pixel } = calibrationValues();
    const distance = Math.hypot(dx, dy);
    if (units === 'reciprocal_angstrom') return distance;
    const mm = units === 'px' ? distance * pixel : distance;
    return camera > 0 ? mm / camera : 0;
  }

  function frameSize() {
    if (state.source?.kind === 'gallery') {
      const pattern = state.gallery.data.pattern;
      return { width: pattern.width_px, height: pattern.height_px };
    }
    if (state.image) return { width: state.image.width, height: state.image.height };
    return null;
  }

  function drawPattern() {
    const size = frameSize();
    if (!size) {
      frame.setContent(
        el('div.stage__placeholder', {
          text: 'Choose a practice pattern or open a micrograph to start picking.',
        }),
      );
      frame.setStatus('');
      return;
    }
    const { width, height } = size;
    const root = svg('svg', {
      viewBox: `0 0 ${width} ${height}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': 'Diffraction pattern',
    });

    if (state.source.kind === 'image') {
      root.append(
        svg('image', {
          href: state.image.source,
          x: 0,
          y: 0,
          width,
          height,
          preserveAspectRatio: 'none',
        }),
      );
    } else {
      drawSimulatedPattern(root, state.gallery.data.pattern);
    }

    if (state.showLattice && state.fit) drawFittedLattice(root, width, height);
    if (state.showCalculated) drawCalculatedPattern(root, width, height);

    const marker = Math.max(width, height) / 140;
    if (state.picks.centre) {
      const [cx, cy] = state.picks.centre;
      root.append(
        svg('circle', {
          cx, cy, r: marker * 1.4,
          fill: 'none',
          stroke: 'var(--accent)',
          'stroke-width': marker / 3,
        }),
        svg('line', {
          x1: cx - marker * 2.4, y1: cy, x2: cx + marker * 2.4, y2: cy,
          stroke: 'var(--accent)', 'stroke-width': marker / 4,
        }),
        svg('line', {
          x1: cx, y1: cy - marker * 2.4, x2: cx, y2: cy + marker * 2.4,
          stroke: 'var(--accent)', 'stroke-width': marker / 4,
        }),
      );
    }

    state.picks.spots.forEach((spot, index) => {
      const dx = state.picks.centre ? spot.x - state.picks.centre[0] : 0;
      const dy = state.picks.centre ? spot.y - state.picks.centre[1] : 0;
      const g = reciprocalRadius(dx, dy);
      const node = svg('circle', {
        cx: spot.x, cy: spot.y, r: marker,
        fill: 'none',
        stroke: 'var(--ok, #17683a)',
        'stroke-width': marker / 3,
      });
      root.append(node);
      root.append(
        svg('text', {
          x: spot.x + marker * 1.8, y: spot.y - marker * 0.6,
          'font-size': marker * 2.4,
          fill: 'var(--ok, #17683a)',
          'paint-order': 'stroke',
          stroke: 'var(--bg-raised)',
          'stroke-width': marker / 2,
          text: String(index + 1),
        }),
      );
      frame.hoverable(node, {
        Spot: index + 1,
        x: spot.x,
        y: spot.y,
        'Distance from beam': Math.hypot(dx, dy),
        '|g| / Å⁻¹': g,
        'd / Å': g > 0 ? 1 / g : null,
      });
    });

    root.addEventListener('click', (event) => {
      const point = eventToImage(event, root, width, height);
      if (!point) return;
      if (!state.picks.centre) {
        state.picks.centre = [point.x, point.y];
        state.pickedCentre = null;
      } else state.picks.spots.push({ x: point.x, y: point.y });
      scheduleFit();
      drawPattern();
    });

    frame.configure({
      toData: (x, y) => ({ x, y }),
      formatCursor: (point) => {
        if (!state.picks.centre) return `${formatNumber(point.x, 0)}, ${formatNumber(point.y, 0)} px`;
        const dx = point.x - state.picks.centre[0];
        const dy = point.y - state.picks.centre[1];
        const g = reciprocalRadius(dx, dy);
        return (
          `${formatNumber(Math.hypot(dx, dy), 1)} px from beam · ` +
          `|g| ${formatNumber(g, 4)} Å⁻¹ · d ${g > 0 ? formatNumber(1 / g, 4) : '∞'} Å`
        );
      },
    });
    frame.setContent(root);
    frame.setStatus(
      state.picks.centre
        ? `Beam marked · ${state.picks.spots.length} spot(s) picked · click to add more`
        : 'Click the transmitted beam first — it is the origin every spot is measured from',
    );
  }

  /**
   * Draw a simulated pattern: dark ground, glowing spots, optional indices.
   *
   * A plate is dark with bright spots in every theme, because that is what a
   * diffraction pattern is; inverting it for a light theme would make the
   * practice pattern look unlike the thing it is practice for. Each reflection
   * is a soft radial glow with a solid core, scaled by the apparent radius the
   * simulation reports, so a strong reflection reads as strong at a glance in
   * the way it does on a real exposure.
   */
  function drawSimulatedPattern(root, pattern) {
    const [cx, cy] = pattern.centre_px;
    const gradientId = 'saed-glow';
    root.append(
      svg('defs', {}, [
        svg('radialGradient', { id: gradientId }, [
          svg('stop', { offset: '0%', 'stop-color': SPOT_COLOUR, 'stop-opacity': '0.95' }),
          svg('stop', { offset: '45%', 'stop-color': SPOT_COLOUR, 'stop-opacity': '0.28' }),
          svg('stop', { offset: '100%', 'stop-color': SPOT_COLOUR, 'stop-opacity': '0' }),
        ]),
      ]),
      svg('rect', {
        x: 0, y: 0, width: pattern.width_px, height: pattern.height_px, fill: '#05070d',
      }),
    );

    // The transmitted beam, drawn unmistakably brighter and larger than any
    // reflection. This is not decoration: the first instruction the panel gives
    // is "click the transmitted beam", and in an earlier draft the beam was
    // barely distinguishable from a strong 200 spot, which makes that
    // instruction impossible to follow. On a real plate the direct beam is
    // orders of magnitude brighter — that is why a beam stop exists.
    //
    // Its size is tied to the nearest reflection rather than fixed, so a dense
    // hexagonal zone does not have its innermost spots swallowed by the glow.
    const nearest = pattern.spots.length
      ? Math.min(...pattern.spots.map((spot) => Math.hypot(spot.x - cx, spot.y - cy)))
      : pattern.width_px / 8;
    const beamGlow = Math.max(28, Math.min(0.5 * nearest, pattern.width_px / 8));
    root.append(
      svg('circle', { cx, cy, r: beamGlow, fill: `url(#${gradientId})`, opacity: 0.95 }),
      svg('circle', { cx, cy, r: beamGlow * 0.55, fill: `url(#${gradientId})`, opacity: 0.95 }),
      svg('circle', { cx, cy, r: beamGlow * 0.2, fill: '#ffffff', opacity: 0.95 }),
    );

    const marker = Math.max(pattern.width_px, pattern.height_px) / 140;
    for (const spot of pattern.spots) {
      const glow = svg('circle', {
        cx: spot.x, cy: spot.y, r: spot.radius_px * 3.0,
        fill: `url(#${gradientId})`,
        // Opacity carries the intensity, because apparent radius alone barely
        // separates reflections whose kinematic intensities differ by tens of
        // percent — which within one zone is the usual case.
        opacity: 0.25 + 0.7 * spot.intensity,
      });
      const core = svg('circle', {
        cx: spot.x, cy: spot.y, r: spot.radius_px * 0.62,
        fill: SPOT_COLOUR,
        opacity: 0.45 + 0.55 * spot.intensity,
      });
      root.append(glow, core);
      frame.hoverable(core, {
        Index: spot.label,
        'd / Å': spot.d_angstrom,
        '|g| / Å⁻¹': spot.g_inv_angstrom,
        'Relative intensity': spot.intensity,
        x: spot.x,
        y: spot.y,
      });
      if (state.showAnswer) {
        root.append(
          svg('text', {
            x: spot.x + spot.radius_px * 1.6,
            y: spot.y - spot.radius_px * 1.2,
            'font-size': marker * 2.1,
            fill: '#9fd2ff',
            text: spot.label,
          }),
        );
      }
    }

    drawScaleBar(root, pattern, marker);
  }

  /**
   * A reciprocal-space scale bar, in inverse angstroms.
   *
   * A diffraction pattern without a scale is a picture. The bar makes the
   * calibration visible rather than merely stored: change the camera length and
   * the bar changes with the pattern, which is the relation the calibration
   * field is trying to teach. Its length is chosen from a 1-2-5 sequence so the
   * label is a round number and the bar lands between a sixth and a third of the
   * frame — long enough to read against, short enough not to cross the pattern.
   */
  function drawScaleBar(root, pattern, marker) {
    const pxPerInvAngstrom = pattern.camera_constant_mm_angstrom / pattern.pixel_size_mm;
    if (!Number.isFinite(pxPerInvAngstrom) || pxPerInvAngstrom <= 0) return;
    const target = pattern.width_px * 0.22;
    const choices = [0.05, 0.1, 0.2, 0.5, 1, 2, 5];
    let best = choices[0];
    for (const value of choices) {
      if (Math.abs(value * pxPerInvAngstrom - target) < Math.abs(best * pxPerInvAngstrom - target)) {
        best = value;
      }
    }
    const length = best * pxPerInvAngstrom;
    const x = pattern.width_px * 0.06;
    const y = pattern.height_px * 0.94;
    root.append(
      svg('line', {
        x1: x, y1: y, x2: x + length, y2: y,
        stroke: SPOT_COLOUR, 'stroke-width': marker * 0.5, 'stroke-opacity': 0.75,
      }),
      svg('text', {
        x, y: y - marker * 1.1,
        'font-size': marker * 2.2,
        fill: SPOT_COLOUR,
        'fill-opacity': 0.75,
        text: `${best} Å⁻¹`,
      }),
    );
  }

  /**
   * The fitted lattice, as two families of lines through the refined centre.
   *
   * Lines rather than dots. A grid of points is a second set of spots to confuse
   * with the first; ruled lines are unmistakably an overlay, and they make the
   * two things a user is checking visible at a glance — whether the rows pass
   * through the spots, and whether their intersection sits on the beam.
   */
  function drawFittedLattice(root, width, height) {
    const data = state.fit.data;
    const [ax, ay] = data.fit.basis[0];
    const [bx, by] = data.fit.basis[1];
    const [cx, cy] = data.centre;
    const reach = Math.ceil(Math.hypot(width, height) / Math.min(
      Math.hypot(ax, ay), Math.hypot(bx, by),
    ));
    const lines = [];
    for (let k = -reach; k <= reach; k += 1) {
      lines.push([
        cx + k * bx - reach * ax, cy + k * by - reach * ay,
        cx + k * bx + reach * ax, cy + k * by + reach * ay,
      ]);
      lines.push([
        cx + k * ax - reach * bx, cy + k * ay - reach * by,
        cx + k * ax + reach * bx, cy + k * ay + reach * by,
      ]);
    }
    const group = svg('g', { 'pointer-events': 'none' });
    const weight = Math.max(width, height) / 900;
    for (const [x1, y1, x2, y2] of lines) {
      group.append(
        svg('line', {
          x1, y1, x2, y2,
          stroke: HALO_COLOUR,
          'stroke-opacity': 0.55,
          'stroke-width': weight * 3,
        }),
      );
    }
    for (const [x1, y1, x2, y2] of lines) {
      group.append(
        svg('line', {
          x1, y1, x2, y2,
          stroke: LATTICE_COLOUR,
          'stroke-opacity': 0.55,
          'stroke-width': weight,
        }),
      );
    }
    // The origin of the fitted lattice, which is where the spots say the beam is.
    const marker = Math.max(width, height) / 110;
    group.append(
      svg('circle', {
        cx, cy, r: marker,
        fill: 'none',
        stroke: LATTICE_COLOUR,
        'stroke-width': marker / 4,
        'stroke-dasharray': `${marker / 2} ${marker / 2}`,
      }),
    );
    root.append(group);
    drawBasisVectors(root, width, height);
  }

  /**
   * The two basis vectors, as labelled arrows from the beam to the spots that
   * generate them.
   *
   * The grid alone shows *that* the picks are consistent; the arrows show
   * *which two* picks are carrying the whole lattice. That matters while
   * adjusting: move the spot an arrow points at and every line in the grid
   * turns with it, so the two picks worth being careful about are the two the
   * arrows are on. Each arrow ends on the picked spot rather than on the ideal
   * node, so the gap between the head and the node is the error in that pick,
   * visible without reading anything.
   */
  function drawBasisVectors(root, width, height) {
    const vectors = state.fit?.data.basis_vectors ?? [];
    if (!vectors.length) return;
    const marker = Math.max(width, height) / 110;
    const group = svg('g', { 'pointer-events': 'none' });
    for (const vector of vectors) {
      const [x1, y1] = vector.from;
      const [x2, y2] = vector.to;
      const dx = x2 - x1;
      const dy = y2 - y1;
      const length = Math.hypot(dx, dy);
      if (!length) continue;
      const ux = dx / length;
      const uy = dy / length;
      // Stop the shaft short of the spot so the arrowhead sits beside it rather
      // than on top of the thing being measured.
      const tipX = x2 - ux * marker * 1.4;
      const tipY = y2 - uy * marker * 1.4;
      // Generous, because the head is the only part of the arrow that says
      // which end is which, and at marker * 1.5 it was lost against the spot.
      const head = marker * 2.6;
      const wing = marker * 1.1;
      const points = [
        `${tipX},${tipY}`,
        `${tipX - ux * head + -uy * wing},${tipY - uy * head + ux * wing}`,
        `${tipX - ux * head - -uy * wing},${tipY - uy * head - ux * wing}`,
      ].join(' ');
      for (const [colour, weight, opacity] of [
        [HALO_COLOUR, marker / 2, 0.6],
        [LATTICE_COLOUR, marker / 5, 1],
      ]) {
        group.append(
          svg('line', {
            x1, y1, x2: tipX, y2: tipY,
            stroke: colour,
            'stroke-width': weight,
            'stroke-opacity': opacity,
            'stroke-linecap': 'round',
          }),
        );
      }
      group.append(svg('polygon', { points, fill: LATTICE_COLOUR }));
      group.append(
        svg('text', {
          x: x1 + dx * 0.55 - uy * marker * 1.6,
          y: y1 + dy * 0.55 + ux * marker * 1.6,
          'font-size': marker * 2.2,
          'font-style': 'italic',
          fill: LATTICE_COLOUR,
          'paint-order': 'stroke',
          stroke: HALO_COLOUR,
          'stroke-width': marker / 2,
          'text-anchor': 'middle',
          text: vector.label,
        }),
      );
    }
    root.append(group);
  }

  /**
   * The pattern a chosen solution predicts, drawn over the one measured.
   *
   * Open rings, so a measured spot showing through the middle of one reads as
   * agreement. Every reflection of the zone is drawn, not only those a spot was
   * picked for: a predicted ring with nothing inside it is either a spot not yet
   * picked or evidence against the candidate, and both are worth seeing.
   */
  function drawCalculatedPattern(root, width, height) {
    const solution = state.solutions[state.selected];
    if (!solution?.overlay?.length) return;
    // The ring must be wider than the spot it is claiming to explain. Drawn at
    // half this size it sat *inside* the bright core and read as part of the
    // spot rather than as a prediction about it, which is the one thing a
    // superimposed pattern must never do.
    const marker = Math.max(width, height) / 90;
    const group = svg('g');
    for (const spot of solution.overlay) {
      group.append(
        svg('circle', {
          cx: spot.x, cy: spot.y, r: marker,
          fill: 'none',
          stroke: HALO_COLOUR,
          'stroke-width': marker / 3,
          'stroke-opacity': 0.5,
        }),
      );
      const node = svg('circle', {
        cx: spot.x, cy: spot.y, r: marker,
        fill: 'none',
        stroke: CALCULATED_COLOUR,
        'stroke-width': marker / 7,
        'stroke-opacity': 0.95,
      });
      group.append(node);
      frame.hoverable(node, {
        Calculated: spot.label,
        'd / Å': spot.d,
        '|g| / Å⁻¹': spot.g,
        x: spot.x,
        y: spot.y,
      });
    }
    root.append(group);
  }

  function eventToImage(event, node, width, height) {
    const rect = node.getBoundingClientRect();
    if (!rect.width || !rect.height) return null;
    const scale = Math.min(rect.width / width, rect.height / height);
    const offsetX = (rect.width - width * scale) / 2;
    const offsetY = (rect.height - height * scale) / 2;
    return {
      x: (event.clientX - rect.left - offsetX) / scale,
      y: (event.clientY - rect.top - offsetY) / scale,
    };
  }

  /* -------------------------------------------------------------- actions */

  async function solve() {
    solveButton.disabled = true;
    solveButton.textContent = 'Indexing…';
    state.solveForm.clearErrors();
    try {
      const request = { ...state.solveForm.values(), picks: state.picks };
      // A practice plate knows its own answer, so the indexing is checked
      // rather than merely reported. The comparison happens in Python, where
      // the symmetry lives: [110] and [101] are the same bcc pattern.
      if (state.source?.kind === 'gallery') {
        request.expected_zone_axis = state.gallery.data.pattern.zone_axis;
      }
      const result = await call('tem.solve_pattern', request);
      state.solution = result;
      state.solutions = result.data.alternatives ?? [];
      state.selected = 0;
      state.accepted = null;
      calculatedButton.hidden = state.solutions.length === 0;
      state.showCalculated = true;
      calculatedButton.setAttribute('aria-pressed', 'true');
      renderResult(details, result, {
        extra: [verdictCard(result), solutionsCard(result), deviationCard(result)].filter(Boolean),
      });
      drawPattern();
      // Carry the answer into the panels below: the axis just indexed is the
      // axis the next move starts from, and retyping it is both tedious and a
      // chance to get it wrong.
      frame.setStatus(
        `Indexed as ${result.data.phase_name} down ${result.data.zone_axis_label}` +
          ` · score ${formatNumber(result.data.score.score, 3)}` +
          ` · ${state.picks.spots.length} picks · accept a solution to plan from it`,
      );
    } catch (error) {
      // The verdict answers one attempt. Leaving the previous "Correct" card on
      // screen beside a failed index reads as if the failure were the correct
      // answer, so it goes; the rest of the previous result stays, because
      // losing it to a mistyped camera constant would be worse.
      for (const card of details.querySelectorAll('.card--verdict')) card.remove();
      reportError(state.solveForm, error);
    } finally {
      solveButton.disabled = false;
      solveButton.textContent = 'Index the pattern';
    }
  }

  /** The answer check, when the pattern came with an answer. */
  function verdictCard(result) {
    const check = result.data.check;
    if (!check) return null;
    return el('section.card.card--verdict', {
      class: check.correct ? 'card--verdict-correct' : 'card--verdict-wrong',
    }, [
      el('div.card__header', {}, [
        el('h2.card__title', {
          text: check.correct ? 'Correct — that is the axis' : 'Not the expected axis',
        }),
      ]),
      el('div.card__body', {}, [
        el('p.summary', {
          text: check.correct
            ? `You indexed this pattern as ${result.data.zone_axis_label}, which is the axis it ` +
              `was built from, to within ${formatNumber(check.deviation_deg, 3)}°. The check is ` +
              'made up to symmetry, because symmetry-equivalent axes give identical patterns.'
            : `You indexed this pattern as ${result.data.zone_axis_label}, but it was built from ` +
              `${check.expected_label} — ${formatNumber(check.deviation_deg, 2)}° away. Check the ` +
              'beam pick first: an error there biases every d-spacing at once. Then the camera ' +
              'constant, then whether every spot you picked really is a reflection.',
        }),
      ]),
    ]);
  }

  /**
   * The candidates, ranked by the score, each one selectable and acceptable.
   *
   * Selecting draws that candidate's calculated pattern over the measured one,
   * which is how the choice should be made; accepting is a separate, deliberate
   * act that carries the answer into the panels below. Keeping them apart means
   * looking at a candidate costs nothing and commits to nothing.
   */
  function solutionsCard(result) {
    const alternatives = result.data.alternatives ?? [];
    if (!alternatives.length) return null;
    const list = el('div.examples');

    const paint = () => {
      list.replaceChildren(
        ...alternatives.map((entry, index) =>
          el(
            'button.example',
            {
              type: 'button',
              class: index === state.selected ? 'example--selected' : null,
              'aria-pressed': String(index === state.selected),
              onclick: () => {
                state.selected = index;
                paint();
                drawPattern();
              },
            },
            [
              el('strong', {
                text: `${entry.phase} ${entry.zone_axis} · score ${formatNumber(entry.score, 3)}`,
              }),
              el('span', {
                text:
                  `${entry.matched_spots} of ${state.picks.spots.length} spots indexed · ` +
                  `d to ${formatNumber(100 * entry.rms_relative_length_deviation, 2)}% · ` +
                  `angles to ${formatNumber(entry.rms_angle_deviation_deg, 2)}°`,
              }),
              scoreBar(entry),
            ],
          ),
        ),
      );
    };
    paint();

    return el('section.card', {}, [
      el('div.card__header', {}, [
        el('h2.card__title', { text: 'Candidate solutions' }),
        el('p.card__subtitle', {
          text:
            'Ranked by the accuracy score. Choosing one draws the pattern it predicts over the ' +
            'one you measured — the honest way to decide is to look at whether they coincide.',
        }),
        el('div.button-row', { style: 'margin-left:auto' }, [
          el('button.button.button--primary', {
            type: 'button',
            text: 'Accept this solution',
            onclick: () => acceptSolution(result),
          }),
        ]),
      ]),
      el('div.card__body', {}, [list]),
    ]);
  }

  /** Three bars: how well lengths, angles and coverage each agree. */
  function scoreBar(entry) {
    const term = (label, value, colour) =>
      el('span.score-term', { title: `${label}: ${formatNumber(value, 3)} of 1` }, [
        el('span.score-term__label', { text: label }),
        el('span.score-term__track', {}, [
          el('span.score-term__fill', {
            style: `width:${Math.max(0, Math.min(1, value)) * 100}%;background:${colour}`,
          }),
        ]),
      ]);
    return el('span.score-bars', {}, [
      term('d', entry.length_agreement, 'var(--accent)'),
      term('angle', entry.angle_agreement, 'var(--teal)'),
      term('spots', entry.coverage_agreement, 'var(--violet)'),
    ]);
  }

  /** The measured-against-calculated evidence, spot by spot and pair by pair. */
  function deviationCard(result) {
    const score = result.data.score;
    if (!score) return null;
    const angles = score.angle_deviations.slice(0, 10);
    return el('section.card', {}, [
      el('div.card__header', {}, [
        el('h2.card__title', { text: 'Where it disagrees' }),
        el('p.card__subtitle', {
          text:
            'The same deviation on every spot is a camera constant; a scatter of them is an ' +
            'indexing error. Angles do not depend on the calibration at all, so an angular ' +
            'disagreement is evidence about the phase.',
        }),
      ]),
      el('div.card__body', {}, [
        el('p.summary', { text: score.describe }),
        angles.length
          ? el('table.table--compact', {}, [
              el('thead', {}, [
                el('tr', {}, [
                  el('th', { text: 'Pair' }),
                  el('th.numeric', { text: 'Measured / °' }),
                  el('th.numeric', { text: 'Calculated / °' }),
                  el('th.numeric', { text: 'Δ / °' }),
                ]),
              ]),
              el(
                'tbody',
                {},
                angles.map((pair) =>
                  el('tr', {}, [
                    el('td', { text: pair.pair }),
                    el('td.numeric', { text: formatNumber(pair.measured_deg, 2) }),
                    el('td.numeric', { text: formatNumber(pair.calculated_deg, 2) }),
                    el('td.numeric', { text: formatNumber(pair.deviation_deg, 2) }),
                  ]),
                ),
              ),
            ])
          : null,
      ]),
    ]);
  }

  /**
   * Accept a candidate: it becomes the orientation everything downstream uses.
   *
   * Explicit rather than automatic. Indexing produces candidates; deciding which
   * one is right is the user's judgement, and a tilt planned from a solution
   * nobody chose is a tilt planned from a guess.
   */
  function acceptSolution(result) {
    const entry = state.solutions[state.selected];
    if (!entry) return;
    state.accepted = entry;
    // The phase carried forward is the value the picker emitted, not the
    // expanded description the result echoes back: a catalogue choice that makes
    // a round trip through `PhaseSpec.to_json` comes back as a full description,
    // which the picker then labels "(edited)" — untrue, and alarming on a phase
    // nobody edited.
    const phase = state.solveForm.values().phase;
    const axis = entry.zone_axis_indices ?? result.data.zone_axis;
    state.atlasForm?.setValues({ phase, current_zone_axis: axis });
    state.tiltForm.setValues({ phase, current_zone_axis: axis });
    frame.setStatus(
      `Accepted ${entry.phase} down ${entry.zone_axis} (score ` +
        `${formatNumber(entry.score, 3)}). It is now the starting orientation for the zone-axis ` +
        'list and the tilt plan below.',
    );
    atlasButton.scrollIntoView({ behavior: 'smooth', block: 'center' });
  }

  async function runAtlas() {
    if (!state.atlasForm) return;
    atlasButton.disabled = true;
    atlasButton.textContent = 'Searching…';
    state.atlasForm.clearErrors();
    try {
      const result = await call('tem.zone_axis_atlas', state.atlasForm.values());
      renderResult(details, result, {
        teaches: state.teaches,
        extra: [destinationsCard(result)],
      });
      state.teaches = null;
    } catch (error) {
      reportError(state.atlasForm, error);
    } finally {
      atlasButton.disabled = false;
      atlasButton.textContent = 'List the zone axes';
    }
  }

  /**
   * The atlas rows as destinations you can act on.
   *
   * The table below states the numbers; this turns each row into the decision
   * it represents. Reading "⟨111⟩ at 35.26°, 36 reflections, six-fold,
   * reachable" and then retyping [111] into a form two sections down is exactly
   * the kind of transcription that puts the wrong indices into a tilt plan.
   */
  function destinationsCard(result) {
    const rows = result.table.rows.filter((row) => row.verdict !== 'current axis');
    return el('section.card', {}, [
      el('div.card__header', {}, [
        el('h2.card__title', { text: 'Choose a destination' }),
        el('p.card__subtitle', {
          text:
            'Nearest first. Choosing one sets it as the tilt target below — the number of ' +
            'reflections is what the trip buys, and the symmetry is what you will recognise ' +
            'when you arrive.',
        }),
      ]),
      el('div.card__body', {}, [
        el(
          'div.examples',
          {},
          rows.map((row) =>
            el(
              'button.example',
              {
                type: 'button',
                class: row.reachable ? 'example--reachable' : 'example--unreachable',
                onclick: () => chooseTarget(row.indices, row.target),
              },
              [
                el('strong', { text: `${row.family} · ${row.target}` }),
                el('span', {
                  text:
                    `${formatNumber(row.angle_deg, 2)}° away · ${row.reflections} reflections · ` +
                    `${row.symmetry} · ${row.family_size} members · ` +
                    (row.reachable
                      ? `reachable, Δα ${formatNumber(row.delta_alpha_deg, 1)}° ` +
                        `Δβ ${formatNumber(row.delta_beta_deg, 1)}°`
                      : 'out of the holder’s range in one move'),
                }),
              ],
            ),
          ),
        ),
      ]),
    ]);
  }

  function chooseTarget(indices, label) {
    if (!indices) return;
    state.tiltForm.setValues({ target_zone_axis: indices });
    frame.setStatus(`Target set to ${label ?? indices.join(' ')} — plan the tilt below.`);
    tiltButton.scrollIntoView({ behavior: 'smooth', block: 'center' });
    tiltButton.focus();
  }

  async function planTilt() {
    tiltButton.disabled = true;
    tiltButton.textContent = 'Planning…';
    state.tiltForm.clearErrors();
    try {
      const result = await call('tem.plan_tilt', state.tiltForm.values());
      renderResult(details, result, {
        teaches: state.teaches,
        extra: [tiltMap(result)],
      });
      state.teaches = null;
    } catch (error) {
      reportError(state.tiltForm, error);
    } finally {
      tiltButton.disabled = false;
      tiltButton.textContent = 'Plan the tilt';
    }
  }

  renderForms();
  renderCentreTool();
  drawPattern();
  if (galleryEntries.length) loadGallery(galleryEntries[0].value);

  return {
    help: () => {
      if (state.mode === 'pick') return solveOperation;
      return atlasOperation ?? tiltOperation;
    },
  };
}

/**
 * The tilt map: candidate destinations drawn in the holder's own coordinates.
 *
 * Alpha across, beta up, the envelope as a rectangle. A table of angles says
 * where each solution is; this says at a glance which ones are comfortably
 * inside the envelope and which are pressed against a stop, which is the
 * difference between a move that works and one that drifts off axis.
 */
function tiltMap(result) {
  const { envelope, start } = result.data;
  const alphaLimit = envelope.alpha_limit_deg;
  const betaLimit = envelope.beta_limit_deg;
  const pad = 1.15;
  const scale = 100 / (Math.max(alphaLimit, betaLimit) * pad);

  const root = svg('svg', {
    viewBox: '-110 -110 220 220',
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': 'Tilt map',
    style: 'max-height:22rem',
  });

  root.append(
    svg('rect', {
      x: -alphaLimit * scale, y: -betaLimit * scale,
      width: 2 * alphaLimit * scale, height: 2 * betaLimit * scale,
      fill: 'none', stroke: 'currentColor', 'stroke-opacity': 0.5, 'stroke-width': 0.8,
      'stroke-dasharray': '3 2',
    }),
    svg('line', { x1: -105, y1: 0, x2: 105, y2: 0, stroke: 'currentColor', 'stroke-opacity': 0.25, 'stroke-width': 0.4 }),
    svg('line', { x1: 0, y1: -105, x2: 0, y2: 105, stroke: 'currentColor', 'stroke-opacity': 0.25, 'stroke-width': 0.4 }),
    svg('text', { x: alphaLimit * scale, y: 8, 'font-size': 5, fill: 'currentColor', 'fill-opacity': 0.6, 'text-anchor': 'end', text: `α ${alphaLimit}°` }),
    svg('text', { x: 2, y: -betaLimit * scale - 2, 'font-size': 5, fill: 'currentColor', 'fill-opacity': 0.6, text: `β ${betaLimit}°` }),
  );

  for (const row of result.table.rows) {
    const x = row.alpha_deg * scale;
    const y = -row.beta_deg * scale;
    const reachable = row.margin_deg >= 0;
    root.append(
      svg('line', {
        x1: start.alpha_deg * scale, y1: -start.beta_deg * scale, x2: x, y2: y,
        stroke: reachable ? 'var(--accent)' : 'var(--danger)',
        'stroke-opacity': 0.35,
        'stroke-width': 0.5,
      }),
      svg('circle', {
        cx: x, cy: y, r: 2.4,
        fill: reachable ? 'var(--accent)' : 'var(--danger)',
        'fill-opacity': reachable ? 0.9 : 0.5,
      }),
      svg('text', {
        x: x + 3.5, y: y + 1.6,
        'font-size': 4.5,
        fill: 'currentColor',
        text: row.member,
      }),
    );
  }

  root.append(
    svg('circle', {
      cx: start.alpha_deg * scale, cy: -start.beta_deg * scale, r: 3,
      fill: 'none', stroke: 'currentColor', 'stroke-width': 0.9,
    }),
  );

  return el('section.card', {}, [
    el('div.card__header', {}, [
      el('h2.card__title', { text: 'Tilt map' }),
      el('p.card__subtitle', {
        text: 'Alpha across, beta up. The dashed rectangle is the holder; open circle is where you are.',
      }),
    ]),
    el('div.card__body', {}, [root]),
  ]);
}
