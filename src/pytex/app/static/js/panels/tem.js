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

import { append, el, formatNumber, markdown, svg } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';
import { call } from '../core/api.js';
import * as log from '../core/logbook.js';

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
// Kikuchi bands, and the one band that is a route rather than a landmark.
//
// The route is deliberately not the amber of the basis vectors below: they are
// drawn in the same region of the plate, and two different claims about the
// crystal must not share a colour.
const KIKUCHI_COLOUR = '#8fd9ff';
const KIKUCHI_ROUTE_COLOUR = '#a3e635';
const HALO_COLOUR = '#05070d';

/*
 * Four things are drawn over the plate and each means something different, so
 * each gets its own colour. They were not all distinct before, and the picture
 * could not be read: the fitted grid and the two basis vectors shared the teal,
 * so an arrow into empty space looked like part of the grid rather than like the
 * claim it is — "this pick generates this lattice".
 *
 * What is clicked is warm; what is computed is cool. Picks and the beam are
 * amber and white — the user's own marks. The fitted grid is teal and the basis
 * arrows gold, both derived. The calculated pattern stays violet.
 */
const PICK_COLOUR = '#ffc46b';
const BEAM_COLOUR = '#ffffff';
const BASIS_COLOUR = '#ffd447';
/** The centre the fit prefers, when it is not the one that was clicked. */
const REFINED_COLOUR = '#7ce4a8';
/** The pick under the pointer or selected in the coordinate table. */
const SELECTED_COLOUR = '#ff7ad9';

/**
 * A serial number for the pattern now on the plate.
 *
 * Opening the same file twice is a deliberate act — the file may have been
 * re-exported — so identity by name alone would treat the second open as a
 * redraw and keep the first one's camera. The serial makes every open distinct.
 */
let patternSerial = 0;
function nextPatternSerial() {
  patternSerial += 1;
  return patternSerial;
}

export function mount(context) {
  const operations = context.manifest.operations.filter((entry) => entry.panel === panel.id);
  const galleryOperation = operations.find((entry) => entry.id === 'tem.gallery_pattern');
  const solveOperation = operations.find((entry) => entry.id === 'tem.solve_pattern');
  const atlasOperation = operations.find((entry) => entry.id === 'tem.zone_axis_atlas');
  const tiltOperation = operations.find((entry) => entry.id === 'tem.plan_tilt');
  const stereogramOperation = operations.find((entry) => entry.id === 'tem.stereogram');
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
    // Which pick the coordinate table is working on: 'centre', or a spot index.
    // The nudge pad and the arrow keys act on it, and it is drawn ringed. Named
    // apart from `selected` below, which is the chosen *solution*: one object,
    // two selections, and sharing the key silently broke the solution picker.
    selectedPick: 'centre',
    // Which pick the pointer is over, in the table or on the plate.
    highlight: null,
    // The live coordinate inputs, so a committed edit can refresh the derived
    // columns without rebuilding the table under the user's cursor.
    pickRows: [],
    nudgeLabel: null,
    // The uploaded micrograph's pixels, for snapping a click to a spot centroid.
    pixels: null,
    // Which pattern is on the plate, and which one the frame was last drawn
    // for. Equal means a redraw of the same picture, which keeps its camera;
    // different means a pattern was just opened, which gets a fresh Fit.
    patternKey: null,
    drawnKey: null,
    showLattice: true,
    showCalculated: true,
    // The Kikuchi bands the accepted solution predicts. There is nothing to
    // draw before a solution is accepted, so the toggle appears with it.
    showKikuchi: false,
    kikuchi: null,
    kikuchiRequest: null,
    kikuchiPending: null,
    // The calibration tool: whether it is taking clicks, the two ends of the
    // measured line, and the length the user says it is.
    calibrate: { active: false, points: [], length: '', unit: 'inv_angstrom' },
    // The stereogram beside the pattern: its last result, the request that
    // produced it, and the two things it can be asked to stop drawing.
    stereo: null,
    stereoRequest: null,
    stereoPending: null,
    showEnvelope: true,
    showPoleLabels: true,
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

  /*
   * Calibrate: measure a known length on the image itself.
   *
   * An image that arrives by email has a scale bar and no recorded camera
   * length, and the camera equation only ever uses one number — the reciprocal
   * space one pixel spans. Drawing a line across something whose length is known
   * measures that number directly. Two known lengths are useful and they answer
   * different questions, so the tool asks which: a reciprocal length (a scale
   * bar in Å⁻¹, or a reflection whose spacing is known) fixes the scale itself;
   * a real length on the plate (cm on a print, mm on the detector) fixes the
   * pixel size, and the camera constant then does the rest as before.
   */
  const calibrateButton = el('button.button', {
    type: 'button',
    text: 'Calibrate',
    'aria-pressed': 'false',
    title: 'Measure a known length on the image to set the scale',
    onclick: () => setCalibrating(!state.calibrate.active),
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

  /*
   * Kikuchi: the bands the accepted orientation predicts, on the same plate.
   *
   * A detector records directions of outgoing electrons, and the spots and the
   * bands are placed in that one angular space by the same reciprocal lattice
   * and the same orientation, so superimposing them mixes nothing. The metrics
   * agree as well: a band is exactly as wide as its own 000→g distance and
   * perpendicular to it, which is a check the user can make by eye — and the
   * reason nothing here needs the diffraction rotation or the parity, neither
   * of which one indexed pattern determines.
   */
  const kikuchiButton = el('button.button', {
    type: 'button',
    text: 'Kikuchi',
    hidden: true,
    'aria-pressed': 'false',
    title: 'Superimpose the Kikuchi bands the accepted solution predicts',
    onclick: () => {
      state.showKikuchi = !state.showKikuchi;
      kikuchiButton.setAttribute('aria-pressed', String(state.showKikuchi));
      if (state.showKikuchi) refreshKikuchi();
      else drawPattern();
    },
  });

  const frame = plotFrame({
    title: 'Pattern',
    toolbar: [
      autoPickButton,
      calibrateButton,
      latticeButton,
      calculatedButton,
      kikuchiButton,
      answerButton,
      el('button.button', {
        type: 'button',
        text: 'Undo pick',
        onclick: () => {
          if (state.picks.spots.length) state.picks.spots.pop();
          else state.picks.centre = null;
          state.selectedPick = state.picks.spots.length ? state.picks.spots.length - 1 : 'centre';
          state.highlight = null;
          renderPickTool();
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
          state.selectedPick = 'centre';
          state.highlight = null;
          calculatedButton.hidden = true;
          resetKikuchi();
          renderPickTool();
          drawPattern();
        },
      }),
    ],
  });

  /*
   * The stereogram, beside the pattern rather than under it.
   *
   * A diffraction pattern is square and the stage is wide, so the pattern alone
   * left half the workspace blank while the question it raises — *where do I go
   * next, and can the holder get there* — was answered only by a table further
   * down the page. The two belong side by side: the pattern says what is on the
   * beam, the stereogram says what else is within reach and in which direction.
   */
  const stereoFrame = plotFrame({
    title: 'Stereogram',
    toolbar: [
      el('button.button', {
        type: 'button',
        text: 'Envelope',
        'aria-pressed': 'true',
        title: 'Outline the poles the holder can bring onto the beam',
        onclick: (event) => {
          state.showEnvelope = !state.showEnvelope;
          event.currentTarget.setAttribute('aria-pressed', String(state.showEnvelope));
          drawStereogram();
        },
      }),
      el('button.button', {
        type: 'button',
        text: 'Labels',
        'aria-pressed': 'true',
        title: 'Name the low-index poles',
        onclick: (event) => {
          state.showPoleLabels = !state.showPoleLabels;
          event.currentTarget.setAttribute('aria-pressed', String(state.showPoleLabels));
          drawStereogram();
        },
      }),
    ],
  });

  const details = el('div');
  context.stage.append(
    el('div.tem-stage', {}, [frame.element, stereoFrame.element]),
    details,
  );

  /* ------------------------------------------------------------ controls */

  const fileInput = el('input', {
    type: 'file',
    accept: 'image/*',
    onchange: (event) => loadImage(event.target.files?.[0]),
  });

  const galleryHost = el('div');
  const centreHost = el('div.centre-tool');
  // Inside the pick tool, and stable across a fit: see `syncFitViews`.
  const centreButtonHost = el('div.button-row');
  const hintHost = el('div');
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

  /*
   * The four steps are one accordion, not four open panels.
   *
   * Open together they are 3400 px of rail against 760 px of window, so every
   * step but the first was reached by scrolling and the step being worked on was
   * rarely the one on screen. `name` makes the browser close the others when one
   * opens — the steps are sequential, so exactly one is ever the current one —
   * and `openStep` moves the accordion on as the workflow advances, which is
   * what the scroll calls below were already trying to do by hand.
   */
  const STEP_GROUP = 'tem-workflow';
  const stepOne = el('details.group.step', { name: STEP_GROUP, open: true }, [
    el('summary', { text: '1 · Open a pattern' }),
  ]);
  const stepTwo = el('details.group.step', { name: STEP_GROUP }, [
    el('summary', { text: '2 · Calibrate and index' }),
  ]);
  const stepThree = el('details.group.step', { name: STEP_GROUP }, [
    el('summary', { text: '3 · Where to go next' }),
  ]);
  const stepFour = el('details.group.step', { name: STEP_GROUP }, [
    el('summary', { text: '4 · Tilt to the target' }),
  ]);

  /** Open a step and put it on screen; a no-op if it is already the open one. */
  function openStep(step) {
    step.open = true;
    step.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }

  context.rail.append(
    append(stepOne, [
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
    append(stepTwo, [
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
    append(stepThree, [
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
    append(stepFour, [el('div.group__body', {}, [tiltHost, tiltButton])]),
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

    // The stereogram is a view of this form, so it follows every edit to it
    // rather than waiting for the plan to be pressed.
    state.tiltForm = buildForm(tiltOperation, {
      initial: tilt,
      onChange: () => {
        scheduleStereogram();
        // The target zone axis chooses the connecting band, so the overlay
        // follows the same control rather than adding one of its own.
        scheduleKikuchi();
      },
    });
    tiltHost.replaceChildren(state.tiltForm.element);
    scheduleStereogram();
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
      state.patternKey = `gallery:${entryId}:${nextPatternSerial()}`;
      state.mode = 'pick';
      state.picks = { centre: null, spots: [] };
      state.pickedCentre = null;
      state.fit = null;
      state.solutions = [];
      state.accepted = null;
      state.selectedPick = 'centre';
      state.highlight = null;
      calculatedButton.hidden = true;
      resetKikuchi();
      state.image = null;
      state.pixels = null;
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

      renderPickTool();
      drawPattern();
      renderResult(details, result, {
        teaches: state.teaches,
        extra: [galleryGuidance(result)],
      });
      state.teaches = null;
      // The pattern is open; picking is what happens next.
      openStep(stepTwo);
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
        state.patternKey = `image:${file.name}:${nextPatternSerial()}`;
        state.pixels = readLuminance(image);
        state.source = { kind: 'image' };
        state.gallery = null;
        state.mode = 'pick';
        state.picks = { centre: null, spots: [] };
        state.pickedCentre = null;
        state.fit = null;
        state.solutions = [];
        state.accepted = null;
        state.selectedPick = 'centre';
        state.highlight = null;
        calculatedButton.hidden = true;
        resetKikuchi();
        answerButton.hidden = true;
        autoPickButton.hidden = true;
        renderPickTool();
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

  /**
   * The micrograph's brightness, one byte per pixel, for centroiding a pick.
   *
   * Read once on load rather than per click: a click is an interaction and a
   * megapixel `getImageData` is not free. It stays in this tab — the panel's
   * promise is that only coordinates are sent, and this is what makes the
   * coordinates worth sending.
   *
   * Returns null if the canvas refuses to give up its pixels, which happens for
   * a cross-origin image; picking then works exactly as before, without the
   * snap.
   */
  function readLuminance(image) {
    try {
      const canvas = document.createElement('canvas');
      canvas.width = image.width;
      canvas.height = image.height;
      const context2d = canvas.getContext('2d', { willReadFrequently: true });
      if (!context2d) return null;
      context2d.drawImage(image, 0, 0);
      const { data } = context2d.getImageData(0, 0, image.width, image.height);
      const luminance = new Uint8ClampedArray(image.width * image.height);
      for (let index = 0; index < luminance.length; index += 1) {
        const at = index * 4;
        // Rec. 601 luma. A diffraction plate is grey, so any sensible weighting
        // agrees; this one is the conventional choice.
        luminance[index] =
          0.299 * data[at] + 0.587 * data[at + 1] + 0.114 * data[at + 2];
      }
      return { width: image.width, height: image.height, data: luminance };
    } catch {
      return null;
    }
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
    state.selectedPick = 'centre';
    state.highlight = null;
    renderPickTool();
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
      syncFitViews();
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
    // Update, never rebuild: the answer arrives while the user is still typing.
    syncFitViews();
    drawPattern();
  }

  /* --------------------------------------------------- the picks as numbers */

  /** Every pick as one addressable thing: `'centre'`, or a spot's index. */
  function pickAt(which) {
    if (which === 'centre') {
      return state.picks.centre
        ? { x: state.picks.centre[0], y: state.picks.centre[1] }
        : null;
    }
    const spot = state.picks.spots[which];
    return spot ? { x: spot.x, y: spot.y } : null;
  }

  /**
   * Put a pick at an exact coordinate, from a typed value or a nudge.
   *
   * One route for every way a coordinate can change, so that a typed number, an
   * arrow key and a click all leave the same state behind and all trigger the
   * same re-fit. `structural` rebuilds the table; a coordinate change does not,
   * because the table holds the input the user is typing into.
   */
  function movePick(which, x, y, { structural = false } = {}) {
    if (!Number.isFinite(x) || !Number.isFinite(y)) return;
    if (which === 'centre') {
      if (!state.picks.centre) return;
      // Remember where the beam was clicked the first time it is moved by hand,
      // so "Undo refinement" has something to go back to.
      if (!state.pickedCentre) state.pickedCentre = [...state.picks.centre];
      state.picks.centre = [x, y];
    } else {
      const spot = state.picks.spots[which];
      if (!spot) return;
      spot.x = x;
      spot.y = y;
    }
    if (structural) renderPickTool();
    else refreshPickReadouts();
    scheduleFit();
    drawPattern();
  }

  /** Move the selected pick by the current step, in picking units. */
  function nudgeSelected(dx, dy) {
    const pick = pickAt(state.selectedPick);
    if (!pick) return;
    movePick(state.selectedPick, pick.x + dx * state.nudgeStep, pick.y + dy * state.nudgeStep);
    const row = pickRowFor(state.selectedPick);
    if (row) {
      row.x.value = formatNumber(pickAt(state.selectedPick).x, 2);
      row.y.value = formatNumber(pickAt(state.selectedPick).y, 2);
    }
  }

  function pickRowFor(which) {
    return state.pickRows.find((row) => row.which === which) ?? null;
  }

  /**
   * Make a pick the one the pad and the arrow keys act on.
   *
   * In place, never by rebuilding the table. Selection happens on focus, and a
   * rebuild on focus destroys the input the user has just clicked into — the
   * field is replaced, the caret goes with it, and typing a coordinate becomes
   * impossible in a way that looks like the keyboard is broken.
   */
  function selectPick(which) {
    state.selectedPick = which;
    for (const row of state.pickRows) {
      row.element.classList.toggle('picks__row--selected', state.selectedPick === row.which);
      row.element.classList.toggle('picks__row--marked', isMarked(row.which));
    }
    if (state.nudgeLabel) state.nudgeLabel.textContent = nudgeLabelText();
    drawPattern();
  }

  function highlightPick(which) {
    state.highlight = which;
    for (const row of state.pickRows) {
      row.element.classList.toggle('picks__row--marked', isMarked(row.which));
    }
    drawPattern();
  }

  function removeSpot(index) {
    state.picks.spots.splice(index, 1);
    if (state.selectedPick === index) state.selectedPick = 'centre';
    else if (typeof state.selectedPick === 'number' && state.selectedPick > index) state.selectedPick -= 1;
    state.highlight = null;
    renderPickTool();
    scheduleFit();
    drawPattern();
  }

  /**
   * Make a spot the transmitted beam, and the old beam a spot.
   *
   * A swap rather than a move, because the commonest way to get this wrong is to
   * click a strong reflection first and realise on seeing the overlay; throwing
   * the old beam away would then lose a perfectly good pick, and it is usually a
   * reflection.
   */
  function promoteToBeam(index) {
    const spot = state.picks.spots[index];
    if (!spot) return;
    const previous = state.picks.centre;
    state.picks.spots.splice(index, 1);
    if (previous) state.picks.spots.splice(index, 0, { x: previous[0], y: previous[1] });
    state.picks.centre = [spot.x, spot.y];
    state.pickedCentre = null;
    state.selectedPick = 'centre';
    renderPickTool();
    scheduleFit();
    drawPattern();
  }

  function adoptRefinedCentre() {
    // Only a centre that was actually solved for is worth adopting; below four
    // spots the "refined" centre is the clicked one.
    if (!state.fit?.data.centre_refined) return;
    if (!state.pickedCentre) state.pickedCentre = [...state.picks.centre];
    state.picks.centre = [...state.fit.data.centre];
    renderPickTool();
    scheduleFit();
    drawPattern();
  }

  function restorePickedCentre() {
    if (!state.pickedCentre) return;
    state.picks.centre = [...state.pickedCentre];
    state.pickedCentre = null;
    renderPickTool();
    scheduleFit();
    drawPattern();
  }

  /**
   * The picks as numbers: every coordinate readable, typable and nudgeable.
   *
   * Clicking is how a pattern is picked and it is not enough on its own. A click
   * carries the precision of a mouse over a spot a few pixels across, there is
   * no way to say "put it exactly there" or to check what was actually stored,
   * and a coordinate that can only be produced by clicking cannot be copied out
   * of a previous session or a published figure. So every pick appears here as
   * two editable numbers, with what they measure — radius, d-spacing, and the
   * angle from the first spot — beside them.
   *
   * **Rebuilt only when the set of rows changes.** Everything that happens while
   * a coordinate is being edited — the fit returning two hundred milliseconds
   * later, the derived columns changing, a hint appearing — updates the existing
   * nodes instead. Replacing the table under an editor destroys the input it is
   * typing into, which looks from the outside like a keyboard that has stopped
   * working, and that is exactly what an earlier draft of this panel did on
   * every fit and on every focus.
   */
  function renderPickTool() {
    state.pickRows = [];
    if (!state.picks.centre) {
      centreHost.replaceChildren(
        el('p.field__help', {
          text:
            'Click the transmitted beam on the pattern to start. Every pick then appears here ' +
            'as a coordinate you can read, type over, or nudge a fraction of a pixel at a time.',
        }),
      );
      return;
    }
    centreHost.replaceChildren(
      el('div.field__label', { text: 'Picked coordinates' }),
      el('p.field__help', {
        text:
          'The beam is not a reflection: it is the origin every spot is measured from, so an ' +
          'error there biases every d-spacing at once while leaving the pattern ' +
          'self-consistent. Select a row to nudge it with the pad or the arrow keys; the ' +
          'overlay re-fits as you go.',
      }),
      pickTable(),
      nudgePad(),
      centreButtonHost,
      hintHost,
      coordinateText(),
    );
    syncFitViews();
  }

  /**
   * Bring everything that depends on the fit up to date, without rebuilding.
   *
   * The fit arrives while the user is still working — that is the point of it
   * being live — so it must not be able to take the panel apart while they are
   * mid-edit.
   */
  function syncFitViews() {
    if (!state.picks.centre) return;
    const fit = state.fit?.data;
    const refined = Boolean(fit?.centre_refined);
    centreButtonHost.replaceChildren(...centreButtons(refined));
    hintHost.replaceChildren(...fitHints(fit));
    refreshPickReadouts();
  }

  /** One row per pick: label, x, y, and what the pair measures. */
  function pickTable() {
    const rows = [
      pickRow('centre', 'Beam'),
      ...state.picks.spots.map((_, index) => pickRow(index, `Spot ${index + 1}`)),
    ];
    return el('div.picks', {}, [
      el('div.picks__head', {}, [
        el('span', { text: 'Pick' }),
        el('span', { text: 'x / px' }),
        el('span', { text: 'y / px' }),
        el('span', { text: '' }),
      ]),
      ...rows,
    ]);
  }

  function pickRow(which, label) {
    const pick = pickAt(which);
    const coordinate = (axis) =>
      el('input.picks__input', {
        type: 'number',
        step: '0.1',
        value: formatNumber(pick[axis], 2),
        'aria-label': `${label} ${axis}`,
        onfocus: () => selectPick(which),
        // `change` rather than `input`: a re-fit per keystroke would fire a
        // request for "5", "51" and "512" on the way to a three-digit
        // coordinate, and the first two are somewhere else entirely.
        onchange: (event) => commitCoordinate(which, axis, event.target),
      });
    const x = coordinate('x');
    const y = coordinate('y');
    const readout = el('span.picks__measure', { text: measureText(which) });
    const element = el(
      'div.picks__row',
      {
        onpointerenter: () => highlightPick(which),
        onpointerleave: () => highlightPick(null),
        onkeydown: (event) => onPickKey(event, which),
      },
      [
        el('button.picks__label', {
          type: 'button',
          text: label,
          title: 'Select this pick for the nudge pad and the arrow keys',
          onclick: () => selectPick(which),
        }),
        x,
        y,
        which === 'centre'
          ? el('span')
          : el('span.picks__actions', {}, [
              el('button.button.button--small', {
                type: 'button',
                text: '⌖',
                title: 'Make this spot the transmitted beam, and the beam a spot',
                'aria-label': `Make spot ${which + 1} the transmitted beam`,
                onclick: () => promoteToBeam(which),
              }),
              el('button.button.button--small', {
                type: 'button',
                text: '×',
                title: 'Remove this pick',
                'aria-label': `Remove spot ${which + 1}`,
                onclick: () => removeSpot(which),
              }),
            ]),
        // What the pair measures, on its own line under the numbers. Beside
        // them it took the width the numbers needed, and a coordinate field too
        // narrow to show its own value is worse than no field at all.
        readout,
      ],
    );
    element.classList.toggle('picks__row--selected', state.selectedPick === which);
    element.classList.toggle('picks__row--marked', isMarked(which));
    if (fitVerdict(which) === 'off') element.classList.add('picks__row--outlier');
    state.pickRows.push({ which, element, x, y, readout });
    return element;
  }

  /**
   * Accept a typed coordinate, or refuse it and put the old one back.
   *
   * Refusing in place matters more than it sounds: an empty or malformed field
   * that silently became `NaN` would move the pick to nowhere, take the overlay
   * with it, and leave the user looking at a blank field wondering what they
   * broke.
   */
  function commitCoordinate(which, axis, input) {
    const pick = pickAt(which);
    if (!pick) return;
    const value = Number(input.value);
    if (!Number.isFinite(value)) {
      input.value = formatNumber(pick[axis], 2);
      return;
    }
    const next = { ...pick, [axis]: value };
    movePick(which, next.x, next.y);
    input.value = formatNumber(value, 2);
  }

  /** Arrow keys nudge the pick whose row has focus, by the current step. */
  function onPickKey(event, which) {
    const steps = {
      ArrowLeft: [-1, 0],
      ArrowRight: [1, 0],
      ArrowUp: [0, -1],
      ArrowDown: [0, 1],
    };
    const step = steps[event.key];
    // A number input spends its own arrow keys on incrementing itself, which is
    // the same gesture by a different route; leave it to do that.
    if (!step || event.target.tagName === 'INPUT') return;
    event.preventDefault();
    state.selectedPick = which;
    nudgeSelected(step[0], step[1]);
  }

  /** Update the derived columns without touching the inputs. */
  function refreshPickReadouts() {
    for (const row of state.pickRows) {
      row.readout.textContent = measureText(row.which);
      row.element.classList.toggle('picks__row--outlier', fitVerdict(row.which) === 'off');
    }
  }

  /** What a pick measures: radius from the beam, d, and angle from spot 1. */
  function measureText(which) {
    if (which === 'centre') {
      const fit = state.fit?.data;
      if (!fit) return 'origin';
      if (!fit.centre_refined) return 'origin · held';
      return `origin · fit says ${formatNumber(fit.centre_shift, 1)} px away`;
    }
    const pick = pickAt(which);
    const centre = state.picks.centre;
    if (!pick || !centre) return '';
    const dx = pick.x - centre[0];
    const dy = pick.y - centre[1];
    const radius = Math.hypot(dx, dy);
    const g = reciprocalRadius(dx, dy);
    const first = state.picks.spots[0];
    let angle = '';
    if (which !== 0 && first) {
      const fx = first.x - centre[0];
      const fy = first.y - centre[1];
      const degrees = angleAtBeamDeg(
        { dx, dy, r: radius },
        { dx: fx, dy: fy, r: Math.hypot(fx, fy) },
      );
      if (degrees !== null) angle = ` · ${formatNumber(degrees, 2)}°`;
    }
    return (
      `${formatNumber(radius, 1)} px · ` +
      `${g > 0 ? formatNumber(1 / g, 4) : '∞'} Å${angle}`
    );
  }

  /**
   * The angle at the beam between two offsets, in degrees.
   *
   * Image-plane arithmetic on two clicked points and nothing more — the same
   * class of measurement as the cursor's distance readout. It is defined once
   * and used by both the overlay card and the coordinate table, so the number
   * beside a row and the number in the card cannot disagree.
   *
   * @param {{dx: number, dy: number, r: number}} row
   * @param {{dx: number, dy: number, r: number}} reference
   * @returns {number|null} Null when either offset has no length to take an
   *   angle from.
   */
  function angleAtBeamDeg(row, reference) {
    const norms = row.r * reference.r;
    if (!norms) return null;
    const cosine = (row.dx * reference.dx + row.dy * reference.dy) / norms;
    return (Math.acos(Math.min(Math.max(cosine, -1), 1)) * 180) / Math.PI;
  }

  /** Whether the fit explains this pick: 'on', 'off', or null if there is no fit. */
  function fitVerdict(which) {
    if (which === 'centre' || !state.fit) return null;
    return state.fit.data.outliers.includes(which + 1) ? 'off' : 'on';
  }

  function nudgeLabelText() {
    return state.selectedPick === 'centre' ? 'Nudging: beam' : `Nudging: spot ${state.selectedPick + 1}`;
  }

  function nudgePad() {
    // The buttons name the selection rather than a fixed pick, because which
    // pick they move changes without the pad being rebuilt.
    const arrow = (glyph, dx, dy, description) =>
      el('button.button.button--small', {
        type: 'button',
        text: glyph,
        title: description,
        'aria-label': description,
        onclick: () => nudgeSelected(dx, dy),
      });
    state.nudgeLabel = el('span.centre-tool__readout', { text: nudgeLabelText() });
    return el('div.centre-tool__pad', {}, [
      state.nudgeLabel,
      arrow('◀', -1, 0, 'Move the selected pick left'),
      el('div.centre-tool__column', {}, [
        arrow('▲', 0, -1, 'Move the selected pick up'),
        arrow('▼', 0, 1, 'Move the selected pick down'),
      ]),
      arrow('▶', 1, 0, 'Move the selected pick right'),
      el(
        'select.centre-tool__step',
        {
          'aria-label': 'Nudge step',
          onchange: (event) => {
            state.nudgeStep = Number(event.target.value);
          },
        },
        [0.1, 0.5, 1, 2, 5, 10].map((step) =>
          el('option', {
            value: String(step),
            text: `${step} px`,
            selected: step === state.nudgeStep,
          }),
        ),
      ),
    ]);
  }

  function centreButtons(refined) {
    // Offering to adopt a centre the fit did not solve for would hand the user a
    // copy of their own click dressed up as a measurement.
    return [
      el('button.button.button--small', {
        type: 'button',
        text: 'Refine beam from the spots',
        disabled: !refined,
        title: refined
          ? 'Move the beam to where the fitted lattice puts it'
          : 'Needs four or more spots: with fewer, the centre cannot be solved for',
        onclick: () => adoptRefinedCentre(),
      }),
      el('button.button.button--small', {
        type: 'button',
        text: 'Undo beam move',
        disabled: !state.pickedCentre,
        title: 'Put the beam back where it was clicked',
        onclick: () => restorePickedCentre(),
      }),
    ];
  }

  function fitHints(fit) {
    if (!fit) {
      return [
        el('p.field__hint', {
          text: 'Pick two spots to lay a lattice over the pattern, and four to solve for the beam.',
        }),
      ];
    }
    return [
      fit.basis_vectors?.length
        ? el('p.field__hint', {
            text:
              'The gold arrows are the two lattice vectors ' +
              fit.basis_vectors
                .map((vector) =>
                  vector.on_a_pick
                    ? `${vector.label} (spot ${vector.spot}, ` +
                      `${formatNumber(vector.length, 1)} px)`
                    : `${vector.label} (drawn dashed: no pick sits on this node)`,
                )
                .join(' and ') +
              '. Every line in the teal grid turns with them, so those are the two picks worth ' +
              'adjusting first.',
          })
        : null,
      fit.fit.notes.length ? el('p.field__hint', { text: fit.fit.notes.join(' ') }) : null,
      fit.outliers.length
        ? el('p.field__hint', {
            text: `Spot(s) ${fit.outliers.join(', ')} do not sit on the fitted lattice.`,
          })
        : null,
    ].filter(Boolean);
  }

  /**
   * The whole pick set as text, to carry out of the session and back into it.
   *
   * A set of picks is a measurement, and a measurement that exists only as
   * clicks in a browser tab cannot be checked, repeated, or handed to anyone.
   * One line per pick, beam first — the same order the table shows and the
   * service receives.
   */
  function coordinateText() {
    // A textarea's value is its content, not an attribute.
    const area = el(
      'textarea.picks__text',
      { rows: 4, spellcheck: 'false', 'aria-label': 'Every pick as text, beam first' },
      [picksAsText()],
    );
    return el('details.picks__io', {}, [
      el('summary', { text: 'Coordinates as text' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'One "x, y" per line, transmitted beam first. Copy it to keep this measurement, ' +
            'or paste a set in and apply it.',
        }),
        area,
        el('div.button-row', {}, [
          el('button.button.button--small', {
            type: 'button',
            text: 'Apply these coordinates',
            onclick: () => applyPicksFromText(area),
          }),
          el('button.button.button--small', {
            type: 'button',
            text: 'Reset to the picks',
            onclick: () => {
              area.value = picksAsText();
            },
          }),
        ]),
      ]),
    ]);
  }

  function picksAsText() {
    if (!state.picks.centre) return '';
    return [state.picks.centre, ...state.picks.spots.map((spot) => [spot.x, spot.y])]
      .map(([x, y]) => `${formatNumber(x, 2)}, ${formatNumber(y, 2)}`)
      .join('\n');
  }

  function applyPicksFromText(area) {
    const points = [];
    for (const line of area.value.split('\n')) {
      const trimmed = line.trim();
      if (!trimmed) continue;
      const parts = trimmed.split(/[,;\s]+/).map(Number);
      if (parts.length < 2 || !parts.slice(0, 2).every(Number.isFinite)) {
        frame.setStatus(`Could not read "${trimmed}" as a coordinate pair. Nothing was changed.`);
        return;
      }
      points.push([parts[0], parts[1]]);
    }
    if (!points.length) {
      frame.setStatus('No coordinates to apply.');
      return;
    }
    state.picks = {
      centre: [points[0][0], points[0][1]],
      spots: points.slice(1).map(([x, y]) => ({ x, y })),
    };
    state.pickedCentre = null;
    state.selectedPick = 'centre';
    state.highlight = null;
    renderPickTool();
    scheduleFit();
    drawPattern();
    frame.setStatus(
      `Beam and ${state.picks.spots.length} spot(s) set from the coordinates you entered.`,
    );
  }

  function calibrationValues() {
    const values = state.solveForm ? state.solveForm.values() : {};
    return {
      units: values.units ?? 'px',
      camera: Number(values.camera_constant_mm_angstrom ?? 180),
      pixel: Number(values.pixel_size_mm ?? 0.05),
      scale: Number(values.reciprocal_per_px_angstrom ?? 0),
    };
  }

  /**
   * Convert a picked coordinate offset into the reciprocal-space radius.
   *
   * The same four cases the service applies, because the cursor readout and the
   * indexed result must be the same measurement: a readout that disagreed with
   * the table below it would be worse than no readout.
   */
  function reciprocalRadius(dx, dy) {
    const { units, camera, pixel, scale } = calibrationValues();
    const distance = Math.hypot(dx, dy);
    if (units === 'reciprocal_angstrom') return distance;
    if (units === 'px_scale') return scale > 0 ? distance * scale : 0;
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
    const outer = svg('svg', {
      viewBox: `0 0 ${width} ${height}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': 'Diffraction pattern',
    });

    /*
     * Everything is drawn inside a group clipped to the image rectangle.
     *
     * A viewBox is a coordinate system, not a boundary: with `meet` the drawing
     * is letterboxed, and anything painted outside the viewBox is painted over
     * the letterbox. The fitted lattice is generated by walking node indices
     * outwards until it certainly covers the diagonal, so it always overshoots —
     * and those overshoots ran across the blank margins beside the pattern,
     * claiming lattice where there is no image. The clip makes the picture's own
     * edge the edge of every overlay, which is what a person reading the pattern
     * assumes it already is.
     */
    const clipId = 'tem-pattern-clip';
    const root = svg('g', { 'clip-path': `url(#${clipId})` });
    outer.append(
      svg('defs', {}, [
        svg('clipPath', { id: clipId }, [svg('rect', { x: 0, y: 0, width, height })]),
      ]),
      root,
    );

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

    // Under the spots: the bands are the background the spots are read against.
    if (state.showKikuchi) drawKikuchiBands(root, width, height);
    if (state.showLattice && state.fit) drawFittedLattice(root, width, height);
    if (state.showCalculated) drawCalculatedPattern(root, width, height);
    if (state.calibrate.active) drawCalibrationLine(root, width, height);

    const marker = Math.max(width, height) / 140;
    if (state.picks.centre) {
      const [cx, cy] = state.picks.centre;
      const beamStroke = isMarked('centre') ? SELECTED_COLOUR : BEAM_COLOUR;
      root.append(
        svg('circle', {
          cx, cy, r: marker * 1.4,
          fill: 'none',
          stroke: beamStroke,
          'stroke-width': marker / 3,
        }),
        svg('line', {
          x1: cx - marker * 2.4, y1: cy, x2: cx + marker * 2.4, y2: cy,
          stroke: beamStroke, 'stroke-width': marker / 4,
        }),
        svg('line', {
          x1: cx, y1: cy - marker * 2.4, x2: cx, y2: cy + marker * 2.4,
          stroke: beamStroke, 'stroke-width': marker / 4,
        }),
      );
      drawRefinedCentre(root, cx, cy, marker);
    }

    state.picks.spots.forEach((spot, index) => {
      const dx = state.picks.centre ? spot.x - state.picks.centre[0] : 0;
      const dy = state.picks.centre ? spot.y - state.picks.centre[1] : 0;
      const g = reciprocalRadius(dx, dy);
      const marked = isMarked(index);
      const colour = marked ? SELECTED_COLOUR : PICK_COLOUR;
      const node = svg('circle', {
        cx: spot.x, cy: spot.y, r: marker,
        fill: 'none',
        stroke: colour,
        'stroke-width': marker / 3,
      });
      root.append(node);
      // A second, wider ring on the pick the coordinate table is working on, so
      // "row 2" and "that spot" are the same thing without having to count.
      if (marked) {
        root.append(
          svg('circle', {
            cx: spot.x, cy: spot.y, r: marker * 2.2,
            fill: 'none',
            stroke: SELECTED_COLOUR,
            'stroke-width': marker / 4,
            'stroke-opacity': 0.8,
            'stroke-dasharray': `${marker / 2} ${marker / 2}`,
          }),
        );
      }
      root.append(
        svg('text', {
          x: spot.x + marker * 1.8, y: spot.y - marker * 0.6,
          'font-size': marker * 2.4,
          fill: colour,
          'paint-order': 'stroke',
          stroke: HALO_COLOUR,
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

    outer.addEventListener('click', (event) => {
      // The frame owns the camera, so the frame converts the pointer. A private
      // copy of this arithmetic is right only while the view is unzoomed and
      // unpanned, which is exactly when nobody is picking carefully.
      const point = frame.pointerToData(event);
      if (!point || !Number.isFinite(point.x) || !Number.isFinite(point.y)) return;
      // The letterbox beside a non-square image is inside the element but
      // outside the picture. A pick there has no pixel under it, so it can only
      // be a mis-click, and accepting one silently put a spot in a place the
      // user could see was empty.
      if (point.x < 0 || point.y < 0 || point.x > width || point.y > height) {
        frame.setStatus('That is outside the image — pick inside the pattern.');
        return;
      }
      // While calibrating, a click is a ruler end rather than a pick. The tool
      // is modal on purpose: the two gestures are the same gesture, and a
      // modifier key nobody is told about is not a discoverable alternative.
      if (state.calibrate.active) {
        if (state.calibrate.points.length >= 2) state.calibrate.points = [];
        state.calibrate.points.push({ x: point.x, y: point.y });
        drawPattern();
        return;
      }
      const snapped = snapToSpot(point);
      let picked;
      if (!state.picks.centre) {
        state.picks.centre = [snapped.x, snapped.y];
        state.pickedCentre = null;
        state.selectedPick = 'centre';
        picked = 'Direct beam';
      } else {
        state.picks.spots.push({ x: snapped.x, y: snapped.y });
        state.selectedPick = state.picks.spots.length - 1;
        picked = `Spot ${state.picks.spots.length}`;
      }
      // The console keeps every coordinate the session produced. A pick is
      // otherwise recoverable only by reading it back off the rail before the
      // next one replaces the reader's attention.
      log.info(
        `${picked} is selected: coordinates are (${formatNumber(snapped.x, 2)}, ` +
          `${formatNumber(snapped.y, 2)}) px.`,
        {
          source: 'tem',
          detail: {
            x: Number(snapped.x.toFixed(2)),
            y: Number(snapped.y.toFixed(2)),
            snapped_px: Number(snapped.moved.toFixed(2)),
          },
        },
      );
      if (snapped.moved > 0.05) {
        frame.setStatus(
          `Picked at ${formatNumber(point.x, 1)}, ${formatNumber(point.y, 1)} and moved ` +
            `${formatNumber(snapped.moved, 1)} px to the spot centroid. Type an exact ` +
            'coordinate in the rail if you want it somewhere else.',
        );
      }
      renderPickTool();
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
    /*
     * The camera survives a redraw of the *same* pattern, and only that.
     *
     * Every pick, nudge and typed coordinate rebuilds this SVG, and without the
     * preserved view it snapped back to Fit each time — so zooming in to place
     * a spot precisely was impossible, which is the one situation where zooming
     * in is the point.
     *
     * But a *newly opened* pattern is a different picture, and inheriting the
     * previous one's zoom and pan put it off camera: open a practice plate,
     * zoom to 400%, then open a micrograph of your own, and the frame holds a
     * ninety-pixel crop of an image it has never shown whole. Nothing errors,
     * so it reads as the panel simply refusing to display the file. Keying the
     * camera to the pattern rather than to the frame is what makes "opened" and
     * "redrawn" different events.
     */
    const samePattern = state.patternKey !== null && state.patternKey === state.drawnKey;
    state.drawnKey = state.patternKey;
    frame.setContent(outer, { preserveViewport: samePattern });
    frame.setOverlay(measurementCard());
    frame.setControls(state.calibrate.active ? calibrationCard() : null);
    if (state.calibrate.active) {
      frame.setStatus(
        state.calibrate.points.length === 2
          ? `Measured ${formatNumber(calibrationLengthPx(), 1)} px — say how long that is below.`
          : 'Calibrating: click the two ends of a length you know.',
      );
      return;
    }
    frame.setStatus(
      state.picks.centre
        ? `Beam marked · ${state.picks.spots.length} spot(s) picked · click to add more`
        : 'Click the transmitted beam first — it is the origin every spot is measured from',
    );
  }

  /** Whether this pick is the selected or hovered one, and so drawn ringed. */
  function isMarked(which) {
    return state.selectedPick === which || state.highlight === which;
  }

  /**
   * The centre the fit prefers, drawn beside the one that was clicked.
   *
   * Only when the two differ, and never *instead of* the clicked beam. The
   * overlay used to draw its grid from the refined centre while the crosshair
   * stayed on the click, with nothing on screen to say they were two different
   * points — so a refinement that had run away looked like a pick that had
   * landed in the wrong place. Two marks and a stated distance is the honest
   * form: it is a proposal until the user adopts it.
   */
  function drawRefinedCentre(root, cx, cy, marker) {
    const data = state.fit?.data;
    if (!data?.centre_refined) return;
    const [rx, ry] = data.centre;
    if (Math.hypot(rx - cx, ry - cy) < 0.25) return;
    const group = svg('g', { 'pointer-events': 'none' });
    group.append(
      svg('line', {
        x1: cx, y1: cy, x2: rx, y2: ry,
        stroke: REFINED_COLOUR,
        'stroke-width': marker / 4,
        'stroke-dasharray': `${marker} ${marker}`,
        'stroke-opacity': 0.9,
      }),
      svg('circle', {
        cx: rx, cy: ry, r: marker * 1.4,
        fill: 'none',
        stroke: REFINED_COLOUR,
        'stroke-width': marker / 3,
        'stroke-dasharray': `${marker / 2} ${marker / 2}`,
      }),
      svg('text', {
        x: rx + marker * 2, y: ry + marker * 3,
        'font-size': marker * 2,
        fill: REFINED_COLOUR,
        'paint-order': 'stroke',
        stroke: HALO_COLOUR,
        'stroke-width': marker / 2,
        text: 'fitted centre',
      }),
    );
    root.append(group);
  }

  /**
   * Move a click to the centre of the spot it landed on, if it landed on one.
   *
   * A click is worth a few pixels; a spot centroid is worth a fraction of one,
   * and every d-spacing in the pattern is measured from these two points. The
   * radius is generous enough to catch a click on the edge of a spot and tight
   * enough that a click on background stays where it was put — a pick on
   * background is a legitimate thing to do, and moving it silently to a
   * neighbouring reflection would be worse than not snapping at all.
   *
   * On a practice plate the true centres are known, so the nearest one is the
   * answer. On a micrograph the intensity-weighted centroid of a small window is
   * used, with the window's own floor subtracted so the background does not drag
   * the answer toward the middle of the box.
   *
   * @returns {{x: number, y: number, moved: number}}
   */
  function snapToSpot(point) {
    const size = frameSize();
    if (!size) return { x: point.x, y: point.y, moved: 0 };
    const radius = Math.max(8, Math.max(size.width, size.height) / 64);
    const centred = state.gallery
      ? nearestSimulatedSpot(point, radius)
      : iteratedCentroid(point, radius);
    if (!centred) return { x: point.x, y: point.y, moved: 0 };
    const moved = Math.hypot(centred.x - point.x, centred.y - point.y);
    // A centroid that has walked further than the window it started in has
    // found a different spot, not this one. Refusing is the safe answer: the
    // pick stays where it was put and the user can see that it did.
    if (moved > radius) return { x: point.x, y: point.y, moved: 0 };
    return { x: centred.x, y: centred.y, moved };
  }

  /**
   * Centre of mass, re-centred on its own answer until it stops moving.
   *
   * One pass is not enough and is worse than it looks. The window is centred on
   * the *click*, so on a peak wider than the window the intensity across it is
   * near-linear and the centre of mass lands close to the middle of the window —
   * that is, back at the click. A snap that returns the click looks like a snap
   * that is working. Re-centring the window on each estimate makes the true
   * peak the fixed point, and three or four passes reach it from anywhere inside
   * the starting window.
   */
  function iteratedCentroid(point, radius) {
    let current = { x: point.x, y: point.y };
    for (let pass = 0; pass < 6; pass += 1) {
      const next = centroidOfPixels(current, radius);
      if (!next) return pass === 0 ? null : current;
      const step = Math.hypot(next.x - current.x, next.y - current.y);
      current = next;
      if (step < 0.05) break;
    }
    return current;
  }

  function nearestSimulatedSpot(point, radius) {
    const pattern = state.gallery.data.pattern;
    let best = null;
    let bestDistance = radius;
    // The transmitted beam is a feature of the plate like any other, and it is
    // the brightest one: a centroid on a real micrograph would find it too.
    const beam = pattern.centre_px;
    const candidates = beam
      ? [{ x: beam[0], y: beam[1] }, ...(pattern.spots ?? [])]
      : (pattern.spots ?? []);
    for (const spot of candidates) {
      const distance = Math.hypot(spot.x - point.x, spot.y - point.y);
      if (distance < bestDistance) {
        bestDistance = distance;
        best = { x: spot.x, y: spot.y };
      }
    }
    return best;
  }

  /**
   * Intensity-weighted centroid of an uploaded image, in a window about a point.
   *
   * The window's minimum is subtracted before weighting. Without that, a uniform
   * background contributes a large, perfectly symmetric weight centred on the
   * *window*, which pulls every centroid toward the click and makes the snap
   * look like it is working when it is doing nothing.
   */
  function centroidOfPixels(point, radius) {
    const pixels = state.pixels;
    if (!pixels) return null;
    const left = Math.max(0, Math.floor(point.x - radius));
    const top = Math.max(0, Math.floor(point.y - radius));
    const right = Math.min(pixels.width - 1, Math.ceil(point.x + radius));
    const bottom = Math.min(pixels.height - 1, Math.ceil(point.y + radius));
    if (right <= left || bottom <= top) return null;
    let floor = Infinity;
    let peak = -Infinity;
    for (let y = top; y <= bottom; y += 1) {
      for (let x = left; x <= right; x += 1) {
        const value = pixels.data[y * pixels.width + x];
        if (value < floor) floor = value;
        if (value > peak) peak = value;
      }
    }
    // Nothing to centre on: a window with no contrast is background, and a pick
    // on background belongs where it was put.
    if (!(peak - floor > 8)) return null;
    let sum = 0;
    let sumX = 0;
    let sumY = 0;
    for (let y = top; y <= bottom; y += 1) {
      for (let x = left; x <= right; x += 1) {
        if (Math.hypot(x - point.x, y - point.y) > radius) continue;
        const weight = pixels.data[y * pixels.width + x] - floor;
        if (weight <= 0) continue;
        sum += weight;
        sumX += weight * x;
        sumY += weight * y;
      }
    }
    if (!(sum > 0)) return null;
    return { x: sumX / sum, y: sumY / sum };
  }

  /**
   * What the picks measure, on the pattern rather than under it.
   *
   * A zone axis is identified from three numbers per spot and nothing else:
   * the d-spacing, its ratio to the first spot's, and the angle between them.
   * Those are the numbers an experienced microscopist reads off a plate against
   * a table of the phase, and they are worth having in view *while* picking —
   * a ratio near 1.000 with a 90° angle says "two ⟨220⟩ at right angles" before
   * any indexing runs, and a ratio that lands on nothing sensible says a pick
   * is on the wrong spot while it is still cheap to undo.
   *
   * Measured, not fitted: these come from the clicked coordinates and the
   * calibration, with no solution involved, so they can be compared with the
   * indexed answer rather than being a restatement of it. The first picked spot
   * is the reference for both the ratio and the angle, which is the convention
   * of every ratio table in the literature.
   */
  function measurementCard() {
    const spots = state.picks.spots;
    if (!state.picks.centre || !spots.length) return null;
    const [cx, cy] = state.picks.centre;
    const rows = spots.map((spot, index) => {
      const dx = spot.x - cx;
      const dy = spot.y - cy;
      const g = reciprocalRadius(dx, dy);
      return { index: index + 1, dx, dy, r: Math.hypot(dx, dy), d: g > 0 ? 1 / g : null };
    });
    const reference = rows[0];

    /** The angle at the beam between this spot and the first, in degrees. */
    const angleToReference = (row) => angleAtBeamDeg(row, reference);

    const cell = (text, className = '') =>
      el(className ? `td.${className}` : 'td', { text });

    return el('div.measure', {}, [
      el('div.measure__title', { text: `Measured picks · ${rows.length} spot(s)` }),
      el('table.measure__table', {}, [
        el('thead', {}, [
          el('tr', {}, [
            el('th', { text: '#' }),
            el('th', { text: 'R / px' }),
            el('th', { text: 'd / Å' }),
            el('th', { text: `d${reference.index}/d`, title: 'Ratio to the first picked spot' }),
            el('th', { text: '∠ / °', title: 'Angle at the beam from the first picked spot' }),
          ]),
        ]),
        el(
          'tbody',
          {},
          rows.map((row) => {
            const angle = angleToReference(row);
            const ratio = row.d && reference.d ? reference.d / row.d : null;
            return el('tr', {}, [
              cell(String(row.index)),
              cell(formatNumber(row.r, 1), 'measure__num'),
              cell(row.d ? formatNumber(row.d, 4) : '—', 'measure__num'),
              cell(ratio ? formatNumber(ratio, 3) : '—', 'measure__num'),
              cell(row.index === reference.index ? '—' : formatNumber(angle, 2), 'measure__num'),
            ]);
          }),
        ),
      ]),
      el('p.measure__note', {
        text: 'From the picked coordinates and the calibration — no solution involved.',
      }),
    ]);
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
    // The origin of the fitted lattice is drawn by `drawRefinedCentre`, with the
    // clicked beam beside it, because the only useful thing to say about it is
    // whether it agrees with the pick.
    root.append(group);
    drawBasisVectors(root, width, height);
  }

  /* --------------------------------------------------------- calibration */

  /** Units a measured line can be given in, and what each one calibrates. */
  const CALIBRATION_UNITS = [
    { id: 'inv_angstrom', label: 'Å⁻¹', reciprocal: true, toAngstrom: (value) => value },
    // 1 nm⁻¹ = 0.1 Å⁻¹, and a published scale bar is as often in nm⁻¹.
    { id: 'inv_nm', label: 'nm⁻¹', reciprocal: true, toAngstrom: (value) => value * 0.1 },
    { id: 'cm', label: 'cm on the plate', reciprocal: false, toMillimetre: (v) => v * 10 },
    { id: 'mm', label: 'mm on the plate', reciprocal: false, toMillimetre: (v) => v },
  ];

  function setCalibrating(active) {
    state.calibrate.active = Boolean(active);
    calibrateButton.setAttribute('aria-pressed', String(state.calibrate.active));
    if (!state.calibrate.active) state.calibrate.points = [];
    drawPattern();
  }

  function calibrationLengthPx() {
    const [first, second] = state.calibrate.points;
    if (!first || !second) return 0;
    return Math.hypot(second.x - first.x, second.y - first.y);
  }

  /**
   * Turn the measured line into a calibration and write it into the index form.
   *
   * Whichever unit was chosen, the result is one of the two numbers the form
   * already has, so nothing downstream learns a new concept: a reciprocal
   * length sets the direct scale and switches the coordinate units to it; a
   * length on the plate sets the pixel size and leaves the camera constant
   * doing its usual job.
   */
  function applyCalibration() {
    const distance = calibrationLengthPx();
    const unit = CALIBRATION_UNITS.find((entry) => entry.id === state.calibrate.unit);
    const length = Number(state.calibrate.length);
    if (!(distance > 0) || !Number.isFinite(length) || length <= 0 || !unit) {
      frame.setStatus('Draw a line across a known length, then say how long it is.');
      return;
    }
    let announcement;
    if (unit.reciprocal) {
      const scale = unit.toAngstrom(length) / distance;
      state.solveForm.setValues({ units: 'px_scale', reciprocal_per_px_angstrom: scale });
      announcement =
        `Calibrated: ${formatNumber(distance, 1)} px across ${length} ${unit.label} means ` +
        `1 px = ${formatNumber(scale, 5)} Å⁻¹. Every spacing below now uses it.`;
      log.notice(announcement, {
        source: 'tem',
        detail: { distance_px: Number(distance.toFixed(2)), reciprocal_per_px_angstrom: scale },
      });
    } else {
      const pixelSize = unit.toMillimetre(length) / distance;
      state.solveForm.setValues({ units: 'px', pixel_size_mm: pixelSize });
      announcement =
        `Calibrated: ${formatNumber(distance, 1)} px across ${length} ${unit.label} means a ` +
        `pixel size of ${formatNumber(pixelSize, 5)} mm, used with the camera constant.`;
      log.notice(announcement, {
        source: 'tem',
        detail: { distance_px: Number(distance.toFixed(2)), pixel_size_mm: pixelSize },
      });
    }
    // Leaving the tool redraws, and the redraw writes the ordinary picking
    // status — so the result of the calibration is announced after it, or it is
    // replaced by "click the transmitted beam" the instant it appears.
    setCalibrating(false);
    scheduleFit();
    refreshPickReadouts();
    frame.setStatus(announcement);
  }

  /** The card shown under the pattern while the calibration tool is on. */
  function calibrationCard() {
    const distance = calibrationLengthPx();
    const lengthInput = el('input.calibrate__length', {
      type: 'number',
      step: 'any',
      min: '0',
      value: state.calibrate.length,
      placeholder: 'known length',
      'aria-label': 'Known length of the drawn line',
      oninput: (event) => {
        state.calibrate.length = event.target.value;
      },
    });
    const unitSelect = el(
      'select.calibrate__unit',
      {
        'aria-label': 'Unit of the known length',
        onchange: (event) => {
          state.calibrate.unit = event.target.value;
        },
      },
      CALIBRATION_UNITS.map((entry) =>
        el('option', {
          value: entry.id,
          text: entry.label,
          selected: entry.id === state.calibrate.unit,
        }),
      ),
    );
    return el('div.calibrate', {}, [
      el('span.calibrate__hint', {
        text:
          state.calibrate.points.length < 2
            ? 'Click the two ends of something whose length you know — a scale bar, or a ' +
              'reflection whose spacing you know.'
            : `Measured ${formatNumber(distance, 1)} px. That length is:`,
      }),
      state.calibrate.points.length === 2 ? lengthInput : null,
      state.calibrate.points.length === 2 ? unitSelect : null,
      state.calibrate.points.length === 2
        ? el('button.button.button--primary', {
            type: 'button',
            text: 'Use this scale',
            onclick: () => applyCalibration(),
          })
        : null,
      state.calibrate.points.length
        ? el('button.button', {
            type: 'button',
            text: 'Redraw',
            onclick: () => {
              state.calibrate.points = [];
              drawPattern();
            },
          })
        : null,
    ]);
  }

  /** The measured line itself, drawn over the pattern with its length. */
  function drawCalibrationLine(root, width, height) {
    const marker = Math.max(width, height) / 140;
    const [first, second] = state.calibrate.points;
    const group = svg('g', { 'pointer-events': 'none' });
    for (const point of state.calibrate.points) {
      group.append(
        svg('circle', {
          cx: point.x, cy: point.y, r: marker,
          fill: 'none',
          stroke: REFINED_COLOUR,
          'stroke-width': marker / 3,
        }),
      );
    }
    if (first && second) {
      group.append(
        svg('line', {
          x1: first.x, y1: first.y, x2: second.x, y2: second.y,
          stroke: HALO_COLOUR,
          'stroke-width': marker,
          'stroke-opacity': 0.6,
        }),
        svg('line', {
          x1: first.x, y1: first.y, x2: second.x, y2: second.y,
          stroke: REFINED_COLOUR,
          'stroke-width': marker / 3,
        }),
        svg('text', {
          x: (first.x + second.x) / 2 + marker,
          y: (first.y + second.y) / 2 - marker,
          'font-size': marker * 2.4,
          fill: REFINED_COLOUR,
          'paint-order': 'stroke',
          stroke: HALO_COLOUR,
          'stroke-width': marker / 2,
          text: `${formatNumber(calibrationLengthPx(), 1)} px`,
        }),
      );
    }
    root.append(group);
  }

  /* ---------------------------------------------------------- stereogram */

  /**
   * Ask Python for the stereogram of whatever the tilt form currently says.
   *
   * The tilt step already holds every input this needs — the phase, the axis on
   * the beam, the target, the stage reading and the holder limits — so the
   * stereogram is a view of that form rather than a second set of controls to
   * keep in step with it. Debounced, because it redraws on every keystroke in
   * those fields.
   */
  /* ----------------------------------------------------- Kikuchi bands */

  /** Forget the bands: they belong to an orientation that no longer stands. */
  function resetKikuchi() {
    clearTimeout(state.kikuchiPending);
    state.kikuchi = null;
    state.kikuchiRequest = null;
    state.showKikuchi = false;
    kikuchiButton.hidden = true;
    kikuchiButton.setAttribute('aria-pressed', 'false');
  }

  function scheduleKikuchi() {
    if (!state.showKikuchi) return;
    clearTimeout(state.kikuchiPending);
    state.kikuchiPending = setTimeout(() => refreshKikuchi(), 180);
  }

  /**
   * Ask for the bands of the current orientation, in this plate's own pixels.
   *
   * The accepted solution when there is one, and otherwise the candidate
   * currently selected — the same rule the calculated spots follow, and for the
   * same reason: deciding between candidates is done by looking, and the bands
   * are one more thing to look at.
   *
   * Everything sent is already on screen: that orientation, the calibration
   * that indexed the pattern, the picked beam, and the image size. Nothing is
   * asked of the microscope — the overlay never leaves the pattern frame, so
   * the diffraction rotation and the parity, which one indexed pattern cannot
   * determine, are not needed and are not guessed.
   */
  function kikuchiSolution() {
    return state.accepted ?? state.solutions[state.selected] ?? null;
  }

  async function refreshKikuchi() {
    const size = frameSize();
    const accepted = kikuchiSolution();
    if (!accepted?.crystal_to_pattern || !state.picks.centre || !size) return;
    const values = state.solveForm.values();
    const request = {
      phase: values.phase,
      orientation: { crystal_to_pattern: accepted.crystal_to_pattern },
      units: values.units,
      camera_constant_mm_angstrom: values.camera_constant_mm_angstrom,
      pixel_size_mm: values.pixel_size_mm,
      reciprocal_per_px_angstrom: values.reciprocal_per_px_angstrom,
      centre_x: state.picks.centre[0],
      centre_y: state.picks.centre[1],
      frame_width: size.width,
      frame_height: size.height,
      target_zone_axis: state.tiltForm?.values().target_zone_axis ?? [0, 0, 0],
    };
    const key = JSON.stringify(request);
    if (key === state.kikuchiRequest && state.kikuchi) {
      drawPattern();
      return;
    }
    try {
      const result = await call('tem.kikuchi_overlay', request);
      state.kikuchiRequest = key;
      state.kikuchi = result;
      drawPattern();
      const data = result.data;
      const route = data.connecting
        ? ` · ${data.connecting.text}` +
          (data.connecting.waypoints.length
            ? ` via ${data.connecting.waypoints.map((way) => way.label).join(', ')}`
            : '')
        : data.connecting_note
          ? ` · ${data.connecting_note}`
          : '';
      frame.setStatus(
        `${data.bands.length} predicted band(s), each as wide as its own 000→g distance` +
          `${route} · positions and widths are geometry; contrast, excess/deficient sides and ` +
          'HOLZ lines are not modelled, and a thin foil may show these spots with no bands at all',
      );
    } catch (error) {
      // The bands are an aid; failing to draw them must not disturb the picking
      // and indexing the user is actually doing.
      state.kikuchi = null;
      state.showKikuchi = false;
      kikuchiButton.setAttribute('aria-pressed', 'false');
      drawPattern();
      frame.setStatus(error?.message ?? String(error));
    }
  }

  /**
   * Draw the bands: fine dotted edges under everything else.
   *
   * Two rules the geometry dictates. Bands are drawn *through* the transmitted
   * beam but labelled well away from it, because at an exact zone axis every
   * band of the zone crosses at 000 — the most crowded and least informative
   * point of the figure, and where the beam marker and the picks live. And the
   * connecting band is drawn distinctly, because it is the one instruction here
   * that survives the missing holder calibration: following a band is a
   * pattern-frame move, where "tilt alpha by +12.3" is not.
   */
  function drawKikuchiBands(root, width, height) {
    const data = state.kikuchi?.data;
    if (!data?.bands?.length) return;
    const group = svg('g', { 'pointer-events': 'none' });
    const weight = Math.max(width, height) / 900;
    const font = Math.max(width, height) / 52;
    for (const band of data.bands) {
      const colour = band.connecting ? KIKUCHI_ROUTE_COLOUR : KIKUCHI_COLOUR;
      const prominence = 0.35 + 0.45 * Math.min(1, Number(band.intensity) || 0);
      const [[x1, y1], [x2, y2]] = band.centre;
      group.append(
        svg('line', {
          x1, y1, x2, y2,
          stroke: colour,
          'stroke-opacity': band.connecting ? 0.85 : prominence * 0.6,
          'stroke-width': band.connecting ? weight * 1.6 : weight * 0.9,
          'stroke-dasharray': `${weight * 6} ${weight * 6}`,
        }),
      );
      for (const edge of band.edges) {
        for (const run of edge) {
          group.append(
            svg('polyline', {
              points: run.map(([x, y]) => `${x},${y}`).join(' '),
              fill: 'none',
              stroke: colour,
              'stroke-opacity': prominence,
              'stroke-width': band.connecting ? weight * 2 : weight * 1.2,
              'stroke-dasharray': `${weight * 2} ${weight * 3}`,
            }),
          );
        }
      }
      const [labelX, labelY] = band.label_at;
      const text = band.connecting && data.connecting ? data.connecting.text : band.label;
      group.append(
        svg('text', {
          x: labelX, y: labelY,
          'font-size': band.connecting ? font * 1.1 : font,
          'font-weight': band.connecting ? '600' : '400',
          fill: colour,
          stroke: HALO_COLOUR,
          'stroke-width': font / 6,
          'paint-order': 'stroke',
          'text-anchor': 'middle',
          'dominant-baseline': 'middle',
          text,
        }),
      );
    }
    root.append(group);
  }

  function scheduleStereogram() {
    clearTimeout(state.stereoPending);
    state.stereoPending = setTimeout(() => refreshStereogram(), 180);
  }

  async function refreshStereogram() {
    if (!stereogramOperation || !state.tiltForm) return;
    const values = state.tiltForm.values();
    const request = {
      phase: values.phase,
      zone_axis: values.current_zone_axis,
      target_zone_axis: values.target_zone_axis,
      alpha_deg: values.alpha_deg,
      beta_deg: values.beta_deg,
      alpha_limit_deg: values.alpha_limit_deg,
      beta_limit_deg: values.beta_limit_deg,
      beam_rotation_deg: values.beam_rotation_deg,
    };
    const key = JSON.stringify(request);
    if (key === state.stereoRequest && state.stereo) return;
    try {
      const result = await call('tem.stereogram', request);
      state.stereoRequest = key;
      state.stereo = result;
      drawStereogram();
    } catch (error) {
      // A stereogram that cannot be drawn must not block picking or indexing,
      // which is what the user is actually doing. It says so and stays empty.
      state.stereo = null;
      stereoFrame.setContent(
        el('div.stage__placeholder', { text: 'No stereogram for these inputs yet.' }),
      );
      stereoFrame.setStatus(error?.message ?? String(error));
    }
  }

  /** Inverse stereographic projection: a point on the drawing to a holder direction. */
  function holderDirectionAt(x, y) {
    const squared = x * x + y * y;
    const scale = 1.0 / (1.0 + squared);
    return [2 * x * scale, 2 * y * scale, (1 - squared) * scale];
  }

  /**
   * The stage reading that brings a holder direction onto the beam.
   *
   * The principal branch of `pytex.tem.navigation.solve_tilts_for_direction`,
   * written out here so the cursor readout answers while the pointer moves
   * rather than a request later. The browser test compares it against the
   * server's own value at every plotted pole, so the two cannot drift.
   */
  function stageAnglesFor([x, y, z]) {
    const rho = Math.hypot(x, z);
    if (rho < 1e-9) return { alpha: y > 0 ? 90 : -90, beta: 0 };
    return {
      alpha: (Math.atan2(y, rho) * 180) / Math.PI,
      beta: (Math.atan2(-x, z) * 180) / Math.PI,
    };
  }

  /** The plotted pole nearest a point on the drawing, and how far away it is. */
  function poleNear(x, y) {
    const axes = state.stereo?.data.axes ?? [];
    let best = null;
    let bestDistance = Infinity;
    for (const entry of axes) {
      const distance = Math.hypot(entry.x - x, entry.y - y);
      if (distance < bestDistance) {
        bestDistance = distance;
        best = entry;
      }
    }
    return best ? { entry: best, distance: bestDistance } : null;
  }

  function drawStereogram() {
    const result = state.stereo;
    if (!result) {
      stereoFrame.setContent(
        el('div.stage__placeholder', {
          text: 'Index a pattern, or set a zone axis below, to draw the stereogram.',
        }),
      );
      stereoFrame.setStatus('');
      return;
    }
    const data = result.data;
    // A little past the unit circle, so a label on a pole at the rim is not cut
    // in half by the edge of the drawing.
    const root = svg('svg', {
      viewBox: '-1.28 -1.28 2.56 2.56',
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': 'Stereogram',
    });
    const point = (x, y) => `${x},${-y}`;

    root.append(
      svg('circle', {
        cx: 0, cy: 0, r: 1,
        fill: 'none',
        stroke: 'currentColor',
        'stroke-opacity': 0.45,
        'stroke-width': 0.006,
      }),
      svg('line', {
        x1: -1, y1: 0, x2: 1, y2: 0,
        stroke: 'currentColor', 'stroke-opacity': 0.15, 'stroke-width': 0.004,
      }),
      svg('line', {
        x1: 0, y1: -1, x2: 0, y2: 1,
        stroke: 'currentColor', 'stroke-opacity': 0.15, 'stroke-width': 0.004,
      }),
      // The two holder axes, named where they point. Without them the drawing is
      // a circle of poles with no way to tell which way the stage moves.
      svg('text', {
        x: 0, y: -1.12, 'font-size': 0.075, fill: 'currentColor', 'fill-opacity': 0.65,
        'text-anchor': 'middle', text: '+α',
      }),
      svg('text', {
        x: -1.1, y: 0.03, 'font-size': 0.075, fill: 'currentColor', 'fill-opacity': 0.65,
        'text-anchor': 'middle', text: '+β',
      }),
    );

    if (state.showEnvelope) {
      const boundary = data.envelope.boundary.map(([x, y]) => point(x, y)).join(' ');
      root.append(
        svg('polygon', {
          points: boundary,
          fill: 'var(--accent)',
          'fill-opacity': 0.08,
          stroke: 'var(--accent)',
          'stroke-opacity': 0.55,
          'stroke-width': 0.006,
          'stroke-dasharray': '0.02 0.014',
        }),
      );
    }

    // The route, before the poles, so a waypoint marker sits on top of its path.
    if (data.path) {
      root.append(
        svg('polyline', {
          points: data.path.points.map(([x, y]) => point(x, y)).join(' '),
          fill: 'none',
          stroke: BASIS_COLOUR,
          'stroke-width': 0.012,
          'stroke-dasharray': '0.03 0.024',
          'stroke-linecap': 'round',
        }),
      );
      for (const waypoint of data.path.waypoints) {
        const node = svg('circle', {
          cx: waypoint.x, cy: -waypoint.y, r: 0.028,
          fill: 'none',
          stroke: BASIS_COLOUR,
          'stroke-width': 0.01,
        });
        root.append(node);
        root.append(
          svg('text', {
            x: waypoint.x + 0.04, y: -waypoint.y + 0.062,
            'font-size': 0.058,
            fill: BASIS_COLOUR,
            'paint-order': 'stroke',
            stroke: 'var(--bg-raised)',
            'stroke-width': 0.022,
            'stroke-linejoin': 'round',
            text: waypoint.label,
          }),
        );
        stereoFrame.hoverable(node, {
          'Stop at': waypoint.label,
          'α / °': waypoint.alpha_deg,
          'β / °': waypoint.beta_deg,
          Reachable: waypoint.reachable,
        });
      }
    }

    for (const entry of data.axes) {
      const radius = entry.labelled ? 0.019 : 0.009;
      const node = svg('circle', {
        cx: entry.x, cy: -entry.y, r: radius,
        fill: entry.reachable ? 'var(--accent)' : 'currentColor',
        'fill-opacity': entry.reachable ? 0.95 : 0.35,
      });
      root.append(node);
      if (state.showPoleLabels && entry.labelled) {
        root.append(
          svg('text', {
            x: entry.x + 0.026, y: -entry.y - 0.022,
            'font-size': 0.052,
            fill: 'currentColor',
            // Poles crowd towards the rim, where labels overlap each other and
            // the ticks. The halo is what keeps a name readable where it does.
            'paint-order': 'stroke',
            stroke: 'var(--bg-raised)',
            'stroke-width': 0.02,
            'stroke-linejoin': 'round',
            text: entry.label,
          }),
        );
      }
      stereoFrame.hoverable(node, {
        Pole: entry.label,
        'From beam / °': entry.angle_from_beam_deg,
        'α / °': entry.alpha_deg,
        'β / °': entry.beta_deg,
        'Δα / °': entry.delta_alpha_deg,
        'Δβ / °': entry.delta_beta_deg,
        Reachable: entry.reachable,
      });
    }

    if (data.target) {
      root.append(
        svg('circle', {
          cx: data.target.x, cy: -data.target.y, r: 0.045,
          fill: 'none',
          stroke: data.target.reachable ? 'var(--ok, #2f9e63)' : 'var(--danger, #d1495b)',
          'stroke-width': 0.014,
        }),
      );
    }

    // The beam last, so nothing is drawn over where you are.
    root.append(
      svg('circle', {
        cx: data.beam.x, cy: -data.beam.y, r: 0.034,
        fill: 'none',
        stroke: 'currentColor',
        'stroke-width': 0.012,
      }),
      svg('line', {
        x1: data.beam.x - 0.06, y1: -data.beam.y, x2: data.beam.x + 0.06, y2: -data.beam.y,
        stroke: 'currentColor', 'stroke-width': 0.006,
      }),
      svg('line', {
        x1: data.beam.x, y1: -data.beam.y - 0.06, x2: data.beam.x, y2: -data.beam.y + 0.06,
        stroke: 'currentColor', 'stroke-width': 0.006,
      }),
    );

    stereoFrame.configure({
      // The drawing's y runs up while the SVG's runs down, so the mapping is
      // where that sign is undone — once, here, rather than at every readout.
      toData: (x, y) => ({ x, y: -y }),
      formatCursor: (position) => {
        const radius = Math.hypot(position.x, position.y);
        if (radius > 1.02) return 'outside the hemisphere';
        const direction = holderDirectionAt(position.x, position.y);
        const { alpha, beta } = stageAnglesFor(direction);
        const rho = (2 * Math.atan(radius) * 180) / Math.PI;
        const near = poleNear(position.x, position.y);
        const pole =
          near && near.distance < 0.045 ? ` · near ${near.entry.label}` : '';
        return (
          `α ${formatNumber(alpha, 1)}° · β ${formatNumber(beta, 1)}° · ` +
          `${formatNumber(rho, 1)}° from the holder axis${pole}`
        );
      },
    });

    stereoFrame.setContent(root, { preserveViewport: true });
    const reachable = data.axes.filter((entry) => entry.reachable).length;
    stereoFrame.setStatus(
      `${data.zone_axis_label} on the beam · ${data.axes.length} poles, ${reachable} reachable` +
        (data.target ? ` · route to ${data.target.label}` : '') +
        ' · hover for the tilt that reaches a point',
    );
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
      // An arrow to a node with no pick on it is a different claim from an arrow
      // to a picked spot — it says "the lattice puts a node here and you have
      // not picked it" — so it is drawn dashed rather than solid. Drawn
      // identically, as it was, it read as a pick that had landed in empty
      // space, which is precisely the thing it is not.
      const dashed = !vector.on_a_pick;
      for (const [colour, weight, opacity] of [
        [HALO_COLOUR, marker / 2, 0.6],
        [BASIS_COLOUR, marker / 5, 1],
      ]) {
        group.append(
          svg('line', {
            x1, y1, x2: tipX, y2: tipY,
            stroke: colour,
            'stroke-width': weight,
            'stroke-opacity': opacity,
            'stroke-linecap': 'round',
            'stroke-dasharray': dashed ? `${marker} ${marker}` : null,
          }),
        );
      }
      group.append(
        svg('polygon', {
          points,
          fill: dashed ? 'none' : BASIS_COLOUR,
          stroke: BASIS_COLOUR,
          'stroke-width': marker / 5,
        }),
      );
      group.append(
        svg('text', {
          x: x1 + dx * 0.55 - uy * marker * 1.6,
          y: y1 + dy * 0.55 + ux * marker * 1.6,
          'font-size': marker * 2.2,
          'font-style': 'italic',
          fill: BASIS_COLOUR,
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
      // A new indexing supersedes the previous orientation, and with it the
      // bands that were drawn from it.
      resetKikuchi();
      calculatedButton.hidden = state.solutions.length === 0;
      // The bands need an orientation, and every candidate carries one, so the
      // toggle arrives with the indexing rather than waiting for an acceptance.
      // Which candidate it draws is the one selected below — looking at its
      // bands is part of judging it, exactly as the calculated spots are.
      kikuchiButton.hidden = !state.solutions[0]?.crystal_to_pattern;
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
                // The bands belong to the candidate being considered, so they
                // follow the selection rather than staying on the last one.
                scheduleKikuchi();
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
    // Already visible from the indexing; accepting only settles which
    // orientation the bands are drawn from.
    kikuchiButton.hidden = !entry.crystal_to_pattern;
    scheduleKikuchi();
    state.tiltForm.setValues({ phase, current_zone_axis: axis });
    // `setValues` is a programmatic write and fires no change event, so the
    // stereogram is told explicitly. Without this the solved axis appeared in
    // the form while the drawing beside it still showed the previous one.
    scheduleStereogram();
    frame.setStatus(
      `Accepted ${entry.phase} down ${entry.zone_axis} (score ` +
        `${formatNumber(entry.score, 3)}). It is now the starting orientation for the zone-axis ` +
        'list and the tilt plan below.',
    );
    openStep(stepThree);
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
    scheduleStereogram();
    frame.setStatus(`Target set to ${label ?? indices.join(' ')} — plan the tilt below.`);
    openStep(stepFour);
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
  renderPickTool();
  drawPattern();
  drawStereogram();
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
