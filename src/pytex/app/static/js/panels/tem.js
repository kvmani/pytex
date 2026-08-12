/**
 * The TEM solver panel: upload, calibrate, pick, index, tilt.
 *
 * The picking surface is the whole point. A diffraction pattern is indexed by
 * clicking on it, and everything a user needs while clicking — which pick is
 * the beam, how far each spot is from it, what d-spacing that corresponds to —
 * has to be visible without leaving the image. So the pattern is the plot, the
 * cursor readout reports the distance from the marked beam in both millimetres
 * and inverse angstroms, and every pick is a hoverable entity like any other.
 *
 * The image never leaves the browser. Only the picked coordinates are sent,
 * which keeps an unpublished micrograph on the machine it was opened on and
 * makes the request small enough to be instant on a slow intranet.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';
import { call } from '../core/api.js';

export const panel = {
  id: 'tem',
  title: 'TEM Solver',
  tagline: 'Pick a pattern, index it, and plan the tilt to the next zone axis.',
};

export function mount(context) {
  const operations = context.manifest.operations.filter((entry) => entry.panel === panel.id);
  const solveOperation = operations.find((entry) => entry.id === 'tem.solve_pattern');
  const tiltOperation = operations.find((entry) => entry.id === 'tem.plan_tilt');
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const state = {
    mode: 'tilt', // 'pick' once an image is loaded
    image: null,
    picks: { centre: null, spots: [] },
    solveForm: null,
    tiltForm: null,
    teaches: null,
    solution: null,
  };

  const frame = plotFrame({
    title: 'Pattern',
    toolbar: [
      el('button.button', {
        type: 'button',
        text: 'Undo pick',
        onclick: () => {
          state.picks.spots.pop();
          drawPattern();
        },
      }),
      el('button.button', {
        type: 'button',
        text: 'Clear picks',
        onclick: () => {
          state.picks = { centre: null, spots: [] };
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

  const solveHost = el('div');
  const tiltHost = el('div');

  const solveButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Index the pattern',
    onclick: () => solve(),
  });
  const tiltButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Plan the tilt',
    onclick: () => planTilt(),
  });

  context.rail.append(
    el('details.group', { open: true }, [
      el('summary', { text: '1 · Pattern and picks' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'Open a pattern image, then click the transmitted beam first and the reflections ' +
            'after. The image stays on this machine; only the coordinates are sent.',
        }),
        fileInput,
        solveHost,
        solveButton,
      ]),
    ]),
    el('details.group', { open: true }, [
      el('summary', { text: '2 · Tilt to the next zone axis' }),
      el('div.group__body', {}, [tiltHost, tiltButton]),
    ]),
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'Tilt plans you can run without a pattern of your own.',
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

  function renderForms(solveInitial = {}, tiltInitial = {}) {
    state.solveForm = buildForm(solveOperation, {
      initial: solveInitial,
      // The picks live on the canvas, not in a text box, so the parameter is
      // hidden from the generated form and supplied at call time.
    });
    hideField(state.solveForm.element, 'picks');
    solveHost.replaceChildren(state.solveForm.element);

    state.tiltForm = buildForm(tiltOperation, { initial: tiltInitial });
    tiltHost.replaceChildren(state.tiltForm.element);
  }

  function hideField(root, name) {
    for (const field of root.querySelectorAll('.field')) {
      const control = field.querySelector(`[id^="ctl-${name}-"]`);
      if (control) field.hidden = true;
    }
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderForms({}, example.request);
    planTilt();
  }

  /* --------------------------------------------------------------- image */

  function loadImage(file) {
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const image = new Image();
      image.onload = () => {
        state.image = { source: reader.result, width: image.width, height: image.height };
        state.mode = 'pick';
        state.picks = { centre: null, spots: [] };
        drawPattern();
      };
      image.src = reader.result;
    };
    reader.readAsDataURL(file);
  }

  /* -------------------------------------------------------------- picking */

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

  function drawPattern() {
    if (!state.image) {
      frame.setContent(
        el('div.stage__placeholder', {
          text: 'Open a pattern image to start picking, or run a tilt example below.',
        }),
      );
      frame.setStatus('');
      return;
    }
    const { width, height } = state.image;
    const root = svg('svg', {
      viewBox: `0 0 ${width} ${height}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': 'Diffraction pattern',
    });
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
      if (!state.picks.centre) state.picks.centre = [point.x, point.y];
      else state.picks.spots.push({ x: point.x, y: point.y });
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
      const result = await call('tem.solve_pattern', {
        ...state.solveForm.values(),
        picks: state.picks,
      });
      state.solution = result;
      renderResult(details, result);
      // Carry the answer into the tilt form: the axis just indexed is the axis
      // the next move starts from, and retyping it is both tedious and a chance
      // to get it wrong.
      state.tiltForm.setValues({
        phase: result.inputs.phase,
        current_zone_axis: result.data.zone_axis,
      });
      frame.setStatus(
        `Indexed as ${result.data.phase_name} down ${result.data.zone_axis_label}` +
          ` · ${state.picks.spots.length} picks`,
      );
    } catch (error) {
      if (!state.solveForm.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      solveButton.disabled = false;
      solveButton.textContent = 'Index the pattern';
    }
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
    } catch (error) {
      if (!state.tiltForm.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      tiltButton.disabled = false;
      tiltButton.textContent = 'Plan the tilt';
    }
  }

  renderForms();
  drawPattern();
  if (examples.length) loadExample(examples[0]);

  return { help: () => (state.mode === 'pick' ? solveOperation : tiltOperation) };
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
