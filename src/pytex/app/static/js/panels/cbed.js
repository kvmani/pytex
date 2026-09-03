/**
 * CBED: the disc pattern, its fringes, and what they measure.
 *
 * Three views over one technique, chosen with the picker in the rail:
 *
 * - **CBED pattern** — the simulated exposure. The discs arrive as one
 *   rasterised intensity image with a stated extent in millimetres, and the
 *   panel draws the outlines, labels and HOLZ rings over it as vectors. That
 *   split is deliberate: the geometry has to stay measurable under the cursor
 *   and hoverable per disc, while the intensity inside a disc is a picture and
 *   pretending otherwise would throw away the fringes that are the measurement.
 * - **Thickness from fringes** — the two-beam inversion, drawn as the straight
 *   line it actually is, so a misassigned fringe order shows up as a point off
 *   the line rather than as a wrong number.
 * - **HOLZ rings** — where the higher-order Laue zones fall, which is how the
 *   repeat *along* the beam gets measured.
 *
 * The greyscale bytes are turned into an image through a canvas rather than
 * being sent as a PNG. It keeps the Python side free of an image encoder, and
 * it is what lets the contrast control re-map the same bytes without another
 * round trip — a CBED user adjusts contrast constantly, because the fringes at
 * the rim of a disc are orders of magnitude fainter than its centre.
 */

import { call } from '../core/api.js';
import { explainer } from '../core/explainer.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import * as log from '../core/logbook.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'cbed',
  title: 'CBED',
  tagline: 'Convergent-beam discs, the fringes inside them, and the thickness they measure.',
};

/** The operations this panel presents, in the order the picker lists them. */
const VIEWS = ['cbed.pattern', 'cbed.thickness_from_fringes', 'cbed.holz_rings'];

/** Side of the SVG the pattern is drawn into, in user units. */
const SIZE = 720;

const DEFAULT_APPEARANCE = Object.freeze({
  gamma: 0.6,
  outlines: true,
  labels: true,
  holz: true,
  invert: false,
});

export function mount(context) {
  const operations = VIEWS.map((id) =>
    context.manifest.operations.find((entry) => entry.id === id),
  ).filter(Boolean);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const state = {
    operation: operations[0],
    result: null,
    teaches: null,
    form: null,
    appearance: { ...DEFAULT_APPEARANCE },
    /** Decoded intensity bytes, kept so contrast re-maps without another call. */
    image: null,
  };

  const frame = plotFrame({ title: 'CBED pattern', units: 'mm', digits: 2 });
  const details = el('div');
  const formHost = el('div');
  const appearanceHost = el('div');

  const viewSelect = el(
    'select',
    {
      'aria-label': 'View',
      onchange: () => {
        state.operation = operations.find((entry) => entry.id === viewSelect.value);
        state.teaches = null;
        state.result = null;
        state.image = null;
        renderControls();
        renderAppearance();
        updateRunLabel();
      },
    },
    operations.map((entry) =>
      el('option', { value: entry.id, text: entry.title, title: entry.summary }),
    ),
  );

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Simulate CBED pattern',
    onclick: () => run(),
  });

  context.rail.append(
    el('div.field', {}, [
      el('label.field__label', { text: 'View' }),
      viewSelect,
      explainer(
        'The pattern itself, the thickness its fringes measure, and the Laue-zone rings that ' +
          'measure the repeat along the beam.',
        { label: 'What each view is' },
      ),
    ]),
    formHost,
    runButton,
    appearanceHost,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        explainer(
          'The first two are the same crystal and thickness at two convergence angles — the ' +
            'one parameter that decides whether CBED can measure a thickness at all.',
          { label: 'What these examples show' },
        ),
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
  context.stage.append(frame.element, details);

  renderControls();
  renderAppearance();
  // Draw immediately, as every other panel does. An empty stage waiting to be
  // pressed teaches nothing, and the default pattern is the one that shows what
  // the convergence-angle control is for.
  run();

  function updateRunLabel(busy = false) {
    const labels = {
      'cbed.pattern': ['Simulate CBED pattern', 'Simulating…'],
      'cbed.thickness_from_fringes': ['Fit the thickness', 'Fitting…'],
      'cbed.holz_rings': ['Find the Laue-zone rings', 'Searching…'],
    };
    const pair = labels[state.operation.id] ?? ['Run', 'Running…'];
    runButton.textContent = busy ? pair[1] : pair[0];
  }

  function renderControls(initial = {}) {
    state.form = buildForm(state.operation, { initial });
    formHost.replaceChildren(state.form.element);
  }

  function renderAppearance() {
    if (state.operation.id !== 'cbed.pattern') {
      appearanceHost.replaceChildren();
      return;
    }
    appearanceHost.replaceChildren(
      appearanceControls(state.appearance, () => {
        if (state.result) draw();
      }),
    );
  }

  function loadExample(example) {
    state.operation = operations.find((entry) => entry.id === example.operation);
    viewSelect.value = state.operation.id;
    state.teaches = example.teaches;
    state.image = null;
    renderControls(example.request);
    renderAppearance();
    updateRunLabel();
    run();
  }

  async function run() {
    runButton.disabled = true;
    updateRunLabel(true);
    state.form.clearErrors();
    try {
      state.result = await call(state.operation.id, state.form.values());
      state.image = null;
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      updateRunLabel(false);
    }
  }

  function draw() {
    if (state.operation.id === 'cbed.pattern') drawPattern();
    else if (state.operation.id === 'cbed.thickness_from_fringes') drawThicknessFit();
    else drawRings();
  }

  /* ------------------------------------------------------------ the pattern */

  function drawPattern() {
    const data = state.result.data;
    const extent = data.image.extent_mm;
    // Data millimetres to SVG user units. The image spans [-extent, +extent] in
    // both axes, and y grows downward on screen while it grows upward in the
    // detector frame — so a sign lives here and nowhere else.
    const scale = SIZE / (2 * extent);
    const toSvgX = (mm) => SIZE / 2 + mm * scale;
    const toSvgY = (mm) => SIZE / 2 - mm * scale;

    if (!state.image) state.image = decodeGray(data.image);
    const href = toDataUrl(state.image, state.appearance);

    const children = [
      svg('image', {
        href,
        x: 0,
        y: 0,
        width: SIZE,
        height: SIZE,
        // The pattern is sampled data, not a photograph: smoothing it would
        // invent fringe contrast between samples that the simulation never
        // computed.
        style: 'image-rendering: pixelated',
      }),
    ];

    if (state.appearance.holz) {
      for (const ring of data.holz_rings ?? []) {
        const radius = ring.radius_mm * scale;
        if (!Number.isFinite(radius) || radius <= 0) continue;
        children.push(
          svg('circle', {
            cx: SIZE / 2,
            cy: SIZE / 2,
            r: radius,
            fill: 'none',
            stroke: 'var(--accent)',
            'stroke-width': 1.2,
            'stroke-dasharray': '6 5',
            opacity: 0.85,
          }),
          svg('text', {
            x: SIZE / 2,
            y: SIZE / 2 - radius - 5,
            'text-anchor': 'middle',
            fill: 'var(--accent)',
            'font-size': 13,
            text: `HOLZ ${ring.order}`,
          }),
        );
      }
    }

    const discRadius = data.disc_radius_mm * scale;
    const marks = [];
    for (const disc of data.discs) {
      const cx = toSvgX(disc.x_mm);
      const cy = toSvgY(disc.y_mm);
      const node = svg('circle', {
        cx,
        cy,
        r: discRadius,
        fill: 'transparent',
        stroke: state.appearance.outlines ? 'var(--accent)' : 'transparent',
        'stroke-width': 1,
        opacity: disc.hkl_label.includes('direct') ? 1 : 0.75,
      });
      // Hoverable whether or not the outline is drawn: hiding the annotation
      // must not take the measurements away with it.
      frame.hoverable(node, {
        Reflection: disc.hkl_label,
        'd / Å': disc.d_angstrom,
        '|g| / Å⁻¹': disc.g_inv_angstrom,
        'Centre x / mm': disc.x_mm,
        'Centre y / mm': disc.y_mm,
        'ξ_g / Å': disc.extinction_distance_angstrom,
        '|F| / Å': disc.structure_factor_amplitude,
        'Mean intensity': disc.mean_intensity,
      });
      marks.push(node);
      if (state.appearance.labels && !disc.hkl_label.includes('direct')) {
        marks.push(
          svg('text', {
            x: cx,
            y: cy - discRadius - 4,
            'text-anchor': 'middle',
            fill: 'var(--accent)',
            'font-size': 12,
            'pointer-events': 'none',
            text: disc.hkl_label,
          }),
        );
      }
    }

    const root = svg(
      'svg',
      { viewBox: `0 0 ${SIZE} ${SIZE}`, width: '100%', height: '100%' },
      [...children, ...marks],
    );
    frame.setContent(root);
    frame.configure({
      toData: (x, y) => ({ x: (x - SIZE / 2) / scale, y: (SIZE / 2 - y) / scale }),
      formatCursor: (point) => {
        const radius = Math.hypot(point.x, point.y);
        const g = radius / data.camera_constant_mm_angstrom;
        return (
          `${formatNumber(point.x, 2)}, ${formatNumber(point.y, 2)} mm · ` +
          `|g| ${formatNumber(g, 4)} Å⁻¹ · d ${g > 0 ? formatNumber(1 / g, 4) : '∞'} Å`
        );
      },
    });
    frame.setStatus(
      `${data.discs.length} discs · ${data.regime === 'kossel-moellenstedt' ? 'separated' : 'overlapping'} ` +
        `· disc radius ${formatNumber(data.disc_radius_mm, 2)} mm · nearest centres ` +
        `${formatNumber(data.nearest_disc_separation_mm, 2)} mm apart`,
    );
    frame.setOverlay(
      el('div.plot__caption', {}, [
        el('strong', { text: `${data.phase_name} ${data.zone_axis_label}` }),
        el('span', {
          text: `${data.method === 'bloch' ? 'Bloch wave' : 'Two-beam'} · ${data.regime}`,
        }),
      ]),
    );
  }

  /* ----------------------------------------------------------- the analyses */

  /**
   * The linearised thickness fit, drawn as the straight line it is.
   *
   * Plotting it matters more than it looks: the method's usual failure is
   * assigning the fringe orders wrongly, and a wrong assignment shows here as
   * points that curve away from the line while the fitted number still looks
   * plausible.
   */
  function drawThicknessFit() {
    const data = state.result.data;
    const points = data.fit_points ?? [];
    const margin = { left: 96, right: 32, top: 32, bottom: 72 };
    const width = 900;
    const height = 520;
    if (!points.length) {
      frame.setContent(el('p.field__help', {
                         text: 'No fringe minima to plot.',
                       }));
      return;
    }

    const xs = points.map((point) => point.inverse_order_squared);
    const ys = points.map((point) => point.s_over_n_squared);
    const xMax = Math.max(...xs, 0) * 1.1;
    const yMax = Math.max(...ys, 0) * 1.1;
    const toX = (value) => margin.left + (value / xMax) * (width - margin.left - margin.right);
    const toY = (value) => height - margin.bottom - (value / yMax) * (height - margin.top - margin.bottom);

    // The fitted line, from the reported thickness and extinction distance
    // rather than from a second regression here: two fits of one dataset are
    // two chances to disagree, and the drawing must show the reported answer.
    const intercept = 1 / (data.thickness_angstrom * data.thickness_angstrom);
    const slope = -1 / (data.extinction_distance_angstrom * data.extinction_distance_angstrom);

    const nodes = [
      axis(margin, width, height),
      svg('line', {
        x1: toX(0),
        y1: toY(intercept),
        x2: toX(xMax),
        y2: toY(intercept + slope * xMax),
        stroke: 'var(--accent)',
        'stroke-width': 2,
      }),
      ...points.map((point) => {
        const node = svg('circle', {
          cx: toX(point.inverse_order_squared),
          cy: toY(point.s_over_n_squared),
          r: 6,
          fill: 'var(--accent)',
        });
        frame.hoverable(node, {
          'Minimum n': point.order,
          's_n / Å⁻¹': point.excitation_error_inv_angstrom,
          '1/n²': point.inverse_order_squared,
          '(s_n/n)² / Å⁻²': point.s_over_n_squared,
        });
        return node;
      }),
      label(width / 2, height - 24, '1 / n²'),
      label(24, height / 2, '(s_n / n)²  /  Å⁻²', -90),
    ];

    frame.setContent(
      svg('svg', { viewBox: `0 0 ${width} ${height}`, width: '100%', height: '100%' }, nodes),
    );
    frame.configure({
      toData: (x, y) => ({
        x: ((x - margin.left) / (width - margin.left - margin.right)) * xMax,
        y: ((height - margin.bottom - y) / (height - margin.top - margin.bottom)) * yMax,
      }),
      formatCursor: (point) => `1/n² ${formatNumber(point.x, 4)} · (s/n)² ${formatNumber(point.y, 7)}`,
    });
    frame.setStatus(
      `Thickness ${formatNumber(data.thickness_nm, 1)} nm · extinction distance ` +
        `${formatNumber(data.extinction_distance_angstrom, 1)} Å · orders assigned from n = ${data.first_order}`,
    );
    frame.setOverlay(null);
  }

  /** The Laue-zone rings, to scale, around the pattern centre. */
  function drawRings() {
    const data = state.result.data;
    const rings = data.rings ?? [];
    const outer = Math.max(...rings.map((ring) => ring.radius_inv_angstrom), 1) * 1.2;
    const scale = SIZE / (2 * outer);

    const nodes = [
      svg('circle', { cx: SIZE / 2, cy: SIZE / 2, r: 4, fill: 'var(--ink)' }),
      ...rings.flatMap((ring) => {
        const radius = ring.radius_inv_angstrom * scale;
        const node = svg('circle', {
          cx: SIZE / 2,
          cy: SIZE / 2,
          r: radius,
          fill: 'none',
          stroke: 'var(--accent)',
          'stroke-width': 2,
        });
        frame.hoverable(node, {
          'Laue zone': ring.order,
          'Radius / Å⁻¹': ring.radius_inv_angstrom,
          'Radius / mm': ring.radius_mm,
        });
        return [
          node,
          label(SIZE / 2, SIZE / 2 - radius - 8, `HOLZ ${ring.order}`),
        ];
      }),
    ];

    frame.setContent(
      svg('svg', { viewBox: `0 0 ${SIZE} ${SIZE}`, width: '100%', height: '100%' }, nodes),
    );
    frame.configure({
      toData: (x, y) => ({ x: (x - SIZE / 2) / scale, y: (SIZE / 2 - y) / scale }),
      formatCursor: (point) => `|G| ${formatNumber(Math.hypot(point.x, point.y), 4)} Å⁻¹`,
    });
    frame.setStatus(
      `Layer spacing H = ${formatNumber(data.layer_spacing_inv_angstrom, 5)} Å⁻¹, which is a ` +
        `${formatNumber(data.real_space_repeat_angstrom, 4)} Å repeat along ${data.zone_axis_label}`,
    );
    frame.setOverlay(null);
  }

  return {
    help: () => state.operation,
  };
}

/* ---------------------------------------------------------------- rendering */

function axis(margin, width, height) {
  return svg('g', { stroke: 'var(--line-strong)', 'stroke-width': 1.5 }, [
    svg('line', {
      x1: margin.left,
      y1: height - margin.bottom,
      x2: width - margin.right,
      y2: height - margin.bottom,
    }),
    svg('line', {
      x1: margin.left,
      y1: margin.top,
      x2: margin.left,
      y2: height - margin.bottom,
    }),
  ]);
}

function label(x, y, text, rotate = 0) {
  return svg('text', {
    x,
    y,
    'text-anchor': 'middle',
    fill: 'var(--ink-muted)',
    'font-size': 14,
    transform: rotate ? `rotate(${rotate} ${x} ${y})` : null,
    text,
  });
}

/**
 * Decode the base64 greyscale payload into bytes.
 *
 * @param {object} image - `{width, height, encoding, data}` from the service.
 * @returns {{width: number, height: number, bytes: Uint8Array}}
 */
function decodeGray(image) {
  if (image.encoding !== 'base64-gray8') {
    log.error(`The CBED image arrived in an unknown encoding: ${image.encoding}.`, {
      source: 'cbed',
    });
    return { width: 1, height: 1, bytes: new Uint8Array(1) };
  }
  const binary = atob(image.data);
  const bytes = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) bytes[index] = binary.charCodeAt(index);
  return { width: image.width, height: image.height, bytes };
}

/**
 * Re-map the intensity bytes through the contrast curve and return a data URL.
 *
 * The gamma is applied here rather than server-side because it is a *viewing*
 * choice: a CBED user moves it constantly to bring out the faint outer fringes,
 * and each move would otherwise be a re-simulation.
 */
function toDataUrl({ width, height, bytes }, appearance) {
  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  const context = canvas.getContext('2d');
  const pixels = context.createImageData(width, height);
  const gamma = Math.max(appearance.gamma, 0.05);
  // 256 entries computed once rather than a pow() per pixel: at 512² that is a
  // quarter of a million calls saved on every contrast change.
  const curve = new Uint8Array(256);
  for (let level = 0; level < 256; level += 1) {
    const value = Math.round(255 * (level / 255) ** gamma);
    curve[level] = appearance.invert ? 255 - value : value;
  }
  for (let index = 0; index < bytes.length; index += 1) {
    const value = curve[bytes[index]];
    const offset = index * 4;
    pixels.data[offset] = value;
    pixels.data[offset + 1] = value;
    pixels.data[offset + 2] = value;
    pixels.data[offset + 3] = 255;
  }
  context.putImageData(pixels, 0, 0);
  return canvas.toDataURL('image/png');
}

function appearanceControls(appearance, onChange) {
  const change = (key) => (event) => {
    const target = event.currentTarget;
    appearance[key] = target.type === 'checkbox' ? target.checked : Number(target.value);
    onChange();
  };

  return el('details.group', { open: true }, [
    el('summary', { text: 'Display' }),
    el('div.group__body', {}, [
      el('div.field', {}, [
        el('label.field__label', { for: 'cbed-gamma', text: 'Contrast (gamma)' }),
        el('input', {
          id: 'cbed-gamma',
          type: 'range',
          min: 0.2,
          max: 2,
          step: 0.05,
          value: appearance.gamma,
          oninput: change('gamma'),
        }),
        explainer(
          'Below 1 lifts the faint outer fringes, which are orders of magnitude weaker than ' +
            'the disc centres. It changes the display only; the numbers in the table are the ' +
            'computed intensities.',
          { label: 'How the contrast curve works' },
        ),
      ]),
      toggle('cbed-outlines', 'Disc outlines', appearance.outlines, change('outlines')),
      toggle('cbed-labels', 'Reflection labels', appearance.labels, change('labels')),
      toggle('cbed-holz', 'HOLZ rings', appearance.holz, change('holz')),
      toggle('cbed-invert', 'Invert (plate-like)', appearance.invert, change('invert')),
    ]),
  ]);
}

function toggle(id, text, checked, onChange) {
  const input = el('input', { id, type: 'checkbox', onchange: onChange });
  input.checked = checked;
  return el('div.checkbox', {}, [input, el('label', { for: id, text })]);
}
