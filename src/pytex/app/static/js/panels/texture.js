/**
 * The texture panel: a polycrystal seen the three ways texture is read.
 *
 * A pole figure is a continuous density and is drawn as one — filled contour
 * bands rather than a scatter of points, because the question asked of it is
 * "how strong, and where", and a scatter answers neither. The inverse pole
 * figure is a scatter, because it genuinely is one point per grain and drawing
 * it as a density would need a second kernel choice on top of the ODF's.
 *
 * Everything is in multiples of a random distribution, and the colour scale is
 * anchored at 1 m.r.d. rather than at the minimum of the data: 1 is where a
 * texture-free material sits, so a figure whose scale started at its own
 * minimum would make an untextured sample look textured. That anchoring is the
 * single most important thing about how a pole figure is coloured, and it is
 * the thing most often got wrong.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';
import { call } from '../core/api.js';

export const panel = {
  id: 'texture',
  title: 'Texture',
  tagline: 'Pole figures, inverse pole figures and ODF sections of a model texture.',
};

/** The views, in the order the picker offers them. */
const VIEWS = [
  'texture.pole_figure',
  'texture.inverse_pole_figure',
  'texture.odf_sections',
];

const VIEW = 100;

/**
 * The m.r.d. colour ramp.
 *
 * Sequential and perceptually ordered, running cool through warm, so that
 * "more" always reads as "warmer" without a legend. It is anchored so that
 * exactly 1 m.r.d. — the random baseline — sits at the pale middle of the ramp:
 * below 1 is cool, above 1 is warm, and an untextured material comes out a flat
 * neutral rather than a dramatic colour that happens to be its own minimum.
 */
const RAMP = [
  [0.0, [49, 84, 140]],
  [0.5, [90, 148, 194]],
  [1.0, [232, 234, 228]],
  [2.0, [233, 178, 84]],
  [4.0, [214, 106, 43]],
  [8.0, [154, 32, 32]],
];

function rampColor(mrd) {
  const value = Math.max(mrd, 0);
  for (let index = 1; index < RAMP.length; index += 1) {
    const [lowStop, low] = RAMP[index - 1];
    const [highStop, high] = RAMP[index];
    if (value <= highStop || index === RAMP.length - 1) {
      const span = highStop - lowStop || 1;
      const t = Math.min(Math.max((value - lowStop) / span, 0), 1);
      const mix = low.map((channel, c) => Math.round(channel + (high[c] - channel) * t));
      return `rgb(${mix[0]} ${mix[1]} ${mix[2]})`;
    }
  }
  return 'rgb(154 32 32)';
}

export function mount(context) {
  const operations = VIEWS.map((id) =>
    context.manifest.operations.find((entry) => entry.id === id),
  ).filter(Boolean);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = { operation: operations[0], result: null, teaches: null, form: null };

  const frame = plotFrame({ title: 'Pole figure', toolbar: [] });
  const legend = el('div.legend');
  const details = el('div');
  const formHost = el('div');

  const viewSelect = el(
    'select',
    {
      'aria-label': 'View',
      onchange: () => {
        state.operation = operations.find((entry) => entry.id === viewSelect.value);
        state.teaches = null;
        renderControls(carryOver());
        run();
      },
    },
    operations.map((entry) =>
      el('option', { value: entry.id, text: entry.title, title: entry.summary }),
    ),
  );

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Build texture',
    onclick: () => run(),
  });

  context.rail.append(
    el('div.field', {}, [
      el('label.field__label', { text: 'View' }),
      viewSelect,
      el('p.field__help', {
        text: 'The same texture seen three ways: where a plane points, which direction is along a specimen axis, and the orientation density itself.',
      }),
    ]),
    formHost,
    runButton,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'Run the random baseline first: it is flat at 1 m.r.d., which is what every other figure here is measured against.',
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

  /**
   * The texture-defining settings, carried across a change of view.
   *
   * All three views are of *the same texture*; only the question differs. A
   * user who has set up a 15-degree brass texture and switches from the pole
   * figure to the ODF means to see that texture's ODF, not the default one.
   */
  function carryOver() {
    if (!state.form) return {};
    const values = state.form.values();
    const carried = {};
    for (const key of ['phase', 'model', 'spread_deg', 'grain_count', 'halfwidth_deg', 'seed']) {
      if (values[key] !== undefined) carried[key] = values[key];
    }
    return carried;
  }

  function renderControls(initial = {}) {
    state.form = buildForm(state.operation, { initial });
    formHost.replaceChildren(state.form.element);
  }

  function loadExample(example) {
    state.operation = operations.find((entry) => entry.id === example.operation);
    viewSelect.value = state.operation.id;
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Building…';
    state.form.clearErrors();
    try {
      const result = await call(state.operation.id, state.form.values());
      state.result = result;
      draw();
      renderResult(details, result, { teaches: state.teaches });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Build texture';
    }
  }

  function draw() {
    const data = state.result.data;
    const kind = state.operation.id;
    frame.element.hidden = false;
    legend.hidden = false;

    if (kind === 'texture.odf_sections') {
      frame.configure({ toData: null, formatCursor: null });
      frame.setContent(renderSections(data, frame));
      renderRampLegend(data.max_mrd);
      frame.setStatus(
        `three φ₂ sections · peak ${formatNumber(data.max_mrd, 2)} m.r.d. · ` +
          `${data.grain_count} grains · φ₁ across, Φ down`,
      );
      return;
    }

    const isPoleFigure = kind === 'texture.pole_figure';
    frame.configure({
      toData: (x, y) => ({ x: x / VIEW, y: -y / VIEW }),
      formatCursor: (point) => {
        const radius = Math.hypot(point.x, point.y);
        if (radius > 1.0001) return 'outside the projection';
        const polar =
          data.projection === 'equal_area'
            ? 2 * Math.asin(Math.min(radius / Math.SQRT2, 1))
            : 2 * Math.atan(radius);
        const azimuth = ((Math.atan2(point.y, point.x) * 180) / Math.PI + 360) % 360;
        return (
          `${formatNumber((polar * 180) / Math.PI, 1)}° from ` +
          `${isPoleFigure ? 'ND' : '[001]'} · ${formatNumber(azimuth, 1)}° azimuth`
        );
      },
    });

    if (isPoleFigure) {
      frame.setContent(renderDensity(data, frame));
      renderRampLegend(data.max_mrd);
      frame.setStatus(
        `${data.pole_label} · peak ${formatNumber(data.max_mrd, 2)} m.r.d. · ` +
          `mean ${formatNumber(data.mean_mrd, 3)} m.r.d. (1 by construction) · ` +
          `${data.grain_count} grains · hover for the intensity at a point`,
      );
    } else {
      frame.setContent(renderScatter(data, frame));
      legend.replaceChildren(
        el('span.legend__item', {}, [
          el('span', {
            text: `One point per grain — ${data.grain_count} grains, folded into the standard triangle.`,
          }),
        ]),
      );
      frame.setStatus(
        `crystal direction along ${data.axis_label} · ${data.grain_count} grains · ` +
          'hover a point for its grain and nearest direction',
      );
    }
  }

  /** The colour key, as bands with their m.r.d. values. */
  function renderRampLegend(maxMrd) {
    const stops = [0.25, 0.5, 1, 2, 4, 8].filter(
      (value, index, all) => value <= Math.max(maxMrd, 2) || all[index - 1] < maxMrd,
    );
    legend.replaceChildren(
      el('span.legend__item', {}, [el('span', { text: 'm.r.d.' })]),
      ...stops.map((value) =>
        el('span.legend__item', {}, [
          el('span.legend__swatch', { style: `background:${rampColor(value)}` }),
          el('span', { text: value === 1 ? '1 (random)' : String(value) }),
        ]),
      ),
    );
  }

  renderControls();
  context.stage.append(frame.element, legend, details);
  if (examples.length) loadExample(examples[0]);

  return { help: () => state.operation };
}

/* -------------------------------------------------------------- rendering */

function discFrame(root, axes) {
  root.append(
    svg('circle', {
      cx: 0,
      cy: 0,
      r: VIEW,
      fill: 'none',
      stroke: 'currentColor',
      'stroke-opacity': 0.6,
      'stroke-width': 1,
    }),
  );
  // Specimen axes, named. A pole figure without them cannot be read: the whole
  // content of "brass" versus "copper" is where the poles sit relative to RD.
  const positions = [
    [axes[0], VIEW * 1.08, 2],
    [axes[1], 0, -VIEW * 1.06],
  ];
  for (const [label, x, y] of positions) {
    root.append(
      svg('text', {
        x,
        y,
        'text-anchor': 'middle',
        'font-size': 7,
        fill: 'currentColor',
        'fill-opacity': 0.65,
        text: label,
      }),
    );
  }
}

/**
 * The pole figure, drawn as a filled density.
 *
 * The grid is equispaced on the sphere, not a raster, so there is no rectangular
 * cell to fill. Each sample is drawn as a disc whose radius is a little over
 * half the grid spacing on the projection, which overlaps its neighbours and
 * leaves no gaps — a Voronoi tessellation would be exact and is not worth the
 * code for a figure whose resolution is set by the kernel anyway.
 */
function renderDensity(data, frame) {
  const root = svg('svg', {
    viewBox: `${-VIEW * 1.16} ${-VIEW * 1.16} ${2.32 * VIEW} ${2.32 * VIEW}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': `${data.pole_label} pole figure`,
  });

  const clipId = `pf-clip-${Math.random().toString(36).slice(2, 9)}`;
  const defs = svg('defs');
  const clip = svg('clipPath', { id: clipId });
  clip.append(svg('circle', { cx: 0, cy: 0, r: VIEW }));
  defs.append(clip);
  root.append(defs);

  const body = svg('g', { 'clip-path': `url(#${clipId})` });
  // A radius derived from the point count rather than a constant: the control
  // that sets grid resolution would otherwise leave gaps at coarse settings and
  // a mush of overlap at fine ones.
  const spacing = (2 * VIEW) / Math.sqrt(Math.max(data.points.length, 1));
  const radius = spacing * 0.78;
  const columns = data.columns;
  for (const point of data.points) {
    const node = svg('circle', {
      cx: point.x * VIEW,
      cy: -point.y * VIEW,
      r: radius,
      fill: rampColor(point.mrd),
    });
    body.append(node);
    frame.hoverable(node, point, columns);
  }
  root.append(body);
  discFrame(root, data.specimen_axes ?? ['RD', 'TD']);
  return root;
}

/** The inverse pole figure: one marker per grain, plus the sector outline. */
function renderScatter(data, frame) {
  const root = svg('svg', {
    viewBox: `${-VIEW * 0.2} ${-VIEW * 0.9} ${VIEW * 1.1} ${VIEW * 1.1}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': `Inverse pole figure along ${data.axis_label}`,
  });

  const vertices = data.sector_vertices ?? [];
  if (vertices.length >= 3) {
    root.append(
      svg('polygon', {
        points: vertices.map(([x, y]) => `${x * VIEW},${-y * VIEW}`).join(' '),
        fill: 'none',
        stroke: 'currentColor',
        'stroke-opacity': 0.65,
        'stroke-width': 1,
      }),
    );
    // The corners are what an IPF is read against, so they are named.
    const corners = ['[001]', '[101]', '[111]'];
    vertices.forEach(([x, y], index) => {
      root.append(
        svg('text', {
          x: x * VIEW + (index === 0 ? -4 : 5),
          y: -y * VIEW + (index === 2 ? -4 : 8),
          'text-anchor': index === 0 ? 'end' : 'start',
          'font-size': 6,
          fill: 'currentColor',
          'fill-opacity': 0.65,
          text: corners[index] ?? '',
        }),
      );
    });
  }

  const columns = data.columns;
  for (const point of data.points) {
    const node = svg('circle', {
      cx: point.x * VIEW,
      cy: -point.y * VIEW,
      r: 1.3,
      fill: 'var(--accent)',
      'fill-opacity': 0.45,
    });
    root.append(node);
    frame.hoverable(node, point, columns);
  }
  return root;
}

/**
 * The three phi-2 sections, side by side.
 *
 * Drawn as three small maps rather than one large one because they are read
 * together — the argument about which component dominates a rolling texture is
 * made by comparing the 45 and 65 degree sections, and putting them on separate
 * screens defeats it.
 */
function renderSections(data, frame) {
  const sections = data.sections ?? [];
  const gap = 14;
  const size = 100;
  const width = sections.length * size + (sections.length - 1) * gap;
  const root = svg('svg', {
    viewBox: `-6 -16 ${width + 12} ${size + 32}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': 'ODF sections',
  });

  sections.forEach((section, index) => {
    const originX = index * (size + gap);
    const phi1 = section.phi1_deg;
    const bigPhi = section.big_phi_deg;
    const cellWidth = size / Math.max(phi1.length, 1);
    const cellHeight = size / Math.max(bigPhi.length, 1);

    root.append(
      svg('text', {
        x: originX + size / 2,
        y: -5,
        'text-anchor': 'middle',
        'font-size': 6,
        fill: 'currentColor',
        text: `φ₂ = ${section.phi2_deg}°`,
      }),
    );

    for (let row = 0; row < phi1.length; row += 1) {
      for (let column = 0; column < bigPhi.length; column += 1) {
        const mrd = section.densities[row][column];
        const node = svg('rect', {
          // phi-1 runs across, Phi down: the standard layout of a Bunge
          // section, and the one every published figure uses.
          x: originX + row * cellWidth,
          y: column * cellHeight,
          width: cellWidth + 0.4,
          height: cellHeight + 0.4,
          fill: rampColor(mrd),
        });
        root.append(node);
        frame.hoverable(node, {
          phi2_deg: section.phi2_deg,
          phi1_deg: phi1[row],
          big_phi_deg: bigPhi[column],
          mrd,
        }, data.columns);
      }
    }

    root.append(
      svg('rect', {
        x: originX,
        y: 0,
        width: size,
        height: size,
        fill: 'none',
        stroke: 'currentColor',
        'stroke-opacity': 0.5,
        'stroke-width': 0.6,
      }),
      svg('text', {
        x: originX + size / 2,
        y: size + 9,
        'text-anchor': 'middle',
        'font-size': 5,
        fill: 'currentColor',
        'fill-opacity': 0.6,
        text: 'φ₁ →',
      }),
    );
  });

  root.append(
    svg('text', {
      x: -3,
      y: size / 2,
      'text-anchor': 'middle',
      'font-size': 5,
      fill: 'currentColor',
      'fill-opacity': 0.6,
      transform: `rotate(-90 -3 ${size / 2})`,
      text: 'Φ →',
    }),
  );
  return root;
}
