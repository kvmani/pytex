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
import { download, renderResult } from '../core/result.js';
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

const DEFAULT_CONTOUR_STYLE = Object.freeze({
  mode: 'filled_lines',
  levelCount: 8,
  customLevels: '',
  scaleMax: 0,
  palette: 'mrd',
  lineColor: '#172033',
  lineWidth: 0.75,
  fillOpacity: 1,
  gridSize: 81,
});

function defaultContourStyle() {
  return { ...DEFAULT_CONTOUR_STYLE };
}

/**
 * The m.r.d. colour ramp.
 *
 * Sequential and perceptually ordered, running cool through warm, so that
 * "more" always reads as "warmer" without a legend. It is anchored so that
 * exactly 1 m.r.d. — the random baseline — sits at the pale middle of the ramp:
 * below 1 is cool, above 1 is warm, and an untextured material comes out a flat
 * neutral rather than a dramatic colour that happens to be its own minimum.
 */
const MRD_RAMP = [
  [0.0, [49, 84, 140]],
  [0.5, [90, 148, 194]],
  [1.0, [232, 234, 228]],
  [2.0, [233, 178, 84]],
  [4.0, [214, 106, 43]],
  [8.0, [154, 32, 32]],
];

const NORMALIZED_RAMPS = {
  viridis: [
    [0, [68, 1, 84]], [0.25, [59, 82, 139]], [0.5, [33, 145, 140]],
    [0.75, [94, 201, 98]], [1, [253, 231, 37]],
  ],
  turbo: [
    [0, [48, 18, 59]], [0.2, [50, 103, 188]], [0.4, [32, 204, 188]],
    [0.6, [164, 252, 60]], [0.8, [245, 125, 21]], [1, [122, 4, 3]],
  ],
};

function interpolateRamp(value, ramp) {
  for (let index = 1; index < ramp.length; index += 1) {
    const [lowStop, low] = ramp[index - 1];
    const [highStop, high] = ramp[index];
    if (value <= highStop || index === ramp.length - 1) {
      const span = highStop - lowStop || 1;
      const t = Math.min(Math.max((value - lowStop) / span, 0), 1);
      const mix = low.map((channel, c) => Math.round(channel + (high[c] - channel) * t));
      return `rgb(${mix[0]} ${mix[1]} ${mix[2]})`;
    }
  }
  return 'rgb(154 32 32)';
}

function displayMaximum(maxMrd, style) {
  return style.scaleMax > 0 ? style.scaleMax : Math.max(maxMrd, 1);
}

function paletteColor(mrd, style, maxMrd) {
  const value = Math.max(Number(mrd) || 0, 0);
  if (style.palette === 'mrd') return interpolateRamp(value, MRD_RAMP);
  const normalized = Math.min(value / displayMaximum(maxMrd, style), 1);
  return interpolateRamp(normalized, NORMALIZED_RAMPS[style.palette] ?? NORMALIZED_RAMPS.viridis);
}

function customContourLevels(text) {
  const values = String(text ?? '')
    .trim()
    .split(/[,;\s]+/)
    .filter(Boolean)
    .map(Number);
  if (!values.length || values.some((value) => !Number.isFinite(value) || value <= 0)) return null;
  const unique = [...new Set(values)].sort((left, right) => left - right);
  return unique.length >= 2 ? unique : null;
}

function contourLevels(maxMrd, style) {
  const custom = customContourLevels(style.customLevels);
  if (custom) return custom;
  const maximum = displayMaximum(maxMrd, style);
  const count = Math.max(2, Math.round(style.levelCount));
  const levels = Array.from({ length: count }, (_, index) => maximum * (index + 1) / (count + 1));
  if (maximum > 1 && !levels.some((level) => Math.abs(level - 1) < maximum / (count * 3))) {
    levels[Math.max(0, Math.min(levels.length - 1, Math.round(count / maximum) - 1))] = 1;
    levels.sort((left, right) => left - right);
  }
  return [...new Set(levels)];
}

function bandValue(value, levels) {
  const index = levels.findIndex((level) => value < level);
  if (index < 0) return levels.at(-1) * 1.08;
  const low = index === 0 ? 0 : levels[index - 1];
  return (low + levels[index]) / 2;
}

function contourStyleControl(style, { disabled, onChange, onReset }) {
  const fieldset = el('fieldset.contour-controls', { disabled });
  const mode = el('select', {
    oninput: (event) => {
      style.mode = event.currentTarget.value;
      onChange();
    },
  }, [
    ['filled_lines', 'Filled + lines'],
    ['filled', 'Filled contours'],
    ['lines', 'Contour lines'],
  ].map(([value, text]) => el('option', { value, text, selected: style.mode === value })));
  const palette = el('select', {
    oninput: (event) => {
      style.palette = event.currentTarget.value;
      onChange();
    },
  }, [
    ['mrd', 'm.r.d. baseline'],
    ['viridis', 'Viridis'],
    ['turbo', 'Turbo'],
  ].map(([value, text]) => el('option', { value, text, selected: style.palette === value })));
  const levelOutput = el('output', { text: String(style.levelCount) });
  const customHint = el('span.field__hint', {
    text: style.customLevels && !customContourLevels(style.customLevels)
      ? 'Enter at least two positive levels, for example 0.5, 1, 2, 4.'
      : 'Optional exact m.r.d. isolines; when set, these replace the automatic count.',
  });
  const lineOutput = el('output', { text: `${style.lineWidth.toFixed(2)} px` });
  const opacityOutput = el('output', { text: style.fillOpacity.toFixed(2) });
  const lineColorOutput = el('output', { text: style.lineColor.toUpperCase() });
  fieldset.append(
    el('label.field', {}, [el('span.field__label', { text: 'Contour display' }), mode]),
    el('label.field', {}, [
      el('span.field__label', { text: 'Automatic levels' }),
      el('span.range-control', {}, [
        el('input', {
          type: 'range', min: 2, max: 20, step: 1, value: style.levelCount,
          oninput: (event) => {
            style.levelCount = Number(event.currentTarget.value);
            levelOutput.textContent = String(style.levelCount);
            onChange();
          },
        }),
        levelOutput,
      ]),
      el('span.field__hint', { text: 'Number of isolines when custom levels are empty.' }),
    ]),
    el('label.field', {}, [
      el('span.field__label', { text: 'Custom levels' }),
      el('input', {
        type: 'text', value: style.customLevels, placeholder: '0.5, 1, 2, 4',
        oninput: (event) => {
          style.customLevels = event.currentTarget.value;
          customHint.textContent = style.customLevels && !customContourLevels(style.customLevels)
            ? 'Enter at least two positive levels, for example 0.5, 1, 2, 4.'
            : 'Optional exact m.r.d. isolines; when set, these replace the automatic count.';
          onChange();
        },
      }),
      customHint,
    ]),
    el('label.field', {}, [
      el('span.field__label', { text: 'Upper colour limit' }),
      el('input', {
        type: 'number', min: 0, step: 0.1, value: style.scaleMax,
        oninput: (event) => {
          style.scaleMax = Math.max(Number(event.currentTarget.value) || 0, 0);
          onChange();
        },
      }),
      el('span.field__hint', { text: '0 uses the data peak; a positive value clips the colour scale.' }),
    ]),
    el('label.field', {}, [el('span.field__label', { text: 'Colour palette' }), palette]),
    el('label.field', {}, [
      el('span.field__label', { text: 'Line colour' }),
      el('span.color-control', {}, [
        el('input', {
          type: 'color', value: style.lineColor,
          oninput: (event) => {
            style.lineColor = event.currentTarget.value;
            lineColorOutput.textContent = style.lineColor.toUpperCase();
            onChange();
          },
        }),
        lineColorOutput,
      ]),
    ]),
    el('label.field', {}, [
      el('span.field__label', { text: 'Line width' }),
      el('span.range-control', {}, [
        el('input', {
          type: 'range', min: 0.25, max: 2.5, step: 0.05, value: style.lineWidth,
          oninput: (event) => {
            style.lineWidth = Number(event.currentTarget.value);
            lineOutput.textContent = `${style.lineWidth.toFixed(2)} px`;
            onChange();
          },
        }), lineOutput,
      ]),
    ]),
    el('label.field', {}, [
      el('span.field__label', { text: 'Fill opacity' }),
      el('span.range-control', {}, [
        el('input', {
          type: 'range', min: 0.1, max: 1, step: 0.05, value: style.fillOpacity,
          oninput: (event) => {
            style.fillOpacity = Number(event.currentTarget.value);
            opacityOutput.textContent = style.fillOpacity.toFixed(2);
            onChange();
          },
        }), opacityOutput,
      ]),
    ]),
    el('label.field', {}, [
      el('span.field__label', { text: 'Display grid' }),
      el('select', {
        oninput: (event) => {
          style.gridSize = Number(event.currentTarget.value);
          onChange({ rebuildGrid: true });
        },
      }, [49, 65, 81, 97, 129].map((value) =>
        el('option', { value, text: `${value} × ${value}`, selected: style.gridSize === value }),
      )),
      el('span.field__hint', { text: 'Display interpolation only; source values and exports stay unchanged.' }),
    ]),
    el('button.button', { type: 'button', text: 'Reset contour properties', onclick: onReset }),
  );
  return el('details.group.appearance', {}, [
    el('summary', { text: 'Contour properties' }),
    el('div.group__body', {}, [
      el('p.field__help', {
        text: disabled
          ? 'Contour properties apply to pole figures and ODF sections; the inverse pole figure is one point per grain.'
          : 'Presentation only. Levels, palette and interpolation do not alter the ODF or exported m.r.d. rows.',
      }),
      fieldset,
    ]),
  ]);
}

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
    contour: defaultContourStyle(),
    contourGrid: null,
    plotNode: null,
  };

  const frame = plotFrame({
    title: 'Pole figure',
    toolbar: [
      el('button.button', {
        type: 'button',
        text: 'SVG',
        title: 'Save the current interactive texture figure as SVG',
        onclick: () => {
          if (!state.plotNode) return;
          const markup = new XMLSerializer().serializeToString(state.plotNode);
          download('pytex-texture-figure.svg', markup, 'image/svg+xml');
        },
      }),
    ],
  });
  const legend = el('div.legend');
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
        state.contourGrid = null;
        renderControls(carryOver());
        renderAppearanceControls();
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
    appearanceHost,
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

  function renderAppearanceControls() {
    const supportsContours = state.operation.id !== 'texture.inverse_pole_figure';
    appearanceHost.replaceChildren(
      contourStyleControl(state.contour, {
        disabled: !supportsContours,
        onChange: ({ rebuildGrid = false } = {}) => {
          if (rebuildGrid) state.contourGrid = null;
          if (state.result) draw(true);
        },
        onReset: () => {
          state.contour = defaultContourStyle();
          state.contourGrid = null;
          renderAppearanceControls();
          if (state.result) draw(true);
        },
      }),
    );
  }

  function loadExample(example) {
    state.operation = operations.find((entry) => entry.id === example.operation);
    viewSelect.value = state.operation.id;
    state.teaches = example.teaches;
    state.contourGrid = null;
    renderControls(example.request);
    renderAppearanceControls();
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Building…';
    state.form.clearErrors();
    try {
      const result = await call(state.operation.id, state.form.values());
      state.result = result;
      state.contourGrid = null;
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

  function draw(preserveViewport = false) {
    const data = state.result.data;
    const kind = state.operation.id;
    frame.element.hidden = false;
    legend.hidden = false;

    if (kind === 'texture.odf_sections') {
      frame.configure({ toData: null, formatCursor: null });
      state.plotNode = renderSections(data, frame, state.contour);
      frame.setContent(state.plotNode, { preserveViewport });
      renderRampLegend(data.max_mrd, state.contour);
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
      state.contourGrid ??= interpolatePoleFigure(data.points, state.contour.gridSize);
      state.plotNode = renderDensity(data, frame, state.contour, state.contourGrid);
      frame.setContent(state.plotNode, { preserveViewport });
      renderRampLegend(data.max_mrd, state.contour);
      frame.setStatus(
        `${data.pole_label} · peak ${formatNumber(data.max_mrd, 2)} m.r.d. · ` +
          `mean ${formatNumber(data.mean_mrd, 3)} m.r.d. (1 by construction) · ` +
          `${data.grain_count} grains · hover for the intensity at a point`,
      );
    } else {
      state.plotNode = renderScatter(data, frame);
      frame.setContent(state.plotNode, { preserveViewport });
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
  function renderRampLegend(maxMrd, style) {
    const stops = contourLevels(maxMrd, style);
    legend.replaceChildren(
      el('span.legend__item', {}, [el('span', { text: 'm.r.d.' })]),
      ...stops.map((value) =>
        el('span.legend__item', {}, [
          el('span.legend__swatch', { style: `background:${paletteColor(value, style, maxMrd)}` }),
          el('span', { text: value === 1 ? '1 (random)' : String(value) }),
        ]),
      ),
    );
  }

  renderControls();
  renderAppearanceControls();
  // The legend is a control, so it rides inside the frame rather than under it:
  // toggling a source and seeing the drawing change must not need a scroll.
  frame.setControls(legend);
  context.stage.append(frame.element, details);
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

/** Display-only inverse-distance interpolation onto the projection disc. */
function interpolatePoleFigure(points, size) {
  const xValues = Array.from({ length: size }, (_, index) => -1 + 2 * index / (size - 1));
  const yValues = [...xValues];
  const neighbours = 6;
  const values = yValues.map((y) => xValues.map((x) => {
    if (x * x + y * y > 1.0001) return Number.NaN;
    const nearest = [];
    for (const point of points) {
      const distance2 = (point.x - x) ** 2 + (point.y - y) ** 2;
      if (distance2 < 1e-12) return point.mrd;
      if (nearest.length < neighbours || distance2 < nearest.at(-1)[0]) {
        nearest.push([distance2, point.mrd]);
        nearest.sort((left, right) => left[0] - right[0]);
        if (nearest.length > neighbours) nearest.pop();
      }
    }
    let weighted = 0;
    let total = 0;
    for (const [distance2, value] of nearest) {
      const weight = 1 / (distance2 * distance2 + 1e-12);
      weighted += weight * value;
      total += weight;
    }
    return total ? weighted / total : Number.NaN;
  }));
  return { xValues, yValues, values };
}

function crossingPoint(a, b, level) {
  const span = b.value - a.value;
  const t = Math.abs(span) < 1e-14 ? 0.5 : Math.min(Math.max((level - a.value) / span, 0), 1);
  return { x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t };
}

/** One SVG path containing every marching-squares segment for one level. */
function contourPath(grid, level, mapPoint) {
  const commands = [];
  for (let row = 0; row < grid.yValues.length - 1; row += 1) {
    for (let column = 0; column < grid.xValues.length - 1; column += 1) {
      const corners = [
        { x: grid.xValues[column], y: grid.yValues[row], value: grid.values[row][column] },
        { x: grid.xValues[column + 1], y: grid.yValues[row], value: grid.values[row][column + 1] },
        { x: grid.xValues[column + 1], y: grid.yValues[row + 1], value: grid.values[row + 1][column + 1] },
        { x: grid.xValues[column], y: grid.yValues[row + 1], value: grid.values[row + 1][column] },
      ];
      if (corners.some((corner) => !Number.isFinite(corner.value))) continue;
      const crossings = [];
      [[0, 1], [1, 2], [2, 3], [3, 0]].forEach(([start, end], edge) => {
        const a = corners[start];
        const b = corners[end];
        if ((a.value < level) !== (b.value < level)) {
          crossings.push({ edge, point: crossingPoint(a, b, level) });
        }
      });
      let pairs = [];
      if (crossings.length === 2) {
        pairs = [[crossings[0], crossings[1]]];
      } else if (crossings.length === 4) {
        const center = corners.reduce((sum, corner) => sum + corner.value, 0) / 4;
        const byEdge = Object.fromEntries(crossings.map((entry) => [entry.edge, entry]));
        pairs = center >= level
          ? [[byEdge[0], byEdge[1]], [byEdge[2], byEdge[3]]]
          : [[byEdge[0], byEdge[3]], [byEdge[1], byEdge[2]]];
      }
      for (const [start, end] of pairs) {
        const a = mapPoint(start.point.x, start.point.y);
        const b = mapPoint(end.point.x, end.point.y);
        commands.push(`M ${a.x} ${a.y} L ${b.x} ${b.y}`);
      }
    }
  }
  return commands.join(' ');
}

function drawContourGrid(parent, grid, style, maxMrd, mapPoint) {
  const levels = contourLevels(maxMrd, style);
  if (style.mode !== 'lines') {
    const fill = svg('g', { 'fill-opacity': style.fillOpacity });
    for (let row = 0; row < grid.yValues.length - 1; row += 1) {
      for (let column = 0; column < grid.xValues.length - 1; column += 1) {
        const samples = [
          grid.values[row][column], grid.values[row][column + 1],
          grid.values[row + 1][column], grid.values[row + 1][column + 1],
        ];
        const finite = samples.filter(Number.isFinite);
        if (!finite.length) continue;
        const value = finite.reduce((sum, sample) => sum + sample, 0) / finite.length;
        const a = mapPoint(grid.xValues[column], grid.yValues[row]);
        const b = mapPoint(grid.xValues[column + 1], grid.yValues[row + 1]);
        fill.append(svg('rect', {
          x: Math.min(a.x, b.x),
          y: Math.min(a.y, b.y),
          width: Math.abs(b.x - a.x) + 0.35,
          height: Math.abs(b.y - a.y) + 0.35,
          fill: paletteColor(bandValue(value, levels), style, maxMrd),
        }));
      }
    }
    parent.append(fill);
  }
  if (style.mode !== 'filled') {
    const lines = svg('g', { 'aria-label': `Contour lines at ${levels.join(', ')} m.r.d.` });
    for (const level of levels) {
      const d = contourPath(grid, level, mapPoint);
      if (!d) continue;
      lines.append(svg('path', {
        d,
        fill: 'none',
        stroke: style.lineColor,
        'stroke-width': style.lineWidth,
        'stroke-linejoin': 'round',
        'vector-effect': 'non-scaling-stroke',
        'data-contour-level': level,
      }));
    }
    parent.append(lines);
  }
  return levels;
}

/**
 * The pole figure, drawn as a filled density.
 *
 * The computed support is equispaced on the sphere rather than on screen. For
 * display only it is inverse-distance interpolated onto a clipped projection
 * grid; marching squares and band quantisation consume that same grid so line
 * and filled contours share exactly the declared levels. Transparent hit
 * regions remain on the computed samples, not the interpolated cells.
 */
function renderDensity(data, frame, style, grid) {
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
  drawContourGrid(body, grid, style, data.max_mrd, (x, y) => ({ x: x * VIEW, y: -y * VIEW }));
  const columns = data.columns;
  const hitRadius = (2 * VIEW) / Math.sqrt(Math.max(data.points.length, 1));
  for (const point of data.points) {
    const node = svg('circle', {
      cx: point.x * VIEW,
      cy: -point.y * VIEW,
      r: hitRadius,
      fill: 'transparent',
      'pointer-events': 'all',
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
function renderSections(data, frame, style) {
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
    const grid = {
      xValues: phi1.map((_, point) => point / Math.max(phi1.length - 1, 1)),
      yValues: bigPhi.map((_, point) => point / Math.max(bigPhi.length - 1, 1)),
      values: bigPhi.map((_, column) => phi1.map((__, row) => section.densities[row][column])),
    };

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

    drawContourGrid(
      root,
      grid,
      style,
      data.max_mrd,
      (x, y) => ({ x: originX + x * size, y: y * size }),
    );

    for (let row = 0; row < phi1.length; row += 1) {
      for (let column = 0; column < bigPhi.length; column += 1) {
        const mrd = section.densities[row][column];
        const node = svg('rect', {
          x: originX + row * cellWidth - cellWidth / 2,
          y: column * cellHeight - cellHeight / 2,
          width: cellWidth,
          height: cellHeight,
          fill: 'transparent',
          'pointer-events': 'all',
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
