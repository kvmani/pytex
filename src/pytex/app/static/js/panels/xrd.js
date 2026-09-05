/** Powder XRD: an indexed diffractogram whose peaks remain inspectable. */

import { call } from '../core/api.js';
import { explainer } from '../core/explainer.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'xrd',
  title: 'XRD',
  tagline: 'Structure-aware powder diffraction with indexed, inspectable peaks.',
};

const WIDTH = 1000;
const HEIGHT = 520;
const MARGIN = { left: 76, right: 24, top: 30, bottom: 68 };
const DEFAULT_APPEARANCE = Object.freeze({
  lineColor: '#2563eb',
  stickColor: '#7c3aed',
  lineWidth: 2,
  fill: true,
  sticks: true,
  labels: true,
  labelThreshold: 0.08,
  yScale: 'linear',
});

/**
 * The four questions this workspace answers, in the order they depend on each
 * other: what would this phase diffract, what is the background of a real
 * scan, how much of a peak's width is the instrument, and what does the scan
 * say the structure actually is.
 */
const VIEWS = [
  'xrd.powder_pattern',
  'xrd.background',
  'xrd.rietveld',
  'xrd.size_strain',
];

/** Which kind of drawing each view needs. */
const VIEW_MODES = {
  'xrd.powder_pattern': 'profile',
  'xrd.background': 'overlay',
  'xrd.rietveld': 'overlay',
  'xrd.size_strain': 'scatter',
};

/** What the run button says while each view is working. */
const VIEW_ACTIONS = {
  'xrd.powder_pattern': ['Simulate XRD pattern', 'Simulating…'],
  'xrd.background': ['Estimate background', 'Estimating…'],
  'xrd.rietveld': ['Refine against the scan', 'Refining…'],
  'xrd.size_strain': ['Separate size and strain', 'Fitting…'],
};

/** Curves drawn by the overlay views, in draw order, with their roles. */
const OVERLAY_SERIES = {
  'xrd.background': [
    { key: 'observed', label: 'Observed', color: '#94a3b8', width: 1 },
    { key: 'background', label: 'Estimated background', color: '#dc2626', width: 2 },
    { key: 'subtracted', label: 'Observed − background', color: '#2563eb', width: 1.4 },
  ],
  'xrd.rietveld': [
    { key: 'observed', label: 'Observed', color: '#dc2626', width: 1 },
    { key: 'calculated', label: 'Calculated', color: '#2563eb', width: 1.6 },
    { key: 'background', label: 'Refined background', color: '#15803d', width: 1.2 },
  ],
};

export function mount(context) {
  const operations = VIEWS
    .map((id) => context.manifest.operations.find((entry) => entry.id === id))
    .filter(Boolean);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = {
    operation: operations[0],
    result: null,
    form: null,
    teaches: null,
    appearance: { ...DEFAULT_APPEARANCE },
  };

  const frame = plotFrame({ title: 'Powder diffractogram' });
  const legend = el('div.legend');
  const details = el('div');
  const formHost = el('div');
  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: VIEW_ACTIONS[state.operation.id][0],
    onclick: () => run(),
  });

  function redraw() {
    if (state.result) draw();
  }

  // The appearance controls only govern the simulated-profile renderer. They
  // stay mounted so the panel keeps one control rail, and hide themselves for
  // the views whose drawing they cannot change -- a visible control that does
  // nothing is worse than an absent one.
  const appearance = appearanceControls(state.appearance, redraw);

  // The vertical scale is not an appearance choice in the way a colour is: on a
  // linear axis scaled to 20,000-count peaks a 150-count background is a line on
  // the axis, so the background view cannot show the thing it exists to show.
  // It therefore lives outside the appearance group and stays available to every
  // view that draws a profile.
  const scaleControl = el('label.field', {}, [
    el('span.field__label', { text: 'Vertical display scale' }),
    el('select', {
      oninput: (event) => {
        state.appearance.yScale = event.currentTarget.value;
        redraw();
      },
    }, [
      el('option', { value: 'linear', text: 'Linear' }),
      el('option', { value: 'sqrt', text: 'Square-root' }),
      el('option', { value: 'log', text: 'Log-like' }),
    ]),
    el('span.field__hint', {
      text: 'Square-root and log-like views reveal the background and the weak peaks without '
        + 'changing a single number.',
    }),
  ]);

  const viewSelect = el('select', {
    oninput: (event) => selectView(event.currentTarget.value),
  }, operations.map((entry) => el('option', {
    value: entry.id,
    text: entry.title,
    selected: entry.id === state.operation.id,
  })));

  function selectView(id) {
    const chosen = operations.find((entry) => entry.id === id);
    if (!chosen || chosen === state.operation) return;
    state.operation = chosen;
    state.result = null;
    state.teaches = null;
    frame.setTitle(chosen.title);
    frame.setStatus(chosen.summary);
    frame.setContent(null);
    legend.replaceChildren();
    details.replaceChildren();
    appearance.hidden = mode() !== 'profile';
    scaleControl.hidden = mode() === 'scatter';
    // A run left in flight by the previous view will decline to touch the
    // button, so the new view has to hand it back itself.
    runButton.disabled = false;
    runButton.textContent = VIEW_ACTIONS[chosen.id][0];
    renderControls();
  }

  function mode() {
    return VIEW_MODES[state.operation.id] ?? 'profile';
  }

  context.rail.append(
    el('label.field', {}, [
      el('span.field__label', { text: 'View' }),
      viewSelect,
      el('span.field__hint', {
        text: 'Simulation first, then the three questions a measured scan raises.',
      }),
    ]),
    formHost,
    runButton,
    scaleControl,
    appearance,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        explainer(
          'These standards isolate extinction rules, wavelength, doublet splitting and hexagonal metrics.',
          { label: 'What these examples show' },
        ),
        el('div.examples', {}, examples.map((example) =>
          el('button.example', { type: 'button', onclick: () => loadExample(example) }, [
            el('strong', { text: example.title }),
            el('span', { text: example.summary }),
          ]),
        )),
      ]),
    ]),
  );
  // The legend is a control, so it rides inside the frame rather than under it:
  // toggling a source and seeing the drawing change must not need a scroll.
  frame.setControls(legend);
  context.stage.append(frame.element, details);

  function renderControls(initial = {}) {
    state.form = buildForm(state.operation, { initial });
    formHost.replaceChildren(state.form.element);
  }

  function loadExample(example) {
    // An example names the operation it exercises, so choosing one switches the
    // view. Without this a background example loaded into the simulation form
    // would silently drop every parameter the simulation does not have.
    const target = operations.find((entry) => entry.id === example.operation);
    if (target && target !== state.operation) {
      state.operation = target;
      viewSelect.value = target.id;
      frame.setTitle(target.title);
      appearance.hidden = mode() !== 'profile';
      scaleControl.hidden = mode() === 'scatter';
      runButton.disabled = false;
      runButton.textContent = VIEW_ACTIONS[target.id][0];
    }
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    // Each run remembers the view it was launched for. A request in flight when
    // the view changes would otherwise resolve into the new view and be drawn by
    // a renderer expecting different keys -- a refinement's profile handed to the
    // background drawing, which fails as a type error rather than as a wrong
    // picture. The launching view is the one allowed to consume the answer.
    const launched = state.operation;
    const [idle, busy] = VIEW_ACTIONS[launched.id];
    runButton.disabled = true;
    runButton.textContent = busy;
    state.form.clearErrors();
    try {
      const result = await call(launched.id, state.form.values());
      if (state.operation !== launched) return;
      state.result = result;
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
    } catch (error) {
      if (state.operation !== launched) return;
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      if (state.operation === launched) {
        runButton.disabled = false;
        runButton.textContent = idle;
      }
    }
  }

  function draw() {
    if (mode() === 'overlay') return drawOverlay();
    if (mode() === 'scatter') return drawScatter();
    return drawProfile();
  }

  function drawOverlay() {
    const data = state.result.data;
    const series = OVERLAY_SERIES[state.operation.id];
    const minimum = data.two_theta_deg[0];
    const maximum = data.two_theta_deg[data.two_theta_deg.length - 1];
    let ceiling = 0;
    for (const entry of series) {
      for (const value of data[entry.key]) if (value > ceiling) ceiling = value;
    }
    ceiling = ceiling > 0 ? ceiling : 1;
    frame.configure({
      toData: (x, y) => {
        const plotWidth = WIDTH - MARGIN.left - MARGIN.right;
        const plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom;
        if (
          x < MARGIN.left || x > WIDTH - MARGIN.right ||
          y < MARGIN.top || y > HEIGHT - MARGIN.bottom
        ) return null;
        return {
          x: minimum + ((x - MARGIN.left) / plotWidth) * (maximum - minimum),
          y: (1 - (y - MARGIN.top) / plotHeight) * ceiling,
        };
      },
      formatCursor: (point) =>
        `${formatNumber(point.x, 3)}° 2θ · ${formatNumber(point.y, 0)} counts`,
    });
    frame.setContent(
      renderOverlay(data, series, ceiling, state.operation.id, state.appearance.yScale),
    );
    legend.replaceChildren(...series.map((entry) => el('span.legend__item', {}, [
      el('span.legend__swatch', { style: `background:${entry.color}` }),
      el('span', { text: entry.label }),
    ])), ...(state.operation.id === 'xrd.rietveld' ? [el('span.legend__item', {}, [
      el('span.legend__swatch', { style: 'background:#7c3aed' }),
      el('span', { text: 'Observed − calculated' }),
    ])] : []));
    frame.setStatus(overlayStatus(data, state.operation.id));
  }

  function drawScatter() {
    const data = state.result.data;
    frame.configure({
      toData: () => null,
      formatCursor: () => '',
    });
    frame.setContent(renderWilliamsonHall(data));
    legend.replaceChildren(
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: 'background:#2563eb' }),
        el('span', { text: 'Reflections' }),
      ]),
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: 'background:#dc2626' }),
        el('span', { text: 'Fitted line' }),
      ]),
    );
    frame.setStatus(
      `D = ${formatNumber(data.crystallite_size_nm, 3)} nm · ` +
      `ε = ${formatNumber(data.microstrain, 5)} · R² = ${formatNumber(data.r_squared, 5)} · ` +
      'the intercept is the size and the slope is the strain',
    );
  }

  function drawProfile() {
    const data = state.result.data;
    const minimum = data.two_theta_deg[0];
    const maximum = data.two_theta_deg[data.two_theta_deg.length - 1];
    frame.configure({
      toData: (x, y) => {
        const plotWidth = WIDTH - MARGIN.left - MARGIN.right;
        const plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom;
        if (
          x < MARGIN.left || x > WIDTH - MARGIN.right ||
          y < MARGIN.top || y > HEIGHT - MARGIN.bottom
        ) return null;
        return {
          x: minimum + ((x - MARGIN.left) / plotWidth) * (maximum - minimum),
          y: 1 - (y - MARGIN.top) / plotHeight,
        };
      },
      formatCursor: (point) =>
        `${formatNumber(point.x, 3)}° 2θ · displayed height ${formatNumber(point.y, 3)}`,
    });
    frame.setContent(renderPattern(data, state.appearance, frame));
    legend.replaceChildren(
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: `background:${state.appearance.lineColor}` }),
        el('span', { text: `${data.phase_name} profile` }),
      ]),
      ...(state.appearance.sticks ? [el('span.legend__item', {}, [
        el('span.legend__swatch', { style: `background:${state.appearance.stickColor}` }),
        el('span', { text: 'Kα1 reflection families' }),
      ])] : []),
    );
    frame.setStatus(
      `${data.reflections.length} indexed Kα1 families · ${data.radiation_name} · ` +
      `λ ${formatNumber(data.wavelength_angstrom, 6)} Å · ` +
      `${state.appearance.yScale} display · hover a peak for d, multiplicity and |F|`,
    );
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);
  // The help drawer must follow the view. Returning the operation captured at
  // mount would document the simulation while a refinement is on screen.
  return { help: () => state.operation };
}

function appearanceControls(style, onChange) {
  const root = el('details.group.appearance');
  const body = el('div.group__body');

  const colorControl = (label, key, help) => {
    const output = el('output', { text: style[key].toUpperCase() });
    const input = el('input', {
      type: 'color',
      value: style[key],
      oninput: (event) => {
        style[key] = event.currentTarget.value;
        output.textContent = style[key].toUpperCase();
        onChange();
      },
    });
    return el('label.field', {}, [
      el('span.field__label', { text: label }),
      el('span.color-control', {}, [input, output]),
      el('span.field__hint', { text: help }),
    ]);
  };

  const checkbox = (label, key, help) => el('label.field', {}, [
    el('span.checkbox', {}, [
      el('input', {
        type: 'checkbox',
        checked: style[key],
        onchange: (event) => {
          style[key] = event.currentTarget.checked;
          onChange();
        },
      }),
      el('span', { text: label }),
    ]),
    el('span.field__hint', { text: help }),
  ]);

  const lineOutput = el('output', { text: `${style.lineWidth.toFixed(1)} px` });
  const thresholdOutput = el('output', { text: `${Math.round(style.labelThreshold * 100)}%` });
  body.append(
    explainer(
      'Display controls redraw the existing profile. They never change peak positions, integrated intensities or exports.',
      { label: 'What these controls change' },
    ),
    colorControl('Profile colour', 'lineColor', 'Line and optional area-fill colour.'),
    colorControl('Reflection-stick colour', 'stickColor', 'Independent colour for indexed families.'),
    el('label.field', {}, [
      el('span.field__label', { text: 'Profile line width' }),
      el('span.range-control', {}, [
        el('input', {
          type: 'range', min: 0.5, max: 5, step: 0.25, value: style.lineWidth,
          oninput: (event) => {
            style.lineWidth = Number(event.currentTarget.value);
            lineOutput.textContent = `${style.lineWidth.toFixed(1)} px`;
            onChange();
          },
        }),
        lineOutput,
      ]),
    ]),
    checkbox('Fill below profile', 'fill', 'A translucent fill makes peak envelopes easier to compare.'),
    checkbox('Show reflection sticks', 'sticks', 'Draw the primary Kα1 families beneath the profile.'),
    checkbox('Label strong peaks', 'labels', 'Labels use canonical crystallographic notation.'),
    el('label.field', {}, [
      el('span.field__label', { text: 'Peak-label threshold' }),
      el('span.range-control', {}, [
        el('input', {
          type: 'range', min: 0, max: 1, step: 0.01, value: style.labelThreshold,
          oninput: (event) => {
            style.labelThreshold = Number(event.currentTarget.value);
            thresholdOutput.textContent = `${Math.round(style.labelThreshold * 100)}%`;
            onChange();
          },
        }),
        thresholdOutput,
      ]),
      el('span.field__hint', { text: 'Relative Kα1 integrated intensity required for a label.' }),
    ]),
    el('button.button', {
      type: 'button',
      text: 'Reset XRD appearance',
      onclick: () => {
        Object.assign(style, DEFAULT_APPEARANCE);
        root.replaceWith(appearanceControls(style, onChange));
        onChange();
      },
    }),
  );
  root.append(el('summary', { text: 'Appearance' }), body);
  return root;
}

function transformedIntensity(value, mode) {
  const safe = Math.max(0, Math.min(1, Number(value) || 0));
  if (mode === 'sqrt') return Math.sqrt(safe);
  if (mode === 'log') return Math.log10(1 + 99 * safe) / 2;
  return safe;
}

function renderPattern(data, appearance, frame) {
  const root = svg('svg', {
    viewBox: `0 0 ${WIDTH} ${HEIGHT}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': `Powder XRD pattern of ${data.phase_name}`,
  });
  const plotWidth = WIDTH - MARGIN.left - MARGIN.right;
  const plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom;
  const minimum = data.two_theta_deg[0];
  const maximum = data.two_theta_deg[data.two_theta_deg.length - 1];
  const xAt = (angle) => MARGIN.left + ((angle - minimum) / (maximum - minimum)) * plotWidth;
  const yAt = (intensity) =>
    MARGIN.top + (1 - transformedIntensity(intensity, appearance.yScale)) * plotHeight;

  for (let value = 0; value <= 1.0001; value += 0.25) {
    const y = yAt(value);
    root.append(
      svg('line', {
        x1: MARGIN.left, y1: y, x2: WIDTH - MARGIN.right, y2: y,
        stroke: 'currentColor', 'stroke-opacity': value === 0 ? 0.5 : 0.1,
        'stroke-width': value === 0 ? 1 : 0.6,
      }),
      svg('text', {
        x: MARGIN.left - 12, y: y + 4, 'text-anchor': 'end', 'font-size': 12,
        fill: 'currentColor', 'fill-opacity': 0.6, text: value.toFixed(2),
      }),
    );
  }
  const tickStep = niceStep(maximum - minimum, 8);
  const firstTick = Math.ceil(minimum / tickStep) * tickStep;
  for (let value = firstTick; value <= maximum + 1e-9; value += tickStep) {
    const x = xAt(value);
    root.append(
      svg('line', {
        x1: x, y1: MARGIN.top, x2: x, y2: HEIGHT - MARGIN.bottom,
        stroke: 'currentColor', 'stroke-opacity': 0.08, 'stroke-width': 0.6,
      }),
      svg('text', {
        x, y: HEIGHT - MARGIN.bottom + 24, 'text-anchor': 'middle', 'font-size': 12,
        fill: 'currentColor', 'fill-opacity': 0.65, text: formatNumber(value, 0),
      }),
    );
  }
  root.append(
    svg('text', {
      x: MARGIN.left + plotWidth / 2, y: HEIGHT - 16, 'text-anchor': 'middle',
      'font-size': 14, fill: 'currentColor', text: '2θ (°)',
    }),
    svg('text', {
      x: 19, y: MARGIN.top + plotHeight / 2, 'text-anchor': 'middle',
      'font-size': 14, fill: 'currentColor',
      transform: `rotate(-90 19 ${MARGIN.top + plotHeight / 2})`,
      text: 'Normalized intensity',
    }),
  );

  const points = data.two_theta_deg.map((angle, index) =>
    [xAt(angle), yAt(data.intensity[index])]);
  const linePath = points.map(([x, y], index) =>
    `${index === 0 ? 'M' : 'L'} ${x.toFixed(2)} ${y.toFixed(2)}`).join(' ');
  if (appearance.fill) {
    root.append(svg('path', {
      d: `${linePath} L ${WIDTH - MARGIN.right} ${HEIGHT - MARGIN.bottom} ` +
        `L ${MARGIN.left} ${HEIGHT - MARGIN.bottom} Z`,
      fill: appearance.lineColor, 'fill-opacity': 0.12, stroke: 'none',
    }));
  }
  root.append(svg('path', {
    d: linePath, fill: 'none', stroke: appearance.lineColor,
    'stroke-width': appearance.lineWidth, 'vector-effect': 'non-scaling-stroke',
  }));

  for (const [index, reflection] of data.reflections.entries()) {
    const x = xAt(reflection.two_theta_deg);
    const indexOnGrid = nearestIndex(data.two_theta_deg, reflection.two_theta_deg);
    const peakY = yAt(data.intensity[indexOnGrid]);
    if (appearance.sticks) {
      root.append(svg('line', {
        x1: x, y1: HEIGHT - MARGIN.bottom,
        x2: x, y2: HEIGHT - MARGIN.bottom - 52 * reflection.relative_intensity,
        stroke: appearance.stickColor, 'stroke-width': 1.2, 'stroke-opacity': 0.82,
      }));
    }
    if (appearance.labels && reflection.relative_intensity >= appearance.labelThreshold) {
      const labelY = Math.max(peakY - 9 - (index % 2) * 13, MARGIN.top + 11);
      root.append(svg('text', {
        x: x + 3, y: labelY, 'font-size': 11, fill: 'currentColor',
        'fill-opacity': 0.75, transform: `rotate(-45 ${x + 3} ${labelY})`,
        text: reflection.hkl_label,
      }));
    }
    const hit = svg('line', {
      x1: x, y1: MARGIN.top, x2: x, y2: HEIGHT - MARGIN.bottom,
      stroke: 'transparent', 'stroke-width': 10,
    });
    root.append(hit);
    frame.hoverable(hit, reflection, data.columns);
  }
  return root;
}

function overlayStatus(data, operationId) {
  const synthetic = data.synthetic
    ? 'generated demonstration scan (not a measurement) · '
    : '';
  if (operationId === 'xrd.background') {
    return (
      `${synthetic}${data.method} · ` +
      `${formatNumber(100 * data.background_fraction, 1)}% of the signal assigned to background`
    );
  }
  return (
    `${synthetic}R_wp ${formatNumber(100 * data.weighted_profile_r_factor, 3)}% · ` +
    `R_exp ${formatNumber(100 * data.expected_r_factor, 3)}% · ` +
    `GoF ${formatNumber(data.goodness_of_fit, 3)} · ` +
    `R_Bragg ${formatNumber(100 * data.bragg_r_factor, 3)}% · ` +
    `Durbin–Watson ${formatNumber(data.durbin_watson, 3)}` +
    (data.converged ? '' : ' · DID NOT CONVERGE')
  );
}

/**
 * Draw several curves against a shared count axis.
 *
 * The refinement view also gets a difference band beneath the profile, because
 * the residual is the most informative single output a refinement has and it is
 * unreadable at the scale of the peaks it came from.
 */
function renderOverlay(data, series, ceiling, operationId, yScale) {
  const withResidual = operationId === 'xrd.rietveld';
  const root = svg('svg', {
    viewBox: `0 0 ${WIDTH} ${HEIGHT}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': withResidual ? 'Rietveld refinement' : 'Estimated background',
  });
  const residualHeight = withResidual ? 108 : 0;
  const plotWidth = WIDTH - MARGIN.left - MARGIN.right;
  const plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom - residualHeight;
  const minimum = data.two_theta_deg[0];
  const maximum = data.two_theta_deg[data.two_theta_deg.length - 1];
  const xAt = (angle) => MARGIN.left + ((angle - minimum) / (maximum - minimum)) * plotWidth;
  const yAt = (counts) =>
    MARGIN.top + (1 - transformedIntensity(counts / ceiling, yScale)) * plotHeight;

  for (let step = 0; step <= 4; step += 1) {
    // Ticks are placed at even fractions of the *drawn* height and labelled with
    // the count they correspond to, so a square-root axis stays readable rather
    // than crowding every label into the bottom fifth.
    const fraction = step / 4;
    const value = ceiling * inverseIntensity(fraction, yScale);
    const y = MARGIN.top + (1 - fraction) * plotHeight;
    root.append(
      svg('line', {
        x1: MARGIN.left, y1: y, x2: WIDTH - MARGIN.right, y2: y,
        stroke: 'currentColor', 'stroke-opacity': step === 0 ? 0.5 : 0.1,
        'stroke-width': step === 0 ? 1 : 0.6,
      }),
      svg('text', {
        x: MARGIN.left - 12, y: y + 4, 'text-anchor': 'end', 'font-size': 12,
        fill: 'currentColor', 'fill-opacity': 0.6, text: formatNumber(value, 0),
      }),
    );
  }
  const tickStep = niceStep(maximum - minimum, 8);
  const firstTick = Math.ceil(minimum / tickStep) * tickStep;
  for (let value = firstTick; value <= maximum + 1e-9; value += tickStep) {
    const x = xAt(value);
    root.append(
      svg('line', {
        x1: x, y1: MARGIN.top, x2: x, y2: HEIGHT - MARGIN.bottom,
        stroke: 'currentColor', 'stroke-opacity': 0.08, 'stroke-width': 0.6,
      }),
      svg('text', {
        x, y: HEIGHT - MARGIN.bottom + 24, 'text-anchor': 'middle', 'font-size': 12,
        fill: 'currentColor', 'fill-opacity': 0.65, text: formatNumber(value, 0),
      }),
    );
  }
  root.append(
    svg('text', {
      x: MARGIN.left + plotWidth / 2, y: HEIGHT - 16, 'text-anchor': 'middle',
      'font-size': 14, fill: 'currentColor', text: '2θ (°)',
    }),
    svg('text', {
      x: 19, y: MARGIN.top + plotHeight / 2, 'text-anchor': 'middle',
      'font-size': 14, fill: 'currentColor',
      transform: `rotate(-90 19 ${MARGIN.top + plotHeight / 2})`,
      text: 'Intensity (counts)',
    }),
  );

  for (const entry of series) {
    const values = data[entry.key];
    if (!values) continue;
    const path = data.two_theta_deg.map((angle, index) =>
      `${index === 0 ? 'M' : 'L'} ${xAt(angle).toFixed(2)} ${yAt(values[index]).toFixed(2)}`,
    ).join(' ');
    root.append(svg('path', {
      d: path, fill: 'none', stroke: entry.color, 'stroke-width': entry.width,
      'vector-effect': 'non-scaling-stroke',
    }));
  }

  if (withResidual && data.residual) {
    // The residual gets its own band and its own symmetric scale. Drawn against
    // the peaks it came from it would be a flat line, and its shape is the whole
    // reason to look at it.
    const top = HEIGHT - MARGIN.bottom - residualHeight + 18;
    const half = (residualHeight - 26) / 2;
    const centre = top + half;
    let extreme = 0;
    for (const value of data.residual) extreme = Math.max(extreme, Math.abs(value));
    extreme = extreme > 0 ? extreme : 1;
    root.append(svg('line', {
      x1: MARGIN.left, y1: centre, x2: WIDTH - MARGIN.right, y2: centre,
      stroke: 'currentColor', 'stroke-opacity': 0.35, 'stroke-width': 1,
    }));
    const path = data.two_theta_deg.map((angle, index) => {
      const y = centre - (data.residual[index] / extreme) * half;
      return `${index === 0 ? 'M' : 'L'} ${xAt(angle).toFixed(2)} ${y.toFixed(2)}`;
    }).join(' ');
    root.append(
      svg('path', {
        d: path, fill: 'none', stroke: '#7c3aed', 'stroke-width': 1,
        'vector-effect': 'non-scaling-stroke',
      }),
      svg('text', {
        x: MARGIN.left, y: top - 4, 'font-size': 11, fill: 'currentColor',
        'fill-opacity': 0.6,
        text: `observed − calculated (±${formatNumber(extreme, 0)} counts)`,
      }),
    );
  }

  for (const reflection of data.reflections ?? []) {
    const x = xAt(reflection.two_theta_deg);
    if (x < MARGIN.left || x > WIDTH - MARGIN.right) continue;
    root.append(svg('line', {
      x1: x, y1: HEIGHT - MARGIN.bottom - residualHeight - 6,
      x2: x, y2: HEIGHT - MARGIN.bottom - residualHeight + 8,
      stroke: '#7c3aed', 'stroke-width': 1.2, 'stroke-opacity': 0.8,
    }));
  }
  return root;
}

/** Draw the Williamson-Hall points and the line fitted through them. */
function renderWilliamsonHall(data) {
  const root = svg('svg', {
    viewBox: `0 0 ${WIDTH} ${HEIGHT}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': 'Williamson-Hall plot',
  });
  const plotWidth = WIDTH - MARGIN.left - MARGIN.right;
  const plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom;
  // The abscissa starts at zero because the intercept *is* the answer: a plot
  // cropped to the data would hide the quantity being read off it.
  const maxX = Math.max(...data.abscissa) * 1.08;
  const maxY = Math.max(...data.ordinate, data.intercept) * 1.15;
  const xAt = (value) => MARGIN.left + (value / maxX) * plotWidth;
  const yAt = (value) => MARGIN.top + (1 - value / maxY) * plotHeight;

  for (let step = 0; step <= 4; step += 1) {
    const value = (maxY * step) / 4;
    const y = yAt(value);
    root.append(
      svg('line', {
        x1: MARGIN.left, y1: y, x2: WIDTH - MARGIN.right, y2: y,
        stroke: 'currentColor', 'stroke-opacity': step === 0 ? 0.5 : 0.1,
        'stroke-width': step === 0 ? 1 : 0.6,
      }),
      svg('text', {
        x: MARGIN.left - 12, y: y + 4, 'text-anchor': 'end', 'font-size': 12,
        fill: 'currentColor', 'fill-opacity': 0.6, text: value.toExponential(2),
      }),
    );
  }
  const tickStep = niceStep(maxX, 6);
  for (let value = 0; value <= maxX + 1e-9; value += tickStep) {
    const x = xAt(value);
    root.append(
      svg('line', {
        x1: x, y1: MARGIN.top, x2: x, y2: HEIGHT - MARGIN.bottom,
        stroke: 'currentColor', 'stroke-opacity': 0.08, 'stroke-width': 0.6,
      }),
      svg('text', {
        x, y: HEIGHT - MARGIN.bottom + 24, 'text-anchor': 'middle', 'font-size': 12,
        fill: 'currentColor', 'fill-opacity': 0.65, text: formatNumber(value, 2),
      }),
    );
  }
  root.append(
    svg('text', {
      x: MARGIN.left + plotWidth / 2, y: HEIGHT - 16, 'text-anchor': 'middle',
      'font-size': 14, fill: 'currentColor', text: '4 sin θ',
    }),
    svg('text', {
      x: 19, y: MARGIN.top + plotHeight / 2, 'text-anchor': 'middle',
      'font-size': 14, fill: 'currentColor',
      transform: `rotate(-90 19 ${MARGIN.top + plotHeight / 2})`,
      text: 'β cos θ (rad)',
    }),
    svg('line', {
      x1: xAt(0), y1: yAt(data.intercept),
      x2: xAt(maxX), y2: yAt(data.intercept + data.slope * maxX),
      stroke: '#dc2626', 'stroke-width': 1.8,
    }),
    svg('circle', {
      cx: xAt(0), cy: yAt(data.intercept), r: 5, fill: 'none',
      stroke: '#dc2626', 'stroke-width': 1.6,
    }),
    svg('text', {
      x: xAt(0) + 12, y: yAt(data.intercept) - 10, 'font-size': 12, fill: '#dc2626',
      text: `intercept = Kλ/D → ${formatNumber(data.crystallite_size_nm, 3)} nm`,
    }),
  );
  for (const [index, x] of data.abscissa.entries()) {
    root.append(svg('circle', {
      cx: xAt(x), cy: yAt(data.ordinate[index]), r: 4.5,
      fill: '#2563eb', 'fill-opacity': 0.85,
    }));
  }
  return root;
}

/** Invert `transformedIntensity`, so a transformed axis can still be labelled. */
function inverseIntensity(value, mode) {
  const safe = Math.max(0, Math.min(1, Number(value) || 0));
  if (mode === 'sqrt') return safe * safe;
  if (mode === 'log') return ((10 ** (2 * safe)) - 1) / 99;
  return safe;
}

function niceStep(range, targetTicks) {
  const rough = range / targetTicks;
  const power = 10 ** Math.floor(Math.log10(Math.max(rough, 1e-9)));
  const fraction = rough / power;
  if (fraction <= 1) return power;
  if (fraction <= 2) return 2 * power;
  if (fraction <= 5) return 5 * power;
  return 10 * power;
}

function nearestIndex(values, target) {
  let low = 0;
  let high = values.length - 1;
  while (low < high) {
    const middle = Math.floor((low + high) / 2);
    if (values[middle] < target) low = middle + 1;
    else high = middle;
  }
  if (low > 0 && Math.abs(values[low - 1] - target) < Math.abs(values[low] - target)) {
    return low - 1;
  }
  return low;
}
