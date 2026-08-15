/** Powder XRD: an indexed diffractogram whose peaks remain inspectable. */

import { call } from '../core/api.js';
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

export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'xrd.powder_pattern',
  );
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = {
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
    text: 'Simulate XRD pattern',
    onclick: () => run(),
  });

  function redraw() {
    if (state.result) draw();
  }

  context.rail.append(
    formHost,
    runButton,
    appearanceControls(state.appearance, redraw),
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'These standards isolate extinction rules, wavelength, doublet splitting and hexagonal metrics.',
        }),
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
    state.form = buildForm(operation, { initial });
    formHost.replaceChildren(state.form.element);
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Simulating…';
    state.form.clearErrors();
    try {
      state.result = await call(operation.id, state.form.values());
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Simulate XRD pattern';
    }
  }

  function draw() {
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
  return { help: () => operation };
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
  const yScale = el('select', {
    oninput: (event) => {
      style.yScale = event.currentTarget.value;
      onChange();
    },
  }, [
    el('option', { value: 'linear', text: 'Linear', selected: style.yScale === 'linear' }),
    el('option', { value: 'sqrt', text: 'Square-root', selected: style.yScale === 'sqrt' }),
    el('option', { value: 'log', text: 'Log-like', selected: style.yScale === 'log' }),
  ]);

  body.append(
    el('p.field__help', {
      text: 'Display controls redraw the existing profile. They never change peak positions, integrated intensities or exports.',
    }),
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
    el('label.field', {}, [
      el('span.field__label', { text: 'Vertical display scale' }),
      yScale,
      el('span.field__hint', { text: 'Square-root and log-like views reveal weak peaks without changing the data.' }),
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
