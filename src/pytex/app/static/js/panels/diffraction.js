/**
 * The diffraction panel: a composite pattern you can interrogate spot by spot.
 *
 * The reason this panel is worth building rather than exporting a figure is the
 * hover. A Kurdjumov-Sachs composite pattern down [001] carries two hundred
 * spots from twenty-five different reciprocal lattices, and the only question a
 * researcher ever asks of it is "what is *that* one?" — indices, spacing,
 * intensity, and above all which variant. Every spot therefore carries its full
 * row, taken from the same table the CSV export writes.
 *
 * The cursor readout reports reciprocal-space radius as well as detector
 * position, because the two are what a pattern is measured in: a plate is read
 * in millimetres and interpreted in inverse angstroms, and the camera constant
 * relating them is exactly the thing people get wrong.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';
import { call } from '../core/api.js';

export const panel = {
  id: 'diffraction',
  title: 'Diffraction',
  tagline: 'Composite SAED of a parent phase and its product variants.',
};

const VIEW = 100;

/**
 * Colours for the sources in a composite pattern.
 *
 * The parent is deliberately the one achromatic entry: it is the reference
 * everything else is measured against, and a reader should never have to ask
 * which lattice is the parent. Variants cycle through a hue wheel, which is
 * honest about there being no meaningful order to them.
 */
const PARENT_COLOR = 'var(--ink)';

function variantColor(index) {
  const hue = (index * 137.508) % 360; // golden angle: adjacent variants stay distinguishable
  return `hsl(${hue.toFixed(1)} 70% 52%)`;
}

export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'diffraction.composite_saed',
  );
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = { result: null, teaches: null, form: null, hidden: new Set() };

  const frame = plotFrame({
    title: 'Composite pattern',
    units: 'Å⁻¹',
    toData: null, // set once a pattern exists and the camera constant is known
    toolbar: [],
  });

  const legend = el('div.legend');
  const details = el('div');
  context.stage.append(frame.element, legend, details);

  const formHost = el('div');
  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Simulate pattern',
    onclick: () => run(),
  });

  context.rail.append(
    formHost,
    runButton,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'Run the first two in order: the same zone axis with 24 variants and then with one.',
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
      const result = await call(operation.id, state.form.values());
      state.result = result;
      state.hidden = new Set();
      draw();
      renderResult(details, result, { teaches: state.teaches });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Simulate pattern';
    }
  }

  function draw() {
    const data = state.result.data;
    const cameraConstant = data.camera_constant_mm_angstrom || 1;
    const scale = (VIEW * 0.94) / (data.detector_radius_mm || 1);

    // The readout is in both currencies at once: millimetres are what a plate is
    // measured in, inverse angstroms are what the measurement means.
    frame.configure({
      toData: (x, y) => ({ x: x / scale, y: -y / scale }),
      formatCursor: (point) => {
        const radiusMm = Math.hypot(point.x, point.y);
        const g = radiusMm / cameraConstant;
        const spacing = g > 1e-9 ? 1 / g : Infinity;
        return (
          `${formatNumber(point.x, 2)}, ${formatNumber(point.y, 2)} mm · ` +
          `|g| ${formatNumber(g, 4)} Å⁻¹ · d ${
            Number.isFinite(spacing) ? formatNumber(spacing, 4) : '∞'
          } Å`
        );
      },
    });

    frame.setContent(renderPattern(data, { scale, frame, hidden: state.hidden }));
    renderLegend(data);
    const shown = data.spots.filter((spot) => !state.hidden.has(sourceKey(spot))).length;
    frame.setStatus(
      `${shown} of ${data.spots.length} spots · camera constant ` +
        `${formatNumber(cameraConstant, 1)} mm·Å · ${formatNumber(data.beam_energy_kev, 0)} kV · ` +
        'hover a spot for its indices, spacing, intensity and variant',
    );
  }

  function renderLegend(data) {
    legend.replaceChildren(
      ...data.sources.map((source, index) => {
        const key = `${source.source}|${source.variant}`;
        const hidden = state.hidden.has(key);
        return el(
          'button.legend__item',
          {
            type: 'button',
            'aria-pressed': String(!hidden),
            title: hidden ? 'Show these spots' : 'Hide these spots',
            onclick: () => {
              if (hidden) state.hidden.delete(key);
              else state.hidden.add(key);
              draw();
            },
          },
          [
            el('span.legend__swatch', {
              style: `background:${source.source === 'parent' ? PARENT_COLOR : variantColor(index)}`,
            }),
            el('span', { text: `${source.label} (${source.spots.length})` }),
          ],
        );
      }),
    );
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);

  return { help: () => operation };
}

function sourceKey(spot) {
  return `${spot.source}|${spot.variant}`;
}

function renderPattern(data, { scale, frame, hidden }) {
  const root = svg('svg', {
    viewBox: `${-VIEW} ${-VIEW} ${2 * VIEW} ${2 * VIEW}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': 'Composite diffraction pattern',
  });

  // Reciprocal-space rings at whole inverse-angstrom radii, so the scale of the
  // pattern is readable without measuring anything.
  const cameraConstant = data.camera_constant_mm_angstrom || 1;
  const maxG = (data.detector_radius_mm || 1) / cameraConstant;
  for (let g = 0.5; g <= maxG; g += 0.5) {
    root.append(
      svg('circle', {
        cx: 0, cy: 0, r: g * cameraConstant * scale,
        fill: 'none',
        stroke: 'currentColor',
        'stroke-opacity': 0.14,
        'stroke-width': 0.3,
        'stroke-dasharray': '2 3',
      }),
      svg('text', {
        x: 0, y: -g * cameraConstant * scale - 1,
        'text-anchor': 'middle',
        'font-size': 3,
        fill: 'currentColor',
        'fill-opacity': 0.45,
        text: `${g.toFixed(1)} Å⁻¹`,
      }),
    );
  }

  const colors = new Map();
  data.sources.forEach((source, index) => {
    colors.set(
      `${source.source}|${source.variant}`,
      source.source === 'parent' ? PARENT_COLOR : variantColor(index),
    );
  });

  const columns = data.columns;
  // Draw weak spots first so a strong reflection is never buried under a faint
  // one it overlaps — the same reason the crystal viewer depth-sorts.
  const ordered = [...data.spots].sort(
    (left, right) => left.relative_intensity - right.relative_intensity,
  );
  for (const spot of ordered) {
    const key = sourceKey(spot);
    if (hidden.has(key)) continue;
    const x = spot.detector_x_mm * scale;
    const y = -spot.detector_y_mm * scale;
    // Radius by the fourth root of intensity: a linear map makes everything
    // below a tenth of the maximum invisible, and those are exactly the
    // superlattice and variant spots the panel exists to show.
    const radius = 0.7 + 2.6 * Math.pow(Math.max(spot.relative_intensity, 0), 0.25);
    const node = svg('circle', {
      cx: x, cy: y, r: radius,
      fill: colors.get(key) ?? 'currentColor',
      'fill-opacity': spot.double_diffraction ? 0.35 : 0.92,
      stroke: spot.double_diffraction ? colors.get(key) ?? 'currentColor' : 'none',
      'stroke-width': 0.4,
      'stroke-dasharray': spot.double_diffraction ? '1 1' : null,
    });
    root.append(node);
    frame.hoverable(node, spot, columns);
  }

  // The direct beam, marked so the centre of the pattern is never ambiguous.
  root.append(
    svg('circle', { cx: 0, cy: 0, r: 1.6, fill: 'none', stroke: 'currentColor', 'stroke-width': 0.5 }),
    svg('line', { x1: -3, y1: 0, x2: 3, y2: 0, stroke: 'currentColor', 'stroke-width': 0.3 }),
    svg('line', { x1: 0, y1: -3, x2: 0, y2: 3, stroke: 'currentColor', 'stroke-width': 0.3 }),
  );
  return root;
}
