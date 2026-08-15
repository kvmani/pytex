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
import {
  markerNode,
  markerRadius,
  markerStyleControl,
  productColor,
  variantMarkerStyle,
} from '../core/visualstyle.js';

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
export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'diffraction.composite_saed',
  );
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = {
    result: null,
    teaches: null,
    form: null,
    hidden: new Set(),
    // Which pattern the legend was built for, so a redraw that only changes
    // which sources are shown updates the buttons instead of replacing them.
    legendFor: null,
    appearance: null,
  };

  const frame = plotFrame({
    title: 'Composite pattern',
    units: 'Å⁻¹',
    toData: null, // set once a pattern exists and the camera constant is known
    toolbar: [],
  });

  const legend = el('div.legend');
  const details = el('div');
  // The legend is a control, so it rides inside the frame rather than under it:
  // toggling a source and seeing the drawing change must not need a scroll.
  frame.setControls(legend);
  context.stage.append(frame.element, details);

  const formHost = el('div');
  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Simulate pattern',
    onclick: () => run(),
  });
  const appearance = markerStyleControl({
    onChange: (style) => {
      state.appearance = style;
      if (state.result) draw();
    },
  });
  state.appearance = appearance.style;

  context.rail.append(
    formHost,
    runButton,
    appearance.element,
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

    frame.setContent(
      renderPattern(data, {
        scale,
        frame,
        hidden: state.hidden,
        appearance: state.appearance,
      }),
    );
    // Built once per pattern, updated in place thereafter. Rebuilding it on
    // every redraw destroys the button that was just pressed, and the browser
    // then moves focus to the body — so a keyboard user who hides a variant
    // loses their place and has to tab back through the whole page.
    if (state.legendFor !== data) buildLegend(data);
    else updateLegend();
    const shown = data.spots.filter((spot) => !state.hidden.has(sourceKey(spot))).length;
    frame.setStatus(
      `${shown} of ${data.spots.length} spots · camera constant ` +
        `${formatNumber(cameraConstant, 1)} mm·Å · ${formatNumber(data.beam_energy_kev, 0)} kV · ` +
        'hover a spot for its indices, spacing, intensity and variant',
    );
  }

  /** Build the legend for a new pattern. Called once per pattern, not per redraw. */
  function buildLegend(data) {
    state.legendFor = data;
    const sourceKeys = data.sources.map((source) => `${source.source}|${source.variant}`);
    const parentSource = data.sources.find((source) => source.source === 'parent');
    const parentKey = parentSource ? `${parentSource.source}|${parentSource.variant}` : null;
    const focusSelect = el('select.legend__focus', {
      'aria-label': 'Focus on one variant while retaining the parent pattern',
      title: 'Show the parent and one selected variant; use the legend chips for finer control.',
      onchange: (event) => {
        const focusKey = event.currentTarget.value;
        if (!focusKey) return;
        state.hidden = new Set(
          sourceKeys.filter((key) => key !== focusKey && key !== parentKey),
        );
        event.currentTarget.value = '';
        draw();
      },
    });
    focusSelect.append(
      el('option', { value: '', text: 'Focus a variant…', selected: true }),
      ...data.sources
        .filter((source) => source.source !== 'parent')
        .map((source) =>
          el('option', {
            value: `${source.source}|${source.variant}`,
            text: `${source.label} + parent`,
          }),
        ),
    );

    const visibilityTools = el('div.legend__toolbar', {}, [
      el('span.legend__guide', { text: 'Display' }),
      el('button.button.button--small', {
        type: 'button',
        text: 'Show all',
        title: 'Restore every source in this composite pattern.',
        onclick: () => {
          state.hidden.clear();
          draw();
        },
      }),
      ...(parentSource
        ? [
            el('button.button.button--small', {
              type: 'button',
              text: 'Parent only',
              title: 'Hide every product variant and retain the parent reference pattern.',
              onclick: () => {
                state.hidden = new Set(sourceKeys.filter((key) => key !== parentKey));
                draw();
              },
            }),
          ]
        : []),
      focusSelect,
      el('span.legend__guide', { text: 'Click a chip to toggle one source.' }),
    ]);

    const items = el(
      'div.legend__items',
      {},
      data.sources.map((source, index) => {
        const key = `${source.source}|${source.variant}`;
        const marker = source.source === 'parent'
          ? { shape: state.appearance.shape, scale: 1 }
          : variantMarkerStyle(index, state.appearance);
        return el(
          'button.legend__item',
          {
            type: 'button',
            dataset: { key, source: source.source, colorIndex: index },
            onclick: () => {
              if (state.hidden.has(key)) state.hidden.delete(key);
              else state.hidden.add(key);
              draw();
            },
          },
          [
            el('span.legend__swatch', {
              dataset: { shape: marker.shape, fill: state.appearance.fill },
              style: `--swatch-color:${
                source.source === 'parent'
                  ? state.appearance.parentColor
                  : productColor(index, state.appearance)
              };--swatch-scale:${marker.scale}`,
            }),
            el('span', { text: `${source.label} (${source.spots.length})` }),
          ],
        );
      }),
    );
    legend.replaceChildren(
      visibilityTools,
      items,
      el('span.legend__item', {}, [
        el('span.legend__swatch', {
          dataset: { shape: 'transmitted' },
          style: '--swatch-color:var(--ink);--swatch-scale:1',
        }),
        el('span', { text: 'Transmitted beam (000)' }),
      ]),
    );
    updateLegend();
  }

  /** Reflect the hidden set onto the existing buttons, without replacing them. */
  function updateLegend() {
    for (const button of legend.querySelectorAll('button.legend__item')) {
      const hidden = state.hidden.has(button.dataset.key);
      button.setAttribute('aria-pressed', String(!hidden));
      button.title = hidden ? 'Show these spots' : 'Hide these spots';
      const swatch = button.querySelector('.legend__swatch');
      const marker = button.dataset.source === 'parent'
        ? { shape: state.appearance.shape, scale: 1 }
        : variantMarkerStyle(Number(button.dataset.colorIndex), state.appearance);
      swatch.dataset.shape = marker.shape;
      swatch.dataset.fill = state.appearance.fill;
      swatch.style.setProperty('--swatch-scale', marker.scale);
      const color = button.dataset.source === 'parent'
        ? state.appearance.parentColor
        : productColor(Number(button.dataset.colorIndex), state.appearance);
      swatch.style.setProperty('--swatch-color', color);
      swatch.style.background = state.appearance.fill === 'outline' ? 'transparent' : color;
    }
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);

  return { help: () => operation };
}

function sourceKey(spot) {
  return `${spot.source}|${spot.variant}`;
}

function renderPattern(data, { scale, frame, hidden, appearance }) {
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

  const sourceStyles = new Map();
  data.sources.forEach((source, index) => {
    const marker = source.source === 'parent'
      ? { shape: appearance.shape, scale: 1 }
      : variantMarkerStyle(index, appearance);
    sourceStyles.set(
      `${source.source}|${source.variant}`,
      {
        color: source.source === 'parent'
          ? appearance.parentColor
          : productColor(index, appearance),
        ...marker,
      },
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
    const sourceStyle = sourceStyles.get(key) ?? {
      color: 'currentColor', shape: appearance.shape, scale: 1,
    };
    const x = spot.detector_x_mm * scale;
    const y = -spot.detector_y_mm * scale;
    // Radius by the fourth root of intensity: a linear map makes everything
    // below a tenth of the maximum invisible, and those are exactly the
    // superlattice and variant spots the panel exists to show.
    const node = markerNode(svg, {
      x,
      y,
      radius: markerRadius(spot.relative_intensity, appearance) * sourceStyle.scale,
      shape: sourceStyle.shape,
      color: sourceStyle.color,
      hollow: appearance.fill === 'outline' || spot.double_diffraction,
      dashed: spot.double_diffraction,
    });
    node.setAttribute('fill-opacity', spot.double_diffraction ? '0.35' : '0.92');
    root.append(node);
    frame.hoverable(node, spot, columns);
  }

  // The direct beam, marked so the centre of the pattern is never ambiguous.
  root.append(
    svg('circle', { cx: 0, cy: 0, r: 2.2, fill: 'none', stroke: 'currentColor', 'stroke-width': 0.75 }),
    svg('circle', { cx: 0, cy: 0, r: 0.8, fill: 'currentColor' }),
    svg('line', { x1: -3.8, y1: 0, x2: -2.5, y2: 0, stroke: 'currentColor', 'stroke-width': 0.45 }),
    svg('line', { x1: 2.5, y1: 0, x2: 3.8, y2: 0, stroke: 'currentColor', 'stroke-width': 0.45 }),
    svg('line', { x1: 0, y1: -3.8, x2: 0, y2: -2.5, stroke: 'currentColor', 'stroke-width': 0.45 }),
    svg('line', { x1: 0, y1: 2.5, x2: 0, y2: 3.8, stroke: 'currentColor', 'stroke-width': 0.45 }),
    svg('text', {
      x: 4.2, y: -3.2, 'font-size': 3.2, fill: 'currentColor',
      'font-weight': 600, text: '(000) transmitted',
    }),
  );
  return root;
}
