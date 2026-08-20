/**
 * EBSD distributions: grain size, misorientation angle, KAM, GROD, any channel.
 *
 * A map answers *where*; this answers *how much of what*, and the two are read
 * together. The histogram is drawn here from the bin edges and counts the
 * service returns, so the picture and the table below it are the same numbers.
 *
 * Two things the drawing insists on. The bars are hoverable, because a bar
 * without its interval and its count is a shape rather than a measurement. And
 * when a random reference exists — it does for the misorientation-angle
 * distribution — it is drawn *over* the bars as a line, because the whole
 * meaning of a measured misorientation distribution is its departure from the
 * random one, and putting them in two figures would leave that comparison to the
 * reader's memory.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import { scanControls, withScan } from '../core/ebsdscan.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'ebsd_distribution',
  title: 'Distributions',
  tagline: 'Grain size, misorientation angle, KAM, GROD and the measured channels.',
};

const WIDTH = 1000;
const HEIGHT = 520;
const MARGIN = { left: 78, right: 26, top: 28, bottom: 66 };

export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'ebsd.distribution',
  );
  const examples = context.manifest.examples.filter(
    (entry) => entry.operation === operation.id,
  );

  const state = { result: null, form: null, teaches: null };
  const frame = plotFrame({ title: 'Distribution' });
  const details = el('div');
  const formHost = el('div');
  const legend = el('div.legend');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Plot the distribution',
    onclick: () => run(),
  });

  const scan = scanControls({ onChange: () => run(), showError: context.showError });

  context.rail.append(
    scan.element,
    formHost,
    runButton,
    examples.length
      ? el('details.group', {}, [
          el('summary', { text: 'Try an example' }),
          el('div.group__body', {}, [
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
        ])
      : null,
  );

  frame.setControls(legend);
  context.stage.append(frame.element, details);

  function renderControls(initial = {}) {
    state.form = buildForm(operation, { initial });
    formHost.replaceChildren(state.form.element);
    for (const field of state.form.element.querySelectorAll('.field')) {
      if (field.querySelector('[id^="ctl-scan_file-"]')) field.hidden = true;
    }
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Counting…';
    state.form.clearErrors();
    try {
      state.result = await call(operation.id, withScan(state.form.values()));
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
      state.teaches = null;
    } catch (error) {
      if (error?.field === 'scan_file') {
        scan.setStatus(error.message);
        context.showError(error);
      } else if (!state.form.showError(error)) {
        context.showError(error);
      } else {
        context.showError(error, { quiet: true });
      }
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Plot the distribution';
    }
  }

  function draw() {
    const data = state.result.data;
    const edges = data.edges;
    const counts = data.counts;
    const total = counts.reduce((sum, value) => sum + value, 0) || 1;
    const fractions = counts.map((value) => value / total);
    const reference = data.reference?.fractions ?? null;
    const peak = Math.max(...fractions, ...(reference ?? [0])) || 1;

    const plotWidth = WIDTH - MARGIN.left - MARGIN.right;
    const plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom;
    const lowest = edges[0];
    const highest = edges[edges.length - 1];
    const span = highest - lowest || 1;
    const x = (value) => MARGIN.left + ((value - lowest) / span) * plotWidth;
    const y = (value) => MARGIN.top + plotHeight * (1 - value / peak);

    const root = svg('svg', {
      viewBox: `0 0 ${WIDTH} ${HEIGHT}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': `${data.label} distribution`,
    });

    // Axes first, so nothing is drawn under a bar it should be over.
    root.append(
      svg('line', {
        x1: MARGIN.left,
        y1: MARGIN.top + plotHeight,
        x2: MARGIN.left + plotWidth,
        y2: MARGIN.top + plotHeight,
        stroke: 'currentColor',
        'stroke-opacity': 0.6,
      }),
      svg('line', {
        x1: MARGIN.left,
        y1: MARGIN.top,
        x2: MARGIN.left,
        y2: MARGIN.top + plotHeight,
        stroke: 'currentColor',
        'stroke-opacity': 0.6,
      }),
    );

    for (let tick = 0; tick <= 5; tick += 1) {
      const value = lowest + (span * tick) / 5;
      root.append(
        svg('text', {
          x: x(value),
          y: MARGIN.top + plotHeight + 22,
          'text-anchor': 'middle',
          'font-size': 15,
          fill: 'currentColor',
          'fill-opacity': 0.75,
          text: formatNumber(value, 3),
        }),
      );
      const height = (peak * tick) / 5;
      root.append(
        svg('text', {
          x: MARGIN.left - 10,
          y: y(height) + 5,
          'text-anchor': 'end',
          'font-size': 15,
          fill: 'currentColor',
          'fill-opacity': 0.75,
          text: formatNumber(100 * height, 1),
        }),
      );
    }

    root.append(
      svg('text', {
        x: MARGIN.left + plotWidth / 2,
        y: HEIGHT - 14,
        'text-anchor': 'middle',
        'font-size': 17,
        fill: 'currentColor',
        text: data.units ? `${data.label} / ${data.units}` : data.label,
      }),
      svg('text', {
        x: 20,
        y: MARGIN.top + plotHeight / 2,
        'text-anchor': 'middle',
        'font-size': 17,
        fill: 'currentColor',
        transform: `rotate(-90 20 ${MARGIN.top + plotHeight / 2})`,
        text: data.weighted_by_length ? 'Length fraction / %' : 'Fraction / %',
      }),
    );

    fractions.forEach((fraction, index) => {
      const left = x(edges[index]);
      const right = x(edges[index + 1]);
      const bar = svg('rect', {
        x: left + 1,
        y: y(fraction),
        width: Math.max(right - left - 2, 1),
        height: Math.max(MARGIN.top + plotHeight - y(fraction), 0),
        fill: 'var(--accent)',
        'fill-opacity': 0.72,
      });
      root.append(bar);
      frame.hoverable(bar, {
        From: edges[index],
        To: edges[index + 1],
        [data.weighted_by_length ? 'Length / µm' : 'Count']: counts[index],
        Fraction: fraction,
      });
    });

    if (reference) {
      // A step line, not a smooth one: it is a histogram of the same bins, and
      // drawing it as a curve would suggest a density estimate nobody made.
      const points = [];
      reference.forEach((fraction, index) => {
        points.push([x(edges[index]), y(fraction)], [x(edges[index + 1]), y(fraction)]);
      });
      root.append(
        svg('polyline', {
          points: points.map(([px, py]) => `${px},${py}`).join(' '),
          fill: 'none',
          stroke: 'var(--warn)',
          'stroke-width': 2.4,
          'stroke-dasharray': '7 5',
        }),
      );
    }

    frame.configure({
      toData: (px, py) => {
        if (px < MARGIN.left || px > MARGIN.left + plotWidth) return null;
        if (py < MARGIN.top || py > MARGIN.top + plotHeight) return null;
        return {
          x: lowest + ((px - MARGIN.left) / plotWidth) * span,
          y: peak * (1 - (py - MARGIN.top) / plotHeight),
        };
      },
      formatCursor: (point) =>
        `${formatNumber(point.x, 4)}${data.units ? ` ${data.units}` : ''} · ` +
        `${formatNumber(100 * point.y, 2)} %`,
    });
    frame.setContent(root);

    legend.replaceChildren(
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: 'background: var(--accent)' }),
        el('span', { text: `${data.label} (${data.statistics.count.toLocaleString()} ${data.population})` }),
      ]),
      ...(reference
        ? [
            el('span.legend__item', {}, [
              el('span.legend__swatch', { style: 'background: var(--warn)' }),
              el('span', { text: data.reference.label }),
            ]),
          ]
        : []),
    );

    const statistics = data.statistics;
    const unit = data.units ? ` ${data.units}` : '';
    frame.setStatus(
      `mean ${formatNumber(statistics.mean, 4)}${unit} · ` +
        `median ${formatNumber(statistics.median, 4)}${unit} · ` +
        `10th to 90th percentile ${formatNumber(statistics.p10, 4)} to ` +
        `${formatNumber(statistics.p90, 4)}${unit} · ` +
        `${statistics.count.toLocaleString()} ${data.population}` +
        (data.weighted_by_length ? ' · weighted by boundary length' : ''),
    );
  }

  renderControls();
  // The first example rather than the bare defaults, as every other panel does:
  // the practice dataset a control defaults to is the bicrystal, which is the
  // right first map and a poor first distribution — two grains make one bar.
  if (examples.length) loadExample(examples[0]);
  else run();

  return { help: () => operation };
}
