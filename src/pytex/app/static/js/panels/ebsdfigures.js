/**
 * Discrete pole and inverse pole figures of an EBSD scan.
 *
 * The scatter a contoured figure is an estimate of. Every drawn point is one
 * measured orientation's pole, projected stereographically onto the upper
 * hemisphere — no kernel, no interpolation, and therefore no width that somebody
 * had to choose.
 *
 * The markers are small and slightly transparent on purpose. Overlap is
 * information here: where the scatter saturates is where the texture is, and a
 * fully opaque marker would throw that away by making one point look the same as
 * fifty. It is the honest opposite of a contour — the density is visible without
 * ever being estimated.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import { activeDataset, adoptForm, scanControls, withScan } from '../core/ebsdscan.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'ebsd_figures',
  title: 'Pole figures',
  tagline: 'Discrete pole and inverse pole figures: the measurement, not a contour of it.',
};

/** Half-width of the drawing area, in viewBox units. */
const VIEW = 320;
const PAD = 46;

export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'ebsd.discrete_figure',
  );
  const examples = context.manifest.examples.filter(
    (entry) => entry.operation === operation.id,
  );

  const state = { result: null, form: null, teaches: null };
  const frame = plotFrame({ title: 'Discrete figure' });
  const details = el('div');
  const formHost = el('div');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Plot the figure',
    onclick: () => run(),
  });

  const scan = scanControls({
    operation,
    onChange: () => run(),
    showError: context.showError,
  });

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
                el('button.example', { type: 'button', onclick: () => loadExample(example, { chosen: true }) }, [
                  el('strong', { text: example.title }),
                  el('span', { text: example.summary }),
                ]),
              ),
            ),
          ]),
        ])
      : null,
  );

  context.stage.append(frame.element, details);

  /**
   * Rebuild the form.
   *
   * `chosen` says whether these values came from a user action. It decides one
   * thing: whether a `dataset` among them may move the whole EBSD workspace.
   * Clicking an example is choosing its dataset; this panel opening itself on
   * an example is not, and the workspace's own choice wins there.
   */
  function renderControls(initial = {}, { chosen = false } = {}) {
    const dataset = activeDataset();
    state.form = buildForm(operation, {
      initial: chosen ? { dataset, ...initial } : { ...initial, dataset },
    });
    formHost.replaceChildren(state.form.element);
    // The scan file and the dataset are workspace-wide and presented beside
    // the form, so their generated controls are hidden rather than shown twice.
    adoptForm(state.form, { adoptDataset: chosen });
  }

  function loadExample(example, { chosen = false } = {}) {
    state.teaches = example.teaches;
    renderControls(example.request, { chosen });
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Projecting…';
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
      runButton.textContent = 'Plot the figure';
    }
  }

  function draw() {
    const data = state.result.data;
    const size = VIEW + PAD;
    const root = svg('svg', {
      viewBox: `${-size} ${-size} ${2 * size} ${2 * size}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': data.label,
    });

    root.append(
      svg('circle', {
        cx: 0,
        cy: 0,
        r: VIEW,
        fill: 'none',
        stroke: 'currentColor',
        'stroke-opacity': 0.55,
        'stroke-width': 1.6,
      }),
      svg('line', {
        x1: -VIEW, y1: 0, x2: VIEW, y2: 0,
        stroke: 'currentColor', 'stroke-opacity': 0.16,
      }),
      svg('line', {
        x1: 0, y1: -VIEW, x2: 0, y2: VIEW,
        stroke: 'currentColor', 'stroke-opacity': 0.16,
      }),
    );

    // The specimen axes, named. A pole figure whose frame is not on it can be
    // read the wrong way round without anything looking wrong.
    if (data.kind === 'pole') {
      for (const [label, x, y] of [
        ['X', VIEW + 26, 6],
        ['Y', 0, -VIEW - 18],
      ]) {
        root.append(
          svg('text', {
            x, y,
            'text-anchor': 'middle',
            'font-size': 22,
            fill: 'currentColor',
            'fill-opacity': 0.65,
            text: label,
          }),
        );
      }
    }

    // Marker size falls with the count so a dense figure stays a scatter rather
    // than becoming a disc, and the transparency lets overlap show density.
    const count = data.points.length || 1;
    const radius = Math.max(1.6, Math.min(5.0, 900 / Math.sqrt(count) / 6));
    for (const point of data.points) {
      root.append(
        svg('circle', {
          cx: point.x * VIEW,
          cy: -point.y * VIEW,
          r: radius,
          fill: 'var(--accent)',
          'fill-opacity': 0.45,
        }),
      );
    }

    frame.configure({
      toData: (x, y) => {
        const rx = x / VIEW;
        const ry = -y / VIEW;
        if (Math.hypot(rx, ry) > 1) return null;
        return { x: rx, y: ry };
      },
      formatCursor: (point) => {
        // The polar angle a projected radius stands for: rho = 2 arctan(r), the
        // inverse of the stereographic projection this figure is drawn in.
        const r = Math.hypot(point.x, point.y);
        const polar = (2 * Math.atan(r) * 180) / Math.PI;
        const azimuth = (Math.atan2(point.y, point.x) * 180) / Math.PI;
        return (
          `${formatNumber(point.x, 3)}, ${formatNumber(point.y, 3)} · ` +
          `${formatNumber(polar, 1)}° from the centre at ${formatNumber(azimuth, 1)}°`
        );
      },
    });
    frame.setContent(root);
    frame.setStatus(
      `${data.drawn_points.toLocaleString()} projected points from ` +
        `${data.measurement_points.toLocaleString()} measurement point(s)` +
        (data.subsampled ? ` of ${data.scan_points.toLocaleString()}` : '') +
        ' · stereographic, upper hemisphere · this is the measurement, not a density estimate',
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
