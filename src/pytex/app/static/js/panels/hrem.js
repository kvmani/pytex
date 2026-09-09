/**
 * HRTEM Simulation and Objective Lens Contrast Transfer Function panel.
 *
 * Two views: a multislice phase-contrast micrograph with its Thon-ring power
 * spectrum, and the objective lens transfer function. They are sub-tabs in the
 * shell's own strip rather than a pair of buttons on the rail, because a view
 * tab and a workspace sub-tab are the same affordance and a strip drawn on the
 * stage would push the plot card past the bottom of it.
 *
 * The CTF view draws the azimuthal band whenever the lens carries a residual
 * non-round aberration. A single radial cut of a lens with astigmatism, coma or
 * trefoil describes one direction in the image; drawing the cut alone would let
 * a reader take a directional number for the resolution of the instrument.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import { explainer } from '../core/explainer.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'tem_hrem',
  title: 'HRTEM Simulation',
  tagline: 'Multislice simulation and CTF transfer for double-corrected TEM.',
};

const VIEWS = [
  {
    id: 'tem.simulate_hrem',
    title: 'Micrograph',
    summary: 'Phase-contrast image and its Thon-ring power spectrum.',
    action: 'Run HRTEM simulation',
  },
  {
    id: 'tem.ctf_calculator',
    title: 'Transfer function',
    summary: 'Objective lens CTF, coherence envelopes and azimuthal anisotropy.',
    action: 'Calculate CTF',
  },
];

// Curve colours. Distinct in hue and in dash pattern, so the plate stays
// readable in monochrome print and to a reader who cannot separate the hues.
const COLOR_TRANSFER = '#2563eb';
const COLOR_ENVELOPE = '#f59e0b';
const COLOR_UNDAMPED = '#94a3b8';
const COLOR_BAND = '#818cf8';
const COLOR_POINT_RES = '#10b981';
const COLOR_INFO_LIMIT = '#ef4444';

export function mount(context) {
  const operations = VIEWS.map((view) => ({
    ...view,
    operation: context.manifest.operations.find((entry) => entry.id === view.id),
  })).filter((view) => view.operation);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const state = {
    view: operations[0],
    results: {},
    forms: {},
    teaches: null,
  };

  const simFrame = plotFrame({ title: 'Simulated Micrograph', units: 'Å', toolbar: [] });
  const fftFrame = plotFrame({ title: 'Power Spectrum (Thon Rings)', units: 'Å⁻¹', toolbar: [] });
  const ctfFrame = plotFrame({
    title: 'Contrast Transfer Function & Envelopes',
    units: 'Å⁻¹',
    toolbar: [],
  });

  // The figures own the stage height and the reading matter follows below it,
  // which is what every other two-figure panel does: a result card inside the
  // figure stage would be drawn over the figures it describes.
  const simStage = el('div.hrem-sim-stage', {}, [simFrame.element, fftFrame.element]);
  const ctfStage = el('div.hrem-ctf-stage', {}, [ctfFrame.element]);
  const legend = el('div.hrem-legend');
  const details = el('div');
  context.stage.append(simStage, ctfStage, legend, details);

  const formHosts = {};
  for (const view of operations) {
    const host = el('div.hrem-form');
    state.forms[view.id] = buildForm(view.operation, { initial: {}, onChange: () => {} });
    host.append(state.forms[view.id].element);
    formHosts[view.id] = host;
  }

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: state.view.action,
    onclick: () => runCurrent(),
  });

  const exampleItems = examples.map((ex) =>
    el(
      'button.button.button--ghost.button--block',
      {
        type: 'button',
        text: ex.title,
        title: ex.summary,
        onclick: () => loadExample(ex),
      },
      [el('small.text-muted', { text: ` ${ex.summary}` })],
    ),
  );

  const exampleGroup = el('details.group', { open: false }, [
    el('summary', { text: 'Try an example' }),
    el('div.group__body', {}, [
      explainer(
        'Preset scenarios demonstrate high-resolution TEM contrast, CTF calibration and the '
          + 'anisotropy a residual aberration imposes.',
        { label: 'HRTEM Presets' },
      ),
      el('div.stack.stack--tight', {}, exampleItems),
    ]),
  ]);

  context.rail.append(...Object.values(formHosts), runButton, exampleGroup);

  const markActiveTab = context.setViews(operations, (id) => selectView(id));

  function selectView(id) {
    const chosen = operations.find((entry) => entry.id === id);
    if (!chosen || chosen === state.view) return;
    state.view = chosen;
    state.teaches = null;
    applyView();
  }

  function applyView() {
    const isSim = state.view.id === 'tem.simulate_hrem';
    markActiveTab(state.view.id);
    simStage.hidden = !isSim;
    ctfStage.hidden = isSim;
    for (const [id, host] of Object.entries(formHosts)) {
      host.hidden = id !== state.view.id;
    }
    runButton.textContent = state.view.action;
    runButton.disabled = false;
    updateLegend();
    updateDetails();
  }

  async function runCurrent() {
    // Each run remembers the view it was launched for: a request still in flight
    // when the reader switches tabs must not be drawn by the other view's
    // renderer, which expects different keys.
    const launched = state.view;
    const form = state.forms[launched.id];
    if (!form) return;
    runButton.disabled = true;
    runButton.textContent = 'Working…';
    form.clearErrors();
    try {
      const result = await call(launched.id, form.values());
      state.results[launched.id] = result;
      if (state.view !== launched) return;
      if (launched.id === 'tem.simulate_hrem') drawSimulation();
      else drawCTF();
      updateLegend();
      updateDetails();
    } catch (error) {
      if (state.view !== launched) return;
      if (!form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      if (state.view === launched) {
        runButton.disabled = false;
        runButton.textContent = launched.action;
      }
    }
  }

  function loadExample(ex) {
    const target = operations.find((entry) => entry.id === ex.operation);
    if (!target) return;
    state.view = target;
    applyView();
    state.forms[ex.operation]?.setValues(ex.request);
    state.teaches = ex.teaches;
    runCurrent();
  }

  function drawSimulation() {
    const res = state.results['tem.simulate_hrem'];
    if (!res || !res.data) return;

    const imgData = res.data.image_png;
    if (imgData) {
      const img = el('img.hrem-image', { src: imgData, alt: 'Simulated HRTEM micrograph' });
      const extent = res.data.extent_angstrom || [10, 10];
      const caption = el('div.text-center.text-muted', {
        text:
          `Field of view: ${formatNumber(extent[0], 1)} × ${formatNumber(extent[1], 1)} Å `
          + `(pixel size: ${formatNumber(res.data.pixel_size_angstrom, 3)} Å/px)`,
      });
      simFrame.setContent(el('div.stack', {}, [img, caption]));
    }

    const psData = res.data.power_spectrum_png;
    if (psData) {
      const psImg = el('img.hrem-image', {
        src: psData,
        alt: '2D FFT Power Spectrum (Thon Rings)',
      });
      const psCaption = el('div.text-center.text-muted', {
        text:
          `Point resolution: ${formatNumber(res.data.point_resolution_angstrom, 2)} Å | `
          + `Information limit: ${formatNumber(res.data.information_limit_angstrom, 2)} Å`,
      });
      fftFrame.setContent(el('div.stack', {}, [psImg, psCaption]));
    }
  }

  function drawCTF() {
    const res = state.results['tem.ctf_calculator'];
    if (!res || !res.data) return;

    const q = res.data.spatial_frequencies || [];
    const transfer = res.data.transfer_function || [];
    const undamped = res.data.undamped_ctf || [];
    const envelope = res.data.total_envelope || [];
    const bandMin = res.data.azimuthal_transfer_min || [];
    const bandMax = res.data.azimuthal_transfer_max || [];
    const anisotropic = Boolean(res.data.has_azimuthal_aberrations);

    if (q.length < 2) return;

    const width = 640;
    const height = 320;
    const margin = { top: 20, right: 30, bottom: 44, left: 58 };
    const innerW = width - margin.left - margin.right;
    const innerH = height - margin.top - margin.bottom;

    const qMax = q[q.length - 1] || 2.0;
    const xScale = (freq) => margin.left + (freq / qMax) * innerW;
    const yScale = (val) => margin.top + innerH / 2 - (val * innerH) / 2.2;

    const chart = svg('svg', {
      class: 'hrem-ctf-chart',
      viewBox: `0 0 ${width} ${height}`,
      width: '100%',
      height: '100%',
      role: 'img',
      'aria-label': 'Contrast transfer function against spatial frequency',
    });

    // Gridlines and tick labels: a transfer plot without a y scale cannot be
    // read off, and the zero line alone does not say where +/-1 is.
    for (const level of [-1, -0.5, 0, 0.5, 1]) {
      chart.append(
        svg('line', {
          x1: margin.left,
          y1: yScale(level),
          x2: margin.left + innerW,
          y2: yScale(level),
          stroke: 'var(--border-color, #ccc)',
          'stroke-width': level === 0 ? '1.2' : '0.5',
          'stroke-dasharray': level === 0 ? '' : '2,4',
        }),
        svg('text', {
          x: margin.left - 8,
          y: yScale(level) + 4,
          'text-anchor': 'end',
          'font-size': '10',
          fill: 'currentColor',
          text: String(level),
        }),
      );
    }
    for (let step = 0; step <= 5; step += 1) {
      const freq = (qMax * step) / 5;
      chart.append(
        svg('text', {
          x: xScale(freq),
          y: margin.top + innerH + 16,
          'text-anchor': 'middle',
          'font-size': '10',
          fill: 'currentColor',
          text: formatNumber(freq, 1),
        }),
      );
    }

    const points = (values) =>
      q.map((freq, i) => `${xScale(freq)},${yScale(values[i])}`).join(' ');

    // The azimuthal band goes down first, so the curves read on top of it.
    if (anisotropic && bandMin.length === q.length && bandMax.length === q.length) {
      const upper = q.map((freq, i) => `${xScale(freq)},${yScale(bandMax[i])}`);
      const lower = q
        .map((freq, i) => `${xScale(freq)},${yScale(bandMin[i])}`)
        .reverse();
      chart.append(
        svg('polygon', {
          points: [...upper, ...lower].join(' '),
          fill: COLOR_BAND,
          'fill-opacity': '0.22',
          stroke: 'none',
        }),
      );
    }

    chart.append(
      svg('polyline', {
        points: points(undamped),
        fill: 'none',
        stroke: COLOR_UNDAMPED,
        'stroke-width': '1.2',
        'stroke-dasharray': '3,3',
      }),
      svg('polyline', {
        points: points(envelope),
        fill: 'none',
        stroke: COLOR_ENVELOPE,
        'stroke-width': '1.5',
        'stroke-dasharray': '4,4',
      }),
      svg('polyline', {
        points: q.map((freq, i) => `${xScale(freq)},${yScale(-envelope[i])}`).join(' '),
        fill: 'none',
        stroke: COLOR_ENVELOPE,
        'stroke-width': '1.5',
        'stroke-dasharray': '4,4',
      }),
      svg('polyline', {
        points: points(transfer),
        fill: 'none',
        stroke: COLOR_TRANSFER,
        'stroke-width': '2',
      }),
    );

    const marker = (value, color, label, anchor) => {
      if (!value || value <= 0 || Number.isNaN(value)) return;
      const qv = 1.0 / value;
      if (qv > qMax) return;
      chart.append(
        svg('line', {
          x1: xScale(qv),
          y1: margin.top,
          x2: xScale(qv),
          y2: margin.top + innerH,
          stroke: color,
          'stroke-width': '1.5',
          'stroke-dasharray': '2,2',
        }),
        svg('text', {
          x: anchor === 'end' ? xScale(qv) - 4 : xScale(qv) + 4,
          y: anchor === 'end' ? margin.top + 28 : margin.top + 14,
          fill: color,
          'font-size': '11',
          'text-anchor': anchor,
          text: label,
        }),
      );
    };

    marker(
      res.data.point_resolution_angstrom,
      COLOR_POINT_RES,
      `d₀ = ${formatNumber(res.data.point_resolution_angstrom, 2)} Å`,
      'start',
    );
    marker(
      res.data.information_limit_angstrom,
      COLOR_INFO_LIMIT,
      `d_info = ${formatNumber(res.data.information_limit_angstrom, 2)} Å`,
      'end',
    );

    chart.append(
      svg('text', {
        x: margin.left + innerW / 2,
        y: height - 6,
        'text-anchor': 'middle',
        'font-size': '12',
        fill: 'currentColor',
        text: 'Spatial frequency q (Å⁻¹)',
      }),
      svg('text', {
        x: 16,
        y: margin.top + innerH / 2,
        'text-anchor': 'middle',
        'font-size': '12',
        fill: 'currentColor',
        transform: `rotate(-90 16 ${margin.top + innerH / 2})`,
        text: 'Transfer T(q)',
      }),
    );

    ctfFrame.setContent(chart);
    ctfFrame.setStatus(
      anisotropic
        ? `Cut at θ = ${formatNumber(res.data.azimuth_deg, 1)}°; shaded band spans every azimuth.`
        : 'Round lens: transfer is the same at every azimuth.',
    );
  }

  function swatch(color, label, dashed = false, filled = false) {
    const mark = svg('svg', { width: '22', height: '10', 'aria-hidden': 'true' }, [
      filled
        ? svg('rect', { x: '1', y: '1', width: '20', height: '8', fill: color, 'fill-opacity': '0.35' })
        : svg('line', {
            x1: '1',
            y1: '5',
            x2: '21',
            y2: '5',
            stroke: color,
            'stroke-width': '2',
            'stroke-dasharray': dashed ? '3,3' : '',
          }),
    ]);
    return el('span.hrem-legend__item', {}, [mark, el('span', { text: label })]);
  }

  function updateLegend() {
    legend.replaceChildren();
    if (state.view.id !== 'tem.ctf_calculator') return;
    const res = state.results['tem.ctf_calculator'];
    if (!res || !res.data) return;

    const items = [
      swatch(COLOR_TRANSFER, 'Damped transfer T(q)'),
      swatch(COLOR_ENVELOPE, 'Coherence envelope ±E(q)', true),
      swatch(COLOR_UNDAMPED, 'Undamped sin χ(q)', true),
      swatch(COLOR_POINT_RES, 'First zero d₀'),
      swatch(COLOR_INFO_LIMIT, 'Information limit'),
    ];
    if (res.data.has_azimuthal_aberrations) {
      items.push(swatch(COLOR_BAND, 'Transfer over all azimuths', false, true));
    }
    legend.append(el('div.hrem-legend__row', {}, items));

    if (res.data.has_azimuthal_aberrations) {
      const best = res.data.azimuthal_point_resolution_best_angstrom;
      const worst = res.data.azimuthal_point_resolution_worst_angstrom;
      const terms = (res.data.residual_aberrations || [])
        .map(
          (term) =>
            `${term.symbol} = ${formatNumber(term.amplitude_angstrom, 1)} Å `
            + `at ${formatNumber(term.azimuth_deg, 1)}°`,
        )
        .join(', ');
      const spread = Number.isFinite(best) && Number.isFinite(worst)
        ? `Point resolution runs from ${formatNumber(best, 2)} Å to ${formatNumber(worst, 2)} Å `
          + `across azimuth, an anisotropy of `
          + `${formatNumber(res.data.resolution_anisotropy_angstrom, 2)} Å, coarsest at `
          + `${formatNumber(res.data.worst_azimuth_deg, 1)}°.`
        : 'No azimuth crosses zero within the sampled passband, so no directional point '
          + 'resolution is quoted.';
      legend.append(
        explainer(`Residual ${terms} makes transfer depend on direction. ${spread}`, {
          label: 'Anisotropic transfer',
        }),
      );
    }
  }

  function updateDetails() {
    const result = state.results[state.view.id];
    if (!result) {
      details.replaceChildren();
      return;
    }
    renderResult(details, result, { teaches: state.teaches });
  }

  applyView();

  if (examples.length > 0) {
    loadExample(examples[0]);
  }
}
