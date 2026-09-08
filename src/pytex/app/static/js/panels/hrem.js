/**
 * HRTEM Simulation and Objective Lens Contrast Transfer Function panel.
 *
 * Implements interactive multislice high-resolution TEM simulations and optical
 * CTF diagnostics for double-corrected transmission electron microscopy.
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

export function mount(context) {
  const simOperation = context.manifest.operations.find(
    (entry) => entry.id === 'tem.simulate_hrem',
  );
  const ctfOperation = context.manifest.operations.find(
    (entry) => entry.id === 'tem.ctf_calculator',
  );
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const state = {
    mode: 'simulation', // 'simulation' | 'ctf'
    simResult: null,
    ctfResult: null,
    simForm: null,
    ctfForm: null,
    teaches: null,
  };

  // Stage frames
  const simFrame = plotFrame({
    title: 'Simulated Micrograph',
    units: 'Å',
    toolbar: [],
  });

  const fftFrame = plotFrame({
    title: 'Power Spectrum (Thon Rings)',
    units: 'Å⁻¹',
    toolbar: [],
  });

  const ctfFrame = plotFrame({
    title: 'Contrast Transfer Function & Envelopes',
    units: 'Å⁻¹',
    toolbar: [],
  });

  const stageWrapper = el('div.hrem-stage');
  const simStage = el('div.hrem-sim-stage', {}, [simFrame.element, fftFrame.element]);
  const ctfStage = el('div.hrem-ctf-stage', {}, [ctfFrame.element]);
  const details = el('div.hrem-details');

  stageWrapper.append(simStage, ctfStage, details);
  context.stage.append(stageWrapper);

  // Rail controls
  const modeSwitchSim = el('button.button', {
    type: 'button',
    text: 'HRTEM Simulation',
    'aria-pressed': 'true',
    onclick: () => setMode('simulation'),
  });

  const modeSwitchCtf = el('button.button', {
    type: 'button',
    text: 'CTF & Aberrations',
    'aria-pressed': 'false',
    onclick: () => setMode('ctf'),
  });

  const modeSelector = el('div.button-group.button-group--block', {}, [
    modeSwitchSim,
    modeSwitchCtf,
  ]);

  const simFormHost = el('div.hrem-sim-form');
  const ctfFormHost = el('div.hrem-ctf-form');

  if (simOperation) {
    state.simForm = buildForm(simOperation, {
      values: {},
      onChange: () => {},
    });
    simFormHost.append(state.simForm.element);
  }

  if (ctfOperation) {
    state.ctfForm = buildForm(ctfOperation, {
      values: {},
      onChange: () => {},
    });
    ctfFormHost.append(state.ctfForm.element);
  }

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Run simulation',
    onclick: () => runCurrent(),
  });

  // Example scenario selector
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
      explainer('Preset scenarios demonstrate high-resolution TEM contrast and CTF calibration.', {
        label: 'HRTEM Presets',
      }),
      el('div.stack.stack--tight', {}, exampleItems),
    ]),
  ]);

  context.rail.append(modeSelector, simFormHost, ctfFormHost, runButton, exampleGroup);

  function setMode(mode) {
    state.mode = mode;
    const isSim = mode === 'simulation';
    modeSwitchSim.setAttribute('aria-pressed', String(isSim));
    modeSwitchCtf.setAttribute('aria-pressed', String(!isSim));
    simStage.style.display = isSim ? 'block' : 'none';
    ctfStage.style.display = isSim ? 'none' : 'block';
    simFormHost.style.display = isSim ? 'block' : 'none';
    ctfFormHost.style.display = isSim ? 'none' : 'block';
    runButton.textContent = isSim ? 'Run HRTEM simulation' : 'Calculate CTF';
    updateDetails();
  }

  function runCurrent() {
    if (state.mode === 'simulation' && state.simForm) {
      const values = state.simForm.read();
      context.spin(runButton, async () => {
        try {
          const result = await call(simOperation.id, values);
          state.simResult = result;
          state.teaches = null;
          drawSimulation();
          updateDetails();
        } catch (err) {
          context.reportError(err);
        }
      });
    } else if (state.mode === 'ctf' && state.ctfForm) {
      const values = state.ctfForm.read();
      context.spin(runButton, async () => {
        try {
          const result = await call(ctfOperation.id, values);
          state.ctfResult = result;
          state.teaches = null;
          drawCTF();
          updateDetails();
        } catch (err) {
          context.reportError(err);
        }
      });
    }
  }

  function loadExample(ex) {
    if (ex.operation === 'tem.simulate_hrem') {
      setMode('simulation');
      if (state.simForm) state.simForm.fill(ex.request);
    } else if (ex.operation === 'tem.ctf_calculator') {
      setMode('ctf');
      if (state.ctfForm) state.ctfForm.fill(ex.request);
    }
    state.teaches = ex.teaches;
    runCurrent();
  }

  function drawSimulation() {
    const res = state.simResult;
    if (!res || !res.data) return;

    // Draw micrograph
    const imgData = res.data.image_png;
    if (imgData) {
      const img = el('img', {
        src: imgData,
        alt: 'Simulated HRTEM micrograph',
        style: 'max-width: 100%; height: auto; border-radius: 4px; display: block; margin: 0 auto;',
      });
      const extent = res.data.extent_angstrom || [10, 10];
      const caption = el('div.text-center.text-muted', {
        text: `Field of view: ${formatNumber(extent[0], 1)} × ${formatNumber(extent[1], 1)} Å (pixel size: ${formatNumber(res.data.pixel_size_angstrom, 3)} Å/px)`,
      });
      simFrame.setContent(el('div.stack', {}, [img, caption]));
    }

    // Draw FFT power spectrum
    const psData = res.data.power_spectrum_png;
    if (psData) {
      const psImg = el('img', {
        src: psData,
        alt: '2D FFT Power Spectrum (Thon Rings)',
        style: 'max-width: 100%; height: auto; border-radius: 4px; display: block; margin: 0 auto;',
      });
      const psCaption = el('div.text-center.text-muted', {
        text: `Point resolution: ${formatNumber(res.data.point_resolution_angstrom, 2)} Å | Information limit: ${formatNumber(res.data.information_limit_angstrom, 2)} Å`,
      });
      fftFrame.setContent(el('div.stack', {}, [psImg, psCaption]));
    }
  }

  function drawCTF() {
    const res = state.ctfResult;
    if (!res || !res.data) return;

    const q = res.data.spatial_frequencies || [];
    const transfer = res.data.transfer_function || [];
    const undamped = res.data.undamped_ctf || [];
    const envelope = res.data.total_envelope || [];

    if (q.length < 2) return;

    const width = 600;
    const height = 300;
    const margin = { top: 20, right: 30, bottom: 40, left: 50 };
    const innerW = width - margin.left - margin.right;
    const innerH = height - margin.top - margin.bottom;

    const qMax = q[q.length - 1] || 2.0;
    const xScale = (freq) => margin.left + (freq / qMax) * innerW;
    const yScale = (val) => margin.top + innerH / 2 - (val * innerH) / 2.2;

    const chart = svg('svg', {
      viewBox: `0 0 ${width} ${height}`,
      width: '100%',
      height: '100%',
      style: 'background: var(--surface-bg, #ffffff); border-radius: 4px;',
    });

    // Zero axis
    chart.append(
      svg('line', {
        x1: margin.left,
        y1: yScale(0),
        x2: margin.left + innerW,
        y2: yScale(0),
        stroke: 'var(--border-color, #ccc)',
        'stroke-width': '1',
      }),
    );

    // Build polyline points
    const transferPoints = q.map((freq, i) => `${xScale(freq)},${yScale(transfer[i])}`).join(' ');
    const undampedPoints = q.map((freq, i) => `${xScale(freq)},${yScale(undamped[i])}`).join(' ');
    const envUpper = q.map((freq, i) => `${xScale(freq)},${yScale(envelope[i])}`).join(' ');
    const envLower = q.map((freq, i) => `${xScale(freq)},${yScale(-envelope[i])}`).join(' ');

    // Undamped curve (dashed gray)
    chart.append(
      svg('polyline', {
        points: undampedPoints,
        fill: 'none',
        stroke: '#94a3b8',
        'stroke-width': '1.2',
        'stroke-dasharray': '3,3',
      }),
    );

    // Envelope curves (amber dashed)
    chart.append(
      svg('polyline', {
        points: envUpper,
        fill: 'none',
        stroke: '#f59e0b',
        'stroke-width': '1.5',
        'stroke-dasharray': '4,4',
      }),
      svg('polyline', {
        points: envLower,
        fill: 'none',
        stroke: '#f59e0b',
        'stroke-width': '1.5',
        'stroke-dasharray': '4,4',
      }),
    );

    // Damped transfer function (solid blue)
    chart.append(
      svg('polyline', {
        points: transferPoints,
        fill: 'none',
        stroke: '#2563eb',
        'stroke-width': '2',
      }),
    );

    // Resolution markers
    const ptRes = res.data.point_resolution_angstrom;
    if (ptRes && ptRes > 0 && !isNaN(ptRes)) {
      const qZero = 1.0 / ptRes;
      if (qZero <= qMax) {
        chart.append(
          svg('line', {
            x1: xScale(qZero),
            y1: margin.top,
            x2: xScale(qZero),
            y2: margin.top + innerH,
            stroke: '#10b981',
            'stroke-width': '1.5',
            'stroke-dasharray': '2,2',
          }),
          svg('text', {
            x: xScale(qZero) + 4,
            y: margin.top + 14,
            fill: '#10b981',
            'font-size': '11',
            text: `d₀ = ${formatNumber(ptRes, 2)} Å`,
          }),
        );
      }
    }

    const infoLim = res.data.information_limit_angstrom;
    if (infoLim && infoLim > 0 && !isNaN(infoLim)) {
      const qInfo = 1.0 / infoLim;
      if (qInfo <= qMax) {
        chart.append(
          svg('line', {
            x1: xScale(qInfo),
            y1: margin.top,
            x2: xScale(qInfo),
            y2: margin.top + innerH,
            stroke: '#ef4444',
            'stroke-width': '1.5',
            'stroke-dasharray': '2,2',
          }),
          svg('text', {
            x: xScale(qInfo) - 4,
            y: margin.top + 28,
            fill: '#ef4444',
            'font-size': '11',
            'text-anchor': 'end',
            text: `d_info = ${formatNumber(infoLim, 2)} Å`,
          }),
        );
      }
    }

    // Axes and labels
    chart.append(
      svg('text', {
        x: margin.left + innerW / 2,
        y: height - 5,
        'text-anchor': 'middle',
        'font-size': '12',
        fill: 'currentColor',
        text: 'Spatial Frequency q (Å⁻¹)',
      }),
      svg('text', {
        x: 15,
        y: margin.top + innerH / 2,
        'text-anchor': 'middle',
        'font-size': '12',
        fill: 'currentColor',
        transform: `rotate(-90 15 ${margin.top + innerH / 2})`,
        text: 'Transfer T(q)',
      }),
    );

    ctfFrame.setContent(chart);
  }

  function updateDetails() {
    details.replaceChildren();
    const result = state.mode === 'simulation' ? state.simResult : state.ctfResult;
    if (result) {
      if (state.teaches) {
        details.append(explainer(state.teaches, { label: 'Lesson' }));
      }
      details.append(renderResult(result, context));
    }
  }

  // Set initial mode
  setMode('simulation');

  // Trigger initial example if available
  if (examples.length > 0) {
    loadExample(examples[0]);
  }
}
