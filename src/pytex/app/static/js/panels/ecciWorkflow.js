/**
 * The ECCI workflow: an EBSD orientation, its on-axis view, and the tilt to a
 * two-beam condition.
 *
 * Three things share the screen. The EBSD Kikuchi pattern (left) is drawn the
 * same way `ebsdkikuchi.js` draws it — bands as the gap between two Kossel-cone
 * edges, zone axes as circles sized by how many bands meet there — because it
 * is the same simulation, just embedded rather than standalone. The on-axis
 * view (right) is a TEM-style zone-axis pattern, drawn with the shared
 * `core/saedplot.js` helper an ordinary SAED simulation uses, because an
 * on-axis BSE detector *is* that geometry: its normal is the beam, whatever the
 * specimen tilt. Under both, tilt and rotation controls call the live
 * re-simulation operation on every move and redraw both pictures, with the
 * current zone-axis proximity always shown beside them.
 *
 * "Solve" finds the stage move that brings a chosen direction exactly onto the
 * beam; "Go to solution" copies it into the tilt/rotation controls so the two
 * patterns can be watched converging onto the two-beam condition the solver
 * predicted, one small step at a time or in a single jump.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import { bandLabelNode, labelAngleDeg } from '../core/kikuchilabel.js';
import { plotFrame } from '../core/plotframe.js';
import { drawSimulatedPattern } from '../core/saedplot.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'ecci_workflow',
  title: 'ECCI workflow',
  tagline: 'From an EBSD orientation: the on-axis view and the tilt to a two-beam condition.',
};

const SCREEN_COLOUR = '#05070d';
const BAND_COLOUR = '#dbeafe';
const TRACE_COLOUR = '#7dd3fc';
const AXIS_COLOUR = '#fbbf24';

export function mount(context) {
  const solveOp = context.manifest.operations.find((entry) => entry.id === 'ecci.solve_workflow');
  const liveOp = context.manifest.operations.find((entry) => entry.id === 'ecci.resimulate');
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const state = {
    result: null,
    form: null,
    teaches: null,
    busy: false,
  };

  const kikuchiFrame = plotFrame({ title: 'EBSD Kikuchi pattern (current stage state)' });
  const onAxisFrame = plotFrame({ title: 'On-axis view (TEM-style, current stage state)' });

  const proximityCard = el('div.measure');
  const solutionCard = el('div.measure');

  const tiltInput = numberField('Stage tilt', 'deg', 0, 89.9, 0.05);
  const rotationInput = numberField('Stage rotation', 'deg', -180, 180, 0.1);
  tiltInput.input.addEventListener('input', () => scheduleLive());
  rotationInput.input.addEventListener('input', () => scheduleLive());

  const goToSolutionButton = el('button.button', {
    type: 'button',
    text: 'Go to solved tilt/rotation',
    disabled: true,
    onclick: () => {
      if (!state.result?.data?.solution) return;
      const solution = state.result.data.solution;
      tiltInput.input.value = String(solution.tilt_deg);
      rotationInput.input.value = String(solution.rotation_deg);
      runLive();
    },
  });

  const formHost = el('div');
  const solveButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Solve the ECCI tilt',
    onclick: () => runSolve(),
  });

  const details = el('div');

  context.rail.append(
    el('details.group', { open: true }, [
      el('summary', { text: 'Orientation, geometry, and target direction' }),
      el('div.group__body', {}, [formHost, solveButton]),
    ]),
    el('details.group', { open: true }, [
      el('summary', { text: 'Stage state (scrub to re-simulate live)' }),
      el('div.group__body', {}, [
        tiltInput.row,
        rotationInput.row,
        goToSolutionButton,
      ]),
    ]),
    el('details.group', { open: true }, [
      el('summary', { text: 'Zone-axis proximity' }),
      el('div.group__body', {}, [proximityCard]),
    ]),
    el('details.group', {}, [
      el('summary', { text: 'Solved stage move' }),
      el('div.group__body', {}, [solutionCard]),
    ]),
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
    details,
  );

  context.stage.append(
    el('div.stage__split', {}, [kikuchiFrame.element, onAxisFrame.element]),
  );

  function numberField(label, units, min, max, step) {
    const input = el('input', { type: 'number', min, max, step, value: '70' });
    const row = el('label.field', {}, [
      el('span.field__label', { text: `${label} (${units})` }),
      input,
    ]);
    return { row, input };
  }

  function renderControls(initial = {}) {
    state.form = buildForm(solveOp, { initial });
    formHost.replaceChildren(state.form.element);
    const values = state.form.values();
    if (typeof values.stage_tilt_deg === 'number') tiltInput.input.value = String(values.stage_tilt_deg);
    if (typeof values.stage_rotation_deg === 'number') {
      rotationInput.input.value = String(values.stage_rotation_deg);
    }
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    runSolve();
  }

  function currentLiveRequest() {
    const base = { ...state.form.values() };
    base.stage_tilt_deg = Number(tiltInput.input.value);
    base.stage_rotation_deg = Number(rotationInput.input.value);
    return base;
  }

  async function runSolve() {
    solveButton.disabled = true;
    solveButton.textContent = 'Solving…';
    state.form.clearErrors();
    try {
      state.result = await call(solveOp.id, state.form.values());
      const values = state.form.values();
      tiltInput.input.value = String(values.stage_tilt_deg);
      rotationInput.input.value = String(values.stage_rotation_deg);
      draw(state.result.data);
      goToSolutionButton.disabled = !state.result.data.solution;
      renderResult(details, state.result, { teaches: state.teaches });
      state.teaches = null;
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      solveButton.disabled = false;
      solveButton.textContent = 'Solve the ECCI tilt';
    }
  }

  let liveTimer = null;
  function scheduleLive() {
    if (liveTimer) window.clearTimeout(liveTimer);
    liveTimer = window.setTimeout(() => runLive(), 120);
  }

  async function runLive() {
    if (!state.form || state.busy) return;
    state.busy = true;
    try {
      const payload = await call(liveOp.id, currentLiveRequest());
      draw(payload.data);
    } catch (error) {
      context.showError(error, { quiet: true });
    } finally {
      state.busy = false;
    }
  }

  function draw(data) {
    if (!data) return;
    drawKikuchi(data.kikuchi);
    drawOnAxis(data.on_axis);
    drawProximity(data);
  }

  function drawKikuchi(kikuchi) {
    if (!kikuchi) {
      kikuchiFrame.setContent(el('div.stage__placeholder', { text: 'Solve or scrub to simulate.' }));
      return;
    }
    const width = kikuchi.width_px;
    const height = kikuchi.height_px;
    const outer = svg('svg', { viewBox: `0 0 ${width} ${height}`, preserveAspectRatio: 'xMidYMid meet' });
    const clipId = 'ecci-kikuchi-clip';
    const root = svg('g', { 'clip-path': `url(#${clipId})` });
    outer.append(
      svg('defs', {}, [svg('clipPath', { id: clipId }, [svg('rect', { x: 0, y: 0, width, height })])]),
      root,
    );
    root.append(svg('rect', { x: 0, y: 0, width, height, fill: SCREEN_COLOUR }));
    const weight = Math.max(width, height) / 900;
    const font = Math.max(width, height) / 48;
    for (const band of kikuchi.bands) {
      const prominence = 0.28 + 0.6 * Math.min(1, Number(band.intensity) || 0);
      const [narrow, far] = band.edges;
      if (narrow.length === 1 && far.length === 1) {
        const outline = [...narrow[0], ...[...far[0]].reverse()];
        root.append(
          svg('polygon', {
            points: outline.map(([x, y]) => `${x},${y}`).join(' '),
            fill: BAND_COLOUR,
            'fill-opacity': 0.1 * prominence,
            stroke: 'none',
          }),
        );
      }
      for (const edge of [narrow, far]) {
        for (const run of edge) {
          root.append(
            svg('polyline', {
              points: run.map(([x, y]) => `${x},${y}`).join(' '),
              fill: 'none',
              stroke: BAND_COLOUR,
              'stroke-opacity': prominence,
              'stroke-width': weight * 1.3,
              'stroke-linecap': 'round',
            }),
          );
        }
      }
      for (const run of band.centre) {
        const node = svg('polyline', {
          points: run.map(([x, y]) => `${x},${y}`).join(' '),
          fill: 'none',
          stroke: TRACE_COLOUR,
          'stroke-opacity': 0.45 * prominence + 0.2,
          'stroke-width': weight,
          'stroke-dasharray': `${weight * 7} ${weight * 7}`,
        });
        root.append(node);
        kikuchiFrame.hoverable(node, {
          Plane: band.label,
          'd / Å': band.d_angstrom,
          'Band width / deg': band.width_deg,
        });
      }
      if (band.centre.length && band.centre[0].length > 2) {
        const run = band.centre[0];
        const at = Math.min(run.length - 1, Math.round(run.length * 0.25));
        root.append(
          bandLabelNode({
            x: run[at][0],
            y: run[at][1],
            angleDeg: labelAngleDeg(run[Math.max(0, at - 2)], run[at]),
            text: band.label,
            fontSize: font,
            colour: TRACE_COLOUR,
            haloColour: SCREEN_COLOUR,
            haloWidth: font / 5,
          }),
        );
      }
    }
    for (const axis of kikuchi.zone_axes) {
      if (!axis.on_detector) continue;
      const radius = weight * (3 + 1.1 * Math.min(axis.order, 8));
      root.append(
        svg('circle', {
          cx: axis.x, cy: axis.y, r: radius, fill: 'none', stroke: AXIS_COLOUR,
          'stroke-opacity': 0.85, 'stroke-width': weight * 1.2,
        }),
        svg('text', {
          x: axis.x, y: axis.y - radius - font * 0.35, 'font-size': font * 0.9,
          fill: AXIS_COLOUR, stroke: SCREEN_COLOUR, 'stroke-width': font / 6,
          'paint-order': 'stroke', 'text-anchor': 'middle', text: axis.label,
        }),
      );
    }
    kikuchiFrame.setContent(outer);
    kikuchiFrame.setStatus(
      `${kikuchi.bands.length} band(s) · λ = ${formatNumber(kikuchi.wavelength_angstrom, 4)} Å`,
    );
  }

  function drawOnAxis(onAxis) {
    if (!onAxis) {
      onAxisFrame.setContent(el('div.stage__placeholder', { text: 'Solve or scrub to simulate.' }));
      return;
    }
    const outer = svg('svg', {
      viewBox: `0 0 ${onAxis.width_px} ${onAxis.height_px}`,
      preserveAspectRatio: 'xMidYMid meet',
    });
    drawSimulatedPattern(outer, onAxis, {
      labels: true,
      hoverable: (node, row) => onAxisFrame.hoverable(node, row),
    });
    onAxisFrame.setContent(outer);
    const zoneLabel = onAxis.nearest_zone_axis.join(' ');
    onAxisFrame.setStatus(
      `Nearest zone [${zoneLabel}], ${formatNumber(onAxis.nearest_zone_axis_deviation_deg, 2)}° off · `
      + `max |excitation error| ${formatNumber(onAxis.max_abs_excitation_error_inv_angstrom, 4)} Å⁻¹`,
    );
  }

  function drawProximity(data) {
    const row = (label, value) =>
      el('tr', {}, [el('th', { text: label }), el('td.measure__num', { text: value })]);
    const proximity = data.proximity;
    const target = data.target;
    proximityCard.replaceChildren(
      el('table.measure__table', {}, [
        el('tbody', {}, [
          row('Nearest zone axis', proximity ? proximity.label : '—'),
          row('Deviation', proximity ? `${formatNumber(proximity.deviation_deg, 2)}°` : '—'),
          row('Target direction', target ? target.label : '—'),
          ...(target && typeof target.angle_from_beam_deg === 'number'
            ? [row('Target off beam by', `${formatNumber(target.angle_from_beam_deg, 2)}°`)]
            : []),
        ]),
      ]),
    );
    if (data.solution) {
      const solution = data.solution;
      solutionCard.replaceChildren(
        el('table.measure__table', {}, [
          el('tbody', {}, [
            row('Solved tilt', `${formatNumber(solution.tilt_deg, 2)}°`),
            row('Solved rotation', `${formatNumber(solution.rotation_deg, 2)}°`),
            row('Move required', `Δtilt ${formatNumber(solution.delta_tilt_deg, 2)}°, Δrotation ${formatNumber(solution.delta_rotation_deg, 2)}°`),
            row('Forward-validated residual', `${formatNumber(solution.residual_deg, 4)}°`),
          ]),
        ]),
      );
    }
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);

  return { help: () => solveOp };
}
