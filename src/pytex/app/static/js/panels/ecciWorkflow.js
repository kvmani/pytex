/**
 * The ECCI workflow: an EBSD orientation, its on-axis view, and the tilt to a
 * two-beam condition.
 *
 * Four things share the screen. The EBSD Kikuchi pattern (top left) is drawn the
 * same way `ebsdkikuchi.js` draws it — bands as the gap between two Kossel-cone
 * edges, zone axes as circles sized by how many bands meet there — because it
 * is the same simulation, just embedded rather than standalone. The on-axis
 * view (right) is a TEM-style zone-axis pattern, drawn with the shared
 * `core/saedplot.js` helper an ordinary SAED simulation uses, because an
 * on-axis BSE detector *is* that geometry: its normal is the beam, whatever the
 * specimen tilt. Under both is the stage console: a schematic of the stage, a
 * tracker showing where the target sits relative to the beam and which way each
 * control moves it, and the tilt and rotation sliders themselves, which call the
 * live re-simulation operation on every move and redraw both pictures above.
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
  // Must equal the `panel=` the ECCI operations and examples are registered
  // under, as every other panel's does. It did not, so this panel's examples
  // never matched their own filter and it opened empty every time, and feature
  // search could not resolve an ECCI operation back to the panel that runs it.
  id: 'ecci',
  title: 'ECCI workflow',
  tagline: 'From an EBSD orientation: the on-axis view and the tilt to a two-beam condition.',
};

const SCREEN_COLOUR = '#05070d';
const BAND_COLOUR = '#dbeafe';
const TRACE_COLOUR = '#7dd3fc';
const AXIS_COLOUR = '#fbbf24';


/*
 * The stage console.
 *
 * The two pictures above say what the beam currently sees. They do not say what
 * to do about it, and that is the question an operator actually has: the target
 * is twelve degrees off the beam, so which way does the tilt knob move it, and
 * how far? Watching a single deviation number while nudging a control is the
 * slow way to find out, and it is what the console replaces.
 *
 * Three things, left to right. A schematic of the stage itself, so the numbers
 * in the controls have a physical picture attached. A target tracker: the beam
 * at the centre of a set of rings, the target as a dot that has to be walked
 * into the middle, and an arrow for each control showing which way it pushes
 * that dot at the current stage state. And the controls, as sliders that
 * re-simulate as they move, so the effect of an input is seen rather than
 * inferred.
 */

const STAGE_INK = '#94a3b8';
const STAGE_BEAM = '#fbbf24';
const STAGE_SPECIMEN = '#38bdf8';
const STAGE_TARGET = '#f472b6';
const STAGE_TILT_ARROW = '#a3e635';
const STAGE_ROT_ARROW = '#c084fc';

/** Degrees of arc the guide arrows are drawn for, so a 1-degree step is visible. */
const GUIDE_SCALE_DEG = 10;

/**
 * A slider and a number box that are the same value.
 *
 * The slider is for exploring — drag it and watch both patterns move — and the
 * box is for arriving, because "54.7356" cannot be reached by dragging. Nudge
 * buttons sit between the two: a known step, repeatable, which is how a real
 * stage is driven.
 */
function stageControl({ label, units, min, max, step, value, steps, onInput }) {
  const slider = el('input.ecci-control__slider', {
    type: 'range', min, max, step, value: String(value),
    'aria-label': `${label} (${units})`,
  });
  const box = el('input.ecci-control__box', {
    type: 'number', min, max, step, value: String(value),
    'aria-label': `${label} (${units}), exact value`,
  });

  let silent = false;
  const publish = (next, from) => {
    if (silent) return;
    const clamped = Math.min(max, Math.max(min, Number(next)));
    if (!Number.isFinite(clamped)) return;
    silent = true;
    if (from !== 'slider') slider.value = String(clamped);
    if (from !== 'box') box.value = String(clamped);
    silent = false;
    onInput(clamped);
  };

  slider.addEventListener('input', () => publish(slider.value, 'slider'));
  box.addEventListener('input', () => publish(box.value, 'box'));

  const nudges = el(
    'div.ecci-control__nudges',
    {},
    steps.flatMap((amount) =>
      [-amount, amount].map((delta) =>
        el('button.ecci-control__nudge', {
          type: 'button',
          text: `${delta > 0 ? '+' : '−'}${Math.abs(delta)}`,
          title: `${delta > 0 ? 'Increase' : 'Decrease'} ${label.toLowerCase()} by ${Math.abs(delta)}°`,
          onclick: () => publish(Number(box.value) + delta, null),
        }),
      ),
    ),
  );

  const element = el('div.ecci-control', {}, [
    el('div.ecci-control__head', {}, [
      el('span.ecci-control__label', { text: label }),
      box,
      el('span.ecci-control__units', { text: units }),
    ]),
    slider,
    nudges,
  ]);

  return {
    element,
    get value() {
      return Number(box.value);
    },
    set(next) {
      silent = true;
      const clamped = Math.min(max, Math.max(min, Number(next)));
      slider.value = String(clamped);
      box.value = String(clamped);
      silent = false;
    },
  };
}

/**
 * The stage, drawn side-on, with a plan view of the rotation beside it.
 *
 * Deliberately a schematic and not a rendering: what has to be legible is the
 * tilt angle between the beam and the specimen normal, and the fact that
 * rotation happens about that normal and therefore does nothing at all when the
 * target already lies on it. A photorealistic holder would communicate neither.
 */
function drawStageSchematic(root, { tiltDeg, rotationDeg }) {
  root.replaceChildren();
  const width = 220;
  const height = 150;
  root.setAttribute('viewBox', `0 0 ${width} ${height}`);

  const pivotX = 78;
  const pivotY = 96;
  const halfLength = 46;
  // The stage tilts about the laboratory x axis, which is perpendicular to this
  // elevation, so a tilt appears here as the specimen bar swinging away from
  // horizontal by exactly the tilt angle.
  const radians = (tiltDeg * Math.PI) / 180;
  const dx = Math.cos(radians) * halfLength;
  const dy = -Math.sin(radians) * halfLength;

  root.append(
    // The beam, always vertical: it is the fixed thing, and the stage moves.
    svg('line', {
      x1: pivotX, y1: 12, x2: pivotX, y2: pivotY - 4,
      stroke: STAGE_BEAM, 'stroke-width': 2.4, 'stroke-linecap': 'round',
    }),
    svg('polygon', {
      points: `${pivotX - 4},${pivotY - 10} ${pivotX + 4},${pivotY - 10} ${pivotX},${pivotY - 2}`,
      fill: STAGE_BEAM,
    }),
    svg('text', {
      x: pivotX + 7, y: 20, 'font-size': 8, fill: STAGE_BEAM, text: 'beam',
    }),
    // The horizontal, so the tilt angle has something to be measured from.
    svg('line', {
      x1: pivotX - halfLength, y1: pivotY, x2: pivotX + halfLength, y2: pivotY,
      stroke: STAGE_INK, 'stroke-width': 0.8, 'stroke-dasharray': '3 3', 'stroke-opacity': 0.6,
    }),
    // The specimen.
    svg('line', {
      x1: pivotX - dx, y1: pivotY - dy, x2: pivotX + dx, y2: pivotY + dy,
      stroke: STAGE_SPECIMEN, 'stroke-width': 5, 'stroke-linecap': 'round',
    }),
    // The specimen normal, which is what rotation turns about.
    svg('line', {
      x1: pivotX, y1: pivotY,
      x2: pivotX + Math.sin(radians) * 34, y2: pivotY - Math.cos(radians) * 34,
      stroke: STAGE_SPECIMEN, 'stroke-width': 1.2, 'stroke-dasharray': '4 3',
    }),
    svg('circle', { cx: pivotX, cy: pivotY, r: 2.6, fill: STAGE_INK }),
    svg('text', {
      x: 6, y: pivotY + 24, 'font-size': 8.5, fill: STAGE_INK,
      text: `tilt ${tiltDeg.toFixed(1)}°`,
    }),
  );

  // The tilt angle itself, as an arc from vertical to the specimen normal.
  const arcRadius = 22;
  const endX = pivotX + Math.sin(radians) * arcRadius;
  const endY = pivotY - Math.cos(radians) * arcRadius;
  root.append(
    svg('path', {
      d: `M ${pivotX} ${pivotY - arcRadius} A ${arcRadius} ${arcRadius} 0 0 1 ${endX} ${endY}`,
      fill: 'none', stroke: STAGE_INK, 'stroke-width': 1,
    }),
  );

  // The plan view: rotation about the specimen normal, seen down that normal.
  const planX = 178;
  const planY = 96;
  const planR = 26;
  const rotRad = (rotationDeg * Math.PI) / 180;
  root.append(
    svg('circle', {
      cx: planX, cy: planY, r: planR, fill: 'none',
      stroke: STAGE_INK, 'stroke-width': 1, 'stroke-opacity': 0.7,
    }),
    svg('line', {
      x1: planX, y1: planY, x2: planX, y2: planY - planR,
      stroke: STAGE_INK, 'stroke-width': 0.8, 'stroke-dasharray': '3 3', 'stroke-opacity': 0.6,
    }),
    svg('line', {
      x1: planX, y1: planY,
      x2: planX + Math.sin(rotRad) * planR, y2: planY - Math.cos(rotRad) * planR,
      stroke: STAGE_SPECIMEN, 'stroke-width': 2.2, 'stroke-linecap': 'round',
    }),
    svg('circle', { cx: planX, cy: planY, r: 2, fill: STAGE_INK }),
    svg('text', {
      x: planX, y: planY + planR + 14, 'font-size': 8.5, fill: STAGE_INK,
      'text-anchor': 'middle', text: `rotation ${rotationDeg.toFixed(1)}°`,
    }),
    svg('text', {
      x: planX, y: 20, 'font-size': 8, fill: STAGE_INK,
      'text-anchor': 'middle', text: 'seen down the normal',
    }),
  );
}

/**
 * Where the target is relative to the beam, and which way each control moves it.
 *
 * The beam is the centre. Rings are drawn every fifteen degrees away from it,
 * so the distance from the middle *is* the deviation the panel reports in
 * words. The target is a dot on that map, and reaching the two-beam condition
 * is exactly the act of walking the dot into the bullseye — which is a thing a
 * person can aim at, where "reduce 12.4 to zero" is not.
 *
 * The two arrows are what makes the controls legible. Each is the displacement
 * the target actually undergoes for a positive move of that control at this
 * stage state, magnified so a degree is visible. Their lengths are not equal
 * and their directions are not fixed: near the pole a degree of rotation moves
 * the target almost nowhere, and the arrow shrinks to show it, which is the
 * single most useful thing to know when a rotation control seems to be doing
 * nothing.
 */
function drawTargetTracker(root, stageView, { angleDeg, targetLabel, zoneLabel }) {
  root.replaceChildren();
  const size = 200;
  const centre = size / 2;
  const radius = 82;
  root.setAttribute('viewBox', `0 0 ${size} ${size}`);

  // Directions are mapped by the sine of their angle from the beam, so the
  // outermost ring is 90 degrees away and the map is the near hemisphere seen
  // down the beam. Equal-angle rings therefore crowd towards the edge exactly
  // as they do on a real stereogram.
  const project = (x, y) => [centre + x * radius, centre - y * radius];

  for (const ring of [15, 30, 45, 60, 75, 90]) {
    const r = Math.sin((ring * Math.PI) / 180) * radius;
    root.append(
      svg('circle', {
        cx: centre, cy: centre, r, fill: 'none', stroke: STAGE_INK,
        'stroke-width': ring === 90 ? 1.2 : 0.6,
        'stroke-opacity': ring === 90 ? 0.8 : 0.35,
      }),
    );
    if (ring < 90) {
      root.append(
        svg('text', {
          x: centre + 2, y: centre - r - 1.5, 'font-size': 6,
          fill: STAGE_INK, 'fill-opacity': 0.75, text: `${ring}°`,
        }),
      );
    }
  }

  // The tilt axis, so the schematic beside this one has a shared reference: the
  // stage tilts about laboratory x, and the target swings perpendicular to it.
  root.append(
    svg('line', {
      x1: centre - radius, y1: centre, x2: centre + radius, y2: centre,
      stroke: STAGE_INK, 'stroke-width': 0.7, 'stroke-dasharray': '5 4', 'stroke-opacity': 0.5,
    }),
    svg('text', {
      x: centre + radius - 2, y: centre - 4, 'font-size': 6, 'text-anchor': 'end',
      fill: STAGE_INK, 'fill-opacity': 0.8, text: 'tilt axis',
    }),
    // The bullseye: the beam, and the condition being aimed for.
    svg('circle', {
      cx: centre, cy: centre, r: 5.5, fill: 'none',
      stroke: STAGE_BEAM, 'stroke-width': 1.1, 'stroke-opacity': 0.9,
    }),
    svg('circle', { cx: centre, cy: centre, r: 1.6, fill: STAGE_BEAM }),
  );

  if (!stageView) {
    root.append(
      svg('text', {
        x: centre, y: size - 8, 'font-size': 8, 'text-anchor': 'middle',
        fill: STAGE_INK, text: 'Solve or move the stage.',
      }),
    );
    return;
  }

  const [tx, ty] = stageView.target_lab;
  const [px, py] = project(tx, ty);

  const arrow = (delta, colour, label) => {
    const [dx, dy] = delta;
    const scale = GUIDE_SCALE_DEG / (stageView.probe_step_deg || 1);
    const [ax, ay] = project(tx + dx * scale, ty + dy * scale);
    const length = Math.hypot(ax - px, ay - py);
    // Below about a pixel the arrow is not a short arrow, it is a smudge; a
    // stated "no effect here" is more useful than a mark too small to read.
    if (length < 1.5) return [];
    const angle = Math.atan2(ay - py, ax - px);
    const head = 4.5;
    return [
      svg('line', {
        x1: px, y1: py, x2: ax, y2: ay, stroke: colour,
        'stroke-width': 1.6, 'stroke-linecap': 'round',
      }),
      svg('polygon', {
        points: [
          `${ax},${ay}`,
          `${ax - head * Math.cos(angle - 0.45)},${ay - head * Math.sin(angle - 0.45)}`,
          `${ax - head * Math.cos(angle + 0.45)},${ay - head * Math.sin(angle + 0.45)}`,
        ].join(' '),
        fill: colour,
      }),
      svg('text', {
        x: ax + 3 * Math.cos(angle), y: ay + 3 * Math.sin(angle) - 2,
        'font-size': 6.5, fill: colour, text: label,
      }),
    ];
  };

  root.append(
    ...arrow(stageView.per_tilt_degree, STAGE_TILT_ARROW, `tilt +${GUIDE_SCALE_DEG}°`),
    ...arrow(stageView.per_rotation_degree, STAGE_ROT_ARROW, `rot +${GUIDE_SCALE_DEG}°`),
    svg('circle', {
      cx: px, cy: py, r: 4.5, fill: STAGE_TARGET,
      stroke: '#0b1220', 'stroke-width': 1,
    }),
    svg('text', {
      x: px + 7, y: py + 3, 'font-size': 7.5, fill: STAGE_TARGET,
      text: targetLabel || 'target',
    }),
    svg('text', {
      x: centre, y: size - 6, 'font-size': 7.5, 'text-anchor': 'middle', fill: STAGE_INK,
      text: zoneLabel ? `beam is on ${zoneLabel}` : '',
    }),
  );

  if (angleDeg !== null && angleDeg !== undefined && angleDeg < 0.05) {
    root.append(
      svg('text', {
        x: centre, y: 16, 'font-size': 9, 'text-anchor': 'middle',
        fill: STAGE_BEAM, 'font-weight': '700', text: 'on axis',
      }),
    );
  }
}

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

  // The stage controls live in the console under the pictures, not in the rail,
  // because their whole value is watching the pictures change as they move and
  // a control on the far side of the screen from its effect is not a live
  // control. One pair of them exists, so there is nothing to keep in step.
  const tiltControl = stageControl({
    label: 'Stage tilt', units: '°', min: 0, max: 89.9, step: 0.05,
    value: 70, steps: [0.5, 5], onInput: () => scheduleLive(),
  });
  const rotationControl = stageControl({
    label: 'Stage rotation', units: '°', min: -180, max: 180, step: 0.1,
    value: 0, steps: [1, 15], onInput: () => scheduleLive(),
  });

  // `svg()` takes attributes literally and has no `.class` shorthand, unlike
  // `el()`; the class has to be an attribute.
  const schematicSvg = svg('svg', { class: 'ecci-console__svg', viewBox: '0 0 220 150' });
  const trackerSvg = svg('svg', { class: 'ecci-console__svg', viewBox: '0 0 200 200' });
  const deviationReadout = el('p.ecci-console__deviation', {
    text: 'Solve, or move the stage, to see where the target sits.',
  });
  const deviationBar = el('div.ecci-console__bar-fill');

  const goToSolutionButton = el('button.button', {
    type: 'button',
    text: 'Go to solved tilt/rotation',
    disabled: true,
    onclick: () => {
      if (!state.result?.data?.solution) return;
      const solution = state.result.data.solution;
      tiltControl.set(solution.tilt_deg);
      rotationControl.set(solution.rotation_deg);
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
    el('div.ecci', {}, [
      el('div.ecci__views', {}, [kikuchiFrame.element, onAxisFrame.element]),
      el('section.plot.ecci-console', {}, [
        el('header.plot__header', {}, [
          el('h2.plot__title', { text: 'Stage — move it and watch both views' }),
        ]),
        el('div.ecci-console__body', {}, [
          el('figure.ecci-console__figure', {}, [
            schematicSvg,
            el('figcaption', { text: 'The stage, side-on, with its rotation seen down the normal' }),
          ]),
          el('figure.ecci-console__figure', {}, [
            trackerSvg,
            el('figcaption', {
              text: 'The target relative to the beam. Arrows show which way each control moves it.',
            }),
          ]),
          el('div.ecci-console__controls', {}, [
            deviationReadout,
            el('div.ecci-console__bar', {}, [deviationBar]),
            tiltControl.element,
            rotationControl.element,
            goToSolutionButton,
          ]),
        ]),
      ]),
    ]),
  );

  function renderControls(initial = {}) {
    state.form = buildForm(solveOp, { initial });
    formHost.replaceChildren(state.form.element);
    const values = state.form.values();
    if (typeof values.stage_tilt_deg === 'number') tiltControl.set(values.stage_tilt_deg);
    if (typeof values.stage_rotation_deg === 'number') rotationControl.set(values.stage_rotation_deg);
    // The console under the pictures *is* the stage state, so the rail must not
    // offer a second pair of boxes for it: they would sit there holding the
    // value the panel opened with while the console showed the value in force,
    // and there would be no way to tell which one Solve was about to use.
    const hiddenGroups = new Set();
    for (const name of ['stage_tilt_deg', 'stage_rotation_deg']) {
      const field = state.form.field(name);
      if (!field) continue;
      field.element.hidden = true;
      const group = field.element.closest('details');
      if (group) hiddenGroups.add(group);
    }
    // A group whose every control has moved to the console is a heading with
    // nothing under it, which reads as a section that failed to load.
    for (const group of hiddenGroups) {
      const remaining = Array.from(group.querySelectorAll('.field')).some((node) => !node.hidden);
      if (!remaining) group.hidden = true;
    }
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    runSolve();
  }

  // The two operations do not take the same parameters: solving accepts
  // `allow_reverse`, which chooses among stage branches, and re-simulating has
  // no branches to choose between and rejects it. Copying the solve form
  // wholesale therefore made every live call fail with an unknown-parameter
  // error, silently, which is why moving the controls appeared to do nothing.
  // Filtering against the operation's own declared parameters fixes that and
  // cannot drift: a parameter added to one and not the other is handled by
  // construction rather than by remembering.
  const liveParameterNames = new Set((liveOp.parameters ?? []).map((entry) => entry.name));

  function currentSolveRequest() {
    return {
      ...state.form.values(),
      stage_tilt_deg: tiltControl.value,
      stage_rotation_deg: rotationControl.value,
    };
  }

  function currentLiveRequest() {
    const values = currentSolveRequest();
    const request = {};
    for (const [key, value] of Object.entries(values)) {
      if (liveParameterNames.has(key)) request[key] = value;
    }
    return request;
  }

  async function runSolve() {
    solveButton.disabled = true;
    solveButton.textContent = 'Solving…';
    state.form.clearErrors();
    try {
      // Solving is always *from where the stage is now*, and where it is now
      // is what the console says.
      state.result = await call(solveOp.id, currentSolveRequest());
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
      // Quiet, because this fires on every slider move and a modal per frame is
      // unusable — but never invisible: the readout the user is already
      // watching says the stage stopped following them.
      context.showError(error, { quiet: true });
      deviationReadout.textContent = `Could not re-simulate: ${error.message ?? error}`;
      deviationReadout.classList.add('is-error');
    } finally {
      state.busy = false;
    }
  }

  function draw(data) {
    if (!data) return;
    drawKikuchi(data.kikuchi);
    drawOnAxis(data.on_axis);
    drawProximity(data);
    drawConsole(data);
  }

  function drawConsole(data) {
    drawStageSchematic(schematicSvg, {
      tiltDeg: tiltControl.value,
      rotationDeg: rotationControl.value,
    });
    deviationReadout.classList.remove('is-error');
    const angle = data.target ? data.target.angle_from_beam_deg : null;
    drawTargetTracker(trackerSvg, data.stage_view ?? null, {
      angleDeg: angle,
      targetLabel: data.target ? data.target.label : '',
      zoneLabel: data.proximity ? data.proximity.label : '',
    });
    if (typeof angle !== 'number') {
      deviationReadout.textContent = 'Solve, or move the stage, to see where the target sits.';
      deviationBar.style.width = '0%';
      return;
    }
    const label = data.target.label;
    deviationReadout.textContent =
      angle < 0.05
        ? `${label} is on the beam — this is the two-beam condition.`
        : `${label} is ${formatNumber(angle, 2)}° off the beam.`;
    deviationReadout.classList.toggle('is-on-axis', angle < 0.05);
    // Full bar at 90 degrees off, empty on axis: the bar shrinking as the dot
    // walks into the bullseye is the same fact told twice, which is what makes
    // a slider worth dragging.
    deviationBar.style.width = `${Math.min(100, (angle / 90) * 100)}%`;
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
