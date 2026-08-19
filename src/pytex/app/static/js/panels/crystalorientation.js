/**
 * The crystal viewer's orientation dock: a pole figure, an inverse pole figure,
 * and the Euler angles of the view, all describing the structure beside them.
 *
 * Why the figures are drawn here rather than fetched
 * --------------------------------------------------
 * A pole figure that only updates when the mouse is released is a pole figure
 * of the *last* view, and the one thing this dock exists to show — that turning
 * the crystal turns its poles — is exactly what that would hide. So the whole
 * crystallographic content arrives once, with the scene, in the Cartesian
 * crystal frame: the symmetry operators, the expanded pole families, and the
 * fundamental sector with its inward edge normals. What happens on each frame
 * is a matrix product and a stereographic divide, which is the same contract
 * the structure itself is drawn under.
 *
 * The numbers are a different matter. Euler angles depend on a convention, and
 * a convention implemented twice is a convention that will eventually disagree
 * with itself, so nothing here decodes a rotation into angles: a debounced call
 * to `crystal.orientation` answers that, and the readout settles a moment after
 * the drag stops while the pictures stay live throughout.
 *
 * The frame
 * ---------
 * The screen is the specimen frame. RD is screen right, TD is screen up, and ND
 * points out of the screen towards the viewer — a right-handed triad, and the
 * one the texture panel already draws its pole figures in. Pole figures in the
 * literature are usually drawn with RD at the top instead; the ninety-degree
 * difference is deliberate, because the point of this dock is that the pole
 * figure is the *same view* of the *same crystal* as the structure next to it.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { applyMatrix, applyTranspose } from '../core/rotation3.js';
import { call } from '../core/api.js';

/** Half-width of each figure's drawing area, in viewBox units. */
const VIEW = 100;

/** How long the drag must be still before the numbers are asked for, in ms. */
const READOUT_DELAY = 220;

/** Positions kept in a fly-by trail. About four seconds of continuous dragging. */
const TRAIL_LENGTH = 90;

/** Angular step below which a trail position is treated as a repeat, in radians. */
const TRAIL_STEP = 0.012;

/* Theme tokens rather than literal colours, so both figures follow the app into
 * dark mode along with everything else. */
const FAMILY_COLORS = ['var(--accent)', 'var(--teal)', 'var(--violet)', 'var(--warn)'];
const AXIS_COLORS = { rd: 'var(--accent)', td: 'var(--teal)', nd: 'var(--violet)' };

/**
 * Build the dock.
 *
 * @param {object} options
 * @param {() => number[]} options.camera - The viewer's current rotation, row-major.
 * @param {(matrix: number[]) => void} options.setCamera - Turn the structure to a rotation.
 * @param {() => object} options.request - The phase and overlay inputs of the current scene.
 * @param {(error: unknown) => void} options.showError - The panel's error channel.
 * @returns {{element: HTMLElement, controls: HTMLElement, setScene: Function,
 *   draw: Function, resetTrail: Function, scheduleReadout: Function}}
 */
export function orientationDock({ camera, setCamera, request, showError }) {
  const state = {
    orientation: null,
    families: new Set(),
    axes: new Set(['nd']),
    trails: { rd: [], td: [], nd: [] },
    showTrail: true,
    showOverlays: true,
    convention: 'bunge',
    fit: null,
    readoutTimer: null,
    readoutToken: 0,
    editing: false,
  };

  const poleFigure = el('div.orient__canvas');
  const inverseFigure = el('div.orient__canvas');
  const familyChips = el('div.orient__chips');
  const axisChips = el('div.orient__chips');
  const presetChips = el('div.orient__chips');
  const angleFields = el('div.orient__angles');
  const readout = el('p.orient__readout', { text: 'Build a structure to read its orientation.' });
  const conventionSelect = el(
    'select',
    {
      'aria-label': 'Euler-angle convention',
      onchange: (event) => {
        state.convention = event.currentTarget.value;
        renderAngleFields();
        scheduleReadout({ immediate: true });
      },
    },
    [],
  );

  const angleInputs = [];

  const trailToggle = el('label.orient__toggle', {}, [
    el('input', {
      type: 'checkbox',
      checked: true,
      onchange: (event) => {
        state.showTrail = event.currentTarget.checked;
        if (!state.showTrail) resetTrail();
        draw();
      },
    }),
    el('span', { text: 'Fly-by trail' }),
  ]);

  const overlayToggle = el('label.orient__toggle', {}, [
    el('input', {
      type: 'checkbox',
      checked: true,
      onchange: (event) => {
        state.showOverlays = event.currentTarget.checked;
        draw();
      },
    }),
    el('span', { text: "The scene's own planes" }),
  ]);

  const setButton = el('button.button.button--primary', {
    type: 'button',
    text: 'Set view',
    title: 'Turn the structure to these Euler angles',
    onclick: () => applyAngles(),
  });

  /*
   * Open on a desktop, closed on a phone. On a narrow window the cards stack
   * below the structure, which puts an already off-screen figure behind a long
   * scroll; starting closed there costs one tap and keeps the structure the
   * whole of what the panel opens with.
   */
  const wide = window.matchMedia?.('(min-width: 48rem)')?.matches ?? true;

  /*
   * Two figures in the stage, the controls in the rail.
   *
   * Splitting them is not a compromise for space; it is where each belongs. The
   * stage of every panel in this application holds pictures and the rail holds
   * the controls that change them, and an orientation entered as three numbers
   * is a control. Keeping the pair together in the stage cost the figures a
   * third of the room and put a form somewhere a reader has learnt not to look
   * for one.
   */
  const element = el('details.orient', { open: wide }, [
    el('summary', { text: 'Orientation figures' }),
    el('div.orient__grid', {}, [
      el('figure.orient__cell', {}, [
        el('figcaption.orient__caption', {}, [
          el('span.orient__title', { text: 'Pole figure' }),
          familyChips,
        ]),
        poleFigure,
        overlayToggle,
      ]),
      el('figure.orient__cell', {}, [
        el('figcaption.orient__caption', {}, [
          el('span.orient__title', { text: 'Inverse pole figure' }),
          axisChips,
        ]),
        inverseFigure,
        trailToggle,
      ]),
    ]),
  ]);

  const controls = el('details.group.orient-controls', { open: true }, [
    el('summary', { text: 'Orientation' }),
    el('div.group__body', {}, [
      el('p.field__help', {
        text:
          'The screen is the specimen frame: RD points right, TD up, and ND out of the screen. '
          + 'Turning the structure and setting these angles are the same act, so the numbers '
          + 'follow the drag and the drag follows the numbers.',
      }),
      el('label.field', {}, [
        el('span.field__label', { text: 'Euler convention' }),
        conventionSelect,
      ]),
      angleFields,
      el('div.orient__actions', {}, [setButton]),
      el('span.field__hint', { text: 'Named ideal orientations' }),
      presetChips,
      readout,
    ]),
  ]);

  /* ------------------------------------------------------------- lifecycle */

  function setScene(scene) {
    state.orientation = scene?.orientation ?? null;
    state.fit = null;
    resetTrail();
    if (!state.orientation) {
      poleFigure.replaceChildren();
      inverseFigure.replaceChildren();
      return;
    }
    const families = state.orientation.pole_families ?? [];
    // Start on the last family of the list — {111} for a cubic phase, the
    // pyramidal set for hexagonal. It is the one a reader recognises fastest,
    // and one family at a time is the readable default; the others are a chip
    // away.
    state.families = new Set(families.length ? [families[families.length - 1].key] : []);
    const conventions = state.orientation.euler_conventions ?? [];
    conventionSelect.replaceChildren(
      ...conventions.map((entry) =>
        el('option', {
          value: entry.key,
          text: entry.label,
          title: entry.help_text,
          selected: entry.key === state.convention,
        }),
      ),
    );
    renderFamilyChips();
    renderAxisChips();
    renderPresetChips();
    renderAngleFields();
    draw();
    scheduleReadout({ immediate: true });
  }

  function resetTrail() {
    state.trails = { rd: [], td: [], nd: [] };
  }

  function draw() {
    if (!state.orientation) return;
    const rotation = camera();
    recordTrail(rotation);
    poleFigure.replaceChildren(renderPoleFigure(state, rotation));
    inverseFigure.replaceChildren(renderInverseFigure(state, rotation));
    // The pictures are live; the numbers settle. Scheduling from here means
    // every route that moves the camera -- drag, zoom, an axis button, a
    // preset -- updates the readout without each having to remember to.
    scheduleReadout();
  }

  /* ---------------------------------------------------------------- chips */

  function chip(label, active, color, onToggle, title) {
    return el('button.orient__chip', {
      type: 'button',
      text: label,
      title,
      'aria-pressed': String(active),
      style: color ? `--chip-color: ${color}` : null,
      onclick: onToggle,
    });
  }

  function renderFamilyChips() {
    const families = state.orientation?.pole_families ?? [];
    familyChips.replaceChildren(
      ...families.map((family, index) =>
        chip(
          family.label,
          state.families.has(family.key),
          FAMILY_COLORS[index % FAMILY_COLORS.length],
          () => {
            if (state.families.has(family.key)) state.families.delete(family.key);
            else state.families.add(family.key);
            renderFamilyChips();
            draw();
          },
          `Show the ${family.label} poles (${family.vectors.length} with their antipodes)`,
        ),
      ),
    );
  }

  function renderAxisChips() {
    const axes = state.orientation?.specimen_axes ?? [];
    axisChips.replaceChildren(
      ...axes.map((axis) =>
        chip(
          axis.label,
          state.axes.has(axis.key),
          AXIS_COLORS[axis.key],
          () => {
            if (state.axes.has(axis.key)) state.axes.delete(axis.key);
            else state.axes.add(axis.key);
            state.trails[axis.key] = [];
            renderAxisChips();
            draw();
          },
          axis.help_text,
        ),
      ),
    );
  }

  function renderPresetChips() {
    const components = state.orientation?.components ?? [];
    presetChips.replaceChildren(
      ...(components.length
        ? components.map((component) =>
            chip(
              component.name,
              false,
              null,
              () => setAngles(component.angles_deg, { convention: 'bunge', apply: true }),
              `${component.label} — Bunge ${component.angles_deg
                .map((value) => `${formatNumber(value, 1)}°`)
                .join(', ')}`,
            ),
          )
        : [
            el('span.orient__hint', {
              text: 'Named ideal orientations are a cubic rolling-texture catalogue, so this phase has none.',
            }),
          ]),
    );
  }

  /* --------------------------------------------------------------- angles */

  function angleNames() {
    return state.convention === 'bunge' ? ['φ₁', 'Φ', 'φ₂'] : ['α', 'β', 'γ'];
  }

  function renderAngleFields() {
    const names = angleNames();
    const values = angleInputs.map((input) => input.value);
    angleInputs.length = 0;
    angleFields.replaceChildren(
      ...names.map((name, index) =>
        el('label.orient__angle', {}, [
          el('span', { text: name }),
          (() => {
            const input = el('input', {
              type: 'number',
              step: '0.1',
              value: values[index] ?? '0',
              'aria-label': `${name} in degrees`,
              onfocus: () => {
                state.editing = true;
              },
              onblur: () => {
                state.editing = false;
              },
              onkeydown: (event) => {
                if (event.key === 'Enter') applyAngles();
              },
            });
            angleInputs.push(input);
            return input;
          })(),
          el('span.orient__unit', { text: '°' }),
        ]),
      ),
    );
  }

  function setAngles(angles, { convention = state.convention, apply = false } = {}) {
    if (convention !== state.convention) {
      state.convention = convention;
      conventionSelect.value = convention;
      renderAngleFields();
    }
    angles.forEach((value, index) => {
      if (angleInputs[index]) angleInputs[index].value = formatNumber(value, 2);
    });
    if (apply) applyAngles();
  }

  async function applyAngles() {
    if (!state.orientation) return;
    const values = angleInputs.map((input) => Number(input.value));
    if (values.some((value) => !Number.isFinite(value))) {
      readout.textContent = 'Enter three angles in degrees.';
      return;
    }
    setButton.disabled = true;
    try {
      const result = await call('crystal.orientation', {
        ...request(),
        euler_convention: state.convention,
        angle1: values[0],
        angle2: values[1],
        angle3: values[2],
        camera_matrix: '',
      });
      resetTrail();
      setCamera(result.data.camera_matrix);
      showReadout(result);
    } catch (error) {
      showError(error);
    } finally {
      setButton.disabled = false;
    }
  }

  /* -------------------------------------------------------------- readout */

  /**
   * Ask Python what the current camera is, once the drag has settled.
   *
   * Debounced rather than throttled: during a drag the answer is stale the
   * moment it arrives, and the only request worth making is the last one. The
   * token guards against an earlier reply landing after a later one and
   * rewriting the angle fields with an orientation the user has already turned
   * away from.
   */
  function scheduleReadout({ immediate = false } = {}) {
    if (!state.orientation) return;
    window.clearTimeout(state.readoutTimer);
    state.readoutTimer = window.setTimeout(() => void requestReadout(), immediate ? 0 : READOUT_DELAY);
  }

  async function requestReadout() {
    const token = (state.readoutToken += 1);
    try {
      const result = await call('crystal.orientation', {
        ...request(),
        euler_convention: state.convention,
        camera_matrix: camera().join(' '),
      });
      if (token !== state.readoutToken) return;
      showReadout(result);
    } catch (error) {
      // A readout is a courtesy, not the picture: a failed one leaves the
      // figures alone and says so in place of the numbers, rather than raising
      // a toast for every frame of a drag that crossed a bad state.
      readout.textContent = error?.message ? `Readout unavailable: ${error.message}` : 'Readout unavailable.';
    }
  }

  function showReadout(result) {
    const data = result?.data;
    if (!data) return;
    if (!state.editing) setAngles(data.euler.angles_deg);
    const axis = data.axis_angle.axis.map((value) => formatNumber(value, 3)).join(' ');
    const directions = data.ipf_points
      .map((point) => `${point.label} ∥ ${point.direction}`)
      .join(' · ');
    readout.replaceChildren(
      el('span', { text: `${formatNumber(data.axis_angle.angle_deg, 2)}° about [${axis}]` }),
      el('span.orient__directions', { text: directions }),
    );
  }

  /* --------------------------------------------------------------- trail */

  function recordTrail(rotation) {
    if (!state.showTrail) return;
    for (const axis of state.orientation.specimen_axes ?? []) {
      if (!state.axes.has(axis.key)) continue;
      const reduced = reduceToSector(state.orientation, applyTranspose(rotation, axis.vector));
      if (!reduced) continue;
      const trail = state.trails[axis.key] ?? (state.trails[axis.key] = []);
      const last = trail[trail.length - 1];
      // Angular, not planar: near a triangle corner the projection compresses,
      // and a planar test would drop every position there and leave the trail
      // with a gap exactly where the crystal is doing something interesting.
      if (last && dot(last, reduced) > Math.cos(TRAIL_STEP)) continue;
      trail.push(reduced);
      if (trail.length > TRAIL_LENGTH) trail.shift();
    }
  }

  return { element, controls, setScene, draw, resetTrail, scheduleReadout };
}

/* ------------------------------------------------------------------ maths */

function dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * Stereographic projection of an upper-hemisphere direction.
 *
 * The browser half of `pytex.core.sphere.project_directions` with
 * `method="stereographic"`: `r = tan(rho / 2)`, which in Cartesian components is
 * a divide by `1 + z`. Only this one formula is implemented here, and only for
 * the upper hemisphere, because that is all the dock draws; anything wanting
 * equal area or a signed lower hemisphere asks Python.
 */
function project(vector) {
  const [x, y, z] = vector;
  if (z < 0) return null;
  const divisor = 1 + z;
  if (divisor < 1e-9) return null;
  return [x / divisor, y / divisor];
}

/**
 * Fold a crystal direction into the fundamental sector.
 *
 * Applies each symmetry operator, and each operator's antipode, until one image
 * satisfies every inward edge normal of the sector. Python supplies both the
 * operators and the normals; the test is a dot product, and the fold is the
 * definition of the standard triangle rather than an approximation of it.
 */
function reduceToSector(orientation, vector) {
  const normals = orientation.sector?.edge_normals ?? [];
  const operators = orientation.operators ?? [];
  if (!normals.length || !operators.length) return null;
  let fallback = null;
  for (const operator of operators) {
    for (const sign of [1, -1]) {
      const image = applyMatrix(operator, vector).map((value) => value * sign);
      if (fallback === null) fallback = image;
      if (normals.every((normal) => dot(image, normal) > -1e-9)) return image;
    }
  }
  return fallback;
}

/* ---------------------------------------------------------------- figures */

function figureRoot(label) {
  return svg('svg', {
    viewBox: `${-VIEW - 22} ${-VIEW - 22} ${2 * (VIEW + 22)} ${2 * (VIEW + 22)}`,
    preserveAspectRatio: 'xMidYMid meet',
    role: 'img',
    'aria-label': label,
  });
}

function marker(x, y, color, { radius = 3.4, shape = 'circle', title = '' } = {}) {
  const node =
    shape === 'diamond'
      ? svg('polygon', {
          points: [
            [x, y - radius * 1.3],
            [x + radius * 1.3, y],
            [x, y + radius * 1.3],
            [x - radius * 1.3, y],
          ]
            .map(([px, py]) => `${px},${py}`)
            .join(' '),
        })
      : svg('circle', { cx: x, cy: y, r: radius });
  node.setAttribute('fill', color);
  node.setAttribute('stroke', 'var(--bg-raised)');
  node.setAttribute('stroke-width', '0.9');
  if (title) node.append(svg('title', { text: title }));
  return node;
}

function renderPoleFigure(state, rotation) {
  const orientation = state.orientation;
  const root = figureRoot('Pole figure of the current view');

  root.append(
    svg('circle', {
      cx: 0,
      cy: 0,
      r: VIEW,
      fill: 'none',
      stroke: 'currentColor',
      'stroke-opacity': 0.55,
      'stroke-width': 1,
    }),
    svg('line', { x1: -VIEW, y1: 0, x2: VIEW, y2: 0, stroke: 'currentColor', 'stroke-opacity': 0.16 }),
    svg('line', { x1: 0, y1: -VIEW, x2: 0, y2: VIEW, stroke: 'currentColor', 'stroke-opacity': 0.16 }),
  );

  // RD right and TD up: the projection plane is the screen, so this figure is
  // the structure beside it seen from the same place.
  for (const [label, x, y] of [
    ['RD', VIEW + 12, 2.5],
    ['TD', 0, -VIEW - 8],
  ]) {
    root.append(
      svg('text', {
        x,
        y,
        'text-anchor': 'middle',
        'font-size': 9,
        fill: 'currentColor',
        'fill-opacity': 0.6,
        text: label,
      }),
    );
  }

  (orientation.pole_families ?? []).forEach((family, index) => {
    if (!state.families.has(family.key)) return;
    const color = FAMILY_COLORS[index % FAMILY_COLORS.length];
    for (const vector of family.vectors) {
      const projected = project(applyMatrix(rotation, vector));
      if (!projected) continue;
      root.append(
        marker(projected[0] * VIEW, -projected[1] * VIEW, color, {
          title: `${family.label} pole`,
        }),
      );
    }
  });

  if (state.showOverlays) {
    for (const pole of orientation.overlay_poles ?? []) {
      const projected = project(applyMatrix(rotation, pole.vector));
      if (!projected) continue;
      const x = projected[0] * VIEW;
      const y = -projected[1] * VIEW;
      root.append(
        marker(x, y, 'var(--ink)', { shape: 'diamond', radius: 3.6, title: `${pole.label} pole` }),
        svg('text', {
          x,
          y: y - 7,
          'text-anchor': 'middle',
          'font-size': 8,
          fill: 'var(--ink)',
          'paint-order': 'stroke',
          stroke: 'var(--bg-raised)',
          'stroke-width': 2,
          text: pole.label,
        }),
      );
    }
    for (const direction of orientation.overlay_directions ?? []) {
      const projected = project(applyMatrix(rotation, direction.vector));
      if (!projected) continue;
      root.append(
        marker(projected[0] * VIEW, -projected[1] * VIEW, 'var(--ok)', {
          radius: 2.6,
          title: `${direction.label} direction`,
        }),
      );
    }
  }

  return root;
}

/**
 * The scale and offset that put the standard triangle in the middle of its box.
 *
 * A fundamental sector occupies a small and system-dependent corner of the
 * projection disc — the cubic triangle spans about 0.41 by 0.37 of it — so a
 * figure drawn at the disc's own scale would be mostly empty. The fit is
 * computed once per scene from the outline Python sent, so the triangle fills
 * its cell whatever the point group.
 */
function sectorFit(orientation) {
  const outline = orientation.sector?.outline ?? [];
  if (outline.length < 3) return { scale: VIEW, x: 0, y: 0 };
  const xs = outline.map((point) => point[0]);
  const ys = outline.map((point) => point[1]);
  const minimum = [Math.min(...xs), Math.min(...ys)];
  const maximum = [Math.max(...xs), Math.max(...ys)];
  const span = Math.max(maximum[0] - minimum[0], maximum[1] - minimum[1], 1e-9);
  return {
    scale: (2 * VIEW * 0.86) / span,
    x: (minimum[0] + maximum[0]) / 2,
    y: (minimum[1] + maximum[1]) / 2,
  };
}

function placeInSector(fit, point) {
  return [(point[0] - fit.x) * fit.scale, -(point[1] - fit.y) * fit.scale];
}

function renderInverseFigure(state, rotation) {
  const orientation = state.orientation;
  const root = figureRoot('Inverse pole figure of the current view');
  const fit = state.fit ?? (state.fit = sectorFit(orientation));
  const outline = orientation.sector?.outline ?? [];

  if (outline.length >= 3) {
    root.append(
      svg('polygon', {
        points: outline.map((point) => placeInSector(fit, point).join(',')).join(' '),
        fill: 'var(--accent-soft)',
        'fill-opacity': 0.45,
        stroke: 'currentColor',
        'stroke-opacity': 0.55,
        'stroke-width': 1.2,
      }),
    );
  }

  for (const corner of orientation.sector?.corners ?? []) {
    const projected = project(corner.vector);
    if (!projected) continue;
    const [x, y] = placeInSector(fit, projected);
    root.append(
      svg('circle', { cx: x, cy: y, r: 1.8, fill: 'currentColor', 'fill-opacity': 0.5 }),
      svg('text', {
        x,
        y: y + (y > 0 ? 13 : -7),
        'text-anchor': 'middle',
        'font-size': 9,
        fill: 'currentColor',
        'fill-opacity': 0.75,
        text: corner.label,
      }),
    );
  }

  for (const axis of orientation.specimen_axes ?? []) {
    if (!state.axes.has(axis.key)) continue;
    const color = AXIS_COLORS[axis.key];
    const trail = state.trails[axis.key] ?? [];
    if (state.showTrail) {
      // Dots, not a polyline. Folding into the sector makes the path jump
      // whenever the direction crosses a symmetry boundary, and a line would
      // draw a chord across the triangle that the crystal never took.
      trail.forEach((position, index) => {
        const projected = project(position);
        if (!projected) return;
        const [x, y] = placeInSector(fit, projected);
        root.append(
          svg('circle', {
            cx: x,
            cy: y,
            r: 1.6,
            fill: color,
            'fill-opacity': (0.08 + 0.5 * (index / Math.max(trail.length - 1, 1))).toFixed(3),
          }),
        );
      });
    }
    const reduced = reduceToSector(orientation, applyTranspose(rotation, axis.vector));
    const projected = reduced ? project(reduced) : null;
    if (!projected) continue;
    const [x, y] = placeInSector(fit, projected);
    root.append(
      marker(x, y, color, { radius: 4.2, title: `${axis.label} in the crystal` }),
      svg('text', {
        x,
        y: y - 8,
        'text-anchor': 'middle',
        'font-size': 9,
        fill: color,
        'paint-order': 'stroke',
        stroke: 'var(--bg-raised)',
        'stroke-width': 2,
        text: axis.label,
      }),
    );
  }

  return root;
}
