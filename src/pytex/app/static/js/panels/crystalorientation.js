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
import { bandLabelNode, labelAngleDeg } from '../core/kikuchilabel.js';
import { applyMatrix, applyTranspose } from '../core/rotation3.js';
import { call } from '../core/api.js';

/** Half-width of each figure's drawing area, in viewBox units. */
const VIEW = 100;

/** Half-width of a figure's viewBox, drawing area plus its label margin. */
const HALF = VIEW + 22;

/** How far the Kikuchi map may be magnified. */
const MAP_ZOOM_LIMITS = { min: 1, max: 16 };

/** Multiplier per wheel notch. A tenth is a step the eye can follow. */
const MAP_ZOOM_STEP = 1.1;

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
    // The pole figure's own fly-by: one entry per recorded moment, each holding
    // the projected position of every pole of every shown family at that moment.
    // Kept per family so switching a chip on does not invent a history it never
    // had.
    poleTrails: {},
    showTrail: true,
    // The Kikuchi map: its result, the request that produced it, and the axis it
    // is centred on as the user typed it.
    kikuchi: null,
    kikuchiRequest: null,
    kikuchiCentre: '0 0 1',
    kikuchiPending: false,
    // What part of the map is on screen: a magnification and the point of the
    // map the figure is centred on, in drawing units. A band is a fraction of a
    // degree wide and its neighbours are degrees apart, so a whole-hemisphere
    // map is a picture of the network rather than something to read indices
    // off; magnifying it is what turns it back into an atlas.
    mapView: { zoom: 1, x: 0, y: 0 },
    mapPan: null,
    showOverlays: true,
    convention: 'bunge',
    fit: null,
    readoutTimer: null,
    readoutToken: 0,
    kikuchiToken: 0,
    editing: false,
  };

  const poleFigure = el('div.orient__canvas');
  const inverseFigure = el('div.orient__canvas');
  const kikuchiFigure = el('div.orient__canvas.orient__canvas--zoomable', {
    title: 'Scroll to magnify the map, drag to move it, double-click to fit it again',
  });
  const kikuchiStatus = el('p.field__hint', { text: 'Centre the map on a zone axis.' });

  /*
   * Wheel to magnify, drag to move, double-click to fit.
   *
   * The listener is not passive, because zooming a figure and scrolling the
   * page past it are different intentions and the page must not do both. It is
   * attached to the container rather than to the SVG, which is replaced on
   * every frame of a drag of the structure.
   *
   * Magnifying about the pointer rather than about the middle is what makes it
   * usable at high zoom: the band under the cursor is the one being examined,
   * so it is the one that must stay still.
   */
  kikuchiFigure.addEventListener(
    'wheel',
    (event) => {
      if (!state.kikuchi) return;
      event.preventDefault();
      const before = mapPointAt(event);
      const factor = MAP_ZOOM_STEP ** (-Math.sign(event.deltaY));
      const zoom = clamp(
        state.mapView.zoom * factor,
        MAP_ZOOM_LIMITS.min,
        MAP_ZOOM_LIMITS.max,
      );
      if (zoom === state.mapView.zoom) return;
      state.mapView.zoom = zoom;
      const after = mapPointAt(event);
      state.mapView.x += before[0] - after[0];
      state.mapView.y += before[1] - after[1];
      clampMapView();
      drawKikuchi();
    },
    { passive: false },
  );

  kikuchiFigure.addEventListener('pointerdown', (event) => {
    if (!state.kikuchi || state.mapView.zoom <= 1) return;
    state.mapPan = { at: mapPointAt(event), id: event.pointerId };
    kikuchiFigure.setPointerCapture(event.pointerId);
  });

  kikuchiFigure.addEventListener('pointermove', (event) => {
    if (!state.mapPan || state.mapPan.id !== event.pointerId) return;
    const now = mapPointAt(event);
    state.mapView.x += state.mapPan.at[0] - now[0];
    state.mapView.y += state.mapPan.at[1] - now[1];
    clampMapView();
    drawKikuchi();
  });

  for (const name of ['pointerup', 'pointercancel', 'pointerleave']) {
    kikuchiFigure.addEventListener(name, (event) => {
      if (state.mapPan?.id !== event.pointerId) return;
      state.mapPan = null;
      if (kikuchiFigure.hasPointerCapture?.(event.pointerId)) {
        kikuchiFigure.releasePointerCapture(event.pointerId);
      }
    });
  }

  kikuchiFigure.addEventListener('dblclick', () => {
    state.mapView = { zoom: 1, x: 0, y: 0 };
    drawKikuchi();
  });

  const kikuchiInput = el('input.orient__axis', {
    type: 'text',
    value: '0 0 1',
    size: 7,
    'aria-label': 'Kikuchi map centre [uvw]',
    title: 'The zone axis at the middle of the map — the direction on the beam',
    onchange: (event) => {
      state.kikuchiCentre = event.currentTarget.value;
      void refreshKikuchi();
    },
  });
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
      /*
       * The Kikuchi map: the band network of the crystal, centred where the
       * user says rather than where the view happens to be.
       *
       * Deliberately *not* turned by the camera. The other two figures are the
       * same view of the same crystal as the structure beside them, and turning
       * with it is their whole point. This one is the atlas: it is fixed to the
       * crystal, and what moves on it is the marker showing where the current
       * view's beam direction falls — which is exactly how a map is read.
       */
      el('figure.orient__cell', {}, [
        el('figcaption.orient__caption', {}, [
          el('span.orient__title', { text: 'Kikuchi map' }),
          el('label.orient__axis-field', {}, [
            el('span', { text: 'centre [uvw]' }),
            kikuchiInput,
          ]),
        ]),
        kikuchiFigure,
        kikuchiStatus,
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
    void refreshKikuchi();
  }

  /* --------------------------------------------------- the map's own view */

  /** Where a pointer event falls in the map's drawing units. */
  function mapPointAt(event) {
    const root = kikuchiFigure.firstElementChild;
    const matrix = root?.getScreenCTM?.();
    if (!matrix) return [state.mapView.x, state.mapView.y];
    // The SVG letterboxes inside its box, so the screen transform is the only
    // honest conversion; deriving one from the bounding rectangle would be
    // wrong by the letterbox margin at every aspect ratio but one.
    const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(matrix.inverse());
    return [point.x, point.y];
  }

  /** Keep the magnified window over the map rather than off its edge. */
  function clampMapView() {
    const span = HALF - HALF / state.mapView.zoom;
    state.mapView.x = clamp(state.mapView.x, -span, span);
    state.mapView.y = clamp(state.mapView.y, -span, span);
  }

  function drawKikuchi() {
    const rotation = camera();
    kikuchiFigure.replaceChildren(renderKikuchiMap(state, rotation));
  }

  function resetTrail() {
    state.trails = { rd: [], td: [], nd: [] };
    state.poleTrails = {};
  }

  function draw() {
    if (!state.orientation) return;
    const rotation = camera();
    recordTrail(rotation);
    poleFigure.replaceChildren(renderPoleFigure(state, rotation));
    inverseFigure.replaceChildren(renderInverseFigure(state, rotation));
    drawKikuchi();
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
    // The residual travels with the label. A general orientation puts no axis on
    // a low-index direction, so "ND ∥ [10-11]" alone would read as an exact
    // statement about a direction that is four degrees away from it.
    const directions = data.ipf_points
      .map((point) => `${point.label} ∥ ${point.direction} ±${formatNumber(point.residual_deg, 1)}°`)
      .join(' · ');
    // Bunge unconditionally, whatever the picker is set to. It is the
    // convention every EBSD file, every ODF section and every published
    // orientation is written in, so it is what a reader needs in front of them
    // even while they are looking at the angles in another one.
    const bunge = (data.euler_bunge?.angles_deg ?? [])
      .map((value) => formatNumber(value, 2))
      .join(', ');
    readout.replaceChildren(
      el('span.orient__bunge', {
        title: 'Bunge ZXZ angles of the view, whichever convention is selected above',
        text: bunge ? `Bunge (φ₁, Φ, φ₂) = ${bunge}°` : '',
      }),
      el('span', { text: `${formatNumber(data.axis_angle.angle_deg, 2)}° about [${axis}]` }),
      el('span.orient__directions', { text: directions }),
    );
  }

  /* --------------------------------------------------------------- trail */

  function recordTrail(rotation) {
    if (!state.showTrail) return;
    recordPoleTrail(rotation);
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

  /**
   * Remember where the poles of each shown family were, for the fly-by.
   *
   * The inverse figure's trail follows three specimen axes; the pole figure has
   * a whole family to follow, so a moment is stored as the family's *set* of
   * projected positions and drawn as a fading cloud. Storing the projection
   * rather than the direction is deliberate here: unlike the inverse figure
   * there is no fold into a sector, so nothing about a past position can change
   * afterwards, and re-projecting ninety stored orientations on every frame
   * would be work done to reach the same answer.
   */
  function recordPoleTrail(rotation) {
    for (const family of state.orientation?.pole_families ?? []) {
      if (!state.families.has(family.key)) {
        delete state.poleTrails[family.key];
        continue;
      }
      const trail = state.poleTrails[family.key] ?? (state.poleTrails[family.key] = []);
      const first = family.vectors[0];
      if (!first) continue;
      // The gate is angular and on one representative pole: the family turns
      // rigidly, so one member's motion is every member's motion.
      const moved = applyMatrix(rotation, first);
      const last = trail[trail.length - 1];
      if (last && dot(last.gate, moved) > Math.cos(TRAIL_STEP)) continue;
      const points = [];
      for (const vector of family.vectors) {
        const projected = project(applyMatrix(rotation, vector));
        if (projected) points.push(projected);
      }
      trail.push({ gate: moved, points });
      if (trail.length > TRAIL_LENGTH) trail.shift();
    }
  }

  /* --------------------------------------------------------- Kikuchi map */

  /**
   * Fetch the band network for the phase, centred where the user asked.
   *
   * Once per phase and centre, not per frame: the map is a property of the
   * crystal, and the only thing the camera changes about it is where the beam
   * marker falls.
   */
  async function refreshKikuchi() {
    const scene = request();
    if (!scene?.phase) return;
    const centre = parseIndices(state.kikuchiCentre);
    if (!centre) {
      kikuchiStatus.textContent = 'Three whole numbers, please — for example 0 0 1.';
      return;
    }
    const payload = { phase: scene.phase, centre_direction: centre.join(' ') };
    const key = JSON.stringify(payload);
    if (key === state.kikuchiRequest && state.kikuchi) return;
    /*
     * The guard above cannot catch a request that is still in flight -- there
     * is no `state.kikuchi` yet to compare against -- so two calls for the same
     * map can be running at once, and the slower one used to land afterwards
     * and refit the view. A user who had zoomed in the meantime watched the map
     * jump back for no reason they could see. A token makes the newest request
     * the only one whose answer is applied.
     */
    const token = (state.kikuchiToken += 1);
    state.kikuchiPending = true;
    kikuchiStatus.textContent = 'Computing the band network…';
    try {
      const computed = await call('crystal.kikuchi_map', payload);
      if (token !== state.kikuchiToken) return;
      // A *new* map is a new picture; a window into the old one means nothing
      // on it, so the view is fitted again. The same map arriving twice is not
      // a new picture, and refitting there would only discard a magnification
      // the user asked for.
      const isNewMap = key !== state.kikuchiRequest;
      state.kikuchi = computed;
      state.kikuchiRequest = key;
      if (isNewMap) state.mapView = { zoom: 1, x: 0, y: 0 };
      const data = state.kikuchi.data;
      kikuchiStatus.textContent =
        `${data.bands.length} bands, ${data.zone_axes.length} zone axes within ` +
        `${data.max_polar_angle_deg}° of ${data.centre_label}. The cross is the ` +
        'direction now on the beam. Scroll to magnify, drag to move, double-click to fit.';
      draw();
    } catch (error) {
      if (token !== state.kikuchiToken) return;
      // A map is an atlas beside the figures, not one of them: if it cannot be
      // computed the other two must carry on turning with the crystal.
      state.kikuchi = null;
      state.kikuchiRequest = null;
      kikuchiStatus.textContent = error?.message
        ? `No map: ${error.message}`
        : 'No map for this phase and centre.';
      kikuchiFigure.replaceChildren();
    } finally {
      if (token === state.kikuchiToken) state.kikuchiPending = false;
    }
  }

  return { element, controls, setScene, draw, resetTrail, scheduleReadout };
}

/* ------------------------------------------------------------------ maths */

/**
 * Read three whole numbers from a typed direction.
 *
 * Accepts `0 0 1`, `0,0,1`, `[0 0 1]` and `[001]` — the last because a
 * crystallographer writing a low-index axis writes it closed up, and refusing
 * that would be pedantry about notation rather than about the crystal. Returns
 * null for anything else, including the zero direction, which is not a
 * direction.
 */
function parseIndices(text) {
  const cleaned = String(text ?? '').replace(/[[\]()]/g, ' ').trim();
  let parts = cleaned.split(/[\s,]+/).filter(Boolean);
  if (parts.length === 1 && /^-?\d{3}$/.test(parts[0]) === false && /^\d{3}$/.test(parts[0])) {
    parts = parts[0].split('');
  }
  if (parts.length !== 3) return null;
  const values = parts.map((part) => Number(part));
  if (values.some((value) => !Number.isInteger(value))) return null;
  if (values.every((value) => value === 0)) return null;
  return values;
}

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

function clamp(value, low, high) {
  return Math.min(Math.max(value, low), high);
}

/**
 * A figure's drawing surface.
 *
 * `view` magnifies: the window is `HALF / zoom` about `(x, y)` in drawing
 * units, which is the whole box at zoom one. Zooming through the viewBox rather
 * than through a transform means every drawn coordinate stays in the map's own
 * units, so nothing downstream has to know whether the figure is magnified.
 */
function figureRoot(label, view = null) {
  const zoom = Math.max(view?.zoom ?? 1, 1);
  const half = HALF / zoom;
  const x = view?.x ?? 0;
  const y = view?.y ?? 0;
  return svg('svg', {
    viewBox: `${x - half} ${y - half} ${2 * half} ${2 * half}`,
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
    /*
     * The fly-by, under the live poles.
     *
     * Dots rather than a polyline, for the same reason the inverse figure uses
     * them: a pole leaving the upper hemisphere simply stops being drawn, and a
     * line would join its last position to wherever the *next* pole of the
     * family happened to be — a path the crystal never took. Older positions
     * fade, so the direction of travel is legible without an arrow.
     */
    if (state.showTrail) {
      const trail = state.poleTrails[family.key] ?? [];
      trail.forEach((moment, position) => {
        const fade = (0.06 + 0.42 * (position / Math.max(trail.length - 1, 1))).toFixed(3);
        for (const point of moment.points) {
          root.append(
            svg('circle', {
              cx: point[0] * VIEW,
              cy: -point[1] * VIEW,
              r: 1.5,
              fill: color,
              'fill-opacity': fade,
            }),
          );
        }
      });
    }
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
 * The Kikuchi map: the crystal's whole band network about a chosen axis.
 *
 * Every band is drawn as its two edges — what a plate actually records — with
 * the plane's trace dashed between them, because the trace is a construction
 * rather than something visible. Zone axes are marked where bands cross, sized
 * by how many bands meet there, since that count is the n-fold symmetry of the
 * pattern seen on arriving.
 *
 * The map does not turn with the camera; the marker on it does. The map is the
 * atlas — fixed to the crystal, centred where the user asked — and the cross
 * shows which direction the current view has on the beam, projected into the
 * map's own frame by the matrix the service sent. Turning the crystal moves the
 * cross across a stationary map, which is how a map is read.
 */
/**
 * A band's indices, written along the band and beside it.
 *
 * Where the name goes is not a free choice. Bands crowd towards the middle of a
 * stereogram — at the centre every band of the zone crosses at once — and they
 * leave the drawing at its rim, so a name is placed part of the way out from
 * whatever the reader currently has in the middle. That "whatever" is the point
 * of taking the view: magnify a corner of the map and the names follow the
 * bands into it instead of staying behind at the rim of a picture nobody is
 * looking at any more.
 *
 * @returns {SVGElement|null} Null when no sampled point of the band falls in a
 *   readable place, which is the honest outcome for a band that barely clips
 *   the visible window.
 */
function bandName(band, { at, inside, view, shrink }) {
  const target = 0.55 * HALF * shrink;
  let best = null;
  for (const run of band.centre) {
    const points = run.filter(inside).map(at);
    for (let index = 1; index < points.length; index += 1) {
      const [x, y] = points[index];
      const miss = Math.abs(Math.hypot(x - view.x, y - view.y) - target);
      if (!best || miss < best.miss) best = { miss, point: points[index], before: points[index - 1] };
    }
  }
  if (!best || best.miss > target) return null;
  const [x, y] = best.point;
  const dx = x - best.before[0];
  const dy = y - best.before[1];
  const length = Math.hypot(dx, dy);
  if (!(length > 0)) return null;
  // Beside the band, along its normal, by about the height of the text: on the
  // line the glyphs would sit on top of the very thing being named.
  const clearance = 5 * shrink;
  return bandLabelNode({
    x: x - (dy / length) * clearance,
    y: y + (dx / length) * clearance,
    angleDeg: labelAngleDeg(best.before, best.point),
    text: band.label,
    fontSize: 7.5 * shrink,
    colour: 'var(--accent)',
    haloColour: 'var(--bg-raised)',
    haloWidth: 1.8 * shrink,
  });
}

function renderKikuchiMap(state, rotation) {
  const data = state.kikuchi?.data;
  const mapView = state.mapView ?? { zoom: 1, x: 0, y: 0 };
  const root = figureRoot('Kikuchi map of the crystal', mapView);
  if (!data) return root;

  // Text and markers are annotation, not geometry: they keep their size on the
  // screen while the map grows under them. Line weights follow the same rule,
  // so a magnified band does not turn into a stripe.
  const zoom = Math.max(mapView.zoom, 1);
  const shrink = 1 / zoom;
  const radius = Math.max(data.boundary_radius ?? 1, 1e-6);
  const scale = VIEW / radius;
  const at = (point) => [point[0] * scale, -point[1] * scale];
  const inside = (point) => Math.hypot(point[0], point[1]) <= radius * 1.001;

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
  );

  const polyline = (run, attrs) => {
    // Clipped to the mapped region rather than to the drawing box: a band that
    // leaves the region is not part of this map, and letting it run to the edge
    // of the cell would claim a wider map than was computed.
    const points = run.filter(inside).map(at);
    if (points.length < 2) return null;
    return svg('polyline', {
      points: points.map(([x, y]) => `${x},${y}`).join(' '),
      fill: 'none',
      ...attrs,
    });
  };

  // How many bands are named. The bands arrive strongest first, so this keeps
  // the ones a reader identifies the map by; magnification buys room, so it
  // buys names — which is the point of being able to magnify at all.
  const named = Math.min(data.bands.length, Math.round(6 * zoom ** 0.9));
  data.bands.forEach((band, index) => {
    // Width carries prominence, so the strong bands read first — the same
    // ordering a real map is legible by.
    const weight = (0.5 + 1.1 * Math.min(1, Number(band.intensity) || 0)) * shrink;
    for (const edge of band.edges) {
      for (const run of edge) {
        root.append(
          polyline(run, {
            stroke: 'var(--ink)',
            'stroke-opacity': 0.55,
            'stroke-width': weight,
            'stroke-linecap': 'round',
          }),
        );
      }
    }
    for (const run of band.centre) {
      root.append(
        polyline(run, {
          stroke: 'var(--accent)',
          'stroke-opacity': 0.4,
          'stroke-width': weight * 0.6,
          'stroke-dasharray': `${3 * shrink} ${3 * shrink}`,
        }),
      );
    }
    const name = index < named ? bandName(band, { at, inside, view: mapView, shrink }) : null;
    if (name) root.append(name);
  });

  for (const axis of data.zone_axes) {
    if (!inside([axis.x, axis.y])) continue;
    const [x, y] = at([axis.x, axis.y]);
    const node = marker(x, y, 'var(--violet)', {
      radius: 2.2 + 0.35 * Math.min(axis.order, 8),
      title: `${axis.label} — ${axis.order} bands meet here`,
    });
    root.append(
      node,
      svg('text', {
        x,
        y: y - 6,
        'text-anchor': 'middle',
        'font-size': 8,
        fill: 'var(--ink)',
        'paint-order': 'stroke',
        stroke: 'var(--bg-raised)',
        'stroke-width': 2,
        text: axis.label,
      }),
    );
  }

  // Where the current view looks: ND in crystal coordinates, carried into the
  // map frame by the matrix the map was built with.
  const view = data.view_matrix;
  if (view && rotation) {
    const beamCrystal = applyTranspose(rotation, [0, 0, 1]);
    const inMap = applyMatrix(view, beamCrystal);
    const projected = project(inMap[2] >= 0 ? inMap : inMap.map((value) => -value));
    if (projected && inside(projected)) {
      const [x, y] = at(projected);
      const arm = 6;
      root.append(
        svg('line', {
          x1: x - arm, y1: y, x2: x + arm, y2: y,
          stroke: 'var(--ok)', 'stroke-width': 1.8,
        }),
        svg('line', {
          x1: x, y1: y - arm, x2: x, y2: y + arm,
          stroke: 'var(--ok)', 'stroke-width': 1.8,
        }),
        svg('circle', {
          cx: x, cy: y, r: 3.2,
          fill: 'none', stroke: 'var(--ok)', 'stroke-width': 1.2,
        }),
        svg('title', { text: 'The direction this view has on the beam' }),
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
