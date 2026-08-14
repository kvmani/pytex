/**
 * The crystal viewer panel.
 *
 * All the crystallography arrived from Python as finished vertices in Cartesian
 * angstrom. What happens here is *viewing*: a rotation matrix accumulated from
 * the drag, an orthographic projection, and a painter's-algorithm depth sort.
 * That is the entire contract of Decision 4 in the architecture record — the
 * browser makes the picture move, and Python decides what the picture is of.
 *
 * The depth sort deserves a note. Atoms, bonds, plane polygons and cell edges
 * are drawn into one list ordered by camera depth, rather than in layers by
 * kind. Layering is easier and wrong: a translucent (111) polygon drawn on top
 * of everything hides the atoms behind it *and* the ones in front, so it becomes
 * impossible to see which atoms the plane actually passes through, which is the
 * whole reason the plane was drawn.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult, download, saveBlob } from '../core/result.js';
import { call } from '../core/api.js';

export const panel = {
  id: 'crystal',
  title: 'Crystal Viewer',
  tagline: 'Turn a structure in 3D with planes and directions superimposed.',
};

/** Half-width of the drawing area in viewBox units. The scene is scaled to fit. */
const VIEW = 100;

export function mount(context) {
  const operations = context.manifest.operations.filter((entry) => entry.panel === panel.id);
  const sceneOperation = operations.find((entry) => entry.id === 'crystal.scene');
  const renderOperation = operations.find((entry) => entry.id === 'crystal.render');
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const camera = { rotation: identity(), zoom: 1, scale: 1, centre: [0, 0, 0] };
  const state = { scene: null, result: null, teaches: null, form: null };

  // The renderer offers both formats and the help explains when each is right,
  // so the choice has to be on the toolbar. Without it the help describes a
  // control that does not exist, and an SVG -- the one format that stays
  // editable -- is unreachable from the viewer that produces the figure.
  const figureFormat = el(
    'select.button',
    { 'aria-label': 'Figure format', title: 'Format for the published figure' },
    [
      el('option', { value: 'png', text: 'PNG 600 dpi' }),
      el('option', { value: 'svg', text: 'SVG' }),
    ],
  );

  const frame = plotFrame({
    title: 'Structure',
    viewport: false,
    units: 'Å',
    digits: 2,
    toData: (x, y) => cameraToCrystal(camera, x, y),
    formatCursor: (point) =>
      `${formatNumber(point.x, 2)}, ${formatNumber(point.y, 2)}, ${formatNumber(point.z, 2)} Å`,
    toolbar: [
      viewButton('a', [1, 0, 0]),
      viewButton('b', [0, 1, 0]),
      viewButton('c', [0, 0, 1]),
      el('button.button', { type: 'button', text: 'Reset', onclick: () => resetCamera() }),
      figureFormat,
      el('button.button', {
        type: 'button',
        text: 'Figure',
        title: 'Render this exact view through the publication renderer',
        onclick: () => exportFigure(),
      }),
    ],
  });

  const details = el('div');
  context.stage.append(frame.element, details);

  /* ------------------------------------------------------------ controls */

  const formHost = el('div');
  const drawButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Build structure',
    onclick: () => run(),
  });

  context.rail.append(
    formHost,
    drawButton,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'Structures chosen so that one rotation makes a crystallographic point obvious.',
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
    state.form = buildForm(sceneOperation, { initial });
    formHost.replaceChildren(state.form.element);
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    drawButton.disabled = true;
    drawButton.textContent = 'Building…';
    state.form.clearErrors();
    try {
      const result = await call('crystal.scene', state.form.values());
      state.result = result;
      state.scene = result.data.scene;
      resetCamera();
      renderResult(details, result, { teaches: state.teaches });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      drawButton.disabled = false;
      drawButton.textContent = 'Build structure';
    }
  }

  async function exportFigure() {
    if (!state.result) return;
    try {
      const result = await call('crystal.render', {
        ...state.result.inputs,
        // Send the camera itself rather than angles derived here: the
        // conversion to elevation and azimuth lives in Python, so the exported
        // figure cannot drift from the on-screen view through two slightly
        // different derivations of the same thing.
        camera_matrix: camera.rotation.join(' '),
        show_legend: true,
        show_frame_indicator: true,
        format: figureFormat.value,
        dpi: 600,
      });
      // Trim the separator a trailing bracket leaves behind: "Halite (NaCl)"
      // otherwise files as `halite-nacl--structure`.
      const slug = state.scene.phase.name.replace(/\W+/g, '-').replace(/^-|-$/g, '').toLowerCase();
      const name = `${slug}-structure`;
      if (result.data.encoding === 'base64') {
        downloadBase64(`${name}.png`, result.data.image, 'image/png');
      } else {
        download(`${name}.svg`, result.data.image, 'image/svg+xml');
      }
      frame.setStatus(`Figure exported: ${Math.round(result.data.bytes / 1024)} kB.`);
    } catch (error) {
      context.showError(error);
    }
  }

  /* -------------------------------------------------------------- camera */

  function viewButton(label, axis) {
    return el('button.button', {
      type: 'button',
      text: `Along ${label}`,
      title: `Look down the ${label} axis`,
      onclick: () => {
        const vectors = state.scene?.axes;
        if (!vectors) return;
        camera.rotation = lookAlong(vectors[axis.indexOf(1)].vector);
        draw();
      },
    });
  }

  function resetCamera() {
    camera.rotation = multiply(rotationX(-1.2), rotationY(0.6));
    camera.zoom = 1;
    if (state.scene) {
      camera.centre = state.scene.centre;
      camera.scale = (VIEW * 0.82) / (state.scene.radius || 1);
    }
    draw();
  }

  /* ------------------------------------------------------------ drawing */

  function draw() {
    if (!state.scene) return;
    const node = renderScene(state.scene, camera, frame);
    frame.setContent(node);
    attachPointer(node);
    const scene = state.scene;
    frame.setStatus(
      `${scene.atoms.length} atoms · ${scene.bonds.length} bonds · ` +
        `${scene.planes.length} planes · ${scene.directions.length} directions · ` +
        `drag to rotate, scroll to zoom`,
    );
  }

  // The drag is watched on the frame, which outlives the drawing inside it.
  //
  // Turning the crystal redraws it, and a redraw builds a new SVG and throws
  // the old one away. Handlers owned by that SVG therefore see the first
  // movement of a drag and nothing after it: the element they belong to no
  // longer exists, and it took its pointer capture with it. The symptom is a
  // crystal that nudges once however far you pull — the whole of "turn it in
  // your hands" failing on the second frame, with nothing in the console.
  //
  // The frame is created once per mount, so capturing the pointer on it also
  // keeps a drag alive when the pointer leaves the picture mid-turn, and takes
  // the listeners with it when the panel is replaced.
  let dragging = null;

  frame.element.addEventListener('pointerdown', (event) => {
    if (!state.scene || !event.target.closest('svg')) return;
    dragging = { x: event.clientX, y: event.clientY };
    frame.element.setPointerCapture(event.pointerId);
    event.preventDefault();
  });
  frame.element.addEventListener('pointermove', (event) => {
    if (!dragging) return;
    const dx = (event.clientX - dragging.x) * 0.01;
    const dy = (event.clientY - dragging.y) * 0.01;
    dragging = { x: event.clientX, y: event.clientY };
    // Rotate about the *screen* axes, not the model's: pre-multiplying keeps
    // "drag right turns right" true no matter how the crystal is already
    // oriented, which is what makes the control feel like a physical object.
    camera.rotation = multiply(multiply(rotationY(dx), rotationX(dy)), camera.rotation);
    draw();
  });
  for (const ending of ['pointerup', 'pointercancel']) {
    frame.element.addEventListener(ending, (event) => {
      dragging = null;
      if (frame.element.hasPointerCapture?.(event.pointerId)) {
        frame.element.releasePointerCapture(event.pointerId);
      }
    });
  }

  function attachPointer(node) {
    node.addEventListener(
      'wheel',
      (event) => {
        event.preventDefault();
        camera.zoom = Math.min(Math.max(camera.zoom * (event.deltaY < 0 ? 1.12 : 0.89), 0.2), 12);
        draw();
      },
      { passive: false },
    );
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);

  return { help: () => sceneOperation ?? renderOperation };
}

/* ------------------------------------------------------------------ scene */

function renderScene(scene, camera, frame) {
  const root = svg('svg', {
    viewBox: `${-VIEW} ${-VIEW} ${2 * VIEW} ${2 * VIEW}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': 'Crystal structure',
  });
  const defs = svg('defs');
  root.append(defs);

  const project = (point) => projectPoint(camera, point);
  const items = [];

  for (const edge of scene.cell_edges) {
    const a = project(edge[0]);
    const b = project(edge[1]);
    items.push({
      depth: (a.depth + b.depth) / 2 - 1e3, // cell edges sit behind matter at equal depth
      node: svg('line', {
        x1: a.x, y1: a.y, x2: b.x, y2: b.y,
        stroke: 'currentColor',
        'stroke-width': 0.45,
        'stroke-opacity': 0.5,
      }),
    });
  }

  for (const plane of scene.planes) {
    const points = plane.vertices.map(project);
    if (points.length < 3) continue;
    const polygon = svg('polygon', {
      points: points.map((point) => `${point.x},${point.y}`).join(' '),
      fill: plane.color,
      'fill-opacity': plane.alpha ?? 0.28,
      stroke: plane.color,
      'stroke-width': 0.5,
    });
    const depth = points.reduce((sum, point) => sum + point.depth, 0) / points.length;
    items.push({ depth, node: polygon, row: { Plane: plane.label ?? 'plane' } });
    if (plane.label) {
      const centroid = points.reduce(
        (acc, point) => ({ x: acc.x + point.x / points.length, y: acc.y + point.y / points.length }),
        { x: 0, y: 0 },
      );
      items.push({
        depth: depth + 0.5,
        node: svg('text', {
          x: centroid.x, y: centroid.y,
          'text-anchor': 'middle',
          'font-size': 4.5,
          fill: plane.color,
          'paint-order': 'stroke',
          stroke: 'var(--bg-raised)',
          'stroke-width': 1.2,
          text: plane.label,
        }),
      });
    }
  }

  for (const bond of scene.bonds) {
    const a = project(bond.start);
    const b = project(bond.end);
    items.push({
      depth: (a.depth + b.depth) / 2,
      node: svg('line', {
        x1: a.x, y1: a.y, x2: b.x, y2: b.y,
        stroke: bond.color,
        'stroke-width': 1.1,
        'stroke-linecap': 'round',
        'stroke-opacity': 0.85,
      }),
      row: { Bond: bond.species, 'Length / Å': formatNumber(bond.length, 4) },
    });
  }

  for (const atom of scene.atoms) {
    const point = project(atom.position);
    const radius = atom.radius * camera.scale * camera.zoom;
    items.push({
      depth: point.depth,
      node: svg('circle', {
        cx: point.x, cy: point.y, r: radius,
        fill: atom.color,
        stroke: 'var(--bg-raised)',
        'stroke-width': 0.35,
      }),
      row: {
        Element: atom.species,
        Site: atom.label ?? '—',
        'x / Å': formatNumber(atom.position[0], 4),
        'y / Å': formatNumber(atom.position[1], 4),
        'z / Å': formatNumber(atom.position[2], 4),
        Occupancy: formatNumber(atom.occupancy, 3),
      },
    });
    if (atom.label) {
      items.push({
        depth: point.depth + 1e-3,
        node: svg('text', {
          x: point.x, y: point.y + radius * 0.35,
          'text-anchor': 'middle',
          'font-size': Math.max(radius * 0.8, 2.5),
          fill: 'var(--ink)',
          text: atom.label,
        }),
      });
    }
  }

  for (const direction of scene.directions) {
    const a = project(direction.start);
    const b = project(direction.end);
    items.push({
      depth: Math.max(a.depth, b.depth) + 1e3, // arrows stay legible in front
      node: svg('g', {}, [
        svg('line', {
          x1: a.x, y1: a.y, x2: b.x, y2: b.y,
          stroke: direction.color,
          'stroke-width': 1.4,
        }),
        arrowHead(a, b, direction.color),
        direction.label
          ? svg('text', {
              x: b.x, y: b.y - 2,
              'text-anchor': 'middle',
              'font-size': 5,
              fill: direction.color,
              'paint-order': 'stroke',
              stroke: 'var(--bg-raised)',
              'stroke-width': 1.4,
              text: direction.label,
            })
          : null,
      ]),
      row: { Direction: direction.label ?? 'direction' },
    });
  }

  items.sort((left, right) => left.depth - right.depth);
  for (const item of items) {
    root.append(item.node);
    if (item.row) frame.hoverable(item.node, item.row);
  }

  root.append(axisGizmo(scene.axes, camera));
  return root;
}

/**
 * The axis gizmo, without which a rotated cubic cell is unreadable.
 *
 * Drawn in the corner at fixed size and rotated with the camera, so it reports
 * where a, b and c now point rather than where they started.
 */
function axisGizmo(axes, camera) {
  const group = svg('g', { transform: `translate(${-VIEW + 18} ${VIEW - 18})` });
  const length = 12;
  for (const axis of axes ?? []) {
    const direction = normalise(axis.vector);
    const rotated = applyMatrix(camera.rotation, direction);
    const x = rotated[0] * length;
    const y = -rotated[1] * length;
    group.append(
      svg('line', {
        x1: 0, y1: 0, x2: x, y2: y,
        stroke: 'var(--ink-muted)',
        'stroke-width': 0.8,
      }),
      svg('text', {
        x: x * 1.25, y: y * 1.25,
        'text-anchor': 'middle',
        'dominant-baseline': 'middle',
        'font-size': 5,
        fill: 'var(--ink-muted)',
        'font-style': 'italic',
        text: axis.label,
      }),
    );
  }
  return group;
}

function arrowHead(from, to, color) {
  const dx = to.x - from.x;
  const dy = to.y - from.y;
  const length = Math.hypot(dx, dy) || 1;
  const ux = dx / length;
  const uy = dy / length;
  const size = 3.2;
  const points = [
    [to.x, to.y],
    [to.x - ux * size - uy * size * 0.45, to.y - uy * size + ux * size * 0.45],
    [to.x - ux * size + uy * size * 0.45, to.y - uy * size - ux * size * 0.45],
  ];
  return svg('polygon', {
    points: points.map(([x, y]) => `${x},${y}`).join(' '),
    fill: color,
  });
}

/* ----------------------------------------------------------- projection */

/**
 * Project one crystal-frame point to view coordinates.
 *
 * Orthographic, not perspective: a perspective divide makes equal lattice
 * spacings look unequal, and someone reading a structure is entitled to trust
 * that two identical rows of atoms are drawn identically.
 */
function projectPoint(camera, point) {
  const centred = [
    point[0] - camera.centre[0],
    point[1] - camera.centre[1],
    point[2] - camera.centre[2],
  ];
  const rotated = applyMatrix(camera.rotation, centred);
  const scale = camera.scale * camera.zoom;
  return { x: rotated[0] * scale, y: -rotated[1] * scale, depth: rotated[2] };
}

/** Invert the projection for the cursor readout, on the plane through the centre. */
function cameraToCrystal(camera, x, y) {
  const scale = camera.scale * camera.zoom;
  if (!scale) return null;
  const view = [x / scale, -y / scale, 0];
  const model = applyTranspose(camera.rotation, view);
  return {
    x: model[0] + camera.centre[0],
    y: model[1] + camera.centre[1],
    z: model[2] + camera.centre[2],
  };
}

/** A camera looking down a crystal direction. */
function lookAlong(vector) {
  const forward = normalise(vector);
  const reference = Math.abs(forward[1]) > 0.95 ? [0, 0, 1] : [0, 1, 0];
  const right = normalise(cross(reference, forward));
  const up = cross(forward, right);
  return [right[0], right[1], right[2], up[0], up[1], up[2], forward[0], forward[1], forward[2]];
}

/* Row-major 3x3 helpers. Small enough to read, and the only maths in the panel. */

function identity() {
  return [1, 0, 0, 0, 1, 0, 0, 0, 1];
}

function multiply(a, b) {
  const out = new Array(9).fill(0);
  for (let row = 0; row < 3; row += 1) {
    for (let column = 0; column < 3; column += 1) {
      let sum = 0;
      for (let k = 0; k < 3; k += 1) sum += a[row * 3 + k] * b[k * 3 + column];
      out[row * 3 + column] = sum;
    }
  }
  return out;
}

function rotationX(angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return [1, 0, 0, 0, c, -s, 0, s, c];
}

function rotationY(angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return [c, 0, s, 0, 1, 0, -s, 0, c];
}

function applyMatrix(m, v) {
  return [
    m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
    m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
    m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
  ];
}

function applyTranspose(m, v) {
  return [
    m[0] * v[0] + m[3] * v[1] + m[6] * v[2],
    m[1] * v[0] + m[4] * v[1] + m[7] * v[2],
    m[2] * v[0] + m[5] * v[1] + m[8] * v[2],
  ];
}

function cross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

function normalise(v) {
  const length = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / length, v[1] / length, v[2] / length];
}

/**
 * Save a base64 payload, through the one save path the shells share.
 *
 * A PNG arrives base64-encoded because JSON has no bytes. Only the decoding
 * belongs here; handing the blob to an anchor directly is what left the desktop
 * shell unable to save a figure while the SVG beside it saved fine.
 */
function downloadBase64(filename, base64, mime) {
  const binary = atob(base64);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i += 1) bytes[i] = binary.charCodeAt(i);
  return saveBlob(filename, new Blob([bytes], { type: mime }));
}
