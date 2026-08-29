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
import {
  applyMatrix,
  applyTranspose,
  identity,
  lookAlong,
  multiply,
  normalise,
  rotationX,
  rotationY,
} from '../core/rotation3.js';
import { orientationDock } from './crystalorientation.js';
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

/* Zoom limits. Reset draws the scene at 100%, and below it is as useful as
 * above: a large supercell seen whole, with room around it, is what "zoom out"
 * is asked for. */
const ZOOM_MIN = 0.2;
const ZOOM_MAX = 12;

const DEFAULT_APPEARANCE = Object.freeze({
  showAtoms: true,
  showBonds: true,
  showCells: true,
  showPlanes: true,
  showDirections: true,
  showLabels: true,
  showGizmo: true,
  atomScale: 1,
  atomOpacity: 1,
  surfaceFinish: 'glossy',
  lightAzimuth: -135,
  lightElevation: 45,
  ambientLight: 0.42,
  diffuseLight: 0.78,
  specularLight: 0.38,
  atomShininess: 26,
  depthCue: 0.18,
  bondColor: '#64748b',
  bondWidth: 1,
  bondOpacity: 0.85,
  cellColor: '#64748b',
  cellWidth: 1,
  cellOpacity: 0.5,
  planeColor: '#0f766e',
  planeOpacity: 0.28,
  directionColor: '#2563eb',
  directionWidth: 1,
  directionOpacity: 0.96,
  annotationScale: 1,
});

export function defaultAppearance() {
  return { ...DEFAULT_APPEARANCE, speciesColors: {} };
}

export function seedSpeciesColors(appearance, scene) {
  for (const atom of scene?.atoms ?? []) {
    if (!appearance.speciesColors[atom.species]) {
      appearance.speciesColors[atom.species] = atom.color;
    }
  }
}

function publicationAppearance(appearance, camera) {
  const screenLight = lightingDirection(appearance);
  const lightDirection = applyTranspose(camera.rotation, [
    screenLight.x,
    -screenLight.y,
    screenLight.z,
  ]);
  return {
    show_atoms: appearance.showAtoms,
    show_bonds: appearance.showBonds,
    show_cells: appearance.showCells,
    show_planes: appearance.showPlanes,
    show_directions: appearance.showDirections,
    show_labels: appearance.showLabels,
    show_gizmo: appearance.showGizmo,
    atom_scale: appearance.atomScale,
    atom_opacity: appearance.atomOpacity,
    surface_finish: appearance.surfaceFinish,
    light_direction: lightDirection,
    light_ambient: appearance.ambientLight,
    light_diffuse: appearance.diffuseLight,
    light_specular: appearance.specularLight,
    atom_shininess: appearance.atomShininess,
    depth_cue_strength: appearance.depthCue,
    species_colors: appearance.speciesColors,
    bond_color: appearance.bondColor,
    bond_width: appearance.bondWidth,
    bond_opacity: appearance.bondOpacity,
    cell_color: appearance.cellColor,
    cell_width: appearance.cellWidth,
    cell_opacity: appearance.cellOpacity,
    plane_color: appearance.planeColor,
    plane_opacity: appearance.planeOpacity,
    direction_color: appearance.directionColor,
    direction_width: appearance.directionWidth,
    direction_opacity: appearance.directionOpacity,
    annotation_scale: appearance.annotationScale,
  };
}

function appearanceSelect(label, hint, appearance, key, options, onChange) {
  return el('label.field', {}, [
    el('span.field__label', { text: label }),
    el(
      'select',
      {
        onchange: (event) => {
          appearance[key] = event.currentTarget.value;
          onChange();
        },
      },
      options.map(([value, text]) =>
        el('option', { value, text, selected: appearance[key] === value }),
      ),
    ),
    el('span.field__hint', { text: hint }),
  ]);
}

function appearanceToggle(label, appearance, key, onChange) {
  return el('label.object-toggle', {}, [
    el('input', {
      type: 'checkbox',
      checked: appearance[key],
      onchange: (event) => {
        appearance[key] = event.currentTarget.checked;
        onChange();
      },
    }),
    el('span', { text: label }),
  ]);
}

function appearanceRange(label, hint, appearance, key, { min, max, step, suffix = '' }, onChange) {
  const output = el('output', { text: `${Number(appearance[key]).toFixed(2)}${suffix}` });
  return el('label.field', {}, [
    el('span.field__label', { text: label }),
    el('span.range-control', {}, [
      el('input', {
        type: 'range', min, max, step, value: appearance[key],
        oninput: (event) => {
          appearance[key] = Number(event.currentTarget.value);
          output.textContent = `${appearance[key].toFixed(2)}${suffix}`;
          onChange();
        },
      }),
      output,
    ]),
    el('span.field__hint', { text: hint }),
  ]);
}

function appearanceColor(label, hint, appearance, key, onChange) {
  const output = el('output', { text: appearance[key].toUpperCase() });
  return el('label.field', {}, [
    el('span.field__label', { text: label }),
    el('span.color-control', {}, [
      el('input', {
        type: 'color', value: appearance[key],
        oninput: (event) => {
          appearance[key] = event.currentTarget.value;
          output.textContent = appearance[key].toUpperCase();
          onChange();
        },
      }),
      output,
    ]),
    el('span.field__hint', { text: hint }),
  ]);
}

function appearanceControl(appearance, scene, { onChange, onReset }) {
  seedSpeciesColors(appearance, scene);
  const species = Object.entries(appearance.speciesColors).sort(([left], [right]) =>
    left.localeCompare(right),
  );
  return el('details.group.appearance', {}, [
    el('summary', { text: 'Object properties' }),
    el('div.group__body', {}, [
      el('p.field__help', {
        text: 'Presentation only: these controls redraw existing geometry. Positions, planes, indices and exported data do not change.',
      }),
      el('div.object-toggle-grid', {}, [
        appearanceToggle('Atoms', appearance, 'showAtoms', onChange),
        appearanceToggle('Bonds', appearance, 'showBonds', onChange),
        appearanceToggle('Cells', appearance, 'showCells', onChange),
        appearanceToggle('Planes', appearance, 'showPlanes', onChange),
        appearanceToggle('Directions', appearance, 'showDirections', onChange),
        appearanceToggle('Labels', appearance, 'showLabels', onChange),
        appearanceToggle('Axis gizmo', appearance, 'showGizmo', onChange),
      ]),
      el('h3.group__subheading', { text: 'Atoms' }),
      appearanceRange('Atom size', 'Scale every atomic radius without moving its centre.', appearance, 'atomScale', {
        min: 0.2, max: 2.5, step: 0.05, suffix: '×',
      }, onChange),
      appearanceRange('Atom opacity', 'Lower opacity reveals bonds and planes inside dense cells.', appearance, 'atomOpacity', {
        min: 0.1, max: 1, step: 0.05,
      }, onChange),
      el('h3.group__subheading', { text: 'Lighting and depth' }),
      appearanceSelect(
        'Surface finish',
        'Flat is diagrammatic; matte and glossy use sphere-normal shading.',
        appearance,
        'surfaceFinish',
        [
          ['glossy', 'Glossy spheres'],
          ['matte', 'Matte spheres'],
          ['flat', 'Flat colour'],
        ],
        onChange,
      ),
      appearanceRange('Light azimuth', 'Moves the studio light around the screen.', appearance, 'lightAzimuth', {
        min: -180, max: 180, step: 5, suffix: '°',
      }, onChange),
      appearanceRange('Light elevation', 'Higher values move the highlight toward each sphere centre.', appearance, 'lightElevation', {
        min: 0, max: 90, step: 5, suffix: '°',
      }, onChange),
      appearanceRange('Ambient light', 'Base illumination retained on the shadowed limb.', appearance, 'ambientLight', {
        min: 0.05, max: 1, step: 0.01,
      }, onChange),
      appearanceRange('Diffuse light', 'Controls broad light-to-shadow modelling of each surface.', appearance, 'diffuseLight', {
        min: 0, max: 1.25, step: 0.01,
      }, onChange),
      appearanceRange('Specular highlight', 'Controls the bright reflection on glossy atoms and bonds.', appearance, 'specularLight', {
        min: 0, max: 1, step: 0.01,
      }, onChange),
      appearanceRange('Highlight sharpness', 'Higher values make a smaller, harder highlight.', appearance, 'atomShininess', {
        min: 2, max: 96, step: 1,
      }, onChange),
      appearanceRange('Depth cue', 'Fades distant atoms so foreground structure separates immediately.', appearance, 'depthCue', {
        min: 0, max: 0.75, step: 0.01,
      }, onChange),
      ...(species.length
        ? [el('div.object-style-list', {}, species.map(([name, color]) => {
            const output = el('output', { text: color.toUpperCase() });
            return el('label.object-style-row', {}, [
              el('span', { text: name }),
              el('input', {
                type: 'color', value: color, 'aria-label': `${name} atom colour`,
                oninput: (event) => {
                  appearance.speciesColors[name] = event.currentTarget.value;
                  output.textContent = event.currentTarget.value.toUpperCase();
                  onChange();
                },
              }),
              output,
            ]);
          }))]
        : [el('p.field__hint', { text: 'Build a structure to edit its per-species colours.' })]),
      el('h3.group__subheading', { text: 'Bonds and cells' }),
      appearanceColor('Bond colour', 'A uniform interactive colour keeps the network legible.', appearance, 'bondColor', onChange),
      appearanceRange('Bond width', 'Scales the bond stroke or cylinder diameter.', appearance, 'bondWidth', {
        min: 0.2, max: 3, step: 0.05, suffix: '×',
      }, onChange),
      appearanceRange('Bond opacity', 'Reduce when bonds obscure atom positions.', appearance, 'bondOpacity', {
        min: 0.05, max: 1, step: 0.05,
      }, onChange),
      appearanceColor('Cell colour', 'Colour of the direct-basis cell edges.', appearance, 'cellColor', onChange),
      appearanceRange('Cell width', 'Scales every cell edge.', appearance, 'cellWidth', {
        min: 0.2, max: 3, step: 0.05, suffix: '×',
      }, onChange),
      appearanceRange('Cell opacity', 'Keep the cell visible without competing with atoms.', appearance, 'cellOpacity', {
        min: 0.05, max: 1, step: 0.05,
      }, onChange),
      el('h3.group__subheading', { text: 'Planes, directions and labels' }),
      appearanceColor('Plane colour', 'Applied to all bounded plane overlays in this view.', appearance, 'planeColor', onChange),
      appearanceRange('Plane opacity', 'Low values reveal which atoms lie before and behind a plane.', appearance, 'planeOpacity', {
        min: 0.02, max: 0.85, step: 0.01,
      }, onChange),
      appearanceColor('Direction colour', 'Applied to direction arrows and their labels.', appearance, 'directionColor', onChange),
      appearanceRange('Direction width', 'Scales arrow strokes and publication linewidth.', appearance, 'directionWidth', {
        min: 0.2, max: 3, step: 0.05, suffix: '×',
      }, onChange),
      appearanceRange('Direction opacity', 'Reduce for dense overlay comparisons.', appearance, 'directionOpacity', {
        min: 0.05, max: 1, step: 0.05,
      }, onChange),
      appearanceRange('Annotation size', 'Scales atom, plane and direction labels together.', appearance, 'annotationScale', {
        min: 0.5, max: 2.5, step: 0.05, suffix: '×',
      }, onChange),
      el('button.button', { type: 'button', text: 'Reset object properties', onclick: onReset }),
    ]),
  ]);
}

export function mount(context) {
  const operations = context.manifest.operations.filter((entry) => entry.panel === panel.id);
  const sceneOperation = operations.find((entry) => entry.id === 'crystal.scene');
  const renderOperation = operations.find((entry) => entry.id === 'crystal.render');
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  // `pan` is a translation in view coordinates, applied after the projection:
  // the 3-D camera has no viewBox to move, so panning is the only way to bring
  // a corner of a large cell to the middle of the picture without turning it.
  const camera = {
    rotation: identity(),
    zoom: 1,
    scale: 1,
    centre: [0, 0, 0],
    pan: { x: 0, y: 0 },
    panTool: false,
  };
  const state = {
    scene: null,
    result: null,
    teaches: null,
    form: null,
    appearance: defaultAppearance(),
  };

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

  /*
   * The viewer owns its own camera, so it also owns the buttons the shared plot
   * frame gives every other panel. They are the same controls in the same order
   * — zoom out, readout, zoom in, pan — because a reader moving between the
   * structure and a pole figure should not have to learn two viewers.
   */
  const zoomReadout = el('output.plot__zoom', { text: '100%', title: 'Current view zoom' });
  const setZoom = (value) => {
    camera.zoom = Math.min(Math.max(value, ZOOM_MIN), ZOOM_MAX);
    zoomReadout.textContent = `${Math.round(camera.zoom * 100)}%`;
    draw();
  };
  const zoomOutButton = el('button.button.button--icon', {
    type: 'button', text: '−', title: 'Zoom out', 'aria-label': 'Zoom out',
    onclick: () => setZoom(camera.zoom / 1.2),
  });
  const zoomInButton = el('button.button.button--icon', {
    type: 'button', text: '+', title: 'Zoom in', 'aria-label': 'Zoom in',
    onclick: () => setZoom(camera.zoom * 1.2),
  });
  const panButton = el('button.button.button--icon', {
    type: 'button',
    text: '✥',
    title: 'Pan tool: drag to move the structure instead of turning it',
    'aria-label': 'Pan tool',
    'aria-pressed': 'false',
    onclick: () => {
      camera.panTool = !camera.panTool;
      panButton.setAttribute('aria-pressed', String(camera.panTool));
      const surface = frame.element.querySelector('svg');
      if (surface) surface.dataset.pan = camera.panTool ? 'tool' : '';
    },
  });

  const frame = plotFrame({
    title: 'Structure',
    viewport: false,
    units: 'Å',
    digits: 2,
    toData: (x, y) => cameraToCrystal(camera, x, y),
    formatCursor: (point) =>
      `${formatNumber(point.x, 2)}, ${formatNumber(point.y, 2)}, ${formatNumber(point.z, 2)} Å`,
    toolbar: [
      zoomOutButton,
      zoomReadout,
      zoomInButton,
      panButton,
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

  /*
   * The orientation dock sits under the structure rather than beside it. The
   * structure is the hero and wants the width; the two figures beside it are
   * square and small, and a strip of three cards that reflows to one column on
   * a narrow window uses the space a wide desktop actually has without
   * squeezing the picture that earns it. It collapses, because a reader who
   * came to look at bonds should be able to put it away.
   */
  const dock = orientationDock({
    camera: () => camera.rotation,
    setCamera: (matrix) => {
      camera.rotation = matrix.slice();
      draw();
    },
    request: () => ({
      phase: state.result?.inputs?.phase,
      poles: state.result?.inputs?.planes ?? [],
    }),
    showError: (error) => context.showError(error),
  });

  const details = el('div');
  context.stage.append(
    el('div.crystal-workspace', {}, [frame.element, dock.element]),
    details,
  );

  /* ------------------------------------------------------------ controls */

  const formHost = el('div');
  const appearanceHost = el('div');
  const drawButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Build structure',
    onclick: () => run(),
  });

  context.rail.append(
    formHost,
    drawButton,
    dock.controls,
    appearanceHost,
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

  function renderAppearanceControls() {
    appearanceHost.replaceChildren(
      appearanceControl(state.appearance, state.scene, {
        onChange: () => draw(),
        onReset: () => {
          state.appearance = defaultAppearance();
          seedSpeciesColors(state.appearance, state.scene);
          renderAppearanceControls();
          draw();
        },
      }),
    );
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
      seedSpeciesColors(state.appearance, state.scene);
      renderAppearanceControls();
      dock.setScene(state.scene);
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
        appearance: publicationAppearance(state.appearance, camera),
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
        dock.resetTrail();
        draw();
      },
    });
  }

  function resetCamera() {
    camera.rotation = multiply(rotationX(-1.2), rotationY(0.6));
    dock.resetTrail();
    camera.zoom = 1;
    camera.pan = { x: 0, y: 0 };
    zoomReadout.textContent = '100%';
    if (state.scene) {
      camera.centre = state.scene.centre;
      camera.scale = (VIEW * 0.82) / (state.scene.radius || 1);
    }
    draw();
  }

  /* ------------------------------------------------------------ drawing */

  function draw() {
    if (!state.scene) return;
    const node = renderScene(state.scene, camera, frame, state.appearance);
    frame.setContent(node);
    attachPointer(node);
    dock.draw();
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
    // Panning and turning are the same gesture told apart by modifier, tool or
    // button, in the same way as on every other plot in the application.
    const panning = event.button === 1 || event.shiftKey || camera.panTool;
    dragging = { x: event.clientX, y: event.clientY, panning };
    frame.element.setPointerCapture(event.pointerId);
    event.preventDefault();
  });
  frame.element.addEventListener('pointermove', (event) => {
    if (!dragging) return;
    const moveX = event.clientX - dragging.x;
    const moveY = event.clientY - dragging.y;
    dragging = { ...dragging, x: event.clientX, y: event.clientY };
    if (dragging.panning) {
      // Screen pixels to view units: the drawing is VIEW-wide in a box the
      // browser letterboxes, so a pixel is worth 2·VIEW / the rendered width.
      const surface = frame.element.querySelector('svg');
      const width = surface?.getBoundingClientRect().width || 1;
      const height = surface?.getBoundingClientRect().height || 1;
      const perPixel = (2 * VIEW) / Math.min(width, height);
      camera.pan = {
        x: camera.pan.x + moveX * perPixel,
        y: camera.pan.y + moveY * perPixel,
      };
      draw();
      return;
    }
    // Rotate about the *screen* axes, not the model's: pre-multiplying keeps
    // "drag right turns right" true no matter how the crystal is already
    // oriented, which is what makes the control feel like a physical object.
    camera.rotation = multiply(
      multiply(rotationY(moveX * 0.01), rotationX(moveY * 0.01)),
      camera.rotation,
    );
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
        setZoom(camera.zoom * (event.deltaY < 0 ? 1.12 : 0.89));
      },
      { passive: false },
    );
  }

  renderControls();
  renderAppearanceControls();
  if (examples.length) loadExample(examples[0]);

  return { help: () => sceneOperation ?? renderOperation };
}

/* ------------------------------------------------------------------ scene */

function clamp(value, minimum = 0, maximum = 1) {
  return Math.min(Math.max(value, minimum), maximum);
}

function rgb(color) {
  const value = String(color).replace('#', '');
  return [0, 2, 4].map((offset) => Number.parseInt(value.slice(offset, offset + 2), 16));
}

function hex(values) {
  return `#${values.map((value) => Math.round(clamp(value, 0, 255)).toString(16).padStart(2, '0')).join('')}`;
}

function mixColor(color, target, amount) {
  const source = rgb(color);
  const destination = rgb(target);
  const fraction = clamp(amount);
  return hex(source.map((value, index) => value + (destination[index] - value) * fraction));
}

function scaleColor(color, factor) {
  return hex(rgb(color).map((value) => value * factor));
}

function lightingDirection(appearance) {
  const azimuth = (appearance.lightAzimuth * Math.PI) / 180;
  const elevation = (appearance.lightElevation * Math.PI) / 180;
  const radial = Math.cos(elevation);
  return {
    x: Math.cos(azimuth) * radial,
    y: Math.sin(azimuth) * radial,
    z: Math.sin(elevation),
  };
}

function spherePaint(defs, cache, color, appearance) {
  if (appearance.surfaceFinish === 'flat') return color;
  const key = [
    color,
    appearance.surfaceFinish,
    appearance.lightAzimuth,
    appearance.lightElevation,
    appearance.ambientLight,
    appearance.diffuseLight,
    appearance.specularLight,
    appearance.atomShininess,
  ].join('|');
  if (cache.has(key)) return `url(#${cache.get(key)})`;

  const light = lightingDirection(appearance);
  const id = `atom-sphere-${cache.size}`;
  const glossy = appearance.surfaceFinish === 'glossy';
  const specular = appearance.specularLight * (glossy ? 1 : 0.12);
  const ambient = appearance.ambientLight;
  const diffuse = appearance.diffuseLight;
  const highlightStop = clamp(0.3 - appearance.atomShininess * 0.0023, 0.07, 0.3);
  const highlight = mixColor(
    scaleColor(color, clamp(ambient + diffuse, 0.12, 1.18)),
    '#ffffff',
    specular * 0.88,
  );
  const lit = scaleColor(color, clamp(ambient + diffuse * 0.82, 0.12, 1.18));
  const middle = scaleColor(color, clamp(ambient + diffuse * 0.42, 0.1, 1.02));
  const limb = scaleColor(color, clamp(ambient + diffuse * 0.12, 0.08, 0.72));
  const gradient = svg('radialGradient', {
    id,
    cx: '50%', cy: '50%', r: '62%',
    fx: `${50 + light.x * 30}%`, fy: `${50 + light.y * 30}%`,
  }, [
    svg('stop', { offset: '0%', 'stop-color': highlight }),
    svg('stop', { offset: `${Math.round(highlightStop * 100)}%`, 'stop-color': lit }),
    svg('stop', { offset: '58%', 'stop-color': middle }),
    svg('stop', { offset: '84%', 'stop-color': limb }),
    svg('stop', { offset: '100%', 'stop-color': scaleColor(limb, 0.68) }),
  ]);
  defs.append(gradient);
  cache.set(key, id);
  return `url(#${id})`;
}

function bondGlyph(a, b, color, width, opacity, appearance) {
  if (appearance.surfaceFinish === 'flat') {
    return svg('line', {
      x1: a.x, y1: a.y, x2: b.x, y2: b.y,
      stroke: color,
      'stroke-width': width,
      'stroke-linecap': 'round',
      'stroke-opacity': opacity,
    });
  }
  const length = Math.hypot(b.x - a.x, b.y - a.y) || 1;
  let nx = -(b.y - a.y) / length;
  let ny = (b.x - a.x) / length;
  const light = lightingDirection(appearance);
  if (nx * light.x + ny * light.y < 0) {
    nx *= -1;
    ny *= -1;
  }
  const shift = width * 0.2;
  const gloss = appearance.surfaceFinish === 'glossy' ? appearance.specularLight : 0.08;
  return svg('g', { 'data-surface': 'cylinder' }, [
    svg('line', {
      x1: a.x, y1: a.y, x2: b.x, y2: b.y,
      stroke: scaleColor(color, clamp(appearance.ambientLight * 0.75, 0.16, 0.62)),
      'stroke-width': width * 1.35,
      'stroke-linecap': 'round',
      'stroke-opacity': opacity,
    }),
    svg('line', {
      x1: a.x, y1: a.y, x2: b.x, y2: b.y,
      stroke: scaleColor(color, clamp(appearance.ambientLight + appearance.diffuseLight * 0.58, 0.2, 1.1)),
      'stroke-width': width,
      'stroke-linecap': 'round',
      'stroke-opacity': opacity,
    }),
    svg('line', {
      x1: a.x + nx * shift, y1: a.y + ny * shift,
      x2: b.x + nx * shift, y2: b.y + ny * shift,
      stroke: '#ffffff',
      'stroke-width': Math.max(width * 0.18, 0.12),
      'stroke-linecap': 'round',
      'stroke-opacity': opacity * gloss * 0.72,
    }),
  ]);
}

export function renderScene(scene, camera, frame, appearance) {
  const root = svg('svg', {
    viewBox: `${-VIEW} ${-VIEW} ${2 * VIEW} ${2 * VIEW}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': 'Crystal structure',
  });
  const defs = svg('defs');
  root.append(defs);

  const project = (point) => projectPoint(camera, point);
  const items = [];
  const gradientCache = new Map();
  const projectedAtoms = (appearance.showAtoms ? scene.atoms : []).map((atom) => ({
    atom,
    point: project(atom.position),
  }));
  const atomDepths = projectedAtoms.map(({ point }) => point.depth);
  const minimumAtomDepth = atomDepths.length ? Math.min(...atomDepths) : 0;
  const atomDepthSpan = atomDepths.length ? Math.max(...atomDepths) - minimumAtomDepth : 0;

  for (const edge of appearance.showCells ? scene.cell_edges : []) {
    const a = project(edge[0]);
    const b = project(edge[1]);
    items.push({
      depth: (a.depth + b.depth) / 2 - 1e3, // cell edges sit behind matter at equal depth
      node: svg('line', {
        x1: a.x, y1: a.y, x2: b.x, y2: b.y,
        stroke: appearance.cellColor,
        'stroke-width': 0.45 * appearance.cellWidth,
        'stroke-opacity': appearance.cellOpacity,
      }),
    });
  }

  for (const plane of appearance.showPlanes ? scene.planes : []) {
    const points = plane.vertices.map(project);
    if (points.length < 3) continue;
    const polygon = svg('polygon', {
      points: points.map((point) => `${point.x},${point.y}`).join(' '),
      fill: appearance.planeColor,
      'fill-opacity': appearance.planeOpacity,
      stroke: appearance.planeColor,
      'stroke-width': 0.5 * appearance.cellWidth,
    });
    const depth = points.reduce((sum, point) => sum + point.depth, 0) / points.length;
    items.push({ depth, node: polygon, row: { Plane: plane.label ?? 'plane' } });
    if (plane.label && appearance.showLabels) {
      const centroid = points.reduce(
        (acc, point) => ({ x: acc.x + point.x / points.length, y: acc.y + point.y / points.length }),
        { x: 0, y: 0 },
      );
      items.push({
        depth: depth + 0.5,
        node: svg('text', {
          x: centroid.x, y: centroid.y,
          'text-anchor': 'middle',
          'font-size': 4.5 * appearance.annotationScale,
          fill: appearance.planeColor,
          'paint-order': 'stroke',
          stroke: 'var(--bg-raised)',
          'stroke-width': 1.2,
          text: plane.label,
        }),
      });
    }
  }

  for (const bond of appearance.showBonds ? scene.bonds : []) {
    const a = project(bond.start);
    const b = project(bond.end);
    const width = 1.1 * appearance.bondWidth;
    items.push({
      depth: (a.depth + b.depth) / 2,
      node: bondGlyph(
        a,
        b,
        appearance.bondColor,
        width,
        appearance.bondOpacity,
        appearance,
      ),
      row: { Bond: bond.species, 'Length / Å': formatNumber(bond.length, 4) },
    });
  }

  for (const { atom, point } of projectedAtoms) {
    const radius = atom.radius * camera.scale * camera.zoom * appearance.atomScale;
    const color = appearance.speciesColors[atom.species] ?? atom.color;
    const nearness = atomDepthSpan > 1e-12
      ? (point.depth - minimumAtomDepth) / atomDepthSpan
      : 1;
    const depthFade = appearance.depthCue * (1 - nearness);
    const depthOpacity = 1 - depthFade;
    const light = lightingDirection(appearance);
    const shadowOpacity = appearance.surfaceFinish === 'flat'
      ? 0
      : (1 - depthFade * 0.6) * (0.055 + appearance.diffuseLight * 0.055);
    const shadow = svg('circle', {
      cx: point.x - light.x * radius * 0.09,
      cy: point.y - light.y * radius * 0.09,
      r: radius * 1.015,
      fill: '#000000',
      'fill-opacity': shadowOpacity,
    });
    const sphere = svg('circle', {
      cx: point.x, cy: point.y, r: radius,
      fill: spherePaint(defs, gradientCache, color, appearance),
      'fill-opacity': 1,
      stroke: scaleColor(color, clamp(appearance.ambientLight * 0.7, 0.12, 0.55)),
      'stroke-opacity': appearance.surfaceFinish === 'flat' ? 0.55 : 0.82,
      'stroke-width': appearance.surfaceFinish === 'flat' ? 0.35 : 0.28,
      'data-surface': appearance.surfaceFinish === 'flat' ? 'disc' : 'sphere',
      'data-depth-opacity': depthOpacity.toFixed(3),
    });
    const depthVeil = svg('circle', {
      cx: point.x, cy: point.y, r: radius,
      fill: 'var(--bg-raised)',
      'fill-opacity': depthFade * 0.62,
      'pointer-events': 'none',
    });
    items.push({
      depth: point.depth,
      node: svg('g', {
        'data-atom': atom.species,
        opacity: appearance.atomOpacity,
      }, [shadow, sphere, depthVeil]),
      row: {
        Element: atom.species,
        Site: atom.label ?? '—',
        'x / Å': formatNumber(atom.position[0], 4),
        'y / Å': formatNumber(atom.position[1], 4),
        'z / Å': formatNumber(atom.position[2], 4),
        Occupancy: formatNumber(atom.occupancy, 3),
      },
    });
    if (atom.label && appearance.showLabels) {
      items.push({
        depth: point.depth + 1e-3,
        node: svg('text', {
          x: point.x, y: point.y + radius * 0.35,
          'text-anchor': 'middle',
          'font-size': Math.max(radius * 0.8, 2.5) * appearance.annotationScale,
          fill: 'var(--ink)',
          text: atom.label,
        }),
      });
    }
  }

  for (const direction of appearance.showDirections ? scene.directions : []) {
    const a = project(direction.start);
    const b = project(direction.end);
    items.push({
      depth: Math.max(a.depth, b.depth) + 1e3, // arrows stay legible in front
      node: svg('g', {}, [
        svg('line', {
          x1: a.x, y1: a.y, x2: b.x, y2: b.y,
          stroke: appearance.directionColor,
          'stroke-width': 1.4 * appearance.directionWidth,
          'stroke-opacity': appearance.directionOpacity,
        }),
        arrowHead(a, b, appearance.directionColor, appearance.directionOpacity),
        direction.label && appearance.showLabels
          ? svg('text', {
              x: b.x, y: b.y - 2,
              'text-anchor': 'middle',
              'font-size': 5 * appearance.annotationScale,
              fill: appearance.directionColor,
              'fill-opacity': appearance.directionOpacity,
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

  if (appearance.showGizmo) root.append(axisGizmo(scene.axes, camera));
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

function arrowHead(from, to, color, opacity = 1) {
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
    'fill-opacity': opacity,
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
  const pan = camera.pan ?? { x: 0, y: 0 };
  return { x: rotated[0] * scale + pan.x, y: -rotated[1] * scale + pan.y, depth: rotated[2] };
}

/** Invert the projection for the cursor readout, on the plane through the centre. */
function cameraToCrystal(camera, x, y) {
  const scale = camera.scale * camera.zoom;
  if (!scale) return null;
  const pan = camera.pan ?? { x: 0, y: 0 };
  const view = [(x - pan.x) / scale, -(y - pan.y) / scale, 0];
  const model = applyTranspose(camera.rotation, view);
  return {
    x: model[0] + camera.centre[0],
    y: model[1] + camera.centre[1],
    z: model[2] + camera.centre[2],
  };
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
