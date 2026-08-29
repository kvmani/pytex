/**
 * EBSD: one map, four independent choices about what it shows.
 *
 * The map arrives as a raster at its native resolution — a crystal map on a
 * regular grid has one measurement per pixel, so there is nothing to
 * interpolate and any smoothing would invent microstructure. The boundary
 * network arrives separately as line segments in map coordinates and is drawn
 * as vectors on top, so it stays sharp under zoom and can be classified into
 * low- and high-angle by colour without touching the pixels underneath.
 *
 * The colour key is drawn here rather than sent as an image: for a scalar
 * colouring it is a gradient bar with the real range on it, and a colour bar
 * without numbers is decoration. A GROD map whose brightest pixel might be two
 * degrees or twenty cannot be read at all.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import { scanControls, withScan } from '../core/ebsdscan.js';
import { offer } from '../core/handoff.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'ebsd_map',
  title: 'IPF map',
  tagline: 'The orientation map: IPF colour, greyed by any channel, boundaries on top.',
};

/*
 * The same panel, opened on a different colouring.
 *
 * GROD and KAM are not other panels — they are this one with one control set
 * differently, and building them as copies would give the workspace three
 * places to fix the same bug. But they are also not *findable* behind a select
 * on a form: "show me the local misorientation" is a thing a user comes to the
 * workspace to do, and a sub-tab is where they will look for it. One
 * implementation, three doors into it, and every control still present inside.
 */
export function mapPanel({ id, title, tagline, colouring }) {
  return {
    panel: { id, title, tagline },
    mount: (context) => mount(context, { colouring }),
  };
}

export const grodPanel = mapPanel({
  id: 'ebsd_grod',
  title: 'GROD',
  tagline: 'How far each point has turned from its own grain: intragranular gradient.',
  colouring: 'grod',
});

export const kamPanel = mapPanel({
  id: 'ebsd_kam',
  title: 'KAM',
  tagline: 'How far each point has turned from its neighbours: local deformation.',
  colouring: 'kam',
});

/** Side of the SVG the map is drawn into, in user units. */
const SIZE = 700;

/**
 * How boundary character is drawn: high-angle black and heavier, low-angle red
 * and finer, as the literature draws them.
 *
 * Fixed colours rather than theme tokens, deliberately. These lines sit on top
 * of the map image, so they need contrast against *the map* — IPF colours, a
 * viridis ramp — and not against the page background. A token that inverted
 * with the theme would put a near-white boundary on a pale IPF grain.
 */
const BOUNDARY_STYLE = {
  high: { stroke: '#111111', width: 1.6 },
  low: { stroke: '#e03131', width: 0.9 },
};

/**
 * The two picked grains, and how they are drawn on the map.
 *
 * Fixed colours for the same reason the boundaries have them: these outlines
 * sit on top of IPF colours and a viridis ramp, so they need contrast against
 * *the map*. Each is drawn over a white halo, which is what makes an outline
 * legible on a pale grain and on a dark one without knowing which it landed on.
 *
 * Parent and child rather than first and second, because that is what the
 * relationship they are being picked for calls them, and the two are not
 * interchangeable: the parent is the phase that transformed.
 */
const PICK_SLOTS = [
  { role: 'parent', label: 'Parent grain', colour: '#f08c00' },
  { role: 'child', label: 'Child grain', colour: '#d6336c' },
];

export function mount(context, { colouring = 'ipf' } = {}) {
  const operation = context.manifest.operations.find((entry) => entry.id === 'ebsd.map');
  // By operation rather than by panel: the whole workspace declares `panel:
  // "ebsd"` in Python, because its operations belong to one subject, and the
  // sub-panels here are views of it rather than separate panels in the manifest.
  const examples = context.manifest.examples.filter(
    (entry) => entry.operation === operation.id,
  );

  // `scan` holds the user's own file once one is opened: {name, text} for a
  // text format, {name, data_base64} for an HDF5 one. It is kept beside the
  // form rather than in it, because a scan is a megabyte of text and a form
  // field is not where a megabyte of text belongs — the generated control for
  // it is hidden and the value supplied at call time, the same arrangement the
  // TEM panel uses for its picks.
  //
  // `picks` holds the grain ids chosen on the map, parent first, and `labels`
  // the per-pixel grain array the click resolves through — decoded once per
  // result rather than per click.
  const state = {
    result: null,
    form: null,
    teaches: null,
    picks: [],
    labels: null,
    picksLayer: null,
    geometry: null,
  };

  const frame = plotFrame({ title: 'Orientation map', units: 'µm', digits: 2 });
  const legend = el('div.legend');
  const picksHost = el('section.picks');
  const details = el('div');
  const formHost = el('div');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Draw the map',
    onclick: () => run(),
  });

  const scan = scanControls({ onChange: () => run(), showError: context.showError });

  context.rail.append(
    scan.element,
    formHost,
    runButton,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'Each dataset is a construction with a known answer, so the numbers on screen can be ' +
            'checked rather than trusted. The GROD and KAM examples are the same microstructure ' +
            'seen two ways.',
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
  frame.setControls(legend);
  context.stage.append(frame.element, picksHost, details);

  // On the frame rather than on the drawing: every redraw replaces the SVG, and
  // a handler owned by it would stop working the first time a control changed.
  // The frame swallows the click that ends a pan before it reaches here, so
  // dragging the map is not also picking whatever the drag finished on.
  frame.element.addEventListener('click', (event) => {
    if (!state.labels || !event.target.closest('svg')) return;
    const point = frame.pointerToData(event);
    if (!point || !(point.grain >= 0)) return;
    pick(point.grain);
  });

  // The sub-tab decides what the map opens on; the control is still there, so a
  // reader who arrived at GROD can switch to KAM without changing tabs.
  renderControls({ colouring });
  // Draw immediately, as every other panel does. A workspace that opens empty
  // and waits to be pressed teaches nothing and looks broken; the first map is
  // also the fastest way to see what the four controls above it do.
  run();

  function renderControls(initial = {}) {
    state.form = buildForm(operation, { initial });
    formHost.replaceChildren(state.form.element);
    // The scan travels beside the form, so its generated control is not shown.
    for (const field of state.form.element.querySelectorAll('.field')) {
      if (field.querySelector('[id^="ctl-scan_file-"]')) field.hidden = true;
    }
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  /* -------------------------------------------------------------- picking */

  /**
   * Choose a grain, or unchoose it.
   *
   * Two slots, filled parent then child. Clicking a grain that is already
   * chosen releases it; clicking a third grain when both slots are full starts
   * a *new* pair with that grain as the parent, rather than silently replacing
   * one of the two and leaving the user to work out which. A pair is a gesture
   * with a beginning, and starting over is what a third click means.
   */
  function pick(grainId) {
    const existing = state.picks.indexOf(grainId);
    if (existing >= 0) state.picks.splice(existing, 1);
    else if (state.picks.length >= PICK_SLOTS.length) state.picks = [grainId];
    else state.picks.push(grainId);
    renderPicks();
  }

  function grainRow(grainId) {
    return (state.result?.data?.grains ?? []).find((row) => row.grain_id === grainId) ?? null;
  }

  /** Redraw the outlines and rebuild the card under the map. */
  function renderPicks() {
    drawPickOutlines();
    picksHost.replaceChildren(...pickCard());
  }

  /**
   * The chosen grains outlined on the map itself.
   *
   * Outlined rather than tinted: a tint would change the colour the map is
   * *for*, and on an IPF map the colour is the measurement. The outline follows
   * the cells that carry the label, so what is highlighted is the grain the
   * segmentation found rather than a circle around where the click landed.
   */
  function drawPickOutlines() {
    if (!state.picksLayer || !state.labels || !state.geometry) return;
    const nodes = [];
    state.picks.forEach((grainId, index) => {
      const slot = PICK_SLOTS[index] ?? PICK_SLOTS[PICK_SLOTS.length - 1];
      const path = grainOutlinePath(state.labels, state.geometry, grainId);
      if (!path) return;
      nodes.push(
        // Six units of halo under three of colour, in a drawing 700 units
        // across: measured on screen rather than guessed, because at the
        // boundaries' own 1.6 the outline of a small grain is a sub-pixel line
        // and the highlight is invisible at the size the map is actually shown.
        svg('path', {
          d: path,
          fill: 'none',
          stroke: '#ffffff',
          'stroke-width': 6,
          'stroke-linecap': 'round',
          'stroke-opacity': 0.85,
        }),
        svg('path', {
          d: path,
          fill: 'none',
          stroke: slot.colour,
          'stroke-width': 3,
          'stroke-linecap': 'round',
        }),
      );
    });
    state.picksLayer.replaceChildren(...nodes);
  }

  /**
   * What has been picked, and what it costs to believe the answer.
   *
   * The grain orientation spread is on the card beside the angles, not left
   * behind in the table: two mean orientations fed into a relationship give a
   * residual of exactly zero for one pair, by construction, and the spread is
   * the only measure of what that zero conceals. It travels with the pair into
   * the Variants workspace for the same reason.
   */
  function pickCard() {
    if (!state.result) return [];
    const slots = PICK_SLOTS.map((slot, index) => {
      const grainId = state.picks[index];
      const row = grainId === undefined ? null : grainRow(grainId);
      return el('div.picks__slot', { dataset: { role: slot.role } }, [
        el('span.picks__swatch', { style: `background:${slot.colour}` }),
        el('div.picks__body', {}, [
          el('strong', { text: row ? `${slot.label}: grain ${row.grain_id}` : slot.label }),
          row
            ? el('span.picks__angles', {
                text:
                  `${formatNumber(row.mean_phi1_deg, 2)}, ${formatNumber(row.mean_Phi_deg, 2)}, ` +
                  `${formatNumber(row.mean_phi2_deg, 2)}° Bunge · ` +
                  `${row.phase_name ?? 'phase not named'}`,
              })
            : el('span.picks__angles', { text: 'Click a grain on the map.' }),
          row
            ? el('span.picks__spread', {
                title:
                  'Grain orientation spread: how far the grain’s own points sit from its ' +
                  'mean. The relationship is computed from the means, so this is the only ' +
                  'measure of what its zero residual conceals.',
                text:
                  `GOS ${formatNumber(row.grain_orientation_spread_deg, 3)}° ` +
                  `over ${row.size} points`,
              })
            : null,
        ]),
        row
          ? el('button.picks__release', {
              type: 'button',
              title: 'Release this grain',
              text: '×',
              onclick: () => pick(row.grain_id),
            })
          : null,
      ]);
    });

    const ready = state.picks.length === PICK_SLOTS.length;
    return [
      el('div.picks__slots', {}, slots),
      el('div.picks__actions', {}, [
        el('button.button.button--primary', {
          type: 'button',
          text: 'Send the pair to Variants',
          disabled: !ready,
          title: ready
            ? 'Open the relationship view with these two grains already in it'
            : 'Pick two grains on the map first',
          onclick: () => send(),
        }),
        el('button.button', {
          type: 'button',
          text: 'Swap',
          disabled: !ready,
          title: 'Exchange parent and child: which is which changes the relationship',
          onclick: () => {
            state.picks.reverse();
            renderPicks();
          },
        }),
        el('button.button', {
          type: 'button',
          text: 'Clear',
          disabled: !state.picks.length,
          onclick: () => {
            state.picks = [];
            renderPicks();
          },
        }),
      ]),
      el('p.field__help', {
        text:
          'Click a grain to pick it, and a second for the other side of the relationship. What ' +
          'travels across is each grain’s symmetry-aware mean orientation, its phase as ' +
          'the scan names it, and its spread — nothing is retyped, and nothing ' +
          'about the scatter is left behind.',
      }),
    ];
  }

  /**
   * Hand the pair to the Variants workspace and follow it there.
   *
   * The offer carries the numbers *and* their provenance — which map, which
   * dataset or file, at which segmentation threshold — because a pair of Euler
   * triples with no history is exactly the unattributable number this
   * repository refuses. The phases go as the scan names them, with the built-in
   * each name resolves to when it resolves to one; where it does not, the
   * receiving panel asks rather than assuming a lattice the measurement never
   * claimed.
   */
  function send() {
    const data = state.result?.data;
    if (!data || state.picks.length !== PICK_SLOTS.length) return;
    const builtins = data.phase_builtins ?? {};
    const grains = state.picks.map((grainId, index) => {
      const row = grainRow(grainId);
      return {
        role: PICK_SLOTS[index].role,
        grain_id: row.grain_id,
        phi1_deg: row.mean_phi1_deg,
        Phi_deg: row.mean_Phi_deg,
        phi2_deg: row.mean_phi2_deg,
        phase_name: row.phase_name ?? null,
        phase_builtin: builtins[row.phase_name] ?? null,
        grain_orientation_spread_deg: row.grain_orientation_spread_deg,
        size: row.size,
      };
    });
    offer('measured-pair', {
      euler_convention: 'bunge',
      grains,
      source: {
        panel: panel.id,
        map: data.dataset?.title ?? 'this map',
        grain_threshold_deg: state.result.inputs?.grain_threshold_deg ?? null,
      },
    });
    context.openPanel('variants');
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Drawing…';
    state.form.clearErrors();
    try {
      state.result = await call(operation.id, withScan(state.form.values()));
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
    } catch (error) {
      // The scan's own control is hidden — the file travels beside the form —
      // so an error naming it has nowhere on the form to land. It belongs on
      // the line under the file button, where the file's name is, and it is
      // worth a toast too: a file that could not be read is the one failure
      // where the panel goes on drawing something else.
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
      runButton.textContent = 'Draw the map';
    }
  }

  function draw() {
    const data = state.result.data;
    state.labels = decodeLabels(data.grain_ids);
    // A pick names a grain, and a redraw at a different threshold is a
    // different segmentation: the same number would then outline a different
    // region. Dropping the picks is the only honest answer.
    if (!state.labels) state.picks = [];
    const [, , width, height] = data.extent_um;
    // The image spans one step beyond the last point centre in each direction,
    // because a pixel is a measurement's area rather than its centre.
    const step = data.step_um;
    const spanX = width + step;
    const spanY = height + step;
    const scale = SIZE / Math.max(spanX, spanY);
    const drawWidth = spanX * scale;
    const drawHeight = spanY * scale;

    const nodes = [
      svg('image', {
        href: toDataUrl(data.image),
        x: 0,
        y: 0,
        width: drawWidth,
        height: drawHeight,
        // One pixel is one measurement. Smoothing would draw orientations
        // between points that were never measured.
        style: 'image-rendering: pixelated',
      }),
    ];

    // Under the boundaries, so a boundary is never hidden by a highlight, and
    // over the image, so the highlight is visible on the grain it names.
    state.picksLayer = svg('g', { 'data-role': 'picks' });
    nodes.push(state.picksLayer);

    for (const line of data.boundaries ?? []) {
      const character =
        line.misorientation_deg >= data.high_angle_threshold_deg ? 'high' : 'low';
      const style = BOUNDARY_STYLE[character];
      const node = svg('line', {
        x1: (line.x1 + step / 2) * scale,
        y1: (line.y1 + step / 2) * scale,
        x2: (line.x2 + step / 2) * scale,
        y2: (line.y2 + step / 2) * scale,
        stroke: style.stroke,
        'stroke-width': style.width,
        'stroke-linecap': 'square',
      });
      frame.hoverable(node, {
        Boundary: character === 'high' ? 'High-angle' : 'Low-angle',
        'Misorientation / °': line.misorientation_deg,
      });
      nodes.push(node);
    }

    frame.setContent(
      svg(
        'svg',
        { viewBox: `0 0 ${drawWidth} ${drawHeight}`, width: '100%', height: '100%' },
        nodes,
      ),
    );
    state.geometry = state.labels
      ? {
          cellWidth: drawWidth / state.labels.width,
          cellHeight: drawHeight / state.labels.height,
        }
      : null;
    frame.configure({
      // The grain travels with the position, so one lookup serves both the
      // readout and the click: the number under the cursor is the number a
      // click would pick, and they cannot disagree.
      toData: (x, y) => ({
        x: x / scale - step / 2,
        y: y / scale - step / 2,
        grain: labelAt(state.labels, state.geometry, x, y),
      }),
      formatCursor: (point) =>
        `${formatNumber(point.x, 2)}, ${formatNumber(point.y, 2)} µm` +
        (point.grain >= 0 ? ` · grain ${point.grain}` : ''),
    });
    renderPicks();

    const summary = data.boundary_summary ?? [];
    frame.setStatus(
      `${data.grain_count} grains · ${(data.boundaries ?? []).length} boundary segments drawn` +
        (summary.length
          ? ` · ${summary
              .map((row) => `${row.character} ${(row.fraction * 100).toFixed(0)}%`)
              .join(', ')}`
          : ''),
    );
    frame.setOverlay(
      el('div.plot__caption', {}, [
        el('strong', { text: data.dataset.title }),
        el('span', { text: `${(data.grid_shape ?? []).join(' × ')} points at ${step} µm step` }),
      ]),
    );
    renderLegend(data);
  }

  function renderLegend(data) {
    legend.replaceChildren();

    if (data.colour_scale) {
      legend.append(colourBar(data.colour_scale));
    } else if (data.colouring === 'ipf') {
      legend.append(
        el('p.field__help', {
          text:
            `Colour is the crystal direction along specimen ${data.ipf_direction}, folded into ` +
            'the fundamental sector. The colour key belongs to the point group: two maps of ' +
            'different symmetries are not colour-comparable.',
        }),
      );
    } else {
      legend.append(
        el('p.field__help', {
          text: 'One arbitrary hue per grain. The colours carry no orientation meaning.',
        }),
      );
    }

    if ((data.boundaries ?? []).length) {
      legend.append(
        el('div.legend__row', {}, [
          swatch(BOUNDARY_STYLE.high.stroke, `High-angle (≥${data.high_angle_threshold_deg}°)`),
          swatch(BOUNDARY_STYLE.low.stroke, `Low-angle (<${data.high_angle_threshold_deg}°)`),
        ]),
      );
    }
    if (data.modulate_by && data.modulate_by !== 'none') {
      legend.append(
        el('p.field__help', {
          text: `Brightness is modulated by ${data.modulate_by.replace(/_/g, ' ')}; darker means worse.`,
        }),
      );
    }
  }

  return { help: () => operation };
}

/* ------------------------------------------------------------------ picking */

/**
 * The per-pixel grain labels, decoded from the payload's base64 int32 array.
 *
 * A click has to resolve to the grain the *segmentation* found, and colour
 * cannot answer that: two grains of one orientation carry exactly the same
 * colour, and a colour-matching pick would join them into one silently. The
 * labels are the segmentation's own answer, aligned cell for cell with the
 * image.
 */
function decodeLabels(payload) {
  if (!payload?.data) return null;
  const binary = atob(payload.data);
  const bytes = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) bytes[index] = binary.charCodeAt(index);
  return {
    width: payload.width,
    height: payload.height,
    // Int32Array over the bytes, little-endian as NumPy wrote them. A DataView
    // read per cell would be correct too and an order of magnitude slower on a
    // map of a quarter of a million points.
    values: new Int32Array(bytes.buffer),
  };
}

/** The label at a point in the drawing's own units, or -1 where there is none. */
function labelAt(labels, geometry, x, y) {
  if (!labels || !geometry) return -1;
  const column = Math.floor(x / geometry.cellWidth);
  const row = Math.floor(y / geometry.cellHeight);
  if (column < 0 || row < 0 || column >= labels.width || row >= labels.height) return -1;
  return labels.values[row * labels.width + column];
}

/**
 * The outline of one grain, as a single SVG path.
 *
 * An edge is drawn wherever a cell carrying the label meets one that does not,
 * which gives the grain's true boundary including its holes — a twin lamella
 * enclosed in a grain outlines as two rings, and drawing a bounding box instead
 * would claim the lamella as part of it.
 *
 * One path rather than a segment element each: a grain in a fine map has
 * thousands of edges, and thousands of nodes is where this stops being
 * interactive.
 */
function grainOutlinePath(labels, geometry, grainId) {
  if (!labels || !geometry) return null;
  const { width, height, values } = labels;
  const { cellWidth, cellHeight } = geometry;
  const parts = [];
  const edge = (x1, y1, x2, y2) => {
    parts.push(`M${x1.toFixed(2)} ${y1.toFixed(2)}L${x2.toFixed(2)} ${y2.toFixed(2)}`);
  };
  for (let row = 0; row < height; row += 1) {
    for (let column = 0; column < width; column += 1) {
      if (values[row * width + column] !== grainId) continue;
      const left = column * cellWidth;
      const top = row * cellHeight;
      const right = left + cellWidth;
      const bottom = top + cellHeight;
      if (column === 0 || values[row * width + column - 1] !== grainId)
        edge(left, top, left, bottom);
      if (column === width - 1 || values[row * width + column + 1] !== grainId)
        edge(right, top, right, bottom);
      if (row === 0 || values[(row - 1) * width + column] !== grainId) edge(left, top, right, top);
      if (row === height - 1 || values[(row + 1) * width + column] !== grainId)
        edge(left, bottom, right, bottom);
    }
  }
  return parts.length ? parts.join('') : null;
}

/* ---------------------------------------------------------------- rendering */

/**
 * A gradient bar with the real range written on it.
 *
 * Built from the same colour maps the service used, sampled here so the bar and
 * the pixels cannot drift apart in appearance — and labelled with the numbers,
 * because an unlabelled colour bar tells a reader nothing about whether a
 * bright pixel is two degrees of misorientation or twenty.
 */
function colourBar(scale) {
  const stops = 32;
  const gradient = Array.from({ length: stops }, (_, index) => {
    const position = index / (stops - 1);
    const [r, g, b] = sampleColourMap(scale.colour_map, position);
    return `rgb(${r},${g},${b}) ${(position * 100).toFixed(1)}%`;
  }).join(', ');

  return el('div.legend__scale', {}, [
    el('span.legend__scale-label', {
      text: `${scale.label}${scale.units ? ` / ${scale.units}` : ''}`,
    }),
    el('span.legend__scale-bar', { style: `background: linear-gradient(to right, ${gradient})` }),
    el('span.legend__scale-ends', {
      text: `${formatNumber(scale.minimum, 3)} — ${formatNumber(scale.maximum, 3)}`,
    }),
  ]);
}

/** The same anchors `pytex.app.services.ebsd` interpolates, kept in step by hand. */
const COLOUR_ANCHORS = {
  viridis: [
    [68, 1, 84],
    [58, 82, 139],
    [32, 145, 140],
    [94, 201, 98],
    [253, 231, 37],
  ],
  magma: [
    [0, 0, 4],
    [81, 18, 124],
    [183, 55, 121],
    [251, 135, 97],
    [252, 253, 191],
  ],
};

function sampleColourMap(name, position) {
  const anchors = COLOUR_ANCHORS[name] ?? COLOUR_ANCHORS.viridis;
  const span = (anchors.length - 1) * Math.min(Math.max(position, 0), 1);
  const lower = Math.min(Math.floor(span), anchors.length - 2);
  const blend = span - lower;
  return anchors[lower].map((value, channel) =>
    Math.round(value + (anchors[lower + 1][channel] - value) * blend),
  );
}

function swatch(colour, text) {
  return el('span.legend__item', {}, [
    el('span.legend__swatch', { style: `background: ${colour}` }),
    el('span', { text }),
  ]);
}

/**
 * Decode the base64 RGB raster into a data URL through a canvas.
 *
 * The bytes are sent rather than a PNG so that Python needs no image encoder,
 * and a canvas is the one place a browser will turn raw pixels into something
 * an `<image>` element accepts.
 */
function toDataUrl(image) {
  const binary = atob(image.data);
  const canvas = document.createElement('canvas');
  canvas.width = image.width;
  canvas.height = image.height;
  const context = canvas.getContext('2d');
  const pixels = context.createImageData(image.width, image.height);
  for (let index = 0; index < image.width * image.height; index += 1) {
    const source = index * 3;
    const target = index * 4;
    pixels.data[target] = binary.charCodeAt(source);
    pixels.data[target + 1] = binary.charCodeAt(source + 1);
    pixels.data[target + 2] = binary.charCodeAt(source + 2);
    pixels.data[target + 3] = 255;
  }
  context.putImageData(pixels, 0, 0);
  return canvas.toDataURL('image/png');
}
