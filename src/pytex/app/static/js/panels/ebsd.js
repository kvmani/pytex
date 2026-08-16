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
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'ebsd',
  title: 'EBSD',
  tagline: 'Orientation maps: IPF, grains, GROD, KAM, greyed by any channel, boundaries on top.',
};

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

export function mount(context) {
  const operation = context.manifest.operations.find((entry) => entry.id === 'ebsd.map');
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const state = { result: null, form: null, teaches: null };

  const frame = plotFrame({ title: 'Orientation map', units: 'µm', digits: 2 });
  const legend = el('div.legend');
  const details = el('div');
  const formHost = el('div');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Draw the map',
    onclick: () => run(),
  });

  context.rail.append(
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
  context.stage.append(frame.element, details);

  renderControls();
  // Draw immediately, as every other panel does. A workspace that opens empty
  // and waits to be pressed teaches nothing and looks broken; the first map is
  // also the fastest way to see what the four controls above it do.
  run();

  function renderControls(initial = {}) {
    state.form = buildForm(operation, { initial });
    formHost.replaceChildren(state.form.element);
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Drawing…';
    state.form.clearErrors();
    try {
      state.result = await call(operation.id, state.form.values());
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Draw the map';
    }
  }

  function draw() {
    const data = state.result.data;
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
    frame.configure({
      toData: (x, y) => ({ x: x / scale - step / 2, y: y / scale - step / 2 }),
      formatCursor: (point) => `${formatNumber(point.x, 2)}, ${formatNumber(point.y, 2)} µm`,
    });

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
