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

/** Scan extensions that hold HDF5 rather than text, and so travel base64. */
const HDF5_SUFFIXES = ['.oh5', '.h5'];

/**
 * A file's bytes as base64, for the JSON field that carries it to the server.
 *
 * Encoded in chunks rather than by spreading the whole array into
 * `String.fromCharCode`: a scan is megabytes, and one argument per byte
 * overflows the call stack somewhere in the low hundreds of thousands.
 */
async function readAsBase64(file) {
  const bytes = new Uint8Array(await file.arrayBuffer());
  const chunk = 0x8000;
  const pieces = [];
  for (let start = 0; start < bytes.length; start += chunk) {
    pieces.push(String.fromCharCode.apply(null, bytes.subarray(start, start + chunk)));
  }
  return btoa(pieces.join(''));
}

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

  // `scan` holds the user's own file once one is opened: {name, text} for a
  // text format, {name, data_base64} for an HDF5 one. It is kept beside the
  // form rather than in it, because a scan is a megabyte of text and a form
  // field is not where a megabyte of text belongs — the generated control for
  // it is hidden and the value supplied at call time, the same arrangement the
  // TEM panel uses for its picks.
  const state = { result: null, form: null, teaches: null, scan: null };

  const frame = plotFrame({ title: 'Orientation map', units: 'µm', digits: 2 });
  const legend = el('div.legend');
  const details = el('div');
  const formHost = el('div');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Draw the map',
    onclick: () => run(),
  });

  const scanStatus = el('p.field__help', {
    text: 'No scan open — the practice dataset chosen below is being analysed.',
  });

  const scanInput = el('input', {
    type: 'file',
    accept: '.ang,.ctf,.oh5,.h5',
    'aria-label': 'Open an EBSD scan file',
    onchange: (event) => openScan(event.target.files?.[0]),
  });

  const closeScanButton = el('button.button', {
    type: 'button',
    text: 'Close the scan',
    hidden: true,
    onclick: () => {
      state.scan = null;
      scanInput.value = '';
      closeScanButton.hidden = true;
      scanStatus.textContent =
        'No scan open — the practice dataset chosen below is being analysed.';
      run();
    },
  });

  context.rail.append(
    el('details.group', { open: true }, [
      el('summary', { text: 'Open a scan' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'An EDAX/TSL .ang, an Oxford/HKL .ctf, or an EDAX OIM .oh5 or .h5 — the last two ' +
            'being one HDF5 format under two extensions. It is read by the same importer a ' +
            'script would call, so the phases, the symmetry and the quality channels come from ' +
            'its own header. While one is open it replaces the practice dataset, and every ' +
            'control below means exactly what it means for a practice map.',
        }),
        el('p.field__help', {
          text:
            'An HDF5 scan saved with its diffraction patterns is far larger than a request can ' +
            'carry — export it without the patterns, or read it with pytex.adapters.read_scan ' +
            'in a script.',
        }),
        scanInput,
        scanStatus,
        closeScanButton,
      ]),
    ]),
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
    // The scan travels beside the form, so its generated control is not shown.
    for (const field of state.form.element.querySelectorAll('.field')) {
      if (field.querySelector('[id^="ctl-scan_file-"]')) field.hidden = true;
    }
  }

  /**
   * Read a scan file in the browser and analyse it.
   *
   * The contents are sent with the next call rather than uploaded to a store:
   * there is no store, the server keeps nothing between requests, and a scan
   * that failed to parse should leave nothing behind to clean up.
   *
   * A .ang or .ctf is text and travels as text. An .oh5 or .h5 is HDF5, so it
   * travels base64-encoded in the same field — one request path for both, at
   * the cost of the 4/3 inflation base64 carries.
   */
  async function openScan(file) {
    if (!file) return;
    scanStatus.textContent = `Reading ${file.name}…`;
    try {
      state.scan = HDF5_SUFFIXES.some((suffix) => file.name.toLowerCase().endsWith(suffix))
        ? { name: file.name, data_base64: await readAsBase64(file) }
        : { name: file.name, text: await file.text() };
      closeScanButton.hidden = false;
      scanStatus.textContent = `${file.name} — ${formatNumber(file.size / 1024, 0)} kB open.`;
      await run();
    } catch (error) {
      state.scan = null;
      closeScanButton.hidden = true;
      scanStatus.textContent = `${file.name} could not be read in the browser.`;
      context.showError(error);
    }
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
      const request = { ...state.form.values() };
      if (state.scan) request.scan_file = state.scan;
      state.result = await call(operation.id, request);
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
    } catch (error) {
      // The scan's own control is hidden — the file travels beside the form —
      // so an error naming it has nowhere on the form to land. It belongs on
      // the line under the file button, where the file's name is, and it is
      // worth a toast too: a file that could not be read is the one failure
      // where the panel goes on drawing something else.
      if (error?.field === 'scan_file') {
        scanStatus.textContent = error.message;
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
