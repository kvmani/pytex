/**
 * The variants panel: the child orientations one parent grain produces.
 *
 * Two views of the same variant set, chosen by which question is being asked.
 *
 * The **pole figure** shows where each variant puts a chosen child plane, in
 * the parent frame, so it is directly comparable with a measured pole figure of
 * the same grain. The reason to draw it here rather than export it is the same
 * reason the diffraction panel exists: seventy-two poles from twenty-four
 * variants all look alike, and the only question ever asked of one is "which
 * variant is *that*?". Every pole therefore carries its full row.
 *
 * The **spectrum** is a table, and gets the plain result view. Its content is
 * numbers, not geometry, and a histogram of ten discrete values would be a
 * worse presentation of them than the table already is.
 *
 * Colour carries packet, not variant. Twenty-four hues are twenty-four hues
 * whatever order they are in, and nobody reads a legend that long; four are
 * distinguishable at a glance, and the packet is the grouping that means
 * something in a micrograph.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult, saveBlob } from '../core/result.js';
import { call } from '../core/api.js';

export const panel = {
  id: 'variants',
  title: 'Variants',
  tagline: 'The child orientations one parent grain produces, and how they differ.',
};

/** Half-width of the drawing area. The unit disc is drawn at radius `VIEW`. */
const VIEW = 100;

/**
 * Packet colours.
 *
 * A small fixed set rather than a generated wheel, because packets are few (4
 * or 6) and are compared with each other constantly; fixed hues mean packet 1
 * is the same colour in every figure a reader puts side by side. Chosen to stay
 * separable for the common forms of colour vision deficiency, and to keep
 * roughly equal lightness so no packet reads as more important than another.
 */
const PACKET_COLORS = [
  'hsl(212 72% 52%)',
  'hsl(28 82% 50%)',
  'hsl(158 62% 38%)',
  'hsl(322 60% 52%)',
  'hsl(266 58% 58%)',
  'hsl(52 72% 42%)',
];

/** The parent's own poles: achromatic, because they are the reference. */
const PARENT_COLOR = 'var(--ink)';

function packetColor(packet) {
  return PACKET_COLORS[(packet - 1) % PACKET_COLORS.length];
}

/**
 * Radius on the unit disc of a pole at `polarDeg` from the projection axis.
 *
 * The service normalises both projections so the rim is radius 1, which is what
 * makes these two closed forms the whole story:
 *
 *   stereographic  r = tan(θ/2)
 *   equal area     r = √2·sin(θ/2)
 *
 * Both are 0 at the centre and 1 at θ = 90°, and differ everywhere between.
 */
function radiusFromPolar(polarDeg, projection) {
  const theta = (polarDeg * Math.PI) / 180;
  return projection === 'equal_area'
    ? Math.SQRT2 * Math.sin(theta / 2)
    : Math.tan(theta / 2);
}

/** The inverse of `radiusFromPolar`, in radians. */
function polarFromRadius(radius, projection) {
  return projection === 'equal_area'
    ? 2 * Math.asin(Math.min(radius / Math.SQRT2, 1))
    : 2 * Math.atan(radius);
}

/**
 * The two views, in the order the picker offers them.
 *
 * Named rather than filtered by panel: `variants.render` also belongs to this
 * panel and is not a view — it is the Figure button. Listing the views is the
 * only way for the panel to say which of its operations a user chooses between,
 * and a new operation added to the panel is then a deliberate addition here
 * rather than a surprise entry in a dropdown.
 */
const VIEWS = ['variants.pole_figure', 'variants.intervariant_misorientations'];

export function mount(context) {
  const operations = VIEWS.map((id) =>
    context.manifest.operations.find((entry) => entry.id === id),
  ).filter(Boolean);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = {
    operation: operations[0],
    result: null,
    teaches: null,
    form: null,
    hiddenPackets: new Set(),
    // Which result the legend was built for, so a redraw that only changes
    // which packets are shown updates the buttons instead of replacing them.
    legendFor: null,
  };

  // SVG first, and so the default: a pole figure is line art and a few dozen
  // markers, where vector is both smaller and more useful. This is the opposite
  // of the crystal viewer, where a sphere mesh makes SVG the expensive choice —
  // and the order here matches what `variants.render` declares and what its
  // help text tells the user to choose, which is the point.
  const formatSelect = el(
    'select.button',
    { 'aria-label': 'Figure format', title: 'File format for the published figure' },
    [
      el('option', { value: 'svg', text: 'SVG' }),
      el('option', { value: 'png', text: 'PNG 600 dpi' }),
    ],
  );

  const publishButton = el('button.button', {
    type: 'button',
    text: 'Figure',
    title: 'Render this figure through the publication renderer',
    onclick: () => publish(),
  });

  const frame = plotFrame({
    title: 'Pole figure',
    toolbar: [formatSelect, publishButton],
  });

  const legend = el('div.legend');
  const details = el('div');
  const formHost = el('div');

  const operationSelect = el(
    'select',
    {
      'aria-label': 'View',
      onchange: () => {
        state.operation = operations.find((entry) => entry.id === operationSelect.value);
        state.teaches = null;
        renderControls();
        run();
      },
    },
    operations.map((entry) =>
      el('option', { value: entry.id, text: entry.title, title: entry.summary }),
    ),
  );

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Show variants',
    onclick: () => run(),
  });

  context.rail.append(
    el('div.field', {}, [
      el('label.field__label', { text: 'View' }),
      operationSelect,
      el('p.field__help', {
        text: 'The pole figure shows where the variants point; the spectrum shows how they differ from each other.',
      }),
    ]),
    formHost,
    runButton,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'Run the first two in order: 24 variants in 4 packets, then the same figure with half as many.',
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

  function isPoleFigure() {
    return state.operation.id === 'variants.pole_figure';
  }

  function renderControls(initial = {}) {
    state.form = buildForm(state.operation, { initial });
    formHost.replaceChildren(state.form.element);
  }

  function loadExample(example) {
    state.operation = operations.find((entry) => entry.id === example.operation);
    operationSelect.value = state.operation.id;
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Computing…';
    state.form.clearErrors();
    try {
      const result = await call(state.operation.id, state.form.values());
      state.result = result;
      state.hiddenPackets = new Set();
      draw();
      renderResult(details, result, { teaches: state.teaches });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Show variants';
    }
  }

  function draw() {
    // The spectrum is a table. Hiding the frame rather than drawing an empty
    // one keeps the stage from carrying a picture of nothing, which reads as a
    // failed render rather than as a view that has no picture.
    const showFigure = isPoleFigure();
    frame.element.hidden = !showFigure;
    legend.hidden = !showFigure;
    if (!showFigure) return;

    const data = state.result.data;
    frame.configure({
      // Screen to data: the disc is drawn at radius VIEW and the y axis points
      // up in the projection and down on screen.
      toData: (x, y) => ({ x: x / VIEW, y: -y / VIEW }),
      formatCursor: (point) => {
        const radius = Math.hypot(point.x, point.y);
        if (radius > 1.0001) return 'outside the projection';
        // Inverting the projection is what makes the readout an instrument
        // rather than a decoration: the reader wants the angle from the
        // projection axis, not the position on the page. Both inverses are of
        // the *normalised* radius the service returns, where the rim is 1.
        const polar = polarFromRadius(radius, data.projection);
        const azimuth = ((Math.atan2(point.y, point.x) * 180) / Math.PI + 360) % 360;
        return (
          `polar ${formatNumber((polar * 180) / Math.PI, 2)}° · ` +
          `azimuth ${formatNumber(azimuth, 2)}°`
        );
      },
    });

    frame.setContent(renderPoleFigure(data, { frame, hidden: state.hiddenPackets }));
    // The legend is built once per result and updated in place thereafter.
    // Rebuilding it on every redraw destroys the button that was just pressed,
    // and the browser moves focus to the body — so a keyboard user who toggles
    // a packet loses their place and has to tab back through the whole page.
    if (state.legendFor !== data) buildLegend(data);
    else updateLegend();
    const shown = data.poles.filter((pole) => !state.hiddenPackets.has(pole.packet)).length;
    frame.setStatus(
      `${shown} of ${data.poles.length} variant poles · ${data.variant_count} variants in ` +
        `${data.packet_count} packets · ${
          data.projection === 'equal_area' ? 'equal-area' : 'stereographic'
        } projection · hover a pole for its variant, packet and indices`,
    );
  }

  /** Build the legend for a new result. Called once per result, not per redraw. */
  function buildLegend(data) {
    state.legendFor = data;
    const entries = Object.entries(data.packet_sizes)
      .map(([packet, size]) => ({ packet: Number(packet), size }))
      .sort((left, right) => left.packet - right.packet);
    legend.replaceChildren(
      ...entries.map(({ packet, size }) =>
        el(
          'button.legend__item',
          {
            type: 'button',
            dataset: { packet: String(packet) },
            onclick: () => {
              if (state.hiddenPackets.has(packet)) state.hiddenPackets.delete(packet);
              else state.hiddenPackets.add(packet);
              draw();
            },
          },
          [
            el('span.legend__swatch', { style: `background:${packetColor(packet)}` }),
            el('span', { text: `Packet ${packet} (${size} variants)` }),
          ],
        ),
      ),
      data.parent_poles.length
        ? el('span.legend__item', {}, [
            el('span.legend__swatch', { style: `background:${PARENT_COLOR}` }),
            el('span', { text: 'Parent' }),
          ])
        : null,
    );
    updateLegend();
  }

  /** Reflect the hidden set onto the existing buttons, without replacing them. */
  function updateLegend() {
    for (const button of legend.querySelectorAll('button.legend__item')) {
      const hidden = state.hiddenPackets.has(Number(button.dataset.packet));
      button.setAttribute('aria-pressed', String(!hidden));
      button.title = hidden ? 'Show this packet' : 'Hide this packet';
    }
  }

  /**
   * Render the same poles through the Python renderer and save the file.
   *
   * The request is the pole figure's own inputs plus a format, so the published
   * figure is of exactly what is on screen — including the projection and the
   * packet plane, which are the two settings it would be easiest to publish a
   * figure that quietly disagrees about.
   */
  async function publish() {
    if (!state.result || !isPoleFigure()) return;
    // The button is held by name rather than found by selector: the format
    // select carries `.button` too, for its styling, so `.plot__toolbar .button`
    // matches the select first and disables the wrong control.
    publishButton.disabled = true;
    publishButton.textContent = 'Rendering…';
    try {
      const rendered = await call('variants.render', {
        ...state.result.inputs,
        format: formatSelect.value,
      });
      const { image, format, encoding } = rendered.data;
      const mime = format === 'svg' ? 'image/svg+xml' : 'image/png';
      const filename = `variant-pole-figure.${format}`;
      if (encoding === 'base64') {
        // Bytes arrive base64-encoded because JSON has none. Only the decoding
        // happens here; the blob goes out through the one save path both shells
        // share, which is what keeps the desktop shell able to write it.
        const binary = atob(image);
        const bytes = new Uint8Array(binary.length);
        for (let index = 0; index < binary.length; index += 1) {
          bytes[index] = binary.charCodeAt(index);
        }
        await saveBlob(filename, new Blob([bytes], { type: mime }));
      } else {
        await saveBlob(filename, new Blob([image], { type: mime }));
      }
    } catch (error) {
      context.showError(error);
    } finally {
      publishButton.disabled = false;
      publishButton.textContent = 'Figure';
    }
  }

  renderControls();
  context.stage.append(frame.element, legend, details);
  if (examples.length) loadExample(examples[0]);

  return { help: () => state.operation };
}

/**
 * Draw the projection.
 *
 * The net is drawn first and faintly: small circles every 30° of polar angle
 * and radii every 30° of azimuth, so a reader can place a pole to within a few
 * degrees without measuring. Poles are drawn on top, parent last, because the
 * parent poles are the reference and must never be buried.
 */
function renderPoleFigure(data, { frame, hidden }) {
  const root = svg('svg', {
    viewBox: `${-VIEW * 1.12} ${-VIEW * 1.12} ${2.24 * VIEW} ${2.24 * VIEW}`,
    preserveAspectRatio: 'xMidYMid meet',
    'aria-label': 'Variant pole figure',
  });

  // Small circles of constant polar angle, at the radius the active projection
  // actually puts them — a grid drawn with the wrong projection is worse than
  // no grid, because it invites measurement.
  for (const polar of [30, 60]) {
    root.append(
      svg('circle', {
        cx: 0,
        cy: 0,
        r: radiusFromPolar(polar, data.projection) * VIEW,
        fill: 'none',
        stroke: 'currentColor',
        'stroke-opacity': 0.16,
        'stroke-width': 0.4,
        'stroke-dasharray': '2 3',
      }),
    );
  }
  for (let azimuth = 0; azimuth < 180; azimuth += 30) {
    const radians = (azimuth * Math.PI) / 180;
    root.append(
      svg('line', {
        x1: -VIEW * Math.cos(radians),
        y1: VIEW * Math.sin(radians),
        x2: VIEW * Math.cos(radians),
        y2: -VIEW * Math.sin(radians),
        stroke: 'currentColor',
        'stroke-opacity': 0.12,
        'stroke-width': 0.4,
      }),
    );
  }

  // The primitive circle, drawn solid: it is the equator, not decoration.
  root.append(
    svg('circle', {
      cx: 0,
      cy: 0,
      r: VIEW,
      fill: 'none',
      stroke: 'currentColor',
      'stroke-opacity': 0.55,
      'stroke-width': 0.9,
    }),
  );

  // Frame axes, named. Without them the figure has no orientation at all and
  // two figures of the same grain cannot be compared.
  for (const [label, x, y] of [
    ['x', VIEW, 0],
    ['y', 0, -VIEW],
  ]) {
    root.append(
      svg('text', {
        x: x * 1.06,
        y: y * 1.06 + 2,
        'text-anchor': 'middle',
        'font-size': 6,
        fill: 'currentColor',
        'fill-opacity': 0.55,
        text: label,
      }),
    );
  }

  const columns = data.columns;
  for (const pole of data.poles) {
    if (hidden.has(pole.packet)) continue;
    const node = svg('circle', {
      cx: pole.x * VIEW,
      cy: -pole.y * VIEW,
      r: 3.2,
      fill: packetColor(pole.packet),
      'fill-opacity': 0.85,
      stroke: 'var(--bg-raised)',
      'stroke-width': 0.6,
    });
    root.append(node);
    frame.hoverable(node, pole, columns);
  }

  // Parent poles as open squares: a different shape as well as a different
  // colour, so the distinction survives being printed in grey.
  for (const pole of data.parent_poles) {
    const size = 3.4;
    const node = svg('rect', {
      x: pole.x * VIEW - size,
      y: -pole.y * VIEW - size,
      width: 2 * size,
      height: 2 * size,
      fill: 'none',
      stroke: PARENT_COLOR,
      'stroke-width': 1.1,
    });
    root.append(node);
    frame.hoverable(node, pole, columns);
  }

  return root;
}
