/**
 * The EBSD Kikuchi simulator: what the camera would see, before it sees it.
 *
 * Every other view in this workspace starts from a scan that has already been
 * indexed. This one starts from the lattice and the geometry, and produces the
 * pattern an indexing routine would have been given — which is what makes it
 * the place to answer the questions a scan cannot: what does this orientation
 * look like, is the camera geometry the one I think it is, and would these two
 * phases even be distinguishable on this screen.
 *
 * Drawn the way a pattern is read
 * -------------------------------
 * A band is the *gap between its two edges*, not a line, so the edges are drawn
 * and the space between them is filled faintly; the plane's own trace is dashed
 * between them because it is a construction rather than something the camera
 * records. Every band is named along itself — the convention `core/kikuchilabel.js`
 * owns — because on a screen where a dozen bands cross, a horizontal caption
 * belongs to whichever band is nearest the text.
 *
 * The zone axes are marked where the bands meet, sized by how many meet there,
 * since that count is the symmetry a reader recognises the pattern by. The
 * pattern centre is marked too: it is the one point on the screen that is a
 * statement about the *camera* rather than about the crystal, and half of EBSD
 * calibration is about where it is.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import { bandLabelNode, labelAngleDeg } from '../core/kikuchilabel.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'ebsd_kikuchi',
  title: 'Kikuchi simulator',
  tagline: 'The pattern a tilted specimen throws onto the EBSD camera, from the lattice.',
};

/* The screen is dark with bright bands whatever the interface theme, because
 * that is what a phosphor screen looks like; inverting it for a light theme
 * would make a simulated pattern look unlike the thing it is a simulation of. */
const SCREEN_COLOUR = '#05070d';
const BAND_COLOUR = '#dbeafe';
const TRACE_COLOUR = '#7dd3fc';
const AXIS_COLOUR = '#fbbf24';
const CENTRE_COLOUR = '#f87171';
const HALO_COLOUR = '#05070d';

export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'ebsd.simulate_kikuchi_pattern',
  );
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const state = {
    result: null,
    form: null,
    teaches: null,
    showLabels: true,
    showAxes: true,
  };

  const labelButton = toggle('Indices', 'Name each band along the band', true, (on) => {
    state.showLabels = on;
    draw();
  });
  const axisButton = toggle('Zone axes', 'Mark where the bands meet', true, (on) => {
    state.showAxes = on;
    draw();
  });

  const frame = plotFrame({
    title: 'Simulated EBSD pattern',
    readout: true,
    toolbar: [labelButton, axisButton],
  });

  const details = el('div');
  const formHost = el('div');
  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Simulate the pattern',
    onclick: () => run(),
  });

  context.rail.append(
    el('details.group', { open: true }, [
      el('summary', { text: 'The crystal, and the camera that sees it' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'The beam is the laboratory z axis and the stage tilt axis is x; the camera sits '
            + 'on −y. Tilting the stage tips the specimen normal towards the camera, and that '
            + 'normal then falls at gnomonic radius tan(90° − tilt + elevation) from the '
            + 'pattern centre — about 0.36 at the usual 70°.',
        }),
        formHost,
        runButton,
      ]),
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
  );

  context.stage.append(frame.element, details);

  function toggle(text, title, pressed, onChange) {
    return el('button.button', {
      type: 'button',
      text,
      title,
      'aria-pressed': String(pressed),
      onclick: (event) => {
        const next = event.currentTarget.getAttribute('aria-pressed') !== 'true';
        event.currentTarget.setAttribute('aria-pressed', String(next));
        onChange(next);
      },
    });
  }

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
    runButton.textContent = 'Simulating…';
    state.form.clearErrors();
    try {
      state.result = await call(operation.id, state.form.values());
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
      state.teaches = null;
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Simulate the pattern';
    }
  }

  function draw() {
    const data = state.result?.data;
    if (!data) {
      frame.setContent(
        el('div.stage__placeholder', { text: 'Choose a phase and a geometry to simulate.' }),
      );
      return;
    }
    const width = data.width_px;
    const height = data.height_px;
    const outer = svg('svg', {
      viewBox: `0 0 ${width} ${height}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': 'Simulated EBSD Kikuchi pattern',
    });
    // Clipped to the screen: a viewBox is a coordinate system rather than a
    // boundary, and a band edge is a conic that runs well past the picture.
    const clipId = 'ebsdkikuchi-clip';
    const root = svg('g', { 'clip-path': `url(#${clipId})` });
    outer.append(
      svg('defs', {}, [
        svg('clipPath', { id: clipId }, [svg('rect', { x: 0, y: 0, width, height })]),
      ]),
      root,
    );
    root.append(svg('rect', { x: 0, y: 0, width, height, fill: SCREEN_COLOUR }));

    const weight = Math.max(width, height) / 900;
    const font = Math.max(width, height) / 44;
    data.bands.forEach((band, index) => {
      drawBand(root, band, {
        index,
        weight,
        font,
        width,
        height,
        hoverable: (node, row) => frame.hoverable(node, row),
      });
    });
    if (state.showAxes) drawZoneAxes(root, data, { weight, font });
    drawPatternCentre(root, data, { weight, font });

    frame.configure({
      toData: (x, y) => ({ x, y }),
      formatCursor: (point) => {
        // Gnomonic radius, because that is the coordinate the geometry is
        // stated in: one unit is one camera distance, so the radius is the
        // tangent of the angle from the camera axis.
        const [cx, cy] = data.pattern_centre_px;
        const radius = Math.hypot(point.x - cx, point.y - cy);
        const gnomonic = (radius * data.pixel_size_mm) / data.camera_length_mm;
        return (
          `${formatNumber(radius, 1)} px from the pattern centre · `
          + `gnomonic ${formatNumber(gnomonic, 3)} · `
          + `${formatNumber((Math.atan(gnomonic) * 180) / Math.PI, 2)}° off the camera axis`
        );
      },
    });
    frame.setContent(outer);
    frame.setReadout(geometryCard(data));
    frame.setStatus(statusLine(data));
  }

  function drawBand(root, band, { index, weight, font, width, height, hoverable }) {
    const prominence = 0.28 + 0.6 * Math.min(1, Number(band.intensity) || 0);
    const [narrow, far] = band.edges;

    // The band itself is the gap between the edges. It is filled only when both
    // edges survived the clip as one run each: anything else and the polygon
    // would close across a stretch of screen the band never covered.
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
      if (hoverable) {
        hoverable(node, {
          Plane: band.label,
          'd / Å': band.d_angstrom,
          'Band width / deg': band.width_deg,
          'Relative intensity': band.intensity,
        });
      }
    }

    if (!state.showLabels) return;
    const anchor = labelAnchor(band.centre, { index, width, height });
    if (!anchor) return;
    root.append(
      bandLabelNode({
        x: anchor.x,
        y: anchor.y,
        angleDeg: labelAngleDeg(anchor.before, anchor.point),
        text: band.label,
        fontSize: font,
        colour: TRACE_COLOUR,
        haloColour: HALO_COLOUR,
        haloWidth: font / 5,
      }),
    );
  }

  /**
   * Where a band's name goes: on the band, inside the picture.
   *
   * Only the samples that are actually on the screen are candidates. The traces
   * arrive clipped to a generous margin — a Kossel conic has to be sampled well
   * past the frame for its shape to survive — so taking a fraction along the
   * whole run would put half the names outside the picture, where the clip eats
   * them silently.
   *
   * The fraction itself is stepped by band index. Parallel planes of successive
   * orders, (022) and (044), share one centre line, and their names would land
   * on top of one another at any fixed fraction.
   */
  function labelAnchor(runs, { index, width, height }) {
    const inset = 0.04 * Math.min(width, height);
    let longest = null;
    for (const run of runs) {
      const visible = run.filter(
        ([x, y]) => x > inset && x < width - inset && y > inset && y < height - inset,
      );
      if (visible.length >= 2 && (!longest || visible.length > longest.length)) longest = visible;
    }
    if (!longest) return null;
    const fraction = 0.2 + 0.12 * (index % 5);
    const at = Math.min(longest.length - 1, Math.max(1, Math.round(longest.length * fraction)));
    // A tangent measured over one sample is noise when the samples crowd; over a
    // twentieth of the visible run it is the direction of the band.
    const span = Math.max(1, Math.round(longest.length / 20));
    const before = longest[Math.max(0, at - span)];
    const point = longest[at];
    return { x: point[0], y: point[1], point, before };
  }

  function drawZoneAxes(root, data, { weight, font }) {
    for (const axis of data.zone_axes) {
      if (!axis.on_detector) continue;
      const radius = weight * (3 + 1.1 * Math.min(axis.order, 8));
      root.append(
        svg('circle', {
          cx: axis.x,
          cy: axis.y,
          r: radius,
          fill: 'none',
          stroke: AXIS_COLOUR,
          'stroke-opacity': 0.85,
          'stroke-width': weight * 1.2,
        }),
        svg('text', {
          x: axis.x,
          y: axis.y - radius - font * 0.35,
          'font-size': font * 0.95,
          fill: AXIS_COLOUR,
          stroke: HALO_COLOUR,
          'stroke-width': font / 6,
          'paint-order': 'stroke',
          'text-anchor': 'middle',
          text: axis.label,
        }),
        svg('title', { text: `${axis.label} — ${axis.order} bands meet here` }),
      );
    }
  }

  function drawPatternCentre(root, data, { weight, font }) {
    const [x, y] = data.pattern_centre_px;
    const arm = font * 0.8;
    root.append(
      svg('line', {
        x1: x - arm, y1: y, x2: x + arm, y2: y,
        stroke: CENTRE_COLOUR, 'stroke-width': weight * 1.4,
      }),
      svg('line', {
        x1: x, y1: y - arm, x2: x, y2: y + arm,
        stroke: CENTRE_COLOUR, 'stroke-width': weight * 1.4,
      }),
      svg('title', { text: 'The pattern centre: where the camera axis meets the screen' }),
    );
  }

  /**
   * The geometry beside the picture.
   *
   * Not decoration: a simulated pattern is only worth anything if the reader
   * knows which geometry it belongs to, and the specimen-normal radius is the
   * one line of arithmetic that checks the whole convention.
   */
  function geometryCard(data) {
    const inputs = state.result?.inputs ?? {};
    const [phi1, phi, phi2] = inputs.euler_deg ?? [0, 0, 0];
    const row = (label, value, title = '') =>
      el('tr', { title }, [el('th', { text: label }), el('td.measure__num', { text: value })]);
    return el('div.measure', {}, [
      el('div.measure__title', { text: 'Geometry' }),
      el('table.measure__table', {}, [
        el('tbody', {}, [
          row(
            'Bunge (φ₁, Φ, φ₂)',
            `${formatNumber(phi1, 1)}, ${formatNumber(phi, 1)}, ${formatNumber(phi2, 1)}°`,
            'Crystal-to-specimen, as an EBSD file reports it',
          ),
          row(
            'Stage tilt',
            `${formatNumber(inputs.sample_tilt_deg ?? 70, 1)}°`,
            'About the laboratory x axis, tipping the specimen normal towards the camera',
          ),
          row(
            'Camera distance z*',
            formatNumber(inputs.detector_distance ?? 0, 3),
            'Specimen-to-screen distance as a fraction of the screen width',
          ),
          row(
            'Field of view',
            `${formatNumber(data.field_of_view_deg, 1)}°`,
            'Across the width of the screen',
          ),
          row(
            'Specimen normal at',
            formatNumber(data.specimen_normal_gnomonic_radius, 3),
            'Gnomonic radius from the pattern centre: tan(90° − tilt + elevation)',
          ),
        ]),
      ]),
    ]);
  }

  function statusLine(data) {
    const axes = data.zone_axes.filter((axis) => axis.on_detector).length;
    return (
      `${data.bands.length} band(s) and ${axes} zone axis/axes on a `
      + `${data.width_px}×${data.height_px} screen · λ = `
      + `${formatNumber(data.wavelength_angstrom, 4)} Å · positions and widths exact, `
      + 'intensities kinematic · hover a band trace for its plane and width'
    );
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);
  else run();

  return { help: () => operation };
}
