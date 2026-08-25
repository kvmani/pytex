/**
 * The SAED simulator: the pattern a crystal *would* give, before any plate exists.
 *
 * The solver next door answers "which zone axis is this exposure". This answers
 * the forward question, and the two are worth having side by side for the same
 * reason a microscopist keeps a book of patterns open at the column: the way to
 * be confident about an indexed answer is to see what the answer predicts and
 * compare it with what is on the screen.
 *
 * Three things this panel insists on:
 *
 * 1. **The orientation is stated, not implied.** Point the crystal by zone axis
 *    or by Bunge Euler angles; either way the readout says which zone is on the
 *    beam, what the roll about it is, and — for an orientation that is not
 *    exactly a zone axis — how many degrees off it the pattern shown really is.
 * 2. **The Kikuchi bands come from the same orientation as the spots.** They are
 *    fetched from the same overlay operation the solver uses, given the matrix
 *    this simulation reports, so the relation worth learning is visible on a
 *    pattern whose answer is known: the band for (hkl) is perpendicular to its
 *    own spot and exactly as wide as that spot is far from the beam.
 * 3. **The plate is drawn by the shared drawing.** `core/saedplot.js` paints
 *    this and the solver's practice plates, because they are the same
 *    calculation and a reader who learns one has learned the other.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber, svg } from '../core/dom.js';
import { plotFrame } from '../core/plotframe.js';
import { renderResult } from '../core/result.js';
import { drawKikuchiBands, drawSimulatedPattern } from '../core/saedplot.js';

export const panel = {
  id: 'tem_simulator',
  title: 'SAED Simulator',
  tagline: 'The pattern a crystal gives down a chosen zone axis, with its Kikuchi bands.',
};

export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'tem.simulate_saed',
  );
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);

  const state = {
    result: null,
    form: null,
    teaches: null,
    showLabels: true,
    showKikuchi: false,
    // The band overlay, and the request that produced it, so switching the
    // toggle back on does not fetch the same answer twice.
    kikuchi: null,
    kikuchiRequest: null,
  };

  const labelButton = el('button.button', {
    type: 'button',
    text: 'Indices',
    title: 'Write each reflection’s index beside it',
    'aria-pressed': 'true',
    onclick: (event) => {
      state.showLabels = !state.showLabels;
      event.currentTarget.setAttribute('aria-pressed', String(state.showLabels));
      draw();
    },
  });

  const kikuchiButton = el('button.button', {
    type: 'button',
    text: 'Kikuchi',
    title: 'Superimpose the bands this orientation predicts',
    'aria-pressed': 'false',
    onclick: (event) => {
      state.showKikuchi = !state.showKikuchi;
      event.currentTarget.setAttribute('aria-pressed', String(state.showKikuchi));
      if (state.showKikuchi) refreshKikuchi();
      else draw();
    },
  });

  const frame = plotFrame({
    title: 'Simulated pattern',
    // The plate is read with the pointer here as much as in the solver: the
    // question asked of a simulated pattern is "how far out is that spot, and
    // what spacing is that", so the readout belongs under the drawing.
    readout: true,
    toolbar: [labelButton, kikuchiButton],
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
      el('summary', { text: 'The crystal, and where it points' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'Choose a phase and a direction to look down it. The pattern is the exact ' +
            'zone-axis section: every reflection sits at the camera constant divided by its ' +
            'own d-spacing, which is the relation a real plate is calibrated by.',
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

  function renderControls(initial = {}) {
    state.form = buildForm(operation, { initial });
    formHost.replaceChildren(state.form.element);
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    // An example that asks for bands turns the toggle on, so what it teaches is
    // what appears rather than what appears after finding a button.
    state.showKikuchi = Boolean(example.request.show_kikuchi);
    kikuchiButton.setAttribute('aria-pressed', String(state.showKikuchi));
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Simulating…';
    state.form.clearErrors();
    try {
      state.result = await call(operation.id, state.form.values());
      // A new pattern is a new orientation, so any bands in hand are about a
      // pattern that is no longer on screen.
      state.kikuchi = null;
      state.kikuchiRequest = null;
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
      state.teaches = null;
      if (state.showKikuchi) refreshKikuchi();
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Simulate the pattern';
    }
  }

  /**
   * Fetch the bands for the orientation now on screen.
   *
   * Through the same operation the solver uses, given the matrix this
   * simulation reports. Nothing about the band geometry is computed here, so
   * the two panels cannot come to disagree about where a band falls.
   */
  async function refreshKikuchi() {
    const data = state.result?.data;
    if (!data) return;
    const pattern = data.pattern;
    const request = {
      phase: data.calibration.phase,
      orientation: { crystal_to_pattern: data.orientation.crystal_to_pattern },
      units: 'px',
      camera_constant_mm_angstrom: data.calibration.camera_constant_mm_angstrom,
      pixel_size_mm: data.calibration.pixel_size_mm,
      centre_x: pattern.centre_px[0],
      centre_y: pattern.centre_px[1],
      frame_width: pattern.width_px,
      frame_height: pattern.height_px,
      accelerating_voltage_kv: state.result.inputs.beam_energy_kev,
    };
    const key = JSON.stringify(request);
    if (key === state.kikuchiRequest && state.kikuchi) {
      draw();
      return;
    }
    try {
      state.kikuchi = await call('tem.kikuchi_overlay', request);
      state.kikuchiRequest = key;
    } catch (error) {
      // Bands that cannot be computed must not take the pattern away with them:
      // the spots are the result, and the overlay is an addition to it.
      state.kikuchi = null;
      state.showKikuchi = false;
      kikuchiButton.setAttribute('aria-pressed', 'false');
      context.showError(error);
    }
    draw();
  }

  function draw() {
    const data = state.result?.data;
    if (!data) {
      frame.setContent(
        el('div.stage__placeholder', { text: 'Choose a phase and a zone axis to simulate.' }),
      );
      return;
    }
    const pattern = data.pattern;
    const { width_px: width, height_px: height } = pattern;
    const outer = svg('svg', {
      viewBox: `0 0 ${width} ${height}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': 'Simulated diffraction pattern',
    });
    // Clipped to the plate for the same reason the solver's is: a viewBox is a
    // coordinate system rather than a boundary, and a band drawn beyond the
    // detector would be painted over the letterbox as though it were data.
    const clipId = 'saedsim-clip';
    const root = svg('g', { 'clip-path': `url(#${clipId})` });
    outer.append(
      svg('defs', {}, [
        svg('clipPath', { id: clipId }, [svg('rect', { x: 0, y: 0, width, height })]),
      ]),
      root,
    );

    drawSimulatedPattern(root, pattern, {
      labels: state.showLabels,
      hoverable: (node, row) => frame.hoverable(node, row),
      gradientId: 'saedsim-glow',
    });
    if (state.showKikuchi && state.kikuchi) {
      drawKikuchiBands(root, state.kikuchi.data, { width, height });
    }

    const scale = data.calibration.scale_px_per_inv_angstrom;
    const centre = pattern.centre_px;
    frame.configure({
      toData: (x, y) => ({ x, y }),
      formatCursor: (point) => {
        const radius = Math.hypot(point.x - centre[0], point.y - centre[1]);
        const g = scale > 0 ? radius / scale : 0;
        return (
          `${formatNumber(radius, 1)} px from 000 · ` +
          `|g| ${formatNumber(g, 4)} Å⁻¹ · d ${g > 0 ? formatNumber(1 / g, 4) : '∞'} Å`
        );
      },
    });
    frame.setContent(outer);
    frame.setReadout(orientationCard(data));
    frame.setStatus(statusLine(data));
  }

  /**
   * What is on the beam, and how the crystal is turned — beside the picture.
   *
   * The zone axis, the Euler angles of the same orientation, the roll about the
   * beam, and the deviation when the orientation asked for was not a zone axis.
   * The last of those is the one that must never be silent: a pattern drawn for
   * an orientation eight degrees off [001] is a picture of [001], and a reader
   * who is not told that will take it for the orientation they typed.
   */
  function orientationCard(data) {
    const orientation = data.orientation;
    const [phi1, phi, phi2] = orientation.euler_bunge_deg;
    const row = (label, value, title = '') =>
      el('tr', { title }, [
        el('th', { text: label }),
        el('td.measure__num', { text: value }),
      ]);
    return el('div.measure', {}, [
      el('div.measure__title', { text: 'Orientation on the beam' }),
      el('table.measure__table', {}, [
        el('tbody', {}, [
          row(
            data.zone_axis_notation_label ?? 'Zone axis',
            data.zone_axis_label,
            'The crystal direction along the beam',
          ),
          data.zone_axis_conversion_note
            ? row(
                'Entered [uvw]',
                data.zone_axis_three_index_label,
                'The same hexagonal direction in the three-index basis used by the input boxes',
              )
            : null,
          row(
            'Bunge (φ₁, Φ, φ₂)',
            `${formatNumber(phi1, 2)}, ${formatNumber(phi, 2)}, ${formatNumber(phi2, 2)}°`,
            'The same orientation, as a measured one would arrive',
          ),
          row(
            'Roll about the beam',
            `${formatNumber(orientation.in_plane_rotation_deg, 2)}°`,
            'One pattern cannot fix this; it is a statement about the plate',
          ),
          orientation.deviation_deg === null
            ? null
            : row(
                'Off that axis by',
                `${formatNumber(orientation.deviation_deg, 2)}°`,
                'How far the orientation given is from the zone actually drawn',
              ),
        ]),
      ]),
      data.zone_axis_conversion_note
        ? el('p.measure__note', { text: data.zone_axis_conversion_note })
        : null,
      orientation.deviation_deg > 0.05
        ? el('p.measure__note', {
            text:
              'The pattern drawn is the exact zone, not the orientation asked for — a spot ' +
              'pattern exists only on a zone axis.',
          })
        : null,
    ]);
  }

  function statusLine(data) {
    const pattern = data.pattern;
    const bands = state.showKikuchi ? state.kikuchi?.data?.bands?.length ?? 0 : 0;
    return (
      `${pattern.spots.length} reflections · ${data.phase_name} down ${data.zone_axis_label} · ` +
      `${formatNumber(data.calibration.scale_px_per_inv_angstrom, 2)} px per Å⁻¹` +
      (state.showKikuchi ? ` · ${bands} Kikuchi band(s)` : '') +
      ' · hover a spot for its index, d and |g|'
    );
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);
  else run();

  return { help: () => operation };
}
