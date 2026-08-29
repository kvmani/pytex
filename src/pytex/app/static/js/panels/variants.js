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
import { claim } from '../core/handoff.js';
import { identity, multiply, rotationX, rotationY } from '../core/rotation3.js';
import {
  compositeScene,
  contactSheetScene,
  measuredCompositeScene,
  screenAlignment,
} from '../core/compositescene.js';
import { defaultAppearance, renderScene } from './crystal.js';

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

/**
 * Phase colours for the composite views.
 *
 * In these views colour carries the **phase**, not the element. It has to: the
 * two phases of an orientation relationship are usually the same element —
 * austenite and ferrite are both iron — so colouring by species would paint
 * both crystals one colour and leave a single blob in which no parallelism can
 * be read. The legend states the change rather than leaving it to be inferred.
 */
const PHASE_COLORS = Object.freeze({
  parent: '#5b7fa6',
  child: '#d97706',
  // The idealized child, when it is shown beside the measured one. Muted,
  // because it is the comparison and not the measurement.
  ideal: '#9ca3af',
});

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
const VIEWS = [
  'variants.pole_figure',
  'variants.intervariant_misorientations',
  'variants.composite_scene',
  'variants.contact_sheet',
  'variants.or_from_grains',
  'variants.measured_composite',
];

/** Which kind of drawing each view needs. `table` has no picture, by design. */
const VIEW_MODES = {
  'variants.pole_figure': 'pole',
  'variants.intervariant_misorientations': 'table',
  'variants.composite_scene': 'composite',
  'variants.contact_sheet': 'sheet',
  'variants.or_from_grains': 'table',
  'variants.measured_composite': 'measured',
};

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
    appearance: defaultAppearance(),
    // Which result the camera was framed for. Reframing on every redraw would
    // make turning the crystal also rescale it.
    framedFor: null,
    // A pair of grains picked on the EBSD map, if the user arrived that way:
    // the seeded control values and the card that says where they came from.
    // Kept rather than applied once, because the two measured-pair views take
    // the same six angles and switching between them rebuilds the form.
    pair: null,
  };

  // One camera drives every crystal on screen -- both crystals of a composite,
  // and all twenty-four panels of a contact sheet. That is the point of having
  // the service place the child in Python: the lock is free and cannot drift,
  // because there is only ever one rotation.
  const camera = {
    rotation: multiply(rotationX(-1.2), rotationY(0.6)),
    zoom: 1,
    scale: 1,
    centre: [0, 0, 0],
    pan: { x: 0, y: 0 },
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

  /** The frame's heading, which has to follow the view rather than the panel. */
  const FRAME_TITLES = {
    pole: 'Pole figure',
    composite: 'Composite scene',
    sheet: 'Variant contact sheet',
    measured: 'The measured pair, in the specimen frame',
  };

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
        text: 'The pole figure shows where the variants point, the spectrum how they differ from each other, and the two composite views what a single variant and the whole family look like as crystals.',
      }),
    ]),
    formHost,
    runButton,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'Run the first two in order: 24 variants in 4 packets, then the same figure with half as many. The composite examples show the same variants as crystals.',
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

  function viewMode() {
    return VIEW_MODES[state.operation.id] ?? 'table';
  }

  function renderControls(initial = {}) {
    // A measured pair survives a change of view, because the two measured-pair
    // views are two questions about *one* pair: switching from the relationship
    // to the picture must not silently answer the second question about the
    // panel's default grains. Anything the caller supplies still wins, so
    // loading an example over a seeded pair replaces it, as it reads.
    const seed = state.pair && takesMeasuredPair() ? state.pair.values : {};
    state.form = buildForm(state.operation, { initial: { ...seed, ...initial } });
    formHost.replaceChildren(state.form.element);
  }

  /** Whether the open view is one that takes two measured orientations. */
  function takesMeasuredPair() {
    return (
      state.operation.id === 'variants.or_from_grains' ||
      state.operation.id === 'variants.measured_composite'
    );
  }

  function loadExample(example) {
    state.operation = operations.find((entry) => entry.id === example.operation);
    operationSelect.value = state.operation.id;
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  /**
   * Which run is the current one.
   *
   * The panel's views answer different questions and return differently shaped
   * results, and switching view starts a new request while the previous one is
   * still in flight. Without a token the slower response wins: the pole figure
   * that was already running lands *after* the composite scene, overwrites it,
   * and the drawing code for the new view is handed the old view's data. That
   * threw where the shapes disagree and, worse, silently drew stale numbers
   * where they happened to agree.
   */
  let runToken = 0;

  async function run() {
    const token = (runToken += 1);
    runButton.disabled = true;
    runButton.textContent = 'Computing…';
    state.form.clearErrors();
    try {
      const result = await call(state.operation.id, state.form.values());
      if (token !== runToken) return;
      state.result = result;
      state.hiddenPackets = new Set();
      frameCamera(result);
      draw();
      renderResult(details, result, { teaches: state.teaches });
      // The verdict and the statement go *above* the catalogue table, because
      // they are the answer and the table is the evidence for it.
      const verdict = verdictCard(result);
      if (verdict) details.prepend(verdict);
      // Above the verdict, because it is the question the verdict answers: the
      // measurement first, then what was computed from it.
      const received = receivedCard();
      if (received) details.prepend(received);
    } catch (error) {
      if (token !== runToken) return;
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      // Only the latest run owns the button; a superseded one must not
      // re-enable it while its successor is still computing.
      if (token === runToken) {
        runButton.disabled = false;
        runButton.textContent = 'Show variants';
      }
    }
  }

  /**
   * Point the camera at a new result, once.
   *
   * The centre and radius come from the service's world extent, which for a
   * contact sheet is taken over *every* variant's placement. Framing each panel
   * to its own contents would rescale them against each other, and a reader
   * comparing orientations across the grid would read that as a difference in
   * the crystallography.
   */
  function frameCamera(result) {
    const world = result?.data?.world;
    if (!world || state.framedFor === result) return;
    state.framedFor = result;
    camera.centre = world.centre;
    camera.scale = (VIEW * 0.8) / (world.radius || 1);
    camera.zoom = 1;
    camera.pan = { x: 0, y: 0 };
    camera.rotation = multiply(rotationX(-1.2), rotationY(0.6));
    // Deliberately *not* seeded with species colours: the renderer prefers
    // `appearance.speciesColors[species]` over the atom's own colour, so a
    // seeded entry for "Fe" would repaint both phases one colour and undo the
    // tint that tells them apart.
    state.appearance.speciesColors = {};
  }

  function draw() {
    // The spectrum and the measured-pair answer are tables. Hiding the frame
    // rather than drawing an empty one keeps the stage from carrying a picture
    // of nothing, which reads as a failed render rather than as a view that has
    // no picture.
    const mode = viewMode();
    frame.element.hidden = mode === 'table';
    legend.hidden = mode === 'table';
    frame.setTitle(FRAME_TITLES[mode] ?? 'Pole figure');
    publishButton.disabled = mode !== 'pole';
    formatSelect.disabled = mode !== 'pole';
    if (mode === 'table') return;
    if (mode === 'composite') return drawComposite();
    if (mode === 'measured') return drawMeasuredComposite();
    if (mode === 'sheet') return drawContactSheet();
    drawPoleFigure();
  }

  function drawPoleFigure() {
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

  /** One variant, both crystals, one depth-sorted drawing. */
  function drawComposite() {
    const data = state.result.data;
    frame.configure({ formatCursor: () => 'drag to turn · scroll to zoom' });
    frame.setContent(
      renderScene(
        compositeScene(data, { parentColor: PHASE_COLORS.parent, childColor: PHASE_COLORS.child }),
        camera,
        frame,
        compositeAppearance(state.appearance),
      ),
    );
    buildPhaseLegend(data);
    const notes = screenAlignment(camera.rotation, data.primitives);
    frame.setStatus(
      `Variant ${data.variant.index} of ${data.variant.count} · ` +
        `${data.parent.label} and ${data.child.label} in the parent frame · ` +
        alignmentStatus(notes),
    );
  }

  /**
   * Two measured grains, where the measurement puts them.
   *
   * The camera needs no locking here: the relative placement came from the data,
   * so turning the view turns one rigid arrangement. The status line carries the
   * clause deviations, because they are what the visible gap between each pair
   * of overlays measures, and a reader should not have to estimate an angle off
   * the screen.
   */
  function drawMeasuredComposite() {
    const data = state.result.data;
    frame.configure({ formatCursor: () => 'drag to turn · scroll to zoom' });
    frame.setContent(
      renderScene(
        measuredCompositeScene(data, {
          parentColor: PHASE_COLORS.parent,
          childColor: PHASE_COLORS.child,
          idealColor: PHASE_COLORS.ideal,
        }),
        camera,
        frame,
        compositeAppearance(state.appearance),
      ),
    );
    buildMeasuredLegend(data);
    const clauses = (data.parallelisms ?? [])
      .map((row) => `${row.parent} ∥ ${row.child} ±${formatNumber(row.deviation_deg, 2)}°`)
      .join(' · ');
    const ideal = data.idealized
      ? ` · idealized child ${formatNumber(data.idealized.turn_deg, 2)}° away`
      : '';
    frame.setStatus(
      `Specimen frame (RD, TD, ND) · ${data.naming.best_label ?? 'no named relationship'} at ` +
        `${formatNumber(data.naming.best_deviation_deg, 2)}° · ${clauses}${ideal}`,
    );
  }

  /**
   * The measured composite's legend.
   *
   * It has one job the catalogue legend does not: saying that each overlay is
   * drawn on both sides, so a reader knows the doubled patch is the point
   * rather than a rendering fault.
   */
  function buildMeasuredLegend(data) {
    if (state.legendFor === data) return;
    state.legendFor = data;
    legend.hidden = false;
    const items = [
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: `background:${PHASE_COLORS.parent}` }),
        el('span', { text: `${data.parent.label} (measured)` }),
      ]),
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: `background:${PHASE_COLORS.child}` }),
        el('span', { text: `${data.child.label} (measured)` }),
      ]),
    ];
    if (data.idealized) {
      items.push(
        el('span.legend__item', {}, [
          el('span.legend__swatch', { style: `background:${PHASE_COLORS.ideal}` }),
          el('span', { text: 'Child as the integer statement would place it' }),
        ]),
      );
    }
    legend.replaceChildren(
      ...items,
      el('span.legend__guide', {
        text:
          'Each parallelism is drawn on both sides, in the two phase colours: the gap between ' +
          'them is the clause deviation, not a rendering fault.',
      }),
    );
  }

  /**
   * Every variant at one camera.
   *
   * The two structures were fetched once; each panel applies its own placement
   * matrix here. Fetching a scene per variant would be the same picture for
   * twenty-four times the traffic, and is exactly what the payload shape exists
   * to avoid.
   */
  function drawContactSheet() {
    const data = state.result.data;
    const appearance = compositeAppearance(state.appearance, { compact: true });
    const sheet = el(
      'div.sheet',
      {},
      data.variants.map((entry) => {
        const plane = entry.parallelisms.find((row) => row.kind === 'plane');
        return el(
          'button.sheet__cell',
          {
            type: 'button',
            dataset: { variant: String(entry.index), packet: String(entry.packet) },
            title: `Open variant ${entry.index} in the composite view`,
            onclick: () => {
              if (dragMoved) return;
              openVariant(entry.index);
            },
          },
          [
            renderScene(
              contactSheetScene(data, entry, {
                parentColor: PHASE_COLORS.parent,
                childColor: PHASE_COLORS.child,
              }),
              camera,
              SILENT_FRAME,
              appearance,
            ),
            el('span.sheet__caption', {}, [
              el('strong', { text: `V${entry.index}` }),
              el('span.sheet__packet', {}, [
                el('span.sheet__swatch', {
                  style: `background:${packetColor(entry.packet)}`,
                }),
                el('span', { text: `packet ${entry.packet}` }),
              ]),
            ]),
            plane
              ? el('span.sheet__pair', { text: `${plane.parent} \u2225 ${plane.child}` })
              : null,
          ],
        );
      }),
    );
    buildPhaseLegend(data);
    frame.setContent(sheet);
    frame.setStatus(
      `${data.variant_count} variants in ${data.packet_count} packets at one camera · ` +
        'drag to turn all of them together · click a panel to open it',
    );
  }

  /**
   * The legend of a composite view: which colour is which phase.
   *
   * Not the packet legend, and not a filter — the two crystals are the figure,
   * so there is nothing to switch off. It exists because colour means something
   * different here than in the pole figure, and an unexplained change of
   * meaning is worse than either meaning alone.
   */
  function buildPhaseLegend(data) {
    if (state.legendFor === data) return;
    state.legendFor = data;
    legend.hidden = false;
    legend.replaceChildren(
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: `background:${PHASE_COLORS.parent}` }),
        el('span', { text: `${data.parent.label} (parent)` }),
      ]),
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: `background:${PHASE_COLORS.child}` }),
        el('span', { text: `${data.child.label} (child)` }),
      ]),
      el('span.legend__guide', {
        text: 'Colour is the phase, not the element: both phases are usually the same element.',
      }),
    );
  }

  /** Switch to the composite view already showing the variant that was clicked. */
  function openVariant(index) {
    const target = operations.find((entry) => entry.id === 'variants.composite_scene');
    if (!target) return;
    const request = { ...state.form.values(), variant: index };
    delete request.packet_plane;
    state.operation = target;
    operationSelect.value = target.id;
    renderControls(request);
    run();
  }

  /**
   * The measured-pair answer, stated before the evidence for it.
   *
   * Four different angles reach this panel, and the one thing the card must not
   * do is let them blur together. Each is shown with the word for what it
   * measures, taken from the service's own `angle_meanings` so the panel and the
   * operation's help cannot drift apart.
   */
  function verdictCard(result) {
    if (state.operation.id !== 'variants.or_from_grains') return null;
    const data = result.data;
    const meanings = data.angle_meanings ?? {};
    const naming = data.naming;
    const statement = data.statement;
    const rows = [
      el('div.verdict__row', {}, [
        el('span.verdict__label', { text: 'Relationship' }),
        el('strong', {
          text: naming.is_conclusive
            ? `${naming.best_label} (conclusive)`
            : `${naming.best_label ?? 'none'} — not conclusive`,
        }),
      ]),
      el('div.verdict__row', {}, [
        el('span.verdict__label', {
          text: 'Catalogue distance',
          title: meanings.catalog,
        }),
        el('span', { text: `${formatNumber(naming.best_deviation_deg, 3)}°` }),
      ]),
      el('div.verdict__row', {}, [
        el('span.verdict__label', {
          text: 'Lead over runner-up',
          title: 'The margin the verdict rests on: it must exceed both the scatter and the misfit.',
        }),
        el('span', { text: `${formatNumber(naming.margin_deg, 3)}°` }),
      ]),
    ];
    if (statement) {
      rows.push(
        el('div.verdict__row', {}, [
          el('span.verdict__label', { text: 'Integer statement' }),
          el('strong', { text: statement.text }),
        ]),
        el('div.verdict__row', {}, [
          el('span.verdict__label', {
            text: 'Cost of writing it',
            title: meanings.rationalization,
          }),
          el('span', { text: `${formatNumber(statement.rationalization_cost_deg, 3)}°` }),
        ]),
      );
    }
    return el('section.verdict', {}, [
      el('h3', { text: 'The answer' }),
      ...rows,
      el('p.verdict__note', {
        text: statement
          ? `The statement is an idealization: it names a nearby exact relationship, not the ` +
            `measurement. ${meanings.rationalization ?? ''}`
          : data.statement_note ??
            'No integer statement was found within the index bound; the rotation stands alone.',
      }),
    ]);
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

  /*
   * Turning the crystals.
   *
   * Watched on the frame rather than on the drawing, because a redraw replaces
   * the drawing: handlers owned by the SVG see the first movement of a drag and
   * nothing after it, and the symptom is a crystal that nudges once however far
   * you pull. This is the same reason, and the same fix, as in the crystal
   * viewer.
   */
  let dragging = null;
  // A contact-sheet panel is a button *and* a drag surface. Capturing the
  // pointer or preventing the default on pointerdown kills the click that opens
  // the variant, so the drag only takes over once the pointer has actually
  // moved -- and the click that follows a real drag is then ignored, which is
  // what stops "turn the sheet" from also meaning "open whichever panel you
  // happened to start on".
  const DRAG_THRESHOLD_PX = 3;
  let dragMoved = false;

  frame.element.addEventListener('pointerdown', (event) => {
    // Cleared for *every* press, including one that starts outside a drawing:
    // a press on a panel's caption must not inherit the "this was a drag" flag
    // from the last turn of the crystals and be swallowed by it.
    dragMoved = false;
    if (!state.result || !is3D()) return;
    if (!event.target.closest('svg')) return;
    dragging = { x: event.clientX, y: event.clientY, id: event.pointerId };
  });
  frame.element.addEventListener('pointermove', (event) => {
    if (!dragging) return;
    const moveX = event.clientX - dragging.x;
    const moveY = event.clientY - dragging.y;
    if (!dragMoved) {
      if (Math.hypot(moveX, moveY) < DRAG_THRESHOLD_PX) return;
      dragMoved = true;
      frame.element.setPointerCapture(event.pointerId);
    }
    event.preventDefault();
    dragging = { ...dragging, x: event.clientX, y: event.clientY };
    // Pre-multiplying rotates about the *screen* axes, so "drag right turns
    // right" stays true however the crystals are already oriented.
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
  frame.element.addEventListener(
    'wheel',
    (event) => {
      if (!state.result || !is3D()) return;
      event.preventDefault();
      camera.zoom = Math.min(Math.max(camera.zoom * (event.deltaY < 0 ? 1.12 : 0.89), 0.2), 12);
      draw();
    },
    { passive: false },
  );

  function is3D() {
    const mode = viewMode();
    return mode === 'composite' || mode === 'sheet' || mode === 'measured';
  }

  /**
   * Take a pair of grains picked on the EBSD map, if one was offered.
   *
   * Claimed while mounting, once: the offer is a gesture the user just made —
   * two grains and a button — and not a setting. Leaving it readable would make
   * every later visit to this workspace re-seed itself from a pick made minutes
   * ago, which is the same silent substitution as a map that quietly analyses
   * the practice dataset instead of the open scan.
   *
   * **The phases are the reason this is not simply "fill the form and run".**
   * An orientation relationship is defined between two *distinct* phases, and a
   * scan names phases without describing them. Where both names resolve to
   * built-in phases and the two differ, the pair is complete and the answer is
   * computed immediately. Where they do not — a single-phase scan, or a name
   * this application does not carry — the angles are seeded and the run is
   * *not* started: the panel's default phases are austenite and ferrite, and
   * computing under them would report a relationship between two phases the
   * measurement never claimed.
   */
  function seedFromPickedPair() {
    const offered = claim('measured-pair');
    if (!offered?.grains || offered.grains.length !== 2) return false;
    const [parent, child] = offered.grains;
    const values = {
      euler_convention: offered.euler_convention ?? 'bunge',
      parent_angle1: parent.phi1_deg,
      parent_angle2: parent.Phi_deg,
      parent_angle3: parent.phi2_deg,
      child_angle1: child.phi1_deg,
      child_angle2: child.Phi_deg,
      child_angle3: child.phi2_deg,
    };
    const distinct =
      parent.phase_builtin && child.phase_builtin && parent.phase_builtin !== child.phase_builtin;
    if (distinct) {
      values.phase = { builtin: parent.phase_builtin };
      values.child_phase = { builtin: child.phase_builtin };
    }
    state.pair = { ...offered, values, phasesResolved: Boolean(distinct) };
    state.operation = operations.find((entry) => entry.id === 'variants.or_from_grains');
    operationSelect.value = state.operation.id;
    state.teaches = null;
    renderControls();
    if (distinct) {
      run();
    } else {
      // Nothing computed, so nothing to put the card above: it goes on its own,
      // where a result would have been. The frame still has to be told this
      // view has no picture, or the stage carries an empty pole figure — which
      // reads as a failed render rather than as a question waiting on a phase.
      draw();
      details.replaceChildren(receivedCard());
    }
    return true;
  }

  /**
   * Where the two orientations came from, and what believing them costs.
   *
   * Shown for as long as the pair is the panel's input, not only on the run
   * that arrived: a user who raises the index bound and runs again is still
   * looking at an answer about two measured grains, and the spread that
   * qualifies it must not scroll away with the first result.
   *
   * The grain orientation spread is on the card because a single pair's
   * residual is zero *by construction* — two mean orientations always fit one
   * rotation exactly — so the spread is the only measure of what that zero
   * conceals. It is the same lesson as the Kearns triad that summed to one
   * whatever the data were.
   */
  function receivedCard() {
    if (!state.pair || !takesMeasuredPair()) return null;
    const [parent, child] = state.pair.grains;
    const spread = Math.max(
      parent.grain_orientation_spread_deg ?? 0,
      child.grain_orientation_spread_deg ?? 0,
    );
    const describe = (grain, role) =>
      el('div.received__grain', {}, [
        el('strong', { text: `${role}: grain ${grain.grain_id}` }),
        el('span', {
          text:
            ` — ${formatNumber(grain.phi1_deg, 2)}, ${formatNumber(grain.Phi_deg, 2)}, ` +
            `${formatNumber(grain.phi2_deg, 2)}° Bunge · ` +
            `${grain.phase_name ?? 'phase not named'}` +
            ` · GOS ${formatNumber(grain.grain_orientation_spread_deg, 3)}°`,
        }),
      ]);
    return el('section.received', {}, [
      el('h3', { text: 'Two grains picked off the map' }),
      describe(parent, 'Parent'),
      describe(child, 'Child'),
      el('p.received__note', {
        text:
          `From ${state.pair.source?.map ?? 'the orientation map'}` +
          (state.pair.source?.grain_threshold_deg
            ? `, segmented at ${state.pair.source.grain_threshold_deg}°`
            : '') +
          '. Each orientation is the grain’s symmetry-aware mean over its own points.',
      }),
      el('p.received__note', {
        text:
          'A mean has no scatter, so this pair fits one rotation exactly and the residual is ' +
          `zero whatever the grains were. The spread is what that zero conceals: up to ` +
          `${formatNumber(spread, 3)}° within a grain. Read every angle below against it.`,
      }),
      state.pair.phasesResolved
        ? null
        : el('p.received__note.received__note--warning', {
            text:
              'The phases were not carried across: the scan names both grains ' +
              `${parent.phase_name ?? 'the same phase'}` +
              (child.phase_name && child.phase_name !== parent.phase_name
                ? ` and ${child.phase_name}`
                : '') +
              ', and a relationship is defined between two distinct phases. Choose the parent ' +
              'and child phases above — the ones shown are this panel’s defaults, not the ' +
              'scan’s — and press Show variants.',
          }),
    ]);
  }

  renderControls();
  // The legend is a control, so it rides inside the frame rather than under it:
  // toggling a source and seeing the drawing change must not need a scroll.
  frame.setControls(legend);
  context.stage.append(frame.element, details);
  // A pair handed over from the map is what the user asked for by pressing the
  // button; the opening example is only what the panel does when nobody asked.
  if (!seedFromPickedPair() && examples.length) loadExample(examples[0]);

  return { help: () => state.operation };
}

/**
 * Draw one variant's two crystals, or the whole family at one camera.
 *
 * Both come from the crystal viewer's renderer, not from a second one: the
 * composition happens in `core/compositescene.js`, which concatenates payloads
 * into the shape that renderer already takes. That matters for more than reuse
 * — a single renderer means a single depth sort, so parent and child atoms
 * occlude each other correctly instead of one crystal being drawn wholly in
 * front of the other.
 */
function compositeAppearance(base, { compact = false } = {}) {
  return {
    ...base,
    showBonds: compact ? false : base.showBonds,
    showLabels: compact ? false : base.showLabels,
    showGizmo: compact ? false : base.showGizmo,
    showCells: base.showCells,
    surfaceFinish: compact ? 'flat' : base.surfaceFinish,
    atomScale: compact ? base.atomScale * 0.9 : base.atomScale,
    annotationScale: base.annotationScale,
  };
}

/** A frame that accepts hover registrations and drops them, for the small panels. */
const SILENT_FRAME = { hoverable() {} };

/**
 * The status line of the composite view: where the parallelism is, right now.
 *
 * A picture cannot answer "am I looking down it yet?" — a plane a couple of
 * degrees off edge-on looks edge-on. The numbers can, so they are on the status
 * line and they update as the crystal turns.
 */
function alignmentStatus(notes) {
  if (!notes.length) return 'no parallel plane or direction is being drawn';
  return notes
    .map(
      (note) =>
        `${note.label} ${formatNumber(note.angleDeg, 1)}° from ${note.aligned}`,
    )
    .join(' · ');
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
