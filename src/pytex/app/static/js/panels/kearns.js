/**
 * The Kearns panel: one texture number, by three routes, read against 1/3.
 *
 * Like the calculator, almost none of this file is about crystallography. The
 * route list, every control, and every help string come from the manifest;
 * Python computes `f`, the orientation tensor, and the tilt profile, and this
 * module lays them out. Adding a `kearns.*` operation in Python makes it appear
 * here with no change to this file.
 *
 * The two panel-specific views exist because the table alone does not answer
 * the questions this number is asked:
 *
 * - **The triad card.** A Kearns parameter is meaningless without its baseline,
 *   so every value is drawn as a bar against 1/3 — the untextured value — with
 *   the bar's zero at 1/3 rather than at 0. A texture-free specimen then reads
 *   as three stubs of nothing, which is what it is. Beneath it sits the closure
 *   check: `f_RD + f_TD + f_ND` is identically 1 for *any* texture, so the sum
 *   tests the measurement and never the material. It is the first thing a
 *   reviewer looks for and the last thing most reports print.
 *
 * - **The tilt profile.** For the two diffraction routes, `f` is the
 *   cos²-weighted mean of the azimuthally averaged pole density `I(phi)`. The
 *   plot draws that density *and* the volume fraction `I sin(phi)` it is
 *   weighted by, on one pair of axes, because the single most misread thing in
 *   the method is that a towering density at low tilt contributes almost
 *   nothing: the ring of crystals at zero tilt has no circumference. Drawing
 *   the two curves together makes that visible instead of arithmetical.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { explainer } from '../core/explainer.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { download, renderResult } from '../core/result.js';
import { call } from '../core/api.js';

export const panel = {
  id: 'kearns',
  title: 'Kearns parameter',
  tagline: 'The basal-texture number f, by three independent routes.',
};

/** The untextured value of f in every direction. The baseline everything is read against. */
const ISOTROPIC = 1 / 3;

/* The profile plot's viewBox.
 *
 * Deliberately landscape, and by roughly the ratio the plot frame's canvas
 * actually has (a wide, short strip beneath the card header). The frame fits a
 * drawing by its limiting dimension, so a squarer viewBox is letterboxed into
 * the middle of the strip and draws the curve at half the size it could be —
 * measured at 40% of the available width before this was widened. */
const PLOT_W = 380;
const PLOT_H = 140;
const MARGIN = { left: 48, right: 18, top: 14, bottom: 32 };

/** How far the triad sum may sit from 1 before the closure check is called a failure. */
const CLOSURE_TOLERANCE = 5e-4;

/**
 * Mount the panel.
 *
 * @param {object} context - `{stage, rail, manifest, showError, setBusy, openHelp}`.
 */
export function mount(context) {
  const operations = context.manifest.operations.filter((entry) => entry.panel === panel.id);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = {
    operation: operations[0],
    form: null,
    teaches: null,
    plotNode: null,
    // Opened pole-figure files, kept on the panel rather than in the form:
    // switching route must not silently drop what the user has opened.
    files: [],
  };

  const chooser = el(
    'select',
    {
      'aria-label': 'Route',
      onchange: () => {
        state.operation = operations.find((entry) => entry.id === chooser.value);
        state.teaches = null;
        renderControls();
      },
    },
    operations.map((operation) =>
      el('option', { value: operation.id, text: operation.title, title: operation.summary }),
    ),
  );

  const routeHelp = el('p.field__help', { text: state.operation?.summary ?? '' });
  const fileHost = el('div.group__body');
  const fileGroup = el('details.group', { open: true, id: 'kearns-files' }, [
    el('summary', { text: 'Open pole figures' }),
    fileHost,
  ]);
  const formHost = el('div');
  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Compute f',
    onclick: () => run(),
  });

  const exampleList = el(
    'div.examples',
    {},
    examples.map((example) =>
      el(
        'button.example',
        { type: 'button', onclick: () => loadExample(example) },
        [el('strong', { text: example.title }), el('span', { text: example.summary })],
      ),
    ),
  );

  context.rail.append(
    el('div.field', {}, [
      el('label.field__label', { text: 'Route' }),
      chooser,
      routeHelp,
    ]),
    fileGroup,
    formHost,
    runButton,
    el('details.group', { open: examples.length > 0 }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        explainer(
          'The three routes on one synthetic specimen whose true f is known. Run them in '
            + 'order and watch the same number come back three different ways.',
          { label: 'What these examples show' },
        ),
        exampleList,
      ]),
    ]),
  );

  /** Whether the selected route reads its data from opened files. */
  function needsFiles() {
    return state.operation.parameters.some((parameter) => parameter.name === 'files');
  }

  function renderControls(initial = {}) {
    chooser.value = state.operation.id;
    routeHelp.textContent = state.operation.summary ?? '';
    state.form = buildForm(state.operation, { initial });
    formHost.replaceChildren(state.form.element);
    renderFileControls();
  }

  /*
   * The file picker, shown only for the routes that read a measurement.
   *
   * Several files at once, because the ODF route needs a *set* — one pole
   * figure cannot constrain the inversion at all. The pole-figure route refuses
   * more than one and says why, rather than quietly integrating the first.
   */
  function renderFileControls() {
    const wanted = needsFiles();
    fileGroup.hidden = !wanted;
    if (!wanted) {
      fileHost.replaceChildren();
      return;
    }
    const input = el('input', {
      type: 'file',
      accept: '.xrdml',
      multiple: true,
      'aria-label': 'Open XRDML pole-figure files',
      onchange: (event) => openFiles([...(event.target.files ?? [])]),
    });
    fileHost.replaceChildren(
      explainer(
        'Panalytical XRDML pole-figure files. Give the plane of each below, in the order '
          + 'they are listed here: the file records the diffraction angle, not the reflection.',
        { label: 'Which files, in which order' },
      ),
      input,
      el('p.field__help', {
        text: state.files.length
          ? `${state.files.length} file(s) open: ${state.files.map((file) => file.name).join(', ')}`
          : 'No file open yet. Choose one or more .xrdml pole-figure files.',
      }),
      state.files.length
        ? el('button.button', {
            type: 'button',
            text: 'Close them',
            onclick: () => {
              state.files = [];
              renderFileControls();
            },
          })
        : null,
    );
  }

  async function openFiles(files) {
    if (!files.length) return;
    try {
      state.files = await Promise.all(
        files.map(async (file) => ({ name: file.name, text: await file.text() })),
      );
      renderFileControls();
      await run();
    } catch (error) {
      context.showError(error);
    }
  }

  function loadExample(example) {
    state.operation = operations.find((entry) => entry.id === example.operation);
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    const operation = state.operation;
    runButton.disabled = true;
    runButton.textContent = 'Computing…';
    state.form.clearErrors();
    try {
      const request = state.form.values();
      if (needsFiles()) {
        if (!state.files.length) {
          context.showError(
            new Error('Open at least one XRDML pole-figure file before computing f.'),
          );
          return;
        }
        request.files = { items: state.files };
      }
      const result = await call(operation.id, request);
      renderResult(context.stage, result, {
        extra: views(result, state),
        teaches: state.teaches,
      });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Compute f';
    }
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);

  return {
    /** The help drawer shows the route currently selected, not the panel. */
    help: () => state.operation,
  };
}

/** The panel-specific cards, in reading order. */
function views(result, state) {
  const nodes = [];
  const triad = triadCard(result);
  if (triad) nodes.push(triad);
  const tensor = tensorCard(result);
  if (tensor) nodes.push(tensor);
  const profile = profileCard(result, state);
  if (profile) nodes.push(profile);
  return nodes;
}

/**
 * The headline: every f as a bar against the untextured 1/3, plus the closure check.
 */
function triadCard(result) {
  const directions = result.data?.directions;
  if (!Array.isArray(directions) || !directions.length) return null;

  const values = directions.map((entry) => Number(entry.f));
  const largest = Math.max(...values, ISOTROPIC * 2);

  const meters = directions.map((entry) => {
    const value = Number(entry.f);
    const ratio = value / ISOTROPIC;
    // The bar runs from 0 to the larger of the data and 2/3, so the 1/3 tick
    // always sits inside the track and a weak texture is not magnified into a
    // strong-looking one by an auto-scaled axis.
    const width = Math.max(0, Math.min(1, value / largest)) * 100;
    const baseline = (ISOTROPIC / largest) * 100;
    return el('div.kearns__row', {}, [
      el('span.kearns__label', { text: `f_${entry.label}` }),
      el('strong.kearns__value', { text: formatNumber(value, 4) }),
      el('div.kearns__track', { title: `${formatNumber(ratio, 2)} times the random value` }, [
        el('div.kearns__bar', { style: `width:${width.toFixed(2)}%` }),
        el('div.kearns__baseline', { style: `left:${baseline.toFixed(2)}%` }),
      ]),
      el('span.kearns__ratio', { text: `${formatNumber(ratio, 2)}× random` }),
    ]);
  });

  return el('section.card', {}, [
    el('div.card__header', {}, [
      el('h2.card__title', { text: 'The Kearns parameter' }),
      el('p.card__subtitle', {
        text: 'The tick on each bar is 1/3, the value of an untextured aggregate.',
      }),
    ]),
    el('div.card__body', {}, [el('div.kearns', {}, meters), closureNote(result), coverageNote(result)]),
  ]);
}

/**
 * The closure check — and, more importantly, what it does *not* prove.
 *
 * The sum rule has content only where the three values were measured
 * independently. A route that builds one pole orientation tensor and reads the
 * triad off its diagonal closes automatically: the tensor averages `c cᵀ` over
 * unit vectors, so its trace is 1 whatever the data were. Reporting that as a
 * passed check would be false reassurance of the exact kind this panel exists
 * to prevent — a basal pole figure truncated at 60° closes to 1.0000 and can
 * still be wrong by half, because the pseudo-norm renormalises over the
 * measured cap alone and that cap is biased towards the pole.
 *
 * So the verdict is written three ways, and the by-construction case points at
 * the coverage diagnostic that genuinely tests it.
 */
function closureNote(result) {
  const sum = result.data?.triad_sum;
  if (sum === null || sum === undefined) {
    return el('p.field__help', {
      text:
        'This route reports one direction per section, so there is no triad to close. '
        + 'Measure all three principal sections to obtain the full set — and then the sum '
        + 'becomes a real check, because those three values are independent measurements.',
    });
  }
  if (result.data?.orientation_tensor) {
    return el('p.notes--info', {}, [
      el('strong', { text: 'Closes by construction. ' }),
      document.createTextNode(
        `The triad sums to ${formatNumber(sum, 4)}, but these three values are the diagonal `
        + 'of one pole orientation tensor whose trace is identically 1. The sum therefore '
        + 'tests the arithmetic, not the measurement, and a route can close perfectly while '
        + 'being badly wrong. Read the coverage below instead.',
      ),
    ]);
  }
  const departure = Math.abs(Number(sum) - 1);
  const closes = departure <= CLOSURE_TOLERANCE;
  return el(closes ? 'p.notes--ok' : 'p.notes--warn', {}, [
    el('strong', { text: closes ? 'Closure check passes. ' : 'Closure check fails. ' }),
    document.createTextNode(
      closes
        ? `The triad sums to ${formatNumber(sum, 4)}. These sections were measured `
          + 'independently, so the sum rule has content here and nothing is unaccounted for.'
        : `The triad sums to ${formatNumber(sum, 4)}, off by ${formatNumber(departure, 4)}. `
          + 'These sections were measured independently and the sum is identically 1 for every '
          + 'texture, so this measures a systematic error in the measurement — an unmeasured '
          + 'tilt range, a wrong random standard, an unbalanced background — and says nothing '
          + 'about the material.',
    ),
  ]);
}

/**
 * How much of the hemisphere was actually integrated.
 *
 * For the pole-figure route this is the diagnostic that matters, and it is the
 * one the closure check is mistaken for. A figure truncated at 60° covers half
 * the hemisphere, and the missing half is the high-tilt half that the `sin φ`
 * weighting makes count most.
 */
function coverageNote(result) {
  const covered = result.data?.diagnostics?.measured_solid_angle_fraction;
  const maxPolar = result.data?.diagnostics?.max_polar_deg;
  if (covered === undefined || maxPolar === undefined) return null;
  const complete = Number(maxPolar) >= 89;
  return el(complete ? 'p.notes--ok' : 'p.notes--warn', {}, [
    el('strong', {
      text: complete ? 'Coverage is complete. ' : 'The figure is incomplete. ',
    }),
    document.createTextNode(
      complete
        ? `Measured to ${formatNumber(maxPolar, 1)}° of tilt, covering the whole hemisphere, `
          + 'so no unmeasured region had to be assumed.'
        : `Measured only to ${formatNumber(maxPolar, 1)}° of tilt, covering `
          + `${formatNumber(Number(covered) * 100, 0)}% of the hemisphere. The result `
          + 'renormalises over the measured cap alone (the Kern-Bergmann pseudo-norm) and so '
          + 'assumes the unmeasured cap resembles it. Because the sin φ weighting makes high '
          + 'tilt count most, that assumption is where this route goes wrong — and the triad '
          + 'still sums to 1 while it does.',
    ),
  ]);
}

/** The orientation tensor and its principal values, when the route determines one. */
function tensorCard(result) {
  const tensor = result.data?.orientation_tensor;
  if (!Array.isArray(tensor) || tensor.length !== 3) return null;
  const principal = result.data?.principal_values;
  const labels = ['RD', 'TD', 'ND'];

  const header = el('tr', {}, [
    el('th', { text: '' }),
    ...labels.map((label) => el('th', { text: label })),
  ]);
  const rows = tensor.map((row, index) =>
    el('tr', {}, [
      el('th', { text: labels[index] }),
      ...row.map((value) =>
        el('td.numeric', {
          text: formatNumber(value, 4),
          title: `${formatNumber(value, 8)}`,
        }),
      ),
    ]),
  );

  return el('section.card', {}, [
    el('div.card__header', {}, [
      el('h2.card__title', { text: 'Pole orientation tensor' }),
      el('p.card__subtitle', {
        text:
          'A = <c cᵀ>, the second moment of the basal-pole directions. f along any direction d '
          + 'is dᵀ A d, so this one object answers for every direction at once — not only the '
          + 'three axes above. Its trace is 1, which is why the triad sums to 1.',
      }),
    ]),
    el('div.card__body', {}, [
      el('div.table-wrap', {}, [
        el('table.result', {}, [el('thead', {}, header), el('tbody', {}, rows)]),
      ]),
      Array.isArray(principal)
        ? el('p.field__help', {
            text:
              `Principal values ${principal.map((value) => formatNumber(value, 4)).join(', ')}. `
              + 'These are f along the texture’s own principal axes, so the largest is the '
              + 'greatest f obtainable in any direction and the smallest the least: together '
              + 'they bound every direction’s value.',
          })
        : null,
    ]),
  ]);
}

/**
 * The tilt profile, with the volume-fraction weighting drawn beside the density.
 */
function profileCard(result, state) {
  const profile = result.data?.profile;
  const polar = profile?.polar_deg;
  const intensity = profile?.intensity;
  if (!Array.isArray(polar) || !Array.isArray(intensity) || polar.length < 2) return null;

  const frame = plotFrame({
    title: 'Basal-pole tilt profile',
    units: '',
    toolbar: [
      el('button.button', {
        type: 'button',
        text: 'SVG',
        title: 'Save the tilt profile as SVG',
        onclick: () => {
          if (!state.plotNode) return;
          const markup = new XMLSerializer().serializeToString(state.plotNode);
          download('pytex-kearns-tilt-profile.svg', markup, 'image/svg+xml');
        },
      }),
    ],
  });

  const node = profilePlot(polar, intensity);
  state.plotNode = node;
  frame.setContent(node);
  frame.setStatus(
    'Blue is the measured pole density I(φ). Amber is the volume fraction I(φ)·sin φ that f '
    + 'actually averages — the two peak at different tilts, and the amber one is the one that '
    + 'counts.',
  );

  return el('section.card', {}, [
    el('div.card__header', {}, [
      el('h2.card__title', { text: 'What f was integrated from' }),
      el('p.card__subtitle', {
        text:
          'f is the cos²φ-weighted mean of the volume fraction. A gap in the middle of the '
          + 'tilt range means the integral is interpolating over something nothing measured.',
      }),
    ]),
    el('div.card__body', {}, [frame.element]),
  ]);
}

/**
 * Draw density and volume fraction against tilt, each on its own normalised scale.
 *
 * Two separate normalisations rather than one shared axis: the curves differ by
 * more than an order of magnitude at low tilt, and a shared axis would flatten
 * the volume fraction into the baseline — hiding the very comparison the plot
 * is drawn for. Each curve is therefore scaled to its own maximum, and the axis
 * is labelled as relative for that reason rather than carrying a false unit.
 */
function profilePlot(polar, intensity) {
  const angles = polar.map(Number);
  const density = intensity.map(Number);
  const fraction = angles.map((angle, index) => density[index] * Math.sin((angle * Math.PI) / 180));
  const maxDensity = Math.max(...density) || 1;
  const maxFraction = Math.max(...fraction) || 1;

  const plotWidth = PLOT_W - MARGIN.left - MARGIN.right;
  const plotHeight = PLOT_H - MARGIN.top - MARGIN.bottom;
  const x = (angle) => MARGIN.left + (angle / 90) * plotWidth;
  const y = (value) => MARGIN.top + (1 - value) * plotHeight;

  const path = (values, max) =>
    values
      .map((value, index) => `${index === 0 ? 'M' : 'L'} ${x(angles[index]).toFixed(2)} ${y(value / max).toFixed(2)}`)
      .join(' ');

  const ticks = [0, 15, 30, 45, 60, 75, 90];
  const gridlines = ticks.map((angle) =>
    svg('line', {
      x1: x(angle), y1: MARGIN.top, x2: x(angle), y2: MARGIN.top + plotHeight,
      stroke: 'currentColor', 'stroke-width': 0.3, opacity: 0.18,
    }),
  );
  const tickLabels = ticks.map((angle) =>
    svg('text', {
      x: x(angle), y: MARGIN.top + plotHeight + 11,
      'text-anchor': 'middle', 'font-size': 7, fill: 'currentColor', text: String(angle),
    }),
  );

  const marks = (values, max, color) =>
    values.map((value, index) =>
      svg('circle', {
        cx: x(angles[index]), cy: y(value / max), r: 1.6, fill: color,
      }, [svg('title', { text: `${angles[index]}°: ${formatNumber(value, 4)}` })]),
    );

  return svg(
    'svg',
    {
      viewBox: `0 0 ${PLOT_W} ${PLOT_H}`,
      xmlns: 'http://www.w3.org/2000/svg',
      role: 'img',
      'aria-label': 'Basal-pole density and volume fraction against tilt angle',
    },
    [
      ...gridlines,
      svg('line', {
        x1: MARGIN.left, y1: MARGIN.top + plotHeight, x2: MARGIN.left + plotWidth,
        y2: MARGIN.top + plotHeight, stroke: 'currentColor', 'stroke-width': 0.6,
      }),
      svg('line', {
        x1: MARGIN.left, y1: MARGIN.top, x2: MARGIN.left, y2: MARGIN.top + plotHeight,
        stroke: 'currentColor', 'stroke-width': 0.6,
      }),
      svg('path', {
        d: path(density, maxDensity), fill: 'none', stroke: '#31548c', 'stroke-width': 1.2,
      }),
      svg('path', {
        d: path(fraction, maxFraction), fill: 'none', stroke: '#d66a2b', 'stroke-width': 1.2,
      }),
      ...marks(density, maxDensity, '#31548c'),
      ...marks(fraction, maxFraction, '#d66a2b'),
      ...tickLabels,
      svg('text', {
        x: MARGIN.left + plotWidth / 2, y: PLOT_H - 6, 'text-anchor': 'middle',
        'font-size': 8, fill: 'currentColor', text: 'Tilt φ from the reference direction (°)',
      }),
      svg('text', {
        x: 10, y: MARGIN.top + plotHeight / 2, 'font-size': 8, fill: 'currentColor',
        'text-anchor': 'middle', transform: `rotate(-90 10 ${MARGIN.top + plotHeight / 2})`,
        text: 'Relative to each curve’s maximum',
      }),
      svg('text', {
        x: MARGIN.left + plotWidth, y: MARGIN.top + 6, 'text-anchor': 'end',
        'font-size': 7, fill: '#31548c', text: 'pole density I(φ)',
      }),
      svg('text', {
        x: MARGIN.left + plotWidth, y: MARGIN.top + 15, 'text-anchor': 'end',
        'font-size': 7, fill: '#d66a2b', text: 'volume fraction I(φ)·sin φ',
      }),
    ],
  );
}
