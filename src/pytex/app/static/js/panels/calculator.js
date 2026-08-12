/**
 * The calculator panel.
 *
 * Almost none of this file is about crystallography, which is the point. The
 * operation list, every control, every help string and every example come from
 * the manifest; this module chooses the layout, runs the call, and hands the
 * result to the shared renderer. Adding a `calc.*` operation in Python makes it
 * appear here with no change to this file.
 *
 * The one panel-specific piece is the angle-matrix view: a table of pairs is the
 * right thing to export and the wrong thing to read when the question is "which
 * of these planes are close to each other". Both come from the same numbers.
 */

import { el, formatNumber } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { renderResult } from '../core/result.js';
import { call } from '../core/api.js';

export const panel = {
  id: 'calculator',
  title: 'Calculator',
  tagline: 'Angles, families, spacings, and orientation relationships.',
};

/**
 * Mount the panel.
 *
 * @param {object} context - `{stage, rail, manifest, showError, setBusy, openHelp}`.
 */
export function mount(context) {
  const operations = context.manifest.operations.filter((entry) => entry.panel === panel.id);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = { operation: operations[0], form: null, teaches: null };

  const chooser = el(
    'select',
    {
      'aria-label': 'Calculation',
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

  const formHost = el('div');
  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Calculate',
    onclick: () => run(),
  });

  const exampleList = el(
    'div.examples',
    {},
    examples.map((example) =>
      el(
        'button.example',
        {
          type: 'button',
          onclick: () => loadExample(example),
        },
        [el('strong', { text: example.title }), el('span', { text: example.summary })],
      ),
    ),
  );

  context.rail.append(
    el('div.field', {}, [
      el('label.field__label', { text: 'Calculation' }),
      chooser,
      el('p.field__help', { text: state.operation?.summary ?? '' }),
    ]),
    formHost,
    runButton,
    el('details.group', { open: examples.length > 0 }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text: 'Complete, runnable cases on the built-in materials. Each one lands you in a working state you can then edit.',
        }),
        exampleList,
      ]),
    ]),
  );

  function renderControls(initial = {}) {
    chooser.value = state.operation.id;
    state.form = buildForm(state.operation, { initial });
    formHost.replaceChildren(state.form.element);
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
    runButton.textContent = 'Calculating…';
    state.form.clearErrors();
    try {
      const result = await call(operation.id, state.form.values());
      renderResult(context.stage, result, {
        extra: matrixViews(result),
        teaches: state.teaches,
      });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      runButton.disabled = false;
      runButton.textContent = 'Calculate';
    }
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);

  return {
    /** The help drawer shows the operation currently selected, not the panel. */
    help: () => state.operation,
  };
}

/**
 * Render an angle matrix when the result carries one.
 *
 * Colour runs from near-parallel to perpendicular through the accent hue, so a
 * block of parallelisms in an orientation-relationship check is visible before
 * a single number is read.
 */
function matrixViews(result) {
  const matrix = result.data?.matrix;
  if (!matrix?.values_deg?.length) return [];

  const header = el(
    'tr',
    {},
    [el('th', { text: '' }), ...matrix.column_labels.map((label) => el('th', { text: label }))],
  );
  const rows = matrix.values_deg.map((values, index) =>
    el('tr', {}, [
      el('th', { text: matrix.row_labels[index] }),
      ...values.map((value) =>
        el('td.numeric', {
          text: formatNumber(value, 2),
          style: `background:${shade(value)}`,
          title: `${formatNumber(value, 4)}°`,
        }),
      ),
    ]),
  );

  return [
    el('section.card', {}, [
      el('div.card__header', {}, [
        el('h2.card__title', { text: 'Angle matrix' }),
        el('p.card__subtitle', { text: 'Degrees. Darker means closer to parallel.' }),
      ]),
      el('div.table-wrap', {}, [
        el('table.result', {}, [el('thead', {}, header), el('tbody', {}, rows)]),
      ]),
    ]),
  ];
}

function shade(angleDeg) {
  const fraction = Math.min(Math.max(Number(angleDeg) / 90, 0), 1);
  const alpha = (0.42 * (1 - fraction)).toFixed(3);
  return `color-mix(in srgb, var(--accent) ${Number(alpha) * 100}%, transparent)`;
}
