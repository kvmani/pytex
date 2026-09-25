/**
 * FIB lamella planning: which grain to cut, where, and at what azimuth.
 *
 * A sub-tab of the EBSD workspace because its input is the scan the workspace
 * has open: the same file, the same practice dataset, chosen once in
 * `core/ebsdscan.js`. Almost nothing here is about crystallography. The form is
 * generated from the `fib.lamella_plan` declaration, the result is drawn by the
 * shared renderer (highlights, warnings, the ranked table, the staged figures),
 * and this module adds the two things only this panel needs: a headline strip
 * that states the answer in the operator's words, and the printable work order,
 * which the server writes as one self-contained HTML page.
 *
 * Why the work order is opened rather than rendered inline
 * -------------------------------------------------------
 * It is the page an operator carries to the FIB, so it is laid out for A4 and
 * for a printer, not for this column. Opening it in its own window and printing
 * from there gives exactly that page; if the browser refuses the window, it is
 * downloaded instead, so the button never silently does nothing.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber } from '../core/dom.js';
import { activeDataset, adoptForm, scanControls, withScan } from '../core/ebsdscan.js';
import { download, renderResult } from '../core/result.js';

export const panel = {
  id: 'fib_lamella',
  title: 'FIB lamella',
  tagline: 'Which grain to cut, where, and at what azimuth for a TEM zone axis.',
};

export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'fib.lamella_plan',
  );
  const examples = context.manifest.examples.filter(
    (entry) => entry.operation === operation.id,
  );

  const state = { result: null, form: null, teaches: null };
  const headline = el('div.summary-cards', { 'data-role': 'fib-headline' });
  const details = el('div');
  const formHost = el('div');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Plan the lamella',
    onclick: () => run(),
  });
  const printButton = el('button.button.button--block', {
    type: 'button',
    text: 'Print work order',
    disabled: true,
    title: 'Open the one-page work order for the FIB operator and print it.',
    onclick: () => printWorkOrder(),
  });
  const saveButton = el('button.button.button--block', {
    type: 'button',
    text: 'Download work order',
    disabled: true,
    title: 'Save the work order as a self-contained HTML page.',
    onclick: () => saveWorkOrder(),
  });

  const scan = scanControls({
    operation,
    onChange: () => run(),
    showError: context.showError,
  });

  context.rail.append(
    scan.element,
    formHost,
    runButton,
    printButton,
    saveButton,
    examples.length
      ? el('details.group', { open: true }, [
          el('summary', { text: 'Try an example' }),
          el('div.group__body', {}, [
            el(
              'div.examples',
              {},
              examples.map((example) =>
                el(
                  'button.example',
                  { type: 'button', onclick: () => loadExample(example, { chosen: true }) },
                  [el('strong', { text: example.title }), el('span', { text: example.summary })],
                ),
              ),
            ),
          ]),
        ])
      : null,
  );

  context.stage.append(headline, details);

  function renderControls(initial = {}, { chosen = false } = {}) {
    const dataset = activeDataset();
    state.form = buildForm(operation, {
      initial: chosen ? { dataset, ...initial } : { ...initial, dataset },
    });
    formHost.replaceChildren(state.form.element);
    adoptForm(state.form, { adoptDataset: chosen });
  }

  function loadExample(example, { chosen = false } = {}) {
    state.teaches = example.teaches;
    renderControls(example.request, { chosen });
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Planning…';
    state.form.clearErrors();
    try {
      state.result = await call(operation.id, withScan(state.form.values()));
      drawHeadline();
      renderResult(details, state.result, { teaches: state.teaches });
      state.teaches = null;
      printButton.disabled = false;
      saveButton.disabled = false;
    } catch (error) {
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
      runButton.textContent = 'Plan the lamella';
    }
  }

  /** The answer in four figures, each with the condition it holds under. */
  function drawHeadline() {
    const report = state.result.data.report;
    const best = report.plans.find((plan) => plan.grain_id === state.result.data.best_grain_id)
      ?? report.plans[0];
    const figures = [
      ['Grain', String(best.grain_id ?? '—'), best.member_text + ' on the beam'],
      [
        'Residual tilt',
        `${formatNumber(best.geometry.eps_deg, 2)} ± ${formatNumber(best.uncertainty.expanded_deg, 2)}°`,
        `the holder must supply this; ${best.feasibility}`,
      ],
      [
        'Lamella normal',
        `${formatNumber(best.geometry.theta_sample_deg, 2)}°`,
        `sample frame; ${formatNumber(best.theta_image_deg, 2)}° in the SEM image`,
      ],
      [
        'FIB rotation',
        `${formatNumber(best.theta_ion_deg, 2)}°`,
        best.chamber.calibrated ? 'from the calibrated chamber' : 'UNCALIBRATED — see warnings',
      ],
    ];
    headline.replaceChildren(
      el(
        'div.summary-headline',
        {},
        figures.map(([label, value, note]) =>
          el('div.summary-figure', {}, [
            el('span.summary-figure__label', { text: label }),
            el('span.summary-figure__value', { text: value }),
            el('span.summary-figure__note', { text: note }),
          ]),
        ),
      ),
    );
  }

  function workOrderName() {
    const grain = state.result?.data?.best_grain_id;
    return `fib_work_order${grain === null || grain === undefined ? '' : `_grain_${grain}`}.html`;
  }

  function saveWorkOrder() {
    if (!state.result) return;
    download(workOrderName(), state.result.data.work_order_html, 'text/html');
  }

  function printWorkOrder() {
    if (!state.result) return;
    const page = window.open('', '_blank');
    if (!page) {
      saveWorkOrder();
      return;
    }
    page.document.open();
    page.document.write(state.result.data.work_order_html);
    page.document.close();
    page.focus();
    page.print();
  }

  renderControls();
  if (examples.length) loadExample(examples[0]);
  else run();

  return { help: () => operation };
}
