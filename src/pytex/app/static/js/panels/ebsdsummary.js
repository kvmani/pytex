/**
 * The EBSD scan summary: what this scan *is*, before anything is read off it.
 *
 * A page of numbers rather than a figure, and deliberately the first sub-tab
 * after the map. Everything on it is a check somebody would otherwise skip:
 * whether the step is fine enough for the microstructure, whether the indexing
 * is good *uniformly* rather than on average, whether the file was read as the
 * material it is, and whether the segmentation that everything below it depends
 * on is sensible.
 *
 * Laid out as sections rather than as one long table because they are answers to
 * different questions, and a reader looking for the step size should not have to
 * read past the grain count to find it.
 */

import { call } from '../core/api.js';
import { buildForm } from '../core/controls.js';
import { el, formatNumber } from '../core/dom.js';
import { scanControls, withScan } from '../core/ebsdscan.js';
import { renderResult } from '../core/result.js';

export const panel = {
  id: 'ebsd_summary',
  title: 'Scan summary',
  tagline: 'Points, grid, phases, indexing quality and microstructure, in one page.',
};

export function mount(context) {
  const operation = context.manifest.operations.find(
    (entry) => entry.id === 'ebsd.scan_summary',
  );
  const examples = context.manifest.examples.filter(
    (entry) => entry.operation === operation.id,
  );

  const state = { result: null, form: null, teaches: null };
  const cards = el('div.summary-cards');
  const details = el('div');
  const formHost = el('div');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Summarise the scan',
    onclick: () => run(),
  });

  const scan = scanControls({ onChange: () => run(), showError: context.showError });

  context.rail.append(
    scan.element,
    formHost,
    runButton,
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

  context.stage.append(cards, details);

  function renderControls(initial = {}) {
    state.form = buildForm(operation, { initial });
    formHost.replaceChildren(state.form.element);
    // The scan travels beside the form, as it does on every EBSD panel.
    for (const field of state.form.element.querySelectorAll('.field')) {
      if (field.querySelector('[id^="ctl-scan_file-"]')) field.hidden = true;
    }
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  async function run() {
    runButton.disabled = true;
    runButton.textContent = 'Reading…';
    state.form.clearErrors();
    try {
      state.result = await call(operation.id, withScan(state.form.values()));
      draw();
      renderResult(details, state.result, { teaches: state.teaches });
      state.teaches = null;
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
      runButton.textContent = 'Summarise the scan';
    }
  }

  /**
   * The summary as sections, with the headline numbers pulled out.
   *
   * The table below carries every row; these cards carry the four a reader
   * checks first and the sentence each of them needs to be read with. A number
   * with no statement of what it depends on — a grain count without its
   * threshold, an indexed fraction without its cut-off — is not a measurement,
   * so the dependence is on the card rather than in a footnote.
   */
  function draw() {
    const data = state.result.data;
    const rows = state.result.table?.rows ?? [];
    const groups = new Map();
    for (const row of rows) {
      if (!groups.has(row.group)) groups.set(row.group, []);
      groups.get(row.group).push(row);
    }

    const headline = [
      ['Points', Number(data.acquisition.point_count).toLocaleString(), 'measured orientations'],
      [
        'Grains',
        Number(data.microstructure.grain_count).toLocaleString(),
        `at a ${formatNumber(data.microstructure.grain_threshold_deg, 1)}° threshold`,
      ],
      data.indexed_fraction === null
        ? null
        : [
            'Indexed',
            `${formatNumber(100 * data.indexed_fraction, 1)} %`,
            `above CI ${formatNumber(data.confidence_threshold, 2)}`,
          ],
      data.microstructure.mean_equivalent_diameter_um === null
        ? null
        : [
            'Mean grain size',
            `${formatNumber(data.microstructure.mean_equivalent_diameter_um, 3)} µm`,
            'equivalent circular diameter',
          ],
    ].filter(Boolean);

    cards.replaceChildren(
      el(
        'div.summary-headline',
        {},
        headline.map(([label, value, note]) =>
          el('div.summary-figure', {}, [
            el('span.summary-figure__label', { text: label }),
            el('span.summary-figure__value', { text: value }),
            el('span.summary-figure__note', { text: note }),
          ]),
        ),
      ),
      ...[...groups.entries()].map(([name, sectionRows]) =>
        el('section.card.summary-section', {}, [
          el('div.card__header', {}, [el('h2.card__title', { text: name })]),
          el('div.card__body', {}, [
            el(
              'dl.summary-list',
              {},
              sectionRows.flatMap((row) => [
                el('dt', { text: row.metric }),
                el('dd', {}, [
                  el('strong', { text: row.value }),
                  el('span.summary-list__note', { text: row.note }),
                ]),
              ]),
            ),
          ]),
        ]),
      ),
    );
  }

  renderControls();
  // The first example rather than the bare defaults, as every other panel does:
  // the practice dataset a control defaults to is the bicrystal, which is the
  // right first map and a poor first distribution — two grains make one bar.
  if (examples.length) loadExample(examples[0]);
  else run();

  return { help: () => operation };
}
