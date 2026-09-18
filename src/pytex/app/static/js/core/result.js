/**
 * Rendering a result, and getting it back out again.
 *
 * Because every operation returns the same shape — title, prose summary,
 * optional table, notes, citations — one renderer serves every panel, and one
 * export path serves every result. A new operation is displayable and
 * exportable the moment it is registered; neither this file nor any panel needs
 * to learn about it.
 *
 * Files are produced by posting the result *back* to the server, which formats
 * the object the user is looking at. Recomputing from the inputs would be tidier
 * to write and wrong to ship: the file must be the numbers on screen, and a
 * second evaluation is a second chance to differ from them.
 */

import { clear, el, formatNumber, markdown } from './dom.js';
import { PRINT_DPI, exportMenu, fileStem } from './imageexport.js';

/**
 * The export formats, as the manifest declares them.
 *
 * Set once at start-up rather than hard-coded here, so a format added in Python
 * appears on every result in every panel without an edit in the browser — the
 * same rule the operations themselves follow. The fallback covers the case of a
 * server older than this file, where three formats is better than none.
 */
let EXPORT_FORMATS = [
  { id: 'csv', label: 'CSV', description: 'One row per entity, at full precision.' },
  { id: 'xlsx', label: 'Excel', description: 'The table, plus a sheet recording the inputs.' },
  { id: 'json', label: 'JSON', description: 'The complete result, reloadable into the app.' },
];

/** Adopt the formats the server publishes. Called once, from the shell. */
export function setExportFormats(formats) {
  if (Array.isArray(formats) && formats.length) EXPORT_FORMATS = formats;
}

/**
 * Render a result into a container.
 *
 * @param {HTMLElement} container
 * @param {object} result - An `AppResult` payload.
 * @param {object} [options]
 * @param {Node[]} [options.extra] - Panel-specific nodes placed above the table.
 * @param {string} [options.teaches] - The "what to notice" line from an example.
 */
export function renderResult(container, result, { extra = [], teaches = null } = {}) {
  clear(container);

  container.append(
    el('section.card.result-head', {}, [
      el('div.card__header', {}, [
        el('h2.card__title', { text: result.title }),
        reportButtons(result),
      ]),
      el('div.card__body', {}, [
        el('p.summary', { text: result.summary }),
        result.highlights?.length ? highlightsList(result.highlights) : null,
        result.warnings?.length ? warningsList(result.warnings) : null,
        teaches ? el('p.teaches', {}, [el('strong', { text: 'What to notice: ' }), teaches]) : null,
        result.notes?.length
          ? el('ul.notes', {}, result.notes.map((note) => el('li', { text: note })))
          : null,
        result.citations?.length
          ? el('p.citations', { text: `Sources: ${result.citations.join('; ')}` })
          : null,
      ]),
    ]),
  );

  for (const node of extra) container.append(node);

  if (result.figures?.length) {
    container.append(
      el('section.card.result-figures', {}, [
        el('div.card__header', {}, [el('h2.card__title', { text: 'Figures' })]),
        el('div.card__body.figure-grid', {}, result.figures.map((figure) => figureCard(figure))),
      ]),
    );
  }

  const sectioned = result.stages?.some((stage) => stage.section);
  if (!sectioned && result.table?.rows?.length) {
    container.append(tableCard(result));
  }

  if (result.stages?.length) {
    if (sectioned) {
      for (const node of sectionCards(result)) container.append(node);
    } else {
      container.append(stagesCard(result.stages));
    }
  }
}

/**
 * The report downloads, offered on every result whether or not it has a table.
 *
 * A result without rows still has prose, figures and provenance worth keeping;
 * before this, such a result had no way out of the page at all.
 */
function reportButtons(result) {
  const formats = EXPORT_FORMATS.filter((format) => ['md', 'zip', 'json'].includes(format.id));
  if (!formats.length) return null;
  return el('div.button-row.result-head__exports', { style: 'margin-left:auto' },
    formats.map((format) => exportButton(result, format.id, format.label, format.description)));
}

/** The answer and the few numbers that decide how far to trust it. */
function highlightsList(highlights) {
  return el('dl.highlights', {}, highlights.flatMap((metric) => {
    const value = typeof metric.value === 'number'
      ? formatNumber(metric.value)
      : formatCell(metric.value);
    return [
      el('dt', { text: metric.label, title: metric.help ?? null }),
      el('dd', {
        text: metric.units ? `${value} ${metric.units}` : value,
        title: metric.help ?? null,
      }),
    ];
  }));
}

function warningsList(warnings) {
  return el('div.result-warnings', { role: 'alert' }, [
    el('strong', {
      text: warnings.length === 1
        ? 'Check this before using the result'
        : `Check these ${warnings.length} points before using the result`,
    }),
    el('ul', {}, warnings.map((warning) => el('li', { text: warning }))),
  ]);
}

/**
 * One server-drawn figure: the picture, what is plotted, and what it implies.
 *
 * Shown as an image rather than inlined, so the figure's own identifiers and
 * styles can never collide with the page's, and so the browser's own "Copy
 * image" and "Save image as" work on it as they do on any picture.
 */
export function figureCard(figure) {
  const url = URL.createObjectURL(new Blob([figure.svg], { type: 'image/svg+xml' }));
  const image = el('img.result-figure__image', {
    src: url,
    alt: `${figure.title}. ${figure.caption ?? ''}`.trim(),
  });
  const widthPx = (figure.width_in ?? 6.4) * 96;
  const heightPx = (figure.height_in ?? 4) * 96;
  const menu = exportMenu({
    figure: () => ({
      markup: figure.svg,
      width: widthPx,
      height: heightPx,
      stem: `pytex-${fileStem(figure.key)}`,
      // Exactly 300 dpi of the size the figure was designed at.
      scale: PRINT_DPI / 96,
    }),
  });
  return el('figure.result-figure', { 'data-figure': figure.key }, [
    el('div.result-figure__header', {}, [
      el('h3.result-figure__title', { text: figure.title }),
      menu,
    ]),
    image,
    el('figcaption.result-figure__caption', {}, [
      figure.caption ? el('p', { text: figure.caption }) : null,
      figure.interpretation
        ? el('p.result-figure__reading', {}, [
            el('strong', { text: 'What it shows: ' }),
            figure.interpretation,
          ])
        : null,
    ]),
  ]);
}

/** The order a report is read in; mirrors `REPORT_SECTIONS` in `pytex.app.results`. */
const REPORT_SECTIONS = [
  ['result', 'Result in detail'],
  ['evidence', 'Evidence'],
  ['diagnostics', 'Diagnostics'],
  ['method', 'Method'],
  ['audit', 'Audit details'],
];

/**
 * The stages grouped as the report is read: result, evidence, diagnostics,
 * method, audit. The final data table belongs to the audit trail. Every stage
 * outside the audit trail starts open, and so does any stage that warns.
 */
function sectionCards(result) {
  const cards = [];
  for (const [section, heading] of REPORT_SECTIONS) {
    const members = result.stages.filter((stage) => stage.section === section);
    const withTable = section === 'audit' && result.table?.rows?.length;
    if (!members.length && !withTable) continue;
    cards.push(el(`section.card.stages.report-section.report-section--${section}`, {
      'data-section': section,
      'data-stages': String(members.length),
    }, [
      el('div.card__header', {}, [el('h2.card__title', { text: heading })]),
      el('div.card__body.stages__list', {}, members.map((stage) =>
        stageSection(stage, section !== 'audit' || stage.status === 'warning'))),
    ]));
    if (withTable) cards.push(tableCard(result));
  }
  const loose = result.stages.filter((stage) => !stage.section);
  if (loose.length) cards.push(stagesCard(loose));
  return cards;
}

/**
 * The intermediate stages behind a result, one disclosure per stage.
 *
 * A final number is only as believable as the steps that produced it, and an
 * analysis can fail at any of them — a missed peak, a wrong assignment, a
 * fit that absorbed an aberration into the cell. Each stage therefore shows
 * its own numbers, its own rows and a note on how to read them. The first
 * stage that reports a warning opens itself, because that is the one a reader
 * has to look at before believing anything below it.
 */
function stagesCard(stages) {
  const firstWarning = stages.findIndex((stage) => stage.status === 'warning');
  return el('section.card.stages', { 'data-stages': String(stages.length) }, [
    el('div.card__header', {}, [
      el('h2.card__title', { text: 'How this result was reached' }),
      el('p.card__subtitle', {
        text: `${stages.length} intermediate stage${stages.length === 1 ? '' : 's'}, in the order `
          + 'they ran. Each is included in the Markdown and Excel exports.',
      }),
    ]),
    el('div.card__body.stages__list', {}, stages.map((stage, index) =>
      stageSection(stage, index === firstWarning)),
    ),
  ]);
}

function stageSection(stage, open) {
  const metrics = stage.metrics?.length
    ? el('dl.stage__metrics', {}, stage.metrics.flatMap((metric) => {
        const value = typeof metric.value === 'number'
          ? formatNumber(metric.value)
          : formatCell(metric.value);
        return [
          el('dt', { text: metric.label, title: metric.help ?? null }),
          el('dd', { text: metric.units ? `${value} ${metric.units}` : value }),
        ];
      }))
    : null;
  const table = stage.table?.rows?.length
    ? el('div.stage__table', {}, [
        stage.table.caption ? el('p.card__subtitle', { text: stage.table.caption }) : null,
        el('div.table-wrap', {}, [buildTable(stage.table.columns, stage.table.rows)]),
      ])
    : null;
  return el(
    `details.stage.stage--${stage.status ?? 'ok'}`,
    { open, 'data-stage': stage.key },
    [
      el('summary.stage__summary', {}, [
        el('span.stage__status', {
          text: stage.status === 'warning' ? '!' : stage.status === 'info' ? 'i' : '✓',
          'aria-label': stage.status ?? 'ok',
        }),
        el('span.stage__title', { text: stage.title }),
      ]),
      el('div.stage__body', {}, [
        el('p.stage__text', { text: stage.summary }),
        metrics,
        stage.figures?.length
          ? el('div.figure-grid', {}, stage.figures.map((figure) => figureCard(figure)))
          : null,
        table,
        stage.explanation
          ? el('p.stage__explanation', {}, [
              el('strong', { text: 'How to read this: ' }),
              stage.explanation,
            ])
          : null,
      ]),
    ],
  );
}

/**
 * How many rows are put on screen before the table is truncated.
 *
 * The table is a *preview* of the export, not the export. A texture ODF is 1083
 * rows and a composite pattern several hundred, and past a couple of hundred
 * the scroll box stops being something anyone reads — it has no search, no
 * sort, and no way to reach row 900 except dragging. The export buttons above
 * it carry every row at full precision, which is what a reader with 1083 rows
 * of questions actually needs.
 *
 * The count is stated in the caption rather than silently applied: a table that
 * quietly shows a subset is worse than one that is too long, because a reader
 * counting rows would get the wrong answer and never know.
 */
const TABLE_PREVIEW_ROWS = 200;

function tableCard(result) {
  const { columns, rows, caption } = result.table;
  const truncated = rows.length > TABLE_PREVIEW_ROWS;
  const shown = truncated ? rows.slice(0, TABLE_PREVIEW_ROWS) : rows;
  const subtitle = truncated
    ? `${caption ? `${caption} ` : ''}Showing the first ${TABLE_PREVIEW_ROWS} of ` +
      `${rows.length} rows; every export below carries all ${rows.length}.`
    : caption;
  return el('section.card', {}, [
    el('div.card__header', {}, [
      el('h2.card__title', { text: 'Data' }),
      subtitle ? el('p.card__subtitle', { text: subtitle }) : null,
      el('div.button-row', { style: 'margin-left:auto' }, [
        ...EXPORT_FORMATS.map((format) =>
          exportButton(result, format.id, format.label, format.description),
        ),
        el('button.button', {
          type: 'button',
          text: 'Copy summary',
          title: 'Copy the prose summary to the clipboard',
          onclick: () => navigator.clipboard?.writeText(describe(result)),
        }),
      ]),
    ]),
    el('div.table-wrap', {}, [buildTable(columns, shown)]),
  ]);
}

function buildTable(columns, rows) {
  const head = el(
    'thead',
    {},
    el(
      'tr',
      {},
      columns.map((column) =>
        el('th', { class: column.numeric ? 'numeric' : null, title: column.help ?? null }, [
          column.label,
          column.units ? el('span.unit', { text: ` / ${column.units}` }) : null,
        ]),
      ),
    ),
  );
  const body = el(
    'tbody',
    {},
    rows.map((row) =>
      el(
        'tr',
        {},
        columns.map((column) => {
          const value = row[column.key];
          const text =
            column.numeric || typeof value === 'number'
              ? formatNumber(value, column.digits)
              : formatCell(value);
          return el('td', { class: column.numeric ? 'numeric' : null, text });
        }),
      ),
    ),
  );
  return el('table.result', {}, [head, body]);
}

function formatCell(value) {
  if (value === null || value === undefined) return '—';
  if (typeof value === 'boolean') return value ? 'yes' : 'no';
  return String(value);
}

/** The prose form: what a user pastes into a lab notebook. */
export function describe(result) {
  const parts = [result.title, '', result.summary];
  if (result.highlights?.length) {
    parts.push('', ...result.highlights.map((metric) => {
      const value = typeof metric.value === 'number'
        ? formatNumber(metric.value)
        : formatCell(metric.value);
      return `${metric.label}: ${value}${metric.units ? ` ${metric.units}` : ''}`;
    }));
  }
  if (result.warnings?.length) {
    parts.push('', ...result.warnings.map((warning) => `Warning: ${warning}`));
  }
  if (result.notes?.length) parts.push('', ...result.notes.map((note) => `Note: ${note}`));
  if (result.citations?.length) parts.push('', `Sources: ${result.citations.join('; ')}`);
  return parts.join('\n');
}

function exportButton(result, format, label, title) {
  const button = el('button.button', {
    type: 'button',
    text: label,
    title,
    onclick: async () => {
      button.disabled = true;
      try {
        await exportResult(result, format);
      } finally {
        button.disabled = false;
      }
    },
  });
  return button;
}

/**
 * Ask the server to format this result and save the reply.
 *
 * The filename comes from `Content-Disposition` so that Python owns naming:
 * a file called `d-spacings-of-nickel-fcc.xlsx` is named the same whether it
 * was produced from the browser, the desktop shell, or a script.
 */
export async function exportResult(result, format) {
  const response = await fetch('/api/export', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ result, format }),
  });
  if (!response.ok) {
    const payload = await response.json().catch(() => ({}));
    throw new Error(payload?.error?.message ?? `Export failed (HTTP ${response.status}).`);
  }
  const disposition = response.headers.get('Content-Disposition') ?? '';
  const match = /filename="([^"]+)"/.exec(disposition);
  const blob = await response.blob();
  return saveBlob(match ? match[1] : `result.${format}`, blob);
}

/** Save text produced in the page itself — a figure, a scene, a report. */
export function download(filename, text, mime) {
  return saveBlob(filename, new Blob([text], { type: mime }));
}

/**
 * Put a blob where the user can find it, by whichever route this shell has.
 *
 * There are two, and the difference is not cosmetic. A browser takes an anchor
 * with a `download` attribute and puts the file in its downloads folder. The
 * embedded web view of the desktop shell accepts that same click and does
 * nothing whatsoever with it — no file, no error — so the desktop shell hands
 * the bytes to Python, which asks where to put them and writes them there.
 *
 * Both routes announce the outcome as a `pytex:saved` event, so the shell can
 * say what happened. Silence after pressing an export button is the failure
 * mode this whole path exists to avoid.
 *
 * @param {string} filename
 * @param {Blob} blob
 * @returns {Promise<string|null>} Where it went, when the shell can say.
 */
export async function saveBlob(filename, blob) {
  const bridge = window.pywebview?.api?.save_file;
  if (bridge) {
    const buffer = await blob.arrayBuffer();
    const bytes = new Uint8Array(buffer);
    let binary = '';
    // In chunks: String.fromCharCode(...bytes) overflows the argument limit on
    // a megabyte-scale figure, which is exactly the case that must not fail.
    for (let index = 0; index < bytes.length; index += 0x8000) {
      binary += String.fromCharCode(...bytes.subarray(index, index + 0x8000));
    }
    const path = await bridge(filename, btoa(binary));
    announceSave(path ? `Saved to ${path}` : 'Not saved.', path);
    return path ?? null;
  }
  const url = URL.createObjectURL(blob);
  const anchor = el('a', { href: url, download: filename });
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
  announceSave(`${filename} downloaded.`, null);
  return null;
}

function announceSave(message, path) {
  document.dispatchEvent(
    new CustomEvent('pytex:saved', { detail: { message, path } }),
  );
}

/** Render help for an operation into the help drawer. */
export function renderHelp(container, operation) {
  clear(container);
  container.append(
    el('h2', { text: operation.title }),
    el('p', { text: operation.summary }),
    ...markdown(operation.help),
  );
  if (operation.documentation) {
    container.append(
      el('aside.docs-link', {}, [
        el('span.docs-link__eyebrow', { text: 'In-App Scientific Documentation' }),
        el('strong', { text: operation.documentation.title }),
        el('p', { text: 'Open the full guide for theory, conventions, worked examples and related APIs.' }),
        el('a.button.button--primary', {
          href: operation.documentation.url,
          target: '_blank',
          rel: 'noopener noreferrer',
          text: 'Read the full guide ↗',
        }),
      ]),
    );
  }
  if (operation.parameters?.length) {
    container.append(
      el('h3', { text: 'Inputs' }),
      el(
        'dl.detail',
        {},
        operation.parameters.flatMap((parameter) => [
          el('dt', { text: parameter.units ? `${parameter.label} (${parameter.units})` : parameter.label }),
          el('dd', { text: parameter.help }),
        ]),
      ),
    );
  }
  if (operation.returns) {
    container.append(el('h3', { text: 'Result' }), el('p', { text: operation.returns }));
  }
  if (operation.citations?.length) {
    container.append(
      el('h3', { text: 'Sources' }),
      el('ul', {}, operation.citations.map((citation) => el('li', { text: citation }))),
    );
  }
}
