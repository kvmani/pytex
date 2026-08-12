/**
 * Rendering a result, and getting it back out again.
 *
 * Because every operation returns the same shape — title, prose summary,
 * optional table, notes, citations — one renderer serves every panel, and one
 * export path serves every result. A new operation is displayable and
 * exportable the moment it is registered; neither this file nor any panel needs
 * to learn about it.
 *
 * CSV is written here rather than fetched from the server for one reason: the
 * rows are already in the browser, and a round trip to re-derive them would be
 * a second computation that could disagree with the one on screen.
 */

import { clear, el, formatNumber, markdown } from './dom.js';

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
    el('section.card', {}, [
      el('div.card__header', {}, [el('h2.card__title', { text: result.title })]),
      el('div.card__body', {}, [
        el('p.summary', { text: result.summary }),
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

  if (result.table?.rows?.length) {
    container.append(tableCard(result));
  }
}

function tableCard(result) {
  const { columns, rows, caption } = result.table;
  return el('section.card', {}, [
    el('div.card__header', {}, [
      el('h2.card__title', { text: 'Data' }),
      caption ? el('p.card__subtitle', { text: caption }) : null,
      el('div.button-row', { style: 'margin-left:auto' }, [
        el('button.button', {
          type: 'button',
          text: 'CSV',
          title: 'Download these rows as CSV',
          onclick: () => downloadCsv(result),
        }),
        el('button.button', {
          type: 'button',
          text: 'JSON',
          title: 'Download the full result, including the inputs that produced it',
          onclick: () => downloadJson(result),
        }),
        el('button.button', {
          type: 'button',
          text: 'Copy summary',
          title: 'Copy the prose summary to the clipboard',
          onclick: () => navigator.clipboard?.writeText(describe(result)),
        }),
      ]),
    ]),
    el('div.table-wrap', {}, [buildTable(columns, rows)]),
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
  if (result.notes?.length) parts.push('', ...result.notes.map((note) => `Note: ${note}`));
  if (result.citations?.length) parts.push('', `Sources: ${result.citations.join('; ')}`);
  return parts.join('\n');
}

/**
 * Write the table as CSV.
 *
 * Full precision, not display precision: the on-screen table rounds so it can be
 * read, and a file that inherited that rounding would be useless for re-plotting.
 */
export function toCsv(result) {
  const { columns, rows } = result.table;
  const header = columns.map((column) =>
    csvCell(column.units ? `${column.label} (${column.units})` : column.label),
  );
  const body = rows.map((row) => columns.map((column) => csvCell(row[column.key])));
  return [header, ...body].map((row) => row.join(',')).join('\r\n');
}

function csvCell(value) {
  if (value === null || value === undefined) return '';
  const text = String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function downloadCsv(result) {
  download(`${slug(result.title)}.csv`, toCsv(result), 'text/csv;charset=utf-8');
}

function downloadJson(result) {
  download(
    `${slug(result.title)}.json`,
    JSON.stringify(result, null, 2),
    'application/json;charset=utf-8',
  );
}

/** Trigger a browser download. The desktop shell intercepts this and writes a file. */
export function download(filename, text, mime) {
  const blob = new Blob([text], { type: mime });
  const url = URL.createObjectURL(blob);
  const anchor = el('a', { href: url, download: filename });
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

function slug(text) {
  return String(text)
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '')
    .slice(0, 60) || 'pytex-result';
}

/** Render help for an operation into the help drawer. */
export function renderHelp(container, operation) {
  clear(container);
  container.append(
    el('h2', { text: operation.title }),
    el('p', { text: operation.summary }),
    ...markdown(operation.help),
  );
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
