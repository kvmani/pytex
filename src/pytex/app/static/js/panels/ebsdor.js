/**
 * The orientation relationship of a measured microstructure.
 *
 * The panel is built around one asymmetry: the *input* is a table and the
 * *answer* is a verdict. Six spin boxes cannot express more than one grain
 * pair, and one pair fits any rotation exactly — so the input here is a grid of
 * rows that can be added, removed, and pasted into wholesale, because what a
 * user has is a column of Euler angles from an indexing run rather than six
 * numbers to retype.
 *
 * The answer is laid out in the order the question is actually asked. The
 * verdict first, because that is what was wanted. Then the evidence for it: how
 * far every catalogued relationship sits from the fit, drawn as bars against
 * the naming tolerance, since "is it Burgers or Shoji-Nishiyama" is a
 * comparison of magnitudes against a threshold and that is what a bar chart
 * with a rule on it says. Then the statement in integers with its price. Then
 * the coincident planes and directions *ranked*, runners-up included, because
 * the winner alone cannot say whether it won clearly or was drawn from a
 * near-tie. Last, one row per measured pair, with the variant it sits on.
 *
 * Four different angles reach this screen — the scatter, the catalogue
 * distance, the rationalization cost and the clause deviation. Every one of
 * them is labelled with the word for what it measures, taken from the service's
 * own `angle_meanings`, so the panel and the operation's help cannot drift.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { buildForm } from '../core/controls.js';
import { renderResult } from '../core/result.js';
import { call } from '../core/api.js';
import { claim, offer } from '../core/handoff.js';

export const panel = {
  id: 'ebsd_or',
  title: 'OR from grains',
  tagline: 'Measured parent and product orientations in; the relationship they show out.',
};

const OPERATION = 'ebsd.or_from_grains';

/** The six columns of one grain pair, in the order they are entered. */
const COLUMNS = [
  { key: 'p1', group: 'parent', label: 'φ₁' },
  { key: 'p2', group: 'parent', label: 'Φ' },
  { key: 'p3', group: 'parent', label: 'φ₂' },
  { key: 'c1', group: 'child', label: 'φ₁' },
  { key: 'c2', group: 'child', label: 'Φ' },
  { key: 'c3', group: 'child', label: 'φ₂' },
];

export function mount(context) {
  const operation = context.manifest.operations.find((entry) => entry.id === OPERATION);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const state = {
    form: null,
    result: null,
    teaches: null,
    /** The grain-pair rows, as strings: what the user typed, not what it parsed to. */
    rows: [],
  };

  const grid = el('div.grains__grid');
  const details = el('div');
  const answer = el('div.answer');
  const formHost = el('div');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Determine the relationship',
    onclick: () => run(),
  });

  context.rail.append(
    el('p.field__help', {
      text:
        'Enter the measured orientations of grains of the two phases — one parent and one ' +
        'product per row — and the panel fits the relationship they show, names it if it ' +
        'is a catalogued one, and states it in integers.',
    }),
    formHost,
    runButton,
    el('details.group', { open: true }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        el('p.field__help', {
          text:
            'The first is the canonical case with exact numbers, so the answer can be ' +
            'checked; the second is the same grains with measurement noise, where the ' +
            'scatter starts to mean something.',
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

  context.stage.append(grainCard(), answer, details);

  /**
   * The controls, with the pair text hidden.
   *
   * The `pairs` parameter is declared as multiline text, which the generic form
   * builder renders as a textarea. That is the right wire format and the wrong
   * control: a user editing six numbers in a row of a table should not be
   * editing whitespace in a paragraph. The generated control is therefore
   * hidden and the grid below writes through to it, so there is exactly one
   * value and the two cannot disagree.
   */
  function renderControls(initial = {}) {
    state.form = buildForm(operation, { initial });
    formHost.replaceChildren(state.form.element);
    const field = state.form.field('pairs');
    if (field) field.element.hidden = true;
    setRows(parseText(state.form.values().pairs ?? ''));
  }

  /** The grain-pair table: the panel's real input control. */
  function grainCard() {
    return el('section.grains', {}, [
      el('div.grains__head', {}, [
        el('h2.grains__title', { text: 'Measured grain pairs' }),
        el('span.grains__note', {
          text:
            'One row per pair, in degrees. Several pairs are worth far more than one: a ' +
            'single pair fits any rotation exactly, so only a set of them has a scatter.',
        }),
      ]),
      el('div.grains__columns', {}, [
        el('span.grains__spacer', { text: '#' }),
        el('span.grains__group', { text: 'Parent grain' }),
        el('span.grains__group', { text: 'Product grain' }),
        el('span.grains__spacer', { text: '' }),
      ]),
      grid,
      el('div.grains__actions', {}, [
        el('button.button', {
          type: 'button',
          text: 'Add a row',
          onclick: () => {
            state.rows.push(['', '', '', '', '', '']);
            renderGrid();
          },
        }),
        el('button.button', {
          type: 'button',
          text: 'Paste rows…',
          title: 'Paste six numbers per line from a spreadsheet or an indexing export',
          onclick: () => pasteRows(),
        }),
        el('button.button', {
          type: 'button',
          text: 'Clear',
          onclick: () => {
            state.rows = [['', '', '', '', '', '']];
            renderGrid();
          },
        }),
      ]),
      el('div.grains__paste', { hidden: true }, [
        el('textarea.grains__pastebox', {
          rows: 5,
          placeholder: '30 40 10   167.5709 58.2280 0.9653',
          'aria-label': 'Rows to paste',
        }),
        el('div.grains__actions', {}, [
          el('button.button.button--primary', {
            type: 'button',
            text: 'Add these rows',
            onclick: () => applyPaste(),
          }),
          el('button.button', {
            type: 'button',
            text: 'Cancel',
            onclick: () => togglePaste(false),
          }),
        ]),
      ]),
    ]);
  }

  function pasteRows() {
    togglePaste(true);
    context.stage.querySelector('.grains__pastebox')?.focus();
  }

  function togglePaste(open) {
    const box = context.stage.querySelector('.grains__paste');
    if (box) box.hidden = !open;
  }

  function applyPaste() {
    const box = context.stage.querySelector('.grains__pastebox');
    const added = parseText(box?.value ?? '');
    if (added.length) {
      // Rows that were left blank are scaffolding, not data: appending under
      // them would leave a hole in the middle of the table.
      const kept = state.rows.filter((row) => row.some((value) => value.trim() !== ''));
      setRows([...kept, ...added]);
    }
    if (box) box.value = '';
    togglePaste(false);
  }

  /** Read the wire format: six numbers a line, `#` comments and blanks ignored. */
  function parseText(text) {
    const rows = [];
    for (const raw of String(text).split('\n')) {
      const line = raw.split('#')[0].trim();
      if (!line) continue;
      const fields = line.replace(/[,\t]/g, ' ').split(/\s+/).filter(Boolean);
      if (fields.length !== 6) continue;
      rows.push(fields);
    }
    return rows;
  }

  function setRows(rows) {
    state.rows = rows.length ? rows.map((row) => [...row]) : [['', '', '', '', '', '']];
    renderGrid();
  }

  /** Write the grid back to the hidden parameter. One value, one source. */
  function syncField() {
    const text = state.rows
      .filter((row) => row.some((value) => String(value).trim() !== ''))
      .map((row) => row.join(' '))
      .join('\n');
    state.form.field('pairs')?.write(text);
  }

  function renderGrid() {
    grid.replaceChildren(
      ...state.rows.map((row, index) =>
        el('div.grains__row', {}, [
          el('span.grains__index', { text: String(index + 1) }),
          ...COLUMNS.map((column, position) =>
            el('input.grains__cell', {
              type: 'text',
              inputmode: 'decimal',
              value: row[position] ?? '',
              dataset: { group: column.group },
              'aria-label': `Pair ${index + 1} ${column.group} ${column.label}`,
              title: `${column.group === 'parent' ? 'Parent' : 'Product'} ${column.label}`,
              placeholder: column.label,
              oninput: (event) => {
                state.rows[index][position] = event.currentTarget.value;
                syncField();
              },
            }),
          ),
          el('button.grains__drop', {
            type: 'button',
            text: '×',
            title: 'Remove this pair',
            'aria-label': `Remove pair ${index + 1}`,
            onclick: () => {
              state.rows.splice(index, 1);
              setRows(state.rows);
            },
          }),
        ]),
      ),
    );
    syncField();
  }

  /**
   * A pair picked on the map, if the user arrived by that gesture.
   *
   * The offer is the same one the map already makes, so a pair picked there can
   * be answered here without retyping. It seeds a *row*: the panel's whole
   * point is that one pair is not an answer, and dropping a picked pair into an
   * otherwise empty table says so by leaving room under it.
   */
  function seedFromPickedPair() {
    const offered = claim('measured-pair');
    if (!offered?.grains || offered.grains.length !== 2) return;
    const [parent, child] = offered.grains;
    const row = [
      parent.phi1_deg,
      parent.Phi_deg,
      parent.phi2_deg,
      child.phi1_deg,
      child.Phi_deg,
      child.phi2_deg,
    ].map((value) => formatNumber(value, 4));
    setRows([row]);
    answer.replaceChildren(
      el('div.notes.notes--info', {}, [
        el('strong', { text: 'Two grains picked off the map. ' }),
        el('span', {
          text:
            'They are in the table as one pair. One pair fits any rotation exactly, so add ' +
            'more before trusting the verdict.',
        }),
      ]),
    );
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    renderControls(example.request);
    run();
  }

  let runToken = 0;

  async function run() {
    const token = (runToken += 1);
    runButton.disabled = true;
    runButton.textContent = 'Fitting…';
    state.form.clearErrors();
    try {
      const result = await call(OPERATION, state.form.values());
      if (token !== runToken) return;
      state.result = result;
      renderResult(details, result, { teaches: state.teaches });
      renderAnswer(result);
    } catch (error) {
      if (token !== runToken) return;
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      if (token === runToken) {
        runButton.disabled = false;
        runButton.textContent = 'Determine the relationship';
      }
    }
  }

  /** The verdict, then the evidence for it. */
  function renderAnswer(result) {
    const data = result.data;
    answer.replaceChildren(
      verdictCard(data),
      el('div.answer__row', {}, [catalogueChart(data), statementCard(data)]),
      coincidenceCard(data),
      pairCard(data),
    );
  }

  function verdictCard(data) {
    const naming = data.naming;
    const fit = data.fit;
    const meanings = data.angle_meanings ?? {};
    const conclusive = naming.is_conclusive;
    return el('section.verdict', { dataset: { conclusive: String(conclusive) } }, [
      el('div.verdict__headline', {}, [
        el('strong', {
          text: conclusive
            ? naming.best_label
            : `${naming.best_label ?? 'No relationship'} — not conclusive`,
        }),
        el('span.verdict__tag', {
          title: conclusive
            ? 'Within tolerance, and leading the runner-up by more than the scatter.'
            : 'Either outside the tolerance, or too close to the runner-up for the data to ' +
              'separate them. Saying so is the honest answer, not a failure.',
          text: conclusive ? 'conclusive' : 'inconclusive',
        }),
      ]),
      el('div.verdict__grid', {}, [
        fact(
          'Fitted rotation',
          `${formatNumber(fit.angle_deg, 2)}° about ${fit.axis_parent.label}`,
          `The same axis is ${fit.axis_child.label} against the product basis. ` +
            `Each label is ${formatNumber(fit.axis_parent.deviation_deg, 2)}° and ` +
            `${formatNumber(fit.axis_child.deviation_deg, 2)}° from the exact axis; ` +
            'these axes are not in general rational in either basis.',
        ),
        fact(
          'Catalogue distance',
          `${formatNumber(naming.best_deviation_deg, 3)}°`,
          meanings.catalog,
        ),
        fact(
          'Scatter',
          `${formatNumber(fit.mean_residual_deg, 3)}° over ${fit.pair_count} pair(s)`,
          fit.pair_count > 1
            ? meanings.residual
            : 'Zero by construction: one pair fits one rotation exactly. Add pairs before ' +
              'reading this number as agreement.',
        ),
        fact(
          'Lead over the runner-up',
          `${formatNumber(naming.margin_deg, 3)}°`,
          'A naming is only as good as this margin compared with the scatter.',
        ),
      ]),
      el('button.button', {
        type: 'button',
        text: 'Open this relationship in the variant wall ▸',
        title:
          'Draw the named relationship as crystals: the parent and every variant it admits',
        disabled: !naming.best,
        onclick: () => handOff(data),
      }),
    ]);
  }

  function fact(label, value, title) {
    return el('div.verdict__row', title ? { title } : {}, [
      el('span.verdict__label', { text: label }),
      el('span.verdict__value', { text: value }),
    ]);
  }

  /**
   * Every catalogued relationship as a bar, with the naming tolerance drawn on.
   *
   * "Is it Burgers or Shoji-Nishiyama" is a comparison of two magnitudes
   * against a threshold, which is a bar chart with a rule on it. The numbers
   * stay in the labels, so nothing has to be measured off the picture.
   */
  function catalogueChart(data) {
    const rows = data.catalog ?? [];
    const tolerance = data.naming.tolerance_deg;
    const largest = Math.max(tolerance * 1.4, ...rows.map((row) => row.deviation_deg), 1);
    const width = 320;
    const rowHeight = 22;
    const left = 132;
    const chart = svg('svg', {
      viewBox: `0 0 ${width} ${rows.length * rowHeight + 18}`,
      preserveAspectRatio: 'xMinYMin meet',
      class: 'ladder',
      'aria-label': 'Distance from the fitted rotation to every catalogued relationship',
    });
    const scale = (value) => left + (value / largest) * (width - left - 40);
    rows.forEach((row, index) => {
      const y = index * rowHeight + 6;
      chart.append(
        svg('text', {
          x: left - 6, y: y + 11,
          'text-anchor': 'end',
          'font-size': 10,
          fill: 'var(--ink)',
          text: row.relationship,
        }),
        svg('rect', {
          x: left, y: y + 3,
          width: Math.max(scale(row.deviation_deg) - left, 1),
          height: 11,
          rx: 2,
          fill: index === 0 ? 'var(--accent)' : 'var(--ink-faint)',
          'fill-opacity': index === 0 ? 0.95 : 0.45,
        }),
        svg('text', {
          x: scale(row.deviation_deg) + 4, y: y + 12,
          'font-size': 9.5,
          fill: 'var(--ink-muted)',
          text: `${formatNumber(row.deviation_deg, 2)}°`,
        }),
      );
    });
    const rule = scale(tolerance);
    chart.append(
      svg('line', {
        x1: rule, y1: 2, x2: rule, y2: rows.length * rowHeight + 2,
        stroke: 'var(--warn, #b45309)',
        'stroke-width': 1,
        'stroke-dasharray': '3 2',
      }),
      svg('text', {
        x: rule + 3, y: rows.length * rowHeight + 12,
        'font-size': 9,
        fill: 'var(--warn, #b45309)',
        text: `tolerance ${formatNumber(tolerance, 1)}°`,
      }),
    );
    return el('section.card', {}, [
      el('h3.card__title', { text: 'How close is every known relationship?' }),
      chart,
      el('p.card__note', {
        text:
          'Distance from the fitted rotation, symmetry-reduced. The dashed rule is the ' +
          'naming tolerance; a bar to its left is close enough to be named, and the gap ' +
          'between the first two bars is what makes the naming safe.',
      }),
    ]);
  }

  function statementCard(data) {
    const statement = data.statement;
    const meanings = data.angle_meanings ?? {};
    if (!statement) {
      return el('section.card', {}, [
        el('h3.card__title', { text: 'Rational statement' }),
        el('p.card__note', {
          text:
            data.statement_note ??
            'No statement in low indices was found; the rotation stands on its own.',
        }),
      ]);
    }
    return el('section.card', {}, [
      el('h3.card__title', { text: 'Stated the way a paper states it' }),
      el('div.statement', {}, [
        el('div.statement__line', { text: `${statement.plane.parent} ∥ ${statement.plane.child}` }),
        el('div.statement__note', {
          title: meanings.clause,
          text: `clause deviation ${formatNumber(statement.plane.deviation_deg, 3)}°`,
        }),
        el('div.statement__line', {
          text: `${statement.direction.parent} ∥ ${statement.direction.child}`,
        }),
        el('div.statement__note', {
          title: meanings.clause,
          text: `clause deviation ${formatNumber(statement.direction.deviation_deg, 3)}°`,
        }),
      ]),
      el('p.card__note', {
        title: meanings.rationalization,
        text:
          `Writing the fit in integers up to ${statement.max_index} costs ` +
          `${formatNumber(statement.cost_deg, 3)}°. Compare that against the scatter: an ` +
          'idealization cheaper than the measurement noise is free, and one dearer than it ' +
          'is a claim about the material rather than a reading of it.',
      }),
    ]);
  }

  /**
   * The coincident planes and directions, ranked.
   *
   * The runners-up are the content. A single chosen pair cannot say whether it
   * won by a mile or was picked out of four candidates within a tenth of a
   * degree of each other, and those are different situations for anyone about
   * to write the relationship down.
   */
  function coincidenceCard(data) {
    const table = (rows, heading) =>
      el('div.coincide', {}, [
        el('h4.coincide__title', { text: heading }),
        el('table.table', {}, [
          el('thead', {}, [
            el('tr', {}, [
              el('th', { text: data.phases.parent }),
              el('th', { text: data.phases.child }),
              el('th', { text: 'Misfit' }),
            ]),
          ]),
          el(
            'tbody',
            {},
            rows.map((row, index) =>
              el('tr', { dataset: { rank: String(index + 1) } }, [
                el('td', { text: row.parent }),
                el('td', { text: row.child }),
                el('td', { text: `${formatNumber(row.deviation_deg, 3)}°` }),
              ]),
            ),
          ),
        ]),
      ]);
    return el('section.card', {}, [
      el('h3.card__title', { text: 'Best coincident planes and directions' }),
      el('div.answer__row', {}, [
        table(data.coincidences.planes ?? [], 'Planes'),
        table(data.coincidences.directions ?? [], 'Directions'),
      ]),
      el('p.card__note', {
        text:
          'Ranked by how far the exact image of the parent object sits from the integer ' +
          'child indices. The runners-up are shown on purpose: a winner that leads by a ' +
          'hundredth of a degree is a choice, not a finding.',
      }),
    ]);
  }

  /** One row per measured pair: its residual, and which variant it sits on. */
  function pairCard(data) {
    const rows = data.pairs ?? [];
    return el('section.card', {}, [
      el('h3.card__title', { text: 'The pairs, and the variant each one sits on' }),
      el('table.table', {}, [
        el('thead', {}, [
          el('tr', {}, [
            el('th', { text: '#' }),
            el('th', { text: 'Parent φ₁ Φ φ₂' }),
            el('th', { text: 'Product φ₁ Φ φ₂' }),
            el('th', { title: data.angle_meanings?.residual, text: 'Residual' }),
            el('th', {
              title:
                'Numbered as the fitted relationship enumerates its own variants, which is ' +
                'not a published variant table.',
              text: data.variant_count ? `Variant (of ${data.variant_count})` : 'Variant',
            }),
            el('th', { text: 'Distance to it' }),
          ]),
        ]),
        el(
          'tbody',
          {},
          rows.map((row) =>
            el('tr', {}, [
              el('td', { text: String(row.pair) }),
              el('td', { text: row.parent_euler.map((v) => formatNumber(v, 2)).join(', ') }),
              el('td', { text: row.child_euler.map((v) => formatNumber(v, 2)).join(', ') }),
              el('td', { text: `${formatNumber(row.residual_deg, 3)}°` }),
              el('td', { text: `V${row.variant}` }),
              el('td', { text: `${formatNumber(row.distance_deg, 3)}°` }),
            ]),
          ),
        ),
      ]),
      el('p.card__note', {
        text:
          'The variant is the nearest one the fitted relationship admits from this parent, ' +
          'and the distance to it is how near. A large distance means the pair does not ' +
          'belong to this parent at all — worth seeing rather than rounded away to an index.',
      }),
    ]);
  }

  /**
   * Hand the named relationship to the variant wall.
   *
   * The natural next click after identifying a relationship is to see it, and
   * the wall already draws exactly that. What crosses is the *name* and the two
   * phases, not the fitted rotation: the wall draws catalogued relationships,
   * and quietly substituting a fitted rotation for the catalogue entry would
   * mean the picture and its label were about different things.
   */
  function handOff(data) {
    // The phases come from the *form*, not from the result's inputs. The
    // service echoes back a resolved phase specification -- lattice, symmetry
    // and nothing else -- and the wall needs to draw atoms, so handing that
    // across arrives as a custom phase with no atomic basis and nothing to
    // draw. What the user chose is a catalogue entry, and that is what travels.
    const values = state.form.values();
    offer('or-catalogue', {
      relationship: data.naming.best,
      phase: values.phase,
      child_phase: values.child_phase,
      source: {
        panel: panel.id,
        pair_count: data.fit.pair_count,
        deviation_deg: data.naming.best_deviation_deg,
      },
    });
    context.openPanel('variants');
  }

  // Last, so that everything the opening run touches -- the run token included
  // -- is initialised before it starts.
  renderControls();
  seedFromPickedPair();
  run();

  return { help: () => operation };
}
