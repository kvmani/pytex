/**
 * Controls built from the manifest, never by hand.
 *
 * A panel does not write a form. It hands this module an operation's manifest
 * entry and gets back a form plus a `values()` reader. Label, units, default,
 * placeholder, validation hint and the inline help popover all come from the
 * same declaration Python validates against, which is what stops help text and
 * behaviour from drifting apart.
 *
 * Parameters are grouped by their declared `group`, and anything marked
 * `advanced` hides behind a disclosure — the "dense but not cluttered" rule
 * applied mechanically instead of panel by panel.
 */

import { append, clear, el, markdown } from './dom.js';
import * as log from './logbook.js';
import { phaseControl } from './phasecontrol.js';

/**
 * Build a form for one operation.
 *
 * @param {object} operation - Manifest entry.
 * @param {object} [options]
 * @param {object} [options.initial] - Starting values, e.g. from an example.
 * @param {Function} [options.onChange] - Called after any edit.
 * @returns {{element: HTMLElement, values: Function, setValues: Function, showError: Function, clearErrors: Function}}
 */
export function buildForm(operation, { initial = {}, onChange = () => {} } = {}) {
  const fields = new Map();
  const groups = new Map();
  const root = el('form.controls', {
    onsubmit: (event) => event.preventDefault(),
  });

  const ungrouped = [];
  const advanced = [];

  for (const parameter of operation.parameters ?? []) {
    const field = buildField(parameter, {
      value: initial[parameter.name] ?? parameter.default,
      onChange,
    });
    fields.set(parameter.name, field);
    if (parameter.advanced) advanced.push(field.element);
    else if (parameter.group) {
      if (!groups.has(parameter.group)) groups.set(parameter.group, []);
      groups.get(parameter.group).push(field.element);
    } else ungrouped.push(field.element);
  }

  append(root, ungrouped);
  for (const [name, nodes] of groups) {
    root.append(disclosure(name, nodes, { open: true }));
  }
  if (advanced.length) root.append(disclosure('Advanced', advanced, { open: false }));

  return {
    element: root,
    /** Current values, ready to send as request parameters. */
    values() {
      const result = {};
      for (const [name, field] of fields) result[name] = field.read();
      return result;
    },
    /**
     * One parameter's control, by name, or `undefined`.
     *
     * For a panel that presents a parameter somewhere of its own — a slider
     * under the picture it changes rather than a box in the rail — and needs to
     * hide the generated control so the two cannot disagree about the value.
     *
     * @param {string} name - The parameter name, as the manifest declares it.
     */
    field(name) {
      return fields.get(name);
    },
    /** Load a set of values, e.g. when an example is chosen. */
    setValues(values) {
      for (const [name, field] of fields) {
        if (name in values) field.write(values[name]);
      }
    },
    /** Put a server error beside the control it names, if it named one. */
    showError(error) {
      this.clearErrors();
      const field = error.field ? fields.get(error.field) : null;
      if (!field) return false;
      field.setError(error.message);
      return true;
    },
    clearErrors() {
      for (const field of fields.values()) field.setError(null);
    },
  };
}

function disclosure(title, nodes, { open }) {
  const body = el('div.group__body', {}, nodes);
  return el('details.group', { open }, [el('summary', { text: title }), body]);
}

function buildField(parameter, { value, onChange }) {
  const errorNode = el('p.field__error', { hidden: true });
  const helpNode = el('p.field__help', {}, markdown(parameter.help));
  helpNode.hidden = true;

  const helpButton = el('button.field__help-button', {
    type: 'button',
    text: '?',
    'aria-expanded': 'false',
    'aria-label': `What is ${parameter.label}?`,
    onclick: () => {
      helpNode.hidden = !helpNode.hidden;
      helpButton.setAttribute('aria-expanded', String(!helpNode.hidden));
    },
  });

  const input = buildInput(parameter, value, onChange);

  // A checkbox puts its label beside the box, so repeating it in the field
  // header would say the same thing twice; the help button moves inline instead.
  const inlineLabel = parameter.kind === 'boolean';
  const label = inlineLabel
    ? null
    : el('label.field__label', { for: input.id }, [
        parameter.label,
        parameter.units ? el('span.field__units', { text: parameter.units }) : null,
        helpButton,
      ]);
  const body = inlineLabel
    ? el('div.field__row', {}, [input.element, helpButton])
    : input.element;

  const element = el('div.field', {}, [label, body, helpNode, errorNode]);

  const field = {
    element,
    read: input.read,
    write: input.write,
    setError(message) {
      errorNode.hidden = !message;
      errorNode.textContent = message ?? '';
    },
  };

  if (input.problem) {
    // Reported on the transition into the bad state only. Logging every
    // keystroke of a mistyped field would push the rest of the session out of a
    // console the user opened to find out what went wrong.
    let reported = false;
    input.element.addEventListener('input', () => {
      const problem = input.problem();
      field.setError(problem);
      if (problem && !reported) {
        log.warning(problem, {
          source: 'controls',
          detail: { field: parameter.name, label: parameter.label },
        });
      }
      reported = Boolean(problem);
    });
  }

  return field;
}

let idCounter = 0;
function nextId(name) {
  idCounter += 1;
  return `ctl-${name}-${idCounter}`;
}

function buildInput(parameter, value, onChange) {
  const id = nextId(parameter.name);
  switch (parameter.kind) {
    case 'number':
    case 'integer':
      return numberInput(parameter, value, onChange, id);
    case 'boolean':
      return booleanInput(parameter, value, onChange, id);
    case 'choice':
      return choiceInput(parameter, value, onChange, id);
    case 'indices':
      return indicesInput(parameter, value, onChange, id, { multi: false });
    case 'indices-list':
      return indicesInput(parameter, value, onChange, id, { multi: true });
    case 'object':
      return objectInput(parameter, value, onChange, id);
    default:
      return textInput(parameter, value, onChange, id);
  }
}

function numberInput(parameter, value, onChange, id) {
  const node = el('input', {
    id,
    type: 'number',
    value: value ?? '',
    step: parameter.kind === 'integer' ? 1 : (parameter.step ?? 'any'),
    min: parameter.minimum,
    max: parameter.maximum,
    oninput: onChange,
  });
  return {
    id,
    element: node,
    read: () => (node.value === '' ? null : Number(node.value)),
    write: (next) => {
      node.value = next ?? '';
    },
    /**
     * Why this cannot be left to the server.
     *
     * A `type="number"` input whose content is not a number reports an *empty*
     * value: the characters are on screen, and `node.value` is `''`. Sent as-is
     * that becomes "this required field is missing", which is a confusing
     * answer to a field the user can see they filled in. The browser does say
     * what happened, through `validity.badInput`, so the check happens here and
     * the message names the real problem.
     *
     * Everything else stays server-side, where the same declaration that
     * generated this control does the validating.
     */
    problem: () => {
      if (!node.validity.badInput) return null;
      return parameter.kind === 'integer'
        ? 'Invalid format of the input: only integers are allowed!'
        : 'Invalid format of the input: only numbers are allowed!';
    },
  };
}

function textInput(parameter, value, onChange, id) {
  const node = parameter.multiline
    ? el('textarea', { id, placeholder: parameter.placeholder, oninput: onChange })
    : el('input', { id, type: 'text', placeholder: parameter.placeholder, oninput: onChange });
  node.value = value ?? '';
  return {
    id,
    element: node,
    read: () => (node.value.trim() === '' ? null : node.value),
    write: (next) => {
      node.value = next ?? '';
    },
  };
}

function booleanInput(parameter, value, onChange, id) {
  const node = el('input', { id, type: 'checkbox', onchange: onChange });
  node.checked = Boolean(value);
  const wrapper = el('div.checkbox', {}, [node, el('span', { text: parameter.label })]);
  return {
    id,
    element: wrapper,
    read: () => node.checked,
    write: (next) => {
      node.checked = Boolean(next);
    },
  };
}

function choiceInput(parameter, value, onChange, id) {
  const node = el(
    'select',
    { id, onchange: onChange },
    (parameter.options ?? []).map((option) =>
      el('option', { value: option.value, text: option.label, title: option.help }),
    ),
  );
  if (value !== null && value !== undefined) node.value = value;
  return {
    id,
    element: node,
    read: () => node.value,
    write: (next) => {
      if (next !== null && next !== undefined) node.value = next;
    },
  };
}

/**
 * Miller indices: one box per index, never a free-text row.
 *
 * Why the free-text row went
 * --------------------------
 * A single box invites `110` for what the user means as `1 1 0`, and it was
 * being typed. The server used to guess that one correctly, which is the worst
 * possible outcome: it taught the habit, and the habit fails silently the
 * moment an index reaches two digits, because `10 10 0` typed as `10100` is
 * `(10, 10, 0)`, `(1, 0, 100)` and `(101, 0, 0)` with nothing to choose
 * between them. A calculation then returns a plausible number for indices
 * nobody entered.
 *
 * So the control no longer has a place to put a run of digits. Each index gets
 * a box of its own, named `h`, `k`, `l` — or `u`, `v`, `w`, or whatever the
 * parameter's label says its indices are called, which Python decides and
 * publishes as `symbols`. A list parameter is a stack of the same rows with an
 * add and a remove button.
 *
 * What is still the server's job
 * ------------------------------
 * Everything numerical. The boxes send what was typed, index by index, and the
 * same validator the API uses decides whether `1 1 0` is a legal row for this
 * parameter. The only judgement made here is "a row you have started must be
 * finished", which the server cannot phrase as well because by then the
 * separate boxes are one list with a hole in it.
 */
function indicesInput(parameter, value, onChange, id, { multi }) {
  const width = parameter.width ?? 3;
  const symbols = normaliseSymbols(parameter.symbols, width);
  const rows = [];
  const rowsNode = el('div.indices__rows');

  const container = el(multi ? 'div.indices.indices--multi' : 'div.indices', {
    role: 'group',
    'aria-label': parameter.label,
  });

  /** One row of `width` boxes, optionally removable. */
  function addRow(rowValue, { focus = false } = {}) {
    const boxes = symbols.map((symbol, position) =>
      el('input.indices__box', {
        // The label points at the first box of the first row, so that clicking
        // the parameter's label lands somewhere sensible.
        id: position === 0 && rows.length === 0 ? id : undefined,
        type: 'text',
        inputmode: 'numeric',
        autocomplete: 'off',
        spellcheck: 'false',
        size: 3,
        placeholder: symbol,
        title: `${parameter.label}: ${symbol}`,
        'aria-label': `${parameter.label}: ${symbol}`,
        oninput: onChange,
      }),
    );
    const cells = splitRow(rowValue, width);
    boxes.forEach((box, position) => {
      box.value = cells[position] ?? '';
    });

    const row = {
      boxes,
      read: () => boxes.map((box) => box.value.trim()),
      isEmpty: () => boxes.every((box) => box.value.trim() === ''),
    };
    const remove = multi
      ? el('button.indices__remove', {
          type: 'button',
          text: '×',
          title: 'Remove this row',
          'aria-label': 'Remove this row',
          onclick: () => {
            const index = rows.indexOf(row);
            if (index >= 0) rows.splice(index, 1);
            row.element.remove();
            if (!rows.length) addRow(null);
            renumber();
            onChange();
          },
        })
      : null;

    row.element = el('div.indices__row', {}, [
      multi ? el('span.indices__ordinal', { text: String(rows.length + 1) }) : null,
      ...interleave(boxes),
      remove,
    ]);
    rows.push(row);
    rowsNode.append(row.element);
    renumber();
    if (focus) boxes[0].focus();
    return row;
  }

  function renumber() {
    if (!multi) return;
    rows.forEach((row, index) => {
      const ordinal = row.element.querySelector('.indices__ordinal');
      if (ordinal) ordinal.textContent = String(index + 1);
    });
  }

  function setRows(next) {
    rows.length = 0;
    clear(rowsNode);
    const incoming = splitRows(next, width, multi);
    if (!incoming.length) incoming.push(null);
    for (const row of incoming) addRow(row);
  }

  setRows(value);
  container.append(rowsNode);
  if (multi) {
    container.append(
      el('button.indices__add', {
        type: 'button',
        text: '+ Add row',
        onclick: () => {
          if (parameter.maxRows && rows.length >= parameter.maxRows) return;
          addRow(null, { focus: true });
        },
      }),
    );
  }

  return {
    id,
    element: container,
    read: () => {
      const filled = rows.filter((row) => !row.isEmpty()).map((row) => row.read());
      if (!filled.length) return null;
      return multi ? filled : filled[0];
    },
    write: (next) => setRows(next),
    /**
     * The one check that belongs here: a half-filled row.
     *
     * Sent as it stands, an empty box arrives as an empty string and the server
     * answers "'' is not a whole number", which names a box the user cannot see
     * a name for. Said here, it names the index that is missing.
     */
    problem: () => {
      const partial = rows.find((row) => !row.isEmpty() && row.read().some((cell) => cell === ''));
      if (!partial) return null;
      const missing = symbols.filter((_, position) => partial.read()[position] === '');
      return `Every index needs a value: ${missing.join(', ')} ${
        missing.length === 1 ? 'is' : 'are'
      } still empty.`;
    },
  };
}

/** Put the separators between the boxes, so a row reads `h k l` on screen. */
function interleave(boxes) {
  return boxes.flatMap((box, position) =>
    position === 0 ? [box] : [el('span.indices__gap', { 'aria-hidden': 'true' }), box],
  );
}

/** `symbols` comes from the manifest; fall back if a build predates it. */
function normaliseSymbols(symbols, width) {
  if (Array.isArray(symbols) && symbols.length === width) return symbols.map(String);
  const defaults = { 3: ['h', 'k', 'l'], 4: ['h', 'k', 'i', 'l'] }[width];
  return defaults ?? Array.from({ length: width }, (_, index) => `i${index + 1}`);
}

/**
 * One row's cells, from whatever a default or an example supplies.
 *
 * A stored value may be `[1, 1, 0]`, or the text `"1 1 0"` an older example
 * carries. Text is split on separators; a run of digits is *not* split, because
 * splitting it here would reintroduce exactly the guess the server now refuses.
 * It lands in the first box instead, where its owner can see it and fix it.
 */
function splitRow(value, width) {
  if (value === null || value === undefined) return [];
  if (Array.isArray(value)) return value.slice(0, width).map((cell) => String(cell));
  return String(value)
    .trim()
    .replace(/^[([{<]|[)\]}>]$/g, '')
    .split(/[\s,]+/)
    .filter((cell) => cell !== '')
    .slice(0, width);
}

/** The rows of a value: one for a single set, any number for a list. */
function splitRows(value, width, multi) {
  if (value === null || value === undefined) return [];
  if (!multi) return [value];
  if (typeof value === 'string') {
    return value
      .replace(/;/g, '\n')
      .split('\n')
      .filter((line) => line.trim() !== '');
  }
  if (!Array.isArray(value)) return [value];
  return value;
}

function objectInput(parameter, value, onChange, id) {
  if (parameter.editor === 'phase') return phaseControl(parameter, value, onChange, id);
  const node = el('textarea', { id, rows: 4, oninput: onChange, spellcheck: 'false' });
  node.value = value ? JSON.stringify(value, null, 2) : '';
  return {
    id,
    element: node,
    read: () => {
      if (node.value.trim() === '') return null;
      try {
        return JSON.parse(node.value);
      } catch {
        // Send it as-is and let the server produce the error message; a
        // client-side parse error here would be a second, differently worded
        // complaint about the same input.
        return node.value;
      }
    },
    write: (next) => {
      node.value = next ? JSON.stringify(next, null, 2) : '';
    },
  };
}

/** Replace a container's contents with a form. */
export function mountForm(container, form) {
  clear(container);
  container.append(form.element);
  return form;
}
