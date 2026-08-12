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

  return {
    element,
    read: input.read,
    write: input.write,
    setError(message) {
      errorNode.hidden = !message;
      errorNode.textContent = message ?? '';
    },
  };
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
 * Miller indices, typed the way people type them.
 *
 * One row for a single set, a textarea for a list — a repeatable row editor was
 * tried and is worse: superimposing a dozen planes is faster as twelve typed
 * lines than as twelve clicks on an "add row" button. Parsing is left to the
 * server, so `1 1 1`, `1,1,1` and `(1 -1 0)` all work and the error message
 * comes from the same validator the API uses.
 */
function indicesInput(parameter, value, onChange, id, { multi }) {
  const placeholder = multi
    ? Array.from({ length: 2 }, () => Array(parameter.width ?? 3).fill('1').join(' ')).join('\n')
    : Array(parameter.width ?? 3).fill('1').join(' ');
  const node = multi
    ? el('textarea', { id, placeholder, rows: 3, oninput: onChange, spellcheck: 'false' })
    : el('input', { id, type: 'text', placeholder, oninput: onChange, spellcheck: 'false' });
  node.value = formatIndices(value, multi);
  return {
    id,
    element: node,
    read: () => (node.value.trim() === '' ? null : node.value),
    write: (next) => {
      node.value = formatIndices(next, multi);
    },
  };
}

function formatIndices(value, multi) {
  if (value === null || value === undefined) return '';
  if (typeof value === 'string') return value;
  if (!multi) return Array.isArray(value) ? value.join(' ') : String(value);
  if (!Array.isArray(value)) return String(value);
  return value.map((row) => (Array.isArray(row) ? row.join(' ') : String(row))).join('\n');
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
