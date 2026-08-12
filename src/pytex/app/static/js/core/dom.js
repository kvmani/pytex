/**
 * Element construction without a framework and without innerHTML.
 *
 * Every node in this application is built here. Nothing is ever assembled from
 * a template string, because service results carry user-supplied phase names
 * and file names, and a single forgotten escape in a template would put them
 * into the document as markup. Building nodes and setting `textContent` makes
 * that class of mistake unavailable rather than merely discouraged.
 */

/**
 * Create an element.
 *
 * @param {string} tag - Tag name, optionally with `.class` suffixes: `div.card.wide`.
 * @param {object} [attrs] - Attributes. `text` sets textContent; `dataset` sets
 *   data attributes; anything starting with `on` binds a listener; `null` and
 *   `undefined` values are skipped so callers can pass optional attributes.
 * @param {(Node|string|null|undefined|Array)[]} [children] - Appended in order.
 * @returns {HTMLElement}
 */
export function el(tag, attrs = {}, children = []) {
  const [name, ...classes] = tag.split('.');
  const node = document.createElement(name);
  if (classes.length) node.className = classes.join(' ');

  for (const [key, value] of Object.entries(attrs)) {
    if (value === null || value === undefined || value === false) continue;
    if (key === 'text') {
      node.textContent = String(value);
    } else if (key === 'class') {
      node.className = node.className ? `${node.className} ${value}` : String(value);
    } else if (key === 'dataset') {
      for (const [dataKey, dataValue] of Object.entries(value)) {
        if (dataValue !== null && dataValue !== undefined) node.dataset[dataKey] = String(dataValue);
      }
    } else if (key.startsWith('on') && typeof value === 'function') {
      node.addEventListener(key.slice(2).toLowerCase(), value);
    } else if (value === true) {
      node.setAttribute(key, '');
    } else {
      node.setAttribute(key, String(value));
    }
  }

  append(node, children);
  return node;
}

/** Create an SVG element. SVG needs its own namespace; `el` would produce an inert node. */
export function svg(tag, attrs = {}, children = []) {
  const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (value === null || value === undefined || value === false) continue;
    if (key === 'text') node.textContent = String(value);
    else if (key.startsWith('on') && typeof value === 'function') {
      node.addEventListener(key.slice(2).toLowerCase(), value);
    } else node.setAttribute(key, String(value));
  }
  append(node, children);
  return node;
}

/** Append children, flattening arrays and skipping empty values. */
export function append(parent, children) {
  const list = Array.isArray(children) ? children : [children];
  for (const child of list) {
    if (child === null || child === undefined || child === false) continue;
    if (Array.isArray(child)) append(parent, child);
    else parent.append(typeof child === 'string' ? document.createTextNode(child) : child);
  }
  return parent;
}

/** Remove every child of a node. */
export function clear(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
  return node;
}

/** Replace a node's children in one step. */
export function replace(node, children) {
  return append(clear(node), children);
}

/**
 * Render a restricted subset of Markdown as elements.
 *
 * Help text is authored in Python docstrings and is Markdown-flavoured. Rather
 * than ship a parser, this handles exactly what the help strings use: blank-line
 * paragraphs, `-` bullet lists, `**bold**`, `*italic*` and `` `code` ``. Anything
 * else renders as literal text, which is the safe failure.
 *
 * @param {string} text
 * @returns {HTMLElement[]}
 */
export function markdown(text) {
  const blocks = String(text ?? '').split(/\n{2,}/);
  const nodes = [];
  for (const block of blocks) {
    const lines = block.split('\n');
    if (lines.every((line) => /^\s*[-*]\s+/.test(line)) && lines.length) {
      nodes.push(el('ul', {}, lines.map((line) => el('li', {}, inline(line.replace(/^\s*[-*]\s+/, ''))))));
    } else {
      nodes.push(el('p', {}, inline(block.replace(/\n/g, ' '))));
    }
  }
  return nodes;
}

function inline(text) {
  const nodes = [];
  const pattern = /(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)/g;
  let index = 0;
  for (const match of text.matchAll(pattern)) {
    if (match.index > index) nodes.push(text.slice(index, match.index));
    const token = match[0];
    if (token.startsWith('**')) nodes.push(el('strong', { text: token.slice(2, -2) }));
    else if (token.startsWith('`')) nodes.push(el('code', { text: token.slice(1, -1) }));
    else nodes.push(el('em', { text: token.slice(1, -1) }));
    index = match.index + token.length;
  }
  if (index < text.length) nodes.push(text.slice(index));
  return nodes;
}

/**
 * Format a number for display.
 *
 * Fixed decimals when the column declares them, otherwise a general format that
 * neither prints sixteen digits of float noise nor silently rounds away a small
 * but meaningful residual — which is exactly the case for an
 * orientation-relationship deviation of a hundredth of a degree.
 */
export function formatNumber(value, digits) {
  if (value === null || value === undefined || value === '') return '—';
  if (typeof value === 'boolean') return value ? 'yes' : 'no';
  const number = Number(value);
  if (!Number.isFinite(number)) return '—';
  if (Number.isInteger(number) && digits === undefined) return String(number);
  if (digits !== undefined) return number.toFixed(digits);
  const magnitude = Math.abs(number);
  if (magnitude !== 0 && (magnitude < 1e-3 || magnitude >= 1e6)) return number.toExponential(4);
  return number.toPrecision(6).replace(/\.?0+$/, '');
}
