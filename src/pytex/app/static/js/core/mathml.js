/**
 * Typeset equations: build MathML elements from a result's equation tree.
 *
 * The tree is made in Python by `pytex.app.mathml` and arrives as plain JSON:
 * `{tag, text?, font?, attrs?, children?}`. The browser renders MathML natively
 * (MathML Core), so an equation needs no script library, no font download and
 * no network — the derivation of a variant identification reads the same on an
 * offline instrument PC as anywhere else.
 *
 * Built element by element with `createElementNS`, never by parsing markup, and
 * only from the whitelist below, which mirrors `ALLOWED_TAGS` and
 * `ALLOWED_ATTRIBUTES` in Python: no string in a result can become a script, a
 * link or an event handler, whatever it contains.
 *
 * Fonts: MathML Core renders only `mathvariant="normal"` reliably, so bold
 * (matrices) and sans-serif (the transpose T) are applied as classes the
 * stylesheet defines, which also keeps them consistent with the HTML report.
 */

const MATHML_NS = 'http://www.w3.org/1998/Math/MathML';

const TAGS = new Set([
  'math', 'mrow', 'mi', 'mn', 'mo', 'mtext', 'mspace', 'msub', 'msup', 'msubsup',
  'mfrac', 'msqrt', 'mphantom', 'mover', 'munder', 'mtable', 'mtr', 'mtd',
]);

const ATTRIBUTES = new Set([
  'display', 'stretchy', 'accent', 'width', 'lspace', 'rspace', 'form', 'largeop',
]);

const LEAVES = new Set(['mi', 'mn', 'mo', 'mtext']);

const FONTS = new Set(['bold', 'upright', 'sans', 'bold-italic']);

/**
 * Build one MathML element from a tree node.
 *
 * An unknown tag is dropped rather than rendered, and says so in the console:
 * a missing symbol is visible to the reader, a silently invented one is not.
 *
 * @param {object} node
 * @returns {Element|null}
 */
export function mathNode(node) {
  if (!node || !TAGS.has(node.tag)) {
    console.warn('Equation node refused:', node?.tag);
    return null;
  }
  if (isBracketMatrix(node)) {
    // The brackets are drawn by the stylesheet: a stretchy bracket needs an
    // OpenType MATH font, which many Linux desktops lack, and there it stays
    // one line high beside a three-line matrix. Mirrors `to_markup` in Python.
    const table = mathNode(node.children[1]);
    table?.classList.add('math-bmatrix');
    return table;
  }
  const element = document.createElementNS(MATHML_NS, node.tag);
  for (const [name, value] of Object.entries(node.attrs ?? {})) {
    if (ATTRIBUTES.has(name)) element.setAttribute(name, String(value));
  }
  if (node.font && FONTS.has(node.font)) {
    element.setAttribute('mathvariant', 'normal');
    element.classList.add(`math-font--${node.font}`);
  }
  if (LEAVES.has(node.tag)) {
    element.textContent = String(node.text ?? '');
    return element;
  }
  for (const child of node.children ?? []) {
    const built = mathNode(child);
    if (built) element.append(built);
  }
  return element;
}

/** A table between stretchy square brackets: `[`, mtable, `]`. */
function isBracketMatrix(node) {
  const children = node.children ?? [];
  return node.tag === 'mrow'
    && children.length === 3
    && children[1]?.tag === 'mtable'
    && children[0]?.tag === 'mo' && children[0]?.text === '['
    && children[2]?.tag === 'mo' && children[2]?.text === ']';
}

/**
 * One displayed equation: its lead sentence, the typeset equation with its
 * number at the right margin, and the note after it.
 *
 * @param {object} equation - `{math, lead?, tag?, note?}` from a result stage.
 * @returns {HTMLElement}
 */
export function equationBlock(equation) {
  const block = document.createElement('div');
  block.className = 'equation';
  if (equation.tag) block.dataset.equation = equation.tag;
  if (equation.lead) {
    const lead = document.createElement('p');
    lead.className = 'equation__lead';
    lead.textContent = equation.lead;
    block.append(lead);
  }
  const line = document.createElement('div');
  line.className = 'equation__line';
  const body = document.createElement('div');
  body.className = 'equation__body';
  const math = mathNode(equation.math);
  if (math) body.append(math);
  line.append(body);
  if (equation.tag) {
    const tag = document.createElement('span');
    tag.className = 'equation__tag';
    tag.textContent = `(${equation.tag})`;
    line.append(tag);
  }
  block.append(line);
  if (equation.note) {
    const note = document.createElement('p');
    note.className = 'equation__note';
    note.textContent = equation.note;
    block.append(note);
  }
  return block;
}
