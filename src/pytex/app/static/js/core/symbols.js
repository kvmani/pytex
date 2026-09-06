/**
 * The registered symbols, as the manifest published them.
 *
 * Why this module exists
 * ----------------------
 * Almost every control in the application is generated from a parameter
 * declaration, and a declaration that names a symbol arrives with the symbol
 * already rendered — `controls.js` never looks anything up. The exception is a
 * control the frontend builds itself, and there is one that matters: the cell
 * editor inside the phase control, whose six boxes are the lattice parameters.
 * Those boxes were labelled with the words `alpha`, `beta`, `gamma`, because
 * there was nowhere to get `α β γ` from.
 *
 * Typing the Greek letters into this file instead would be the mistake
 * `pytex.core.symbols` exists to prevent: a glyph in a JavaScript file is
 * unreviewable against `docs/standards/terminology_and_symbol_registry.md`, and
 * drifts from the Python that draws the same quantity on a figure axis. So the
 * table is served with the manifest and read from here, and there remains
 * exactly one registry.
 *
 * @see pytex.core.symbols
 */

/** @type {Record<string, {text: string, latex: string, meaning: string}>} */
let table = {};

/**
 * Load the table from the manifest. Called once at startup by `main.js`.
 *
 * @param {Record<string, object>} symbols - `manifest.symbols`, or nothing.
 */
export function setSymbols(symbols) {
  table = symbols && typeof symbols === 'object' ? symbols : {};
}

/**
 * The display form of a registered symbol.
 *
 * Falls back to the name itself, which is the readable failure: an older server
 * that publishes no table, or a name this build does not know, then shows
 * `alpha` — which is what the control showed before this module existed — and
 * never an empty label or the word `undefined`.
 *
 * @param {string} name - Registered symbol name, e.g. `"alpha"`.
 * @returns {string}
 */
export function symbolText(name) {
  return table[name]?.text ?? name;
}

/**
 * The one-line meaning of a registered symbol, for a tooltip or an aria label.
 *
 * Falls back to the name for the same reason `symbolText` does. A symbol on its
 * own is not self-explanatory, so anywhere a symbol is shown this is what says
 * what it means.
 *
 * @param {string} name
 * @returns {string}
 */
export function symbolMeaning(name) {
  return table[name]?.meaning ?? name;
}
