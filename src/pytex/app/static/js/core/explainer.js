/**
 * Explanatory prose that is available rather than present.
 *
 * The problem this exists for
 * ---------------------------
 * A panel's own prose — what a view shows, what an example set is for, what a
 * group of controls does not change — was written as a paragraph pinned above
 * the controls it describes. Each one is two to five lines of rail that a user
 * reads once and then scrolls past on every visit afterwards, and there are
 * around forty of them across the workspaces. Added up, they were pushing the
 * controls themselves below the fold on panels that would otherwise fit.
 *
 * Deleting the prose is not the answer: it is the part of the interface that
 * teaches, and this application is meant to teach. So the prose stays, exactly
 * as written, one click away instead of always open.
 *
 * Why `<details>` and not a popover
 * ---------------------------------
 * The disclosure is the browser's own. It is keyboard-reachable, it announces
 * its expanded state to a screen reader, it survives a print, and it costs no
 * state in this module — which is what the field-level `?` popover in
 * `controls.js` had to write by hand and what a second hand-written popover
 * here would have to write again. Closed, it is one muted line; open, it is the
 * paragraph that used to be there.
 *
 * Where the *field* help lives
 * ----------------------------
 * Not here. A parameter's help comes from the manifest and is rendered by
 * `controls.js` behind the `?` beside its label, because it is declared next to
 * the Python that validates the parameter. This module is for the prose a panel
 * writes about itself, which has no parameter to hang from.
 */

import { el, markdown } from './dom.js';

/**
 * A collapsed explainer carrying one passage of panel prose.
 *
 * @param {string} text - The prose. Markdown is rendered, as it is in field help.
 * @param {object} [options]
 * @param {string} [options.label] - The closed line's wording. Say what the
 *   passage is *about*, not "help": "What these examples show" reads as a
 *   promise, "Help" reads as a manual.
 * @param {boolean} [options.open] - Start expanded. For the rare passage a user
 *   must read before the controls make sense at all; the default is closed,
 *   which is the point of the module.
 * @returns {HTMLElement} A `<details class="explainer">`.
 */
export function explainer(text, { label = 'What this shows', open = false } = {}) {
  return el('details.explainer', { open }, [
    el('summary.explainer__summary', {}, [
      el('span.explainer__mark', { text: 'ⓘ', 'aria-hidden': 'true' }),
      el('span', { text: label }),
    ]),
    el('div.explainer__body', {}, markdown(text)),
  ]);
}
