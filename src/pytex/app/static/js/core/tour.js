/**
 * The welcome, and the short tour behind it.
 *
 * The problem this solves
 * -----------------------
 * A colleague opens the workbench for the first time and sees seven tabs, a
 * control rail, and a figure. Nothing on screen says that the search box
 * indexes every operation and every worked example, that each control carries
 * the same help the API documentation does, or that the strip along the bottom
 * is a record of the session they can read afterwards. Those are the three
 * things that change how the application feels, and all three were invisible.
 *
 * The rules it follows
 * --------------------
 * - **Skip is on every step**, not only the first, and it is a real button
 *   rather than a small grey link. Somebody who came here to do a calculation
 *   at four in the afternoon should be able to leave in one click.
 * - **Skipping is remembered**, so the greeting is a first-visit event and not
 *   a toll paid on every load. A demonstration machine can ask for it every
 *   time through `tour.show_every_visit`.
 * - **It never blocks the application.** Escape closes it, clicking the
 *   backdrop closes it, and it is dismissed by anything that matters more.
 * - **It points at real things.** Each step highlights the element it is
 *   talking about, so the tour teaches the page rather than a picture of it.
 *
 * Whether it runs at all is the deployment's decision, published on
 * `/api/experience`.
 */

import { clear, el } from './dom.js';
import * as log from './logbook.js';

/** Remembers that this browser has been greeted. */
const SEEN_KEY = 'pytex-tour-seen';

/**
 * Build the tour.
 *
 * @param {object} options
 * @param {object} options.config - The `tour` block from `/api/experience`.
 * @param {Array} options.steps - `{title, body, target, note}` per step; `target`
 *   is a CSS selector, and a step whose target is absent is skipped rather than
 *   pointing at nothing.
 * @param {string} options.version - The running version, so a returning user is
 *   greeted again after an upgrade rather than never again.
 * @returns {{start: Function, close: Function, shouldGreet: Function}}
 */
export function createTour({ config, steps, version }) {
  // A step whose subject is not on screen is dropped rather than shown pointing
  // at nothing. Presence in the DOM is not enough: a deployment that turned the
  // feedback form off leaves its button in the markup with `hidden` set, and a
  // step introducing an invisible button is worse than one step fewer.
  const usable = steps.filter((step) => !step.target || isOnScreen(step.target));
  let index = 0;
  let highlighted = null;

  const backdrop = el('div.tour', {
    hidden: true,
    role: 'dialog',
    'aria-modal': 'true',
    'aria-label': 'Welcome to PyTex',
    onclick: (event) => {
      if (event.target === backdrop) finish('dismissed');
    },
  });
  /*
   * The card is a sibling of the backdrop, not a child of it, and that is
   * load-bearing rather than stylistic.
   *
   * The backdrop is positioned with a z-index, which makes it a stacking
   * context: every z-index inside it is resolved against the backdrop's own
   * layer, not against the page. The highlighted element has to sit *above*
   * the backdrop to stay bright, and a card nested inside the backdrop can
   * therefore never be above the highlight, however large a number it is
   * given. The observed symptom was a Next button that could not be clicked
   * on exactly the steps that point at the stage and the rail — the two
   * biggest things on screen. As siblings, backdrop, highlight and card are
   * three layers in one context and order as written: 60, 61, 62.
   */
  const card = el('div.tour__card');
  document.body.append(backdrop, card);

  function highlight(selector) {
    if (highlighted) highlighted.classList.remove('tour-target');
    highlighted = selector ? document.querySelector(selector) : null;
    if (highlighted) {
      highlighted.classList.add('tour-target');
      highlighted.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    }
  }

  function render() {
    const step = usable[index];
    if (!step) {
      finish('completed');
      return;
    }
    highlight(step.target);
    clear(card);
    card.append(
      el('p.tour__progress', { text: `Step ${index + 1} of ${usable.length}` }),
      el('h2.tour__title', { text: step.title }),
      ...step.body.map((paragraph) => el('p.tour__body', { text: paragraph })),
      step.note ? el('p.tour__note', { text: step.note }) : null,
      el('div.tour__actions', {}, [
        // Skip first in the source order so it is the first thing a keyboard
        // user reaches, and last in the visual order so it is not the button
        // the eye lands on. Leaving must be easy, not encouraged.
        el('button.button.tour__skip', {
          type: 'button',
          text: index === 0 ? 'Skip the tour' : 'Skip the rest',
          onclick: () => finish('skipped'),
        }),
        el('span.tour__spacer'),
        index > 0
          ? el('button.button', {
              type: 'button',
              text: 'Back',
              onclick: () => {
                index -= 1;
                render();
              },
            })
          : null,
        el('button.button.button--primary', {
          type: 'button',
          text: index === 0 ? 'Show me around' : index === usable.length - 1 ? 'Done' : 'Next',
          onclick: () => {
            index += 1;
            render();
          },
        }),
      ]),
    );
    card.querySelector('.button--primary')?.focus();
  }

  function finish(reason) {
    highlight(null);
    backdrop.hidden = true;
    card.hidden = true;
    remember(version);
    log.info(`Welcome tour ${reason}.`, { source: 'app', detail: { step: index + 1, reason } });
  }

  function start() {
    if (!usable.length) return;
    index = 0;
    backdrop.hidden = false;
    card.hidden = false;
    render();
  }

  document.addEventListener('keydown', (event) => {
    if (backdrop.hidden) return;
    if (event.key === 'Escape') {
      event.preventDefault();
      finish('dismissed');
    } else if (event.key === 'ArrowRight') {
      index += 1;
      render();
    } else if (event.key === 'ArrowLeft' && index > 0) {
      index -= 1;
      render();
    }
  });

  return {
    start,
    close: () => finish('dismissed'),
    /**
     * Whether this visitor should be greeted without being asked.
     *
     * False when the deployment turned the tour off, and false when this
     * browser has already seen this version of it — the greeting is a
     * first-visit event, not a toll on every load. The masthead's Help panel
     * can still call `start()` at any time, which is what makes remembering
     * safe rather than final.
     */
    shouldGreet() {
      if (!config?.enabled) return false;
      if (config.show_every_visit) return true;
      return recall() !== version;
    },
  };
}

function isOnScreen(selector) {
  const node = document.querySelector(selector);
  // `getClientRects` rather than `offsetParent`: the latter is also null for a
  // `position: fixed` element, and the console strip is one, so that test would
  // have silently dropped the step about the message log.
  return Boolean(node) && node.getClientRects().length > 0;
}

function recall() {
  try {
    return localStorage.getItem(SEEN_KEY);
  } catch {
    // A webview with storage disabled greets every time. That is the right way
    // round: showing a skippable welcome twice is a smaller fault than never
    // showing it at all.
    return null;
  }
}

function remember(version) {
  try {
    localStorage.setItem(SEEN_KEY, version);
  } catch {
    // As above.
  }
}
