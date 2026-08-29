/**
 * One panel's answer, offered to another panel as its input.
 *
 * Why this exists
 * ---------------
 * Until now no panel in this application seeded another. Every workspace took
 * its input from its own controls, and the one thing that crossed a boundary —
 * the open scan — is shared *state*, not a request: every EBSD panel reads it
 * when it runs, and none of them hands it to anybody.
 *
 * Two grains picked off a map are a different shape of thing. They are a
 * one-shot request from the panel the user is looking at to a panel in another
 * workspace: *take these, and answer them*. The alternative is what the grain
 * table replaced at a smaller scale — reading six numbers off one screen and
 * typing them into another — and it is the hand transcription this repository
 * refuses everywhere else.
 *
 * The design, and why it is this small
 * ------------------------------------
 * An offer is a value under a key, claimed exactly once. Nothing subscribes and
 * nothing is notified: the offering panel navigates to the receiving one
 * immediately afterwards, so the receiver reads the offer while mounting, which
 * is the only moment it could act on it anyway. A notification channel would
 * add a second way for a panel to change under the user — panels that redraw
 * themselves because something happened in a workspace that is not on screen —
 * and that is a class of surprise worth not having.
 *
 * **Claimed once, deliberately.** A hand-off is a gesture, not a setting: the
 * user picked two grains and pressed a button, and that request is answered
 * once. Leaving it readable would make every later visit to the receiving panel
 * silently re-seed itself from a gesture made minutes ago — the same failure as
 * a map that quietly re-analyses the practice dataset, in the other direction.
 * Use `peek` where a panel needs to know an offer exists without consuming it.
 *
 * The payload is a plain object and never a live reference to another panel's
 * state, so an offer cannot change between being made and being claimed.
 */

/** Open offers, by key. At most one per key: a second offer replaces the first. */
const offers = new Map();

/**
 * Offer a payload to whichever panel claims this key next.
 *
 * @param {string} key - What is being offered, e.g. `'measured-pair'`.
 * @param {object} payload - Plain data. Structured-cloned, so the offering
 *   panel may go on mutating its own state without changing what was offered.
 */
export function offer(key, payload) {
  offers.set(key, structuredClone(payload));
}

/** Read an offer without consuming it, or null. */
export function peek(key) {
  return offers.get(key) ?? null;
}

/** Take an offer, leaving nothing behind. Returns null when there is none. */
export function claim(key) {
  const payload = offers.get(key) ?? null;
  offers.delete(key);
  return payload;
}

/** Withdraw an offer that was never claimed. */
export function withdraw(key) {
  offers.delete(key);
}
