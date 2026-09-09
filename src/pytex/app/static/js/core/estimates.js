/**
 * How long an operation has taken here before.
 *
 * Why this exists
 * ---------------
 * Most operations cannot count their own work — a single matrix solve or one
 * pass of NumPy has no loop to report from — so no honest percentage exists for
 * them. What does exist, after the second run, is a measurement of how long
 * this operation takes *on this machine, with inputs of the kind this user
 * gives it*, and that is a far better basis for "about twenty seconds left"
 * than any figure the code could invent.
 *
 * So this keeps a short history of durations per operation and reports their
 * median. Everything derived from it is presented as an estimate, never as a
 * measured fraction — the two are drawn differently and labelled differently,
 * because a reader who cannot tell them apart has been misled about which
 * number is a measurement.
 *
 * The history is per-browser and per-machine, which is the right scope: an
 * estimate learned from a workstation is wrong on a laptop, and both are wrong
 * for the next person on a shared intranet server.
 */

const STORAGE_KEY = 'pytex.operation-durations.v1';

/** Durations kept per operation. Enough to be robust to one slow outlier. */
const HISTORY_LENGTH = 5;

/**
 * Runs before an estimate is offered.
 *
 * One run is a sample, not an estimate: the first run of an operation is often
 * the slowest of any, because it is the one that imports the scientific stack.
 */
const MINIMUM_SAMPLES = 2;

/** Below this an operation is not worth a progress bar at all. */
export const SLOW_ENOUGH_TO_REPORT_MS = 900;

let cache = null;

function load() {
  if (cache) return cache;
  cache = {};
  try {
    const raw = window.localStorage.getItem(STORAGE_KEY);
    if (raw) {
      const parsed = JSON.parse(raw);
      if (parsed && typeof parsed === 'object') cache = parsed;
    }
  } catch {
    // A private window, cleared site data, or a browser refusing storage. An
    // estimate is a convenience; losing it must never break a calculation.
  }
  return cache;
}

function save() {
  try {
    window.localStorage.setItem(STORAGE_KEY, JSON.stringify(cache ?? {}));
  } catch {
    // As above: storage is best-effort.
  }
}

/**
 * Record how long one call took.
 *
 * @param {string} operation - Operation id.
 * @param {number} durationMs - Wall-clock duration of the call.
 */
export function recordDuration(operation, durationMs) {
  if (!operation || !Number.isFinite(durationMs) || durationMs <= 0) return;
  const store = load();
  const history = Array.isArray(store[operation]) ? store[operation] : [];
  history.push(Math.round(durationMs));
  store[operation] = history.slice(-HISTORY_LENGTH);
  save();
}

/**
 * The expected duration of this operation, in milliseconds, or `null`.
 *
 * The median rather than the mean: one run that happened while the machine was
 * busy elsewhere should not shift every estimate that follows it.
 *
 * @param {string} operation - Operation id.
 * @returns {{ms: number, samples: number} | null} `null` until there is enough
 *   history to say anything, which is the honest answer for a first run.
 */
export function expectedDuration(operation) {
  const history = load()[operation];
  if (!Array.isArray(history) || history.length < MINIMUM_SAMPLES) return null;
  const sorted = [...history].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  const ms =
    sorted.length % 2 === 1 ? sorted[middle] : Math.round((sorted[middle - 1] + sorted[middle]) / 2);
  return { ms, samples: sorted.length };
}

/** Forget every learned duration. Used by the tests and by a reset control. */
export function forgetDurations() {
  cache = {};
  save();
}
