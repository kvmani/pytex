/**
 * The console at the bottom of the window: one stream of what the app is doing.
 *
 * Every part of the application reports here — the panels when a spot is picked
 * or a control is rejected, `api.js` when a call goes out, and Python itself,
 * whose records ride back on each call envelope and are polled from
 * `/api/log` for the events that belong to no call. The point is that there is
 * exactly one place to look. Before this existed the same session's story was
 * split between a toast that disappeared, a strip that counted calls, and a
 * terminal a desktop user never opens.
 *
 * What a good record looks like
 * -----------------------------
 * A finished sentence, addressed to the researcher, with the numbers in it:
 *
 *   `Spot 2 selected: coordinates are (452.48, 709.30) px.`
 *   `Invalid format of the input: only integers are allowed!`
 *   `50% progress. ETA: 2 min 30 sec`
 *
 * Not `handleClick fired`, and not `error`. If a record does not tell the user
 * something they could act on or would want to remember, it is `debug`.
 *
 * Levels
 * ------
 * `progress`, `debug`, `info`, `notice` ("Important"), `success`, `warning`,
 * `error`, `critical` — the same tokens and the same ordering Python uses in
 * `pytex.app.logbook`, because a console that graded severity differently from
 * the code emitting it would be actively misleading.
 *
 * Progress
 * --------
 * Ticks that share a `task` replace one another in place rather than stacking,
 * so a long simulation occupies one line that counts up instead of a thousand
 * lines that bury everything else.
 */

import { clear, el } from './dom.js';
import { expectedDuration, recordDuration, SLOW_ENOUGH_TO_REPORT_MS } from './estimates.js';

/** Severity order, lowest first. Mirrors `pytex.app.logbook.LogLevel`. */
export const LEVELS = {
  progress: { label: 'Progress', severity: 5, mark: '◔' },
  debug: { label: 'Debug', severity: 10, mark: '·' },
  info: { label: 'Info', severity: 20, mark: 'i' },
  notice: { label: 'Important', severity: 25, mark: '!' },
  success: { label: 'Success', severity: 27, mark: '✓' },
  warning: { label: 'Warning', severity: 30, mark: '⚠' },
  error: { label: 'Error', severity: 40, mark: '✕' },
  critical: { label: 'Critical', severity: 50, mark: '✕' },
};

/** Severity thresholds offered in the filter, coarsest last. */
const THRESHOLDS = [
  { value: 0, label: 'Everything' },
  { value: 20, label: 'Info and above' },
  { value: 25, label: 'Important and above' },
  { value: 30, label: 'Warnings and above' },
  { value: 40, label: 'Errors only' },
];

/**
 * How many entries the console keeps.
 *
 * A day-long session must not grow without limit, and nobody scrolls back
 * through more than this. The Python buffer is smaller still; when the console
 * notices the server dropped records it says so rather than leaving a silent
 * gap.
 */
const CAPACITY = 1000;

/** Milliseconds between polls of `/api/log` while nothing is running. */
const POLL_INTERVAL_MS = 2500;

/**
 * Milliseconds between polls while a call is in flight.
 *
 * The idle rate is set by how soon a background message should appear, which is
 * "within a couple of seconds". A progress bar has a stricter requirement: at
 * 2.5 seconds a bar advances in visible jumps and reads as stalled between
 * them. The endpoint returns a few hundred bytes and only the records since the
 * last cursor, so the faster rate costs little and only while someone is
 * actually waiting on it.
 */
const BUSY_POLL_INTERVAL_MS = 700;

const state = {
  entries: [],
  /** Entries keyed by `task`, so a progress tick can find its predecessor. */
  byTask: new Map(),
  clientSequence: 0,
  /** Highest server sequence merged, so a poll asks only for what is missing. */
  serverSequence: 0,
  threshold: 0,
  query: '',
  open: false,
  follow: true,
  /** Counts since the console was last opened, shown on the collapsed bar. */
  unseen: { warning: 0, error: 0 },
  /** Calls currently in flight, keyed by call id. */
  busy: new Map(),
  /**
   * The run currently being shown on the bar, or `null` between runs.
   *
   * `fraction` is null until the operation reports one. That distinction is the
   * whole design: a measured fraction and a fraction inferred from how long
   * this operation took last time are drawn and worded differently, because a
   * reader who cannot tell them apart has been told a guess is a measurement.
   */
  run: null,
  /** Ticks the elapsed time while a run is in flight. */
  ticker: null,
  dom: null,
  poller: null,
};

/**
 * Note that a call is in flight, so the collapsed bar can say so.
 *
 * A call's *records* arrive only when it returns — that is what makes them
 * consistent with its result — so without this the bar would sit on the
 * previous message for the whole of a slow simulation and look stalled.
 *
 * @param {string|number} id - Any value unique to this call.
 * @param {string} label - What to call it on screen.
 */
export function beginCall(id, label, operation = null) {
  state.busy.set(id, label);
  state.run = {
    id,
    label,
    operation,
    startedAt: Date.now(),
    fraction: null,
    stage: null,
    etaSeconds: null,
    // What this operation has taken here before, so a bar exists from the first
    // second even for work that cannot count itself. Null on a first run, which
    // is the honest answer rather than a made-up duration.
    expected: operation ? expectedDuration(operation) : null,
  };
  startTicking();
  render();
}

/** Note that a call has returned, whatever its outcome. */
export function endCall(id) {
  state.busy.delete(id);
  if (state.run?.id === id) {
    const elapsed = Date.now() - state.run.startedAt;
    // Only calls slow enough to have needed a bar teach the estimate. A
    // sub-second call's duration is mostly transport, and letting it into the
    // history would drag the median for the runs that actually make a user wait.
    if (state.run.operation && elapsed >= SLOW_ENOUGH_TO_REPORT_MS) {
      recordDuration(state.run.operation, elapsed);
    }
    state.run = null;
  }
  if (state.busy.size === 0) stopTicking();
  render();
}

/* ------------------------------------------------------------ the progress bar */

/**
 * Keep elapsed time and the countdown moving between records.
 *
 * Server ticks arrive at most every few hundred milliseconds and are polled at
 * a slower rate still, so without this the elapsed time would advance in visible
 * jumps and read as a stall between them.
 */
function startTicking() {
  if (state.ticker) return;
  state.ticker = setInterval(() => renderProgress(), 250);
}

function stopTicking() {
  if (state.ticker) clearInterval(state.ticker);
  state.ticker = null;
}

/**
 * What the bar should currently show.
 *
 * Three cases, and they are deliberately distinguishable on screen:
 *
 * - **measured** — the operation reported a fraction of its own work. The
 *   percentage is real and the remaining time is extrapolated from the rate
 *   observed during this run.
 * - **estimated** — it reported nothing, but this operation has run here before.
 *   The bar shows elapsed against the median of those runs, worded as an
 *   estimate, and stops short of full when it overruns rather than sitting at
 *   100% while the work continues.
 * - **unknown** — a first run of something that cannot count itself. Elapsed
 *   time only, and an indeterminate track: there is nothing else true to say.
 */
export function progressView(now = Date.now()) {
  const run = state.run;
  if (!run) return null;
  const elapsedS = Math.max((now - run.startedAt) / 1000, 0);
  if (run.fraction !== null) {
    return {
      kind: 'measured',
      fraction: run.fraction,
      elapsedS,
      remainingS: run.etaSeconds,
      label: run.stage ?? run.label,
    };
  }
  if (run.expected) {
    const expectedS = run.expected.ms / 1000;
    const raw = expectedS > 0 ? elapsedS / expectedS : 0;
    // Held below full while it overruns: a bar that reaches 100% and then keeps
    // going is the exact thing that makes people distrust progress bars.
    const fraction = Math.min(raw, 0.98);
    return {
      kind: 'estimated',
      fraction,
      elapsedS,
      remainingS: raw < 1 ? Math.max(expectedS - elapsedS, 0) : null,
      samples: run.expected.samples,
      label: run.label,
    };
  }
  return { kind: 'unknown', fraction: null, elapsedS, remainingS: null, label: run.label };
}

function renderProgress() {
  if (!state.dom?.progress) return;
  const { progress, progressFill, progressLabel, progressTimes } = state.dom;
  const view = progressView();
  if (!view) {
    progress.hidden = true;
    return;
  }
  progress.hidden = false;
  progress.dataset.kind = view.kind;
  const percent = view.fraction === null ? null : Math.round(view.fraction * 100);
  progressFill.style.width = view.fraction === null ? '100%' : `${percent}%`;
  progress.setAttribute('aria-valuenow', percent === null ? '' : String(percent));
  progress.setAttribute(
    'aria-valuetext',
    percent === null
      ? `${view.label}, ${formatDuration(view.elapsedS)} elapsed, remaining time unknown`
      : `${view.label}, ${percent}% ${view.kind === 'estimated' ? 'estimated' : 'complete'}`,
  );

  progressLabel.textContent =
    percent === null
      ? view.label
      : view.kind === 'estimated'
        ? `${view.label} — about ${percent}%`
        : `${view.label} — ${percent}%`;

  const parts = [`${formatDuration(view.elapsedS)} elapsed`];
  if (view.remainingS !== null && view.remainingS !== undefined) {
    parts.push(`about ${formatDuration(view.remainingS)} left`);
  } else if (view.kind === 'measured') {
    parts.push('finishing');
  } else {
    parts.push('time remaining unknown');
  }
  if (view.kind === 'estimated') {
    parts.push(`estimated from ${view.samples} previous run${view.samples === 1 ? '' : 's'}`);
  }
  progressTimes.textContent = parts.join(' · ');
}

/* ------------------------------------------------------------------ emitting */

/**
 * Add one record.
 *
 * @param {string} level - A key of {@link LEVELS}.
 * @param {string} message - A finished sentence for the user.
 * @param {object} [options]
 * @param {string} [options.source] - Panel id, operation id, or subsystem.
 * @param {object} [options.detail] - Structured extras shown when expanded.
 * @param {string} [options.task] - Groups ticks that should replace each other.
 * @param {number} [options.progress] - Completion fraction in `[0, 1]`.
 * @param {number} [options.etaSeconds] - Estimated remaining seconds.
 * @returns {object} The stored entry.
 */
export function record(level, message, options = {}) {
  const entry = {
    key: `c${++state.clientSequence}`,
    sequence: null,
    time: Date.now(),
    level: LEVELS[level] ? level : 'info',
    message: String(message),
    source: options.source ?? 'app',
    detail: options.detail ?? null,
    task: options.task ?? null,
    progress: typeof options.progress === 'number' ? options.progress : null,
    etaSeconds: typeof options.etaSeconds === 'number' ? options.etaSeconds : null,
  };
  push(entry);
  document.dispatchEvent(new CustomEvent('pytex:log', { detail: entry }));
  return entry;
}

const emitter = (level) => (message, options) => record(level, message, options);

/** Emit a record only a developer wants. */
export const debug = emitter('debug');
/** Emit an ordinary narration record. */
export const info = emitter('info');
/** Emit an important-but-not-wrong record. */
export const notice = emitter('notice');
/** Emit a record for something that completed as intended. */
export const success = emitter('success');
/** Emit a record for something suspect that did not stop the work. */
export const warning = emitter('warning');
/** Emit a record for work that failed. */
export const error = emitter('error');
/** Emit a record for a failure that leaves the application unusable. */
export const critical = emitter('critical');

/**
 * Emit a progress tick.
 *
 * Ticks sharing `task` replace one another, so the console shows one line that
 * counts up. Pass `etaSeconds` whenever an estimate exists: "50% progress" tells
 * a user how far along the work is, and "ETA: 2 min 30 sec" tells them whether
 * to wait for it, which is the question they actually have.
 *
 * @param {string} task
 * @param {number} fraction - Clamped to `[0, 1]`.
 * @param {object} [options]
 * @param {string} [options.message] - Overrides the generated sentence.
 * @param {number} [options.etaSeconds]
 * @param {string} [options.source]
 */
export function progress(task, fraction, options = {}) {
  const clamped = Math.min(Math.max(Number(fraction) || 0, 0), 1);
  const eta = typeof options.etaSeconds === 'number' ? options.etaSeconds : null;
  const message =
    options.message ??
    `${Math.round(clamped * 100)}% progress.${eta === null ? '' : ` ETA: ${formatDuration(eta)}`}`;
  return record('progress', message, {
    source: options.source,
    detail: options.detail,
    task,
    progress: clamped,
    etaSeconds: eta ?? undefined,
  });
}

/**
 * Merge records that came from Python.
 *
 * Server records carry a monotonic `sequence`; merging on it means a record
 * delivered both on a call envelope and by the next poll appears once.
 *
 * @param {object[]} records - Wire-form records, oldest first.
 */
export function ingest(records) {
  for (const wire of records ?? []) {
    // A progress tick from the operation now running is the measured half of
    // the bar. Ticks for anything else are still logged; they simply do not
    // move a bar that is reporting a different run.
    if (
      state.run &&
      wire.level === 'progress' &&
      typeof wire.progress === 'number' &&
      // A tick emitted before this run began belongs to a previous run of the
      // same operation, replayed from the buffer. Applying it would open a new
      // run at the percentage the last one ended on.
      (typeof wire.time !== 'number' || wire.time * 1000 >= state.run.startedAt - 1000) &&
      (wire.task === state.run.operation || wire.task === null || wire.task === undefined)
    ) {
      state.run.fraction = Math.min(Math.max(wire.progress, 0), 1);
      state.run.stage = wire.detail?.stage ?? state.run.stage;
      state.run.etaSeconds =
        typeof wire.eta_seconds === 'number' ? wire.eta_seconds : null;
    }
    if (typeof wire.sequence === 'number') {
      if (wire.sequence <= state.serverSequence) continue;
      state.serverSequence = wire.sequence;
    }
    push({
      key: `s${wire.sequence}`,
      sequence: wire.sequence ?? null,
      time: (wire.time ?? Date.now() / 1000) * 1000,
      level: LEVELS[wire.level] ? wire.level : 'info',
      message: String(wire.message ?? ''),
      source: wire.source ?? 'server',
      detail: wire.detail ?? null,
      task: wire.task ?? null,
      progress: typeof wire.progress === 'number' ? wire.progress : null,
      etaSeconds: typeof wire.eta_seconds === 'number' ? wire.eta_seconds : null,
    });
  }
  render();
}

function push(entry) {
  if (entry.task) {
    const previous = state.byTask.get(entry.task);
    if (previous) {
      // Replace in place: the tick keeps its original position in the stream so
      // a running task does not jump to the bottom on every update, which makes
      // the surrounding messages unreadable.
      const index = state.entries.indexOf(previous);
      if (index >= 0) {
        state.entries[index] = { ...entry, key: previous.key };
        state.byTask.set(entry.task, state.entries[index]);
        render();
        return;
      }
    }
    state.byTask.set(entry.task, entry);
  }
  state.entries.push(entry);
  if (state.entries.length > CAPACITY) {
    const dropped = state.entries.splice(0, state.entries.length - CAPACITY);
    for (const stale of dropped) if (stale.task) state.byTask.delete(stale.task);
  }
  if (!state.open) {
    if (entry.level === 'warning') state.unseen.warning += 1;
    else if (entry.level === 'error' || entry.level === 'critical') state.unseen.error += 1;
  }
  render();
}

/* ---------------------------------------------------------------- formatting */

/** Render a duration the way `pytex.app.logbook.format_duration` does. */
export function formatDuration(seconds) {
  if (!Number.isFinite(seconds) || seconds < 0) return 'unknown';
  const whole = Math.round(seconds);
  if (whole < 60) return `${whole} sec`;
  if (whole < 3600) return `${Math.floor(whole / 60)} min ${String(whole % 60).padStart(2, '0')} sec`;
  return `${Math.floor(whole / 3600)} hr ${String(Math.floor((whole % 3600) / 60)).padStart(2, '0')} min`;
}

function clockText(milliseconds) {
  const when = new Date(milliseconds);
  const pad = (value) => String(value).padStart(2, '0');
  return `${pad(when.getHours())}:${pad(when.getMinutes())}:${pad(when.getSeconds())}`;
}

function detailText(detail) {
  return Object.entries(detail)
    .map(([key, value]) => `${key}=${typeof value === 'object' ? JSON.stringify(value) : value}`)
    .join(', ');
}

/** The whole console as plain text, for the Copy button. */
export function asText() {
  return state.entries
    .map((entry) => {
      const head = `${clockText(entry.time)}  ${LEVELS[entry.level].label.toUpperCase().padEnd(9)} [${entry.source}] ${entry.message}`;
      return entry.detail ? `${head} (${detailText(entry.detail)})` : head;
    })
    .join('\n');
}

/* ------------------------------------------------------------------ mounting */

/**
 * Build the console and take ownership of `root`.
 *
 * @param {HTMLElement} root - The `<section>` reserved for it in the page.
 * @returns {{open: Function, close: Function}}
 */
export function mountConsole(root) {
  // The bar lives above the message-log toggle and spans the window, because a
  // person waiting should not have to know which panel they are in to find out
  // how long is left. It is the one place in the shell every operation reports
  // to, which is why it is here and not in each panel.
  const progressFill = el('div.progress__fill');
  const progressLabel = el('span.progress__label', { text: '' });
  const progressTimes = el('span.progress__times', { text: '' });
  const progress = el(
    'div.progress',
    {
      hidden: true,
      role: 'progressbar',
      'aria-valuemin': '0',
      'aria-valuemax': '100',
      'aria-label': 'Progress of the running calculation',
    },
    [
      el('div.progress__track', {}, [progressFill]),
      el('div.progress__text', {}, [progressLabel, progressTimes]),
    ],
  );

  const summary = el('span.console__summary', { text: 'Ready', role: 'status', 'aria-live': 'polite' });
  const counts = el('span.console__counts', { text: 'No messages yet' });
  const indicator = el('span.console__indicator', { 'aria-hidden': 'true' });
  const stream = el('ol.console__stream', { id: 'console-stream' });
  const empty = el('p.console__empty', { text: 'Messages appear here as you work.' });

  const toggle = el(
    'button.console__toggle',
    {
      type: 'button',
      id: 'console-toggle',
      'aria-expanded': 'false',
      'aria-controls': 'console-panel',
      'aria-label': 'Open the message log',
      onclick: () => setOpen(!state.open),
    },
    [indicator, summary, counts, el('span.console__chevron', { text: '⌃', 'aria-hidden': 'true' })],
  );

  const thresholdSelect = el(
    'select.console__filter',
    {
      id: 'console-threshold',
      'aria-label': 'Minimum severity shown',
      onchange: (event) => {
        state.threshold = Number(event.target.value);
        render();
      },
    },
    THRESHOLDS.map((option) => el('option', { value: option.value, text: option.label })),
  );

  const querySearch = el('input.console__search', {
    type: 'search',
    id: 'console-search',
    placeholder: 'Filter messages…',
    'aria-label': 'Filter messages by text',
    oninput: (event) => {
      state.query = event.target.value.trim().toLowerCase();
      render();
    },
  });

  const followButton = el('button.button', {
    type: 'button',
    id: 'console-follow',
    text: 'Follow',
    title: 'Scroll to the newest message as it arrives',
    'aria-pressed': 'true',
    onclick: () => {
      state.follow = !state.follow;
      followButton.setAttribute('aria-pressed', String(state.follow));
      if (state.follow) scrollToEnd();
    },
  });

  const panel = el('div.console__panel', { id: 'console-panel', hidden: true }, [
    el('header.console__header', {}, [
      el('div.console__title', {}, [
        el('strong', { text: 'Message log' }),
        el('span', { text: 'Everything the workbench is doing, newest last' }),
      ]),
      el('div.console__tools', {}, [
        thresholdSelect,
        querySearch,
        followButton,
        el('button.button', {
          type: 'button',
          id: 'console-copy',
          text: 'Copy',
          title: 'Copy the whole log as text',
          onclick: copyToClipboard,
        }),
        el('button.button', {
          type: 'button',
          id: 'console-clear',
          text: 'Clear',
          title: 'Discard the messages shown here',
          onclick: () => {
            state.entries = [];
            state.byTask.clear();
            info('Log cleared.', { source: 'app' });
          },
        }),
      ]),
    ]),
    stream,
    empty,
  ]);

  clear(root);
  root.append(progress, toggle, panel);
  state.dom = {
    root,
    summary,
    counts,
    indicator,
    stream,
    empty,
    panel,
    toggle,
    progress,
    progressFill,
    progressLabel,
    progressTimes,
  };
  render();
  startPolling();
  return { open: () => setOpen(true), close: () => setOpen(false) };
}

function setOpen(open) {
  state.open = open;
  if (open) state.unseen = { warning: 0, error: 0 };
  const { panel, toggle } = state.dom;
  panel.hidden = !open;
  toggle.setAttribute('aria-expanded', String(open));
  toggle.setAttribute('aria-label', `${open ? 'Close' : 'Open'} the message log`);
  render();
  if (open && state.follow) scrollToEnd();
}

function visible() {
  return state.entries.filter((entry) => {
    if (LEVELS[entry.level].severity < state.threshold) return false;
    if (!state.query) return true;
    return `${entry.message} ${entry.source}`.toLowerCase().includes(state.query);
  });
}

function render() {
  if (!state.dom) return;
  const { summary, counts, indicator, stream, empty } = state.dom;
  renderProgress();

  const latest = state.entries[state.entries.length - 1];
  const running = state.busy.size > 0;
  summary.textContent = running
    ? `Running ${[...state.busy.values()].at(-1)}…`
    : latest
      ? latest.message
      : 'Ready';
  summary.className = `console__summary console__summary--${running ? 'progress' : latest ? latest.level : 'info'}`;
  indicator.className = `console__indicator${running || latest?.level === 'progress' ? ' console__indicator--busy' : ''}`;

  const parts = [];
  if (state.unseen.error) parts.push(`${state.unseen.error} error${state.unseen.error === 1 ? '' : 's'}`);
  if (state.unseen.warning) parts.push(`${state.unseen.warning} warning${state.unseen.warning === 1 ? '' : 's'}`);
  counts.textContent = parts.length
    ? parts.join(' · ')
    : state.entries.length
      ? `${state.entries.length} message${state.entries.length === 1 ? '' : 's'}`
      : 'No messages yet';
  counts.classList.toggle('console__counts--alert', state.unseen.error > 0);

  if (!state.open) return;
  const shown = visible();
  empty.hidden = shown.length > 0;
  empty.textContent = state.entries.length
    ? 'No message matches this filter.'
    : 'Messages appear here as you work.';
  const atEnd = stream.scrollTop + stream.clientHeight >= stream.scrollHeight - 24;
  clear(stream);
  for (const entry of shown) stream.append(entryNode(entry));
  if (state.follow || atEnd) scrollToEnd();
}

function entryNode(entry) {
  const meta = LEVELS[entry.level];
  return el(`li.console__entry.console__entry--${entry.level}`, { dataset: { key: entry.key } }, [
    el('time.console__time', { text: clockText(entry.time), datetime: new Date(entry.time).toISOString() }),
    el('span.console__mark', { text: meta.mark, title: meta.label, 'aria-label': meta.label }),
    el('span.console__source', { text: entry.source, title: `Reported by ${entry.source}` }),
    el('span.console__message', {}, [
      el('span', { text: entry.message }),
      entry.progress === null
        ? null
        : el('span.console__bar', { role: 'progressbar', 'aria-valuenow': Math.round(entry.progress * 100) }, [
            el('span.console__bar-fill', { style: `width: ${(entry.progress * 100).toFixed(1)}%` }),
          ]),
      entry.detail ? el('small.console__detail', { text: detailText(entry.detail) }) : null,
    ]),
  ]);
}

function scrollToEnd() {
  const { stream } = state.dom ?? {};
  if (stream) stream.scrollTop = stream.scrollHeight;
}

async function copyToClipboard() {
  const text = asText();
  try {
    await navigator.clipboard.writeText(text);
    success(`Copied ${state.entries.length} log messages to the clipboard.`, { source: 'app' });
  } catch {
    // A webview without clipboard permission is a real deployment, and losing
    // the log to a permission prompt nobody sees would be worse than saying so.
    warning('The clipboard is not available in this window; select the log text to copy it.', {
      source: 'app',
    });
  }
}

/* ------------------------------------------------------------------- polling */

/**
 * Poll for the Python-side records that ride on no call envelope.
 *
 * Paused while the page is hidden: a minimised workbench does not need to keep
 * asking, and the sequence cursor means the first poll after it returns catches
 * up in one request rather than losing what it missed.
 */
function startPolling() {
  if (state.poller) return;
  // The first poll seeds the cursor from the server's head instead of reading
  // from it. The buffer is process-wide and outlives any one page, so replaying
  // it would open every reload with hundreds of records from earlier sessions —
  // and from other people, on a shared intranet server. The console narrates
  // *this* session; a call's own records still arrive on its envelope.
  let seeded = false;
  const tick = async () => {
    if (document.hidden) return;
    try {
      const response = await fetch(`/api/log?since=${state.serverSequence}`, {
        headers: { Accept: 'application/json' },
      });
      if (!response.ok) return;
      const payload = await response.json();
      if (!seeded) {
        seeded = true;
        state.serverSequence = Math.max(state.serverSequence, payload.latest ?? 0);
        return;
      }
      ingest(payload.records ?? []);
    } catch {
      // The console must never be the thing that reports the server is down —
      // the call the user actually made will do that, with a message about
      // their work rather than about log transport.
    }
  };
  // Rescheduled after each poll rather than fixed, so the rate can follow
  // whether anything is running without tearing down the cursor.
  const schedule = () => {
    state.poller = setTimeout(async () => {
      await tick();
      schedule();
    }, state.busy.size > 0 ? BUSY_POLL_INTERVAL_MS : POLL_INTERVAL_MS);
  };
  schedule();
  tick();
}

/** Stop polling. Used by tests; the running app polls for its whole life. */
export function stopPolling() {
  if (state.poller) clearTimeout(state.poller);
  state.poller = null;
}

/** Every entry currently held, oldest first. Used by tests. */
export function entries() {
  return state.entries.slice();
}
