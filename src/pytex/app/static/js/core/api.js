/**
 * The one way the frontend talks to Python.
 *
 * Every panel calls `call(operation, params)`. There is no second path, no
 * per-panel fetch, and no place where a URL is spelled out twice. That matters
 * beyond tidiness: the desktop shell and the web shell reach the same endpoint,
 * so a panel cannot accidentally work in one and not the other.
 *
 * Errors arrive as the server's own envelope — a code, a sentence written for
 * the user, and usually a hint. `ServiceCallError` carries all three, so the
 * caller can put the message beside the control that caused it rather than
 * showing a stack trace nobody can act on.
 *
 * Every envelope also carries the log records that call produced, and they are
 * handed to the console here. That is why no panel has to remember to report
 * what it just ran: reporting is a property of the one path calls take, not a
 * courtesy each caller extends. What this module logs *itself* is only what the
 * server cannot know it should — a request that never arrived, and a reply that
 * was not the agreed envelope.
 */

import * as log from './logbook.js';

/** An error the server reported deliberately, with a message meant for a person. */
export class ServiceCallError extends Error {
  constructor({ code, message, hint, details }) {
    super(message);
    this.name = 'ServiceCallError';
    this.code = code ?? 'service.error';
    this.hint = hint ?? null;
    this.details = details ?? {};
  }

  /** The parameter this error is about, when the server named one. */
  get field() {
    return this.details.field ?? null;
  }
}

/** Fetch the manifest the whole interface is built from. */
export async function fetchManifest() {
  const response = await fetch('/api/manifest', { headers: { Accept: 'application/json' } });
  if (!response.ok) {
    throw new ServiceCallError({
      code: 'manifest.unavailable',
      message: `The server did not return a manifest (HTTP ${response.status}).`,
      hint: 'Reload the page; if it persists, restart the server.',
    });
  }
  return response.json();
}

/**
 * Ask the server which shell this page is running in.
 *
 * Only file handling differs between them, and the page asks rather than
 * sniffs for a window object: a shell that changes how it saves should have to
 * say so in one place on the Python side.
 *
 * @returns {Promise<object>} `{shell, can_write_local_files, ...}`; a web shell
 *   is assumed if the route is unavailable, because that is the capability set
 *   that needs no cooperation from the host.
 */
export async function fetchShell() {
  try {
    const response = await fetch('/api/shell', { headers: { Accept: 'application/json' } });
    if (!response.ok) throw new Error(String(response.status));
    SHELL = await response.json();
  } catch {
    SHELL = { shell: 'web', can_write_local_files: false, can_read_local_paths: false };
  }
  return SHELL;
}

/**
 * Ask the server how this deployment wants to greet its users.
 *
 * Two things the *operator* decides, not the code: whether the feedback form is
 * offered and what it says, and whether a first-time visitor is shown the tour.
 * A deployment that cannot answer — an older server, a route behind a proxy
 * that drops it — gets a page with neither, which is the safe failure: an
 * unexplained modal on startup is worse than no tour.
 *
 * @returns {Promise<object>} `{feedback: {...}, tour: {...}}`.
 */
export async function fetchExperience() {
  try {
    const response = await fetch('/api/experience', { headers: { Accept: 'application/json' } });
    if (!response.ok) throw new Error(String(response.status));
    return await response.json();
  } catch {
    return { feedback: { enabled: false }, tour: { enabled: false } };
  }
}

/**
 * Send one feedback or feature-request submission.
 *
 * Not routed through `call()`, because this is not an operation: it computes
 * nothing, it is not in the manifest, and it must keep working in a build where
 * the registry failed to load — which is exactly the build somebody most wants
 * to complain about.
 *
 * @param {object} submission - message, category, and the optional fields.
 * @returns {Promise<object>} `{receipt, acknowledgement}`.
 * @throws {ServiceCallError} The server refused it, with a message for a person.
 */
export async function sendFeedback(submission) {
  let response;
  try {
    response = await fetch('/api/feedback', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
      body: JSON.stringify(submission),
    });
  } catch (error) {
    throw new ServiceCallError({
      code: 'feedback.unreachable',
      message: 'PyTex could not reach the server to send that note.',
      hint: 'Check that the server is still running, then try again.',
      details: { cause: String(error) },
    });
  }
  let payload;
  try {
    payload = await response.json();
  } catch {
    throw new ServiceCallError({
      code: 'feedback.malformed',
      message: `The server answered HTTP ${response.status} without a usable body.`,
    });
  }
  if (!response.ok || payload.ok === false) {
    throw new ServiceCallError(payload.error ?? { message: 'The note was not recorded.' });
  }
  return payload;
}

let SHELL = { shell: 'web', can_write_local_files: false, can_read_local_paths: false };
let CALL_SEQUENCE = 0;

/** What the running shell can do, as last reported by {@link fetchShell}. */
export function shell() {
  return SHELL;
}

/**
 * Invoke one operation.
 *
 * @param {string} operation - Identifier from the manifest.
 * @param {object} [params] - Parameters, validated server-side.
 * @returns {Promise<object>} The result object.
 * @throws {ServiceCallError} On any deliberate failure.
 */
export async function call(operation, params = {}) {
  const callId = ++CALL_SEQUENCE;
  log.beginCall(callId, OPERATION_TITLES.get(operation) ?? operation);
  try {
    return await invoke(operation, params);
  } finally {
    log.endCall(callId);
  }
}

/** Operation id to human title, so an in-flight call can be named on the bar. */
const OPERATION_TITLES = new Map();

/**
 * Teach this module the titles from the manifest.
 *
 * Called once at start-up. Without it a call in flight is announced by its
 * dotted id, which is the name the code uses rather than the name on the button
 * the user pressed.
 *
 * @param {object[]} operations - Manifest operation entries.
 */
export function setOperationTitles(operations) {
  OPERATION_TITLES.clear();
  for (const operation of operations ?? []) OPERATION_TITLES.set(operation.id, operation.title);
}

async function invoke(operation, params) {
  let response;
  try {
    response = await fetch('/api/call', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
      body: JSON.stringify({ operation, params }),
    });
  } catch (cause) {
    // The one failure with no server-side record, because nothing reached the
    // server. If the console stayed silent here it would show a call starting
    // and then nothing at all, which reads as a hang rather than a disconnect.
    log.critical('The PyTex server could not be reached.', {
      source: operation,
      detail: { cause: String(cause) },
    });
    throw new ServiceCallError({
      code: 'network.unreachable',
      message: 'The PyTex server could not be reached.',
      hint: 'Check that it is still running, then try again.',
      details: { cause: String(cause) },
    });
  }

  let payload;
  try {
    payload = await response.json();
  } catch {
    log.critical(`The server returned a response that is not JSON (HTTP ${response.status}).`, {
      source: operation,
    });
    throw new ServiceCallError({
      code: 'response.malformed',
      message: `The server returned a response that is not JSON (HTTP ${response.status}).`,
    });
  }

  // Both envelopes carry the call's narration, so a failed call reports what it
  // managed to do before it failed rather than only that it failed.
  log.ingest(payload.log);

  if (!payload.ok) throw new ServiceCallError(payload.error ?? {});
  return payload.result;
}
