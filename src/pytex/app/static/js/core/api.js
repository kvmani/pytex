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
 */

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

let SHELL = { shell: 'web', can_write_local_files: false, can_read_local_paths: false };

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
  let response;
  try {
    response = await fetch('/api/call', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
      body: JSON.stringify({ operation, params }),
    });
  } catch (cause) {
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
    throw new ServiceCallError({
      code: 'response.malformed',
      message: `The server returned a response that is not JSON (HTTP ${response.status}).`,
    });
  }

  if (!payload.ok) throw new ServiceCallError(payload.error ?? {});
  return payload.result;
}
