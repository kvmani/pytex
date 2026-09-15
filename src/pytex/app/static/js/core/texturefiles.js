/**
 * The pole-figure files the Texture workspace has open, shared by its panels.
 *
 * There is one answer to "which measurement am I analysing": the XRDML files the
 * user opened. The texture analysis, the measured pole-figure view and the
 * Kearns pole-figure and ODF routes all read that one set, so opening the files
 * once is enough, and two panels can never quietly analyse different data.
 *
 * The files never leave the browser except as the body of the request that
 * analyses them; nothing is stored on the server.
 */

/** The open files, as `{name, text}`, in the order they were chosen. */
let items = [];

const listeners = new Set();

function notify() {
  for (const listener of listeners) listener(items);
}

/** The open files. The order matters: it is the order planes are assigned in. */
export function openedFiles() {
  return items;
}

/**
 * Read browser `File` objects and make them the open set.
 *
 * @param {File[]} files
 * @returns {Promise<object[]>}
 */
export async function openFiles(files) {
  items = await Promise.all(
    [...files].map(async (file) => ({ name: file.name, text: await file.text() })),
  );
  notify();
  return items;
}

/** Close every open file. */
export function closeFiles() {
  items = [];
  notify();
}

/**
 * Be told when the open set changes, whichever panel changed it.
 *
 * @param {(items: object[]) => void} listener
 * @returns {() => void} Unsubscribe.
 */
export function onFilesChanged(listener) {
  listeners.add(listener);
  return () => listeners.delete(listener);
}

/** One line saying what is open, for every panel's file control. */
export function describeOpenedFiles(idle) {
  return items.length
    ? `${items.length} file(s) open: ${items.map((item) => item.name).join(', ')}`
    : idle;
}
