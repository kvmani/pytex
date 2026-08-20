/**
 * The scan the EBSD workspace is working on, shared by all of its panels.
 *
 * Why this is a module and not panel state
 * ----------------------------------------
 * The EBSD workspace is six views of one thing — the map, GROD, KAM, the
 * summary, the distributions, the discrete figures — and a scan opened in one of
 * them is open in all of them. Keeping the file in whichever panel happened to
 * receive it would mean a user who opened their data, looked at the summary, and
 * came back to the map would find the map analysing the practice dataset again,
 * silently, which is the worst answer available.
 *
 * So the scan belongs to the session rather than to a panel. There is one of it,
 * every panel reads it when it runs, and the control that opens it is the same
 * control wherever it appears.
 *
 * The file never leaves the browser except as the body of the request that
 * analyses it: there is no store on the server, nothing is kept between
 * requests, and a scan that failed to parse leaves nothing behind to clean up.
 */

import { el, formatNumber } from './dom.js';

/** Extensions that hold HDF5 rather than text, and so travel base64-encoded. */
const HDF5_SUFFIXES = ['.oh5', '.h5'];

/** The open scan: `{name, text}` or `{name, data_base64}`, or null. */
let current = null;

/** What the panel that opened it says about it, so every panel can say it too. */
let note = '';

export function openScanPayload() {
  return current;
}

export function openScanNote() {
  return note;
}

export function clearScan() {
  current = null;
  note = '';
}

/**
 * A file's bytes as base64, for the JSON field that carries it to the server.
 *
 * Encoded in chunks rather than by spreading the whole array into
 * `String.fromCharCode`: a scan is megabytes, and one argument per byte
 * overflows the call stack somewhere in the low hundreds of thousands.
 */
async function readAsBase64(file) {
  const bytes = new Uint8Array(await file.arrayBuffer());
  const chunk = 0x8000;
  const pieces = [];
  for (let start = 0; start < bytes.length; start += chunk) {
    pieces.push(String.fromCharCode.apply(null, bytes.subarray(start, start + chunk)));
  }
  return btoa(pieces.join(''));
}

/**
 * Build the "open a scan" control group for one panel's rail.
 *
 * Every EBSD panel gets one, and they are all the same control over the same
 * scan: opening a file in the distributions view is opening it for the map.
 *
 * @param {object} options
 * @param {() => void} options.onChange - Re-run the panel; called after a scan
 *   is opened and after one is closed.
 * @param {(error: unknown) => void} options.showError
 * @returns {{element: HTMLElement, setStatus: Function}}
 */
export function scanControls({ onChange, showError }) {
  const idle = 'No scan open — the practice dataset chosen below is being analysed.';
  const status = el('p.field__help', { text: note || idle });

  const input = el('input', {
    type: 'file',
    accept: '.ang,.ctf,.oh5,.h5',
    'aria-label': 'Open an EBSD scan file',
    onchange: (event) => openScan(event.target.files?.[0]),
  });

  const closeButton = el('button.button', {
    type: 'button',
    text: 'Close the scan',
    hidden: current === null,
    onclick: () => {
      clearScan();
      input.value = '';
      closeButton.hidden = true;
      status.textContent = idle;
      onChange();
    },
  });

  async function openScan(file) {
    if (!file) return;
    status.textContent = `Reading ${file.name}…`;
    try {
      current = HDF5_SUFFIXES.some((suffix) => file.name.toLowerCase().endsWith(suffix))
        ? { name: file.name, data_base64: await readAsBase64(file) }
        : { name: file.name, text: await file.text() };
      note = `${file.name} — ${formatNumber(file.size / 1024, 0)} kB open.`;
      closeButton.hidden = false;
      status.textContent = note;
      await onChange();
    } catch (error) {
      clearScan();
      closeButton.hidden = true;
      status.textContent = `${file.name} could not be read in the browser.`;
      showError(error);
    }
  }

  const element = el('details.group', { open: true }, [
    el('summary', { text: 'Open a scan' }),
    el('div.group__body', {}, [
      el('p.field__help', {
        text:
          'An EDAX/TSL .ang, an Oxford/HKL .ctf, or an EDAX OIM .oh5 or .h5 — the last two ' +
          'being one HDF5 format under two extensions. It is read by the same importer a ' +
          'script would call, so the phases, the symmetry and the quality channels come from ' +
          'its own header. While one is open it replaces the practice dataset in every view ' +
          'of this workspace, not only in this one.',
      }),
      el('p.field__help', {
        text:
          'An HDF5 scan saved with its diffraction patterns is far larger than a request can ' +
          'carry — export it without the patterns, or read it with pytex.adapters.read_scan ' +
          'in a script.',
      }),
      input,
      status,
      closeButton,
    ]),
  ]);

  return {
    element,
    /** Report something about the scan in place, such as a read failure. */
    setStatus(text) {
      status.textContent = text ?? (note || idle);
    },
  };
}

/**
 * Add the open scan to a request, if there is one.
 *
 * Every EBSD operation takes the same `scan_file` field, so this is the one
 * place that knows a request may carry a file at all.
 */
export function withScan(request) {
  return current ? { ...request, scan_file: current } : { ...request };
}
