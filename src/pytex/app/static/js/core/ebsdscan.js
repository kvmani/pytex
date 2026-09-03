/**
 * What the EBSD workspace is looking at, shared by all of its panels.
 *
 * There is exactly one answer to "which map am I analysing" — a file the user
 * opened, or one of the practice datasets — and it lives here rather than in
 * whichever panel happens to be on screen.
 *
 * Why this is a module and not panel state
 * ----------------------------------------
 * The EBSD workspace is six views of one thing: the map, GROD, KAM, the
 * summary, the distributions, the discrete figures. A user who chooses a
 * dataset on the map, looks at the KAM, and comes back expects all three to
 * have been about the same microstructure.
 *
 * The opened *file* was already shared for that reason. The dataset *choice*
 * was not: every panel's form carried its own copy of the generated `dataset`
 * picker, so the IPF map could be showing the bicrystal while the grain-size
 * distribution beside it counted the polycrystal's twelve grains, with nothing
 * on screen saying they disagreed. Two views of "the same" scan reporting
 * different microstructures is the worst answer available, and it was reachable
 * in two clicks.
 *
 * So both belong to the session. One control chooses, every panel reads it when
 * it runs, and the generated per-panel copies are hidden by `adoptForm` so they
 * cannot disagree with it.
 *
 * An example may still name a dataset — several are written to teach one
 * construction — and loading one moves the *workspace* onto that dataset rather
 * than only the panel that ran it. That is the same rule stated from the other
 * side: there is one answer, and whatever last set it is the answer.
 *
 * The file never leaves the browser except as the body of the request that
 * analyses it: there is no store on the server, nothing is kept between
 * requests, and a scan that failed to parse leaves nothing behind to clean up.
 */

import { el, formatNumber } from './dom.js';
import { explainer } from './explainer.js';

/** Extensions that hold HDF5 rather than text, and so travel base64-encoded. */
const HDF5_SUFFIXES = ['.oh5', '.h5'];

/** The open scan: `{name, text}` or `{name, data_base64}`, or null. */
let current = null;

/** What the panel that opened it says about it, so every panel can say it too. */
let note = '';

/**
 * The practice dataset every panel analyses while no file is open.
 *
 * Seeded from the manifest the first time a picker is built, so the default
 * lives in Python — `pytex.app.ebsd_gallery.DEFAULT_ENTRY_ID`, the equiaxed
 * polycrystal — and is not restated here where it could drift from it.
 */
let dataset = null;

/**
 * Every picker built so far, so choosing in one moves the others.
 *
 * Panels are rebuilt as the user moves between them, so this accumulates
 * pickers whose rail has long since been replaced. They are dropped on the next
 * pass rather than by an unmount hook no panel currently has: a `Set` of
 * detached `<select>` elements is a leak, and one that grows for as long as
 * somebody keeps clicking between the EBSD views.
 */
const pickers = new Set();

/** The dataset the workspace is analysing. */
export function activeDataset() {
  return dataset;
}

/**
 * Move the whole workspace onto a dataset.
 *
 * @param {string} value - A gallery entry id, as the manifest lists them.
 * @returns {boolean} Whether this changed anything.
 */
export function setActiveDataset(value) {
  if (!value || value === dataset) return false;
  dataset = value;
  for (const picker of pickers) {
    if (picker.isConnected) picker.value = value;
    else pickers.delete(picker);
  }
  return true;
}

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
 * Build the "which scan" control group for one panel's rail.
 *
 * Every EBSD panel gets one, and they are all the same control over the same
 * choice: opening a file in the distributions view is opening it for the map,
 * and picking a practice dataset here picks it everywhere.
 *
 * @param {object} options
 * @param {object} options.operation - Manifest entry, read for the `dataset`
 *   parameter's options and default. Taking them from the manifest rather than
 *   listing them here is what keeps the gallery a Python decision.
 * @param {() => void} options.onChange - Re-run the panel; called after a scan
 *   is opened, after one is closed, and after the dataset is changed.
 * @param {(error: unknown) => void} options.showError
 * @returns {{element: HTMLElement, setStatus: Function}}
 */
export function scanControls({ operation, onChange, showError }) {
  const idle = 'No scan open — the practice dataset below is being analysed.';
  const status = el('p.field__help', { text: note || idle });
  const picker = datasetPicker(operation, onChange);

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
      picker.setFileOpen(false);
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
      picker.setFileOpen(true);
      await onChange();
    } catch (error) {
      clearScan();
      closeButton.hidden = true;
      status.textContent = `${file.name} could not be read in the browser.`;
      picker.setFileOpen(false);
      showError(error);
    }
  }

  const element = el('details.group', { open: true }, [
    el('summary', { text: 'Open a scan' }),
    el('div.group__body', {}, [
      explainer(
        'An EDAX/TSL .ang, an Oxford/HKL .ctf, or an EDAX OIM .oh5 or .h5 — the last two ' +
          'being one HDF5 format under two extensions. It is read by the same importer a ' +
          'script would call, so the phases, the symmetry and the quality channels come from ' +
          'its own header. While one is open it replaces the practice dataset in every view ' +
          'of this workspace, not only in this one.',
        { label: 'Which formats, and what opening one changes' },
      ),
      explainer(
        'An HDF5 scan saved with its diffraction patterns is far larger than a request can ' +
          'carry — export it without the patterns, or read it with pytex.adapters.read_scan ' +
          'in a script.',
        { label: 'If the file is too large' },
      ),
      input,
      status,
      closeButton,
      picker.element,
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
 * The workspace's dataset picker, wired to the shared choice.
 *
 * Built from the manifest's own options, and disabled while a file is open —
 * disabled rather than hidden, because a control that vanishes leaves a user
 * wondering where the datasets went, while a greyed one with a reason beside it
 * says what happened. The `scan_file` field is what the server prefers when
 * both are present, and this control saying so is what makes that legible.
 */
function datasetPicker(operation, onChange) {
  const parameter = (operation?.parameters ?? []).find((entry) => entry.name === 'dataset');
  const options = parameter?.options ?? [];
  if (dataset === null) dataset = parameter?.default ?? options[0]?.value ?? null;

  const select = el(
    'select',
    {
      'aria-label': 'Practice dataset',
      onchange: () => {
        if (setActiveDataset(select.value)) onChange();
      },
    },
    options.map((option) =>
      el('option', { value: option.value, text: option.label, title: option.help }),
    ),
  );
  if (dataset !== null) select.value = dataset;
  // A panel opened *while* a file is already open must show the picker held,
  // not enabled and contradicting the map beside it. The handlers below only
  // fire on a change, so the initial state has to be set here.
  select.disabled = current !== null;
  pickers.add(select);

  const held = el('p.field__help', {
    text: 'Your own scan is open, so it is what every view analyses. Close it to go back to the practice datasets.',
    hidden: current === null,
  });

  const element = el('div.field.field--compact', {}, [
    el('label.field__label', { text: 'Dataset' }),
    select,
    held,
  ]);

  return {
    element,
    /** Grey the picker out while a file of the user's own is the live scan. */
    setFileOpen(open) {
      select.disabled = open;
      held.hidden = !open;
    },
  };
}

/**
 * Take a freshly built form under this module's management.
 *
 * Two of the generated controls are presented elsewhere: the scan file by the
 * opener above, and the dataset by the picker beside it. Both are hidden here
 * rather than in each panel, so a panel cannot forget one and end up with two
 * controls for one value that disagree.
 *
 * The dataset is also *read back*, but only when the caller says the values
 * were chosen. A form seeded from an example carries whatever dataset that
 * example names, and a user clicking an example is choosing its dataset — for
 * the workspace, not for one panel. A panel *opening itself* on an example is
 * not choosing anything, and adopting that would drag the workspace back onto
 * whichever dataset the first example happens to name every time the user
 * looked at that panel. Three of the four panels open themselves this way, so
 * the distinction is not hypothetical.
 *
 * @param {object} form - The value returned by `buildForm`.
 * @param {object} [options]
 * @param {boolean} [options.adoptDataset] - Whether the form's dataset should
 *   become the workspace's. False by default: silently moving every other view
 *   is not something to do unless asked.
 */
export function adoptForm(form, { adoptDataset = false } = {}) {
  for (const field of form.element.querySelectorAll('.field')) {
    if (field.querySelector('[id^="ctl-scan_file-"]')) field.hidden = true;
    if (field.querySelector('[id^="ctl-dataset-"]')) field.hidden = true;
  }
  const field = form.field('dataset');
  if (!field) return;
  const seeded = field.read();
  if (adoptDataset && seeded) setActiveDataset(seeded);
  else if (dataset !== null) field.write(dataset);
}

/**
 * Add the open scan to a request, if there is one.
 *
 * Every EBSD operation takes the same `scan_file` and `dataset` fields, so this
 * is the one place that knows a request may carry a file at all, and the one
 * place that decides which practice dataset a request names.
 *
 * The dataset is written unconditionally rather than left to the form. A panel
 * whose generated control is hidden has no other way to say which map it means,
 * and a request that omitted it would silently fall back to the server's
 * default — which is the same value today and need not be tomorrow.
 */
export function withScan(request) {
  const next = { ...request };
  if (dataset !== null) next.dataset = dataset;
  if (current) next.scan_file = current;
  return next;
}
