/**
 * Powder diffraction pattern file management for the XRD workspace.
 *
 * Provides a shared pattern file loader across the three XRD analysis
 * views: Background Estimation, Lattice Parameter Determination, and
 * Rietveld Refinement.
 *
 * Supported formats:
 * - .xy / .dat / .csv / .txt: Two-column ASCII files (2θ and intensity/counts)
 * - .xrdml: PANalytical XML powder diffractograms (1D continuous/stepped scans)
 */

import { el, formatNumber } from './dom.js';
import { explainer } from './explainer.js';

/** The open pattern: `{name, text}` or `{name, data_base64}`, or null. */
let current = null;

/** Informational note about the currently loaded pattern file. */
let note = '';

export function openPatternPayload() {
  return current;
}

export function openPatternNote() {
  return note;
}

export function clearPattern() {
  current = null;
  note = '';
}

/**
 * Build the "Open an experimental pattern" control group for the XRD rail.
 *
 * @param {object} options
 * @param {() => Promise<void>|void} options.onChange - Called after a pattern is loaded or closed.
 * @param {(error: unknown) => void} options.showError - Error handler.
 * @returns {{element: HTMLElement, setStatus: (text: string) => void}}
 */
export function patternControls({ onChange, showError }) {
  const idle = 'No pattern file open — a synthetic demonstration scan is analysed.';
  const status = el('p.field__help', { text: note || idle });

  const input = el('input', {
    type: 'file',
    accept: '.xy,.xrdml,.csv,.dat,.txt',
    'aria-label': 'Open an experimental XRD pattern file',
    onchange: (event) => openPattern(event.target.files?.[0]),
  });

  const closeButton = el('button.button', {
    type: 'button',
    text: 'Close pattern (revert to demo scan)',
    hidden: current === null,
    onclick: () => {
      clearPattern();
      input.value = '';
      closeButton.hidden = true;
      status.textContent = idle;
      onChange();
    },
  });

  async function openPattern(file) {
    if (!file) return;
    status.textContent = `Reading ${file.name}…`;
    try {
      const text = await file.text();
      current = { name: file.name, text };
      note = `${file.name} — ${formatNumber(file.size / 1024, 1)} kB loaded.`;
      closeButton.hidden = false;
      status.textContent = note;
      await onChange();
    } catch (error) {
      clearPattern();
      closeButton.hidden = true;
      status.textContent = `${file.name} could not be read in the browser.`;
      showError(error);
    }
  }

  const element = el('details.group', { open: true }, [
    el('summary', { text: 'Experimental pattern' }),
    el('div.group__body', {}, [
      explainer(
        'Load an experimental diffractogram (.xy, .xrdml, .csv, .dat) for background ' +
          'estimation, lattice parameter determination, or Rietveld refinement. ' +
          'For .xrdml files, the 2θ axis, counts/intensities, and laboratory X-ray ' +
          'wavelength (e.g. Cu Kα doublet) are read directly from the instrument metadata.',
        { label: 'Supported pattern formats' },
      ),
      input,
      status,
      closeButton,
    ]),
  ]);

  return {
    element,
    setStatus(text) {
      status.textContent = text ?? (note || idle);
    },
  };
}

/**
 * Manage form fields for pattern upload and data source.
 * Hides raw object parameters (`scan_file`) and updates data source defaults.
 *
 * @param {object} form - The value returned by `buildForm`.
 */
export function adoptForm(form) {
  const scanFileField = form.field?.('scan_file');
  if (scanFileField?.element) {
    scanFileField.element.hidden = true;
  }
  for (const field of form.element.querySelectorAll('.field')) {
    if (field.querySelector('[id^="ctl-scan_file-"]')) {
      field.hidden = true;
    }
  }
}

/**
 * Inject the open pattern into a request payload if present.
 *
 * @param {object} request - Outgoing service request payload.
 * @returns {object} Updated request payload.
 */
export function withPattern(request) {
  const next = { ...request };
  if (current) {
    next.scan_file = current;
    next.data_source = 'file';
  }
  return next;
}
