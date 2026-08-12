/**
 * The frame every plot in the application is mounted into.
 *
 * Two behaviours are mandatory on every plot, so they live here once rather than
 * being re-implemented per panel (see Decision 8 in the architecture record):
 *
 * 1. **A live cursor readout.** The bottom corner reports the pointer position
 *    in the plot's own physical units, updating as the pointer moves. A plot
 *    whose axes carry units but whose cursor does not is a picture rather than
 *    an instrument.
 * 2. **Hover detail on every drawn entity.** Anything the panel marks as
 *    hoverable shows its full row — for a diffraction spot that is the indices,
 *    d, |g|, relative intensity, and which phase and variant produced it. The
 *    payload is the same row the CSV export writes, so screen and file cannot
 *    disagree.
 *
 * The panel supplies a mapping from screen coordinates to data coordinates,
 * because only the panel knows its own projection; everything else is here.
 */

import { clear, el, formatNumber } from './dom.js';

/**
 * Create a plot frame.
 *
 * @param {object} options
 * @param {string} options.title
 * @param {string} [options.units] - Unit suffix for the cursor readout, e.g. `Å⁻¹`.
 * @param {number} [options.digits] - Decimals in the readout.
 * @param {Function} [options.toData] - `(x, y) => ({x, y, ...})` in data units.
 *   Receives coordinates in the SVG's own viewBox space.
 * @param {Function} [options.formatCursor] - Overrides the default readout text.
 * @param {Node[]} [options.toolbar] - Buttons shown in the frame header.
 * @returns {{element: HTMLElement, setContent: Function, setStatus: Function, hoverable: Function}}
 */
export function plotFrame({
  title,
  units = '',
  digits = 4,
  toData = null,
  formatCursor = null,
  toolbar = [],
} = {}) {
  // Held in one object so `configure` can replace them after construction: a
  // panel often cannot describe its coordinates until it has a result — the
  // diffraction panel needs the camera constant before it can report inverse
  // angstroms — and rebuilding the whole frame to learn that would throw away
  // the toolbar and the mounted content.
  const mapping = { toData, formatCursor, units, digits };
  const canvas = el('div.plot__canvas');
  const cursor = el('output.plot__cursor', { text: '' });
  const detail = el('div.plot__detail', { hidden: true });
  const status = el('p.plot__status', { text: '' });

  const element = el('figure.plot', {}, [
    el('figcaption.plot__header', {}, [
      el('h2.plot__title', { text: title }),
      el('div.plot__toolbar', {}, toolbar),
    ]),
    el('div.plot__stage', {}, [canvas, detail, cursor]),
    status,
  ]);

  function pointerToViewBox(event, svgNode) {
    const box = svgNode.viewBox.baseVal;
    const rect = svgNode.getBoundingClientRect();
    if (!rect.width || !rect.height) return null;
    // preserveAspectRatio="xMidYMid meet" letterboxes the drawing, so the
    // scale is the smaller of the two and the offsets recentre it. Getting this
    // wrong is invisible on a square plot and badly wrong on any other.
    const scale = Math.min(rect.width / box.width, rect.height / box.height);
    const offsetX = (rect.width - box.width * scale) / 2;
    const offsetY = (rect.height - box.height * scale) / 2;
    return {
      x: (event.clientX - rect.left - offsetX) / scale + box.x,
      y: (event.clientY - rect.top - offsetY) / scale + box.y,
    };
  }

  function attachCursor(svgNode) {
    svgNode.addEventListener('pointermove', (event) => {
      if (!mapping.toData) return;
      const point = pointerToViewBox(event, svgNode);
      if (!point) return;
      const data = mapping.toData(point.x, point.y);
      if (!data) {
        cursor.textContent = '';
        return;
      }
      cursor.textContent = mapping.formatCursor
        ? mapping.formatCursor(data)
        : `${formatNumber(data.x, mapping.digits)}, ${formatNumber(data.y, mapping.digits)} ${
            mapping.units
          }`.trim();
    });
    svgNode.addEventListener('pointerleave', () => {
      cursor.textContent = '';
    });
  }

  return {
    element,

    /**
     * Change how the cursor position is interpreted, after construction.
     *
     * @param {object} next - Any of `toData`, `formatCursor`, `units`, `digits`.
     */
    configure(next) {
      Object.assign(mapping, next);
    },

    /** Put an SVG (or any node) on the stage and wire the cursor readout to it. */
    setContent(node) {
      clear(canvas);
      canvas.append(node);
      if (node instanceof SVGSVGElement) attachCursor(node);
      return node;
    },

    /** One line under the plot: counts, scale, what is being shown. */
    setStatus(text) {
      status.textContent = text ?? '';
    },

    /**
     * Mark a drawn element as carrying detail.
     *
     * @param {Element} node - The drawn entity, typically an SVG shape.
     * @param {object} row - What to show. Keys are labelled from `columns`
     *   where given, otherwise from the key itself.
     * @param {object[]} [columns] - Column descriptors from the result table,
     *   so the hover card and the exported CSV agree on labels and units.
     */
    hoverable(node, row, columns = null) {
      node.classList.add('plot__hit');
      const show = () => {
        clear(detail);
        detail.append(detailCard(row, columns));
        detail.hidden = false;
      };
      const hide = () => {
        detail.hidden = true;
      };
      node.addEventListener('pointerenter', show);
      node.addEventListener('pointerleave', hide);
      node.addEventListener('focus', show);
      node.addEventListener('blur', hide);
      node.setAttribute('tabindex', '0');
      return node;
    },
  };
}

function detailCard(row, columns) {
  const entries = columns
    ? columns.map((column) => [column.label, row[column.key], column.units, column.digits])
    : Object.entries(row).map(([key, value]) => [key, value, null, undefined]);
  return el(
    'dl.detail',
    {},
    entries.flatMap(([label, value, units, digits]) => {
      if (value === null || value === undefined) return [];
      const text =
        typeof value === 'number' ? formatNumber(value, digits) : String(value);
      return [
        el('dt', { text: label }),
        el('dd', { text: units ? `${text} ${units}` : text }),
      ];
    }),
  );
}
