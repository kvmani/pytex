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
 * 3. **One viewport language.** The wheel zooms about the pointer, Shift-drag
 *    or middle-drag pans, and Fit restores the complete figure. The viewBox is
 *    the camera, so cursor coordinates remain correct after either operation.
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
 * @param {boolean} [options.viewport] - Shared 2-D zoom/pan controls. Disable
 *   for a panel such as the crystal viewer that owns a 3-D camera.
 * @returns {{element: HTMLElement, setContent: Function, setStatus: Function, hoverable: Function}}
 */
export function plotFrame({
  title,
  units = '',
  digits = 4,
  toData = null,
  formatCursor = null,
  toolbar = [],
  viewport = true,
} = {}) {
  // Held in one object so `configure` can replace them after construction: a
  // panel often cannot describe its coordinates until it has a result — the
  // diffraction panel needs the camera constant before it can report inverse
  // angstroms — and rebuilding the whole frame to learn that would throw away
  // the toolbar and the mounted content.
  const mapping = { toData, formatCursor, units, digits };
  const view = {
    svg: null,
    base: null,
    current: null,
    drag: null,
    minZoom: 1,
    maxZoom: 24,
  };
  const canvas = el('div.plot__canvas');
  const cursor = el('output.plot__cursor', { text: '' });
  const detail = el('div.plot__detail', { hidden: true });
  const status = el('p.plot__status', { text: '' });
  const zoomReadout = el('output.plot__zoom', { text: '100%', title: 'Current plot zoom' });

  const viewportToolbar = viewport
    ? [
        el('button.button.button--icon', {
          type: 'button', text: '−', title: 'Zoom out', 'aria-label': 'Zoom out',
          onclick: () => zoomBy(1 / 1.35),
        }),
        zoomReadout,
        el('button.button.button--icon', {
          type: 'button', text: '+', title: 'Zoom in', 'aria-label': 'Zoom in',
          onclick: () => zoomBy(1.35),
        }),
        el('button.button', {
          type: 'button', text: 'Fit', title: 'Fit the complete plot',
          onclick: () => fitView(),
        }),
      ]
    : [];

  const element = el('figure.plot', {}, [
    el('figcaption.plot__header', {}, [
      el('h2.plot__title', { text: title }),
      el('div.plot__toolbar', {}, [...viewportToolbar, ...toolbar]),
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

  function copyBox(box) {
    return { x: box.x, y: box.y, width: box.width, height: box.height };
  }

  function setBox(box) {
    if (!view.svg || !view.base) return;
    view.current = copyBox(box);
    view.svg.setAttribute(
      'viewBox',
      `${view.current.x} ${view.current.y} ${view.current.width} ${view.current.height}`,
    );
    const zoom = view.base.width / view.current.width;
    zoomReadout.textContent = `${Math.round(zoom * 100)}%`;
  }

  function boundedBox(box) {
    const base = view.base;
    if (!base) return box;
    const marginX = base.width * 0.55;
    const marginY = base.height * 0.55;
    return {
      ...box,
      x: Math.min(Math.max(box.x, base.x - marginX), base.x + base.width + marginX - box.width),
      y: Math.min(Math.max(box.y, base.y - marginY), base.y + base.height + marginY - box.height),
    };
  }

  function zoomBy(multiplier, anchor = null) {
    if (!view.svg || !view.base || !view.current) return;
    const currentZoom = view.base.width / view.current.width;
    const nextZoom = Math.min(
      Math.max(currentZoom * multiplier, view.minZoom),
      view.maxZoom,
    );
    const ratio = currentZoom / nextZoom;
    const point = anchor ?? {
      x: view.current.x + view.current.width / 2,
      y: view.current.y + view.current.height / 2,
    };
    setBox(boundedBox({
      x: point.x - (point.x - view.current.x) * ratio,
      y: point.y - (point.y - view.current.y) * ratio,
      width: view.current.width * ratio,
      height: view.current.height * ratio,
    }));
  }

  function fitView() {
    if (!view.base) return;
    setBox(view.base);
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

  function attachViewport(svgNode) {
    svgNode.classList.add('plot__surface');
    svgNode.setAttribute('data-viewport-help', 'Scroll to zoom; Shift-drag or middle-drag to pan');
    svgNode.addEventListener(
      'wheel',
      (event) => {
        const point = pointerToViewBox(event, svgNode);
        if (!point) return;
        event.preventDefault();
        zoomBy(Math.exp(-event.deltaY * 0.0015), point);
      },
      { passive: false },
    );
    svgNode.addEventListener('pointerdown', (event) => {
      if (!(event.button === 1 || (event.button === 0 && event.shiftKey))) return;
      const point = pointerToViewBox(event, svgNode);
      if (!point) return;
      view.drag = { pointerId: event.pointerId, point };
      svgNode.dataset.panning = 'true';
      svgNode.setPointerCapture(event.pointerId);
      event.preventDefault();
    });
    svgNode.addEventListener('pointermove', (event) => {
      if (!view.drag || view.drag.pointerId !== event.pointerId || !view.current) return;
      const point = pointerToViewBox(event, svgNode);
      if (!point) return;
      const dx = view.drag.point.x - point.x;
      const dy = view.drag.point.y - point.y;
      setBox(boundedBox({ ...view.current, x: view.current.x + dx, y: view.current.y + dy }));
      view.drag.point = pointerToViewBox(event, svgNode) ?? point;
    });
    const endPan = (event) => {
      if (!view.drag || view.drag.pointerId !== event.pointerId) return;
      view.drag = null;
      delete svgNode.dataset.panning;
      if (svgNode.hasPointerCapture?.(event.pointerId)) svgNode.releasePointerCapture(event.pointerId);
    };
    svgNode.addEventListener('pointerup', endPan);
    svgNode.addEventListener('pointercancel', endPan);
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

    /** Put an SVG (or any node) on the stage and wire the instrument interactions to it. */
    setContent(node, { preserveViewport = false } = {}) {
      const oldBase = view.base;
      const oldCurrent = view.current;
      clear(canvas);
      canvas.append(node);
      if (node instanceof SVGSVGElement) {
        attachCursor(node);
        if (viewport) {
          const box = node.viewBox.baseVal;
          view.svg = node;
          view.base = copyBox(box);
          if (preserveViewport && oldBase && oldCurrent) {
            const zoom = oldBase.width / oldCurrent.width;
            const centerX = (oldCurrent.x + oldCurrent.width / 2 - oldBase.x) / oldBase.width;
            const centerY = (oldCurrent.y + oldCurrent.height / 2 - oldBase.y) / oldBase.height;
            const width = view.base.width / Math.min(Math.max(zoom, view.minZoom), view.maxZoom);
            const height = view.base.height * (width / view.base.width);
            setBox(boundedBox({
              x: view.base.x + centerX * view.base.width - width / 2,
              y: view.base.y + centerY * view.base.height - height / 2,
              width,
              height,
            }));
          } else {
            fitView();
          }
          attachViewport(node);
        }
      } else {
        view.svg = null;
        view.base = null;
        view.current = null;
      }
      return node;
    },

    /** Restore the complete drawing; useful to bespoke panel toolbars too. */
    fitView,

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
