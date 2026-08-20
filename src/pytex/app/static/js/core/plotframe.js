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
 * 3. **One viewport language.** The wheel zooms about the pointer, Shift-drag,
 *    middle-drag or the pan tool moves the camera, and Fit restores the
 *    complete figure. The viewBox is the camera, so cursor coordinates remain
 *    correct after either operation. Zoom runs below 100% as well as above it,
 *    because "show me the whole figure and its surroundings" is as common a
 *    request as "show me this spot closely".
 * 4. **The figure's own controls travel with it.** Anything that changes what
 *    the plot shows — a legend that toggles a packet, a variant picker — is
 *    mounted in the frame's control strip rather than as a sibling below it, so
 *    it stays on screen with the figure instead of being pushed under the fold
 *    by the result tables.
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
 * @param {boolean} [options.readout] - Give the frame a readout bar under the
 *   drawing, and move the cursor readout into it. For a panel whose numbers are
 *   read *while* the picture is read — the TEM plate is measured with the
 *   pointer — where a card floating over the drawing covers the data it is
 *   reporting on.
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
  readout = false,
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
    // Below 1 as well as above it: Fit shows the drawing at exactly 100%, and a
    // reader who wants to see a wide pattern's tails, or to leave room around a
    // pole figure before exporting it, has nowhere to go if 100% is the floor.
    minZoom: 0.2,
    maxZoom: 24,
    // The pan tool: a sticky mode for pointers that have no middle button and
    // for anyone who should not have to know that Shift-drag is a thing.
    panTool: false,
    // Whether the pointer moved during the drag now ending, so the click it
    // synthesises can be swallowed rather than read as a pick.
    dragged: false,
  };
  const canvas = el('div.plot__canvas');
  const cursor = el('output.plot__cursor', { text: '' });
  // A panel-owned readout pinned to the top-left of the drawing: the TEM panel
  // puts the measurements taken off the pattern there, where they are read
  // against the pattern itself rather than in a table below it.
  const overlay = el('div.plot__overlay', { hidden: true });
  const detail = el('div.plot__detail', { hidden: true });
  const status = el('p.plot__status', { text: '' });
  const controls = el('div.plot__controls', { hidden: true });
  const zoomReadout = el('output.plot__zoom', { text: '100%', title: 'Current plot zoom' });

  const panButton = el('button.button.button--icon', {
    type: 'button',
    text: '✥',
    title: 'Pan tool: drag the figure with the left button',
    'aria-label': 'Pan tool',
    'aria-pressed': 'false',
    onclick: () => setPanTool(!view.panTool),
  });

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
        panButton,
        el('button.button', {
          type: 'button', text: 'Fit', title: 'Fit the complete plot',
          onclick: () => fitView(),
        }),
      ]
    : [];

  /*
   * The readout bar: the panel's numbers and the live cursor, under the drawing.
   *
   * A card over the figure is right when it annotates a picture that is looked
   * at, and wrong when it annotates a picture that is *worked on*. On the TEM
   * plate the measurements are of the spots the user is clicking, so the card
   * covered exactly the region being measured, and the corner it hid was a
   * corner of the data. Below the drawing it is never in the way, it is always
   * legible without a hover, and it costs the figure only its own height.
   */
  const readoutSlot = el('div.plot__readout-panel');
  const readoutBar = el('div.plot__readout', { hidden: !readout }, [
    readoutSlot,
    readout
      ? el('div.plot__readout-cursor', {}, [
          el('div.measure__title', { text: 'Under the pointer' }),
          cursor,
        ])
      : null,
  ]);
  if (readout) {
    cursor.classList.add('plot__cursor--bar');
    // In the bar the readout has a permanent home, so an empty one reads as
    // broken rather than as "the pointer is elsewhere". It rests on a dash.
    cursor.textContent = '—';
  }

  const element = el('figure.plot', {}, [
    el('figcaption.plot__header', {}, [
      el('h2.plot__title', { text: title }),
      el('div.plot__toolbar', {}, [...viewportToolbar, ...toolbar]),
    ]),
    el('div.plot__stage', {}, [canvas, overlay, detail, readout ? null : cursor]),
    readoutBar,
    controls,
    status,
  ]);

  function setPanTool(active) {
    view.panTool = Boolean(active);
    panButton.setAttribute('aria-pressed', String(view.panTool));
    panButton.title = view.panTool
      ? 'Pan tool on: drag the figure; click again for the cursor'
      : 'Pan tool: drag the figure with the left button';
    if (view.svg) view.svg.dataset.pan = view.panTool ? 'tool' : '';
  }

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

  /**
   * Keep the camera near the drawing.
   *
   * Stated as a bound on the camera's *centre* rather than on its edges. The
   * edge form has no answer once the zoom is below 100%: the viewport is then
   * wider than the drawing plus its margins, the lower bound crosses the upper
   * one, and the figure either snaps to a corner or refuses to pan at all. The
   * centre always has an interval to live in, at any zoom, so panning while
   * zoomed out behaves like panning while zoomed in.
   */
  function boundedBox(box) {
    const base = view.base;
    if (!base) return box;
    const axis = (start, size, baseStart, baseSize) => {
      const margin = baseSize * 0.45;
      const centre = start + size / 2;
      const bounded = Math.min(
        Math.max(centre, baseStart - margin),
        baseStart + baseSize + margin,
      );
      return bounded - size / 2;
    };
    return {
      ...box,
      x: axis(box.x, box.width, base.x, base.width),
      y: axis(box.y, box.height, base.y, base.height),
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

  /*
   * The readout sits *behind* the drawing, and lifts while the pointer is on it.
   *
   * A card pinned over the top-left corner reports on the picture — and covers
   * a corner of it, which on a diffraction pattern is a corner of the data. The
   * picture wins by default: the card is painted under the drawing, so the
   * opaque part of the figure masks it and only the letterbox margin shows it.
   * Bringing the pointer onto its rectangle raises it in full, and leaving
   * restores the unobstructed view.
   *
   * The hover cannot be the card's own `:hover`, because behind an opaque
   * drawing it never receives the pointer — so the rectangle is tested against
   * the pointer on the surface above it. Nothing about clicking changes: the
   * card has always been `pointer-events: none`, and a click in that corner has
   * always gone to the figure.
   */
  function attachOverlayReveal(svgNode) {
    const raise = (event) => {
      if (overlay.hidden) return;
      const box = overlay.getBoundingClientRect();
      if (!box.width) return;
      const inside =
        event.clientX >= box.left &&
        event.clientX <= box.right &&
        event.clientY >= box.top &&
        event.clientY <= box.bottom;
      overlay.classList.toggle('plot__overlay--raised', inside);
    };
    svgNode.addEventListener('pointermove', raise);
    svgNode.addEventListener('pointerleave', () => {
      overlay.classList.remove('plot__overlay--raised');
    });
  }

  function attachCursor(svgNode) {
    svgNode.addEventListener('pointermove', (event) => {
      if (!mapping.toData) return;
      const point = pointerToViewBox(event, svgNode);
      if (!point) return;
      const data = mapping.toData(point.x, point.y);
      if (!data) {
        cursor.textContent = readout ? '—' : '';
        return;
      }
      cursor.textContent = mapping.formatCursor
        ? mapping.formatCursor(data)
        : `${formatNumber(data.x, mapping.digits)}, ${formatNumber(data.y, mapping.digits)} ${
            mapping.units
          }`.trim();
    });
    svgNode.addEventListener('pointerleave', () => {
      cursor.textContent = readout ? '—' : '';
    });
  }

  function attachViewport(svgNode) {
    svgNode.classList.add('plot__surface');
    svgNode.setAttribute(
      'data-viewport-help',
      'Scroll to zoom; Shift-drag, middle-drag or the pan tool to pan',
    );
    svgNode.dataset.pan = view.panTool ? 'tool' : '';
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
      const panning =
        event.button === 1 || (event.button === 0 && (event.shiftKey || view.panTool));
      if (!panning) return;
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
      if (dx || dy) view.dragged = true;
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
    // A drag on a panel that picks by clicking — the TEM pattern is picked
    // entirely by clicking — ends with a click event on the same node, so
    // moving the view would silently drop a pick where the drag finished. The
    // capture phase is the only place to stop it before the panel's own
    // handler, which was attached to the same element first.
    svgNode.addEventListener(
      'click',
      (event) => {
        if (!view.panTool && !view.dragged) return;
        view.dragged = false;
        event.stopPropagation();
        event.preventDefault();
      },
      true,
    );
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

    /**
     * Where a pointer event lands, in the panel's own data coordinates.
     *
     * The camera is the `viewBox`, so this is the *only* correct conversion once
     * zoom and pan exist: a panel that measures from the element's bounding box
     * alone is right at the base view and wrong everywhere else. The TEM panel
     * had exactly that private copy, and every pick made after a pan landed off
     * by the camera offset — visibly, since the marker is drawn at the stored
     * coordinate rather than under the cursor.
     *
     * Returns null when there is no drawing, so a caller can ignore the event
     * rather than pick at NaN.
     *
     * @param {PointerEvent|MouseEvent} event
     * @returns {object|null} The result of `toData`, or viewBox coordinates if
     *   the panel supplied no mapping.
     */
    pointerToData(event) {
      const node = view.svg ?? canvas.firstElementChild;
      if (!(node instanceof SVGSVGElement)) return null;
      const point = pointerToViewBox(event, node);
      if (!point) return null;
      return mapping.toData ? mapping.toData(point.x, point.y) : point;
    },

    /** Put an SVG (or any node) on the stage and wire the instrument interactions to it. */
    setContent(node, { preserveViewport = false } = {}) {
      const oldBase = view.base;
      const oldCurrent = view.current;
      clear(canvas);
      canvas.append(node);
      if (node instanceof SVGSVGElement) {
        attachCursor(node);
        attachOverlayReveal(node);
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

    /**
     * Put a panel's own readout in the top-left of the drawing.
     *
     * Non-interactive and outside the SVG, so it neither scales with the zoom
     * nor swallows a click meant for the picture underneath it.
     *
     * @param {Node|null} node - Pass null to clear it.
     */
    setOverlay(node) {
      if (node) overlay.replaceChildren(node);
      else overlay.replaceChildren();
      overlay.hidden = overlay.childElementCount === 0;
      return overlay;
    },

    /**
     * Put the panel's own numbers in the readout bar, under the drawing.
     *
     * Only frames built with `readout: true` have one. Unlike `setOverlay`,
     * nothing here is ever covered by the figure, so it is where a measurement
     * that is read *while* picking belongs.
     *
     * @param {Node|null} node - Pass null to clear it.
     */
    setReadout(node) {
      if (node) readoutSlot.replaceChildren(node);
      else readoutSlot.replaceChildren();
      return readoutSlot;
    },

    /**
     * Mount the figure's own controls inside the frame, under the drawing.
     *
     * For a legend that toggles what is drawn, or a variant picker: these are
     * part of the instrument, and a reader who has to scroll away from the
     * figure to reach them cannot see the effect of what they just pressed.
     * The strip scrolls internally if it is long, so it can never push the
     * drawing off the screen.
     *
     * @param {...Node} nodes
     */
    setControls(...nodes) {
      controls.replaceChildren(...nodes.filter(Boolean));
      controls.hidden = controls.childElementCount === 0;
      return controls;
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
