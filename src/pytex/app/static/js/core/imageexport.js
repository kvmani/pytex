/**
 * Getting a picture out of the page: SVG, PNG and the clipboard, for every figure.
 *
 * Every figure in the application is SVG — a plot drawn by a panel, a figure
 * drawn by the server — and some carry a raster inside them, such as a
 * simulated micrograph at one pixel per simulated pixel. Three ways out are
 * offered for all of them, written once here:
 *
 * - **SVG**: the complete drawing, with the page's styles written onto each
 *   element so the file looks the same outside the page as inside it.
 * - **PNG**: the SVG rasterized at print resolution. Where the drawing embeds a
 *   raster, the scale is raised until that raster is reproduced at least one
 *   output pixel per source pixel, so an export is never coarser than the data.
 * - **Copy**: the PNG on the clipboard. Browsers grant this only to a secure
 *   origin (HTTPS or localhost) and some embedded web views never grant it, so
 *   a refusal is reported and the download is offered instead: copying is a
 *   convenience, never the only way to get the picture.
 *
 * A native raster — the simulation output itself — is saved by its bytes,
 * never redrawn, so its pixel count is exactly the simulation's.
 *
 * Everything happens in the page: no network request, no external library, so
 * it works on an intranet host with no internet access.
 */

import { saveBlob } from './result.js';

/** Resolution a vector figure is rasterized at, in dots per inch of its drawn size. */
export const PRINT_DPI = 300;

/** No exported raster side is larger than this; browsers cap canvas size near 16k. */
const MAX_SIDE_PX = 8192;

/**
 * The long side of a PNG exported from an on-screen plot is at least this.
 *
 * A plot drawn in a narrow column is a few hundred CSS pixels wide, and 300 dpi
 * of *that* is a thumbnail. A server-drawn figure has a designed physical size
 * and is exported at exactly 300 dpi of it instead (it passes its own scale).
 */
export const MIN_LONG_SIDE_PX = 2400;

/** CSS properties that decide how an SVG element looks, copied onto the export. */
const STYLE_PROPERTIES = [
  'fill', 'fill-opacity', 'fill-rule', 'stroke', 'stroke-width', 'stroke-opacity',
  'stroke-dasharray', 'stroke-dashoffset', 'stroke-linecap', 'stroke-linejoin',
  'opacity', 'visibility', 'display', 'font-family', 'font-size', 'font-weight',
  'font-style', 'text-anchor', 'dominant-baseline', 'letter-spacing', 'paint-order',
  'image-rendering', 'shape-rendering', 'vector-effect', 'mix-blend-mode',
];

const SVG_NS = 'http://www.w3.org/2000/svg';
const XLINK_NS = 'http://www.w3.org/1999/xlink';

/**
 * Serialize a live SVG node as a standalone document.
 *
 * The page styles SVG through classes and CSS variables that do not exist
 * outside the page, so the computed value of every visual property is written
 * onto each element. The camera is reset to the complete figure, and the
 * background the figure is drawn on is painted behind it, so a light-on-dark
 * theme does not export as invisible strokes on a transparent sheet.
 *
 * @param {SVGSVGElement} node
 * @param {object} [options]
 * @param {{x:number,y:number,width:number,height:number}} [options.viewBox] -
 *   The region to export; defaults to the node's current viewBox.
 * @param {string} [options.background] - CSS colour painted behind the drawing.
 * @returns {{markup: string, width: number, height: number}} Markup and the
 *   size, in CSS pixels, the figure occupies on screen.
 */
export function serializeSvg(node, { viewBox = null, background = null } = {}) {
  const { clone, width, height, box } = styledClone(node, viewBox);
  clone.setAttribute('xmlns', SVG_NS);
  clone.setAttribute('xmlns:xlink', XLINK_NS);
  const paint = background ?? backgroundOf(node);
  if (paint && box) {
    const rectNode = document.createElementNS(SVG_NS, 'rect');
    rectNode.setAttribute('x', String(box.x));
    rectNode.setAttribute('y', String(box.y));
    rectNode.setAttribute('width', String(box.width));
    rectNode.setAttribute('height', String(box.height));
    rectNode.setAttribute('fill', paint);
    clone.insertBefore(rectNode, clone.firstChild);
  }
  const markup = new XMLSerializer().serializeToString(clone);
  return { markup, width, height };
}

/** A clone of an SVG with its styles inlined and its camera set, and its drawn size. */
function styledClone(node, viewBox = null) {
  const clone = node.cloneNode(true);
  inlineStyles(node, clone);
  const box = viewBox ?? node.viewBox?.baseVal ?? null;
  const rect = node.getBoundingClientRect();
  let width = rect.width || 800;
  let height = rect.height || 600;
  if (box && box.width > 0 && box.height > 0) {
    clone.setAttribute('viewBox', `${box.x} ${box.y} ${box.width} ${box.height}`);
    // The on-screen box is letterboxed; the export has the drawing's own aspect.
    const aspect = box.width / box.height;
    if (width / height > aspect) width = height * aspect;
    else height = width / aspect;
  }
  clone.setAttribute('width', String(Math.round(width)));
  clone.setAttribute('height', String(Math.round(height)));
  clone.removeAttribute('class');
  return { clone, width, height, box };
}

/**
 * Serialize whatever a plot stage holds as one standalone SVG.
 *
 * A stage is usually one SVG, but several panels lay out a small multiple —
 * EBSD maps beside their key, pole figures side by side, a CBED disc beside
 * its rocking curve — as a block of separate SVGs and canvases. Exporting only
 * the first would silently drop the rest, so every drawing in the block is
 * placed in one document at the position it occupies on screen. HTML text
 * between them (a legend written as HTML) is not drawn: an SVG containing
 * HTML cannot be rasterized by a browser without tainting the canvas.
 *
 * @param {Element} root - The stage's content node.
 * @param {object} [options]
 * @param {object} [options.viewBox] - Camera for a single-SVG stage.
 * @returns {{markup: string, width: number, height: number}|null}
 */
export function serializeFigure(root, { viewBox = null } = {}) {
  if (!root) return null;
  if (root instanceof SVGSVGElement) return serializeSvg(root, { viewBox });
  const drawings = [...root.querySelectorAll('svg, canvas, img')].filter(
    (node) => !node.parentElement?.closest('svg') && isVisible(node),
  );
  if (!drawings.length) return null;
  if (drawings.length === 1 && drawings[0] instanceof SVGSVGElement) {
    return serializeSvg(drawings[0], { viewBox });
  }
  const frame = root.getBoundingClientRect();
  const document_ = document.createElementNS(SVG_NS, 'svg');
  document_.setAttribute('xmlns', SVG_NS);
  document_.setAttribute('xmlns:xlink', XLINK_NS);
  let right = 0;
  let bottom = 0;
  const placed = [];
  for (const node of drawings) {
    const rect = node.getBoundingClientRect();
    const x = rect.left - frame.left;
    const y = rect.top - frame.top;
    if (node instanceof SVGSVGElement) {
      const { clone, width, height } = styledClone(node);
      // Re-centre the letterboxed drawing inside the box it had on screen.
      clone.setAttribute('x', String(x + (rect.width - width) / 2));
      clone.setAttribute('y', String(y + (rect.height - height) / 2));
      placed.push(clone);
    } else {
      const href = node instanceof HTMLCanvasElement ? safeCanvasUrl(node) : node.src;
      if (!href || !href.startsWith('data:')) continue;
      const image = document.createElementNS(SVG_NS, 'image');
      image.setAttribute('href', href);
      image.setAttribute('x', String(x));
      image.setAttribute('y', String(y));
      image.setAttribute('width', String(rect.width));
      image.setAttribute('height', String(rect.height));
      image.setAttribute('preserveAspectRatio', 'none');
      placed.push(image);
    }
    right = Math.max(right, x + rect.width);
    bottom = Math.max(bottom, y + rect.height);
  }
  const width = Math.max(1, Math.round(right));
  const height = Math.max(1, Math.round(bottom));
  document_.setAttribute('viewBox', `0 0 ${width} ${height}`);
  document_.setAttribute('width', String(width));
  document_.setAttribute('height', String(height));
  const paint = document.createElementNS(SVG_NS, 'rect');
  paint.setAttribute('width', String(width));
  paint.setAttribute('height', String(height));
  paint.setAttribute('fill', backgroundOf(root));
  document_.append(paint, ...placed);
  return { markup: new XMLSerializer().serializeToString(document_), width, height };
}

function isVisible(node) {
  const rect = node.getBoundingClientRect();
  return rect.width > 0 && rect.height > 0;
}

function safeCanvasUrl(canvas) {
  try {
    return canvas.toDataURL('image/png');
  } catch {
    return null;
  }
}

function inlineStyles(source, target) {
  if (source.nodeType !== 1) return;
  const computed = window.getComputedStyle(source);
  const declarations = [];
  for (const property of STYLE_PROPERTIES) {
    const value = computed.getPropertyValue(property);
    if (value) declarations.push(`${property}:${value}`);
  }
  if (declarations.length) target.setAttribute('style', declarations.join(';'));
  // Interactive state is page behaviour, not figure content.
  target.removeAttribute('tabindex');
  const sourceChildren = source.children;
  const targetChildren = target.children;
  for (let index = 0; index < sourceChildren.length; index += 1) {
    inlineStyles(sourceChildren[index], targetChildren[index]);
  }
}

/** The first opaque background behind a node, so the export is drawn on what the reader saw. */
function backgroundOf(node) {
  let current = node;
  while (current && current.nodeType === 1) {
    const colour = window.getComputedStyle(current).backgroundColor;
    if (colour && colour !== 'transparent' && !/rgba\(.*,\s*0\)$/.test(colour)) return colour;
    current = current.parentElement;
  }
  return '#ffffff';
}

/**
 * Load SVG markup as an image element.
 *
 * Every embedded raster must be a data URL for this to work: an SVG loaded as
 * an image may not fetch anything, and a referenced file silently draws as
 * nothing. The panels embed their rasters as data URLs for exactly this reason.
 */
function loadSvgImage(markup) {
  return new Promise((resolve, reject) => {
    const image = new Image();
    // A data URL rather than an object URL: object URLs are reserved for
    // saving files, which goes through `saveBlob` alone.
    image.onload = () => resolve(image);
    image.onerror = () => reject(new Error('The figure could not be rasterized.'));
    image.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(markup)}`;
  });
}

/**
 * The scale at which every raster embedded in the markup keeps its own pixels.
 *
 * @param {string} markup
 * @param {number} widthPx - The drawn width the scale multiplies.
 */
async function nativeScaleOf(markup, widthPx) {
  const parsed = new DOMParser().parseFromString(markup, 'image/svg+xml');
  const root = parsed.documentElement;
  const box = root.viewBox?.baseVal;
  const images = [...parsed.getElementsByTagNameNS(SVG_NS, 'image')];
  if (!images.length || !box || !(box.width > 0)) return 1;
  let scale = 1;
  for (const node of images) {
    const href = node.getAttribute('href') ?? node.getAttributeNS(XLINK_NS, 'href');
    const drawnWidth = Number(node.getAttribute('width'));
    if (!href || !(drawnWidth > 0)) continue;
    const size = await naturalSize(href).catch(() => null);
    if (!size) continue;
    // Output pixels per viewBox unit needed so this raster is at least 1:1.
    const unitsToPx = size.width / drawnWidth;
    scale = Math.max(scale, (unitsToPx * box.width) / widthPx);
  }
  return scale;
}

function naturalSize(href) {
  return new Promise((resolve, reject) => {
    const image = new Image();
    image.onload = () => resolve({ width: image.naturalWidth, height: image.naturalHeight });
    image.onerror = reject;
    image.src = href;
  });
}

/**
 * Rasterize SVG markup to a PNG blob.
 *
 * @param {string} markup - A standalone SVG document.
 * @param {object} options
 * @param {number} options.width - Drawn width, in CSS pixels.
 * @param {number} options.height - Drawn height, in CSS pixels.
 * @param {number} [options.scale] - Output pixels per CSS pixel. Defaults to
 *   print resolution for a drawing of this size, raised until the long side is
 *   at least `MIN_LONG_SIDE_PX` and any embedded raster is at native resolution.
 * @param {string} [options.background] - Painted under the drawing; PNG has no
 *   notion of a page, and a transparent plot pasted into a dark document
 *   disappears.
 * @returns {Promise<{blob: Blob, width: number, height: number}>}
 */
export async function svgToPng(markup, { width, height, scale = null, background = '#ffffff' }) {
  const image = await loadSvgImage(markup);
  let factor = scale ?? Math.max(
    PRINT_DPI / 96,
    MIN_LONG_SIDE_PX / Math.max(width, height, 1),
    await nativeScaleOf(markup, width),
  );
  factor = Math.min(factor, MAX_SIDE_PX / Math.max(width, height));
  const canvas = document.createElement('canvas');
  canvas.width = Math.max(1, Math.round(width * factor));
  canvas.height = Math.max(1, Math.round(height * factor));
  const context = canvas.getContext('2d');
  if (background) {
    context.fillStyle = background;
    context.fillRect(0, 0, canvas.width, canvas.height);
  }
  context.drawImage(image, 0, 0, canvas.width, canvas.height);
  const blob = await new Promise((resolve, reject) => {
    canvas.toBlob((result) => (result ? resolve(result) : reject(new Error(
      'The browser refused to encode the PNG.',
    ))), 'image/png');
  });
  return { blob, width: canvas.width, height: canvas.height };
}

/**
 * Decode a data URL (or bare base64) to a blob without redrawing it.
 *
 * @param {string} data - `data:image/png;base64,...` or raw base64.
 * @param {string} [mime]
 */
export function dataUrlToBlob(data, mime = 'image/png') {
  const match = /^data:([^;,]+)(;base64)?,(.*)$/s.exec(data);
  const type = match ? match[1] : mime;
  const payload = match ? match[3] : data;
  const isBase64 = match ? Boolean(match[2]) : true;
  if (!isBase64) return new Blob([decodeURIComponent(payload)], { type });
  const binary = atob(payload);
  const bytes = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) bytes[index] = binary.charCodeAt(index);
  return new Blob([bytes], { type });
}

/**
 * Put a PNG on the clipboard, or say plainly why it cannot be.
 *
 * @param {Blob|Promise<Blob>} blob
 * @returns {Promise<boolean>} Whether the clipboard accepted it.
 */
export async function copyPng(blob) {
  const clipboard = navigator.clipboard;
  if (!clipboard?.write || typeof window.ClipboardItem !== 'function') {
    announce(
      'This window cannot put images on the clipboard (browsers allow it only over HTTPS or '
        + 'on localhost). Use Download PNG instead.',
      'warning',
    );
    return false;
  }
  try {
    // A promise, not a blob: Safari requires the item before any await.
    await clipboard.write([new window.ClipboardItem({ 'image/png': Promise.resolve(blob) })]);
    announce('Figure copied to the clipboard as PNG.', 'success');
    return true;
  } catch (error) {
    announce(
      `The clipboard refused the image (${error?.message ?? 'permission denied'}). `
        + 'Use Download PNG instead.',
      'warning',
    );
    return false;
  }
}

function announce(message, level) {
  document.dispatchEvent(new CustomEvent('pytex:saved', { detail: { message, path: null, level } }));
}

/** A file-name-safe stem: lower case, hyphens, no punctuation. */
export function fileStem(text, fallback = 'pytex-figure') {
  const cleaned = String(text ?? '')
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '')
    .slice(0, 60);
  return cleaned || fallback;
}

/**
 * Save a figure given as markup, in the requested form.
 *
 * @param {object} figure
 * @param {string} figure.markup - Standalone SVG.
 * @param {number} figure.width - Drawn width in CSS pixels.
 * @param {number} figure.height - Drawn height in CSS pixels.
 * @param {string} figure.stem - File stem.
 * @param {number} [figure.scale] - Output pixels per CSS pixel for PNG.
 * @param {'svg'|'png'|'copy'} format
 */
export async function exportMarkup(figure, format) {
  if (format === 'svg') {
    return saveBlob(`${figure.stem}.svg`, new Blob([figure.markup], { type: 'image/svg+xml' }));
  }
  if (format === 'copy') {
    // Started before any await so the clipboard still sees the user's gesture.
    const pending = svgToPng(figure.markup, figure).then((result) => result.blob);
    return copyPng(pending);
  }
  const { blob } = await svgToPng(figure.markup, figure);
  return saveBlob(`${figure.stem}.png`, blob);
}

/**
 * A compact Save menu for a figure: PNG, SVG, Copy, and any native rasters.
 *
 * @param {object} options
 * @param {Function} options.figure - `() => ({markup, width, height, stem})`,
 *   evaluated when a choice is made, so the menu always saves the current
 *   drawing; returns null when there is nothing to save.
 * @param {Function} [options.rasters] - `() => [{label, data, filename, width, height}]`,
 *   native-resolution images saved byte for byte.
 * @param {string} [options.label]
 * @returns {HTMLElement}
 */
export function exportMenu({ figure, rasters = () => [], label = 'Save' }) {
  const menu = document.createElement('details');
  menu.className = 'export-menu';
  const summary = document.createElement('summary');
  summary.className = 'button';
  summary.textContent = label;
  summary.title = 'Download this figure as PNG or SVG, or copy it';
  const list = document.createElement('div');
  list.className = 'export-menu__list';
  list.setAttribute('role', 'menu');
  menu.append(summary, list);

  const item = (text, title, action) => {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'export-menu__item';
    button.setAttribute('role', 'menuitem');
    button.textContent = text;
    button.title = title;
    button.addEventListener('click', async () => {
      menu.open = false;
      button.disabled = true;
      try {
        await action();
      } catch (error) {
        announce(`Export failed: ${error?.message ?? error}`, 'error');
      } finally {
        button.disabled = false;
      }
    });
    return button;
  };

  const withFigure = (format) => async () => {
    const current = figure();
    if (!current) {
      announce('There is no figure to save yet.', 'warning');
      return;
    }
    await exportMarkup(current, format);
  };

  menu.addEventListener('toggle', () => {
    if (!menu.open) return;
    list.replaceChildren(
      item('Download PNG', `Raster at ${PRINT_DPI} dpi, or at the native resolution of any `
        + 'image it contains, whichever is finer', withFigure('png')),
      item('Download SVG', 'Vector drawing, editable in Inkscape or Illustrator', withFigure('svg')),
      item('Copy image', 'Put the PNG on the clipboard', withFigure('copy')),
      ...rasters().map((raster) =>
        item(
          raster.label,
          `${raster.width} × ${raster.height} px, exactly as computed; no resampling`,
          () => saveBlob(raster.filename, dataUrlToBlob(raster.data)),
        )),
    );
  });
  // A click anywhere else closes the menu, as a menu should.
  document.addEventListener('click', (event) => {
    if (menu.open && !menu.contains(event.target)) menu.open = false;
  });
  return menu;
}
