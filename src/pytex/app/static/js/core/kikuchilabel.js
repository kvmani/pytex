/**
 * Writing a band's indices along the band.
 *
 * A Kikuchi band is a line, and the thing that identifies it is which line it
 * is. A horizontal label beside a steeply running band belongs — visually — to
 * whatever else is near it, and on a pattern where six bands cross within a
 * few tens of pixels that is a real ambiguity, not a matter of taste: the
 * reader has to guess which line the text is for. Written *along* the band the
 * text is unambiguous by construction, which is how bands are annotated on a
 * published pattern and on every hand-marked plate.
 *
 * Two rules make it readable rather than merely aligned:
 *
 * - **A line has no direction.** The angle is folded into the half-turn
 *   `(-90, 90]`, so a band running up-left is labelled the same way as the
 *   identical band running down-right, and no label is ever upside down.
 * - **The text sits beside the band, not on it.** The offset is applied along
 *   the band's normal, which the caller already knows, so the line stays
 *   visible under its own name.
 *
 * Both the simulated-plate overlay and the crystal's stereographic map use
 * this, because the annotation convention should not depend on which figure a
 * band happens to appear in.
 */

import { svg } from './dom.js';

/**
 * The angle a label should take to run along a line, in degrees.
 *
 * @param {number[]} from - `[x, y]` in the drawing's own coordinates.
 * @param {number[]} to - `[x, y]`; the other end of the run being labelled.
 * @returns {number} Degrees clockwise from the x-axis, folded into `(-90, 90]`
 *   so that reading the text never requires turning the page over.
 */
export function labelAngleDeg(from, to) {
  const dx = Number(to[0]) - Number(from[0]);
  const dy = Number(to[1]) - Number(from[1]);
  if (!Number.isFinite(dx) || !Number.isFinite(dy) || (dx === 0 && dy === 0)) return 0;
  const angle = (Math.atan2(dy, dx) * 180) / Math.PI;
  if (angle > 90) return angle - 180;
  if (angle <= -90) return angle + 180;
  return angle;
}

/**
 * A text node written along a band, haloed so it survives a busy background.
 *
 * @param {object} options
 * @param {number} options.x - Anchor, in the drawing's coordinates.
 * @param {number} options.y
 * @param {number} options.angleDeg - Normally from {@link labelAngleDeg}.
 * @param {string} options.text
 * @param {number} options.fontSize
 * @param {string} options.colour
 * @param {string} options.haloColour - Drawn under the glyphs by `paint-order`.
 * @param {number} [options.haloWidth] - Defaults to a sixth of the font size.
 * @param {string} [options.weight] - SVG `font-weight`.
 * @returns {SVGElement}
 */
export function bandLabelNode({
  x,
  y,
  angleDeg,
  text,
  fontSize,
  colour,
  haloColour,
  haloWidth = fontSize / 6,
  weight = '400',
}) {
  return svg('text', {
    x,
    y,
    // Rotating about the anchor keeps the placement decision — which is about
    // where on the band the name belongs — separate from the orientation.
    transform: `rotate(${angleDeg.toFixed(2)} ${x} ${y})`,
    'font-size': fontSize,
    'font-weight': weight,
    fill: colour,
    stroke: haloColour,
    'stroke-width': haloWidth,
    'paint-order': 'stroke',
    'text-anchor': 'middle',
    'dominant-baseline': 'middle',
    text,
  });
}
