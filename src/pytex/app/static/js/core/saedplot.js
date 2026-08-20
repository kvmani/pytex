/**
 * Drawing a simulated diffraction plate: the spots, the scale, the bands.
 *
 * Two panels draw the same picture for opposite reasons. The solver draws a
 * practice plate to be *indexed* — clicked on, measured, and compared with an
 * answer it hides until asked. The simulator draws a plate to be *read* — the
 * reference a real exposure is held against. They must look identical, because
 * they are the same calculation and a reader who learns one plate has learned
 * the other; so the drawing lives here once rather than twice.
 *
 * Everything takes plain data — the `pattern` payload a Python service returns,
 * and the `data` of a Kikuchi overlay — and returns nothing. Panel state stays
 * in the panel: what to label, what to make hoverable, and what colour scheme
 * the surrounding interface has are passed in, not reached for.
 */

import { svg } from './dom.js';

/** Colour of a simulated reflection. A plate is monochrome; so is this. */
export const SPOT_COLOUR = '#eaf2ff';

/** The plate itself: dark in every interface theme, because that is what it is. */
export const PLATE_COLOUR = '#05070d';

/*
 * Overlay colours are fixed, not theme tokens, and each is drawn over a dark
 * halo.
 *
 * The plate is always dark whatever the interface theme, so `var(--teal)` and
 * `var(--violet)` — which resolve to deep, saturated values in light mode —
 * would put a dark line on a near-black ground for half the users. And an
 * uploaded micrograph may be a light-ground print, where the reverse happens.
 * A bright stroke over a dark one is legible on both, which is the only way to
 * be right about a background this code does not control.
 */
export const KIKUCHI_COLOUR = '#8fd9ff';
/** The one band that is a route rather than a landmark. */
export const KIKUCHI_ROUTE_COLOUR = '#a3e635';
export const HALO_COLOUR = '#05070d';
/** Index labels drawn beside a simulated spot. */
export const LABEL_COLOUR = '#9fd2ff';

/**
 * Draw a simulated pattern: dark ground, glowing spots, optional indices.
 *
 * A plate is dark with bright spots in every theme, because that is what a
 * diffraction pattern is; inverting it for a light theme would make a practice
 * pattern look unlike the thing it is practice for. Each reflection is a soft
 * radial glow with a solid core, scaled by the apparent radius the simulation
 * reports, so a strong reflection reads as strong at a glance in the way it does
 * on a real exposure.
 *
 * @param {SVGElement} root - Group to draw into, in pattern pixel coordinates.
 * @param {object} pattern - The `pattern` payload of a simulation result.
 * @param {object} [options]
 * @param {boolean} [options.labels] - Write each reflection's index beside it.
 * @param {Function} [options.hoverable] - `(node, row) => void`, normally the
 *   frame's own `hoverable`, so a spot carries its full row under the pointer.
 * @param {boolean} [options.scaleBar] - Draw the reciprocal-space scale bar.
 * @param {string} [options.gradientId] - Unique within the document; two
 *   patterns drawn on one page must not share a gradient id.
 */
export function drawSimulatedPattern(root, pattern, options = {}) {
  const {
    labels = false,
    hoverable = null,
    scaleBar = true,
    gradientId = 'saed-glow',
  } = options;
  const [cx, cy] = pattern.centre_px;
  root.append(
    svg('defs', {}, [
      svg('radialGradient', { id: gradientId }, [
        svg('stop', { offset: '0%', 'stop-color': SPOT_COLOUR, 'stop-opacity': '0.95' }),
        svg('stop', { offset: '45%', 'stop-color': SPOT_COLOUR, 'stop-opacity': '0.28' }),
        svg('stop', { offset: '100%', 'stop-color': SPOT_COLOUR, 'stop-opacity': '0' }),
      ]),
    ]),
    svg('rect', {
      x: 0,
      y: 0,
      width: pattern.width_px,
      height: pattern.height_px,
      fill: PLATE_COLOUR,
    }),
  );

  // The transmitted beam, drawn unmistakably brighter and larger than any
  // reflection. This is not decoration: the solver's first instruction is
  // "click the transmitted beam", and in an earlier draft the beam was barely
  // distinguishable from a strong 200 spot, which makes that instruction
  // impossible to follow. On a real plate the direct beam is orders of
  // magnitude brighter — that is why a beam stop exists.
  //
  // Its size is tied to the nearest reflection rather than fixed, so a dense
  // hexagonal zone does not have its innermost spots swallowed by the glow.
  const nearest = pattern.spots.length
    ? Math.min(...pattern.spots.map((spot) => Math.hypot(spot.x - cx, spot.y - cy)))
    : pattern.width_px / 8;
  const beamGlow = Math.max(28, Math.min(0.5 * nearest, pattern.width_px / 8));
  root.append(
    svg('circle', { cx, cy, r: beamGlow, fill: `url(#${gradientId})`, opacity: 0.95 }),
    svg('circle', { cx, cy, r: beamGlow * 0.55, fill: `url(#${gradientId})`, opacity: 0.95 }),
    svg('circle', { cx, cy, r: beamGlow * 0.2, fill: '#ffffff', opacity: 0.95 }),
  );

  const marker = Math.max(pattern.width_px, pattern.height_px) / 140;
  for (const spot of pattern.spots) {
    const glow = svg('circle', {
      cx: spot.x,
      cy: spot.y,
      r: spot.radius_px * 3.0,
      fill: `url(#${gradientId})`,
      // Opacity carries the intensity, because apparent radius alone barely
      // separates reflections whose kinematic intensities differ by tens of
      // percent — which within one zone is the usual case.
      opacity: 0.25 + 0.7 * spot.intensity,
    });
    const core = svg('circle', {
      cx: spot.x,
      cy: spot.y,
      r: spot.radius_px * 0.62,
      fill: SPOT_COLOUR,
      opacity: 0.45 + 0.55 * spot.intensity,
    });
    root.append(glow, core);
    if (hoverable) {
      hoverable(core, {
        Index: spot.label,
        'd / Å': spot.d_angstrom,
        '|g| / Å⁻¹': spot.g_inv_angstrom,
        'Relative intensity': spot.intensity,
        x: spot.x,
        y: spot.y,
      });
    }
    if (labels) {
      root.append(
        svg('text', {
          x: spot.x + spot.radius_px * 1.6,
          y: spot.y - spot.radius_px * 1.2,
          'font-size': marker * 2.1,
          fill: LABEL_COLOUR,
          text: spot.label,
        }),
      );
    }
  }

  if (scaleBar) drawScaleBar(root, pattern, marker);
}

/**
 * A reciprocal-space scale bar, in inverse angstroms.
 *
 * A diffraction pattern without a scale is a picture. The bar makes the
 * calibration visible rather than merely stored: change the camera length and
 * the bar changes with the pattern, which is the relation the calibration is.
 */
export function drawScaleBar(root, pattern, marker) {
  const pxPerInvAngstrom = pattern.camera_constant_mm_angstrom / pattern.pixel_size_mm;
  if (!Number.isFinite(pxPerInvAngstrom) || pxPerInvAngstrom <= 0) return;
  const target = pattern.width_px * 0.22;
  const choices = [0.05, 0.1, 0.2, 0.5, 1, 2, 5];
  let best = choices[0];
  for (const value of choices) {
    if (Math.abs(value * pxPerInvAngstrom - target) < Math.abs(best * pxPerInvAngstrom - target)) {
      best = value;
    }
  }
  const length = best * pxPerInvAngstrom;
  const x = pattern.width_px * 0.06;
  const y = pattern.height_px * 0.94;
  root.append(
    svg('line', {
      x1: x,
      y1: y,
      x2: x + length,
      y2: y,
      stroke: SPOT_COLOUR,
      'stroke-width': marker * 0.5,
      'stroke-opacity': 0.75,
    }),
    svg('text', {
      x,
      y: y - marker * 1.1,
      'font-size': marker * 2.2,
      fill: SPOT_COLOUR,
      'fill-opacity': 0.75,
      text: `${best} Å⁻¹`,
    }),
  );
}

/**
 * Draw the Kikuchi bands of an overlay result over the pattern.
 *
 * Each band is its centre line — dashed, because the centre line is a
 * construction rather than something a plate shows — between its two edges,
 * which are what is actually recorded. Prominence sets the opacity, so the
 * strong bands read first, and the band that is a *route* to a named target is
 * drawn in its own colour and labelled with the instruction rather than with
 * its indices.
 *
 * @param {SVGElement} root
 * @param {object} data - The `data` of a `tem.kikuchi_overlay` result.
 * @param {object} frame - `{width, height}` of the pattern, in its own pixels.
 */
export function drawKikuchiBands(root, data, { width, height }) {
  if (!data?.bands?.length) return;
  const group = svg('g', { 'pointer-events': 'none' });
  const weight = Math.max(width, height) / 900;
  const font = Math.max(width, height) / 52;
  for (const band of data.bands) {
    const colour = band.connecting ? KIKUCHI_ROUTE_COLOUR : KIKUCHI_COLOUR;
    const prominence = 0.35 + 0.45 * Math.min(1, Number(band.intensity) || 0);
    const [[x1, y1], [x2, y2]] = band.centre;
    group.append(
      svg('line', {
        x1,
        y1,
        x2,
        y2,
        stroke: colour,
        'stroke-opacity': band.connecting ? 0.85 : prominence * 0.6,
        'stroke-width': band.connecting ? weight * 1.6 : weight * 0.9,
        'stroke-dasharray': `${weight * 6} ${weight * 6}`,
      }),
    );
    for (const edge of band.edges) {
      for (const run of edge) {
        group.append(
          svg('polyline', {
            points: run.map(([x, y]) => `${x},${y}`).join(' '),
            fill: 'none',
            stroke: colour,
            'stroke-opacity': prominence,
            'stroke-width': band.connecting ? weight * 2 : weight * 1.2,
            'stroke-dasharray': `${weight * 2} ${weight * 3}`,
          }),
        );
      }
    }
    const [labelX, labelY] = band.label_at;
    const text = band.connecting && data.connecting ? data.connecting.text : band.label;
    group.append(
      svg('text', {
        x: labelX,
        y: labelY,
        'font-size': band.connecting ? font * 1.1 : font,
        'font-weight': band.connecting ? '600' : '400',
        fill: colour,
        stroke: HALO_COLOUR,
        'stroke-width': font / 6,
        'paint-order': 'stroke',
        'text-anchor': 'middle',
        'dominant-baseline': 'middle',
        text,
      }),
    );
  }
  root.append(group);
}
