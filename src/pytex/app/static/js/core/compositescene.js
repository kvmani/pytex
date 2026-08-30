/**
 * Composing several crystal scenes, and the OR overlays, into one drawable scene.
 *
 * The crystal viewer's renderer draws **one** scene payload and depth-sorts
 * everything in it into a single list. That is exactly what a two-crystal
 * figure needs — atoms of the parent must occlude atoms of the child and the
 * other way about — so this module does not add a renderer. It concatenates
 * payloads into the shape that renderer already takes, which is the whole of
 * the composition layer.
 *
 * Two placements arrive from the service and are handled differently on
 * purpose:
 *
 * - `variants.composite_scene` sends both crystals **already placed** in one
 *   world frame, so composing them is concatenation and nothing else.
 * - `variants.contact_sheet` sends each structure **once, in its own crystal
 *   frame**, plus a 3x3 placement matrix per variant, because twenty-four
 *   placed copies of both crystals would be tens of megabytes to say what a
 *   matrix multiply says exactly. `placeScene` is that multiply.
 *
 * Nothing here decides any crystallography. Every matrix was computed in
 * Python; applying one is the same arithmetic the camera already does.
 */

import { applyTranspose } from './rotation3.js';

/** `M v + t`, for a row-major 3x3 given as three rows. */
function mapPoint(matrix, translation, point) {
  return [
    matrix[0][0] * point[0] + matrix[0][1] * point[1] + matrix[0][2] * point[2] + translation[0],
    matrix[1][0] * point[0] + matrix[1][1] * point[1] + matrix[1][2] * point[2] + translation[1],
    matrix[2][0] * point[0] + matrix[2][1] * point[1] + matrix[2][2] * point[2] + translation[2],
  ];
}

/** `M v` — a direction, which does not pick up the translation. */
function mapVector(matrix, vector) {
  return mapPoint(matrix, [0, 0, 0], vector);
}

/**
 * A scene payload with a rigid placement baked into every coordinate.
 *
 * The placement is a rotation, so plane normals map the same way points do; the
 * inverse-transpose that a general affine map would need is the matrix itself
 * here, and writing it out would suggest a generality this does not have.
 */
export function placeScene(scene, matrix, translation = [0, 0, 0]) {
  if (!scene) return scene;
  const point = (value) => mapPoint(matrix, translation, value);
  return {
    ...scene,
    atoms: (scene.atoms ?? []).map((atom) => ({ ...atom, position: point(atom.position) })),
    bonds: (scene.bonds ?? []).map((bond) => ({
      ...bond,
      start: point(bond.start),
      end: point(bond.end),
    })),
    cell_edges: (scene.cell_edges ?? []).map((edge) => edge.map(point)),
    planes: (scene.planes ?? []).map((plane) => ({
      ...plane,
      vertices: plane.vertices.map(point),
      normal: mapVector(matrix, plane.normal),
    })),
    directions: (scene.directions ?? []).map((direction) => ({
      ...direction,
      start: point(direction.start),
      end: point(direction.end),
    })),
  };
}

/**
 * Recolour every atom and bond of a scene, so the two crystals can be told apart.
 *
 * This is not decoration. The phases of an orientation relationship are usually
 * the *same element* — austenite and ferrite are both iron, beta and alpha
 * zirconium are both zirconium — so colouring by species paints both crystals
 * one colour and the composite becomes a single orange blob in which no
 * parallelism can be read. In these views colour therefore carries the **phase**
 * rather than the element, and the legend says so, because a viewer that
 * silently changed what a colour meant would be worse than one that could not
 * tell the crystals apart.
 */
export function tintScene(scene, color, { opacity = 1, scale = 1, cellOpacity = 1 } = {}) {
  if (!scene) return scene;
  return {
    ...scene,
    atoms: (scene.atoms ?? []).map((atom) => ({ ...atom, color, opacity, scale })),
    bonds: (scene.bonds ?? []).map((bond) => ({
      ...bond,
      color,
      start_color: color,
      end_color: color,
      opacity,
    })),
    // The cell frame is tinted too, and survives the ghosting at full strength:
    // in the lattice-only detail level it is the only thing left of the
    // crystal, and a ghost with no frame is an absence rather than a crystal.
    cell_edges: (scene.cell_edges ?? []).map((edge) => {
      const [start, end] = Array.isArray(edge) ? edge : [edge.start, edge.end];
      return { start, end, color, opacity: cellOpacity };
    }),
  };
}

/**
 * The same scene with its atoms and bonds dropped: the cell frame alone.
 *
 * Twelve panels of two ball-and-stick crystals is a wall of coloured dots in
 * which no plane can be seen, and the question a variant wall is asked — where
 * does this crystal point — is not answered by atoms. Lattice-only is
 * therefore the default at wall size, and atoms are one click away.
 */
export function frameOnly(scene) {
  if (!scene) return scene;
  return { ...scene, atoms: [], bonds: [] };
}


/**
 * The world-frame OR overlays in the shape the renderer already draws.
 *
 * An arrow is a direction and a patch is a plane; the service emits the same
 * keys, so this is a rename rather than a conversion. Keeping the names apart
 * on the wire is deliberate — a scene's planes belong to a crystal, while these
 * belong to the relationship between two of them.
 */
export function primitiveOverlays(
  primitives,
  { offset = null, planeColor = null, directionColor = null, alpha = null, label = true } = {},
) {
  const shift = offset
    ? (point) => [point[0] + offset[0], point[1] + offset[1], point[2] + offset[2]]
    : (point) => point;
  return {
    planes: (primitives?.patches ?? []).map((patch) => ({
      vertices: patch.vertices.map(shift),
      normal: patch.normal,
      color: planeColor ?? patch.color,
      alpha: alpha ?? patch.alpha,
      label: label ? patch.label : null,
    })),
    directions: (primitives?.arrows ?? []).map((arrow) => ({
      start: shift(arrow.tail),
      end: shift(arrow.head),
      color: directionColor ?? arrow.color,
      label: label ? arrow.label : null,
    })),
  };
}

/** Concatenate two overlay bundles. */
function joinOverlays(...bundles) {
  return {
    planes: bundles.flatMap((bundle) => bundle?.planes ?? []),
    directions: bundles.flatMap((bundle) => bundle?.directions ?? []),
  };
}

/**
 * Concatenate scenes and overlays into one payload for the crystal renderer.
 *
 * `extent` is the world centre and radius the service computed over everything
 * that can be shown — for a contact sheet, over *every* variant's placement, so
 * the camera does not reframe as the reader steps through the panels and
 * mistake a change of framing for a change of orientation.
 */
export function mergeScenes(
  scenes,
  { overlays = null, extent = null, axes = null, triads = null } = {},
) {
  const parts = scenes.filter(Boolean);
  const merged = {
    atoms: parts.flatMap((scene) => scene.atoms ?? []),
    bonds: parts.flatMap((scene) => scene.bonds ?? []),
    cell_edges: parts.flatMap((scene) => scene.cell_edges ?? []),
    planes: [...parts.flatMap((scene) => scene.planes ?? []), ...(overlays?.planes ?? [])],
    directions: [
      ...parts.flatMap((scene) => scene.directions ?? []),
      ...(overlays?.directions ?? []),
    ],
    axes: axes ?? parts[0]?.axes ?? [],
    triads: triads ?? null,
    centre: extent?.centre ?? parts[0]?.centre ?? [0, 0, 0],
    radius: extent?.radius ?? parts[0]?.radius ?? 1,
    bounds: extent?.bounds ?? parts[0]?.bounds ?? null,
  };
  return merged;
}

/**
 * How each crystal is drawn in a two-phase figure, at one of two weights.
 *
 * The parent is the reference and the child is the subject, so they are not
 * drawn alike: ghosting the parent is what lets the shared plane be seen
 * through it. Ghosting is a *style*, not a claim — both crystals are the same
 * structures at the same scale — so the numbers under the panel are unaffected
 * by it and the legend says which is which.
 */
const GHOST = Object.freeze({ opacity: 0.3, scale: 0.62, cellOpacity: 0.85 });
const SOLID = Object.freeze({ opacity: 1, scale: 1, cellOpacity: 1 });

/**
 * Both triads of a two-phase scene, each in its phase's colour.
 *
 * `frames` is what the service computed: the parent axes in the world frame,
 * and the child axes *after* its placement, so the triad reports where the
 * child's own a1, a2 and c now point rather than where they started.
 */
function frameTriads(frames, { parentColor, childColor, parentLabel, childLabel }) {
  if (!frames) return null;
  return [
    { axes: frames.parent, color: parentColor, label: parentLabel, anchor: 'left' },
    { axes: frames.child, color: childColor, label: childLabel, anchor: 'right' },
  ];
}

/**
 * One panel of the variant wall: the parent and one variant's child, in one frame.
 *
 * The overlays are drawn **twice** when the crystals stand apart — once at the
 * parent and once translated onto the child. That is exact rather than
 * decorative: the two objects are parallel by construction, and a parallel
 * plane is unchanged by a translation, so the second copy is the same plane
 * drawn where the second crystal is. Drawing it once would put the whole
 * statement on one of the two crystals and leave the other unmarked, which is
 * precisely the figure that cannot be read.
 */
export function variantPanelScene(
  data,
  entry,
  {
    parentColor,
    childColor,
    parentLabel = null,
    childLabel = null,
    ghostParent = true,
    showAtoms = true,
    extent = null,
    planeColor = null,
    directionColor = null,
    planeAlpha = null,
  } = {},
) {
  const dress = (scene) => (showAtoms ? scene : frameOnly(scene));
  const parent = tintScene(
    dress(data.parent?.scene),
    parentColor,
    ghostParent ? GHOST : SOLID,
  );
  const child = tintScene(
    dress(placeScene(data.child?.scene, entry.child_matrix, entry.translation)),
    childColor,
    SOLID,
  );
  const translation = entry.translation ?? [0, 0, 0];
  const apart = translation.some((value) => Math.abs(value) > 1e-9);
  const style = { planeColor, directionColor, alpha: planeAlpha };
  const overlays = joinOverlays(
    primitiveOverlays(entry.primitives, style),
    apart
      ? primitiveOverlays(entry.primitives, { ...style, offset: translation, label: false })
      : null,
  );
  return mergeScenes([parent, child], {
    overlays,
    extent: extent ?? data.world,
    axes: data.parent?.scene?.axes,
    triads: frameTriads(entry.frames, { parentColor, childColor, parentLabel, childLabel }),
  });
}

/** The composite scene of `variants.composite_scene`, whose parts arrive placed. */
export function compositeScene(
  data,
  {
    parentColor,
    childColor,
    parentLabel = null,
    childLabel = null,
    ghostParent = false,
    showAtoms = true,
    planeColor = null,
    directionColor = null,
    planeAlpha = null,
  } = {},
) {
  const dress = (scene) => (showAtoms ? scene : frameOnly(scene));
  const translation = data.variant?.translation ?? [0, 0, 0];
  const apart = translation.some((value) => Math.abs(value) > 1e-9);
  const style = { planeColor, directionColor, alpha: planeAlpha };
  return mergeScenes(
    [
      tintScene(dress(data.parent?.scene), parentColor, ghostParent ? GHOST : SOLID),
      tintScene(dress(data.child?.scene), childColor, SOLID),
    ],
    {
      overlays: joinOverlays(
        primitiveOverlays(data.primitives, style),
        apart
          ? primitiveOverlays(data.primitives, { ...style, offset: translation, label: false })
          : null,
      ),
      extent: data.world,
      axes: data.parent?.scene?.axes,
      triads: frameTriads(data.variant?.frames, {
        parentColor,
        childColor,
        parentLabel,
        childLabel,
      }),
    },
  );
}

/**
 * The parent alone, as the wall's reference panel.
 *
 * Drawn solid and at the same camera and framing as every variant panel, so
 * "no rotation" is on screen beside the twelve rotations and the comparison is
 * made by the eye rather than from memory.
 */
export function parentReferenceScene(
  data,
  { parentColor, parentLabel = null, showAtoms = true, extent = null } = {},
) {
  const scene = showAtoms ? data.parent?.scene : frameOnly(data.parent?.scene);
  const frames = data.variants?.[0]?.frames;
  return mergeScenes([tintScene(scene, parentColor, SOLID)], {
    extent: extent ?? data.world,
    axes: data.parent?.scene?.axes,
    triads: frames
      ? [{ axes: frames.parent, color: parentColor, label: parentLabel, anchor: 'left' }]
      : null,
  });
}

/**
 * The composite of two *measured* grains, drawn in the specimen frame.
 *
 * Three differences from a catalogue composite, all of them the measurement's
 * doing rather than presentation choices.
 *
 * The world frame is the **specimen** frame the data arrived in, so the triad
 * is RD/TD/ND and comes from the service rather than from either crystal's
 * axes: a triad still labelled a, b, c would invite the picture to be read in
 * the wrong frame.
 *
 * Both structures arrive in their own crystal frames with one placement matrix
 * each, so the idealized child costs a matrix rather than a second copy of the
 * crystal.
 *
 * The overlays are already in the specimen frame and are drawn on **both**
 * sides, so the visible gap between a parent object and its child partner is
 * the clause deviation.
 */
export function measuredCompositeScene(data, { parentColor, childColor, idealColor } = {}) {
  const parts = [
    tintScene(
      placeScene(data.parent?.scene, data.parent?.matrix, data.parent?.translation),
      parentColor,
    ),
    tintScene(
      placeScene(data.child?.scene, data.child?.matrix, data.child?.translation),
      childColor,
    ),
  ];
  if (data.idealized) {
    parts.push(
      tintScene(
        placeScene(data.child?.scene, data.idealized.child_matrix, data.idealized.translation),
        idealColor,
      ),
    );
  }
  return mergeScenes(parts, {
    overlays: primitiveOverlays(data.primitives),
    extent: data.world,
    axes: data.world_axes ?? [],
  });
}

/** Degrees between two unit-ish vectors, taking the acute value. */
function acuteAngleDeg(left, right) {
  const dot =
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
  const norm =
    Math.hypot(left[0], left[1], left[2]) * Math.hypot(right[0], right[1], right[2]);
  if (norm < 1e-12) return Number.NaN;
  return (Math.acos(Math.min(Math.abs(dot) / norm, 1)) * 180) / Math.PI;
}

/**
 * What the OR overlays are doing relative to the screen, right now.
 *
 * The reader's question while turning a composite is "am I looking down it
 * yet?", and a picture cannot answer it: a plane a few degrees off edge-on
 * looks edge-on. The camera is a rotation from world to screen, so the world
 * vector pointing at the viewer is its transpose applied to the screen z axis,
 * and the rest is one angle per overlay.
 *
 * A plane is reported by how far it is from **edge-on**, because that is the
 * orientation in which a plane parallelism is legible — both lattices sit on
 * one line. A direction is reported by how far it is from pointing **at the
 * viewer**, which is when a direction parallelism is legible. The two are
 * ninety degrees apart, and labelling them identically would be the kind of
 * quiet error this readout exists to prevent.
 */
export function screenAlignment(cameraRotation, primitives) {
  const viewing = applyTranspose(cameraRotation, [0, 0, 1]);
  const notes = [];
  for (const patch of primitives?.patches ?? []) {
    const fromNormal = acuteAngleDeg(patch.normal, viewing);
    if (Number.isNaN(fromNormal)) continue;
    notes.push({
      kind: 'plane',
      label: patch.label ?? 'plane',
      // 0 when the plane contains the viewing direction, i.e. edge-on.
      angleDeg: Math.abs(90 - fromNormal),
      aligned: 'edge-on',
    });
  }
  for (const arrow of primitives?.arrows ?? []) {
    const vector = [
      arrow.head[0] - arrow.tail[0],
      arrow.head[1] - arrow.tail[1],
      arrow.head[2] - arrow.tail[2],
    ];
    const angle = acuteAngleDeg(vector, viewing);
    if (Number.isNaN(angle)) continue;
    notes.push({ kind: 'direction', label: arrow.label ?? 'direction', angleDeg: angle, aligned: 'end-on' });
  }
  return notes;
}
