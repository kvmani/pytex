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
export function tintScene(scene, color) {
  if (!scene) return scene;
  return {
    ...scene,
    atoms: (scene.atoms ?? []).map((atom) => ({ ...atom, color })),
    bonds: (scene.bonds ?? []).map((bond) => ({
      ...bond,
      color,
      start_color: color,
      end_color: color,
    })),
  };
}


/**
 * The world-frame OR overlays in the shape the renderer already draws.
 *
 * An arrow is a direction and a patch is a plane; the service emits the same
 * keys, so this is a rename rather than a conversion. Keeping the names apart
 * on the wire is deliberate — a scene's planes belong to a crystal, while these
 * belong to the relationship between two of them.
 */
export function primitiveOverlays(primitives) {
  return {
    planes: (primitives?.patches ?? []).map((patch) => ({
      vertices: patch.vertices,
      normal: patch.normal,
      color: patch.color,
      alpha: patch.alpha,
      label: patch.label,
    })),
    directions: (primitives?.arrows ?? []).map((arrow) => ({
      start: arrow.tail,
      end: arrow.head,
      color: arrow.color,
      label: arrow.label,
    })),
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
export function mergeScenes(scenes, { overlays = null, extent = null, axes = null } = {}) {
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
    centre: extent?.centre ?? parts[0]?.centre ?? [0, 0, 0],
    radius: extent?.radius ?? parts[0]?.radius ?? 1,
    bounds: extent?.bounds ?? parts[0]?.bounds ?? null,
  };
  return merged;
}

/** The composite scene of `variants.composite_scene`, whose parts arrive placed. */
export function compositeScene(data, { parentColor, childColor }) {
  return mergeScenes(
    [tintScene(data.parent?.scene, parentColor), tintScene(data.child?.scene, childColor)],
    {
      overlays: primitiveOverlays(data.primitives),
      extent: data.world,
      axes: data.parent?.scene?.axes,
    },
  );
}

/** One panel of `variants.contact_sheet`: the parent, plus this variant's child. */
export function contactSheetScene(data, entry, { parentColor, childColor }) {
  const child = placeScene(data.child?.scene, entry.child_matrix, entry.translation);
  return mergeScenes(
    [tintScene(data.parent?.scene, parentColor), tintScene(child, childColor)],
    {
      overlays: primitiveOverlays(entry.primitives),
      extent: data.world,
      axes: data.parent?.scene?.axes,
    },
  );
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
