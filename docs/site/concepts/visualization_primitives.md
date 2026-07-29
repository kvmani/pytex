# Visualization Primitives And Scene Composition

PyTex draws crystallographic geometry from a small set of **renderer-independent primitives** so
that every figure -- from a single vector to two crystals in an orientation relationship -- is
assembled from the same vocabulary with minimal code. This page defines that vocabulary, the
placement model that positions primitives in a shared world frame, and the composition surfaces that
build real-space 3D and stereographic figures on top of them.

The primitives live in `pytex.plotting.primitives`; the composite world scene lives in
`pytex.plotting.scene3d`; the stereographic bridge lives in `pytex.plotting.spherical`. All are
re-exported from the top-level `pytex` namespace.

## Why A Primitive Layer

The single-crystal renderer (`build_crystal_scene` / `plot_crystal_structure_3d`) is complete and
VESTA-class, but it models exactly one crystal in its own Cartesian frame. Composite scientific
figures need more:

- a **bare 3D vector**, a **Miller direction**, a **plane pole**, or a **reference-frame triad** as a
  first-class element that is not tied to an atomistic structure;
- a **placement** so one crystal (or primitive group) can sit at an `Orientation`/`Rotation`
  relative to another;
- a **world scene** that holds several placed crystals plus loose primitives and renders them
  together with globally correct depth sorting.

The primitive layer supplies all three without duplicating the proven single-crystal renderer: it
reuses the same lit-mesh accumulation, so a placed crystal renders identically to one built in its
own frame.

## The Primitive Vocabulary

Every primitive is an immutable dataclass carrying only geometry (world Cartesian, angstrom), colors,
and annotation intent. Nothing is matplotlib-specific, so the same scene can drive any future
backend.

| Primitive | Represents |
| --- | --- |
| `Arrow3D` | a 3D vector drawn tail → head (a bare vector, a Miller direction, a plane pole, a triad axis) |
| `PolyLine3D` | an open or closed polyline (unit-cell edges, projected traces) |
| `PlanePatch3D` | a planar polygon with an outward normal (a lattice plane, slip plane, habit plane) |
| `PointCloud3D` | point markers (lattice nodes, atom sites, poles) |
| `Label3D` | a text label anchored at a world point |
| `AxisTriad3D` | an orthonormal axis triad (specimen RD/TD/ND, crystal a/b/c), expanded to three arrows plus labels |

A `PrimitiveScene3D` is an immutable bag of these; `merge` composes two scenes and `transformed`
places an entire scene at once. `render_primitive_scene_3d` draws a scene standalone (a lone vector,
a triad, a unit-cell wireframe), and the composite world renderer reuses the identical drawing logic.

### Builders From Crystallographic Objects

So the crystallographic surface and the drawing surface share one language, builders turn canonical
objects into primitives:

- `vector_arrow(vector, ...)` -- a bare Cartesian vector.
- `direction_arrow(CrystalDirection, ...)` -- a Miller direction, auto-labelled `[uvw]`.
- `plane_normal_arrow(CrystalPlane, ...)` -- a plane pole.
- `crystal_plane_patch(CrystalPlane, ...)` -- a translucent lattice-plane sheet, auto-labelled `(hkl)`.
- `reference_frame_triad(frame_or_basis, ...)` -- a frame gizmo (crystal `a/b/c` from a direct basis,
  or specimen `x/y/z` from a `ReferenceFrame`).
- `unit_cell_polylines(Phase | Lattice, ...)` -- the twelve edges of a (super)cell.
- `lattice_point_cloud(Phase | Lattice, ...)` -- Bravais lattice nodes over a supercell block.

## Placement: `Transform3D`

A `Transform3D` maps geometry into the world frame as
$\mathbf{x}_\text{world} = \mathbf{T}(\mathbf{x}_\text{local}) = \mathbf{R}\,\mathbf{x}_\text{local} + \mathbf{t}$. Construct it from a
`Rotation`, an `Orientation` (`from_orientation`, which reproduces
`Orientation.map_crystal_vector` so a crystal drawn in the sample frame lands where the orientation
says), or an explicit matrix. It provides `apply_points`, `apply_vector`, `apply_normal` (the
covariant inverse-transpose for plane normals), `compose`, and `inverse`.

Crystal placement requires a **rigid** transform (a rotation plus a translation) so atom spheres and
bond cylinders stay undistorted; `CrystalScene.transformed` enforces this.

## Composite World Scenes

`WorldScene3D` is an immutable composite of `PlacedCrystal`s (each a `CrystalScene` plus a
`Transform3D`) and a `PrimitiveScene3D`. Grow it functionally with `add_crystal` and
`add_primitives`. `render_world_scene_3d` accumulates every placed crystal's atom, bond, and
polyhedron faces into **one** depth-sorted collection, so crystals occlude one another correctly from
any viewing angle; lattice frames, plane patches, direction arrows, and the loose primitives are
drawn on top.

### Two Crystals In An Orientation Relationship

`WorldScene3D.from_orientation_relationship(relationship, ...)` is the minimal-code entry point for
the canonical composite figure. It places the parent crystal in the world frame and the child by
`relationship.parent_to_child_rotation.inverse()`, which makes the relationship's parallel planes and
directions coincide in world coordinates -- the geometric statement of the relationship, shown
directly. With the defaults the parallel directions are drawn as arrows and the parallel planes as
translucent patches, so the alignment is visible.

The correctness of this placement is checked by the worked example
{doc}`../examples/generated/visualization`: after placement, the Kurdjumov-Sachs parallel directions
have direction cosine 1.

## VESTA-Class Single-Crystal Rendering

The atomistic renderer targets — and is ledgered against — VESTA, the reference desktop application
for publication crystal graphics (see the {doc}`../validation/vesta_parity_matrix`). One keyword
switches the whole visual system:

- `render_style="ball_and_stick"` (default) — lit spheres and two-tone cylinders in one globally
  depth-sorted mesh.
- `render_style="space_filling"` — Slater atomic radii at full scale, bonds suppressed.
- `render_style="stick"` — uniform thin cylinders with matching atom caps.
- `render_style="wireframe"` — line bonds and marker atoms (no lit mesh).
- `render_style="polyhedral"` — coordination polyhedra for every eligible species.

VESTA-signature behaviors are automatic or one parameter away:

- **Partial occupancy**: sites sharing one position (mixed species) render as azimuthal pie sectors
  of one sphere; occupancy below one leaves a vacancy sector in the theme `vacancy_color`.
- **Atom labels**: `atom_label_mode="species"` or `"site"`.
- **Vectors on atoms** (magnetic moments, displacements): `site_vectors={site_label: vector}` draws
  the arrow on every periodic copy of that site.
- **Depth cueing**: the `depth_cue_strength` theme key fades distant faces toward the background
  along the view direction.
- **Measurement**: `CrystalScene.bond_lengths_angstrom()` and `bond_length_summary()` replace
  click-driven distance readout with scriptable, regression-tested numbers.

## Stereographic And IPF Projection Of The Same Primitives

`plot_stereographic_vectors(vectors, ...)` projects arbitrary world-frame directions as poles and/or
great-circle traces on a Wulff net. Because the inputs are plain Cartesian vectors, directions from
two crystals placed in a common frame with a `Transform3D` overlay on one stereogram -- the
projection analog of a composite 3D scene. Crystallographic wrappers (`plot_crystal_directions`,
`plot_crystal_planes`) and the IPF colour key (`IPFColorKey`, `plot_ipf_key`) remain the phase-aware
surfaces for stereographic and inverse-pole-figure work.

## See Also

- {doc}`miller_planes_directions` -- the `CrystalDirection` and `CrystalPlane` inputs to the builders.
- {doc}`reference_frames_and_conventions` -- the frame semantics that `Transform3D` places between.
- {doc}`orientation_texture` -- orientations, rotations, and orientation relationships.
- {doc}`../examples/generated/visualization` -- the executable checks of the placement geometry.
- The visualization style guide in `docs/standards/visualization_style_guide.md`.
