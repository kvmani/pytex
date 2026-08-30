# Crystal Visualization Geometry

PyTex treats crystal visualization as a rendering layer over the canonical crystallographic model, not as a second private geometry system.

## Atomic Coordinates

Atomic sites are stored in fractional coordinates relative to the unit cell. For a direct basis matrix

$$
\mathbf{A} = [\mathbf{a}\ \mathbf{b}\ \mathbf{c}]
$$

a site at fractional coordinates $\mathbf{f}$ is rendered at Cartesian crystal-space position

$$
\mathbf{r} = \mathbf{A}\mathbf{f}
$$

Supercell repetition uses explicit integer translations in the direct-lattice basis:

$$
\mathbf{r}_{ijk} = \mathbf{A}(\mathbf{f} + [i, j, k]^T)
$$

## Lattice, Cell Overlays, Planes, And Directions

The current viewer draws lattice edges from the repeated supercell box. Optional unit-cell overlays are generated directly from translated copies of the direct-basis parallelepiped. For a cell-corner anchor $\mathbf{f}_0$ and integer span vector $\mathbf{s}$, the overlay vertices are

$$
\mathbf{r}_{\mathrm{cell}} = \mathbf{A}(\mathbf{f}_0 + \boldsymbol{\delta})
$$

where $\boldsymbol{\delta}$ ranges over the eight parallelepiped corners implied by $\mathbf{s}$.

## Plane overlays are cut by the cell

A plane overlay is drawn as the polygon in which the lattice plane meets the cell box, not as a
patch of arbitrary size placed at the plane's orientation. In fractional coordinates a lattice
plane of a family $(hkl)$ is exactly

$$h x_1 + k x_2 + l x_3 = m, \qquad m \in \mathbb{Z},$$

so the polygon is obtained by intersecting that plane with the twelve edges of the box and ordering
the intersections about the plane normal. Members with fewer than three distinct intersections —
those touching only an edge or a corner — are degenerate and are rejected.

Which $m$ is drawn matters, because the members of a family are not congruent in general. PyTex
takes the member of **largest cross-sectional area** through the box, breaking ties toward the box
centre and then toward the larger offset. For a cubic cell and $(110)$ this is the diagonal
rectangle through two opposite edges; for $(100)$, where the two members are congruent faces, it is
the far face.

A direction drawn alongside a plane is clipped to the *polygon*, giving the chord through its
centroid. This is a statement of geometry rather than of style: the direction of an orientation
relationship lies in its plane, and a chord exhibits that, while a segment of arbitrary length
anchored at the origin does not.

For hexagonal-axis lattices, PyTex can also render an auxiliary hexagonal prism for teaching and visual interpretation. Let the direct basis be

$$
\mathbf{A} = [\mathbf{a}\ \mathbf{b}\ \mathbf{c}]
$$

with $|\mathbf{a}| = |\mathbf{b}|$, $\alpha = \beta = 90^{\circ}$, and $\gamma = 120^{\circ}$. From a lattice-point anchor $\mathbf{r}_0$, the basal hexagon is formed from the six vertices

$$
\mathbf{r}_0 \pm \mathbf{a},\qquad
\mathbf{r}_0 \pm \mathbf{b},\qquad
\mathbf{r}_0 \pm (\mathbf{a}+\mathbf{b})
$$

ordered cyclically in the basal plane and extruded by integer multiples of $\mathbf{c}$.

The prism is also the region plane overlays are clipped to when it is drawn, so a basal plane
appears as the hexagon rather than as the rhombus inside it, and its axis is placed through an
atomic column so that the six corner columns are occupied — with the axis on the cell origin, a
phase whose sites sit at $(1/3, 2/3)$ draws a prism with empty corners.

This hexagonal prism is an explicit visualization overlay, not a replacement for the canonical cell semantics of the PyTex data model. It is useful because it exposes the sixfold basal symmetry more directly than the direct-basis parallelepiped, but users should not confuse it with the primitive cell stored by the library.

Plane overlays are defined from Miller indices $(hkl)$ and the reciprocal-space normal

$$
\mathbf{n}_{hkl} = h \mathbf{a}^{*} + k \mathbf{b}^{*} + l \mathbf{c}^{*}
$$

The viewer computes the polygon formed by intersecting that plane with the repeated cell bounding box. This yields a bounded geometric overlay suitable for teaching and publication figures while preserving the exact lattice and reciprocal conventions from the core model.

Direction overlays are defined from a crystallographic direction vector and an explicit fractional anchor point inside the repeated cell. If the direct-basis fractional ray is

$$
\mathbf{f}(t) = \mathbf{f}_0 + t\mathbf{u}
$$

PyTex clips the ray to the repeated-cell volume and renders the bounded segment in Cartesian crystal coordinates. This keeps direction graphics tied to the same lattice semantics as the structure itself.

## Scientific Labels

Plane and direction annotations use a shared Miller-notation formatter. Negative indices are rendered with bar notation, for example

$$
[11\bar{2}0], \qquad (11\bar{2}1)
$$

so graphics, documentation, and future figure exporters can share one scientific text convention rather than encoding ad hoc minus-sign styles per subsystem.

## View Control

Camera angles are rendering choices, not crystallographic semantics. When PyTex accepts a `CrystalDirection` as a view target, the direction is first interpreted in crystal coordinates and only then converted into visualization-camera angles. This keeps the distinction between scientific direction meaning and graphics control explicit.

## Current Limits

- The bond model is heuristic and based on covalent-radius proximity rather than exhaustive chemistry rules.
- The current backend is a static Matplotlib 3D renderer.
- Slab views are currently implemented as deterministic geometric filtering rather than a full clipping-volume engine.
