# Miller Planes And Directions

`MillerPlane`, `MillerDirection`, `MillerPlaneSet`, and `MillerDirectionSet` are the first-class
PyTex surface for crystallographic plane and direction indices tied to a `Phase`.

These objects exist because the semantic meaning of a plane or direction is not just an integer
triplet. The same indices only become scientifically meaningful once the direct basis, reciprocal
basis, symmetry operators, and hexagonal convention policy are fixed by the owning phase.

## Core Semantics

- `MillerPlane` and `MillerPlaneSet` store exact three-index reciprocal-space plane indices
  `(h, k, l)` as `int64`.
- `MillerDirection` and `MillerDirectionSet` store exact three-index direct-space direction indices
  `(u, v, w)` as `int64`.
- Constructors do not silently reduce indices. PyTex preserves exact input order information, so
  `(2, 0, 0)` remains distinct from `(1, 0, 0)` until a family or canonicalization utility is
  requested explicitly.
- Cartesian vectors are always `float64`.
- Reciprocal vectors use the crystallographer convention with no `2\pi`, so reciprocal units are
  `1/angstrom`.
- `d_spacing_angstrom` and `d_spacings_angstrom()` always return angstroms.
- Angular outputs are always radians.

## Hexagonal And Trigonal Conventions

PyTex stores three-index internal indices and exposes vectorized conversion helpers for the
four-index hexagonal forms:

- planes: `(h, k, l) <-> (h, k, i, l)` with `i = -(h + k)`
- directions: `(u, v, w) <-> (U, V, T, W)` with `U + V + T = 0`

Use:

- `MillerPlane.from_hkil(...)`
- `MillerDirection.from_UVTW(...)`
- `plane_hkl_to_hkil_array(...)`
- `plane_hkil_to_hkl_array(...)`
- `direction_uvw_to_uvtw_array(...)`
- `direction_uvtw_to_uvw_array(...)`

The conversion helpers are vectorized and use integer arithmetic throughout.

## Antipodal And Family Semantics

PyTex separates exact stored indices from family semantics:

- `reduce_indices(...)` performs row-wise integer GCD reduction.
- `canonicalize_sign(...)` makes the first non-zero entry positive.
- `antipodal_keys(...)` reduces and sign-canonicalizes.

Default family rules:

- planes are always treated antipodally
- directions default to antipodal family handling, but direction-angle and direction-family calls
  expose `antipodal=` when the sign distinction matters

This means plane families such as `{100}` deduplicate to three representatives under cubic proper
symmetry, while direction families such as `<100>` can be treated either antipodally or as signed
axes depending on the operation.

## Vectorized Operations

`MillerPlaneSet` and `MillerDirectionSet` are the preferred batch surface when you have `n`
indices:

- `MillerPlaneSet.reciprocal_vectors_cartesian() -> (n, 3)`
- `MillerPlaneSet.normals_cartesian() -> (n, 3)`
- `MillerPlaneSet.d_spacings_angstrom() -> (n,)`
- `MillerDirectionSet.direct_vectors_cartesian() -> (n, 3)`
- `MillerDirectionSet.unit_vectors_cartesian() -> (n, 3)`
- `angle_plane_plane_rad(...)`
- `angle_dir_dir_rad(...)`
- `angle_dir_plane_normal_rad(...)`
- `angle_dir_plane_inclination_rad(...)`
- `project_directions_onto_planes(...) -> ((n, 3), (n,))`

Projection returns both the projected cartesian vectors and a boolean degenerate mask. Rows whose
projection collapses to zero are reported through the mask rather than failing the whole batch.

## Symmetry Expansion

`symmetry_equivalent_indices(...)` applies `Phase.symmetry` to non-normalized lattice vectors and
then converts back to integer indices with explicit rounding tolerance checks.

Two public levels are exposed:

- scalar objects return a single `MillerPlaneSet` or `MillerDirectionSet`
- batch objects provide a vectorized padded integer surface:
  `symmetry_equivalent_indices(...) -> (indices, mask)`

The padded `(n, m, 3)` plus mask representation keeps the core math vectorized even when each input
row expands to a variable family size.

## Compatibility With Existing PyTex Objects

The Miller surface complements the older scalar crystallographic helpers rather than replacing them:

- `MillerPlane.from_crystal_plane(...)`
- `MillerDirection.from_zone_axis(...)`
- `MillerPlane.to_crystal_plane()`
- `MillerDirection.to_zone_axis()`

Use `CrystalPlane` and `ZoneAxis` when you need the older scalar workflow surface. Use Miller
objects when you need vectorized family math, explicit batch APIs, or JSON-contract round trips.

## Related Pages

- {doc}`core_model`
- {doc}`reference_frames_and_conventions`
- {doc}`orientation_texture`
- {doc}`../workflows/vectorized_miller_workflows`
- {doc}`../workflows/orix_kikuchipy_interop`
- {doc}`../theory/index`

## References

### Normative

- {doc}`../architecture/canonical_data_model`
- {doc}`../standards/notation_and_conventions`
- {doc}`../standards/hexagonal_and_trigonal_conventions`

### Informative

- MTEX documentation: <https://mtex-toolbox.github.io/>
- orix vector and Miller documentation: <https://orix.readthedocs.io/>
- KikuchiPy documentation: <https://kikuchipy.org/>
