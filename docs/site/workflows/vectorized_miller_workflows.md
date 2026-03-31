# Vectorized Miller Workflows

This page shows the intended batch workflow for first-class Miller planes and directions in PyTex.

The guiding idea is simple: keep indices attached to a `Phase`, keep exact integer input intact,
and use vectorized set APIs for geometry, families, and interop.

## Create Plane And Direction Sets

```python
import numpy as np

from pytex import MillerDirectionSet, MillerPlaneSet

planes = MillerPlaneSet.from_hkl(
    np.array(
        [
            [1, 0, 0],
            [1, 1, 0],
            [1, 1, 1],
        ],
        dtype=np.int64,
    ),
    phase=phase,
)

directions = MillerDirectionSet.from_uvw(
    np.array(
        [
            [1, 0, 0],
            [1, 1, 0],
            [1, 1, 1],
        ],
        dtype=np.int64,
    ),
    phase=phase,
)
```

Hexagonal forms can be created directly:

```python
basal = MillerPlaneSet.from_hkil([[0, 0, 0, 1], [1, 0, -1, 0]], phase=phase)
slip = MillerDirectionSet.from_UVTW([[2, -1, -1, 0], [1, 0, -1, 0]], phase=phase)
```

## Cartesian Vectors, Normals, And d-Spacings

```python
reciprocal_vectors = planes.reciprocal_vectors_cartesian()
plane_normals = planes.normals_cartesian()
d_spacings = planes.d_spacings_angstrom()

direct_vectors = directions.direct_vectors_cartesian()
direction_units = directions.unit_vectors_cartesian()
```

Units:

- reciprocal vectors: `1/angstrom`
- direct vectors: angstrom
- `d_spacings`: angstrom
- angles: radians

## Family Canonicalization

Exact indices are preserved by construction. Family reduction is explicit:

```python
reduced = directions.reduce_indices()
canonical = directions.canonical_indices(antipodal=True)
unique_directions, inverse = directions.unique(antipodal=True)
```

For planes, `unique()` always uses antipodal family semantics.

## Symmetry Expansion

```python
families = planes.symmetry_equivalents()
equivalent_indices, mask = planes.symmetry_equivalent_indices()
```

`families` is a tuple of `MillerPlaneSet` objects. `equivalent_indices` is the vectorized padded
surface with shape `(n, m, 3)` and a matching boolean mask.

This split is deliberate:

- tuple-returning helpers are ergonomic in interactive work
- padded indices plus mask keep the core symmetry math vectorized for large batches

## Angle And Projection Operations

```python
plane_angles = planes.angle_matrix_rad()
direction_angles = directions.angle_matrix_rad(antipodal=True)

dir_plane_normal = angle_dir_plane_normal_rad(directions, planes)
dir_plane_inclination = angle_dir_plane_inclination_rad(directions, planes)
projected, degenerate = project_directions_onto_planes(directions, planes)
```

`project_directions_onto_planes(...)` returns cartesian projected vectors and a per-row degenerate
mask. This makes mixed-validity batches safe for downstream reporting.

## Compatibility With Older Scalar Objects

```python
from pytex import CrystalPlane, ZoneAxis

plane = MillerPlane.from_crystal_plane(
    CrystalPlane.from_miller_bravais((1, 0, -1, 0), phase=phase)
)
axis = MillerDirection.from_zone_axis(ZoneAxis([1, 0, 0], phase=phase))

legacy_plane = plane.to_crystal_plane()
legacy_axis = axis.to_zone_axis()
```

## Validation And Limits

- zero triplets are rejected at construction time
- invalid hexagonal four-index constraints raise `ValueError`
- phase mismatches raise `ValueError`
- symmetry back-conversion to integer indices is tolerance-checked
- plane family operations are always antipodal

## Related Pages

- {doc}`../concepts/miller_planes_directions`
- {doc}`orix_kikuchipy_interop`
- {doc}`ebsd_to_texture_outputs`

## References

### Normative

- `../../standards/notation_and_conventions.md`
- `../../architecture/canonical_data_model.md`

### Informative

- MTEX documentation: <https://mtex-toolbox.github.io/>
