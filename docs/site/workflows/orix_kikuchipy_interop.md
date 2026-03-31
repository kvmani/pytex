# orix And KikuchiPy Interoperability

PyTex keeps ORIX and KikuchiPy behind adapters so external-tool semantics do not leak into the
stable public model. This page records the intended bridge for rotations, EBSD maps, and the new
first-class Miller plane and direction surface.

## What PyTex Preserves

Across the adapter boundary, PyTex keeps the following explicit:

- the owning `Phase`
- the phase symmetry and crystal frame
- scalar-first quaternion ordering `(w, x, y, z)`
- crystallographer reciprocal-space convention with no `2\pi`

## Miller Bridges

Use the dedicated Miller adapter functions when you need `orix.vector.Miller` objects:

```python
from pytex import (
    MillerDirection,
    MillerPlane,
    from_orix_miller,
    to_orix_miller_direction,
    to_orix_miller_plane,
    to_orix_phase,
)

plane = MillerPlane.from_hkl([1, 1, 1], phase=phase)
direction = MillerDirection.from_uvw([1, 1, 0], phase=phase)

orix_phase = to_orix_phase(phase)
orix_plane = to_orix_miller_plane(plane)
orix_direction = to_orix_miller_direction(direction)

plane_round_trip = from_orix_miller(orix_plane, phase=phase)
direction_round_trip = from_orix_miller(orix_direction, phase=phase)
```

`from_orix_miller(...)` uses the ORIX `coordinate_format` when available so planes and directions
are not misclassified even though ORIX can expose multiple equivalent index representations from the
same cartesian storage.

## Rotation And Orientation Bridges

The existing rotation and orientation adapters remain the preferred path for EBSD map rotations:

```python
from pytex import from_orix_orientation, from_orix_rotation, to_orix_orientation, to_orix_rotation
```

PyTex and ORIX both use scalar-first quaternions for these adapters. PyTex tests keep that
ordering explicit because quaternion component order is a common interop failure mode.

## KikuchiPy CrystalMap Flow

KikuchiPy indexing results commonly expose an ORIX `CrystalMap` or ORIX-backed `xmap` carrying:

- ORIX rotations
- ORIX phases
- map coordinates and scan geometry

The PyTex workflow is:

1. normalize the EBSD payload with `normalize_ebsd(...)`
2. convert rotations through the existing ORIX rotation or orientation adapters
3. use the phase-normalized `CrystalMap` surface for KAM, grains, PF, IPF, or ODF outputs
4. use the Miller adapter only when an external ORIX or KikuchiPy step needs explicit plane or
   direction vectors

## Why The Adapter Is Explicit

PyTex does not expose raw ORIX objects as stable public state. The adapter boundary is deliberate:

- PyTex keeps one canonical frame and symmetry model
- optional dependencies remain optional
- JSON contracts and manifests stay reconstructible without hidden ORIX state

## Related Pages

- {doc}`../concepts/miller_planes_directions`
- {doc}`vectorized_miller_workflows`
- {doc}`ebsd_import_normalization`
- {doc}`ebsd_to_texture_outputs`

## References

### Normative

- `../../standards/data_contracts_and_manifests.md`
- `../../architecture/canonical_data_model.md`

### Informative

- ORIX documentation: <https://orix.readthedocs.io/>
- KikuchiPy documentation: <https://kikuchipy.org/>
