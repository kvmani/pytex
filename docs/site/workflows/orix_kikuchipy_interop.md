# orix And KikuchiPy Interoperability

PyTex keeps ORIX and KikuchiPy behind adapters so external-tool semantics do not leak into the
stable public model.

## What PyTex Preserves

- the owning `Phase`
- the phase symmetry and crystal frame
- scalar-first quaternion ordering `(w, x, y, z)`
- crystallographer reciprocal-space convention with no `2π`
- explicit phase tables and `phase_id` arrays for multiphase crystal maps

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
```

## Rotation And Orientation Bridges

Use `to_orix_rotation(...)`, `from_orix_rotation(...)`, `to_orix_orientation(...)`, and
`from_orix_orientation(...)` for quaternion-backed rotation transfer.

## Multiphase `xmap` Flow

When a KikuchiPy or ORIX crystal map exposes `phases` plus `phase_id` values, PyTex now preserves
that phase-resolved structure directly:

```python
dataset = normalize_ebsd(
    xmap_or_signal,
    frames=(crystal_frame, specimen_frame, map_frame),
    phases={1: ferrite_phase, 2: austenite_phase},
)

phase_table = dataset.crystal_map.resolved_phase_entries
phase_ids = dataset.crystal_map.phase_id_array
```

Texture extraction from such a map requires an explicit `phase=` selector:

```python
odf = dataset.crystal_map.to_odf(phase="ferrite")
```

## Why The Adapter Boundary Is Explicit

- PyTex keeps one canonical frame and symmetry model
- optional dependencies remain optional
- JSON contracts and manifests remain reconstructible without hidden ORIX state

## Related Pages

- {doc}`../concepts/miller_planes_directions`
- {doc}`ebsd_import_normalization`
- {doc}`ebsd_to_texture_outputs`
