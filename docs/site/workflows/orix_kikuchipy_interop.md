# orix And KikuchiPy Interoperability

PyTex keeps ORIX and KikuchiPy behind adapters so external-tool semantics do not leak into the
stable public model.

## What Is Executable Today

The current hardening phase treats this workflow as a narrow interoperability boundary, not as a
claim that PyTex is a drop-in replacement for either external package.

Today PyTex explicitly supports:

- Miller plane and direction transfer through dedicated ORIX adapter helpers
- scalar-first quaternion transfer for rotations and orientations
- point-group and phase preservation at the adapter boundary
- multiphase crystal-map normalization when the incoming payload exposes a phase table plus
  per-point `phase_id` values

These paths are covered by executable tests in:

- `tests/unit/test_orix_miller_adapter.py`
- `tests/unit/test_plotting_and_adapters.py`

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

## What Is Not Being Claimed

- PyTex does not claim broad live-package parity across all ORIX or KikuchiPy workflow surfaces.
- Support for an adapter function does not imply validation of every higher-level workflow that can
  be built on top of it.
- KikuchiPy detector-pattern semantics and package-version-specific behavior are still outside this
  page's validated scope.

## Current Limits

- The validated ORIX boundary is strongest on Miller, rotation, orientation, symmetry, and phase
  transfer, not on the full ORIX workflow stack.
- KikuchiPy-facing normalization is validated primarily through PyTex canonicalization behavior and
  xmap-like payload normalization, not through a full dependency-pinned external environment matrix.
- If you need release-specific package interoperability guarantees, treat the tests above as the
  current controlling evidence and add environment-pinned integration coverage before making
  stronger claims.

## Related Pages

- {doc}`../concepts/miller_planes_directions`
- {doc}`ebsd_import_normalization`
- {doc}`ebsd_to_texture_outputs`
