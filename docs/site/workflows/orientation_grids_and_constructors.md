# Orientation Grids And Constructors

PyTex provides a first MTEX-inspired orientation utility surface for deterministic orientation
grids, named specimen directions, and common orientation constructors. This is a foundational
implementation, not an MTEX parity claim.

## Orientation Grids

Use `OrientationSet.from_so2_grid(...)` for rotations about one specimen axis:

```python
orientations = OrientationSet.from_so2_grid(
    "ND",
    5.0,
    specimen_frame=specimen,
    phase=phase,
)
```

Use `OrientationSet.from_regular_so3_grid(...)` for a deterministic Bunge Euler grid:

```python
orientations = OrientationSet.from_regular_so3_grid(
    15.0,
    specimen_frame=specimen,
    phase=phase,
)
```

Use `OrientationSet.from_equispaced_so3_grid(...)` for a deterministic lower-discrepancy
quaternion grid over SO(3):

```python
orientations = OrientationSet.from_equispaced_so3_grid(
    15.0,
    specimen_frame=specimen,
    phase=phase,
    reduce_to_fundamental_region=True,
)
```

The regular SO(3) grid is explicitly a Bunge-grid construction. The equispaced SO(3) grid uses a
deterministic quaternion sequence, then optionally reduces and deduplicates orientations through
the active crystal symmetry.

## Named Specimen Directions

The helper `specimen_direction_vector(...)` accepts the canonical EBSD-style aliases:

| Alias | Vector |
| --- | --- |
| `RD`, `x` | `[1, 0, 0]` |
| `TD`, `y` | `[0, 1, 0]` |
| `ND`, `z` | `[0, 0, 1]` |

These aliases are accepted by the SO2 grid, IPF color helpers, and Miller notation constructors.

## Orientation Constructors

The scalar `Orientation` surface now mirrors the existing batch constructors:

```python
ori = Orientation.from_euler(35.0, 25.0, 10.0, specimen_frame=specimen, phase=phase)
ori = Orientation.from_axis_angle("ND", np.pi / 2.0, specimen_frame=specimen, phase=phase)
ori = Orientation.from_matrix(matrix, specimen_frame=specimen, phase=phase)
ori = Orientation.from_quaternion(quaternion, specimen_frame=specimen, phase=phase)
```

For ideal orientation notation, use `Orientation.from_miller(...)`:

```python
cube = Orientation.from_miller(
    [0, 0, 1],
    [1, 0, 0],
    specimen_frame=specimen,
    phase=phase,
    specimen_plane_normal="ND",
    specimen_direction="RD",
)
```

This maps the crystal plane normal to `ND` and the crystal direction to `RD`. The implementation
uses the same plane/direction backend as `OrientationSet.from_plane_direction(...)`.

## Current Limits

- MTEX validation is intentionally paused for this slice.
- `from_equispaced_so3_grid(...)` is deterministic and symmetry-aware, but it is not a reproduction
  of MTEX `equispacedSO3Grid`.
- Vendor-specific IPF color keys and ideal-component catalogs remain planned follow-up work.

## Next MTEX-Inspired Utilities

- general `orientation.map`-style two-vector/two-vector constructors
- fiber generation from a crystal direction and specimen direction
- low-index interpretation helpers similar in spirit to `round2Miller`
- named ideal texture components such as Cube, Goss, Brass, Copper, and S
- orientation scatter subsampling helpers for PF, IPF, and ODF-section plotting
- documented vendor-style IPF color-key variants

## Related Material

- {doc}`ipf_colors`
- {doc}`../concepts/orientation_texture`
- `docs/testing/mtex_parity_matrix.md`
