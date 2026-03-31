# EBSD To Texture Outputs (ODF, PF, And IPF)

PyTex now exposes one coherent path from orientation construction and EBSD normalization to
texture outputs:

- `OrientationSet.from_plane_direction(...)`
- `normalize_ebsd(...)`
- `index_hough(...)`
- `refine_orientations(...)`
- `CrystalMap.to_odf(...)`
- `CrystalMap.pole_figure(...)`
- `CrystalMap.inverse_pole_figure(...)`
- `CrystalMap.texture_report(...)`
- `plot_ipf_map(...)`
- `plot_kam_map(...)`

## Orientation Construction From Plane/Direction Pairs

The plane/direction constructor builds a crystal-to-specimen orientation from one crystal plane and
one crystal direction. By default:

- the crystal plane normal maps to specimen `Z`
- the in-plane crystal direction maps to specimen `X`

If the supplied crystal direction is not exactly in the plane, PyTex projects it into the plane
before constructing the right-handed basis.

## Normalized EBSD Surface

`normalize_ebsd(...)` accepts either a KikuchiPy signal carrying `xmap` or a crystal-map-like
object with quaternion rotations and map coordinates. The returned object is always a
`NormalizedEBSDDataset`.

For multiphase imports, provide a `phases` lookup keyed by phase id or phase name so PyTex can
attach real `Phase` objects while preserving per-point phase assignments.

## Producing Texture Outputs

```python
crystal_map = dataset.crystal_map

odf = crystal_map.to_odf()
pole_figure = crystal_map.pole_figure([1, 0, 0])
inverse_pole_figure = crystal_map.inverse_pole_figure("z")
```

On a multiphase map, texture extraction requires an explicit phase selector:

```python
ferrite_report = crystal_map.texture_report(
    phase="ferrite",
    poles=([1, 0, 0], [1, 1, 1]),
    sample_directions=("x", "z"),
    plot=True,
)
```

## Plotting The Map Surface

```python
from pytex import plot_ipf_map, plot_kam_map

ipf_figure = plot_ipf_map(crystal_map, direction="z")
kam_figure = plot_kam_map(crystal_map, order=1, symmetry_aware=True)
```

Both plotting calls accept optional grain-boundary overlays.

## Convention Notes

- `Orientation` and `OrientationSet` remain crystal-to-specimen mappings
- quaternion order remains `w, x, y, z`
- axis-angle outputs use radians
- Rodrigues coordinates use `n tan(ω / 2)` with the `ω = π` singularity kept explicit

## Related Material

- {doc}`ebsd_import_normalization`
- {doc}`ebsd_kam`
- {doc}`texture_odf_inversion`
- [../../tex/algorithms/orientation_representations_and_plane_direction_construction.tex](../../tex/algorithms/orientation_representations_and_plane_direction_construction.tex)
- [../../tex/algorithms/multiphase_ebsd_graph_workflows.tex](../../tex/algorithms/multiphase_ebsd_graph_workflows.tex)
