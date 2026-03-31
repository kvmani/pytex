# EBSD To Texture Outputs (ODF, PF, And IPF)

PyTex now exposes one coherent path from orientation construction and EBSD normalization to texture outputs:

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

The important point is not just convenience. Each step keeps the PyTex semantic model intact: crystal-to-specimen mapping direction, phase attachment, crystal symmetry, specimen frame, and provenance-bearing normalization boundaries.

## Orientation Construction From Plane/Direction Pairs

The new plane/direction constructor builds a crystal-to-specimen orientation from one crystal plane and one crystal direction. By default:

- the crystal plane normal maps to specimen `Z`
- the in-plane crystal direction maps to specimen `X`

If the supplied crystal direction is not exactly in the plane, PyTex projects it into the plane before constructing the right-handed basis. The same projection rule is used for custom specimen reference directions.

```python
import numpy as np

from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)

orientations = OrientationSet.from_plane_direction(
    plane=np.array(
        [
            [0, 0, 1],
            [1, 1, 1],
        ]
    ),
    direction=np.array(
        [
            [1, 0, 0],
            [1, -1, 0],
        ]
    ),
    specimen_frame=specimen,
    phase=phase,
)
```

For theory, singularities, and the exact basis-construction rule, see [../../tex/algorithms/orientation_representations_and_plane_direction_construction.tex](../../tex/algorithms/orientation_representations_and_plane_direction_construction.tex).

## Normalizing A KikuchiPy-Like EBSD Surface

`normalize_ebsd(...)` accepts either a KikuchiPy signal carrying `xmap` or a crystal-map-like object with quaternion rotations and map coordinates. The returned object is always a `NormalizedEBSDDataset`, so downstream texture code stays inside PyTex types.

```python
class XMapStub:
    def __init__(self, quaternions):
        self.rotations = quaternions
        self.x = np.array([0.0, 1.0])
        self.y = np.array([0.0, 0.0])
        self.shape = (1, 2)
        self.dx = 1.0
        self.dy = 1.0
        self.source_file = "indexed.h5"


class SignalStub:
    def __init__(self, xmap):
        self.xmap = xmap


dataset = normalize_ebsd(
    SignalStub(XMapStub(orientations.quaternions)),
    frames=(crystal, specimen, specimen),
    phase=phase,
)
```

If live KikuchiPy indexing or refinement is available, `index_hough(...)` and `refine_orientations(...)` wrap those calls and return a `KikuchiPyWorkflowResult` containing:

- the normalized PyTex dataset
- an `ExperimentManifest`
- the returned `xmap`
- optional index or band arrays from Hough indexing

## Producing Texture Outputs From A Crystal Map

Once the EBSD surface is normalized, `CrystalMap` now exposes texture outputs directly.

```python
crystal_map = dataset.crystal_map

odf = crystal_map.to_odf()
pole_figure = crystal_map.pole_figure([1, 0, 0])
inverse_pole_figure = crystal_map.inverse_pole_figure("z")

report = crystal_map.texture_report(
    poles=([1, 0, 0], [1, 1, 1]),
    sample_directions=("x", "z"),
    plot=True,
)
```

`texture_report(...)` is the shortest route when you want one ODF plus one or more PF/IPF outputs in a single call. The returned `TextureReport` keeps the objects grouped without collapsing them into an ad hoc dictionary.

## Plotting The Map And Texture Surfaces

`plot_ipf_map(...)` colors the EBSD map from the current `IPFColorKey` implementation, while `plot_kam_map(...)` uses the existing KAM calculation and renders the result as a scalar image.

```python
from pytex import plot_ipf_map, plot_kam_map

ipf_map_figure = plot_ipf_map(crystal_map, direction="z")
kam_figure = plot_kam_map(crystal_map, order=1, symmetry_aware=True)
```

Both plotting calls accept optional grain-boundary overlays through either a `GrainSegmentation` or a `GrainBoundaryNetwork`.

## Convention Notes

- `Orientation` and `OrientationSet` remain crystal-to-specimen mappings.
- Quaternion order remains `w, x, y, z`.
- `RotationSet.to_axes_angles()` returns angles in radians.
- `RotationSet.to_rodrigues(frank=False)` uses the three-component Rodrigues vector `n tan(omega / 2)`.
- `RotationSet.to_rodrigues(frank=True)` returns Rodrigues--Frank coordinates `(n_x, n_y, n_z, tan(omega / 2))`.
- Exact `180°` rotations are singular in Rodrigues space; PyTex keeps that singularity explicit in the Rodrigues--Frank scale term.

## Related Material

- {doc}`ebsd_import_normalization`
- {doc}`ebsd_kam`
- {doc}`texture_odf_inversion`
- {doc}`../tutorials/notebooks`
- [../../tex/algorithms/orientation_representations_and_plane_direction_construction.tex](../../tex/algorithms/orientation_representations_and_plane_direction_construction.tex)

## References

### Normative

- `docs/standards/notation_and_conventions.md`
- `docs/testing/mtex_parity_matrix.md`
- `docs/testing/strategy.md`

### Informative

- Bunge, *Texture Analysis in Materials Science: Mathematical Methods*.
- Rowenhorst et al., *Modelling and Simulation in Materials Science and Engineering* (2015), DOI: <https://doi.org/10.1088/0965-0393/23/8/083501>.
- KikuchiPy documentation: [Indexing patterns and orientations](https://kikuchipy.org/en/stable/user/guide/indexing.html).
