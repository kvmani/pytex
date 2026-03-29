# XRDML Texture Import And ODF Reconstruction

PyTex now includes a stable XRDML import boundary for pole-figure-style texture data.

Current entry points:

- `XRDMLPoleFigureMeasurement`
- `read_xrdml_pole_figure(...)`
- `load_xrdml_pole_figure(...)`
- `invert_xrdml_pole_figures(...)`
- `OrientationSet.from_bunge_grid(...)`
- `ODF.evaluate_pole_density(...)`

The adapter lives under `pytex.adapters` because XRDML is a vendor XML format. The public output of the adapter is still a PyTex `PoleFigure`, not a raw vendor payload.

## Important Semantics

- PyTex reads `Phi` and `Psi` positions from PANalytical-style pole-figure XRDML scans.
- The imported measurement is converted into explicit specimen directions on the unit sphere.
- The crystallographic pole is **not** inferred from the file. Callers must provide the corresponding `CrystalPlane` explicitly when converting imported data into a `PoleFigure`.
- Current ODF reconstruction remains the existing discrete kernel inversion model. PyTex does not yet claim full harmonic MTEX-class XRD texture inversion parity.

## Import Example

```python
from pytex import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    load_xrdml_pole_figure,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)

pole = CrystalPlane(MillerIndex([1, 1, 3], phase=phase), phase=phase)
pole_figure = load_xrdml_pole_figure(
    "fixtures/xrdml/polefig_Ge113_xrayutilities.xrdml.bz2",
    pole=pole,
    specimen_frame=specimen,
)
```

## Reconstruction Example

```python
from pytex import KernelSpec, OrientationSet, invert_xrdml_pole_figures

dictionary = OrientationSet.from_bunge_grid(
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
    phase=phase,
    phi1_step_deg=15.0,
    big_phi_step_deg=15.0,
    phi2_step_deg=15.0,
)

report = invert_xrdml_pole_figures(
    [
        "measurements/pf_100.xrdml",
        "measurements/pf_111.xrdml",
    ],
    poles=[
        CrystalPlane(MillerIndex([1, 0, 0], phase=phase), phase=phase),
        CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase),
    ],
    specimen_frame=specimen,
    orientation_dictionary=dictionary,
    kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=10.0),
    include_symmetry_family=False,
)
```

## Synthetic Benchmarking Support

`ODF.evaluate_pole_density(...)` exists so imported or synthetic measurement grids can be compared against the current discrete ODF model without detouring through plotting code. This is used internally for the synthetic XRDML inversion benchmark and is also useful for teaching or debugging the current inversion doctrine.

## Current Limits

- the implemented parser targets pole-figure-style `Phi` and `Psi` scans, not every XRDML measurement family
- PyTex does not currently infer reflection metadata from free-text file comments
- the reconstruction path is dictionary-based and kernel-based, not harmonic
- intensity normalization and defocusing corrections remain the caller's responsibility unless encoded in the supplied data

## Related Material

- {doc}`texture_odf_inversion`
- {doc}`stereographic_projections`
- `docs/testing/mtex_parity_matrix.md`

## References

### Normative

- `docs/standards/data_contracts_and_manifests.md`
- `docs/testing/mtex_parity_matrix.md`
- `benchmarks/texture/xrdml_import_benchmark_manifest.json`

### Informative

- MTEX documentation: [Import Pole Figure Data](https://mtex-toolbox.github.io/PoleFigureImport.html)
- MTEX documentation: [loadPoleFigure_xrdml](https://mtex-toolbox.github.io/loadPoleFigure_xrdml.html)
- MTEX documentation: [ODF Reconstruction from X-Ray Diffraction Data of an Al Alloy Rolled Sheet](https://mtex-toolbox.github.io/ExAlODFReconstruction.html)
- xrayutilities repository: [dkriegner/xrayutilities](https://github.com/dkriegner/xrayutilities)
