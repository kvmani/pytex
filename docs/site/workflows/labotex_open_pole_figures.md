# LaboTex Pole-Figure Import And Open Multi-Material Validation

PyTex now supports direct import of LaboTex `PPF` and `EPF` pole-figure files through the adapter boundary.

Current entry points:

- `LaboTexPoleFigureDescriptor`
- `LaboTexPoleFigureMeasurement`
- `read_labotex_pole_figures(...)`
- `load_labotex_pole_figures(...)`
- `invert_labotex_pole_figures(...)`

This surface is useful for open measured texture data that is distributed in legacy pole-figure formats rather than modern vendor XML. The adapter normalizes those files onto the same `PoleFigure` and dictionary-based `ODF` surfaces used by the rest of PyTex.

## Current Scope

- multi-pole-figure `PPF` and `EPF` files with explicit `alpha` / `beta` scan geometry
- import of the crystallographic pole family encoded inside the file header
- conversion into PyTex `PoleFigure` objects on explicit specimen directions
- direct handoff into `ODF.invert_pole_figures(...)` through `invert_labotex_pole_figures(...)`

`POW` defocussing-correction files are not yet a stable import surface.

## Format Semantics

- `alpha` is treated as the pole-figure tilt angle from the specimen normal
- `beta` is treated as the azimuthal angle around the specimen normal
- each file may contain multiple pole figures, each with its own `hkl` indices and 2-theta value
- intensity normalization remains explicit and caller-controlled, just as in the XRDML path

## Example

```python
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    load_labotex_pole_figures,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
lattice = Lattice(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
phase = Phase("copper_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)

pole_figures = load_labotex_pole_figures(
    "Exercise2-Cu.PPF",
    phase=phase,
    specimen_frame=specimen,
    intensity_normalization="max",
)
```

## Open Validation Program

PyTex includes a helper script for a stronger open-data validation pass:

```bash
python scripts/run_open_labotex_validation.py
```

That script downloads the public LaboTex exercise archive into the ignored `inspection_outputs/`
tree, then runs measured-vs-reconstructed pole-figure comparisons for:

- Cu (`Exercise2-Cu.PPF`)
- Fe (`Exercise4-Fe.epf`)
- Al (`Exercise5-Al.epf`)

For each dataset it generates:

- measured pole-figure contour plots
- reconstructed pole-figure contour plots on the same sampled directions
- absolute-difference pole-figure contour plots
- ODF section plots
- a JSON summary with inversion diagnostics

These files are intentionally not bundled in-repo because the source site exposes them as public
downloads but does not state redistribution terms explicitly.

## Interpretation Notes

- The imported LaboTex datasets are useful for external validation of the software path, not as a blanket claim of full parity with LaboTex or MTEX.
- Reconstruction quality depends strongly on the number of measured poles and on the orientation dictionary used for inversion.
- A low residual is useful evidence that the import geometry, pole indexing, and inversion path are coherent. It is not by itself a proof of full physical texture fidelity.

## Related Material

- {doc}`texture_odf_inversion`
- {doc}`xrdml_texture_import`
- `../../testing/mtex_parity_matrix.md`

## References

### Normative

- `../../standards/data_contracts_and_manifests.md`
- `../../testing/mtex_parity_matrix.md`

### Informative

- LaboTex data-format overview: [Supported data formats](https://labosoft.com.pl/formats.html)
- LaboTex exercises page: [Download](https://labosoft.com.pl/download.html)
