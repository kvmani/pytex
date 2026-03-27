# Powder XRD Generation

PyTex now includes a structure-aware powder XRD workflow built on the same lattice, phase, and reciprocal-space semantics used elsewhere in the library.

![Powder XRD Example](../../figures/powder_xrd_demo.svg)

## Scope

- configurable wavelength through `RadiationSpec`
- reflection enumeration from the canonical lattice model
- `d`-spacing and `2\theta` computation from Bragg's law
- approximate intensity estimation from crystal structure and multiplicity
- optional Gaussian broadening into a continuous powder spectrum
- runtime plotting through the shared YAML style system

## Scientific Model

For a reflection with spacing \(d_{hkl}\) and wavelength \(\lambda\), PyTex applies Bragg's law

$$
2 d_{hkl} \sin \theta = \lambda,
$$

then reports the observable angle \(2\theta\). The current implementation computes

- the reciprocal-lattice vector magnitude \(||\mathbf{g}_{hkl}||\)
- \(d_{hkl} = 1 / ||\mathbf{g}_{hkl}||\)
- \(2\theta = 2 \arcsin(\lambda / 2d_{hkl})\)

The current intensity model is intentionally modest. It uses:

- multiplicity inferred from the phase point-group symmetry
- a simple x-ray structure-factor proxy based on atomic numbers and fractional coordinates
- a Lorentz-polarization factor

This is a useful foundational spectrum model, but it is not yet a fully calibrated scattering-factor implementation.

## Example

```python
import numpy as np

from pytex import (
    AtomicSite,
    Lattice,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
    UnitCell,
    generate_xrd_pattern,
    plot_xrd_pattern,
)
from pytex.core.conventions import FrameDomain, Handedness

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
lattice = Lattice(3.615, 3.615, 3.615, 90.0, 90.0, 90.0, crystal_frame=crystal)
unit_cell = UnitCell(
    lattice=lattice,
    sites=(
        AtomicSite("Cu1", "Cu", np.array([0.0, 0.0, 0.0])),
        AtomicSite("Cu2", "Cu", np.array([0.0, 0.5, 0.5])),
        AtomicSite("Cu3", "Cu", np.array([0.5, 0.0, 0.5])),
        AtomicSite("Cu4", "Cu", np.array([0.5, 0.5, 0.0])),
    ),
)
phase = Phase(
    "Cu",
    lattice=lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=unit_cell,
)

pattern = generate_xrd_pattern(
    phase,
    radiation=RadiationSpec.cu_ka(),
    two_theta_range_deg=(20.0, 120.0),
    resolution_deg=0.02,
    max_index=6,
    broadening_fwhm_deg=0.18,
)
figure = plot_xrd_pattern(pattern, theme="journal")
figure.savefig("cu_powder_xrd.png", dpi=200)
```

## Interpretation Notes

- `PowderReflection` is the reflection-level object carrying `d` spacing, `2\theta`, multiplicity, and intensity metadata.
- `PowderPattern` is the broadened spectrum object carrying the reflection list plus grid-sampled intensity.
- The current intensity surface is suitable for teaching, method prototyping, and structure-sensitive inspection, but it is not yet a full Rietveld-grade scattering engine.

## Current Limits

- no tabulated atomic form factors yet
- no preferred-orientation, absorption, or instrument-function model yet
- no profile families beyond the current Gaussian broadening path

## Related Material

- {doc}`../concepts/technical_glossary_and_symbols`
- {doc}`saed_generation`
- {doc}`style_customization`
- [../../tex/algorithms/powder_xrd_and_saed.tex](../../tex/algorithms/powder_xrd_and_saed.tex)

## References

### Normative

- `../../standards/reference_canon.md`
- `../../standards/notation_and_conventions.md`

### Informative

- `../../testing/diffraction_validation_matrix.md`
