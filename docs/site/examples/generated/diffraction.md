<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Diffraction geometry

Powder scattering angles derived from PyTex interplanar spacings via Bragg's law, checked against a standard reference reflection position.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Ni(111) powder reflection angle for Cu K-alpha1

You are calibrating or interpreting a powder pattern and need to predict where the Ni(111) peak should appear with a copper source. PyTex supplies the interplanar spacing from the lattice metric; Bragg's law then gives the scattering angle. The result should land on the textbook Ni(111) position near 44.5 degrees for Cu K-alpha1.

**Symbols**

- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.
- $\lambda$ &mdash; Radiation wavelength.
- $\theta$ &mdash; Bragg half-angle.
- $2\theta$ &mdash; Powder-diffraction scattering angle.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerPlane,
    Phase,
    RadiationSpec,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
cu_ka1 = RadiationSpec.cu_ka().wavelength_angstrom
```

:::

**Compute**

```python
d_111 = MillerPlane.from_hkl([1, 1, 1], phase=nickel).d_spacing_angstrom
theta = np.arcsin(cu_ka1 / (2.0 * d_111))
result = float(np.degrees(2.0 * theta))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-ni-111-two-theta` | 44.496 | 44.496 | deg | 1.43e-04 | 5e-03 | ✅ pass |

**Why this value**: d_111 = 3.52387 / sqrt(3) = 2.03451 angstrom; with lambda = 1.5406 angstrom, 2*theta = 2*arcsin(lambda / (2 d)) = 44.50 degrees, matching standard Ni powder data.

**Citation**: ICDD PDF 04-0850 (nickel); Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed.

**See also**: {doc}`Powder XRD generation <../../workflows/xrd_generation>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`
