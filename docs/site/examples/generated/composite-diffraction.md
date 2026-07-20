<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Composite OR diffraction

Numerical cornerstones of composite orientation-relationship SAED simulation: the relativistic electron wavelength against the standard 200 kV value, and the exactness of the Kurdjumov-Sachs child-zone mapping.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Relativistic electron wavelength at 200 kV

Every kinematic TEM computation starts from the electron wavelength, which fixes the Ewald-sphere radius k = 1/lambda and hence every excitation error. The relativistic formula lambda = h / sqrt(2 m0 e V (1 + e V / (2 m0 c^2))) must reproduce the standard tabulated value at a 200 kV accelerating voltage.

**Symbols**

- $\lambda$ &mdash; Radiation wavelength.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.diffraction.kinematic import electron_wavelength_angstrom
```

:::

**Compute**

```python
result = electron_wavelength_angstrom(200.0)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `composite-electron-wavelength-200kv` | 0.02508 | 0.02508 | angstrom | 6.60e-07 | 5e-06 | ✅ pass |

**Why this value**: The standard relativistic electron wavelength at 200 kV is 2.508 pm = 0.02508 angstrom.

**Citation**: De Graef, Introduction to Conventional Transmission Electron Microscopy, Cambridge University Press, 2003, Table 2.2.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## KS maps the parent [0 1 -1] zone exactly onto a <1 1 1> child zone

The Kurdjumov-Sachs relationship is defined by the parallelism <-1 0 1>_fcc || <-1 -1 1>_bcc. When the composite SAED simulator maps a parent [0 1 -1] zone axis (a member of the <-1 0 1> family) through all 24 variants, at least one variant's child zone axis must land exactly on a <1 1 1>-type direction: the minimal angular deviation between mapped and rational child zones over the variants is zero.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import simulate_composite_saed

parent_frame = ReferenceFrame(
    "parent_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
child_frame = ReferenceFrame(
    "child_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
austenite = Phase(
    "austenite",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
    crystal_frame=parent_frame,
)
martensite = Phase(
    "martensite",
    lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
    crystal_frame=child_frame,
)
ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=martensite
)
```

:::

**Compute**

```python
zone = ZoneAxis(np.array([0, 1, -1]), phase=austenite)
composite = simulate_composite_saed(ks, zone, include_parent=False)
result = min(
    pattern.nearest_zone_axis.deviation_deg
    for pattern in composite.variant_patterns
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `composite-ks-exact-child-zone` | 0.0000 | 0.0000 | deg | 9.00e-15 | 1e-09 | ✅ pass |

**Why this value**: The defining KS direction parallelism makes the mapped child zone rational, so the deviation of the best variant is exactly 0 degrees.

**Citation**: Kurdjumov and Sachs, Z. Physik 64 (1930) 325; Morito et al., Acta Materialia 51 (2003) 1789 (variant conventions).

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`
