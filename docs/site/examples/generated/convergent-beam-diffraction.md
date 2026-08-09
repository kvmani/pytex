<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Convergent-beam diffraction

The absolute scale of the two-beam extinction distance, checked against a published table, and the fringe analysis that measures a foil thickness without needing one.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Extinction distances of aluminium at 100 kV match the published table

Every dynamical quantity in a CBED pattern is measured in units of the extinction distance, so its absolute scale has to be right - and a wrong scale is invisible in the geometry, showing only as fringes at the wrong spacing. Two things set that scale: the Mott-Bethe conversion of the X-ray form factor into an electron scattering factor in angstrom, and the relativistic factor gamma = 1 + E/m0c^2, which is 1.20 at 100 kV. Aluminium is the calibration case because the fitted scattering-factor parametrization is most accurate for light elements; for heavy elements the same calculation is only indicative, which is why CBED practice measures the extinction distance rather than tabulating it.

**Symbols**

- $\xi_{g}$ &mdash; Two-beam extinction distance of reflection g; the depth period of the intensity exchange between the transmitted and diffracted beams.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    AtomicSite,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SpaceGroupSpec,
    SymmetrySpec,
    UnitCell,
    extinction_distance_angstrom,
    thickness_from_fringe_minima,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
aluminium_lattice = Lattice(4.0495, 4.0495, 4.0495, 90.0, 90.0, 90.0, crystal_frame=crystal)
aluminium = Phase(
    name="aluminium-fcc",
    lattice=aluminium_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=aluminium_lattice,
        sites=tuple(
            AtomicSite(
                label=f"Al{index}",
                species="Al",
                fractional_coordinates=np.asarray(position, dtype=float),
            )
            for index, position in enumerate(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
            )
        ),
    ),
    space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=crystal),
)
```

:::

**Compute**

```python
result = extinction_distance_angstrom(
    aluminium, [(1, 1, 1), (2, 0, 0), (2, 2, 0)], beam_energy_kev=100.0
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-cbed-aluminium-extinction-distances-at-100kv` | [555.2, 663.9, 1062.5] | [556.0, 673.0, 1057.0] | angstrom | 9.13e+00 | 1e+01 | ✅ pass |

**Why this value**: Published two-beam extinction distances for aluminium at 100 kV, 556, 673 and 1057 angstrom for {111}, {200} and {220}. The tolerance of 10 angstrom is about 1.5 percent, the accuracy the fitted scattering-factor parametrization supports for a light element.

**Citation**: Williams and Carter, Transmission Electron Microscopy, 2nd ed. (Springer, 2009), Table 23.1.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`, {doc}`Diffraction API <../../api/index>`

## The Kelly plot recovers both the thickness and the extinction distance from fringe positions

Reading the dark fringes off a single CBED disc gives the local foil thickness - and, from the same straight-line fit, the extinction distance, so the thickness does not inherit the error of a tabulated constant. Here the fringe positions are generated from the two-beam relation for a chosen thickness of 2000 angstrom and extinction distance of 500 angstrom, and the fit is asked to recover both. Because the input is the closed-form relation rather than a simulation, this tests the inversion itself and nothing else.

**Symbols**

- $t$ &mdash; Foil thickness along the beam.
- $\xi_{g}$ &mdash; Two-beam extinction distance of reflection g; the depth period of the intensity exchange between the transmitted and diffracted beams.
- $s_{g}$ &mdash; Excitation error: deviation of reflection g from the exact Bragg condition.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    AtomicSite,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SpaceGroupSpec,
    SymmetrySpec,
    UnitCell,
    extinction_distance_angstrom,
    thickness_from_fringe_minima,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
aluminium_lattice = Lattice(4.0495, 4.0495, 4.0495, 90.0, 90.0, 90.0, crystal_frame=crystal)
aluminium = Phase(
    name="aluminium-fcc",
    lattice=aluminium_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=aluminium_lattice,
        sites=tuple(
            AtomicSite(
                label=f"Al{index}",
                species="Al",
                fractional_coordinates=np.asarray(position, dtype=float),
            )
            for index, position in enumerate(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
            )
        ),
    ),
    space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=crystal),
)
```

:::

**Compute**

```python
thickness, extinction = 2000.0, 500.0
orders = np.arange(5, 11, dtype=float)
minima = np.sqrt((orders / thickness) ** 2 - extinction**-2)
report = thickness_from_fringe_minima(minima, first_order=5)
result = np.array([
    report.thickness_angstrom,
    report.extinction_distance_angstrom,
])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-cbed-two-beam-thickness-inverts-the-fringe-relation` | [2000.000000, 500.000000] | [2000.000000, 500.000000] | angstrom | 4.55e-13 | 1e-06 | ✅ pass |

**Why this value**: An exact inversion of the two-beam minimum condition t s_eff,n = n with s_eff^2 = s^2 + xi^-2, which rearranges to (s_n/n)^2 = 1/t^2 - (1/xi^2)(1/n^2). The generated minima lie on that line by construction, so a least-squares fit returns the intercept 1/t^2 and the slope -1/xi^2 to machine precision.

**Citation**: Kelly, Jostsons, Blake and Napier, Physica Status Solidi (a) 31 (1975) 771-780, for the linearization; Williams and Carter, Transmission Electron Microscopy, 2nd ed. (2009), Chapter 23.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`, {doc}`Diffraction API <../../api/index>`
