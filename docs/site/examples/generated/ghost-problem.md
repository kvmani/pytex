<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# The ghost problem

What diffraction pole figures cannot determine: an asymmetric texture still gives a pole set closed under negation, and excluding the odd harmonic degrees that centrosymmetry annihilates discards nearly half the basis.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## An asymmetric texture still gives a centrosymmetric pole figure

Build a deliberately one-sided orientation population - 200 orientations whose first Euler angle is confined to 0-40 degrees, with nothing symmetric about it - and generate its {111} pole figure. Every pole in the result has its antipode also present. The centrosymmetry is a property of the measurement, not of the specimen: a lattice plane has no sense and Friedel's law makes +h and -h scatter identically. This is the root cause of the ghost problem, and the reason no amount of pole-figure data can recover the odd part of an ODF. The example returns the fraction of the first 300 poles whose antipode is in the set, which must be 1.

**Symbols**

- $P_{\mathbf{h}}(\mathbf{y})$ &mdash; Pole density of plane family h along specimen direction y.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.texture.models import PoleFigure

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
phase = Phase(
    name="nickel",
    lattice=Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
```

:::

**Compute**

```python
angles = np.linspace(0.0, 40.0, 200)
orientations = OrientationSet.from_euler_angles(
    np.column_stack(
        [angles, 25.0 * np.ones_like(angles), np.zeros_like(angles)]
    ),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)
figure = PoleFigure.from_orientations(orientations, pole)
directions = figure.sample_directions
poles = np.asarray(
    getattr(directions, 'values', directions)
)
paired = sum(
    1
    for row in poles[:300]
    if np.min(np.linalg.norm(poles + row, axis=1)) < 1e-9
)
result = paired / 300.0
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ghost-pole-figure-is-centrosymmetric` | 1.000000 | 1.000000 | &mdash; | exact | exact | ✅ pass |

**Why this value**: Analytic: a plane normal enters a pole figure as an axis and Friedel's law equates +h with -h, so the pole set is closed under negation for any orientation distribution whatsoever.

**Citation**: Matthies, On the reproducibility of the orientation distribution function of texture samples from pole figures (ghost phenomena), Phys. Status Solidi B 92 (1979) K135-K138.

**See also**: {doc}`The ghost problem <../../theory/ghost_problem_and_odd_harmonics>`, {doc}`Harmonic ODF reconstruction <../../theory/harmonic_odf_reconstruction>`

## Excluding odd degrees discards nearly half the harmonic basis

Count the generalized spherical harmonic terms retained and discarded when odd degrees are excluded, at bandlimit 22. Because pole figures annihilate every odd degree, those coefficients are not poorly determined but wholly undetermined, and the count says how much of the ODF a diffraction measurement is silent about: 15147 of 32407 terms, or 46.7 percent, tending to one half as the bandlimit grows. The example returns the discarded fraction.

**Symbols**

- $\ell$ &mdash; Degree of a generalized spherical harmonic term.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.texture.harmonics import _enumerate_terms
```

:::

**Compute**

```python
even = len(
    _enumerate_terms(degree_bandlimit=22, even_degrees_only=True)
)
every = len(
    _enumerate_terms(degree_bandlimit=22, even_degrees_only=False)
)
result = (every - even) / every
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ghost-odd-degrees-are-half-the-harmonic-basis` | 0.46740 | 0.46739 | &mdash; | 9.02e-06 | 1e-04 | ✅ pass |

**Why this value**: Exact term count: sum over degrees 0..22 of (2l+1)^2 is 32407, of which the even degrees contribute 17260, leaving 15147 discarded. The fraction tends to 1/2 as the bandlimit grows.

**Citation**: Bunge, Texture Analysis in Materials Science: Mathematical Methods (Butterworths 1969) - the generalized spherical harmonic expansion.

**See also**: {doc}`The ghost problem <../../theory/ghost_problem_and_odd_harmonics>`, {doc}`Harmonic ODF reconstruction <../../theory/harmonic_odf_reconstruction>`
