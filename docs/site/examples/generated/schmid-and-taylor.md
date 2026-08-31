<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Schmid and Taylor plasticity factors

Slip geometry against exact answers: eight fcc systems share a Schmid factor of 1/sqrt(6) under [001] tension, the cube orientation's full-constraint Taylor factor is exactly sqrt(6), and a random fcc texture averages Taylor's 1938 value of 3.06.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Eight fcc systems share a Schmid factor of exactly 1/sqrt(6)

Resolve [001] tension onto the twelve {111}<110> systems of a cube-oriented fcc grain. The magnitudes take only two values: eight systems at exactly 1/sqrt(6) = 0.408248 and four at zero. The eightfold degeneracy is why a cube-oriented grain has no single preferred slip system, and it is the origin of the Taylor ambiguity - many different five-system combinations accommodate the same strain at the same cost. The example returns the largest magnitude.

**Symbols**

- $m$ &mdash; Schmid factor, cos(phi) cos(lambda); bounded by 1/2.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Lattice,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.properties import fcc_octahedral_slip, taylor_factors

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
phase = Phase(
    name="nickel",
    lattice=Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
family = fcc_octahedral_slip(phase)

cube = Orientation(
    rotation=Rotation.identity(),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)
```

:::

**Compute**

```python
factors = np.asarray(
    family.schmid_factors(cube, (0.0, 0.0, 1.0))
).ravel()
result = float(np.abs(factors).max())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `plasticity-fcc-cube-schmid-factor` | 0.408248 | 0.408248 | &mdash; | < 1e-12 | 1e-09 | ✅ pass |

**Why this value**: Analytic: for [001] tension the {111}<110> systems give m = (t.n)(t.d) = (1/sqrt(3))(1/sqrt(2)) = 1/sqrt(6) on the eight systems that are stressed at all.

**Citation**: Schmid and Boas, Kristallplastizitaet (1935); Kocks, Tome and Wenk, Texture and Anisotropy (CUP 1998).

**See also**: {doc}`Schmid factors and the Taylor factor <../../theory/schmid_and_taylor_plasticity>`

## The [001] fcc Taylor factor is exactly sqrt(6)

Solve the full-constraint Taylor problem for a cube-oriented fcc grain in uniaxial tension along [001]. Minimising total slip subject to matching the five independent components of the imposed deviatoric strain gives exactly sqrt(6) = 2.449490. Reproducing a closed form to six figures exercises the linear-programming formulation, the five-constraint reduction (the sixth component is implied because both sides are traceless), and the signed-slip doubling that lets non-negative variables represent either shear sense.

**Symbols**

- $M$ &mdash; Full-constraint Taylor factor: minimum total slip per unit equivalent strain.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Lattice,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.properties import fcc_octahedral_slip, taylor_factors

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
phase = Phase(
    name="nickel",
    lattice=Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
family = fcc_octahedral_slip(phase)

cube = Orientation(
    rotation=Rotation.identity(),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)
```

:::

**Compute**

```python
result = float(
    taylor_factors(family, cube, tension_axis=(0.0, 0.0, 1.0))
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `plasticity-fcc-cube-taylor-factor` | 2.449490 | 2.449490 | &mdash; | < 1e-12 | 1e-06 | ✅ pass |

**Why this value**: Analytic: the full-constraint Taylor factor of the cube orientation for {111}<110> slip in uniaxial tension is sqrt(6).

**Citation**: Taylor, Plastic strain in metals, J. Inst. Metals 62 (1938) 307-324; Bishop and Hill, Phil. Mag. 42 (1951) 414-427.

**See also**: {doc}`Schmid factors and the Taylor factor <../../theory/schmid_and_taylor_plasticity>`

## A random fcc texture has an average Taylor factor near 3.06

Average the full-constraint Taylor factor over 2000 Haar-random orientations of an fcc aggregate in uniaxial tension, recovering Taylor's 1938 result of about 3.06. This is the number that converts a single-crystal critical resolved shear stress into a polycrystal flow stress. The spread matters as much as the mean: over a random texture M runs from roughly 2.29 to 3.67, so the hardest orientation is about 60 percent harder than the softest, which is why a textured sheet can differ substantially in flow stress from a random one at the same composition and grain size.

**Symbols**

- $M$ &mdash; Full-constraint Taylor factor: minimum total slip per unit equivalent strain.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Lattice,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.properties import fcc_octahedral_slip, taylor_factors

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
phase = Phase(
    name="nickel",
    lattice=Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
family = fcc_octahedral_slip(phase)

cube = Orientation(
    rotation=Rotation.identity(),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)

rng = np.random.default_rng(5)
u1, u2, u3 = rng.random(2000), rng.random(2000), rng.random(2000)
quaternions = np.stack(
    [
        np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
        np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
        np.sqrt(u1) * np.sin(2 * np.pi * u3),
        np.sqrt(u1) * np.cos(2 * np.pi * u3),
    ],
    axis=-1,
)
orientations = OrientationSet.from_quaternions(
    quaternions,
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)
```

:::

**Compute**

```python
factors = np.asarray(
    taylor_factors(
        family, orientations, tension_axis=(0.0, 0.0, 1.0)
    )
)
result = float(factors.mean())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `plasticity-random-fcc-taylor-factor` | 3.0546 | 3.0600 | &mdash; | 5.42e-03 | 3e-02 | ✅ pass |

**Why this value**: Taylor's 1938 value for a randomly oriented fcc aggregate deforming by {111}<110> slip in uniaxial tension, M ~ 3.06. Tolerance covers the Monte-Carlo standard error at n = 2000 (about 0.009).

**Citation**: Taylor, Plastic strain in metals, J. Inst. Metals 62 (1938) 307-324.

**See also**: {doc}`Schmid factors and the Taylor factor <../../theory/schmid_and_taylor_plasticity>`
