<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Elastic anisotropy and homogenization

Cubic elasticity against closed forms: [110] and [112] are exactly equally stiff, the Voigt and Reuss bulk moduli of a cubic aggregate are identical so the entire bound gap is in the shear modulus, and a numerically homogenized random aggregate reproduces both shear bounds.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## [110] and [112] are exactly equally stiff in any cubic crystal

Evaluate Young's modulus of copper along [100], [110], [112] and [111]. For cubic symmetry the direction enters only through J = n1^2 n2^2 + n2^2 n3^2 + n3^2 n1^2, which is 0 along <100>, 1/3 along <111>, and exactly 1/4 along both [110] and [112]. Those two directions therefore have identical stiffness - not approximately, and not only for copper. The example returns all four moduli so the extremes and the coincidence are checked together.

**Symbols**

- $E(\hat{\mathbf{n}})$ &mdash; Young's modulus along the crystal direction n.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.properties import StiffnessTensor

# Copper single-crystal stiffness, GPa.
C11, C12, C44 = 168.4, 121.4, 75.4
stiffness = StiffnessTensor.cubic(C11, C12, C44)
compliance = stiffness.compliance()
voigt = np.asarray(compliance.voigt_matrix())
S11, S12, S44 = voigt[0, 0], voigt[0, 1], voigt[3, 3]


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)
```

:::

**Compute**

```python
directions = [
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 2.0],
    [1.0, 1.0, 1.0],
]
result = np.array(
    [float(compliance.youngs_modulus(unit(d))) for d in directions]
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `elastic-cubic-youngs-modulus-110-equals-112` | [66.6888, 130.3376, 130.3376, 191.1497] | [66.6888, 130.3376, 130.3376, 191.1497] | GPa | 4.91e-05 | 1e-03 | ✅ pass |

**Why this value**: Closed form 1/E = S11 - 2(S11 - S12 - S44/2) J with J = 0, 1/4, 1/4, 1/3 respectively; S11, S12, S44 from the analytic cubic inverse of (C11, C12, C44).

**Citation**: Nye, Physical Properties of Crystals (OUP); Simmons and Wang, Single Crystal Elastic Constants (MIT Press 1971) for the copper constants.

**See also**: {doc}`Elastic anisotropy and homogenization <../../theory/elastic_anisotropy_and_homogenization>`

## The Voigt and Reuss bulk moduli of a cubic aggregate are equal

Compute the Voigt and Reuss bulk moduli of a randomly textured cubic aggregate from the closed forms and take their difference. It is exactly zero: a cubic crystal responds isotropically to hydrostatic pressure, so the uniform-stress and uniform-strain assumptions cannot disagree about dilatation. The whole Voigt-Reuss gap therefore lives in the shear modulus, and reporting a 'Hill bulk modulus' for a cubic aggregate implies an uncertainty that does not exist. The example returns the difference K_V - K_R, which must vanish.

**Symbols**

- $K_{V}, K_{R}$ &mdash; Voigt and Reuss bulk moduli of an aggregate.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.properties import StiffnessTensor

# Copper single-crystal stiffness, GPa.
C11, C12, C44 = 168.4, 121.4, 75.4
stiffness = StiffnessTensor.cubic(C11, C12, C44)
compliance = stiffness.compliance()
voigt = np.asarray(compliance.voigt_matrix())
S11, S12, S44 = voigt[0, 0], voigt[0, 1], voigt[3, 3]


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)
```

:::

**Compute**

```python
k_voigt = (C11 + 2.0 * C12) / 3.0
k_reuss = 1.0 / (3.0 * (S11 + 2.0 * S12))
result = float(k_voigt - k_reuss)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `elastic-cubic-voigt-reuss-bulk-moduli-coincide` | -8.527e-14 | 0.000e+00 | GPa | 8.53e-14 | 1e-09 | ✅ pass |

**Why this value**: Analytic identity: K_V = (C11 + 2 C12)/3 and K_R = 1/(3(S11 + 2 S12)) are equal for cubic symmetry, since S11 + 2 S12 = 3/(C11 + 2 C12).

**Citation**: Hill, The elastic behaviour of a crystalline aggregate, Proc. Phys. Soc. A 65 (1952) 349-354, DOI 10.1088/0370-1298/65/5/307.

**See also**: {doc}`Elastic anisotropy and homogenization <../../theory/elastic_anisotropy_and_homogenization>`

## A homogenized random aggregate reproduces the Voigt and Reuss shear moduli

Homogenize copper over 40000 Haar-random orientations under the Voigt and Reuss schemes and read the aggregate shear modulus back as C44 of the averaged tensor, comparing with the closed forms mu_V = (C11 - C12 + 3 C44)/5 and mu_R = 5/(4(S11 - S12) + 3 S44). This exercises the rank-four rotation, the weighted average, and the compliance inversion together. The residual is finite-sample texture in the random orientation set, not an error in the averaging.

**Symbols**

- $\mu_{V}, \mu_{R}$ &mdash; Voigt and Reuss aggregate shear moduli.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.properties import StiffnessTensor

# Copper single-crystal stiffness, GPa.
C11, C12, C44 = 168.4, 121.4, 75.4
stiffness = StiffnessTensor.cubic(C11, C12, C44)
compliance = stiffness.compliance()
voigt = np.asarray(compliance.voigt_matrix())
S11, S12, S44 = voigt[0, 0], voigt[0, 1], voigt[3, 3]


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)

from pytex.core import (
    FrameDomain,
    OrientationSet,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.properties import homogenize_elastic

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)

# Haar-uniform orientations (Shoemake), so the aggregate is texture-free.
rng = np.random.default_rng(17)
u1, u2, u3 = rng.random(40000), rng.random(40000), rng.random(40000)
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
moduli = []
for scheme in ('voigt', 'reuss'):
    aggregate = homogenize_elastic(
        stiffness, orientations, scheme=scheme
    )
    moduli.append(float(np.asarray(aggregate.voigt_matrix())[3, 3]))
result = np.array(moduli)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `elastic-random-aggregate-matches-voigt-reuss-closed-form` | [54.6314, 40.0264] | [54.6400, 40.0339] | GPa | 8.56e-03 | 5e-02 | ✅ pass |

**Why this value**: Closed forms for a randomly textured cubic aggregate: mu_V = (C11 - C12 + 3 C44)/5 = 54.6400 GPa and mu_R = 5/(4(S11 - S12) + 3 S44) = 40.0339 GPa. Tolerance covers the finite-sample texture of 40000 random orientations.

**Citation**: Hill, Proc. Phys. Soc. A 65 (1952) 349-354; Simmons and Wang, Single Crystal Elastic Constants (MIT Press 1971).

**See also**: {doc}`Elastic anisotropy and homogenization <../../theory/elastic_anisotropy_and_homogenization>`
