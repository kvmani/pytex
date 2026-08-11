<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Inverse-pole-figure colouring

What an IPF colour actually is, checked against hand-derived values: the sector corners colour to exact primaries, a direction on the [001]-[111] edge colours to exactly (1, 0, 3/4) by the closed form, and every symmetric equivalent shares one colour.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## The cubic sector corners colour to exactly red, green and blue

Colour the three corners of the cubic standard triangle - [001], [101] and [111]. Each corner has barycentric weights equal to a standard basis vector, so after the saturation power and the max-channel renormalization it must return exactly one primary. This is the identity that makes an IPF legend readable without consulting a lookup table, and it is what fails first if the sector corners or the colour basis are mis-ordered.

**Symbols**

- $\boldsymbol{\beta}$ &mdash; Barycentric weights of a direction in the fundamental-sector corner basis.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    IPFColorKey,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
key = IPFColorKey(crystal_symmetry=symmetry, specimen_direction="z")
```

:::

**Compute**

```python
corners = np.array(
    [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]
)
result = key.colors_from_crystal_directions(corners)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ipf-cubic-sector-corners-are-primaries` | [1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000] | [1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000] | &mdash; | 3.88e-30 | 1e-12 | ✅ pass |

**Why this value**: Analytic identity: a sector corner has barycentric weights equal to a standard basis vector, which the colour map carries to the corresponding primary for any saturation exponent.

**Citation**: Nolze and Hielscher, Orientations - perfectly colored, J. Appl. Cryst. 49 (2016) 1786-1802, DOI 10.1107/S1600576716012942.

**See also**: {doc}`IPF colour keys <../../theory/ipf_color_keys>`, {doc}`IPF colour workflow <../../workflows/ipf_colors>`

## The [113] direction colours to exactly (1, 0, 3/4)

Colour a direction lying on the [001]-[111] edge of the cubic triangle, where the whole chain can be followed by hand. For [113]/sqrt(11) the closed-form weights are beta = (2, 0, sqrt(3))/sqrt(11); the largest is beta_1, so the colour is (1, 0, sqrt(3)/2) raised to the power 1/gamma_s = 2, giving exactly (1, 0, 3/4). Agreement here exercises the symmetry reduction, the barycentric solve, the saturation power and the renormalization together, against a number derived rather than recorded.

**Symbols**

- $\boldsymbol{\beta}$ &mdash; Barycentric weights of a direction in the fundamental-sector corner basis.
- $\gamma_{s}$ &mdash; IPF saturation parameter; channels are raised to the power 1/gamma_s.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    IPFColorKey,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
key = IPFColorKey(crystal_symmetry=symmetry, specimen_direction="z")
```

:::

**Compute**

```python
direction = np.array([[1.0, 1.0, 3.0]])
direction = direction / np.linalg.norm(direction)
result = key.colors_from_crystal_directions(direction)[0]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ipf-cubic-closed-form-colour-113` | [1.000000, 0.000000, 0.750000] | [1.000000, 0.000000, 0.750000] | &mdash; | 9.45e-30 | 1e-12 | ✅ pass |

**Why this value**: Closed form: beta = (dz - dx, sqrt(2)(dx - dy), sqrt(3) dy) evaluated at (1,1,3)/sqrt(11) gives (2, 0, sqrt(3))/sqrt(11), and (beta / max beta) ** 2 = (1, 0, 3/4) exactly.

**Citation**: International Tables for Crystallography, Vol. A - the m-3m asymmetric unit that fixes the sector corners.

**See also**: {doc}`IPF colour keys <../../theory/ipf_color_keys>`, {doc}`IPF colour workflow <../../workflows/ipf_colors>`

## All 24 cubic equivalents of a direction take one colour

Generate every symmetric equivalent of a general direction under m-3m and colour them all. The spread across the orbit must be zero: symmetry-equivalent directions are the same physical direction, so a colouring that separated them would paint contrast where there is no crystallography, and the picture would depend on which equivalent index a file happened to store. The example reports the maximum channel spread over the orbit, which is bounded by the rotation arithmetic rather than by a tolerance in the colouring.

**Symbols**

- $\boldsymbol{\beta}$ &mdash; Barycentric weights of a direction in the fundamental-sector corner basis.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    IPFColorKey,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
key = IPFColorKey(crystal_symmetry=symmetry, specimen_direction="z")
```

:::

**Compute**

```python
direction = np.array([0.3, 0.1, 0.9])
direction = direction / np.linalg.norm(direction)
orbit = symmetry.equivalent_vectors(direction)
orbit = np.asarray(
    orbit.values if hasattr(orbit, 'values') else orbit
).reshape(-1, 3)
colors = key.colors_from_crystal_directions(orbit)
result = float(np.abs(colors - colors[0]).max())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ipf-symmetric-equivalents-share-one-colour` | 7.08e-15 | 0.00e+00 | &mdash; | 7.08e-15 | 1e-12 | ✅ pass |

**Why this value**: Analytic identity: the colour is a function of the symmetry-reduced direction alone, so it is constant on a symmetry orbit by construction. The expected spread is exactly zero.

**Citation**: Nolze and Hielscher, Orientations - perfectly colored, J. Appl. Cryst. 49 (2016) 1786-1802, DOI 10.1107/S1600576716012942.

**See also**: {doc}`IPF colour keys <../../theory/ipf_color_keys>`, {doc}`IPF colour workflow <../../workflows/ipf_colors>`
