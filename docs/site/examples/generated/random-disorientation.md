<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Random disorientation baseline

The null hypothesis every MDF claim is measured against: the exact (1 - cos w)/pi density with no symmetry, the exact cubic maximum 62.7994 degrees from the Rodrigues zone vertex, and the 2.2 percent of random cubic boundaries that are low-angle by chance.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## With no symmetry the mean disorientation angle is pi/2 + 2/pi

Sample disorientation angles between independent uniformly (Haar) distributed orientations with triclinic symmetry, so no reduction occurs, and compare the mean against the exact value. The Haar measure in axis-angle coordinates carries the factor (1 - cos w), which is the volume of the shell of rotations at angle w, so the density is (1 - cos w)/pi and its first moment integrates to pi/2 + 2/pi = 126.4756 degrees. This is the strongest available check on the sampler because the target is analytic rather than another simulation.

**Symbols**

- $\omega$ &mdash; Rotation (disorientation) angle of a misorientation.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import FrameDomain, ReferenceFrame, SymmetrySpec
from pytex.core.misorientation_distribution import (
    random_disorientation_angles_deg,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))


def sampled_angles(point_group, total, chunk=20000, seed0=100):
    # The reduction materialises n*|G|^2 matrices, so a large baseline is
    # generated in chunks and concatenated; the samples are independent.
    symmetry = SymmetrySpec.from_point_group(
        point_group, reference_frame=crystal
    )
    parts = [
        random_disorientation_angles_deg(symmetry, chunk, seed=seed0 + i)
        for i in range(total // chunk)
    ]
    return np.concatenate(parts)
```

:::

**Compute**

```python
angles = sampled_angles('1', 200000, seed0=200)
result = float(angles.mean())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `mdf-triclinic-mean-disorientation-angle` | 126.3581 | 126.4756 | deg | 1.18e-01 | 2e-01 | ✅ pass |

**Why this value**: Analytic: integral of w*(1 - cos w)/pi over [0, pi] equals pi/2 + 2/pi radians = 126.4756 degrees. Tolerance covers the Monte-Carlo standard error at n = 2e5.

**Citation**: Morawiec, Orientations and Rotations (Springer 2004) - the invariant measure on SO(3) in axis-angle coordinates.

**See also**: {doc}`Random disorientation baseline <../../theory/random_disorientation_baseline>`, {doc}`Orientation space and disorientation <../../theory/orientation_space_and_disorientation>`

## The largest cubic disorientation is 2*arctan(sqrt(23-16*sqrt2))

The maximum cubic disorientation is not a sampling outcome but a property of the Rodrigues fundamental zone, which for cubic symmetry is the cube |rho_i| <= sqrt(2)-1 intersected with the octahedron sum|rho_i| <= 1. The angle increases with |rho|, so the maximum sits at the vertex farthest from the origin, (sqrt2-1, sqrt2-1, 3-2sqrt2), which meets both constraints with equality. Its magnitude is sqrt(23 - 16 sqrt2), giving 62.7994 degrees about <1, 1, sqrt2-1>. This example evaluates that closed form; a sampled maximum converges to it only slowly from below, which is why the exact value is the one to quote.

**Symbols**

- $\boldsymbol{\rho}$ &mdash; Rodrigues vector n tan(omega/2); the chart in which the cubic fundamental zone is a cube intersected with an octahedron.
- $\omega$ &mdash; Rotation (disorientation) angle of a misorientation.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import FrameDomain, ReferenceFrame, SymmetrySpec
from pytex.core.misorientation_distribution import (
    random_disorientation_angles_deg,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))


def sampled_angles(point_group, total, chunk=20000, seed0=100):
    # The reduction materialises n*|G|^2 matrices, so a large baseline is
    # generated in chunks and concatenated; the samples are independent.
    symmetry = SymmetrySpec.from_point_group(
        point_group, reference_frame=crystal
    )
    parts = [
        random_disorientation_angles_deg(symmetry, chunk, seed=seed0 + i)
        for i in range(total // chunk)
    ]
    return np.concatenate(parts)
```

:::

**Compute**

```python
rho_max = np.array(
    [np.sqrt(2) - 1, np.sqrt(2) - 1, 3 - 2 * np.sqrt(2)]
)
result = float(
    np.degrees(2 * np.arctan(np.linalg.norm(rho_max)))
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `mdf-cubic-maximum-disorientation-angle` | 62.7994 | 62.7994 | deg | 2.96e-05 | 1e-03 | ✅ pass |

**Why this value**: Analytic: |rho| = sqrt(23 - 16*sqrt(2)) at the cubic fundamental-zone vertex, so omega = 2*arctan(sqrt(23 - 16*sqrt(2))) = 62.79943 degrees.

**Citation**: Mackenzie, Second paper on statistics associated with the random disorientation of cubes, Biometrika 45 (1958) 229-240, DOI 10.1093/biomet/45.1-2.229.

**See also**: {doc}`Random disorientation baseline <../../theory/random_disorientation_baseline>`, {doc}`Orientation space and disorientation <../../theory/orientation_space_and_disorientation>`

## 2.2 percent of random cubic boundaries are low-angle by chance

Compute the fraction of random cubic disorientations falling below the conventional 15 degree low-angle threshold. This is the null hypothesis a low-angle-boundary fraction has to beat: a texture-free aggregate already delivers about 2.2 percent low-angle boundaries from geometry alone, so a map reporting 3 percent has demonstrated essentially nothing while one reporting 30 percent has. The value is recomputed rather than tabulated because it depends on the point group.

**Symbols**

- $\omega$ &mdash; Rotation (disorientation) angle of a misorientation.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import FrameDomain, ReferenceFrame, SymmetrySpec
from pytex.core.misorientation_distribution import (
    random_disorientation_angles_deg,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))


def sampled_angles(point_group, total, chunk=20000, seed0=100):
    # The reduction materialises n*|G|^2 matrices, so a large baseline is
    # generated in chunks and concatenated; the samples are independent.
    symmetry = SymmetrySpec.from_point_group(
        point_group, reference_frame=crystal
    )
    parts = [
        random_disorientation_angles_deg(symmetry, chunk, seed=seed0 + i)
        for i in range(total // chunk)
    ]
    return np.concatenate(parts)
```

:::

**Compute**

```python
angles = sampled_angles('m-3m', 300000, seed0=300)
result = float((angles < 15.0).mean())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `mdf-cubic-random-low-angle-fraction` | 0.0223 | 0.0223 | &mdash; | 3.33e-05 | 2e-03 | ✅ pass |

**Why this value**: Monte-Carlo estimate over 3e5 independent Haar-random misorientations reduced by the 24 proper cubic rotations; cross-checked against an independently written quaternion implementation. Tolerance covers the binomial standard error and the seed.

**Citation**: Randle and Engler, Introduction to Texture Analysis - the conventional 15 degree low-angle boundary threshold.

**See also**: {doc}`Random disorientation baseline <../../theory/random_disorientation_baseline>`, {doc}`Orientation space and disorientation <../../theory/orientation_space_and_disorientation>`
