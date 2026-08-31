<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Directional statistics and mean axes

Averaging axes rather than vectors: the orientation tensor has unit trace, its eigenvalues take exact values at the girdle and cluster limits, and it recovers a fibre axis from randomly signed data where the vector resultant fails outright.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## A girdle gives eigenvalues (0, 1/2, 1/2) and a cluster (0, 0, 1)

Build two limiting direction sets - one spread uniformly around a great circle in the xy-plane, one with every axis parallel to z - and read the eigenvalues of their orientation tensors. The closed forms are exact: a girdle has <cos^2> = <sin^2> = 1/2 around the circle so its eigenvalues are (0, 1/2, 1/2), and a perfect cluster gives (0, 0, 1). Together with the uniform case (1/3, 1/3, 1/3) these are the three corners of the eigenvalue triangle, and they let three numbers classify a fabric without contouring anything. The example returns both eigenvalue triples.

**Symbols**

- $\boldsymbol{\Theta}$ &mdash; Orientation tensor (1/n) sum v v^T; the second moment of a direction set.
- $\lambda_{1} \le \lambda_{2} \le \lambda_{3}$ &mdash; Eigenvalues of the orientation tensor; non-negative and summing to one.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.core import crystal_frame
from pytex.core.sphere import SphericalVectorSet

frame = crystal_frame()


def axes(values, antipodal=True):
    return SphericalVectorSet.from_vectors(
        values, reference_frame=frame, antipodal=antipodal
    )
```

:::

**Compute**

```python
rng = np.random.default_rng(4)
angle = rng.random(100000) * 2.0 * np.pi
girdle = np.stack(
    [np.cos(angle), np.sin(angle), np.zeros_like(angle)], axis=-1
)
cluster = np.tile([0.0, 0.0, 1.0], (1000, 1))
result = np.array(
    [
        np.linalg.eigvalsh(np.asarray(axes(girdle).orientation_tensor())),
        np.linalg.eigvalsh(np.asarray(axes(cluster).orientation_tensor())),
    ]
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `directional-orientation-tensor-limiting-eigenvalues` | [0.0000, 0.4996, 0.5004, 0.0000, 0.0000, 1.0000] | [0.0000, 0.5000, 0.5000, 0.0000, 0.0000, 1.0000] | &mdash; | 3.60e-04 | 5e-03 | ✅ pass |

**Why this value**: Analytic: a uniform girdle in a plane has second moments <cos^2> = <sin^2> = 1/2 and zero out of plane, giving (0, 1/2, 1/2); identical parallel axes give (0, 0, 1).

**Citation**: Woodcock, Specification of fabric shapes using an eigenvalue method, GSA Bulletin 88 (1977) 1231-1236.

**See also**: {doc}`Directional statistics and mean axes <../../theory/directional_statistics_and_mean_axes>`

## The orientation tensor of unit vectors has trace exactly one

Compute the trace of the orientation tensor for a set of random unit directions. Because each vector is normalized, the diagonal sums to the mean of |v|^2, which is exactly one - so the eigenvalues always sum to one and live on a triangle. It is a free check on any implementation of the tensor, and it is what makes the three eigenvalues directly comparable between datasets of different size.

**Symbols**

- $\boldsymbol{\Theta}$ &mdash; Orientation tensor (1/n) sum v v^T; the second moment of a direction set.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.core import crystal_frame
from pytex.core.sphere import SphericalVectorSet

frame = crystal_frame()


def axes(values, antipodal=True):
    return SphericalVectorSet.from_vectors(
        values, reference_frame=frame, antipodal=antipodal
    )
```

:::

**Compute**

```python
rng = np.random.default_rng(4)
values = rng.normal(size=(50000, 3))
values = values / np.linalg.norm(values, axis=1, keepdims=True)
tensor = np.asarray(axes(values).orientation_tensor())
result = float(np.trace(tensor))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `directional-orientation-tensor-unit-trace` | 1.000000000000 | 1.000000000000 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: Analytic identity: tr((1/n) sum v v^T) = (1/n) sum |v|^2 = 1 for unit vectors, independent of the distribution.

**Citation**: Fisher, Lewis and Embleton, Statistical Analysis of Spherical Data (CUP 1987).

**See also**: {doc}`Directional statistics and mean axes <../../theory/directional_statistics_and_mean_axes>`

## The tensor recovers z from randomly signed axes; the resultant does not

Take 3000 axes tightly clustered about z - an unambiguous fibre - and give each an independent random sign, which is what a real axial measurement delivers. The normalized resultant comes out wrong in sign and several degrees off axis, and would change if the signs were redrawn. The orientation tensor is blind to sign, because (-v)(-v)^T = v v^T, so its principal eigenvector recovers z. The example returns the absolute z-component of the tensor mean, which must be 1.

**Symbols**

- $\boldsymbol{\Theta}$ &mdash; Orientation tensor (1/n) sum v v^T; the second moment of a direction set.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.core import crystal_frame
from pytex.core.sphere import SphericalVectorSet

frame = crystal_frame()


def axes(values, antipodal=True):
    return SphericalVectorSet.from_vectors(
        values, reference_frame=frame, antipodal=antipodal
    )
```

:::

**Compute**

```python
rng = np.random.default_rng(4)
tight = np.stack(
    [
        0.2 * rng.normal(size=3000),
        0.2 * rng.normal(size=3000),
        np.ones(3000),
    ],
    axis=-1,
)
tight = tight / np.linalg.norm(tight, axis=1, keepdims=True)
signed = tight * rng.choice([-1.0, 1.0], size=(3000, 1))
result = float(abs(axes(signed).mean_direction()[2]))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `directional-mean-axis-of-randomly-signed-axes` | 0.999982 | 1.000000 | &mdash; | < 5e-05 | 5e-03 | ✅ pass |

**Why this value**: Analytic: the axes are drawn about z, so the principal eigenvector of the orientation tensor is z and its z-component is 1. Random signs cannot affect it, since the tensor is invariant under v -> -v.

**Citation**: Mardia and Jupp, Directional Statistics (Wiley 2000) - axial data and the failure of the resultant.

**See also**: {doc}`Directional statistics and mean axes <../../theory/directional_statistics_and_mean_axes>`
