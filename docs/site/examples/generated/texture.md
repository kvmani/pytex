<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Texture kernels

Analytic identities of the SO(3) kernel surface: normalization (A_0 = 1) and the halfwidth definition, computed live.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## The Gaussian SO(3) kernel is normalized and honors its halfwidth

Construct a Gaussian (Gauss-Weierstrass) kernel with a 10 degree halfwidth and verify the two defining identities: the zeroth Chebyshev coefficient equals one (the kernel integrates to one over SO(3) with the normalized Haar measure), and the density at the halfwidth equals half the peak density.

**Symbols**

- $\psi(\omega)$ &mdash; SO(3) radial kernel density as a function of the rotation angle.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import GaussianSO3Kernel
```

:::

**Compute**

```python
kernel = GaussianSO3Kernel(10.0)
a0 = float(kernel.chebyshev_coefficients(0)[0])
ratio = float(
    kernel.evaluate(np.array([np.deg2rad(10.0)]))[0]
    / kernel.evaluate(np.array([0.0]))[0]
)
result = np.array([a0, ratio])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `texture-gaussian-kernel-normalization-and-halfwidth` | [1.000000, 0.500000] | [1.000000, 0.500000] | &mdash; | 7.48e-14 | 1e-06 | ✅ pass |

**Why this value**: A_0 = 1 is the SO(3) normalization identity for character expansions, and psi(halfwidth) = psi(0)/2 is the definition of the kernel halfwidth (both analytic identities).

**Citation**: Bunge, Texture Analysis in Materials Science (1982), harmonic expansion of ODF kernels.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Texture API <../../api/index>`
