<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Texture kernels

Analytic identities of the SO(3) kernel surface - normalization (A_0 = 1) and the halfwidth definition - together with the m.r.d. scale on which pole densities are reported, plus an exactly identifiable named-component mixture, all computed live.

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
| `texture-gaussian-kernel-normalization-and-halfwidth` | [1.000000, 0.500000] | [1.000000, 0.500000] | &mdash; | 7.47e-14 | 1e-06 | ✅ pass |

**Why this value**: A_0 = 1 is the SO(3) normalization identity for character expansions, and psi(halfwidth) = psi(0)/2 is the definition of the kernel halfwidth (both analytic identities).

**Citation**: Bunge, Texture Analysis in Materials Science (1982), harmonic expansion of ODF kernels.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Texture API <../../api/index>`

## A uniform ODF gives a pole figure flat at one m.r.d.

Build an ODF from an equispaced grid over the cubic fundamental region with equal weights - a texture-free aggregate - and evaluate its {111} pole density along three unrelated specimen directions. A discrete ODF's evaluate_pole_density returns a kernel-weighted response, whose peak is one rather than whose integral is one, so converting it to multiples of a random distribution means dividing by random_pole_density: the response a random texture produces. Skipping that division is a scale error of about two orders of magnitude, not a small one.

**Symbols**

- $P_{\mathbf{h}}(\mathbf{y})$ &mdash; Pole density of plane normal h along specimen direction y, in multiples of a random distribution.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    ODF,
    CrystalPlane,
    FrameDomain,
    KernelSpec,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    TextureComponent,
    fit_odf_components,
    random_pole_density,
)

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
dictionary = OrientationSet.from_equispaced_so3_grid(
    10.0,
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
    phase=phase,
)
kernel = KernelSpec(name='de_la_vallee_poussin', halfwidth_deg=15.0)
odf = ODF.from_orientations(dictionary, kernel=kernel)
directions = np.array(
    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0] / np.sqrt(3.0)]
)
response = np.asarray(odf.evaluate_pole_density(pole, directions))
result = response / random_pole_density(kernel)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `texture-uniform-odf-pole-density-is-one-mrd` | [1.000013, 1.000248, 0.999713] | [1.000000, 1.000000, 1.000000] | m.r.d. | 2.87e-04 | 1e-03 | ✅ pass |

**Why this value**: Analytic identity: a uniform orientation distribution maps poles uniformly onto the sphere, so every pole density equals one multiple of a random distribution by the definition of the m.r.d. scale. The residual deviation is the finite orientation grid, not the scale.

**Citation**: Bunge, Texture Analysis in Materials Science (1982), Sec. 4 - normalization of the ODF and of pole figures to multiples of a random distribution.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Texture API <../../api/index>`

## A named ODF fit recovers a 70/30 cube-Goss mixture

A reconstructed sheet-texture ODF appears dominated by cube and Goss peaks. Fit those named components plus a random term and retain the density residual, rather than integrating independent angular windows that may overlap. Here the ODF is constructed algebraically from 70% cube and 30% Goss kernels; zero-weight dictionary points make the inverse basis identifiable without changing the known mixture.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    ODF,
    CrystalPlane,
    FrameDomain,
    KernelSpec,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    TextureComponent,
    fit_odf_components,
    random_pole_density,
)

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
components = (
    TextureComponent('cube', (0.0, 0.0, 0.0), '{001}<100>'),
    TextureComponent('goss', (0.0, 45.0, 0.0), '{011}<100>'),
)
euler = np.array([
    [0.0, 0.0, 0.0], [0.0, 45.0, 0.0], [15.0, 20.0, 10.0],
    [30.0, 55.0, 15.0], [70.0, 35.0, 40.0], [10.0, 75.0, 80.0],
])
support = OrientationSet.from_euler_angles(
    euler, crystal_frame=crystal, specimen_frame=specimen,
    symmetry=symmetry, phase=phase, degrees=True,
)
odf = ODF.from_orientations(
    support, weights=[0.7, 0.3, 0.0, 0.0, 0.0, 0.0],
    kernel=KernelSpec(halfwidth_deg=12.0),
)
fit = fit_odf_components(odf, components)
result = np.array([fit.fraction_for('cube'), fit.fraction_for('goss'), fit.random_fraction])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `texture-odf-named-component-fit-recovers-exact-mixture` | [0.700000000, 0.300000000, 0.000000000] | [0.700000000, 0.300000000, 0.000000000] | volume fraction | 1.64e-16 | 1e-09 | ✅ pass |

**Why this value**: The target density is defined as the exact linear combination 0.7*K(g,cube) + 0.3*K(g,Goss). The component design uses those same two normalized kernel columns plus a constant random-density column, so the unique full-rank constrained solution is [0.7, 0.3, 0.0].

**Citation**: Bunge, Texture Analysis in Materials Science (1982), component method and normalized ODFs; SciPy SLSQP documentation for constrained least squares.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Texture API <../../api/index>`
