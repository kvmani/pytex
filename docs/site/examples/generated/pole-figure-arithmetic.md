<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Pole-figure arithmetic

Exact identities behind comparing two pole figures: the multiples-of-random scale, resampling onto a shared support, and the additivity of pole densities.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## A pole figure in m.r.d. has unit mean density, and its deviation from random sums to zero

Resample a two-component texture onto an equal-area grid with the default m.r.d. normalization, then verify the two identities that define the scale: the solid-angle-weighted mean density is exactly one, and the weighted integral of the deviation from random is exactly zero. Neither holds for a figure normalized by its maximum or by its sum, which is why those scales cannot be compared between measurements.

**Symbols**

- $P_{hkl}(\mathbf{y})$ &mdash; Pole density of the family in multiples of a random distribution.
- $w_i$ &mdash; Solid-angle integration weight of sampled direction i; the weights sum to one.
- $\Delta P$ &mdash; Signed difference of two pole densities on a shared support.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    S2Grid,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
specimen = ReferenceFrame(
    name="specimen", domain=FrameDomain.SPECIMEN, axes=("x", "y", "z"), handedness=Handedness.RIGHT
)
phase = Phase(
    name="fcc",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
pole = CrystalPlane(MillerIndex((1, 1, 1), phase=phase), phase=phase)
grid = S2Grid.equispaced(10.0, reference_frame=specimen, hemisphere="upper", antipodal=True)
```

:::

**Compute**

```python
orientations = OrientationSet.from_euler_angles(
    np.array([[0.0, 0.0, 0.0], [35.0, 20.0, 10.0]]),
    specimen_frame=specimen,
    phase=phase,
)
figure = PoleFigure.from_orientations(orientations, pole).on_grid(
    grid, halfwidth_deg=15.0
)
mean_density = float(np.sum(grid.weights * figure.intensities))
deviation = figure - 1.0
deviation_integral = float(np.sum(grid.weights * deviation.values))
result = np.array([mean_density, deviation_integral])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `texture-pole-figure-mrd-unit-mean-density` | [1.000000000000, -0.000000000000] | [1.000000000000, 0.000000000000] | m.r.d. | 1.11e-16 | 1e-12 | ✅ pass |

**Why this value**: Definition of the multiples-of-random scale: sum_i w_i P_i = 1 with solid-angle weights summing to one. The deviation identity follows immediately, since sum_i w_i (P_i - 1) = 1 - 1 = 0.

**Citation**: Bunge, Texture Analysis in Materials Science (1982), normalization of pole figures to multiples of a random distribution.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Texture API <../../api/index>`

## Resampling preserves a constant field, and pole densities add

Resample a pole figure that is 2.5 m.r.d. everywhere onto a coarser grid, and add two normalized figures. The interpolating estimator is a weighted mean, so a constant field passes through it unchanged for any kernel halfwidth — the partition-of-unity property that distinguishes it from the summing estimator used for pole clouds. Densities then add pointwise, so the mean of a sum of two m.r.d. figures is exactly two.

**Symbols**

- $P_{hkl}(\mathbf{y})$ &mdash; Pole density of the family in multiples of a random distribution.
- $w_i$ &mdash; Solid-angle integration weight of sampled direction i; the weights sum to one.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    S2Grid,
    SymmetrySpec,
)

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
specimen = ReferenceFrame(
    name="specimen", domain=FrameDomain.SPECIMEN, axes=("x", "y", "z"), handedness=Handedness.RIGHT
)
phase = Phase(
    name="fcc",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
pole = CrystalPlane(MillerIndex((1, 1, 1), phase=phase), phase=phase)
grid = S2Grid.equispaced(10.0, reference_frame=specimen, hemisphere="upper", antipodal=True)
```

:::

**Compute**

```python
flat = PoleFigure(
    pole=pole,
    sample_directions=grid.vectors.values,
    intensities=np.full(len(grid), 2.5),
    specimen_frame=specimen,
    antipodal=True,
    sampling='sampled_density',
)
coarse = S2Grid.equispaced(
    18.0, reference_frame=specimen, hemisphere='upper', antipodal=True
)
resampled = flat.on_grid(coarse, halfwidth_deg=7.0, normalize=False)
constant_field = float(np.max(np.abs(resampled.intensities - 2.5)))

single = PoleFigure.from_orientations(
    OrientationSet.from_euler_angles(
        np.zeros((1, 3)), specimen_frame=specimen, phase=phase
    ),
    pole,
).on_grid(grid, halfwidth_deg=15.0)
total = single + single
sum_mean = float(np.sum(grid.weights * total.intensities))
result = np.array([constant_field, sum_mean])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `texture-pole-figure-resampling-and-addition-identities` | [0.000000000000, 2.000000000000] | [0.000000000000, 2.000000000000] | m.r.d. | 1.33e-15 | 1e-12 | ✅ pass |

**Why this value**: Partition of unity of the Nadaraya-Watson estimator: sum_i K_i f_i / sum_i K_i = c whenever every f_i = c, for any kernel and any query direction. Linearity of the solid-angle mean then gives mean(P + P) = 2 for a figure of unit mean.

**Citation**: Nadaraya (1964), Theory of Probability and its Applications 9:141, and Watson (1964), Sankhya A 26:359, for the weighted-mean estimator.

**See also**: {doc}`Orientations and texture <../../concepts/orientation_texture>`, {doc}`Texture API <../../api/index>`
