<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Pole-figure raster sampling

Why a measured pole figure needs solid-angle weights: the unweighted mean of cos^2 over a tilt raster is exactly 1/2 against a true spherical mean of 1/3, a 50 percent bias that survives halving the raster step, while the weighted mean converges.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## An unweighted raster mean of cos^2 is 1/2, not 1/3

Average cos^2(psi) over a 5 degree tilt/rotation raster without weights. Because the raster is uniform in polar angle rather than in solid angle, the result is the elementary integral of cos^2 over [0, pi/2] divided by pi/2, which is exactly 1/2 - while the true spherical mean is 1/3. The unweighted answer is 3/2 times the correct one, a 50 percent error, and it is a bias rather than a discretisation error: the companion example shows it is unchanged by refining the raster.

**Symbols**

- $\psi$ &mdash; Polar (tilt) angle of a specimen direction from the specimen +Z axis.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.core.sphere import raster_solid_angle_weights


def raster_polar_angles(step_deg):
    polar = np.arange(0.0, 90.0 + step_deg, step_deg)
    azimuth = np.arange(0.0, 360.0, step_deg)
    grid, _ = np.meshgrid(polar, azimuth, indexing="ij")
    return grid.ravel()
```

:::

**Compute**

```python
polar = raster_polar_angles(5.0)
field = np.cos(np.radians(polar)) ** 2
result = float(field.mean())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `pole-figure-raster-unweighted-mean-is-biased` | 0.500000 | 0.500000 | &mdash; | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: Analytic: an unweighted mean over a raster uniform in psi evaluates the integral of cos^2(psi) d(psi) over [0, pi/2] divided by pi/2, which is exactly 1/2. The correct spherical mean is 1/3.

**Citation**: Randle and Engler, Introduction to Texture Analysis - measured pole-figure rasters and their sampling geometry.

**See also**: {doc}`Pole-figure arithmetic and the m.r.d. scale <../../theory/pole_figure_arithmetic_and_mrd>`

## Halving the raster step leaves the unweighted error at exactly 50 percent

Repeat the unweighted average at a 2.5 degree step, half the previous one. The answer is again exactly 1/2. This is the signature of a biased estimator: a discretisation error would halve, while a bias from integrating against the wrong measure is independent of resolution. Buying a finer scan cannot fix it; only weighting can. The example returns the ratio of the fine-step mean to the coarse-step mean, which must be exactly 1.

**Symbols**

- $\psi$ &mdash; Polar (tilt) angle of a specimen direction from the specimen +Z axis.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.core.sphere import raster_solid_angle_weights


def raster_polar_angles(step_deg):
    polar = np.arange(0.0, 90.0 + step_deg, step_deg)
    azimuth = np.arange(0.0, 360.0, step_deg)
    grid, _ = np.meshgrid(polar, azimuth, indexing="ij")
    return grid.ravel()
```

:::

**Compute**

```python
coarse = np.cos(np.radians(raster_polar_angles(5.0))) ** 2
fine = np.cos(np.radians(raster_polar_angles(2.5))) ** 2
result = float(fine.mean() / coarse.mean())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `pole-figure-raster-bias-survives-refinement` | 1.000000000 | 1.000000000 | &mdash; | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: Analytic: both unweighted means equal 1/2 exactly, independent of the raster step, so their ratio is 1. A discretisation error would not behave this way.

**Citation**: Bunge, Texture Analysis in Materials Science (1982) - normalization of pole figures as spherical integrals.

**See also**: {doc}`Pole-figure arithmetic and the m.r.d. scale <../../theory/pole_figure_arithmetic_and_mrd>`

## Solid-angle weights recover the spherical mean, and refine correctly

Average the same field using raster_solid_angle_weights, which give each ring the solid angle of the band midway to its neighbours, cos(lower) - cos(upper), shared among its points. The weighted mean approaches the exact 1/3, and unlike the unweighted one it improves with resolution: the error falls from 4.1 percent at 5 degrees to 2.1 percent at 2.5 degrees. The example returns both weighted means so the convergence is visible.

**Symbols**

- $w_i$ &mdash; Solid-angle integration weight of sampled direction i; weights sum to 1 over the sampled region.
- $\psi$ &mdash; Polar (tilt) angle of a specimen direction from the specimen +Z axis.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.core.sphere import raster_solid_angle_weights


def raster_polar_angles(step_deg):
    polar = np.arange(0.0, 90.0 + step_deg, step_deg)
    azimuth = np.arange(0.0, 360.0, step_deg)
    grid, _ = np.meshgrid(polar, azimuth, indexing="ij")
    return grid.ravel()
```

:::

**Compute**

```python
means = []
for step in (5.0, 2.5):
    polar = raster_polar_angles(step)
    weights = raster_solid_angle_weights(polar)
    field = np.cos(np.radians(polar)) ** 2
    means.append(float((weights * field).sum()))
result = np.array(means)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `pole-figure-raster-weighted-mean-converges` | [0.31960, 0.32627] | [0.31960, 0.32627] | &mdash; | 3.46e-06 | 1e-04 | ✅ pass |

**Why this value**: Converging on the exact spherical mean 1/3 = 0.33333 from below, with the ring-band quadrature error halving as the step halves.

**Citation**: Bunge, Texture Analysis in Materials Science (1982) - spherical integration of pole densities.

**See also**: {doc}`Pole-figure arithmetic and the m.r.d. scale <../../theory/pole_figure_arithmetic_and_mrd>`

## A random standard recovers the defocusing curve and correction order exactly

A synthetic, explicitly labelled random standard has raw ring intensities 110, 90 and 60 counts at 0, 30 and 60 degrees, with a known 10-count background. Subtracting that background and normalizing to the zero-tilt ring must give defocusing factors 1, 0.8 and 0.5. A specimen point measured as 12 counts at 60 degrees with a 2-count background then corrects to (12 - 2) / 0.5 = 20. This last value distinguishes the physical order from the incorrect 12 / 0.5 - 2 = 22.

**Symbols**

- $d_i$ &mdash; Random-standard defocusing factor at specimen tilt i.
- $b$ &mdash; Declared constant background in the same intensity units.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CRYSTAL_FRAME,
    SPECIMEN_FRAME,
    Lattice,
    Phase,
    PoleFigure,
    SymmetrySpec,
    defocus_from_random_standard,
    spherical_angles_to_directions,
)
from pytex.core.lattice import CrystalPlane, MillerIndex

lattice = Lattice(3.52, 3.52, 3.52, 90.0, 90.0, 90.0, crystal_frame=CRYSTAL_FRAME)
phase = Phase(
    "synthetic nickel reference",
    lattice=lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=CRYSTAL_FRAME),
    crystal_frame=CRYSTAL_FRAME,
)
pole = CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase)
tilt = np.repeat([0.0, 30.0, 60.0], 4)
azimuth = np.tile([0.0, 90.0, 180.0, 270.0], 3)
standard = PoleFigure(
    pole=pole,
    sample_directions=spherical_angles_to_directions(tilt, azimuth),
    intensities=np.repeat([110.0, 90.0, 60.0], 4),
    specimen_frame=SPECIMEN_FRAME,
    antipodal=True,
    sampling="sampled_density",
)
target = PoleFigure(
    pole=pole,
    sample_directions=spherical_angles_to_directions([60.0], [0.0]),
    intensities=np.array([12.0]),
    specimen_frame=SPECIMEN_FRAME,
    antipodal=True,
    sampling="sampled_density",
)
```

:::

**Compute**

```python
calibration = defocus_from_random_standard(
    standard, background=10.0, synthetic=True
)
corrected = calibration.correction_spec(
    target, background=2.0, missing_intensity_policy="raise"
).apply(target)
result = np.concatenate((calibration.defocus_factors, corrected.intensities))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `pole-figure-random-standard-defocus-calibration` | [1.000000, 0.800000, 0.500000, 20.000000] | [1.000000, 0.800000, 0.500000, 20.000000] | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: Direct arithmetic independent of PyTex: ([110, 90, 60] - 10) / (110 - 10) = [1, 0.8, 0.5], and the corrected specimen value is (12 - 2) / 0.5 = 20.

**Citation**: MTEX, ODF Reconstruction from X-Ray Diffraction Data of an Al Alloy Rolled Sheet, Background and Defocusing Correction: pf = (pf - pf_bg) ./ pf_def.

**See also**: {doc}`Pole-figure correction and residual auditing <../../theory/foundation_feature_priorities>`, {doc}`Texture reconstruction correction workflow <../../workflows/foundation_feature_priorities>`, {doc}`Pole-figure arithmetic and the m.r.d. scale <../../theory/pole_figure_arithmetic_and_mrd>`
