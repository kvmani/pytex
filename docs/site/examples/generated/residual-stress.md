<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Residual stress by the sin^2(psi) method

The diffraction elastic constants checked against the isotropic formula, the cubic Reuss formula worked by hand and the root of Kroener's cubic; the sin^2(psi) law's slope, strain-free tilt and principal stresses checked against closed-form algebra on exact data; and an end-to-end evaluation of noisy simulated scans checked against the stress they were generated with.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## The isotropic diffraction elastic constant of a steel

For an elastically isotropic solid Hooke's law projected on the scattering vector gives 1/2 S2 = (1 + nu)/E. For a steel with E = 210 GPa and nu = 0.28 that is 1.28 / 210 GPa = 6.0952 x 10^-3 / GPa = 6.0952 / TPa: a 100 MPa stress along the scattering vector strains the planes by 6.1 x 10^-4 through this term.

**Symbols**

- $\tfrac{1}{2}S_{2}$ &mdash; Second diffraction elastic constant of the reflection.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.diffraction.xrd_residual_stress import DiffractionElasticConstants
```

:::

**Compute**

```python
result = DiffractionElasticConstants.isotropic(210.0, 0.28).half_s2_per_tpa
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `stress-isotropic-half-s2` | 6.0952 | 6.0952 | TPa^-1 | 3.81e-05 | 1e-04 | ✅ pass |

**Why this value**: 1.28 / 210 = 0.0060952 per GPa = 6.0952 per TPa, by hand.

**Citation**: Noyan & Cohen, Residual Stress, Springer (1987), doi:10.1007/978-1-4613-9570-6.

**See also**: {doc}`Residual stress by the sin^2(psi) method <../../theory/residual_stress_sin2psi>`

## The Reuss constant 1/2 S2 of ferrite (211), by the cubic formula

Under the Reuss model a cubic reflection has 1/2 S2 = S11 - S12 - 3 S0 Gamma with S0 = S11 - S12 - S44/2 and Gamma = (h^2 k^2 + k^2 l^2 + l^2 h^2)/(h^2 + k^2 + l^2)^2. For ferrite, C11 = 231.4, C12 = 134.7 and C44 = 116.4 GPa (Simmons & Wang). With (C11 - C12)(C11 + 2 C12) = 96.7 x 500.8 = 48427.4 GPa^2: S11 = 366.1 / 48427.4 = 7.5598 /TPa, S12 = -134.7 / 48427.4 = -2.7815 /TPa and S44 = 1/116.4 = 8.5911 /TPa, so S0 = 6.0458 /TPa. For (211), Gamma = (4 + 1 + 4)/36 = 1/4, and 1/2 S2 = 10.3413 - 3 x 6.0458 / 4 = 5.8069 /TPa. The library averages the single-crystal compliance over rotations about the plane normal in Mandel form; it must land on the same number.

**Symbols**

- $\tfrac{1}{2}S_{2}$ &mdash; Second diffraction elastic constant of the reflection.
- $\Gamma$ &mdash; Cubic orientation parameter of a reflection.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.diffraction.xrd_residual_stress import (
    DiffractionElasticConstants, single_crystal_stiffness,
)
```

:::

**Compute**

```python
result = DiffractionElasticConstants.from_single_crystal(
    single_crystal_stiffness('fe_bcc'), [2, 1, 1], model='reuss'
).half_s2_per_tpa
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `stress-reuss-ferrite-211` | 5.8070 | 5.8069 | TPa^-1 | 6.49e-05 | 2e-04 | ✅ pass |

**Why this value**: By hand from the cubic Reuss formula: S11 - S12 = 10.3413, 3 S0 Gamma = 4.5344, 1/2 S2 = 5.8069 /TPa.

**Citation**: Welzel et al., J. Appl. Cryst. 38 (2005) 1, doi:10.1107/S0021889804029516; Simmons & Wang, Single Crystal Elastic Constants, MIT Press (1971).

**See also**: {doc}`Residual stress by the sin^2(psi) method <../../theory/residual_stress_sin2psi>`

## Kroener's self-consistent shear modulus of a ferrite aggregate

For cubic crystals Kroener's self-consistent shear modulus is the positive root of G^3 + (5 C11 + 4 C12)/8 G^2 - C44 (7 C11 - 4 C12)/8 G - C44 (C11 - C12)(C11 + 2 C12)/8 = 0. For ferrite the coefficients are 211.975, -15728.55 and -704618.09 (GPa units), and the positive root is G = 82.448 GPa, between the Reuss (74.47) and Voigt (89.18) shear moduli. The library does not solve this cubic: it iterates the Eshelby-sphere concentration tensor in the general (any crystal system) scheme, which must converge to the root.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.diffraction.xrd_residual_stress import (
    _kroener_grain_compliance, _to_mandel, single_crystal_stiffness,
)
```

:::

**Compute**

```python
stiffness = np.asarray(single_crystal_stiffness('fe_bcc').tensor)
_, bulk, shear = _kroener_grain_compliance(_to_mandel(stiffness))
result = shear
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `stress-kroener-shear-modulus-ferrite` | 82.448 | 82.448 | GPa | 3.63e-04 | 1e-03 | ✅ pass |

**Why this value**: Substituting G = 82.448 into the cubic: 560 454.5 + 1 440 936.7 - 1 296 787.5 - 704 618.1 = -14.4, against a derivative 3G^2 + 2 alpha G + beta = 39 618 per GPa, so the root is 82.448 + 0.0004 GPa.

**Citation**: Kroener, Z. Physik 151 (1958) 504, doi:10.1007/BF01337948.

**See also**: {doc}`Residual stress by the sin^2(psi) method <../../theory/residual_stress_sin2psi>`

## The slope of d against sin^2(psi) is d0 1/2 S2 sigma_phi

Under plane stress d = d0 [1 + 1/2 S2 sigma_phi sin^2(psi) + S1 (sigma_11 + sigma_22)], which is exactly linear in sin^2(psi). For a uniaxial sigma_11 = -300 MPa measured at phi = 0 with 1/2 S2 = 5.8 /TPa and d0 = 1.1702 angstrom, the slope is 1.1702 x 5.8 x 10^-6 x (-300) = -2.03615 x 10^-3 angstrom. Exact peak positions are generated from the fundamental equation, and the per-azimuth regression must return that slope.

**Symbols**

- $d_{0}$ &mdash; Stress-free interplanar spacing.
- $\tfrac{1}{2}S_{2}$ &mdash; Second diffraction elastic constant of the reflection.
- $\sigma_{\varphi}$ &mdash; Normal stress along azimuth phi in the surface.
- $\psi$ &mdash; Tilt of the scattering vector from the surface normal.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.diffraction.xrd_residual_stress import (
    DiffractionElasticConstants, StressPeak, determine_residual_stress,
    measurement_direction,
)
WAVELENGTH = 2.2897
def exact_peaks(stress, dec, d0, phis, psis):
    peaks = []
    for phi in phis:
        for psi in psis:
            m = measurement_direction(phi, psi)
            strain = 1e-6 * (dec.half_s2_per_tpa * float(m @ stress @ m)
                             + dec.s1_per_tpa * float(np.trace(stress)))
            two_theta = math.degrees(2 * math.asin(WAVELENGTH / (2 * d0 * (1 + strain))))
            peaks.append(StressPeak(phi_deg=phi, psi_deg=psi, two_theta_deg=two_theta,
                                    two_theta_uncertainty_deg=1e-3, method='given'))
    return peaks
```

:::

**Compute**

```python
dec = DiffractionElasticConstants(s1_per_tpa=-1.25, half_s2_per_tpa=5.8)
stress = np.diag([-300.0, 0.0, 0.0])
peaks = exact_peaks(stress, dec, 1.1702, (0.0,), (0.0, 15.0, 25.0, 35.0, 45.0))
fit = determine_residual_stress(peaks, wavelength_angstrom=WAVELENGTH,
                                d0_angstrom=1.1702, dec=dec)
result = 1000.0 * fit.regressions[0].slope_angstrom
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `stress-sin2psi-slope-identity` | -2.03615 | -2.03615 | mÅ | 2.00e-06 | 2e-05 | ✅ pass |

**Why this value**: 1.1702 x 5.8e-6 x (-300) = -2.03615e-3 angstrom, by hand.

**Citation**: Macherauch & Mueller, Z. angew. Phys. 13 (1961) 305.

**See also**: {doc}`Residual stress by the sin^2(psi) method <../../theory/residual_stress_sin2psi>`, {doc}`Residual stress: the algorithm and the report <../../algorithms/residual_stress_sin2psi>`

## The strain-free tilt of an equibiaxial stress

Under an equibiaxial stress the sin^2(psi) line crosses d = d0 at sin^2(psi*) = -2 S1 / (1/2 S2), independent of the stress. For isotropic constants that is 2 nu / (1 + nu); with nu = 0.28, 0.56 / 1.28 = 0.4375 (psi* = 41.4 degrees). The tilt is read off the fitted line as (d0 - intercept) / slope.

**Symbols**

- $S_{1}$ &mdash; First diffraction elastic constant of the reflection.
- $\tfrac{1}{2}S_{2}$ &mdash; Second diffraction elastic constant of the reflection.
- $\psi$ &mdash; Tilt of the scattering vector from the surface normal.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.diffraction.xrd_residual_stress import (
    DiffractionElasticConstants, StressPeak, determine_residual_stress,
    measurement_direction,
)
WAVELENGTH = 2.2897
def exact_peaks(stress, dec, d0, phis, psis):
    peaks = []
    for phi in phis:
        for psi in psis:
            m = measurement_direction(phi, psi)
            strain = 1e-6 * (dec.half_s2_per_tpa * float(m @ stress @ m)
                             + dec.s1_per_tpa * float(np.trace(stress)))
            two_theta = math.degrees(2 * math.asin(WAVELENGTH / (2 * d0 * (1 + strain))))
            peaks.append(StressPeak(phi_deg=phi, psi_deg=psi, two_theta_deg=two_theta,
                                    two_theta_uncertainty_deg=1e-3, method='given'))
    return peaks
```

:::

**Compute**

```python
dec = DiffractionElasticConstants.isotropic(210.0, 0.28)
stress = np.diag([-400.0, -400.0, 0.0])
peaks = exact_peaks(stress, dec, 1.1702, (0.0, 45.0, 90.0),
                    (0.0, 15.0, 25.0, 35.0, 45.0))
fit = determine_residual_stress(peaks, wavelength_angstrom=WAVELENGTH,
                                d0_angstrom=1.1702, dec=dec)
line = fit.regressions[0]
result = (1.1702 - line.intercept_angstrom) / line.slope_angstrom
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `stress-strain-free-direction` | 0.4375 | 0.4375 | &mdash; | < 1e-06 | 1e-04 | ✅ pass |

**Why this value**: 2 x 0.28 / 1.28 = 0.4375, by hand.

**Citation**: Noyan & Cohen, Residual Stress, Springer (1987), doi:10.1007/978-1-4613-9570-6.

**See also**: {doc}`Residual stress by the sin^2(psi) method <../../theory/residual_stress_sin2psi>`

## The larger in-plane principal stress from three azimuths

Three azimuths, 0, 45 and 90 degrees, fix sigma_11, sigma_22 and sigma_12, and the in-plane principal stresses are (sigma_11 + sigma_22)/2 +/- sqrt(((sigma_11 - sigma_22)/2)^2 + sigma_12^2). For sigma_11 = -350, sigma_22 = -150 and sigma_12 = 60 MPa the larger is -250 + sqrt(100^2 + 60^2) = -250 + 116.619 = -133.381 MPa.

**Symbols**

- $\varphi$ &mdash; Azimuth of the tilt plane, from S1 towards S2.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.diffraction.xrd_residual_stress import (
    DiffractionElasticConstants, StressPeak, determine_residual_stress,
    measurement_direction,
)
WAVELENGTH = 2.2897
def exact_peaks(stress, dec, d0, phis, psis):
    peaks = []
    for phi in phis:
        for psi in psis:
            m = measurement_direction(phi, psi)
            strain = 1e-6 * (dec.half_s2_per_tpa * float(m @ stress @ m)
                             + dec.s1_per_tpa * float(np.trace(stress)))
            two_theta = math.degrees(2 * math.asin(WAVELENGTH / (2 * d0 * (1 + strain))))
            peaks.append(StressPeak(phi_deg=phi, psi_deg=psi, two_theta_deg=two_theta,
                                    two_theta_uncertainty_deg=1e-3, method='given'))
    return peaks
```

:::

**Compute**

```python
dec = DiffractionElasticConstants(s1_per_tpa=-1.23, half_s2_per_tpa=5.69)
stress = np.array([[-350.0, 60.0, 0.0], [60.0, -150.0, 0.0], [0.0, 0.0, 0.0]])
peaks = exact_peaks(stress, dec, 1.1702, (0.0, 45.0, 90.0),
                    (-45.0, -30.0, -15.0, 0.0, 15.0, 30.0, 45.0))
fit = determine_residual_stress(peaks, wavelength_angstrom=WAVELENGTH,
                                d0_angstrom=1.1702, dec=dec)
result = fit.tensor.in_plane_principal()['sigma_I_mpa']
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `stress-in-plane-principal` | -133.381 | -133.381 | MPa | 3.79e-05 | 1e-03 | ✅ pass |

**Why this value**: -250 + sqrt(13600) = -250 + 116.619 = -133.381 MPa, by hand.

**Citation**: Noyan & Cohen, Residual Stress, Springer (1987), doi:10.1007/978-1-4613-9570-6.

**See also**: {doc}`Residual stress by the sin^2(psi) method <../../theory/residual_stress_sin2psi>`, {doc}`Residual stress: the algorithm and the report <../../algorithms/residual_stress_sin2psi>`

## sigma_11 of a shot-peened ferrite, from simulated Cr K-alpha scans

Twenty-one scans of ferrite (211) with Cr K-alpha - three azimuths, seven tilts of both signs - generated from sigma_11 = -350, sigma_22 = -150, sigma_12 = 60 MPa through the Kroener constants, as K-alpha1/K-alpha2 pseudo-Voigt doublets broadened as 1/cos(psi), shaped by the LPA factor and given Poisson noise. Each peak is located by the doublet profile fit after the LPA correction, and the tensor fitted to all 21 positions. The answer must be the stress the data were generated with, to within the statistical uncertainty of a few MPa.

**Symbols**

- $\psi$ &mdash; Tilt of the scattering vector from the surface normal.
- $\varphi$ &mdash; Azimuth of the tilt plane, from S1 towards S2.


:::{dropdown} Setup (imports and object construction)

```python
import math
from pytex.diffraction.xrd_residual_stress import (
    DiffractionElasticConstants, residual_stress_pipeline,
    simulate_sin2psi_measurement, single_crystal_stiffness,
)
```

:::

**Compute**

```python
d0 = 2.8665 / math.sqrt(6.0)
dec = DiffractionElasticConstants.from_single_crystal(
    single_crystal_stiffness('fe_bcc'), [2, 1, 1], model='kroener')
scans = simulate_sin2psi_measurement(
    d0_angstrom=d0, dec=dec, seed=1,
    stress_mpa={'sigma_11': -350.0, 'sigma_22': -150.0, 'sigma_12': 60.0})
fit = residual_stress_pipeline(scans, d0_angstrom=d0, dec=dec, window_deg=8.0)
result = fit.tensor.component('sigma_11')[0]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `stress-end-to-end-ferrite` | -350.1 | -350.0 | MPa | 1.25e-01 | 8e+00 | ✅ pass |

**Why this value**: The generating stress. The tolerance is about four statistical standard uncertainties of the profile-fit route on this data set.

**Citation**: Macherauch & Mueller, Z. angew. Phys. 13 (1961) 305.

**See also**: {doc}`Residual stress by the sin^2(psi) method <../../theory/residual_stress_sin2psi>`, {doc}`Residual stress: the algorithm and the report <../../algorithms/residual_stress_sin2psi>`
