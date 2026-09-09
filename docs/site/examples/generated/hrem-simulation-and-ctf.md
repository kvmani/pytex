<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# HRTEM simulation and contrast transfer function optics

Optics and contrast transfer of high-resolution transmission electron microscopy: relativistic electron wavelength, Scherzer defocus and resolution, analytical zero-crossing of the CTF, chromatic aberration damping reduction in double-corrected TEM, phase contrast sign inversion in Negative Cs Imaging (NCSI), and the directional splitting of point resolution that a residual two-fold astigmatism imposes on a corrected lens.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Relativistic wavelength, Scherzer defocus, and point resolution at 300 kV

In conventional high-resolution transmission electron microscopy (HRTEM), Scherzer defocus balances the phase shifts induced by the spherical aberration of the objective lens against the defocus term to maximize the bandwidth of constant phase contrast. For an accelerating voltage E_0 = 300 keV and third-order spherical aberration C_s = 1.0 mm, compute the relativistic electron wavelength lambda (in pm), the Scherzer underfocus Delta f_Sch = -1.2 sqrt(C_s lambda) (in nm), and the Scherzer point resolution d_Sch = 0.64 (C_s lambda^3)^(1/4) (in Å).

**Symbols**

- $\lambda$ &mdash; Relativistic electron radiation wavelength.
- $\Delta f$ &mdash; Defocus of the objective lens; negative values correspond to underfocus.
- $C_{s}$ &mdash; Third-order spherical aberration coefficient of the objective lens.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
    relativistic_wavelength_angstrom,
)
```

:::

**Compute**

```python
optics = MicroscopeAberrations.conventional_tem(energy_kev=300.0, cs_mm=1.0)
result = [
    round(optics.wavelength_angstrom * 1e3, 4),
    round(abs(optics.scherzer_defocus_angstrom) * 0.1, 3),
    round(optics.scherzer_resolution_angstrom, 3),
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hrem-scherzer-optics-300kv` | [19.688, 53.245, 1.892] | [19.688, 53.245, 1.892] | &mdash; | < 1e-05 | 1e-03 | ✅ pass |

**Why this value**: Analytical relativistic de Broglie wavelength for 300 keV electrons (lambda = 1.96875 pm), Scherzer defocus Delta f_Sch = -53.245 nm, and Scherzer point resolution d_Sch = 1.892 Å.

**Citation**: Scherzer, O. (1949). The theoretical resolution limit of the electron microscope. J. Appl. Phys. 20, 20-29; Williams & Carter, Transmission Electron Microscopy, 2nd ed. (Springer, 2009), Chapter 28.

**See also**: {doc}`HRTEM multislice and CTF theory <../../theory/hrem_multislice_and_ctf>`, {doc}`Diffraction API <../../api/index>`

## Analytical zero crossing of the CTF at the Scherzer passband edge

Under Scherzer underfocus Delta f_Sch = -1.2 sqrt(C_s lambda), the wave aberration phase chi(q) = pi Delta f lambda q^2 + 0.5 pi C_s lambda^3 q^4 possesses an exact non-trivial root at spatial frequency q_cross = sqrt(2.4 / sqrt(C_s lambda^3)). At this frequency, chi(q_cross) = 0 and sin(chi(q_cross)) vanishes identically to floating-point precision, marking the physical edge of the primary broad phase-contrast passband before rapid spatial frequency oscillations begin.

**Symbols**

- $\lambda$ &mdash; Relativistic electron radiation wavelength.
- $\Delta f$ &mdash; Defocus of the objective lens; negative values correspond to underfocus.
- $C_{s}$ &mdash; Third-order spherical aberration coefficient of the objective lens.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
    relativistic_wavelength_angstrom,
)
```

:::

**Compute**

```python
optics = MicroscopeAberrations.conventional_tem(energy_kev=300.0, cs_mm=1.0)
lam = optics.wavelength_angstrom
cs = optics.cs_angstrom
q_cross = np.sqrt(2.4 / np.sqrt(cs * (lam**3)))
chi = float(optics.wave_aberration(np.array([q_cross]))[0])
result = float(abs(np.sin(chi)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hrem-ctf-scherzer-zero-crossing` | < 1e-12 | 0.00e+00 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: An exact analytical identity: Delta f + 0.5 Cs lambda^2 q^2 = 0 implies chi(q_cross) = 0 and sin(chi) = 0 to machine precision.

**Citation**: Reimer, L. & Kohl, H. Transmission Electron Microscopy: Physics of Image Formation, 5th ed. (Springer, 2008), Chapter 6.

**See also**: {doc}`HRTEM multislice and CTF theory <../../theory/hrem_multislice_and_ctf>`, {doc}`Diffraction API <../../api/index>`

## Double-corrected information limit extension by chromatic aberration damping reduction

In a conventional 300 kV TEM with focal spread Delta = 35 Å (dominated by Cc and energy spread Delta E), chromatic damping envelope E_c(q) = exp(-0.5 pi^2 lambda^2 Delta^2 q^4) rapidly suppresses high spatial frequencies. In a double-corrected instrument where chromatic aberration is corrected to Delta = 5 Å, the 1/e^2 temporal information cutoff frequency q_info = sqrt(2 / (pi lambda Delta)) scales as 1 / sqrt(Delta). The ratio of double-corrected to uncorrected cutoff frequency is an exact algebraic ratio sqrt(35 / 5) = sqrt(7).

**Symbols**

- $\lambda$ &mdash; Relativistic electron radiation wavelength.
- $C_{c}$ &mdash; Chromatic aberration coefficient of the objective lens.
- $\Delta f$ &mdash; Defocus of the objective lens; negative values correspond to underfocus.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
    relativistic_wavelength_angstrom,
)
```

:::

**Compute**

```python
uncorr = MicroscopeAberrations(energy_kev=300.0, cs_mm=1.0, focal_spread_angstrom=35.0)
corr = MicroscopeAberrations(energy_kev=300.0, cs_mm=0.0, focal_spread_angstrom=5.0)
q_uncorr = np.sqrt(2.0 / (np.pi * uncorr.wavelength_angstrom * 35.0))
q_corr = np.sqrt(2.0 / (np.pi * corr.wavelength_angstrom * 5.0))
ratio = q_corr / q_uncorr
expected_ratio = np.sqrt(35.0 / 5.0)
result = float(abs(ratio - expected_ratio))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hrem-double-corrected-information-limit` | 0.00e+00 | 0.00e+00 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: An algebraic scaling identity: for temporal coherence damping, the 1/e^2 cutoff frequency scales inversely with sqrt(Delta). The ratio of cutoff frequencies between Delta = 35 Å and Delta = 5 Å is sqrt(7) exactly.

**Citation**: Rose, H. (2009). Historical aspects of aberration correction. J. Electron Microsc. 58, 77-85; Haider, M. et al. (1998). Electron microscopy image enhanced. Nature 392, 768-769.

**See also**: {doc}`HRTEM multislice and CTF theory <../../theory/hrem_multislice_and_ctf>`, {doc}`Diffraction API <../../api/index>`

## Negative Cs Imaging (NCSI) atomic column contrast inversion

Under conventional Scherzer conditions (positive Cs, underfocus Delta f < 0), sin(chi) < 0 across the passband, causing atomic columns to appear as dark minima on a bright background in bright-field HRTEM. Under Negative Spherical Aberration Imaging (NCSI, Cs < 0, overfocus Delta f > 0), sin(chi) > 0 across the passband, causing atomic columns to appear as bright maxima on a dark background. Verify this contrast inversion for an isolated gold atom: the center intensity is below the mean background (< 1.0) for Scherzer and above the mean background (> 1.0) for NCSI.

**Symbols**

- $C_{s}$ &mdash; Third-order spherical aberration coefficient of the objective lens.
- $\Delta f$ &mdash; Defocus of the objective lens; negative values correspond to underfocus.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
    relativistic_wavelength_angstrom,
)
```

:::

**Compute**

```python
snap = AtomicSnapshot(
    species=('Au',),
    positions=np.array([[2.0, 2.0, 1.0]]),
    cell=np.diag([4.0, 4.0, 4.0]),
)
scherzer_res = pure_python_phase_object_simulation(
    snap, MicroscopeAberrations.conventional_tem(energy_kev=300.0, cs_mm=1.0), sampling_angstrom=0.1
)
ncsi_res = pure_python_phase_object_simulation(
    snap, MicroscopeAberrations.ncsi(energy_kev=300.0), sampling_angstrom=0.1
)
ny, nx = scherzer_res.image.shape
cx, cy = nx // 2, ny // 2
scherzer_dark = float(scherzer_res.image[cy, cx] < 1.0)
ncsi_bright = float(ncsi_res.image[cy, cx] > 1.0)
result = [scherzer_dark, ncsi_bright]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hrem-ncsi-contrast-inversion` | [1, 1] | [1, 1] | &mdash; | exact | exact | ✅ pass |

**Why this value**: Physical principle of NCSI: phase contrast sign reversal between conventional underfocus HRTEM (dark atoms) and NCSI (bright atoms).

**Citation**: Jia, C. L., Lentzen, M. & Urban, K. (2003). Atomic-resolution imaging of oxygen in perovskite ceramics. Science 299, 870-873; Urban, K. W. (2008). Studying microstructure with aberration-corrected transmission electron microscopy. Science 321, 506-510.

**See also**: {doc}`HRTEM multislice and CTF theory <../../theory/hrem_multislice_and_ctf>`, {doc}`Diffraction API <../../api/index>`

## Two-fold astigmatism acts as a defocus offset along its own azimuth

Two-fold astigmatism enters the wave aberration as pi lambda q^2 C_12 cos(2(theta - phi_12)). Along its own azimuth, theta = phi_12, the cosine equals +1 and the term is algebraically indistinguishable from adding C_12 to the defocus; ninety degrees away the cosine equals -1 and it subtracts the same amount. For a 300 kV lens at Delta f = -50 Å with Cs = 1 um and C_12 = 20 Å at phi_12 = 30 degrees, the point resolution of the cut along phi_12 must therefore equal that of a round lens at Delta f = -30 Å, and the cut across it that of a round lens at Delta f = -70 Å. Compute the difference between each cut and its equivalent round lens.

**Symbols**

- $C_{12}$ &mdash; Two-fold astigmatism amplitude of the objective lens.
- $\varphi_{12}$ &mdash; Azimuth of the two-fold astigmatism axis in the back focal plane.
- $\theta_{q}$ &mdash; Azimuth in the back focal plane at which a transfer profile is cut.
- $\Delta f$ &mdash; Defocus of the objective lens; negative values correspond to underfocus.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
    relativistic_wavelength_angstrom,
)
```

:::

**Compute**

```python
astigmatic = MicroscopeAberrations(
    energy_kev=300.0,
    defocus_angstrom=-50.0,
    cs_mm=0.001,
    astigmatism_angstrom=20.0,
    astigmatism_angle_deg=30.0,
)
def round_lens(defocus):
    return MicroscopeAberrations(energy_kev=300.0, defocus_angstrom=defocus, cs_mm=0.001)
def d0(lens, azimuth=0.0):
    return lens.evaluate_ctf_1d(2.5, 4000, azimuth_deg=azimuth).point_resolution_angstrom
result = [
    abs(d0(astigmatic, 30.0) - d0(round_lens(-30.0))),
    abs(d0(astigmatic, 120.0) - d0(round_lens(-70.0))),
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hrem-astigmatism-equivalent-defocus` | [0.00e+00, 0.00e+00] | [0.00e+00, 0.00e+00] | Å | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: An exact algebraic identity of the wave aberration function: at theta = phi_12 the C_12 term reduces to pi lambda q^2 C_12, which is the defocus term with Delta f -> Delta f + C_12, so the two lenses have identical chi and hence identical zero crossings.

**Citation**: Krivanek, O. L., Dellby, N. & Lupini, A. R. (1999). Towards sub-A electron beams. Ultramicroscopy 78, 1-11; Kirkland, E. J. Advanced Computing in Electron Microscopy, 2nd ed. (Springer, 2010), Chapter 3.

**See also**: {doc}`HRTEM multislice and CTF theory <../../theory/hrem_multislice_and_ctf>`, {doc}`Diffraction API <../../api/index>`

## Orthogonal split of point resolution under two-fold astigmatism

Because the C_12 term varies as cos(2 theta), its period in azimuth is 180 degrees and its two extremes lie 90 degrees apart. A lens is therefore resolved best along one direction and worst along the perpendicular one, and the separation between those two directions is fixed by the multiplicity of the aberration rather than by its size. Sample the transfer function of the same 300 kV lens over 180 azimuths and report the angular separation, modulo 180 degrees, between the finest and the coarsest point resolution.

**Symbols**

- $C_{12}$ &mdash; Two-fold astigmatism amplitude of the objective lens.
- $\varphi_{12}$ &mdash; Azimuth of the two-fold astigmatism axis in the back focal plane.
- $\theta_{q}$ &mdash; Azimuth in the back focal plane at which a transfer profile is cut.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.diffraction.hrem import (
    AtomicSnapshot,
    DoubleCorrectionMode,
    MicroscopeAberrations,
    pure_python_phase_object_simulation,
    relativistic_wavelength_angstrom,
)
```

:::

**Compute**

```python
lens = MicroscopeAberrations(
    energy_kev=300.0,
    defocus_angstrom=-50.0,
    cs_mm=0.001,
    astigmatism_angstrom=20.0,
    astigmatism_angle_deg=30.0,
)
band = lens.evaluate_ctf_azimuthal(2.5, 4000, 180)
resolutions = band.point_resolution_angstrom
best = float(band.azimuths_deg[int(np.nanargmin(resolutions))])
result = float(abs(band.worst_azimuth_deg - best) % 180.0)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hrem-astigmatic-resolution-split` | 90.0 | 90.0 | deg | < 1e-02 | 1e+00 | ✅ pass |

**Why this value**: The azimuthal multiplicity of two-fold astigmatism: cos(2 theta) has period 180 degrees, so its maximum and minimum are separated by exactly 90 degrees. The tolerance is one degree, the spacing of the 180-point azimuthal grid.

**Citation**: Krivanek, O. L., Dellby, N. & Lupini, A. R. (1999). Towards sub-A electron beams. Ultramicroscopy 78, 1-11.

**See also**: {doc}`HRTEM multislice and CTF theory <../../theory/hrem_multislice_and_ctf>`, {doc}`Diffraction API <../../api/index>`
