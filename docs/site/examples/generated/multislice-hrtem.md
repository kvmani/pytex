<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Multislice HRTEM simulation

The multislice engine against closed forms and a second method: the absolute scale of the projected potential and the mean inner potential, the band limit, defocus as Fresnel propagation, agreement with Bloch waves on the same potential, and focal integration against Frank's envelope.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## The projected potential of one atom integrates to (h²/2πm₀e) f_e(0)

The independent-atom potential is the Fourier transform of the electron scattering factor scaled by h^2/(2 pi m0 e) = 2 pi a0 e = 47.878 V Å^2, so the zero-frequency coefficient of a projected potential - its integral over the plane - is that constant times f_e(0). Place one silicon atom anywhere in a periodic 8 x 7 Å cell, build its projected potential from the Lobato-Van Dyck parametrization, and integrate it over the cell in V Å^3.

**Symbols**

- $v_{n}$ &mdash; Projected potential of slice n, in V Å.
- $f_{e}(s)$ &mdash; Electron atomic scattering factor in ångström.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.app.phases import phase_from_request
from pytex.diffraction.hrem import AtomicSnapshot, MicroscopeAberrations
from pytex.diffraction.multislice import (
    MultisliceGrid,
    SlicedPotential,
    TemporalCoherence,
    hrtem_image,
    multislice,
    parametrized_electron_scattering_factor,
    periodic_slab,
    slice_potential,
)
```

:::

**Compute**

```python
grid = MultisliceGrid((140, 160), (8.0, 7.0))
v = slice_potential(['Si'], np.array([[2.3, 4.1]]), grid, 'lobato')
dx, dy = grid.sampling_angstrom
result = float(np.sum(v) * dx * dy)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `multislice-projected-potential-integral` | 279.414 | 279.414 | V Å^3 | 5.03e-05 | 1e-03 | ✅ pass |

**Why this value**: 2 pi a0 e x f_e(0) with a0 = 0.529177 Å, e^2/(4 pi eps0) = 14.39965 eV Å (CODATA) and f_e(0) = 2 sum(a_i) = 5.8360 Å from the Lobato-Van Dyck coefficients of Si: 47.8776 x 5.8360 = 279.414 V Å^3.

**Citation**: Kirkland, E. J. (2010). Advanced Computing in Electron Microscopy, 2nd ed., eq. 5.9; Lobato, I. & Van Dyck, D. (2014). Acta Cryst. A70, 636-649.

**See also**: {doc}`Multislice HRTEM theory <../../theory/multislice_hrtem>`, {doc}`Multislice HRTEM algorithm <../../algorithms/multislice_hrtem>`

## Independent-atom mean inner potential of silicon

Averaged over the crystal, the potential of eight silicon atoms per cubic cell of edge a = 5.43102 Å is V0 = (h^2/2 pi m0 e) 8 f_e(0) / a^3 in the independent-atom model. Build a periodic silicon slab along [001], slice it, and report the mean inner potential the sliced potential carries, in volts.

**Symbols**

- $V_{0}$ &mdash; Mean inner potential of a crystal, in volts.
- $f_{e}(s)$ &mdash; Electron atomic scattering factor in ångström.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.app.phases import phase_from_request
from pytex.diffraction.hrem import AtomicSnapshot, MicroscopeAberrations
from pytex.diffraction.multislice import (
    MultisliceGrid,
    SlicedPotential,
    TemporalCoherence,
    hrtem_image,
    multislice,
    parametrized_electron_scattering_factor,
    periodic_slab,
    slice_potential,
)
```

:::

**Compute**

```python
phase = phase_from_request({'builtin': 'si_diamond'})[1]
snap, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1), beam_repeats=3)
grid = MultisliceGrid.from_sampling((snap.cell[0, 0], snap.cell[1, 1]), 0.1)
sliced = SlicedPotential.from_snapshot(snap, grid, phase.lattice.a / 4)
result = sliced.mean_inner_potential_volt()
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `multislice-silicon-mean-inner-potential` | 13.954 | 13.954 | V | 1.56e-04 | 1e-03 | ✅ pass |

**Why this value**: 47.8776 V Å^2 x 8 x 5.8360 Å / (5.43102 Å)^3 = 13.954 V. The independent-atom model ignores bonding, which lowers the measured value to about 12 V: the gap is a known limitation of the model, not of the multislice.

**Citation**: Kirkland (2010), sec. 5.4; Gajdardziska-Josifovska, M. et al. (1993). Ultramicroscopy 50, 285-299 (measured mean inner potentials).

**See also**: {doc}`Multislice HRTEM theory <../../theory/multislice_hrtem>`, {doc}`Multislice HRTEM algorithm <../../algorithms/multislice_hrtem>`

## The largest scattering angle a 0.05 Å grid represents at 200 kV

Multislice band-limits the transmission function and the propagator to two thirds of the Nyquist frequency 1/(2 dx) so the product of wave and transmission function cannot alias. The largest scattering semi-angle the calculation represents is then alpha_max = lambda g_max. Report it in mrad for a 0.05 Å pixel at 200 kV.

**Symbols**

- $g_{\max}$ &mdash; Band limit of the multislice grid, two thirds of Nyquist.
- $\lambda$ &mdash; Relativistic electron wavelength.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.app.phases import phase_from_request
from pytex.diffraction.hrem import AtomicSnapshot, MicroscopeAberrations
from pytex.diffraction.multislice import (
    MultisliceGrid,
    SlicedPotential,
    TemporalCoherence,
    hrtem_image,
    multislice,
    parametrized_electron_scattering_factor,
    periodic_slab,
    slice_potential,
)
```

:::

**Compute**

```python
grid = MultisliceGrid.from_sampling((10.0, 10.0), 0.05)
result = grid.max_scattering_angle_mrad(200.0)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `multislice-band-limit-scattering-angle` | 167.196 | 167.196 | mrad | 3.97e-04 | 1e-03 | ✅ pass |

**Why this value**: alpha_max = lambda (2/3) / (2 dx) = 0.0250793 Å x 6.6667 Å^-1 = 167.196 mrad, with the relativistic wavelength at 200 kV.

**Citation**: Kirkland (2010), sec. 6.8, the two-thirds band-limit rule.

**See also**: {doc}`Multislice HRTEM theory <../../theory/multislice_hrtem>`, {doc}`Multislice HRTEM algorithm <../../algorithms/multislice_hrtem>`

## Imaging at defocus is Fresnel propagation of the exit wave

With no spherical aberration and full coherence, the objective lens multiplies the exit-wave spectrum by exp(-i pi lambda Delta f g^2), which is exactly the free-space propagator over a distance Delta f. A focal series therefore needs one multislice run. Image a silicon [110] slab at Delta f = -150 Å and report the largest difference from |exit wave propagated by -150 Å|^2.

**Symbols**

- $\Delta f$ &mdash; Objective lens defocus; negative is underfocus.
- $\lambda$ &mdash; Relativistic electron wavelength.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.app.phases import phase_from_request
from pytex.diffraction.hrem import AtomicSnapshot, MicroscopeAberrations
from pytex.diffraction.multislice import (
    MultisliceGrid,
    SlicedPotential,
    TemporalCoherence,
    hrtem_image,
    multislice,
    parametrized_electron_scattering_factor,
    periodic_slab,
    slice_potential,
)
```

:::

**Compute**

```python
phase = phase_from_request({'builtin': 'si_diamond'})[1]
snap, _, _ = periodic_slab(phase, (1, 1, 0), (1, 1), thickness_angstrom=30.0)
wave = multislice(snap, 200.0, sampling_angstrom=0.1, slice_thickness_angstrom=1.0)
lens = MicroscopeAberrations(energy_kev=200.0, defocus_angstrom=-150.0, cs_mm=0.0,
                             focal_spread_angstrom=0.0, convergence_semiangle_mrad=0.0)
g2 = wave.grid.frequency_magnitude() ** 2
free = np.fft.ifft2(np.fft.fft2(wave.wave())
                    * np.exp(-1j * math.pi * wave.wavelength_angstrom * -150.0 * g2))
result = float(np.max(np.abs(wave.image(lens) - np.abs(free) ** 2)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `multislice-defocus-is-fresnel-propagation` | < 1e-12 | 0.0e+00 | &mdash; | < 1e-12 | 1e-10 | ✅ pass |

**Why this value**: Identity of the defocus term of the wave aberration, pi lambda Delta f g^2, with the phase of the Fresnel propagator over Delta f: the two images are the same function.

**Citation**: Kirkland (2010), sec. 3.3 and eq. 6.92.

**See also**: {doc}`Multislice HRTEM theory <../../theory/multislice_hrtem>`, {doc}`Multislice HRTEM algorithm <../../algorithms/multislice_hrtem>`

## Multislice and Bloch waves agree on the same crystal potential

Multislice and the Bloch-wave method are two exact solutions of the same high-energy Schrodinger equation. Give both the same Mott-Bethe potential for silicon along [001] at 200 kV, keep only the zero-order Laue zone (the projected potential of one period, here split into eight equal slices so that the splitting error is negligible), propagate through 37 periods (201 Å), and report the largest difference between the two methods' intensities of the transmitted beam, 220 and 400.

**Symbols**

- $\Delta z$ &mdash; Slice thickness of the multislice calculation.
- $v_{n}$ &mdash; Projected potential of slice n, in V Å.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.app.phases import phase_from_request
from pytex.diffraction.hrem import AtomicSnapshot, MicroscopeAberrations
from pytex.diffraction.multislice import (
    MultisliceGrid,
    SlicedPotential,
    TemporalCoherence,
    hrtem_image,
    multislice,
    parametrized_electron_scattering_factor,
    periodic_slab,
    slice_potential,
)
```

:::

**Compute**

```python
import scipy.fft
from pytex.core.lattice import ZoneAxis
from pytex.diffraction.dynamical import beam_set_for_zone, solve_bloch_waves
from pytex.diffraction.hrem import relativistic_interaction_parameter_inv_v_angstrom
from pytex.diffraction.multislice import antialias_aperture, fresnel_propagator
phase = phase_from_request({'builtin': 'si_diamond'})[1]
a = phase.lattice.a
cell, _, _ = periodic_slab(phase, (0, 0, 1), (1, 1))
grid = MultisliceGrid((64, 64), (a, a))
v = slice_potential(cell.species, cell.positions[:, :2], grid, 'mott_bethe')
sigma = relativistic_interaction_parameter_inv_v_angstrom(200.0)
tau = scipy.fft.ifft2(scipy.fft.fft2(np.exp(1j * sigma * v / 8))
                      * antialias_aperture(grid))
propagator = fresnel_propagator(grid, 200.0, a / 8)
psi = np.ones(grid.shape, dtype=complex)
for _ in range(37 * 8):
    psi = scipy.fft.ifft2(scipy.fft.fft2(psi * tau) * propagator)
spectrum = scipy.fft.fft2(psi, norm='forward')
hkl = [(0, 0, 0), (2, 2, 0), (4, 0, 0)]
ms = np.array([abs(spectrum[k, h]) ** 2 for h, k, _ in hkl])
beams = beam_set_for_zone(phase, ZoneAxis(indices=(0, 0, 1), phase=phase),
                          beam_energy_kev=200.0, max_index=24,
                          g_max_inv_angstrom=3.9,
                          max_excitation_error_inv_angstrom=0.6)
bloch = solve_bloch_waves(beams, [[0.0, 0.0]], thickness_angstrom=37 * a)
bw = np.array([float(bloch.intensity_of(h)[0]) for h in hkl])
result = float(np.max(np.abs(ms - bw)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `multislice-agrees-with-bloch-waves` | 0.0004 | 0.0000 | &mdash; | 4.41e-04 | 2e-03 | ✅ pass |

**Why this value**: Two exact solutions of one equation on one potential must agree. The tolerance covers the Bloch solver's finite beam set and its cos(theta) normalisation. With one slice per period instead of eight, the multislice differs by about 1 % at 108 Å: the splitting error of a 5.4 Å slice, which is why slices are kept thin.

**Citation**: Kirkland (2010), ch. 6 (multislice) and ch. 7 (Bloch waves); Self, P. G. et al. (1983). Ultramicroscopy 11, 35-52.

**See also**: {doc}`Multislice HRTEM theory <../../theory/multislice_hrtem>`, {doc}`Multislice HRTEM algorithm <../../algorithms/multislice_hrtem>`

## Integrating the focal spread reproduces Frank's temporal envelope

For a linear (weak) fringe the incoherent average over a Gaussian defocus spread of standard deviation Delta multiplies its contrast by Frank's envelope E_c = exp(-pi^2 lambda^2 Delta^2 g^4 / 2). Image a weak 1 Å^-1 fringe at 300 kV and Delta = 30 Å by explicit Gauss-Hermite focal integration and report the fringe contrast relative to the coherent contrast.

**Symbols**

- $E_{c}$ &mdash; Temporal coherence (focal spread) envelope.
- $\Delta$ &mdash; Focal spread: standard deviation of the defocus distribution.
- $\lambda$ &mdash; Relativistic electron wavelength.


:::{dropdown} Setup (imports and object construction)

```python
import math
import numpy as np
from pytex.app.phases import phase_from_request
from pytex.diffraction.hrem import AtomicSnapshot, MicroscopeAberrations
from pytex.diffraction.multislice import (
    MultisliceGrid,
    SlicedPotential,
    TemporalCoherence,
    hrtem_image,
    multislice,
    parametrized_electron_scattering_factor,
    periodic_slab,
    slice_potential,
)
```

:::

**Compute**

```python
grid = MultisliceGrid((64, 200), (10.0, 6.4))
x = np.arange(200) * 0.05
eps = 1e-4
wave = np.broadcast_to(1.0 + eps * np.cos(2 * math.pi * x), (64, 200))
lens = MicroscopeAberrations(energy_kev=300.0, defocus_angstrom=0.0, cs_mm=0.0,
                             focal_spread_angstrom=30.0, convergence_semiangle_mrad=0.0)
image = hrtem_image(wave, grid, lens,
                    temporal_coherence=TemporalCoherence.FOCAL_INTEGRATION)
fringe = 2 * abs(np.fft.fft2(image)[0, 10]) / image.size
result = float(fringe / (2 * eps))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `multislice-focal-integration-reproduces-frank-envelope` | 0.17881 | 0.17881 | &mdash; | 2.25e-06 | 1e-04 | ✅ pass |

**Why this value**: exp(-0.5 pi^2 lambda^2 Delta^2 g^4) with lambda = 0.0196875 Å (300 kV), Delta = 30 Å, g = 1 Å^-1: exp(-1.72146) = 0.17881.

**Citation**: Frank, J. (1973). The envelope of electron microscopic transfer functions for partially coherent illumination. Optik 38, 519-536.

**See also**: {doc}`Multislice HRTEM theory <../../theory/multislice_hrtem>`, {doc}`Multislice HRTEM algorithm <../../algorithms/multislice_hrtem>`
