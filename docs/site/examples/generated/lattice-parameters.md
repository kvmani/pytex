<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Precise lattice-parameter determination

The extrapolation functions and the drift-column identity checked against closed-form algebra, the hexagonal quadratic form checked against the textbook expression it must reduce to, and an end-to-end determination through a deliberately displaced specimen checked against the pinned fixture cell the pattern was generated from - alongside the same data solved by averaging, which leaves the displacement in the answer.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## The Nelson-Riley function at 2*theta = 120 degrees

Before trusting an extrapolation you should be able to evaluate its function by hand. Nelson and Riley (1945) define f = (cos^2(theta)/sin(theta) + cos^2(theta)/theta)/2 with theta in radians. At theta = 60 degrees = pi/3 radians, cos(theta) = 1/2 and sin(theta) = sqrt(3)/2, so the two terms are (1/4)/(sqrt(3)/2) = 0.2886751 and (1/4)/(pi/3) = 0.2387324, and their mean is 0.2637038.

**Symbols**

- $\theta$ &mdash; Bragg half-angle.
- $f(\theta)$ &mdash; Extrapolation function of the systematic-error term.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.diffraction.xrd_lattice_parameter import extrapolation_values
```

:::

**Compute**

```python
result = float(
    extrapolation_values([120.0], function='nelson_riley')[0]
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `lattice-nelson-riley-at-sixty-degrees` | 0.2637038 | 0.2637038 | &mdash; | 2.54e-08 | 1e-06 | ✅ pass |

**Why this value**: By hand at theta = pi/3: cos^2(theta)/sin(theta) = 0.25 / 0.8660254 = 0.2886751; cos^2(theta)/theta = 0.25 / 1.0471976 = 0.2387324; the mean is 0.2637038.

**Citation**: Nelson & Riley, Proc. Phys. Soc. 57 (1945) 160, doi:10.1088/0959-5309/57/3/302.

**See also**: {doc}`Precise lattice-parameter determination <../../theory/precise_lattice_parameter_determination>`

## Every extrapolation function vanishes at theta = 90 degrees

The reason extrapolation removes an aberration is that every admissible function tends to zero at back-reflection, where the geometric aberrations themselves vanish. At theta = 90 degrees, cos(theta) = 0 exactly, so cos^2(theta)/sin(theta), cos^2(theta) and cos(theta)/sin(theta) are all zero, and the Nelson-Riley mean of the first with cos^2(theta)/theta is zero too. Summing all four at 2*theta = 180 degrees must give exactly zero.

**Symbols**

- $\theta$ &mdash; Bragg half-angle.
- $f(\theta)$ &mdash; Extrapolation function of the systematic-error term.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.diffraction.xrd_lattice_parameter import (
    EXTRAPOLATION_FUNCTIONS,
    extrapolation_values,
)
```

:::

**Compute**

```python
result = float(
    sum(
        abs(extrapolation_values([180.0], function=name)[0])
        for name in EXTRAPOLATION_FUNCTIONS
    )
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `lattice-extrapolation-vanishes-at-backscatter` | < 1e-12 | 0.00e+00 | &mdash; | < 1e-12 | 1e-10 | ✅ pass |

**Why this value**: cos(90 degrees) = 0, so every function built from cos(theta) is zero there, and the sum of their absolute values is zero to floating-point precision. Note that theta = 90 degrees is not a singular point for any of them: only sin(theta) appears in a denominator, and sin(90 degrees) = 1.

**Citation**: Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 11.

**See also**: {doc}`Precise lattice-parameter determination <../../theory/precise_lattice_parameter_determination>`

## The Bradley-Jay drift column is Cohen's sin^2(2*theta)/4

The identity that unifies the classical graphical extrapolation with the least-squares treatment: the design column of the systematic-error term is sin^2(theta) f(theta), with the same f the graphical method plots a against. Taking the Bradley-Jay function f = cos^2(theta) must therefore reproduce Cohen's classical drift column exactly, because sin^2(theta) cos^2(theta) = (2 sin(theta) cos(theta))^2 / 4 = sin^2(2 theta) / 4. Checked here as the largest absolute difference over a wide angular range, which must be zero to machine precision.

**Symbols**

- $\theta$ &mdash; Bragg half-angle.
- $f(\theta)$ &mdash; Extrapolation function of the systematic-error term.
- $D$ &mdash; Systematic-error (drift) coefficient refined with the cell.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.diffraction.xrd_lattice_parameter import extrapolation_values

angles = np.linspace(10.0, 170.0, 321)
theta = np.deg2rad(0.5 * angles)
```

:::

**Compute**

```python
column = np.square(np.sin(theta)) * extrapolation_values(
    angles, function='bradley_jay'
)
cohen = np.square(np.sin(2.0 * theta)) / 4.0
result = float(np.max(np.abs(column - cohen)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `lattice-cohen-drift-column-identity` | < 1e-15 | 0.00e+00 | &mdash; | < 1e-15 | 1e-15 | ✅ pass |

**Why this value**: sin^2(theta) cos^2(theta) = (sin(2 theta) / 2)^2 = sin^2(2 theta) / 4 is a double-angle identity, so the two columns are the same function and their difference is zero.

**Citation**: Cohen, Rev. Sci. Instrum. 6 (1935) 68, doi:10.1063/1.1751937.

**See also**: {doc}`Precise lattice-parameter determination <../../theory/precise_lattice_parameter_determination>`

## The hexagonal quadratic form falls out of the metric tensor

PyTex determines a cell by least squares on the components of the reciprocal metric tensor, with the crystal system entering as a constraint rather than as per-system algebra. For a hexagonal cell a* = b* and gamma* = 60 degrees, so G*_12 = a* b* cos(gamma*) = G*_11 / 2, and the quadratic form h^T G* h must collapse to the textbook 1/d^2 = (4/3)(h^2 + hk + k^2)/a^2 + l^2/c^2. Checked on zirconium (a = 3.232, c = 5.147 angstrom) for the (2 1 3) reflection, where the textbook value is (4/3)(4 + 2 + 1)/3.232^2 + 9/5.147^2 = 0.8935798 + 0.3396494 = 1.2332292 inverse square angstrom.

**Symbols**

- $\mathbf{G}^{*}$ &mdash; Reciprocal metric tensor.
- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import FrameDomain, Handedness, Lattice, ReferenceFrame

crystal = ReferenceFrame(
    name='crystal',
    domain=FrameDomain.CRYSTAL,
    axes=('a', 'b', 'c'),
    handedness=Handedness.RIGHT,
)
# Zirconium (HCP), lattice parameters from the pinned PyTex fixture corpus.
zirconium = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal)
```

:::

**Compute**

```python
reciprocal = zirconium.reciprocal_metric_tensor()
indices = np.array([2.0, 1.0, 3.0])
result = float(indices @ reciprocal @ indices)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `lattice-hexagonal-quadratic-form` | 1.2332292 | 1.2332292 | 1/angstrom^2 | 1.76e-08 | 1e-06 | ✅ pass |

**Why this value**: 1/d^2 = (4/3)(h^2 + hk + k^2)/a^2 + l^2/c^2 with (hkl) = (213), a = 3.232 and c = 5.147 angstrom: (4/3)(7)/10.445824 + 9/26.491609 = 0.8935798 + 0.3396494 = 1.2332292 inverse square angstrom.

**Citation**: Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Appendix; International Tables for Crystallography Vol. A.

**See also**: {doc}`Precise lattice-parameter determination <../../theory/precise_lattice_parameter_determination>`, {doc}`Powder XRD and SAED theory <../../theory/powder_xrd_and_saed>`

## Cohen least squares recovers the nickel cell through a displaced specimen

The end-to-end determination, on a pattern generated from the pinned nickel fixture and then displaced by a known 100 micrometres. The peaks are detected, fitted, indexed and solved by weighted least squares in sin^2(theta) with a drift term against cos^2(theta)/sin(theta) - the exact angular form of a specimen displacement. The answer must be the fixture's own lattice parameter, 3.52387 angstrom, because that is the cell the pattern was generated from.

**Symbols**

- $a$ &mdash; Cubic or hexagonal lattice parameter.
- $f(\theta)$ &mdash; Extrapolation function of the systematic-error term.
- $D$ &mdash; Systematic-error (drift) coefficient refined with the cell.
- $\mathbf{G}^{*}$ &mdash; Reciprocal metric tensor.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import FrameDomain, Handedness, ReferenceFrame
from pytex.core.fixtures import get_phase_fixture
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_corrections import specimen_displacement_shift_deg
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_lattice_parameter import (
    determine_lattice_parameters_from_pattern,
)
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC) from the pinned fixture corpus: a = 3.52387 angstrom.
nickel = get_phase_fixture("ni_fcc").load_phase(crystal_frame=crystal)
radiation = RadiationSpec.cu_ka_doublet()

pattern = generate_xrd_pattern(
    nickel,
    radiation=radiation,
    two_theta_range_deg=(25.0, 150.0),
    resolution_deg=0.01,
    broadening_fwhm_deg=0.12,
    profile="pseudo_voigt",
    max_index=6,
)
axis = np.asarray(pattern.two_theta_grid_deg)
counts = np.random.default_rng(5).poisson(
    np.asarray(pattern.intensity_grid) / pattern.intensity_grid.max() * 30000.0 + 150.0
).astype(float)
# A 100 micrometre specimen displacement, which is an ordinary preparation
# error, injected so the determination has a known aberration to remove.
displaced = MeasuredPowderPattern(
    name="nickel, specimen 100 um off axis",
    two_theta_deg=axis
    + specimen_displacement_shift_deg(
        axis, displacement_mm=0.10, goniometer_radius_mm=240.0
    ),
    intensity=counts,
    radiation=radiation,
    synthetic=True,
)
instrument = InstrumentBroadening.ideal(0.12)
```

:::

**Compute**

```python
determination, _ = determine_lattice_parameters_from_pattern(
    displaced,
    nickel,
    method='cohen',
    extrapolation='cos_squared_over_sin',
    instrument=instrument,
)
result = float(determination.a)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `lattice-cohen-recovers-the-nickel-cell` | 3.52387 | 3.52387 | angstrom | 3.38e-07 | 2e-05 | ✅ pass |

**Why this value**: The pinned fixture fixtures/phases/ni_fcc records a = 3.52387 angstrom (COD 9008476, Wyckoff, Crystal Structures Vol. 1). The pattern was generated from that cell, so the determination must return it once the injected displacement is extrapolated away.

**Citation**: Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 11; Cohen, Rev. Sci. Instrum. 6 (1935) 68, doi:10.1063/1.1751937.

**See also**: {doc}`Precise lattice-parameter determination <../../theory/precise_lattice_parameter_determination>`, {doc}`Powder XRD and SAED theory <../../theory/powder_xrd_and_saed>`

## Averaging over reflections leaves the displacement in the answer

The same peaks, solved by the intuitive method: compute a from each reflection and take the mean. Because Delta d / d = -cot(theta) Delta theta makes a fixed angular error a theta-dependent spacing error, averaging divides the random scatter by sqrt(N) and leaves the systematic part untouched. The relative error is therefore of order the fractional spacing error the displacement produces - about 4 parts in 10^4 - rather than the 1 part in 10^7 the matched extrapolation reaches. Reported as the ratio of the two errors, which must be large.

**Symbols**

- $a$ &mdash; Cubic or hexagonal lattice parameter.
- $\theta$ &mdash; Bragg half-angle.
- $f(\theta)$ &mdash; Extrapolation function of the systematic-error term.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import FrameDomain, Handedness, ReferenceFrame
from pytex.core.fixtures import get_phase_fixture
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_corrections import specimen_displacement_shift_deg
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_lattice_parameter import (
    determine_lattice_parameters_from_pattern,
)
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC) from the pinned fixture corpus: a = 3.52387 angstrom.
nickel = get_phase_fixture("ni_fcc").load_phase(crystal_frame=crystal)
radiation = RadiationSpec.cu_ka_doublet()

pattern = generate_xrd_pattern(
    nickel,
    radiation=radiation,
    two_theta_range_deg=(25.0, 150.0),
    resolution_deg=0.01,
    broadening_fwhm_deg=0.12,
    profile="pseudo_voigt",
    max_index=6,
)
axis = np.asarray(pattern.two_theta_grid_deg)
counts = np.random.default_rng(5).poisson(
    np.asarray(pattern.intensity_grid) / pattern.intensity_grid.max() * 30000.0 + 150.0
).astype(float)
# A 100 micrometre specimen displacement, which is an ordinary preparation
# error, injected so the determination has a known aberration to remove.
displaced = MeasuredPowderPattern(
    name="nickel, specimen 100 um off axis",
    two_theta_deg=axis
    + specimen_displacement_shift_deg(
        axis, displacement_mm=0.10, goniometer_radius_mm=240.0
    ),
    intensity=counts,
    radiation=radiation,
    synthetic=True,
)
instrument = InstrumentBroadening.ideal(0.12)
```

:::

**Compute**

```python
naive, _ = determine_lattice_parameters_from_pattern(
    displaced, nickel, method='average', extrapolation='none',
    instrument=instrument,
)
matched, _ = determine_lattice_parameters_from_pattern(
    displaced, nickel, method='cohen',
    extrapolation='cos_squared_over_sin', instrument=instrument,
)
truth = nickel.lattice.a
result = float(abs(naive.a - truth) / abs(matched.a - truth))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `lattice-averaging-fails-by-a-measured-amount` | 4173 | 4000 | &mdash; | 1.73e+02 | 4e+03 | ✅ pass |

**Why this value**: A 100 micrometre displacement at R = 240 mm shifts the lowest reflection by about 47 millidegrees. At 2*theta = 44.5 degrees that is Delta d / d = cot(theta) Delta theta = 2.46 x (8.1 x 10^-4) / 2 ~ 4 x 10^-4 in the averaged answer, against a matched extrapolation limited only by counting statistics near 10^-7. The ratio is therefore of order 10^3 to 10^4; the tolerance spans that range because the numerator is set by physics and the denominator by the noise realisation.

**Citation**: Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 11 (Precise Parameter Measurements).

**See also**: {doc}`Precise lattice-parameter determination <../../theory/precise_lattice_parameter_determination>`
