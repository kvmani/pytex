<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Phase identification from a powder pattern

A pattern generated from a known fixture ranked against three candidates chosen to be wrong in three different ways; a cell dilation imposed by the example and recovered by the refinement; the algebraic identity that makes that refinement safe rather than a way of flattering any candidate; and the contract a comparison of several uploaded structures depends on, that an impossible candidate is scored rather than dropped.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## A pattern generated from nickel ranks nickel first

The elementary check on any identification: generate a powder pattern from a known structure, hand the ranking that structure along with three plausible competitors, and require that it comes back first. The answer is fixed by the fixture the pattern was generated from, before the calculation starts. The competitors are not chosen to be easy: copper is face-centred cubic like nickel and differs only in cell size, ferrite is cubic of a similar size but body-centred, and halite is face-centred cubic with a two-species basis. Between them they exercise all three ways a candidate can be wrong - the wrong cell dimension, the wrong centring, and the wrong basis. The reported value is the rank of the true phase, which must be 1.

**Symbols**

- $\theta$ &mdash; Bragg half-angle.
- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_phase_identification import identify_phase_from_pattern

def scan(phase, seed=7):
    radiation = RadiationSpec.cu_ka()
    pattern = generate_xrd_pattern(
        phase,
        radiation=radiation,
        two_theta_range_deg=(25.0, 140.0),
        resolution_deg=0.01,
        broadening_fwhm_deg=0.12,
        profile='pseudo_voigt',
    )
    profile = np.asarray(pattern.intensity_grid, dtype=float)
    counts = profile / profile.max() * 30000.0 + 150.0
    return MeasuredPowderPattern(
        name='synthetic specimen',
        two_theta_deg=np.asarray(pattern.two_theta_grid_deg, dtype=float),
        intensity=np.random.default_rng(seed).poisson(counts).astype(float),
        radiation=radiation,
        synthetic=True,
    )
```

:::

**Compute**

```python
truth = builtin_phase('ni_fcc').to_phase()
candidates = {
    key: builtin_phase(key).to_phase()
    for key in ('ni_fcc', 'cu_fcc', 'fe_bcc', 'nacl')
}
report, _ = identify_phase_from_pattern(scan(truth), candidates)
result = float(
    1 + [item.phase_name for item in report].index('ni_fcc')
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `phase-id-ranking-returns-the-generating-phase` | 1 | 1 | &mdash; | exact | exact | ✅ pass |

**Why this value**: The pattern was generated from the pinned ni_fcc fixture, so the identity of the specimen is known independently of the ranking. A correct ranking places it first.

**Citation**: Hanawalt, Rinn & Frevel, Ind. Eng. Chem. Anal. Ed. 10 (1938) 457, doi:10.1021/ac50125a001.

**See also**: {doc}`Phase identification from powder patterns <../../theory/phase_identification_from_powder_patterns>`, {doc}`Phase identification algorithm <../../algorithms/phase_identification>`, {doc}`Powder XRD and SAED theory <../../theory/powder_xrd_and_saed>`

## The refined cell dilation recovers a dilation the example imposed

A CIF records the cell of somebody else's specimen. Yours is a solid solution, or at another temperature, or stressed, and its cell differs by a fraction of a per cent - which by Delta(2*theta) = 2 e tan(theta) displaces a back-reflection line by far more than any sensible matching tolerance. The ranking therefore refines one uniform cell dilation per candidate before indexing, and reports it.

Here the specimen's cell is dilated by exactly 1.0040 relative to the tabulated nickel fixture, and the tabulated fixture is offered as the candidate. The refinement must recover the factor that was imposed. Note what this also demonstrates: the *candidate* is the undilated tabulated cell, exactly as it would arrive from a CIF.

**Symbols**

- $\theta$ &mdash; Bragg half-angle.
- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_phase_identification import identify_phase_from_pattern

def scan(phase, seed=7):
    radiation = RadiationSpec.cu_ka()
    pattern = generate_xrd_pattern(
        phase,
        radiation=radiation,
        two_theta_range_deg=(25.0, 140.0),
        resolution_deg=0.01,
        broadening_fwhm_deg=0.12,
        profile='pseudo_voigt',
    )
    profile = np.asarray(pattern.intensity_grid, dtype=float)
    counts = profile / profile.max() * 30000.0 + 150.0
    return MeasuredPowderPattern(
        name='synthetic specimen',
        two_theta_deg=np.asarray(pattern.two_theta_grid_deg, dtype=float),
        intensity=np.random.default_rng(seed).poisson(counts).astype(float),
        radiation=radiation,
        synthetic=True,
    )
from dataclasses import replace

def dilated(phase, scale, name):
    lattice = replace(
        phase.lattice,
        a=phase.lattice.a * scale,
        b=phase.lattice.b * scale,
        c=phase.lattice.c * scale,
    )
    return replace(
        phase,
        lattice=lattice,
        unit_cell=replace(phase.unit_cell, lattice=lattice),
        name=name,
    )
```

:::

**Compute**

```python
tabulated = builtin_phase('ni_fcc').to_phase()
specimen = dilated(tabulated, 1.0040, 'nickel solid solution')
report, _ = identify_phase_from_pattern(
    scan(specimen), {'tabulated nickel': tabulated}
)
result = float(report.best.cell_scale)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `phase-id-refinement-recovers-an-imposed-cell-dilation` | 1.0040 | 1.0040 | &mdash; | < 5e-06 | 5e-04 | ✅ pass |

**Why this value**: The dilation is imposed by this example, so the value to recover is set before the calculation runs. The tolerance is the resolution of the scale grid searched over the default two per cent range, which is 1.0e-4, widened to 5.0e-4 to absorb the counting noise on the fitted peak positions.

**Citation**: Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 11 - Delta d / d = -cot(theta) Delta(theta), the relation that makes a small cell error a large high-angle displacement.

**See also**: {doc}`Phase identification from powder patterns <../../theory/phase_identification_from_powder_patterns>`, {doc}`Phase identification algorithm <../../algorithms/phase_identification>`

## A uniform dilation leaves every ratio of d spacings unchanged

This is the algebra that makes the cell-scale refinement safe rather than a way of flattering any candidate at all. For a lattice scaled uniformly by s, every interplanar spacing becomes s*d_hkl, so the ratio d_hkl / d_h'k'l' is unchanged exactly, for every pair, whatever the structure. Those ratios are precisely what indexing tests: a candidate whose relative line positions are wrong is wrong at every scale, and one a scale factor can rescue is the right structure with the wrong cell size.

The check dilates nickel by 1.05 - far beyond anything the refinement would search - and reports the largest absolute deviation between the two sets of spacing ratios. The families are matched by their Miller indices rather than by position in the list, because a dilated cell moves every line to lower angle and so pulls an extra one into a fixed angular window; comparing the two lists elementwise would be comparing different reflections. The deviation must be zero to machine precision.

**Symbols**

- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_powder_reflections
from dataclasses import replace

def dilated(phase, scale, name):
    lattice = replace(
        phase.lattice,
        a=phase.lattice.a * scale,
        b=phase.lattice.b * scale,
        c=phase.lattice.c * scale,
    )
    return replace(
        phase,
        lattice=lattice,
        unit_cell=replace(phase.unit_cell, lattice=lattice),
        name=name,
    )
```

:::

**Compute**

```python
radiation = RadiationSpec.cu_ka()
phase = builtin_phase('ni_fcc').to_phase()

def spacings(target):
    lines = generate_powder_reflections(
        target, radiation=radiation, two_theta_range_deg=(20.0, 150.0)
    )
    return {
        tuple(int(value) for value in item.miller_indices): item.d_spacing_angstrom
        for item in lines
    }

plain = spacings(phase)
stretched = spacings(dilated(phase, 1.05, 'stretched'))
shared = sorted(set(plain) & set(stretched))
reference = shared[0]
result = float(
    max(
        abs(plain[hkl] / plain[reference] - stretched[hkl] / stretched[reference])
        for hkl in shared
    )
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `phase-id-uniform-dilation-preserves-spacing-ratios` | < 1e-12 | 0.0e+00 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: Algebra: d_hkl -> s d_hkl under a uniform dilation, so every ratio d_hkl / d_h'k'l' is identically unchanged. The tolerance is floating-point round-off, not a physical allowance.

**Citation**: Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 2 - the plane-spacing equations, in which the cell edges enter as an overall scale.

**See also**: {doc}`Phase identification from powder patterns <../../theory/phase_identification_from_powder_patterns>`, {doc}`Phase identification algorithm <../../algorithms/phase_identification>`

## A candidate that cannot be indexed is scored, not dropped

A user who opens five CIF files and gets back four rows cannot tell which one was discarded or why, and a comparison that aborts because one candidate is impossible wastes the other four. So a candidate whose lines cannot be matched at all - here a cell shrunk until every d spacing falls below lambda/2, for which Bragg's law has no solution at any angle - is recorded with a stated reason and a score of zero, and the ranking of the rest proceeds.

Five candidates are offered, one of them impossible. The reported value is the number of rows returned, which must be five.

**Symbols**

- $\theta$ &mdash; Bragg half-angle.
- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_phase_identification import identify_phase_from_pattern

def scan(phase, seed=7):
    radiation = RadiationSpec.cu_ka()
    pattern = generate_xrd_pattern(
        phase,
        radiation=radiation,
        two_theta_range_deg=(25.0, 140.0),
        resolution_deg=0.01,
        broadening_fwhm_deg=0.12,
        profile='pseudo_voigt',
    )
    profile = np.asarray(pattern.intensity_grid, dtype=float)
    counts = profile / profile.max() * 30000.0 + 150.0
    return MeasuredPowderPattern(
        name='synthetic specimen',
        two_theta_deg=np.asarray(pattern.two_theta_grid_deg, dtype=float),
        intensity=np.random.default_rng(seed).poisson(counts).astype(float),
        radiation=radiation,
        synthetic=True,
    )
from dataclasses import replace

def dilated(phase, scale, name):
    lattice = replace(
        phase.lattice,
        a=phase.lattice.a * scale,
        b=phase.lattice.b * scale,
        c=phase.lattice.c * scale,
    )
    return replace(
        phase,
        lattice=lattice,
        unit_cell=replace(phase.unit_cell, lattice=lattice),
        name=name,
    )
```

:::

**Compute**

```python
truth = builtin_phase('ni_fcc').to_phase()
candidates = {
    key: builtin_phase(key).to_phase()
    for key in ('ni_fcc', 'cu_fcc', 'fe_bcc', 'nacl')
}
candidates['impossible'] = dilated(truth, 0.30, 'impossible')
report, _ = identify_phase_from_pattern(scan(truth), candidates)
rejected = [item for item in report if item.indexing is None]
result = float(len(report) if len(rejected) == 1 else -1)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `phase-id-every-candidate-offered-is-ranked` | 5 | 5 | &mdash; | exact | exact | ✅ pass |

**Why this value**: Five candidates are offered and exactly one of them - the cell shrunk by 0.30, whose largest d spacing is below lambda/2 for Cu K-alpha - can produce no reflection at any angle. A ranking that returned four rows, or that raised, would fail this check; the result is set to -1 unless exactly one candidate was rejected, so a run that indexed the impossible cell would fail too.

**Citation**: Cullity & Stock, Elements of X-Ray Diffraction, 3rd ed., Prentice Hall (2001), Ch. 3 - lambda = 2 d sin(theta) has no solution for d < lambda/2.

**See also**: {doc}`Phase identification from powder patterns <../../theory/phase_identification_from_powder_patterns>`, {doc}`Phase identification algorithm <../../algorithms/phase_identification>`
