<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# The ghost problem, and its correction

What diffraction pole figures cannot determine, and what positivity can recover of it: an asymmetric texture still gives a pole set closed under negation, excluding the odd harmonic degrees discards nearly half the basis, a cubic material has no odd term below degree 9, and the correction removes the negative density without moving the fit.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## An asymmetric texture still gives a centrosymmetric pole figure

Build a deliberately one-sided orientation population - 200 orientations whose first Euler angle is confined to 0-40 degrees, with nothing symmetric about it - and generate its {111} pole figure. Every pole in the result has its antipode also present. The centrosymmetry is a property of the measurement, not of the specimen: a lattice plane has no sense and Friedel's law makes +h and -h scatter identically. This is the root cause of the ghost problem, and the reason no amount of pole-figure data can recover the odd part of an ODF. The example returns the fraction of the first 300 poles whose antipode is in the set, which must be 1.

**Symbols**

- $P_{\mathbf{h}}(\mathbf{y})$ &mdash; Pole density of plane family h along specimen direction y.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    Lattice,
    MillerIndex,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.texture.models import PoleFigure

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
angles = np.linspace(0.0, 40.0, 200)
orientations = OrientationSet.from_euler_angles(
    np.column_stack(
        [angles, 25.0 * np.ones_like(angles), np.zeros_like(angles)]
    ),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)
figure = PoleFigure.from_orientations(orientations, pole)
directions = figure.sample_directions
poles = np.asarray(
    getattr(directions, 'values', directions)
)
paired = sum(
    1
    for row in poles[:300]
    if np.min(np.linalg.norm(poles + row, axis=1)) < 1e-9
)
result = paired / 300.0
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ghost-pole-figure-is-centrosymmetric` | 1.000000 | 1.000000 | &mdash; | exact | exact | ✅ pass |

**Why this value**: Analytic: a plane normal enters a pole figure as an axis and Friedel's law equates +h with -h, so the pole set is closed under negation for any orientation distribution whatsoever.

**Citation**: Matthies, On the reproducibility of the orientation distribution function of texture samples from pole figures (ghost phenomena), Phys. Status Solidi B 92 (1979) K135-K138.

**See also**: {doc}`The ghost problem <../../theory/ghost_problem_and_odd_harmonics>`, {doc}`Harmonic ODF reconstruction <../../theory/harmonic_odf_reconstruction>`

## Excluding odd degrees discards nearly half the harmonic basis

Count the generalized spherical harmonic terms retained and discarded when odd degrees are excluded, at bandlimit 22. Because pole figures annihilate every odd degree, those coefficients are not poorly determined but wholly undetermined, and the count says how much of the ODF a diffraction measurement is silent about: 15147 of 32407 terms, or 46.7 percent, tending to one half as the bandlimit grows. The example returns the discarded fraction.

**Symbols**

- $\ell$ &mdash; Degree of a generalized spherical harmonic term.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.texture.harmonics import _enumerate_terms
```

:::

**Compute**

```python
even = len(
    _enumerate_terms(degree_bandlimit=22, even_degrees_only=True)
)
every = len(
    _enumerate_terms(degree_bandlimit=22, even_degrees_only=False)
)
result = (every - even) / every
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ghost-odd-degrees-are-half-the-harmonic-basis` | 0.46740 | 0.46739 | &mdash; | 9.02e-06 | 1e-04 | ✅ pass |

**Why this value**: Exact term count: sum over degrees 0..22 of (2l+1)^2 is 32407, of which the even degrees contribute 17260, leaving 15147 discarded. The fraction tends to 1/2 as the bandlimit grows.

**Citation**: Bunge, Texture Analysis in Materials Science: Mathematical Methods (Butterworths 1969) - the generalized spherical harmonic expansion.

**See also**: {doc}`The ghost problem <../../theory/ghost_problem_and_odd_harmonics>`, {doc}`Harmonic ODF reconstruction <../../theory/harmonic_odf_reconstruction>`

## A cubic material has no odd harmonic to correct below degree 9

Ghost correction supplies the odd part of an ODF that a pole figure cannot measure. Whether there is an odd part to supply at all is a question about the symmetry, not about the data: the crystal rotation group admits odd-degree terms only where it has an odd-degree invariant. Counting those invariants by character theory - the group average of the SO(3) character - gives the classical answer for the cubic rotation group 432: nothing at degrees 1, 3, 5 or 7, and the first odd invariant at degree 9. A cubic ODF expanded to degree 6 or 8 therefore has no ghost part, and PyTex's correction reports that rather than a correction of size zero. Lower symmetries admit odd terms much earlier - degree 7 for hexagonal 622, degree 3 for orthorhombic 222.

**Symbols**

- $\ell$ &mdash; Degree of a generalized spherical harmonic term.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import FrameDomain, ReferenceFrame, SymmetrySpec

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
cubic = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)


def invariant_count(symmetry, degree):
    """Dimension of the degree-l invariant subspace of a rotation group.

    Character theory: the number of invariants is the group average of the
    character of the degree-l representation of SO(3),
    chi_l(theta) = sin((l + 1/2) theta) / sin(theta / 2).
    """

    operators = np.asarray(symmetry.operators)
    cosines = np.clip((np.trace(operators, axis1=1, axis2=2) - 1.0) / 2.0, -1.0, 1.0)
    theta = np.arccos(cosines)
    half = np.sin(theta / 2.0)
    identity = np.abs(half) < 1e-12
    chi = np.where(
        identity,
        2.0 * degree + 1.0,
        np.sin((degree + 0.5) * theta) / np.where(identity, 1.0, half),
    )
    return float(np.mean(chi))
```

:::

**Compute**

```python
odd_degrees = range(1, 16, 2)
result = min(
    degree for degree in odd_degrees if invariant_count(cubic, degree) > 0.5
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ghost-cubic-first-odd-invariant-is-degree-nine` | 9 | 9 | &mdash; | exact | exact | ✅ pass |

**Why this value**: Standard result of cubic harmonic analysis: the cubic rotation group has invariants at degrees 0, 4, 6, 8, 9, 10, ... and the lowest odd-degree cubic harmonic is degree 9. Independently derivable from the Molien series of the octahedral group.

**Citation**: Bunge, Texture Analysis in Materials Science: Mathematical Methods (Butterworths 1969), tables of cubic symmetric generalized spherical harmonics.

**See also**: {doc}`The ghost problem <../../theory/ghost_problem_and_odd_harmonics>`, {doc}`Harmonic ODF reconstruction <../../theory/harmonic_odf_reconstruction>`

## Ghost correction removes the negative density the even part leaves

The even-degree half of an ODF is all a pole figure determines, and on its own it is not a density: here it falls to about -0.46 multiples of random over part of orientation space, which is the classical ghost artefact. Correction holds that even part fixed - it is what the data say - and adds the smallest odd part that makes the whole non-negative. Because the added part is the smallest one that works, the constraint ends up exactly tight: the minimum density of the corrected distribution sits on zero rather than comfortably above it. The maximum rises at the same time, from 3.75 to 4.22 m.r.d., which is the other half of the ghost signature: the even-only solution pays for its false lobes by depressing the true peak.

**Symbols**

- $\ell$ &mdash; Degree of a generalized spherical harmonic term.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    HarmonicODF,
    KernelSpec,
    Lattice,
    MillerIndex,
    ODF,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    SymmetrySpec,
    random_pole_density,
)
from pytex.diffraction.stereonets import spherical_angles_to_directions

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
# Orthorhombic, because 222 admits odd-degree terms from degree 3; a cubic
# material has no ghost part below degree 9 and nothing to demonstrate here.
symmetry = SymmetrySpec.from_point_group("222", reference_frame=crystal)
phase = Phase(
    name="orthorhombic-demo",
    lattice=Lattice(3.0, 4.0, 5.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
# A single broad component: broad enough that a degree-4 expansion represents
# it, so the demonstration is of the ghost problem and not of truncation.
truth = ODF(
    orientations=OrientationSet.from_euler_angles(
        np.array([[35.0, 45.0, 20.0]]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    ),
    weights=np.array([1.0]),
    kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=50.0),
)
polar, azimuth = np.meshgrid(
    np.arange(0.0, 91.0, 15.0), np.arange(0.0, 360.0, 15.0), indexing="ij"
)
directions = spherical_angles_to_directions(polar, azimuth).reshape(-1, 3)
scale = random_pole_density(truth.kernel, antipodal=True)


def measured(indices):
    pole = CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase)
    return PoleFigure(
        pole=pole,
        sample_directions=directions,
        intensities=truth.evaluate_pole_density(pole, directions, antipodal=True) / scale,
        specimen_frame=specimen,
        antipodal=True,
        includes_symmetry_family=True,
        sampling="sampled_density",
    )


pole_figures = [measured(indices) for indices in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0])]
report = HarmonicODF.invert_pole_figures(
    pole_figures,
    degree_bandlimit=4,
    regularization=1e-6,
    pole_kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=10.0),
    phi1_step_deg=15.0,
    big_phi_step_deg=15.0,
    phi2_step_deg=15.0,
    ghost_correction=True,
)
correction = report.ghost_correction
```

:::

**Compute**

```python
result = float(np.min(correction.odf.quadrature_densities))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ghost-correction-restores-a-non-negative-density` | -0.0000 | 0.0000 | m.r.d. | 3.52e-05 | 1e-03 | ✅ pass |

**Why this value**: Analytic: an orientation distribution is a probability density and cannot be negative, and the correction minimizes the norm of the odd part subject to that constraint. A minimum-norm feasible point lies on the boundary of the feasible set whenever the unconstrained point is infeasible, so the minimum density is zero to within the quadrature resolution.

**Citation**: Dahms and Bunge, The iterative series-expansion method for quantitative texture analysis. I. General outline, J. Appl. Cryst. 22 (1989) 439-447.

**See also**: {doc}`The ghost problem <../../theory/ghost_problem_and_odd_harmonics>`, {doc}`Harmonic ODF reconstruction <../../theory/harmonic_odf_reconstruction>`

## The odd part a correction adds is invisible to the pole figures

A correction that improved the density by moving the fit would be spending data agreement it has no right to spend, so this is the check that matters. Under Friedel's law the forward operator annihilates every odd-degree harmonic exactly, and adding one therefore leaves every predicted pole density where it was. The example returns the largest change over all four figures and every measured direction, against intensities that reach 2.4 m.r.d.; it is nonzero only because the odd basis is orthonormalized on a discrete quadrature, so its orthogonality to the even part is exact only in the continuum. Refining the quadrature drives it to zero.

**Symbols**

- $P_{\mathbf{h}}(\mathbf{y})$ &mdash; Pole density of plane family h along specimen direction y.
- $\ell$ &mdash; Degree of a generalized spherical harmonic term.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CrystalPlane,
    FrameDomain,
    HarmonicODF,
    KernelSpec,
    Lattice,
    MillerIndex,
    ODF,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    SymmetrySpec,
    random_pole_density,
)
from pytex.diffraction.stereonets import spherical_angles_to_directions

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"))
# Orthorhombic, because 222 admits odd-degree terms from degree 3; a cubic
# material has no ghost part below degree 9 and nothing to demonstrate here.
symmetry = SymmetrySpec.from_point_group("222", reference_frame=crystal)
phase = Phase(
    name="orthorhombic-demo",
    lattice=Lattice(3.0, 4.0, 5.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=symmetry,
    crystal_frame=crystal,
)
# A single broad component: broad enough that a degree-4 expansion represents
# it, so the demonstration is of the ghost problem and not of truncation.
truth = ODF(
    orientations=OrientationSet.from_euler_angles(
        np.array([[35.0, 45.0, 20.0]]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    ),
    weights=np.array([1.0]),
    kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=50.0),
)
polar, azimuth = np.meshgrid(
    np.arange(0.0, 91.0, 15.0), np.arange(0.0, 360.0, 15.0), indexing="ij"
)
directions = spherical_angles_to_directions(polar, azimuth).reshape(-1, 3)
scale = random_pole_density(truth.kernel, antipodal=True)


def measured(indices):
    pole = CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase)
    return PoleFigure(
        pole=pole,
        sample_directions=directions,
        intensities=truth.evaluate_pole_density(pole, directions, antipodal=True) / scale,
        specimen_frame=specimen,
        antipodal=True,
        includes_symmetry_family=True,
        sampling="sampled_density",
    )


pole_figures = [measured(indices) for indices in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0])]
report = HarmonicODF.invert_pole_figures(
    pole_figures,
    degree_bandlimit=4,
    regularization=1e-6,
    pole_kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=10.0),
    phi1_step_deg=15.0,
    big_phi_step_deg=15.0,
    phi2_step_deg=15.0,
    ghost_correction=True,
)
correction = report.ghost_correction
```

:::

**Compute**

```python
result = float(correction.pole_figure_max_change)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `ghost-correction-leaves-the-measured-fit-untouched` | 0.000878 | 0.000000 | m.r.d. | 8.78e-04 | 5e-03 | ✅ pass |

**Why this value**: Analytic: a Friedel-symmetric pole figure is the integral of the ODF over a kernel even under h -> -h, and an odd-degree generalized spherical harmonic integrates to zero against an even kernel. The change is therefore exactly zero in the continuum; the tolerance is the quadrature discretization error at a 15 degree Bunge step.

**Citation**: Matthies, On the reproducibility of the orientation distribution function of texture samples from pole figures (ghost phenomena), Phys. Status Solidi B 92 (1979) K135-K138.

**See also**: {doc}`The ghost problem <../../theory/ghost_problem_and_odd_harmonics>`, {doc}`Harmonic ODF reconstruction <../../theory/harmonic_odf_reconstruction>`
