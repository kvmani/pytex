<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Dynamical CBED and symmetry determination

The exact limits that calibrate a many-beam calculation, the HOLZ degeneracy that makes voltage calibration mandatory, and the diffraction-group construction that determines a point group including its centre of symmetry.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Two beams reduce the Bloch-wave solver to the closed form exactly

A many-beam dynamical calculation has no independent standard to be checked against, so the one calibration available is its own limiting case. Restricted to a single reflection, the coupled system has the closed-form solution I_g = sin^2(pi t s_eff) / (xi_g s_eff)^2 with s_eff^2 = s^2 + xi_g^-2. Reproducing it to machine precision pins three conventions simultaneously: the diagonal 2 s_g, the off-diagonal scale |nu_g| = 1 / xi_g, and the factor i pi in the propagator. Any one of them wrong yields a rocking curve of the right general shape and the wrong fringe spacing, which is exactly the error that survives a plausibility check.

**Symbols**

- $\xi_{g}$ &mdash; Two-beam extinction distance of reflection g; the depth period of the intensity exchange between the transmitted and diffracted beams.
- $s_{g}$ &mdash; Excitation error: deviation of reflection g from the exact Bragg condition.
- $\nu_{g}$ &mdash; Complex Fourier coefficient of the scaled lattice potential; the off-diagonal element of the dynamical structure matrix, with |nu_g| = 1 / xi_g.
- $t$ &mdash; Foil thickness along the beam.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    AbsorptionModel,
    AtomicSite,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SpaceGroupSpec,
    SymmetrySpec,
    UnitCell,
    ZoneAxis,
    beam_set_for_zone,
    beam_set_from_indices,
    extinction_distance_angstrom,
    holz_line_pattern,
    solve_bloch_waves,
    two_beam_rocking_curve,
)
from pytex.diffraction.kinematic import electron_wavelength_angstrom

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
nickel_lattice = Lattice(3.5239, 3.5239, 3.5239, 90.0, 90.0, 90.0, crystal_frame=crystal)
nickel = Phase(
    name="nickel-fcc",
    lattice=nickel_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=nickel_lattice,
        sites=tuple(
            AtomicSite(
                label=f"Ni{index}",
                species="Ni",
                fractional_coordinates=np.asarray(position, dtype=float),
            )
            for index, position in enumerate(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
            )
        ),
    ),
    space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=crystal),
)
zone = ZoneAxis(np.array([0, 0, 1]), phase=nickel)
wavelength = electron_wavelength_angstrom(200.0)
```

:::

**Compute**

```python
beams = beam_set_from_indices(nickel, zone, [[2, 2, 0]])
g_zone = beams.g_zone[1]
in_plane = float(np.linalg.norm(g_zone[:2]))
zero_tilt = g_zone[2] - 0.5 * wavelength * beams.g_magnitude_inv_angstrom[1] ** 2
targets = np.array([-0.01, -0.003, 0.0, 0.003, 0.01])
tilts = ((zero_tilt - targets) / in_plane)[:, None] * (g_zone[:2] / in_plane)[None, :]
solution = solve_bloch_waves(beams, tilts, thickness_angstrom=800.0)
closed_form = two_beam_rocking_curve(
    targets,
    thickness_angstrom=800.0,
    extinction_distance_angstrom=float(
        extinction_distance_angstrom(nickel, [[2, 2, 0]])[0]
    ),
)
result = float(np.max(np.abs(solution.intensity_of([2, 2, 0]) - closed_form)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-dynamical-two-beam-limit-of-the-many-beam-solver` | 1.14e-15 | 0.00e+00 | &mdash; | 1.14e-15 | 1e-12 | ✅ pass |

**Why this value**: An analytic identity, not a measurement: the two-beam structure matrix is s I + B with B traceless, and the exponential of a traceless 2x2 matrix is cos(pi s_eff t) I + i sin(pi s_eff t) B / s_eff, which gives the Howie-Whelan expression exactly. The deviation must therefore be zero to floating-point rounding.

**Citation**: Howie and Whelan, Proceedings of the Royal Society A 263 (1961) 217-237; Williams and Carter, Transmission Electron Microscopy, 2nd ed. (Springer, 2009), Chapter 23.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`, {doc}`Diffraction API <../../api/index>`

## Without absorption the coupled beams sum to exactly one

The elastic structure matrix is Hermitian, so the propagator exp(i pi A t) is unitary and the beam intensities sum to one at every thickness and every incident direction. This is the only exact global check available on a many-beam calculation, and it is the one that catches the classic implementation error: obtaining the Bloch-wave excitation amplitudes by projection rather than by solving. The eigenvectors of a complex matrix are not orthogonal, so a projection gives rocking curves of the right shape with the wrong contrast - wrong in a way that only this identity reveals.

**Symbols**

- $\nu_{g}$ &mdash; Complex Fourier coefficient of the scaled lattice potential; the off-diagonal element of the dynamical structure matrix, with |nu_g| = 1 / xi_g.
- $s_{g}$ &mdash; Excitation error: deviation of reflection g from the exact Bragg condition.
- $t$ &mdash; Foil thickness along the beam.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    AbsorptionModel,
    AtomicSite,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SpaceGroupSpec,
    SymmetrySpec,
    UnitCell,
    ZoneAxis,
    beam_set_for_zone,
    beam_set_from_indices,
    extinction_distance_angstrom,
    holz_line_pattern,
    solve_bloch_waves,
    two_beam_rocking_curve,
)
from pytex.diffraction.kinematic import electron_wavelength_angstrom

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
nickel_lattice = Lattice(3.5239, 3.5239, 3.5239, 90.0, 90.0, 90.0, crystal_frame=crystal)
nickel = Phase(
    name="nickel-fcc",
    lattice=nickel_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=nickel_lattice,
        sites=tuple(
            AtomicSite(
                label=f"Ni{index}",
                species="Ni",
                fractional_coordinates=np.asarray(position, dtype=float),
            )
            for index, position in enumerate(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
            )
        ),
    ),
    space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=crystal),
)
zone = ZoneAxis(np.array([0, 0, 1]), phase=nickel)
wavelength = electron_wavelength_angstrom(200.0)
```

:::

**Compute**

```python
beams = beam_set_for_zone(nickel, zone, convergence_semi_angle_mrad=6.0)
tilts = np.array([[0.0, 0.0], [3e-3, -2e-3], [-4e-3, 1e-3], [2e-3, 5e-3]])
deviations = [
    float(np.max(np.abs(
        solve_bloch_waves(beams, tilts, thickness_angstrom=t).total_intensity - 1.0
    )))
    for t in (100.0, 700.0, 2500.0)
]
result = float(max(deviations))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-dynamical-intensity-is-conserved-without-absorption` | 2.11e-15 | 0.00e+00 | &mdash; | 2.11e-15 | 1e-12 | ✅ pass |

**Why this value**: Unitarity of the propagator of a Hermitian generator: sum_g |psi_g|^2 is conserved exactly. The expected value is zero by theorem, with the tolerance set by floating-point accumulation over the beam set rather than by any physical uncertainty.

**Citation**: Hirsch, Howie, Nicholson, Pashley and Whelan, Electron Microscopy of Thin Crystals, 2nd ed. (Krieger, 1977), Chapter 10.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`, {doc}`Diffraction API <../../api/index>`

## A HOLZ line cannot separate lattice strain from accelerating voltage

HOLZ line positions are the sharpest lattice-parameter measurement a convergent-beam pattern offers, and this is the reason the measurement begins with a calibration rather than a specimen. Scaling the lattice by 1 + eps shrinks every g by the same factor, and the line offset d_g = (g_z - lambda |g|^2 / 2) / |g_perp| then depends on eps and on lambda through the same term with opposite signs. A fractional change in lattice parameter and a fractional change in wavelength therefore cancel exactly, at every reflection simultaneously - so a lattice parameter quoted from an uncalibrated microscope is a measurement of its high-tension supply.

**Symbols**

- $\lambda$ &mdash; Radiation wavelength.
- $s_{g}$ &mdash; Excitation error: deviation of reflection g from the exact Bragg condition.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    AbsorptionModel,
    AtomicSite,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SpaceGroupSpec,
    SymmetrySpec,
    UnitCell,
    ZoneAxis,
    beam_set_for_zone,
    beam_set_from_indices,
    extinction_distance_angstrom,
    holz_line_pattern,
    solve_bloch_waves,
    two_beam_rocking_curve,
)
from pytex.diffraction.kinematic import electron_wavelength_angstrom

crystal = ReferenceFrame(
    name="crystal", domain=FrameDomain.CRYSTAL, axes=("a", "b", "c"), handedness=Handedness.RIGHT
)
nickel_lattice = Lattice(3.5239, 3.5239, 3.5239, 90.0, 90.0, 90.0, crystal_frame=crystal)
nickel = Phase(
    name="nickel-fcc",
    lattice=nickel_lattice,
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
    unit_cell=UnitCell(
        lattice=nickel_lattice,
        sites=tuple(
            AtomicSite(
                label=f"Ni{index}",
                species="Ni",
                fractional_coordinates=np.asarray(position, dtype=float),
            )
            for index, position in enumerate(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
            )
        ),
    ),
    space_group=SpaceGroupSpec(symbol="Fm-3m", number=225, reference_frame=crystal),
)
zone = ZoneAxis(np.array([0, 0, 1]), phase=nickel)
wavelength = electron_wavelength_angstrom(200.0)
```

:::

**Compute**

```python
lines = holz_line_pattern(
    nickel,
    zone,
    convergence_semi_angle_mrad=8.0,
    max_index=24,
    g_max_inv_angstrom=6.0,
)
strain = 1e-3
result = float(max(
    abs(
        line.offset_at(
            lattice_strain=strain,
            wavelength_angstrom=wavelength * (1.0 + strain),
        )
        - line.offset_rad
    )
    for line in lines.lines
))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-holz-strain-and-wavelength-are-exactly-degenerate` | 1.56e-17 | 0.00e+00 | radian | 1.56e-17 | 1e-15 | ✅ pass |

**Why this value**: An exact cancellation in the closed form: substituting lambda -> lambda (1 + eps) into d_g(eps, lambda) recovers d_g(0, lambda) identically, for every reflection. The expected value is zero by algebra, and the tolerance is floating-point rounding.

**Citation**: Jones, Rackham and Steeds, Proceedings of the Royal Society A 354 (1977) 197-222, for HOLZ line lattice-parameter determination; Williams and Carter, Transmission Electron Microscopy, 2nd ed. (2009), Chapter 21.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`, {doc}`Diffraction API <../../api/index>`

## The diffraction-group construction yields Buxton's 31 groups

The 31 diffraction groups are usually quoted as a table. PyTex derives them instead: each crystal-point-group operator is classified by its action on the beam direction and contributes its transverse restriction, tagged with the reciprocity flag when it reverses the beam. That map lands in a subgroup of (plane point group) x Z2, and scanning all 32 crystallographic point groups over their characteristic beam directions must realize exactly the 31 subgroups Buxton, Eades, Steeds and Rackham enumerated. Reaching 30 or 32 would mean the construction or the stored operators are wrong - a check no transcribed table can perform on itself.

:::{dropdown} Setup (imports and object construction)

```python
from pytex import (
    SymmetryObservations,
    determine_point_group,
    diffraction_group_for,
    diffraction_group_symbols,
)
```

:::

**Compute**

```python
result = len(diffraction_group_symbols())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-groups-construction-yields-buxtons-thirty-one` | 31 | 31 | &mdash; | exact | exact | ✅ pass |

**Why this value**: The published count of diffraction groups: 10 with no reciprocity element, 10 direct products with Z2 (suffix 1_R), and 11 graphs of a surjection onto Z2. An exact integer, so the tolerance is zero.

**Citation**: Buxton, Eades, Steeds and Rackham, Philosophical Transactions of the Royal Society A 281 (1976) 171-194.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`, {doc}`Diffraction API <../../api/index>`

## The plus-minus-g observation splits the 32 point groups into 21 and 11

This is the arithmetic of the whole technique. Friedel's law makes kinematic diffraction blind to a centre of symmetry, so a selected-area pattern determines only the Laue class: 11 possibilities where there are 32 point groups. The diffraction-group element 2_R requires an operator acting as -1 on the beam direction and as -1 on the transverse plane, which is the inversion and nothing else - so observing whether the +g and -g discs are related by a two-fold recovers exactly the distinction Friedel's law destroyed, partitioning the 32 point groups into the 21 acentric and the 11 centric ones.

:::{dropdown} Setup (imports and object construction)

```python
from pytex import (
    SymmetryObservations,
    determine_point_group,
    diffraction_group_for,
    diffraction_group_symbols,
)
```

:::

**Compute**

```python
acentric = determine_point_group(SymmetryObservations(friedel_pair_two_fold=False))
centric = determine_point_group(SymmetryObservations(friedel_pair_two_fold=True))
result = [len(acentric.point_groups), len(centric.point_groups)]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-groups-friedel-observation-splits-the-point-groups` | [21, 11] | [21, 11] | &mdash; | exact | exact | ✅ pass |

**Why this value**: Of the 32 crystallographic point groups, 11 contain the inversion - the Laue classes - and 21 do not. Exact integers from the International Tables, so the tolerance is zero.

**Citation**: International Tables for Crystallography, Volume A, Chapter 10, for the 32 point groups and the 11 Laue classes; Buxton, Eades, Steeds and Rackham, Philosophical Transactions of the Royal Society A 281 (1976) 171-194, for the 2_R correspondence.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`, {doc}`Diffraction API <../../api/index>`

## Zincblende down [001] gives 4_Rmm_R, and no centre of symmetry

The textbook case. Down a four-fold zone the centrosymmetric cubic group m-3m gives diffraction group 4mm1_R, four-fold in both the bright-field disc and the whole pattern; zincblende -43m gives 4_Rmm_R, four-fold in the disc but only two-fold in the pattern. That difference in whole-pattern symmetry is what a CBED exposure of gallium arsenide reads, and it is the observation that determines the absence of a centre of symmetry - which no kinematic pattern can do. The check below counts the whole-pattern operations: 4 for -43m against 8 for m-3m, that is 2mm against 4mm.

:::{dropdown} Setup (imports and object construction)

```python
from pytex import (
    SymmetryObservations,
    determine_point_group,
    diffraction_group_for,
    diffraction_group_symbols,
)
```

:::

**Compute**

```python
polar = diffraction_group_for('-43m', [0, 0, 1])
centric = diffraction_group_for('m-3m', [0, 0, 1])
result = [
    float(polar.symbol == '4_Rmm_R'),
    float(polar.bright_field_symbol == '4mm'),
    float(polar.whole_pattern_symbol == '2mm'),
    float(polar.has_friedel_symmetry),
    float(centric.symbol == '4mm1_R'),
    float(centric.whole_pattern_symbol == '4mm'),
    float(centric.has_friedel_symmetry),
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `diffraction-groups-zincblende-down-001-loses-the-two-fold` | [1, 1, 1, 0, 1, 1, 1] | [1, 1, 1, 0, 1, 1, 1] | &mdash; | exact | exact | ✅ pass |

**Why this value**: The published diffraction-group assignments for the cubic acentric and centric groups viewed down a four-fold axis: -43m gives 4_Rmm_R with bright-field 4mm over whole-pattern 2mm and no 2_R, while m-3m gives 4mm1_R with 4mm in both and 2_R present. Each entry is a boolean agreement, so the tolerance is zero.

**Citation**: Buxton, Eades, Steeds and Rackham, Philosophical Transactions of the Royal Society A 281 (1976) 171-194, Tables 2 and 3; Williams and Carter, Transmission Electron Microscopy, 2nd ed. (2009), Chapter 21.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`, {doc}`Diffraction API <../../api/index>`
