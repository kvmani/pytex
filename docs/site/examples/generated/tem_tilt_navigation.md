<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# TEM tilt navigation

Holder tilts that bring a target zone axis onto the electron beam: analytic interzonal travel for the standard cubic transitions, the closed-form solid angle a double-tilt holder commands, the cost of an uncalibrated diffraction rotation, and the group-order counts that decide whether a single indexed pattern leaves a real ambiguity.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Crystal travel from [001] to [011] in a cubic crystal

You are down the [001] zone of an FCC metal and want [011]. The engine solves the holder angles and plans the path; the crystal travel along that path must equal the interplanar angle between the two zone axes, which for a cubic crystal is an analytic quantity independent of lattice parameter. This is the basic sanity check on the whole chain: orientation, closed-form solution, forward validation, and geodesic path planning.

**Symbols**

- $[uvw]$ &mdash; Lattice direction brought parallel to the beam.
- $\alpha$ &mdash; Holder tilt about the rod axis.
- $\beta$ &mdash; Holder tilt about the cradle axis carried in the rod.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CurrentState,
    DoubleTiltStage,
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    RectangularEnvelope,
    ReferenceFrame,
    StagePosition,
    SymmetrySpec,
    ZoneAxis,
    plan_tilt_to_zone_axis,
)
from pytex.tem.reconstruction import HOLDER_FRAME

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
# A generous holder, so the geometry rather than the envelope decides.
wide_stage = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))
# The crystal sits with [001] along the beam at zero tilt.
aligned = CurrentState.from_orientation(
    Orientation.from_matrix(
        np.eye(3), specimen_frame=HOLDER_FRAME, phase=nickel, crystal_frame=crystal
    ),
    StagePosition(0.0, 0.0),
    current_zone_axis=ZoneAxis([0, 0, 1], phase=nickel),
)
```

:::

**Compute**

```python
report = plan_tilt_to_zone_axis(
    aligned, ZoneAxis([0, 1, 1], phase=nickel), wide_stage
)
result = float(report.best().path.total_travel_deg)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `tem-tilt-001-to-011-travel` | 45.000 | 45.000 | deg | < 1e-12 | 1e-03 | ✅ pass |

**Why this value**: The angle between [001] and [011] in a cubic lattice is arccos(1/sqrt(2)) = 45 degrees exactly, from the dot product of the two directions divided by their lengths. Independent of the lattice parameter.

**Citation**: Edington, J. W., Practical Electron Microscopy in Materials Science, Macmillan; standard cubic interzonal-angle tables.

**See also**: {doc}`TEM tilt navigation foundation <../../architecture/tem_tilt_navigation_foundation>`, {doc}`TEM tilt navigation notebook <../../tutorials/notebooks/24_tem_tilt_navigation>`

## Crystal travel from [001] to [111] in a cubic crystal

The [001] to [111] hop is the one every TEM course teaches, and it is the longest of the common cubic transitions — far enough that a typical +/-30 degree double-tilt holder cannot make it without help from a symmetry equivalent. Here a wide holder is used so that the geometry alone is tested.

**Symbols**

- $[uvw]$ &mdash; Lattice direction brought parallel to the beam.
- $\theta$ &mdash; Angle between the current and target zone axes.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CurrentState,
    DoubleTiltStage,
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    RectangularEnvelope,
    ReferenceFrame,
    StagePosition,
    SymmetrySpec,
    ZoneAxis,
    plan_tilt_to_zone_axis,
)
from pytex.tem.reconstruction import HOLDER_FRAME

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
# A generous holder, so the geometry rather than the envelope decides.
wide_stage = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))
# The crystal sits with [001] along the beam at zero tilt.
aligned = CurrentState.from_orientation(
    Orientation.from_matrix(
        np.eye(3), specimen_frame=HOLDER_FRAME, phase=nickel, crystal_frame=crystal
    ),
    StagePosition(0.0, 0.0),
    current_zone_axis=ZoneAxis([0, 0, 1], phase=nickel),
)
```

:::

**Compute**

```python
report = plan_tilt_to_zone_axis(
    aligned, ZoneAxis([1, 1, 1], phase=nickel), wide_stage
)
result = float(report.best().path.total_travel_deg)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `tem-tilt-001-to-111-travel` | 54.7356 | 54.7356 | deg | < 1e-12 | 1e-03 | ✅ pass |

**Why this value**: arccos(1/sqrt(3)) = 54.7356 degrees, the analytic angle between <001> and <111> in a cubic lattice, and the standard tetrahedral-angle complement.

**Citation**: Williams, D. B. and Carter, C. B., Transmission Electron Microscopy, Springer, DOI: 10.1007/978-0-387-76501-3.

**See also**: {doc}`TEM tilt navigation foundation <../../architecture/tem_tilt_navigation_foundation>`, {doc}`TEM tilt navigation notebook <../../tutorials/notebooks/24_tem_tilt_navigation>`

## Solid angle a +/-30 degree double-tilt holder reaches

Before asking whether a particular zone axis is reachable, it is worth knowing how much of orientation space the holder commands at all. Because the beam direction in holder coordinates is a spherical coordinate system whose pole is the beta axis, the Jacobian is cos(alpha) and the integral is elementary. The answer — a little over eight percent of all directions — is why symmetry equivalents matter so much in practice.

**Symbols**

- $\Omega$ &mdash; Solid angle of beam directions a holder can reach.
- $\alpha$ &mdash; Holder tilt about the rod axis.
- $\beta$ &mdash; Holder tilt about the cradle axis carried in the rod.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CurrentState,
    DoubleTiltStage,
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    RectangularEnvelope,
    ReferenceFrame,
    StagePosition,
    SymmetrySpec,
    ZoneAxis,
    plan_tilt_to_zone_axis,
)
from pytex.tem.reconstruction import HOLDER_FRAME

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
# A generous holder, so the geometry rather than the envelope decides.
wide_stage = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))
# The crystal sits with [001] along the beam at zero tilt.
aligned = CurrentState.from_orientation(
    Orientation.from_matrix(
        np.eye(3), specimen_frame=HOLDER_FRAME, phase=nickel, crystal_frame=crystal
    ),
    StagePosition(0.0, 0.0),
    current_zone_axis=ZoneAxis([0, 0, 1], phase=nickel),
)
```

:::

**Compute**

```python
envelope = RectangularEnvelope(-30.0, 30.0, -30.0, 30.0)
result = float(envelope.accessible_solid_angle_sr())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `tem-holder-accessible-solid-angle` | 1.04720 | 1.04720 | sr | < 1e-12 | 1e-09 | ✅ pass |

**Why this value**: Omega = (beta_max - beta_min) * (sin alpha_max - sin alpha_min) = (pi/3) * (2 sin 30 deg) = pi/3 = 1.04720 sr, integrating the cos(alpha) Jacobian of the beam-direction map over the tilt rectangle.

**Citation**: Derived in section 10.2 of docs/architecture/tem_tilt_navigation_foundation.md; standard spherical-measure result.

**See also**: {doc}`TEM tilt navigation foundation <../../architecture/tem_tilt_navigation_foundation>`

## Angular miss from a 5 degree diffraction-rotation error over a 90 degree hop

The diffraction rotation is not recorded by instrument metadata and must be calibrated. This example quantifies what an uncalibrated value costs: the miss is 2 asin(sin(dphi/2) sin(theta)), which grows with the length of the hop. The same 5 degree error costs 0.44 degrees over a 5 degree hop and the full 5 degrees over a 90 degree one — which is the argument for routing a long excursion through intermediate zones and re-indexing at each.

**Symbols**

- $\theta$ &mdash; Angle between the current and target zone axes.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.tem.calibration import residual_from_rotation_error_deg
```

:::

**Compute**

```python
result = residual_from_rotation_error_deg(5.0, 90.0)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `tem-diffraction-rotation-residual` | 5.000000 | 5.000000 | deg | < 1e-12 | 1e-09 | ✅ pass |

**Why this value**: For theta = 90 degrees, sin(theta) = 1, so the expression reduces to 2 asin(sin(dphi/2)) = dphi exactly. The residual equals the calibration error itself when the target is perpendicular to the current zone axis.

**Citation**: Derived in section 8.2 of docs/architecture/tem_tilt_navigation_foundation.md and verified numerically over 3000 random orientations.

**See also**: {doc}`TEM tilt navigation foundation <../../architecture/tem_tilt_navigation_foundation>`

## Order of the observation stabilizer for cubic m-3m down [001]

A single indexed SAED pattern determines the orientation only up to the rotations of the Laue class that map the zone plane to itself. Counting them answers the question that decides whether the classical 180-degree ambiguity matters: for cubic m-3m down [001] the stabilizer is the group 422, of order 8, and every one of its operators is already a crystal symmetry — so the ambiguity is entirely absorbed and nothing is left undetermined.

**Symbols**

- $[uvw]$ &mdash; Lattice direction brought parallel to the beam.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.tem.ambiguity import observation_stabilizer
```

:::

**Compute**

```python
result = len(observation_stabilizer("m-3m", [0.0, 0.0, 1.0]))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `tem-observation-stabilizer-cubic-001` | 8 | 8 | operators | exact | exact | ✅ pass |

**Why this value**: The rotations of m-3m fixing the [001] axis line form the point group 422: the identity, three rotations about [001] (90, 180, 270 degrees), and four two-fold rotations about the in-plane <100> and <110> axes. Order 8, from International Tables Volume A.

**Citation**: Hahn, Th. (ed.), International Tables for Crystallography, Volume A, IUCr/Springer, DOI: 10.1107/97809553602060000100.

**See also**: {doc}`TEM tilt navigation foundation <../../architecture/tem_tilt_navigation_foundation>`

## Number of symmetry-equivalent targets for a general cubic direction

The user asks for one zone axis; the crystal offers a whole orbit, every member of which gives an identical diffraction pattern. The choice among them is therefore free, and the engine takes the cheapest reachable one — which is frequently the difference between a target being reachable and not. For a general direction in a cubic crystal the orbit has one member per proper operator, in both senses.

**Symbols**

- $[uvw]$ &mdash; Lattice direction brought parallel to the beam.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CurrentState,
    DoubleTiltStage,
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    Phase,
    RectangularEnvelope,
    ReferenceFrame,
    StagePosition,
    SymmetrySpec,
    ZoneAxis,
    plan_tilt_to_zone_axis,
)
from pytex.tem.reconstruction import HOLDER_FRAME

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
# Nickel (FCC), lattice parameter from the pinned PyTex fixture corpus.
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
# A generous holder, so the geometry rather than the envelope decides.
wide_stage = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))
# The crystal sits with [001] along the beam at zero tilt.
aligned = CurrentState.from_orientation(
    Orientation.from_matrix(
        np.eye(3), specimen_frame=HOLDER_FRAME, phase=nickel, crystal_frame=crystal
    ),
    StagePosition(0.0, 0.0),
    current_zone_axis=ZoneAxis([0, 0, 1], phase=nickel),
)
```

:::

**Compute**

```python
report = plan_tilt_to_zone_axis(
    aligned, ZoneAxis([1, 3, 5], phase=nickel), wide_stage, include_paths=False
)
result = int(report.orbit_size)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `tem-symmetry-orbit-multiplicity` | 48 | 48 | directions | exact | exact | ✅ pass |

**Why this value**: The proper point group of m-3m is 432, of order 24. A general direction has a trivial stabilizer, so its orbit has 24 members; counting both senses of each gives 48 distinct directions.

**Citation**: Hahn, Th. (ed.), International Tables for Crystallography, Volume A, IUCr/Springer, DOI: 10.1107/97809553602060000100.

**See also**: {doc}`TEM tilt navigation foundation <../../architecture/tem_tilt_navigation_foundation>`, {doc}`TEM tilt navigation notebook <../../tutorials/notebooks/24_tem_tilt_navigation>`

## Crystal orientation from an indexed pattern: the identity case

The case that fixes every sign in the chain from indexing to orientation. A crystal whose [001] is on the beam at zero tilt, recorded with a zero diffraction rotation and an identity crystal-to-pattern rotation, must come out at the identity orientation. Anything transposed or composed in the wrong order moves the answer off it, so this is the cheapest possible guard on the composition U = R_stage^T F Rz(phi_D) R, which is what turns a solved SAED pattern plus the holder tilts into a reportable orientation.

**Symbols**

- $\varphi_D$ &mdash; Diffraction rotation of the recorded pattern.
- $(\varphi_1, \Phi, \varphi_2)$ &mdash; Bunge Euler angles of a crystal orientation.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    DoubleTiltStage,
    FrameDomain,
    Handedness,
    IndexedPatternObservation,
    Lattice,
    Phase,
    RectangularEnvelope,
    ReferenceFrame,
    StageCalibration,
    StagePosition,
    SymmetrySpec,
    ZoneAxis,
    orientation_from_indexed_pattern,
    orientation_from_indexed_patterns,
    solve_tilts_for_direction,
)
from pytex.tem.stage import rotation_z

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
wide_stage = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))


def pattern_rotation_for(orientation, position, diffraction_rotation_deg):
    """The crystal-to-pattern rotation indexing would report, obtained by
    inverting U = R_stage^T Rz(phi_D) R so that R = Rz(-phi_D) R_stage U."""
    stage_matrix = wide_stage.rotation_matrix(position.alpha_deg, position.beta_deg)
    return rotation_z(np.deg2rad(-diffraction_rotation_deg)) @ stage_matrix @ orientation


def position_for(orientation, zone):
    """Stage angles putting a zone axis on the beam under this orientation."""
    return StagePosition(*solve_tilts_for_direction(orientation @ zone.unit_vector)[0])
```

:::

**Compute**

```python
stage = DoubleTiltStage(
    calibration=StageCalibration(diffraction_rotation_deg=0.0)
)
indexed = orientation_from_indexed_pattern(
    np.eye(3),
    ZoneAxis([0, 0, 1], phase=nickel),
    StagePosition(0.0, 0.0),
    stage,
)
result = float(np.max(np.abs(indexed.matrix - np.eye(3))))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `tem-indexed-orientation-identity` | 0.00e+00 | 0.00e+00 | dimensionless | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: At zero tilt the stage rotation is the identity; at zero diffraction rotation so is Rz; and an unmirrored pattern has identity parity. The composition therefore reduces to the crystal-to-pattern rotation itself, which is the identity by construction, so the deviation is exactly zero.

**Citation**: Composition derived in section 5 of docs/architecture/tem_tilt_navigation_foundation.md.

**See also**: {doc}`TEM tilt navigation foundation <../../architecture/tem_tilt_navigation_foundation>`, {doc}`TEM tilt navigation notebook <../../tutorials/notebooks/24_tem_tilt_navigation>`

## Diffraction rotation recovered from two indexed patterns

The diffraction rotation is the one constant this subsystem needs and the instrument does not report. This example shows it need not be supplied at all: two patterns indexed at two stage positions determine the crystal orientation *and* that constant together. A rotation of 37 degrees is planted in synthetic patterns and recovered from them, with no calibration given to the stage — which is what makes subsequent single-pattern orientations trustworthy rather than inherited.

**Symbols**

- $\varphi_D$ &mdash; Diffraction rotation of the recorded pattern.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    DoubleTiltStage,
    FrameDomain,
    Handedness,
    IndexedPatternObservation,
    Lattice,
    Phase,
    RectangularEnvelope,
    ReferenceFrame,
    StageCalibration,
    StagePosition,
    SymmetrySpec,
    ZoneAxis,
    orientation_from_indexed_pattern,
    orientation_from_indexed_patterns,
    solve_tilts_for_direction,
)
from pytex.tem.stage import rotation_z

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
nickel = Phase(
    "nickel-fcc",
    lattice=Lattice(3.52387, 3.52387, 3.52387, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
wide_stage = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))


def pattern_rotation_for(orientation, position, diffraction_rotation_deg):
    """The crystal-to-pattern rotation indexing would report, obtained by
    inverting U = R_stage^T Rz(phi_D) R so that R = Rz(-phi_D) R_stage U."""
    stage_matrix = wide_stage.rotation_matrix(position.alpha_deg, position.beta_deg)
    return rotation_z(np.deg2rad(-diffraction_rotation_deg)) @ stage_matrix @ orientation


def position_for(orientation, zone):
    """Stage angles putting a zone axis on the beam under this orientation."""
    return StagePosition(*solve_tilts_for_direction(orientation @ zone.unit_vector)[0])
```

:::

**Compute**

```python
truth = np.linalg.qr(np.random.default_rng(23).normal(size=(3, 3)))[0]
if np.linalg.det(truth) < 0:
    truth[:, 0] *= -1

observations = []
for indices in ([0, 0, 1], [0, 1, 1]):
    zone = ZoneAxis(indices, phase=nickel)
    position = position_for(truth, zone)
    observations.append(
        IndexedPatternObservation(
            pattern_rotation_for(truth, position, 37.0), zone, position
        )
    )

# wide_stage carries no diffraction-rotation calibration at all.
fit = orientation_from_indexed_patterns(observations, wide_stage)
result = float(fit.diffraction_rotation_deg)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `tem-self-calibrated-diffraction-rotation` | 37.000000 | 37.000000 | deg | < 1e-12 | 1e-09 | ✅ pass |

**Why this value**: The value planted in the synthetic patterns. Recovery is exact rather than fitted: once the zone axes fix the orientation, the residual R_stage U R^T is by construction a pure rotation about the beam axis whose angle is the diffraction rotation, so the value is read off directly.

**Citation**: Britton, T. B. et al., Materials Characterization 117 (2016) 113-126, DOI: 10.1016/j.matchar.2016.04.008, on why this constant must be measured rather than assumed.

**See also**: {doc}`TEM tilt navigation foundation <../../architecture/tem_tilt_navigation_foundation>`, {doc}`TEM tilt navigation notebook <../../tutorials/notebooks/24_tem_tilt_navigation>`
