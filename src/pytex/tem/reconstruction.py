"""Determining the crystal-to-holder orientation from what the microscope shows.

Everything in `pytex.tem` reduces to inverting one rigid-body equation,

    v_lab = R_stage(alpha, beta) U v_crystal,                              (M)

for the single unknown rotation ``U``. This module builds ``U`` from the data an
operator actually has, and — as importantly — records *which* data were used,
because what can and cannot be resolved depends entirely on that.

The specimen/holder collapse
----------------------------
A textbook treatment would insert a specimen frame between crystal and holder,
joined by a mounting rotation. **Those two rotations are not separately
identifiable from diffraction data, ever**: no observable in this problem depends
on either alone, only on their product. Carrying them separately manufactures an
unobservable degree of freedom, which is the classic route to a sign error that
hides for months. This module therefore carries the single crystal-to-holder
rotation ``U``, and — because the holder frame is declared to *be* the
specimen-domain frame for TEM work — ``U`` is a plain `pytex.core.orientation.Orientation`
with no new concept required.

The three modes, and why the obvious one is not the recommended one
------------------------------------------------------------------
**Mode A** (one indexed pattern plus a calibrated diffraction rotation) is one
line of algebra and is what most people reach for. It is also the mode that
carries the whole instrumental ambiguity: if the diffraction rotation is wrong by
180 degrees, both tilt angles come out negated, the calculation reports a clean
zero residual, and the operator drives the specimen exactly the wrong way.

**Mode B** (two indexed zone axes at two stage positions) needs **no** diffraction
rotation and no parity bit at all: it uses only zone-axis identities and stage
readouts, and solves the resulting two-vector attitude problem exactly. Since an
operator who chased Kikuchi bands to get here has almost always visited a second
zone, this is the recommended path. It also supplies a calibration-free
consistency test — the interzonal angle is a crystallographic invariant, so
comparing it with the angle the stage model predicts indicts a wrong sign
convention, a wrong indexing, or a bent specimen at zero experimental cost.

**Mode C** is Mode A after running the two-excursion calibration in
`pytex.tem.calibration`, which measures the diffraction rotation and the parity
directly.

See ``docs/architecture/tem_tilt_navigation_foundation.md`` section 5.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._angles import angle_between_unit_vectors_rad
from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.frame_catalog import specimen_frame as catalog_specimen_frame
from pytex.core.lattice import Phase, ZoneAxis
from pytex.core.orientation import Orientation, Rotation
from pytex.core.provenance import ProvenanceRecord
from pytex.tem.ambiguity import AmbiguityReport, analyze_ambiguity
from pytex.tem.stage import (
    BEAM_AXIS_LABORATORY,
    StageModel,
    StagePosition,
    rotation_z,
)

__all__ = [
    "HOLDER_FRAME",
    "CurrentState",
    "ReconstructionMode",
    "ZoneAxisObservation",
]

#: Interzonal-angle residual above which a two-zone reconstruction is indicted.
#:
#: Two degrees is comfortably above the error a well-indexed zone axis carries
#: (a few tenths of a degree) and comfortably below the discrepancy a reversed
#: sign convention or a mis-indexed zone produces (tens of degrees), so the test
#: separates the two cleanly rather than firing on ordinary noise.
_CONSISTENCY_WARNING_DEG = 2.0

#: The holder frame, declared as the specimen-domain frame for TEM work.
#:
#: Per `docs/standards/notation_and_conventions.md` no subsystem may invent a new
#: frame domain. The holder is rigidly attached to the specimen, so it *is* the
#: specimen frame here; naming it explicitly keeps the declaration visible
#: instead of implied.
HOLDER_FRAME = catalog_specimen_frame(
    name="holder",
    description=(
        "Rigidly attached to the specimen and to the holder cartridge: x along the "
        "holder rod (the alpha axis), y along the nominal beta axis, z completing the "
        "right-handed set. Coincides with the microscope frame at zero tilt."
    ),
)


class ReconstructionMode(StrEnum):
    """How the crystal-to-holder orientation was obtained.

    The mode is recorded rather than re-inferred downstream, because the
    ambiguity content differs radically between them and a consumer must not
    have to guess.
    """

    #: One indexed pattern plus a calibrated diffraction rotation.
    SINGLE_PATTERN = "single_pattern"
    #: Two indexed zone axes at two stage positions; needs no rotation calibration.
    TWO_ZONE_AXES = "two_zone_axes"
    #: Least squares over three or more indexed zones.
    MULTI_ZONE_FIT = "multi_zone_fit"
    #: Supplied directly; synthetic, teaching, or externally determined.
    KNOWN_ORIENTATION = "known_orientation"


@dataclass(frozen=True, slots=True)
class ZoneAxisObservation:
    """One indexed zone axis, recorded at a known stage position.

    Purpose
    -------
    The atomic datum of the two-zone and multi-zone reconstruction paths: *this
    crystallographic direction was along the beam when the stage read these two
    angles*. That is all the recommended path needs — no pattern rotation, no
    parity, no detector model.

    Attributes
    ----------
    zone_axis : ZoneAxis
        The indexed zone axis. Its sense may be wrong; the reconstruction tests
        both and reports which fits.
    position : StagePosition
    weight : float
        Relative confidence, used by the multi-zone least-squares fit. A zone
        indexed from a clean pattern with many spots deserves more weight than
        one scraped from three faint reflections.
    label : str
        Free-form identifier carried into reports.
    """

    zone_axis: ZoneAxis
    position: StagePosition
    weight: float = 1.0
    label: str = ""

    def __post_init__(self) -> None:
        if not math.isfinite(self.weight) or self.weight <= 0.0:
            raise ValueError("ZoneAxisObservation.weight must be finite and positive.")

    @property
    def unit_vector(self) -> np.ndarray:
        """The zone axis as a Cartesian crystal-frame unit vector."""

        return np.asarray(self.zone_axis.unit_vector, dtype=np.float64)


def _kabsch_rotation(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Proper rotation best carrying ``source`` rows onto ``target`` rows.

    Standard orthogonal Procrustes with the determinant correction that keeps
    the result a rotation rather than a reflection.
    """

    correlation = target.T @ source
    u, _, vt = np.linalg.svd(correlation)
    correction = np.diag([1.0, 1.0, float(np.sign(np.linalg.det(u @ vt)))])
    return np.asarray(u @ correction @ vt, dtype=np.float64)


def _triad(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """Right-handed orthonormal triad from two non-parallel vectors."""

    e1 = normalize_vector(first)
    residual = second - float(np.dot(second, e1)) * e1
    norm = float(np.linalg.norm(residual))
    if norm < 1e-10:
        raise ValueError(
            "The two zone axes are parallel (or antiparallel) to within 1e-10, so they "
            "cannot determine an orientation. Index a second zone that is well "
            "separated from the first — a few tens of degrees is ample."
        )
    e2 = residual / norm
    return np.column_stack([e1, e2, np.cross(e1, e2)])


@dataclass(frozen=True, slots=True)
class CurrentState:
    """The reconstructed crystal-to-holder orientation and how it was obtained.

    Purpose
    -------
    The input to every tilt calculation: where the crystal currently sits
    relative to the holder, at which stage position, with what confidence, and
    with what left undetermined. Constructed through one of the three classmethod
    paths rather than directly, so that the reconstruction mode and its ambiguity
    content are always recorded together with the number.

    Attributes
    ----------
    orientation : Orientation
        The crystal-to-holder rotation ``U``. Its specimen frame is
        :data:`HOLDER_FRAME`.
    position : StagePosition
        Where the stage read when the reconstruction was made.
    mode : ReconstructionMode
    current_zone_axis : ZoneAxis or None
        The zone axis that was along the beam, when one was indexed.
    ambiguity : AmbiguityReport
    consistency_residual_deg : float or None
        For the two- and multi-zone paths: the discrepancy between the
        crystallographic interzonal angle and the angle the stage model
        predicts. A calibration-free indictment of a wrong sign convention,
        a wrong indexing, or a bent specimen. ``None`` when not applicable.
    orientation_uncertainty_deg : float
        One-sigma scale of the orientation estimate.
    notes : tuple of str
    provenance : ProvenanceRecord or None
    """

    orientation: Orientation
    position: StagePosition
    mode: ReconstructionMode
    current_zone_axis: ZoneAxis | None
    ambiguity: AmbiguityReport
    consistency_residual_deg: float | None = None
    orientation_uncertainty_deg: float = 0.5
    notes: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "notes", tuple(self.notes))
        if self.orientation_uncertainty_deg < 0.0:
            raise ValueError("orientation_uncertainty_deg must be non-negative.")

    @property
    def phase(self) -> Phase:
        """The phase this state describes."""

        if self.orientation.phase is None:
            raise ValueError("CurrentState.orientation must carry a phase.")
        return self.orientation.phase

    @property
    def matrix(self) -> np.ndarray:
        """The crystal-to-holder rotation ``U`` as a 3x3 matrix."""

        return np.asarray(self.orientation.as_matrix(), dtype=np.float64)

    def beam_direction_crystal(self, stage: StageModel) -> np.ndarray:
        """The beam direction in crystal coordinates at the current position.

        Purpose
        -------
        What the operator is currently looking down, expressed
        crystallographically. Used to seed path planning, to label the
        stereographic figure, and to check that the reconstruction reproduces
        the zone axis that was actually indexed.
        """

        beam_holder = stage.beam_direction(self.position.alpha_deg, self.position.beta_deg)
        return np.asarray(self.matrix.T @ beam_holder, dtype=np.float64)

    # ------------------------------------------------------------------ #
    # Constructors
    # ------------------------------------------------------------------ #

    @classmethod
    def from_orientation(
        cls,
        orientation: Orientation,
        position: StagePosition,
        *,
        current_zone_axis: ZoneAxis | None = None,
        orientation_uncertainty_deg: float = 0.0,
        provenance: ProvenanceRecord | None = None,
    ) -> CurrentState:
        """Build from a directly known crystal-to-holder orientation.

        Purpose
        -------
        For synthetic studies, teaching, ground-truth tests, and the case where
        an orientation came from an independent measurement. No reconstruction is
        performed, so no reconstruction ambiguity is introduced.

        Parameters
        ----------
        orientation : Orientation
            Crystal-to-holder. Its specimen frame should be :data:`HOLDER_FRAME`.
        position : StagePosition
        current_zone_axis : ZoneAxis, optional
            Recorded for reporting; not used in the algebra.
        orientation_uncertainty_deg : float, default 0
            Zero is honest here: a supplied orientation is exact by construction.

        Returns
        -------
        CurrentState
        """

        if orientation.phase is None:
            raise ValueError(
                "CurrentState.from_orientation requires an Orientation carrying a phase; "
                "the phase supplies the symmetry that decides which targets are "
                "equivalent."
            )
        zone_cartesian = (
            current_zone_axis.unit_vector
            if current_zone_axis is not None
            else np.array([0.0, 0.0, 1.0])
        )
        return cls(
            orientation=orientation,
            position=position,
            mode=ReconstructionMode.KNOWN_ORIENTATION,
            current_zone_axis=current_zone_axis,
            ambiguity=analyze_ambiguity(
                orientation.phase,
                zone_cartesian,
                rotation_calibrated=True,
                reconstruction_note=(
                    "Orientation supplied directly, so no reconstruction ambiguity "
                    "arises; any residual ambiguity is crystallographic only."
                ),
            ),
            orientation_uncertainty_deg=orientation_uncertainty_deg,
            notes=("Orientation supplied directly rather than reconstructed.",),
            provenance=provenance,
        )

    @classmethod
    def from_two_zone_axes(
        cls,
        first: ZoneAxisObservation,
        second: ZoneAxisObservation,
        stage: StageModel,
        *,
        resolve_senses: bool = True,
        orientation_uncertainty_deg: float = 0.5,
        provenance: ProvenanceRecord | None = None,
    ) -> CurrentState:
        """Reconstruct from two indexed zone axes — the recommended path.

        Purpose
        -------
        Determines the crystal-to-holder orientation using **only** zone-axis
        identities and stage readouts. No diffraction rotation, no parity bit, no
        detector model: the single hardest calibration constant in the problem is
        simply not required.

        When to use
        -----------
        Whenever two zones have been indexed, which after any real session of
        Kikuchi-band chasing is almost always. Prefer it over
        :meth:`from_pattern_solution` even when a rotation calibration exists,
        because it is immune to that calibration being stale.

        Algorithm
        ---------
        By equation (M) the beam direction in holder coordinates at each position
        is ``U n_i``. Two non-parallel correspondences determine a rotation
        uniquely (Wahba's problem), solved here by orthogonal Procrustes on the
        two right-handed triads. When ``resolve_senses`` is set, all four sign
        combinations of the two zone axes are tried and the one whose interzonal
        angle matches the stage-predicted angle is chosen — the sense of an
        indexed zone axis is not determined by a single pattern, so testing it is
        not optional.

        Parameters
        ----------
        first, second : ZoneAxisObservation
            Must be for the same phase and must not be parallel.
        stage : StageModel
        resolve_senses : bool, default True
            Test all four sign combinations. Disable only when the senses are
            independently known.
        orientation_uncertainty_deg : float, default 0.5
            One-sigma scale, propagated into solution uncertainties.

        Returns
        -------
        CurrentState
            With :attr:`consistency_residual_deg` set — inspect it. A value well
            above the indexing error means one of the inputs is wrong.

        Notes
        -----
        A residual two-fold ambiguity survives: flipping *both* zone-axis senses
        admits a second proper rotation, related by a 180 degree turn about
        ``n1 x n2``. That ambiguity is **harmless when the two-fold is itself a
        crystal symmetry**, which this method checks and records.
        """

        if first.zone_axis.phase != second.zone_axis.phase:
            raise ValueError(
                "Both zone-axis observations must belong to the same phase."
            )
        phase = first.zone_axis.phase
        beam_first = stage.beam_direction(
            first.position.alpha_deg, first.position.beta_deg
        )
        beam_second = stage.beam_direction(
            second.position.alpha_deg, second.position.beta_deg
        )
        stage_angle_deg = math.degrees(
            math.acos(max(-1.0, min(1.0, float(np.dot(beam_first, beam_second)))))
        )

        sign_options: tuple[tuple[float, float], ...] = (
            ((1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0))
            if resolve_senses
            else ((1.0, 1.0),)
        )

        best: tuple[float, np.ndarray, tuple[float, float]] | None = None
        for sign_first, sign_second in sign_options:
            n1 = sign_first * first.unit_vector
            n2 = sign_second * second.unit_vector
            crystal_angle_deg = math.degrees(
                math.acos(max(-1.0, min(1.0, float(np.dot(n1, n2)))))
            )
            residual = abs(crystal_angle_deg - stage_angle_deg)
            if best is None or residual < best[0]:
                source = _triad(n1, n2)
                target = _triad(beam_first, beam_second)
                best = (residual, target @ source.T, (sign_first, sign_second))

        assert best is not None
        residual_deg, matrix, signs = best

        orientation = Orientation.from_matrix(
            matrix,
            specimen_frame=HOLDER_FRAME,
            phase=phase,
            crystal_frame=phase.crystal_frame,
        )
        flip_axis = np.cross(first.unit_vector, second.unit_vector)
        harmless = _two_fold_is_symmetry(phase, flip_axis)
        notes = []
        if residual_deg > _CONSISTENCY_WARNING_DEG:
            notes.append(
                f"WARNING: the interzonal-angle residual is {residual_deg:.2f} deg, far "
                "beyond indexing error. The crystallographic angle between the two "
                "indexed zones disagrees with the angle the stage model predicts "
                "between the two positions, so at least one input is wrong. In order of "
                "likelihood: a reversed stage sign convention, a mis-indexed zone, "
                "swapped stage readings, or a bent specimen. This reconstruction should "
                "not be used until that is resolved."
            )
        notes.extend(
            [
                "Reconstructed from two indexed zone axes; no diffraction-rotation or "
                "parity calibration was used or needed.",
                f"Sense assignment chosen: first {signs[0]:+.0f}, second {signs[1]:+.0f}, "
                f"on an interzonal-angle residual of {residual_deg:.3f} deg "
                f"(crystallographic angle versus stage-predicted angle).",
            ]
        )
        if harmless:
            notes.append(
                "The residual two-fold about the common normal of the two zone axes is "
                "itself a crystal symmetry, so the two-zone sign ambiguity is harmless "
                "here and no experimental discrimination is needed."
            )
        else:
            notes.append(
                "The residual two-fold about the common normal of the two zone axes is "
                "NOT a crystal symmetry, so a second orientation hypothesis exists; "
                "resolve it by a known small tilt excursion or by indexing a third zone."
            )
        return cls(
            orientation=orientation,
            position=second.position,
            mode=ReconstructionMode.TWO_ZONE_AXES,
            current_zone_axis=second.zone_axis,
            ambiguity=analyze_ambiguity(
                phase,
                signs[1] * second.unit_vector,
                rotation_calibrated=True,
                reconstruction_note=(
                    "Two-zone reconstruction: the diffraction rotation was not used, so "
                    "the instrumental-rotation ambiguity does not arise."
                ),
            ),
            consistency_residual_deg=residual_deg,
            orientation_uncertainty_deg=orientation_uncertainty_deg,
            notes=tuple(notes),
            provenance=provenance,
        )

    @classmethod
    def from_pattern_solution(
        cls,
        crystal_to_pattern: Rotation | ArrayLike,
        position: StagePosition,
        stage: StageModel,
        phase: Phase,
        *,
        zone_axis: ZoneAxis | None = None,
        orientation_uncertainty_deg: float = 0.5,
        provenance: ProvenanceRecord | None = None,
    ) -> CurrentState:
        """Reconstruct from one indexed pattern plus a calibrated pattern rotation.

        Purpose
        -------
        The single-pattern path. Composes the crystal-to-pattern rotation that
        pattern indexing supplies with the calibrated diffraction rotation and
        the inverse stage rotation:

            U = R_stage(alpha, beta)^T F Rz(phi_D) R_pattern<-crystal.

        When to use
        -----------
        When only one zone has been indexed **and** a current diffraction-rotation
        calibration exists for the camera length in use. Prefer
        :meth:`from_two_zone_axes` otherwise — and note that this method raises
        rather than guessing when the calibration is missing, because a guessed
        pattern azimuth is exactly the error that sends an operator tilting
        backwards with a clean-looking residual.

        Parameters
        ----------
        crystal_to_pattern : Rotation or array_like
            Maps crystal Cartesian vectors into the pattern frame, whose ``z``
            axis is the zone axis toward the viewer. This is precisely
            ``pytex.diffraction.solving.PatternSolution.orientation``.
        position : StagePosition
        stage : StageModel
            Supplies the kinematics and the calibration.
        phase : Phase
        zone_axis : ZoneAxis, optional
            The indexed zone axis, for reporting and for the ambiguity analysis.

        Returns
        -------
        CurrentState

        Raises
        ------
        ValueError
            If the stage calibration has no diffraction rotation, or if the
            resulting orientation is improper — which is the free, exact
            self-check that catches a wrong parity setting.
        """

        calibration = stage.calibration
        if not calibration.is_rotation_calibrated:
            raise ValueError(
                "Single-pattern reconstruction needs a calibrated diffraction rotation, "
                "and this StageCalibration has none. Either run the two-excursion "
                "calibration (pytex.tem.calibration.calibrate_from_tilt_excursions) or "
                "use CurrentState.from_two_zone_axes, which needs no rotation "
                "calibration at all."
            )
        matrix_pattern = (
            crystal_to_pattern.as_matrix()
            if isinstance(crystal_to_pattern, Rotation)
            else as_float_array(crystal_to_pattern, shape=(3, 3))
        )
        phi = math.radians(float(calibration.diffraction_rotation_deg or 0.0))
        parity = (
            np.diag([1.0, -1.0, 1.0]).astype(np.float64)
            if calibration.pattern_is_mirrored
            else np.eye(3, dtype=np.float64)
        )
        stage_matrix = stage.rotation_matrix(position.alpha_deg, position.beta_deg)
        matrix = stage_matrix.T @ parity @ rotation_z(phi) @ matrix_pattern

        determinant = float(np.linalg.det(matrix))
        if determinant < 0.0:
            raise ValueError(
                "The reconstructed crystal-to-holder rotation is improper "
                f"(determinant {determinant:+.3f}), which means the parity of the stored "
                "pattern does not match the StageCalibration. Flip "
                "StageCalibration.pattern_is_mirrored, or re-run the two-excursion "
                "calibration, which measures the parity from the sign of the "
                "alpha-versus-beta motion azimuths."
            )

        orientation = Orientation.from_matrix(
            matrix,
            specimen_frame=HOLDER_FRAME,
            phase=phase,
            crystal_frame=phase.crystal_frame,
        )
        zone_cartesian = (
            zone_axis.unit_vector
            if zone_axis is not None
            else matrix.T @ stage_matrix.T @ BEAM_AXIS_LABORATORY
        )
        return cls(
            orientation=orientation,
            position=position,
            mode=ReconstructionMode.SINGLE_PATTERN,
            current_zone_axis=zone_axis,
            ambiguity=analyze_ambiguity(
                phase,
                zone_cartesian,
                rotation_calibrated=True,
                reconstruction_note=(
                    "Single-pattern reconstruction using a diffraction rotation of "
                    f"{calibration.diffraction_rotation_deg:+.2f} deg. An error of "
                    "180 deg in that constant would negate both tilt angles while still "
                    "reporting a zero residual, so confirm it with a known small tilt "
                    "excursion before a long excursion."
                ),
            ),
            orientation_uncertainty_deg=orientation_uncertainty_deg,
            notes=(
                "Reconstructed from a single indexed pattern; the answer depends on the "
                "diffraction-rotation calibration being current for this camera length.",
            ),
            provenance=provenance,
        )

    @classmethod
    def from_zone_axes(
        cls,
        observations: tuple[ZoneAxisObservation, ...],
        stage: StageModel,
        *,
        orientation_uncertainty_deg: float | None = None,
        provenance: ProvenanceRecord | None = None,
    ) -> CurrentState:
        """Least-squares reconstruction over three or more indexed zone axes.

        Purpose
        -------
        The over-determined path. With more correspondences than the two a
        rotation needs, the fit both averages out indexing error and *measures*
        it: the residual scatter is an empirical uncertainty rather than an
        assumed one, and a systematically large value bounds the specimen bending
        that would otherwise silently violate the rigid-crystal assumption.

        Parameters
        ----------
        observations : tuple of ZoneAxisObservation
            Three or more, for one phase. Weights are honoured.
        stage : StageModel
        orientation_uncertainty_deg : float, optional
            Override the measured scatter. Leave unset to use the fit residual,
            which is the honest figure.

        Returns
        -------
        CurrentState
        """

        if len(observations) < 3:
            raise ValueError(
                "from_zone_axes needs at least three observations; use "
                "from_two_zone_axes for exactly two."
            )
        phase = observations[0].zone_axis.phase
        # Compared rather than collected into a set: Phase is a rich value object
        # and is deliberately not hashable, so equality is the available test.
        if any(observation.zone_axis.phase != phase for observation in observations):
            raise ValueError("All zone-axis observations must belong to one phase.")

        beams = np.stack(
            [
                stage.beam_direction(
                    observation.position.alpha_deg, observation.position.beta_deg
                )
                for observation in observations
            ]
        )
        axes = np.stack([observation.unit_vector for observation in observations])
        weights = np.asarray(
            [observation.weight for observation in observations], dtype=np.float64
        )

        # Resolve each sense against the first by comparing crystallographic and
        # stage-predicted angles: the sense of an indexed zone axis is not
        # determined by its own pattern.
        signs = np.ones(len(observations), dtype=np.float64)
        for index in range(1, len(observations)):
            stage_angle = float(np.dot(beams[0], beams[index]))
            crystal_angle = float(np.dot(axes[0], axes[index]))
            if abs(crystal_angle - stage_angle) > abs(-crystal_angle - stage_angle):
                signs[index] = -1.0
        axes = axes * signs[:, None]

        matrix = _kabsch_rotation(axes * weights[:, None], beams * weights[:, None])
        predicted = np.einsum("ij,nj->ni", matrix, axes)
        # Not arccos of the dot product: a fit this close to exact is precisely
        # where that loses half its digits, and the scatter reported here is
        # meant to be believable at zero. See pytex.core._angles.
        residuals_deg = np.degrees(angle_between_unit_vectors_rad(predicted, beams))
        measured_scatter = float(np.sqrt(np.mean(residuals_deg**2)))

        orientation = Orientation.from_matrix(
            matrix,
            specimen_frame=HOLDER_FRAME,
            phase=phase,
            crystal_frame=phase.crystal_frame,
        )
        return cls(
            orientation=orientation,
            position=observations[-1].position,
            mode=ReconstructionMode.MULTI_ZONE_FIT,
            current_zone_axis=observations[-1].zone_axis,
            ambiguity=analyze_ambiguity(
                phase,
                signs[-1] * observations[-1].unit_vector,
                rotation_calibrated=True,
                reconstruction_note=(
                    f"Least-squares fit over {len(observations)} indexed zone axes; the "
                    "diffraction rotation was not used."
                ),
            ),
            consistency_residual_deg=float(np.max(residuals_deg)),
            orientation_uncertainty_deg=(
                measured_scatter
                if orientation_uncertainty_deg is None
                else orientation_uncertainty_deg
            ),
            notes=(
                f"Fitted over {len(observations)} zone axes; root-mean-square residual "
                f"{measured_scatter:.3f} deg, worst {float(np.max(residuals_deg)):.3f} deg. "
                "A systematically large residual bounds specimen bending, which is the "
                "failure mode of the rigid-crystal assumption.",
            ),
            provenance=provenance,
        )

    # ------------------------------------------------------------------ #
    # Reporting
    # ------------------------------------------------------------------ #

    def describe(self) -> str:
        """Convention-explicit prose: what was used, what came out, what is open."""

        mode_text = {
            ReconstructionMode.SINGLE_PATTERN: (
                "one indexed pattern plus a calibrated diffraction rotation"
            ),
            ReconstructionMode.TWO_ZONE_AXES: (
                "two indexed zone axes at two stage positions, needing no "
                "diffraction-rotation calibration"
            ),
            ReconstructionMode.MULTI_ZONE_FIT: "a least-squares fit over several indexed zones",
            ReconstructionMode.KNOWN_ORIENTATION: "a directly supplied orientation",
        }[self.mode]
        head = (
            f"Crystal-to-holder orientation for {self.phase.name} determined from "
            f"{mode_text}, with the stage at {self.position}. The holder frame is the "
            "specimen-domain frame, so this orientation is a standard crystal-to-specimen "
            "mapping and the mounting rotation — which is not separately observable — is "
            "absorbed into it."
        )
        consistency = ""
        if self.consistency_residual_deg is not None:
            consistency = (
                f" Interzonal-angle consistency residual "
                f"{self.consistency_residual_deg:.3f} deg: the crystallographic angle "
                "between the indexed zones against the angle the stage model predicts, "
                "a check that uses no calibration at all."
            )
        uncertainty = (
            f" Orientation uncertainty {self.orientation_uncertainty_deg:.2f} deg "
            "(one sigma)."
        )
        return " ".join(
            [head + consistency + uncertainty, *self.notes, self.ambiguity.describe()]
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "phase": self.phase.name,
            "mode": self.mode.value,
            "stage_position": {
                "alpha_deg": self.position.alpha_deg,
                "beta_deg": self.position.beta_deg,
            },
            "crystal_to_holder_matrix": self.matrix.tolist(),
            "current_zone_axis": (
                None
                if self.current_zone_axis is None
                else [int(value) for value in self.current_zone_axis.indices]
            ),
            "consistency_residual_deg": self.consistency_residual_deg,
            "orientation_uncertainty_deg": self.orientation_uncertainty_deg,
            "ambiguity": self.ambiguity.to_json_dict(),
            "notes": list(self.notes),
        }


def _two_fold_is_symmetry(phase: Phase, axis: ArrayLike) -> bool:
    """Whether the 180 degree rotation about ``axis`` is a crystal symmetry.

    The check that decides whether the two-zone sign ambiguity matters. For a
    cubic crystal with zones ``[001]`` and ``[110]`` the axis is ``<1-10>`` and
    the answer is yes, so the ambiguity vanishes without any experiment.
    """

    vector = np.asarray(axis, dtype=np.float64)
    norm = float(np.linalg.norm(vector))
    if norm < 1e-12:
        return False
    unit = vector / norm
    outer = np.outer(unit, unit)
    two_fold = 2.0 * outer - np.eye(3)
    operators = np.asarray(phase.symmetry.operators, dtype=np.float64)
    return bool(np.any(np.all(np.abs(operators - two_fold) < 1e-8, axis=(1, 2))))
