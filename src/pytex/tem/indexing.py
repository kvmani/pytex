"""From an indexed SAED pattern to the crystal orientation.

Pattern indexing answers *which zone axis am I looking down*. This module
answers the next question, which is the one texture and microstructure work
actually needs: *given that indexing and the holder tilts, what is the crystal
orientation?*

The composition
---------------
`pytex.diffraction.solving.PatternSolution.orientation` is the crystal-to-pattern
rotation ``R_pattern<-crystal``, carrying crystal Cartesian vectors into the
stored pattern's own axes. Chaining it with the diffraction rotation and the
inverse stage rotation gives the crystal-to-holder orientation,

    U = R_stage(alpha, beta)^T  F  Rz(phi_D)  R_pattern<-crystal,

which — because the holder frame is declared to be the specimen-domain frame for
TEM work — *is* an ordinary `pytex.core.orientation.Orientation`, reportable as
Bunge Euler angles and comparable with EBSD orientations without further
conversion.

Why one pattern is not enough, and what to do about it
------------------------------------------------------
That composition needs ``phi_D``, the diffraction rotation. It is not in the file
metadata, it drifts with the lens history, and an error in it is **not** absorbed
by crystal symmetry: it rotates the reported orientation bodily about the beam
axis. A 180-degree error yields a perfectly self-consistent orientation that is
wrong by 180 degrees about the beam. So :func:`orientation_from_indexed_pattern`
refuses to guess, and raises when no calibration is available.

The way out is :func:`orientation_from_indexed_patterns`, which takes **two or
more** indexed patterns at different stage positions and determines the
orientation *and* ``phi_D`` together, from the data alone. Two zone axes fix the
orientation by themselves; each pattern's in-plane indexing then over-determines
the diffraction rotation, and the spread across patterns is a direct measure of
whether the whole chain hangs together.

That last quantity is the useful one. Writing ``M_i = R_stage_i U R_i^T``, the
composition demands that every ``M_i`` be a rotation **about the beam axis** and
nothing else. How far ``M_i`` tilts the beam axis is therefore a residual that no
choice of ``phi_D`` can absorb, and it indicts — in order of likelihood — a
mirrored stored pattern, a mis-indexed reflection, a reversed stage sign
convention, or a bent specimen.

See ``docs/architecture/tem_tilt_navigation_foundation.md`` sections 5 and 8.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Phase, ZoneAxis
from pytex.core.orientation import Orientation, Rotation
from pytex.core.provenance import ProvenanceRecord
from pytex.tem.ambiguity import AmbiguityReport, analyze_ambiguity
from pytex.tem.reconstruction import HOLDER_FRAME, _kabsch_rotation
from pytex.tem.stage import (
    BEAM_AXIS_LABORATORY,
    StageCalibration,
    StageModel,
    StagePosition,
    rotation_z,
)

__all__ = [
    "INDEXED_ORIENTATION_SCHEMA",
    "IndexedOrientation",
    "IndexedPatternObservation",
    "MultiPatternOrientation",
    "orientation_from_indexed_pattern",
    "orientation_from_indexed_patterns",
    "orientations_from_pattern_report",
]

#: Schema identifier of the indexed-orientation payload.
INDEXED_ORIENTATION_SCHEMA = "pytex.indexed_orientation/1"

#: Beam-axis deviation, in degrees, above which a multi-pattern fit is indicted.
#:
#: The composition requires every per-pattern residual rotation to be about the
#: beam axis alone. A degree is comfortably above the error a well-indexed
#: pattern contributes and comfortably below what a mirrored pattern or a
#: reversed sign convention produces, so the test separates them cleanly.
BEAM_DEVIATION_WARNING_DEG = 1.0


def _pattern_matrix(rotation: Rotation | np.ndarray) -> np.ndarray:
    """The crystal-to-pattern rotation as a plain, validated 3x3 matrix.

    Rejects an improper matrix. Pattern indexing builds right-handed triads from
    both the observed and the calculated vectors, so its output is always a
    proper rotation; a reflection here means the caller composed a mirror into it
    by hand, which would silently corrupt the orientation rather than being
    caught downstream.

    Note that this is *not* how a mirrored stored pattern presents. That case
    still yields a proper matrix — the in-plane mirror is accompanied by a
    reversal of the derived pattern normal, and the pair is a rotation — and it
    is detected instead by the beam-axis deviation in
    :func:`orientation_from_indexed_patterns`.
    """

    matrix = (
        np.asarray(rotation.as_matrix(), dtype=np.float64)
        if isinstance(rotation, Rotation)
        else as_float_array(rotation, shape=(3, 3))
    )
    determinant = float(np.linalg.det(matrix))
    if determinant <= 0.0:
        raise ValueError(
            f"The crystal-to-pattern matrix is improper (determinant {determinant:+.3f}); "
            "it must be a proper rotation. Pattern indexing always returns one, so a "
            "reflection here means a mirror was composed in by hand. If the stored "
            "pattern is mirrored, record that in StageCalibration.pattern_is_mirrored "
            "rather than folding it into the matrix."
        )
    return matrix


@dataclass(frozen=True, slots=True)
class IndexedPatternObservation:
    """One indexed pattern recorded at a known stage position.

    Purpose
    -------
    The input datum of the multi-pattern path: *this pattern, indexed this way,
    was recorded with the stage reading these two angles.* Everything needed to
    determine an orientation without any prior calibration, provided there are
    two of them.

    Attributes
    ----------
    crystal_to_pattern : Rotation or np.ndarray
        Carries crystal Cartesian vectors into the stored pattern's axes. This
        is `PatternSolution.orientation`.
    zone_axis : ZoneAxis
        The indexed zone axis. Its sense may be wrong; the fit tests both.
    position : StagePosition
    label : str
        Free-form identifier carried into reports.
    """

    crystal_to_pattern: Rotation | np.ndarray
    zone_axis: ZoneAxis
    position: StagePosition
    label: str = ""

    @property
    def matrix(self) -> np.ndarray:
        """The crystal-to-pattern rotation as a 3x3 matrix."""

        return _pattern_matrix(self.crystal_to_pattern)

    @property
    def unit_vector(self) -> np.ndarray:
        """The indexed zone axis as a Cartesian crystal-frame unit vector."""

        return np.asarray(self.zone_axis.unit_vector, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class IndexedOrientation:
    """The crystal orientation implied by one indexed pattern and a stage position.

    Purpose
    -------
    What a user wants out of indexing when the goal is microstructure rather than
    phase identification: an `Orientation` object, in Bunge Euler angles,
    comparable with EBSD data and usable in every texture surface in the library.

    The object also carries what the single pattern could **not** determine, so
    that a caller cannot mistake a reported orientation for a uniquely determined
    one.

    Attributes
    ----------
    orientation : Orientation
        Crystal-to-holder, and therefore crystal-to-specimen: the holder frame is
        the specimen-domain frame for TEM work.
    zone_axis : ZoneAxis
        The zone axis that was on the beam.
    position : StagePosition
    ambiguity : AmbiguityReport
        What the observation leaves undetermined. Inspect
        :attr:`AmbiguityReport.is_unique` before treating the orientation as the
        answer.
    equivalent_orientations : tuple of Orientation
        One per ambiguity family beyond the first: the physically distinct
        orientations the same pattern also admits. Empty when the orientation is
        uniquely determined, which for a centrosymmetric crystal it usually is.
    zone_axis_residual_deg : float
        Angle between the reported orientation's image of the zone axis and the
        beam. A forward validation of the composition, not a restatement of it.
    diffraction_rotation_deg : float
        The value used, for the record.
    notes : tuple of str
    provenance : ProvenanceRecord or None
    """

    orientation: Orientation
    zone_axis: ZoneAxis
    position: StagePosition
    ambiguity: AmbiguityReport
    equivalent_orientations: tuple[Orientation, ...] = ()
    zone_axis_residual_deg: float = 0.0
    diffraction_rotation_deg: float = 0.0
    notes: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "notes", tuple(self.notes))
        object.__setattr__(
            self, "equivalent_orientations", tuple(self.equivalent_orientations)
        )

    @property
    def phase(self) -> Phase:
        """The phase this orientation belongs to."""

        if self.orientation.phase is None:  # pragma: no cover - constructed with one
            raise ValueError("IndexedOrientation.orientation must carry a phase.")
        return self.orientation.phase

    @property
    def matrix(self) -> np.ndarray:
        """The crystal-to-specimen rotation matrix."""

        return np.asarray(self.orientation.as_matrix(), dtype=np.float64)

    @property
    def euler_bunge_deg(self) -> tuple[float, float, float]:
        """Bunge ``(phi1, Phi, phi2)`` Euler angles in degrees.

        The representation EBSD and texture work quote, so an orientation
        measured in the TEM can be compared with one measured by EBSD without a
        convention conversion in between.
        """

        return self.orientation.rotation.to_bunge_euler(degrees=True)

    @property
    def is_unique(self) -> bool:
        """Whether the pattern determines this orientation without alternatives."""

        return self.ambiguity.is_unique

    def in_frame(self, specimen_frame: ReferenceFrame) -> Orientation:
        """Re-express the orientation against a differently-named specimen frame.

        Purpose
        -------
        The holder frame and a rolling-geometry ``RD/TD/ND`` frame are the same
        specimen domain under different labels. When a TEM orientation is to be
        pooled with sheet-texture data, the frame it declares should say so.

        This **relabels** and does not rotate: it asserts that the target frame's
        axes coincide with the holder's. If the specimen was mounted with a known
        rotation relative to the sheet axes, compose that rotation yourself
        first — silently absorbing it here is exactly the kind of hidden
        convention this library exists to prevent.
        """

        return Orientation.from_matrix(
            self.matrix,
            specimen_frame=specimen_frame,
            phase=self.phase,
            crystal_frame=self.phase.crystal_frame,
        )

    def describe(self) -> str:
        """Convention-explicit prose: the orientation, and what is not determined."""

        phi1, capital_phi, phi2 = self.euler_bunge_deg
        head = (
            f"Crystal orientation of {self.phase.name} from the indexed pattern down "
            f"{self.zone_axis.indices.tolist()} with the stage at {self.position}: "
            f"Bunge (phi1, Phi, phi2) = ({phi1:.2f}, {capital_phi:.2f}, {phi2:.2f}) deg. "
            "The convention is crystal-to-specimen, and the specimen frame is the "
            "holder frame, so this orientation is directly comparable with an EBSD "
            "measurement expressed in the same frame."
        )
        validation = (
            f" Forward validation: the reported orientation places the indexed zone "
            f"axis {self.zone_axis_residual_deg:.4f} deg from the beam."
        )
        calibration = (
            f" Diffraction rotation used: {self.diffraction_rotation_deg:+.2f} deg. "
            "An error in that constant rotates this orientation bodily about the beam "
            "axis and is not absorbed by crystal symmetry."
        )
        alternatives = ""
        if self.equivalent_orientations:
            alternatives = (
                f" {len(self.equivalent_orientations)} alternative orientation(s) are "
                "equally consistent with this pattern and are reported alongside; the "
                "pattern does not choose between them."
            )
        notes = (" " + " ".join(self.notes)) if self.notes else ""
        return head + validation + calibration + alternatives + notes + " " + (
            self.ambiguity.describe()
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        phi1, capital_phi, phi2 = self.euler_bunge_deg
        return {
            "schema": INDEXED_ORIENTATION_SCHEMA,
            "phase": self.phase.name,
            "zone_axis": [int(value) for value in self.zone_axis.indices],
            "stage_position": {
                "alpha_deg": self.position.alpha_deg,
                "beta_deg": self.position.beta_deg,
            },
            "orientation_matrix": self.matrix.tolist(),
            "euler_bunge_deg": {"phi1": phi1, "Phi": capital_phi, "phi2": phi2},
            "zone_axis_residual_deg": self.zone_axis_residual_deg,
            "diffraction_rotation_deg": self.diffraction_rotation_deg,
            "is_unique": self.is_unique,
            "equivalent_orientation_matrices": [
                np.asarray(orientation.as_matrix()).tolist()
                for orientation in self.equivalent_orientations
            ],
            "ambiguity": self.ambiguity.to_json_dict(),
            "notes": list(self.notes),
        }


def _zone_axis_residual_deg(
    stage: StageModel, matrix: np.ndarray, zone_axis: ZoneAxis, position: StagePosition
) -> float:
    """Angle between the orientation's image of the zone axis and the beam.

    Recomputed through the stage model from the crystal direction, so a mistake
    in the composition cannot validate itself.
    """

    achieved = (
        stage.rotation_matrix(position.alpha_deg, position.beta_deg)
        @ matrix
        @ np.asarray(zone_axis.unit_vector, dtype=np.float64)
    )
    cosine = abs(float(np.dot(normalize_vector(achieved), BEAM_AXIS_LABORATORY)))
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))


def orientation_from_indexed_pattern(
    crystal_to_pattern: Rotation | np.ndarray,
    zone_axis: ZoneAxis,
    position: StagePosition,
    stage: StageModel,
    *,
    specimen_frame: ReferenceFrame | None = None,
    provenance: ProvenanceRecord | None = None,
) -> IndexedOrientation:
    """Crystal orientation from one indexed pattern and the holder tilts.

    Purpose
    -------
    Turns the output of SAED indexing into an orientation. Given the
    crystal-to-pattern rotation that indexing produced, the zone axis it found,
    and the stage angles at which the pattern was recorded, this returns an
    `Orientation` — in Bunge Euler angles, in the crystal-to-specimen convention
    the rest of the library uses, and therefore directly comparable with EBSD.

    When to use
    -----------
    After `pytex.diffraction.solving.solve_saed_pattern` (or
    `index_saed_pattern`) has produced a solution you trust, when the goal is the
    orientation rather than the phase or the indices.

    **A diffraction-rotation calibration is required**, and this raises without
    one rather than guessing. If you have no calibration, use
    :func:`orientation_from_indexed_patterns` with a second pattern at a
    different stage position: two patterns determine the orientation *and* the
    calibration together.

    Parameters
    ----------
    crystal_to_pattern : Rotation or array_like
        `PatternSolution.orientation`: crystal Cartesian vectors into the stored
        pattern's axes.
    zone_axis : ZoneAxis
        `PatternSolution.zone_axis`. Supplies the phase and is used for forward
        validation.
    position : StagePosition
        The holder alpha and beta at which the pattern was recorded.
    stage : StageModel
        Supplies the kinematics and the diffraction-rotation calibration.
    specimen_frame : ReferenceFrame, optional
        Declare the orientation against a differently-named specimen-domain
        frame, for example a rolling ``RD/TD/ND`` frame. Relabels only; see
        :meth:`IndexedOrientation.in_frame`.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    IndexedOrientation
        Carrying the orientation, its Euler angles, the forward-validated
        residual, and the ambiguity the single pattern leaves.

    Raises
    ------
    ValueError
        If the stage has no diffraction-rotation calibration, or if the
        resulting orientation is improper — the free parity self-check.

    Notes
    -----
    A single pattern determines the orientation only up to the rotations of the
    Laue class that map the zone plane to itself. For a centrosymmetric crystal
    every one of those is a crystal symmetry and nothing is lost; for ten of the
    thirty-two point groups it is not, and the alternatives are reported in
    :attr:`IndexedOrientation.equivalent_orientations`.

    Examples
    --------
    A crystal whose ``[001]`` is on the beam at zero tilt, with a zero
    diffraction rotation and an identity pattern rotation, is at the identity
    orientation — the trivial case that fixes the sign conventions of the whole
    composition.
    """

    calibration: StageCalibration = stage.calibration
    if not calibration.is_rotation_calibrated:
        raise ValueError(
            "Reporting a crystal orientation from a single indexed pattern needs a "
            "calibrated diffraction rotation, and this StageCalibration has none. "
            "An assumed value rotates the reported orientation bodily about the beam "
            "axis, and crystal symmetry does not absorb that error. Either run "
            "pytex.tem.calibration.calibrate_from_tilt_excursions, or call "
            "orientation_from_indexed_patterns with a second pattern at a different "
            "stage position, which determines the orientation and the calibration "
            "together."
        )

    phase = zone_axis.phase
    rotation_deg = float(calibration.diffraction_rotation_deg or 0.0)
    parity = (
        np.diag([1.0, -1.0, 1.0]).astype(np.float64)
        if calibration.pattern_is_mirrored
        else np.eye(3, dtype=np.float64)
    )
    stage_matrix = stage.rotation_matrix(position.alpha_deg, position.beta_deg)
    matrix = (
        stage_matrix.T
        @ parity
        @ rotation_z(math.radians(rotation_deg))
        @ _pattern_matrix(crystal_to_pattern)
    )

    determinant = float(np.linalg.det(matrix))
    if determinant < 0.0:
        raise ValueError(
            "The orientation implied by this pattern is improper "
            f"(determinant {determinant:+.3f}), which means the parity recorded in the "
            "StageCalibration does not match the stored pattern. Flip "
            "StageCalibration.pattern_is_mirrored, or re-run the two-excursion "
            "calibration, which measures the parity from the sign of the "
            "alpha-versus-beta motion azimuths."
        )

    target_frame = specimen_frame or HOLDER_FRAME
    orientation = Orientation.from_matrix(
        matrix,
        specimen_frame=target_frame,
        phase=phase,
        crystal_frame=phase.crystal_frame,
    )
    ambiguity = analyze_ambiguity(
        phase,
        zone_axis.unit_vector,
        rotation_calibrated=True,
        reconstruction_note=(
            "Orientation derived from a single indexed pattern using a diffraction "
            f"rotation of {rotation_deg:+.2f} deg."
        ),
    )
    equivalents = tuple(
        Orientation.from_matrix(
            matrix @ family.operator,
            specimen_frame=target_frame,
            phase=phase,
            crystal_frame=phase.crystal_frame,
        )
        for family in ambiguity.families[1:]
    )

    return IndexedOrientation(
        orientation=orientation,
        zone_axis=zone_axis,
        position=position,
        ambiguity=ambiguity,
        equivalent_orientations=equivalents,
        zone_axis_residual_deg=_zone_axis_residual_deg(
            stage, matrix, zone_axis, position
        ),
        diffraction_rotation_deg=rotation_deg,
        notes=(
            "Single-pattern orientation: the answer inherits the diffraction-rotation "
            "calibration, so confirm that constant is current for this camera length "
            "before quoting the orientation.",
        ),
        provenance=provenance,
    )


def orientations_from_pattern_report(
    report: Any,
    position: StagePosition,
    stage: StageModel,
    *,
    max_solutions: int = 3,
    specimen_frame: ReferenceFrame | None = None,
) -> tuple[IndexedOrientation, ...]:
    """Orientations for the ranked solutions of a pattern-solution report.

    Purpose
    -------
    A solved pattern may admit more than one indexing — a different phase, or a
    different zone. Each implies a different orientation, and presenting only the
    best would hide that the *indexing* was the uncertain step rather than the
    orientation arithmetic.

    Parameters
    ----------
    report : PatternSolutionReport
        From `pytex.diffraction.solving.solve_saed_pattern`. Typed loosely to
        keep `pytex.tem` free of a hard import back into the solver.
    position : StagePosition
    stage : StageModel
    max_solutions : int, default 3
    specimen_frame : ReferenceFrame, optional

    Returns
    -------
    tuple of IndexedOrientation
        In the report's ranking order. Empty when the pattern was not solved.
    """

    solutions = getattr(report, "solutions", ())
    return tuple(
        orientation_from_indexed_pattern(
            solution.orientation,
            solution.zone_axis,
            position,
            stage,
            specimen_frame=specimen_frame,
        )
        for solution in solutions[:max_solutions]
    )


@dataclass(frozen=True, slots=True)
class MultiPatternOrientation:
    """Orientation and diffraction rotation determined jointly from several patterns.

    Purpose
    -------
    The self-calibrating result. Two or more indexed patterns at different stage
    positions determine the crystal orientation *and* the diffraction rotation
    from the data alone, so no prior calibration is needed and the answer cannot
    inherit a stale one.

    Attributes
    ----------
    orientation : Orientation
        Crystal-to-specimen, with the holder frame as the specimen frame.
    diffraction_rotation_deg : float
        Determined, not assumed. Mean over the supplied patterns.
    diffraction_rotation_scatter_deg : float
        Spread of the per-pattern estimates. Small means the patterns agree on
        one instrument constant, which is a strong check that the whole chain is
        right.
    beam_deviation_deg : float
        The residual no choice of diffraction rotation can absorb: how far each
        pattern's implied correction tilts the beam axis. **This is the number to
        read.** A large value indicts a mirrored stored pattern, a mis-indexed
        reflection, a reversed stage sign convention, or a bent specimen.
    interzonal_residual_deg : float
        Crystallographic interzonal angle against the angle the stage model
        predicts; a second, independent consistency check that uses no
        calibration at all.
    zone_axis_signs : tuple of int
        The sense chosen for each indexed zone axis. A single pattern does not
        determine it; the fit does.
    ambiguity : AmbiguityReport
    is_consistent : bool
    observations : tuple of IndexedPatternObservation
    notes : tuple of str
    provenance : ProvenanceRecord or None
    """

    orientation: Orientation
    diffraction_rotation_deg: float
    diffraction_rotation_scatter_deg: float
    beam_deviation_deg: float
    interzonal_residual_deg: float
    zone_axis_signs: tuple[int, ...]
    ambiguity: AmbiguityReport
    is_consistent: bool
    observations: tuple[IndexedPatternObservation, ...]
    notes: tuple[str, ...] = ()
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "notes", tuple(self.notes))
        object.__setattr__(self, "observations", tuple(self.observations))
        object.__setattr__(self, "zone_axis_signs", tuple(self.zone_axis_signs))

    @property
    def phase(self) -> Phase:
        """The phase these patterns belong to."""

        if self.orientation.phase is None:  # pragma: no cover
            raise ValueError("MultiPatternOrientation.orientation must carry a phase.")
        return self.orientation.phase

    @property
    def matrix(self) -> np.ndarray:
        """The crystal-to-specimen rotation matrix."""

        return np.asarray(self.orientation.as_matrix(), dtype=np.float64)

    @property
    def euler_bunge_deg(self) -> tuple[float, float, float]:
        """Bunge ``(phi1, Phi, phi2)`` Euler angles in degrees."""

        return self.orientation.rotation.to_bunge_euler(degrees=True)

    def as_calibration(self, base: StageCalibration | None = None) -> StageCalibration:
        """The determined diffraction rotation, as a reusable calibration.

        Purpose
        -------
        Closes the loop: having measured the instrument constant from two
        patterns, subsequent single patterns can be turned into orientations
        directly. Refuses when the fit was inconsistent, because propagating a
        calibration that the data already contradict is worse than having none.
        """

        if not self.is_consistent:
            raise ValueError(
                f"This fit is not self-consistent (beam deviation "
                f"{self.beam_deviation_deg:.2f} deg), so its diffraction rotation must "
                "not be adopted as a calibration. Resolve the inconsistency first: the "
                "usual causes are a mirrored stored pattern, a mis-indexed reflection, "
                "a reversed stage sign convention, or a bent specimen."
            )
        template = base or StageCalibration()
        return StageCalibration(
            axes=template.axes,
            alpha_sign=template.alpha_sign,
            beta_sign=template.beta_sign,
            alpha_zero_deg=template.alpha_zero_deg,
            beta_zero_deg=template.beta_zero_deg,
            diffraction_rotation_deg=self.diffraction_rotation_deg,
            pattern_is_mirrored=template.pattern_is_mirrored,
            camera_length_mm=template.camera_length_mm,
            accelerating_voltage_kv=template.accelerating_voltage_kv,
            backlash_deg=template.backlash_deg,
            angular_uncertainty_deg=template.angular_uncertainty_deg,
            notes=(
                *template.notes,
                f"Diffraction rotation determined from {len(self.observations)} indexed "
                f"patterns, scatter {self.diffraction_rotation_scatter_deg:.2f} deg.",
            ),
        )

    def describe(self) -> str:
        """Convention-explicit prose: the orientation, the constant, and the checks."""

        phi1, capital_phi, phi2 = self.euler_bunge_deg
        head = (
            f"Crystal orientation of {self.phase.name} determined from "
            f"{len(self.observations)} indexed patterns with no prior calibration: "
            f"Bunge (phi1, Phi, phi2) = ({phi1:.2f}, {capital_phi:.2f}, {phi2:.2f}) deg, "
            "crystal-to-specimen with the holder as the specimen frame."
        )
        calibration = (
            f" The diffraction rotation was *determined* rather than assumed: "
            f"{self.diffraction_rotation_deg:+.2f} deg, with a scatter of "
            f"{self.diffraction_rotation_scatter_deg:.2f} deg across the patterns."
        )
        checks = (
            f" Consistency: the beam-axis deviation is {self.beam_deviation_deg:.3f} deg "
            "— the part of the residual that no diffraction rotation can absorb — and "
            f"the interzonal-angle residual is {self.interzonal_residual_deg:.3f} deg, "
            "which uses no calibration at all."
        )
        verdict = (
            " Both checks pass, so the orientation and the instrument constant can be "
            "relied on."
            if self.is_consistent
            else (
                " THE FIT IS NOT SELF-CONSISTENT. In order of likelihood: the stored "
                "pattern is mirrored, a reflection is mis-indexed, a stage sign "
                "convention is reversed, or the specimen is bent. Resolve that before "
                "quoting this orientation."
            )
        )
        signs = " Zone-axis senses chosen: " + ", ".join(
            f"{observation.label or index + 1}: {sign:+d}"
            for index, (observation, sign) in enumerate(
                zip(self.observations, self.zone_axis_signs, strict=True)
            )
        ) + "."
        notes = (" " + " ".join(self.notes)) if self.notes else ""
        return head + calibration + checks + verdict + signs + notes

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        phi1, capital_phi, phi2 = self.euler_bunge_deg
        return {
            "schema": INDEXED_ORIENTATION_SCHEMA,
            "phase": self.phase.name,
            "pattern_count": len(self.observations),
            "orientation_matrix": self.matrix.tolist(),
            "euler_bunge_deg": {"phi1": phi1, "Phi": capital_phi, "phi2": phi2},
            "diffraction_rotation_deg": self.diffraction_rotation_deg,
            "diffraction_rotation_scatter_deg": self.diffraction_rotation_scatter_deg,
            "beam_deviation_deg": self.beam_deviation_deg,
            "interzonal_residual_deg": self.interzonal_residual_deg,
            "zone_axis_signs": [int(sign) for sign in self.zone_axis_signs],
            "is_consistent": self.is_consistent,
            "ambiguity": self.ambiguity.to_json_dict(),
            "notes": list(self.notes),
        }


def _beam_axis_deviation_deg(matrix: np.ndarray) -> float:
    """How far a residual rotation tilts the beam axis, in degrees.

    The composition demands that the per-pattern residual be a rotation about the
    beam axis and nothing else, so this is the component of the mismatch that no
    choice of diffraction rotation can remove.
    """

    image = matrix @ BEAM_AXIS_LABORATORY
    cosine = float(np.dot(normalize_vector(image), BEAM_AXIS_LABORATORY))
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))


def _rotation_about_beam_deg(matrix: np.ndarray) -> float:
    """The beam-axis rotation angle of a residual, in degrees."""

    return float(math.degrees(math.atan2(float(matrix[1, 0]), float(matrix[0, 0]))))


def _circular_mean_and_scatter_deg(angles_deg: Sequence[float]) -> tuple[float, float]:
    """Mean and spread of angles, respecting wrap-around at 360 degrees.

    An arithmetic mean of 359 and 1 degrees is 180, which would be a spectacular
    way to report a calibration; the circular mean gives 0.
    """

    radians = np.deg2rad(np.asarray(angles_deg, dtype=np.float64))
    mean = math.degrees(
        math.atan2(float(np.mean(np.sin(radians))), float(np.mean(np.cos(radians))))
    )
    deviations = np.degrees(
        np.angle(np.exp(1j * (radians - math.radians(mean))))
    )
    scatter = float(np.sqrt(np.mean(np.square(deviations)))) if len(angles_deg) > 1 else 0.0
    return float(mean), scatter


def orientation_from_indexed_patterns(
    observations: Sequence[IndexedPatternObservation],
    stage: StageModel,
    *,
    consistency_tolerance_deg: float = BEAM_DEVIATION_WARNING_DEG,
    specimen_frame: ReferenceFrame | None = None,
    provenance: ProvenanceRecord | None = None,
) -> MultiPatternOrientation:
    """Orientation *and* diffraction rotation from two or more indexed patterns.

    Purpose
    -------
    The self-calibrating path, and the one to prefer. Two indexed patterns
    recorded at different stage positions determine the crystal orientation and
    the instrument's diffraction rotation **together**, from the data alone — so
    the answer needs no prior calibration and cannot inherit a stale one.

    When to use
    -----------
    Whenever a second pattern is available, which after any real session is
    almost always. Use it once at the start of a session to establish the
    diffraction rotation (:meth:`MultiPatternOrientation.as_calibration`), after
    which single patterns can be converted directly by
    :func:`orientation_from_indexed_pattern`.

    Algorithm
    ---------
    The zone axes alone fix the orientation: by the master equation
    ``U n_i = b_H(alpha_i, beta_i)``, and two non-parallel correspondences
    determine a rotation (Wahba's problem). A single pattern does not fix the
    *sense* of its zone axis, so all sign combinations are tried and the one
    minimizing the residual is taken.

    With ``U`` in hand, each pattern's in-plane indexing over-determines the
    diffraction rotation: ``M_i = R_stage_i U R_i^T`` must be a rotation about
    the beam axis, whose angle *is* the diffraction rotation. Two quantities come
    out of that, and both are reported: the rotation angle, whose scatter across
    patterns measures agreement; and how far ``M_i`` tilts the beam axis, which
    no choice of rotation can absorb and which therefore indicts the inputs.

    Parameters
    ----------
    observations : sequence of IndexedPatternObservation
        Two or more, for one phase, at stage positions whose beam directions are
        not parallel.
    stage : StageModel
        Its diffraction-rotation calibration is **not** used and need not exist.
    consistency_tolerance_deg : float, default 1
        Beam-axis deviation above which the fit is flagged inconsistent.
    specimen_frame : ReferenceFrame, optional
    provenance : ProvenanceRecord, optional

    Returns
    -------
    MultiPatternOrientation

    Raises
    ------
    ValueError
        Fewer than two observations, mixed phases, or parallel zone axes.

    Examples
    --------
    Two patterns down ``[001]`` and ``[011]`` of the same cubic grain, recorded
    45 degrees apart on the stage, fix the orientation exactly and return the
    instrument's diffraction rotation as a by-product.
    """

    if len(observations) < 2:
        raise ValueError(
            "Determining an orientation and a diffraction rotation together needs at "
            "least two indexed patterns at different stage positions. With one "
            "pattern, use orientation_from_indexed_pattern and supply a calibration."
        )
    phase = observations[0].zone_axis.phase
    if any(observation.zone_axis.phase != phase for observation in observations):
        raise ValueError("All indexed patterns must belong to one phase.")

    beams = np.stack(
        [
            stage.beam_direction(
                observation.position.alpha_deg, observation.position.beta_deg
            )
            for observation in observations
        ]
    )
    axes = np.stack([observation.unit_vector for observation in observations])

    # A single pattern does not fix the sense of its zone axis. Fix the first
    # arbitrarily (a global flip is a different, separately reported ambiguity)
    # and choose the rest by the interzonal angle, which is a crystallographic
    # invariant and needs no calibration.
    signs = np.ones(len(observations), dtype=np.float64)
    for index in range(1, len(observations)):
        stage_cosine = float(np.dot(beams[0], beams[index]))
        crystal_cosine = float(np.dot(axes[0], axes[index]))
        if abs(crystal_cosine - stage_cosine) > abs(-crystal_cosine - stage_cosine):
            signs[index] = -1.0
    signed_axes = axes * signs[:, None]

    if float(np.linalg.norm(np.cross(signed_axes[0], signed_axes[1]))) < 1e-8:
        raise ValueError(
            "The first two zone axes are parallel, so they cannot determine an "
            "orientation. Use patterns from zones that are well separated — a few tens "
            "of degrees is ample."
        )

    matrix = _kabsch_rotation(signed_axes, beams)

    interzonal_residual = float(
        np.max(
            [
                abs(
                    math.degrees(
                        math.acos(
                            max(-1.0, min(1.0, float(np.dot(signed_axes[0], signed_axes[i]))))
                        )
                    )
                    - math.degrees(
                        math.acos(max(-1.0, min(1.0, float(np.dot(beams[0], beams[i])))))
                    )
                )
                for i in range(1, len(observations))
            ]
        )
    )

    rotations_deg: list[float] = []
    deviations_deg: list[float] = []
    for observation, beam_index in zip(observations, range(len(observations)), strict=True):
        stage_matrix = stage.rotation_matrix(
            observation.position.alpha_deg, observation.position.beta_deg
        )
        residual = stage_matrix @ matrix @ observation.matrix.T
        rotations_deg.append(_rotation_about_beam_deg(residual))
        deviations_deg.append(_beam_axis_deviation_deg(residual))
        del beam_index

    rotation_deg, scatter_deg = _circular_mean_and_scatter_deg(rotations_deg)
    beam_deviation = float(np.max(deviations_deg))
    is_consistent = bool(
        beam_deviation <= consistency_tolerance_deg
        and interzonal_residual <= max(2.0, consistency_tolerance_deg)
    )

    target_frame = specimen_frame or HOLDER_FRAME
    orientation = Orientation.from_matrix(
        matrix,
        specimen_frame=target_frame,
        phase=phase,
        crystal_frame=phase.crystal_frame,
    )
    notes = [
        "The diffraction rotation was determined from the data, not supplied, so this "
        "orientation does not depend on any stored calibration.",
    ]
    if not is_consistent:
        notes.append(
            f"The beam-axis deviation of {beam_deviation:.2f} deg exceeds the "
            f"{consistency_tolerance_deg:.2f} deg tolerance. That component of the "
            "residual is orthogonal to the diffraction rotation, so no value of that "
            "constant would remove it."
        )

    return MultiPatternOrientation(
        orientation=orientation,
        diffraction_rotation_deg=rotation_deg,
        diffraction_rotation_scatter_deg=scatter_deg,
        beam_deviation_deg=beam_deviation,
        interzonal_residual_deg=interzonal_residual,
        zone_axis_signs=tuple(int(sign) for sign in signs),
        ambiguity=analyze_ambiguity(
            phase,
            signed_axes[0],
            rotation_calibrated=True,
            reconstruction_note=(
                f"Orientation and diffraction rotation determined jointly from "
                f"{len(observations)} indexed patterns."
            ),
        ),
        is_consistent=is_consistent,
        observations=tuple(observations),
        notes=tuple(notes),
        provenance=provenance,
    )
