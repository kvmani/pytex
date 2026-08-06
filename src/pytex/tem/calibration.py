"""Measuring the constants the instrument does not tell you.

The honest summary of TEM tilt calibration is this: **the instrument reports the
two numbers you already knew and nothing about the transformation between them
and the diffraction pattern.** Stage angles, camera length, voltage, pixel size
and magnification all arrive in the file metadata. The diffraction rotation, the
parity of the stored image, the readout sign conventions, the axis geometry and
the backlash do not, and every one of them changes the answer.

The two-excursion procedure
---------------------------
The core measurement, and the cheapest decisive test in the whole system. At a
position where a Kikuchi pattern is visible: record a reference, apply a known
positive alpha of 5-10 degrees and record how a tracked feature moved, return,
then do the same for beta.

A rigid crystal rotation carries a Kikuchi pole at the pattern centre to
``R z_lab``. At zero tilt, ``Rx(da) z = (0, -sin da, cos da)`` — a displacement
along ``-y_lab`` — and ``Ry(db) z = (sin db, 0, cos db)`` — along ``+x_lab``. So
with the pattern-to-laboratory relation ``F Rz(phi_D)``:

    phi_D = -psi_beta,     psi_alpha - psi_beta = -90 deg  (unmirrored)
                                                 +90 deg  (mirrored)

Three results from two exposures: the diffraction rotation numerically, the
parity from a sign, and a redundancy check — the two azimuths must be 90 degrees
apart, and if they are not, either the axes are non-orthogonal or the tracked
feature was misidentified, so the procedure diagnoses itself.

The displacement *magnitude* is a bonus: a pole must move by exactly the applied
angle, so the same two exposures calibrate the angular scale of the Kikuchi
pattern and check the camera length.

See ``docs/architecture/tem_tilt_navigation_foundation.md`` section 9.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from pytex.core.provenance import ProvenanceRecord
from pytex.tem.reconstruction import CurrentState, ZoneAxisObservation
from pytex.tem.stage import (
    GeneralStageAxes,
    StageCalibration,
    StageModel,
)

__all__ = [
    "CalibrationResult",
    "TiltExcursionObservation",
    "calibrate_from_tilt_excursions",
    "fit_stage_and_orientation",
    "predicted_excursion_azimuth_deg",
    "residual_from_rotation_error_deg",
]


def residual_from_rotation_error_deg(
    rotation_error_deg: float, hop_angle_deg: float
) -> float:
    """Angular miss caused by an error in the diffraction rotation.

    Purpose
    -------
    The most operationally useful formula in the package, and the reason
    multi-hop routing exists. An error ``dphi`` in the diffraction rotation is a
    rotation about the beam axis; following the tilts it implies lands the
    specimen not on the target but

        2 asin( sin(dphi/2) sin(theta) )

    away from it, where ``theta`` is the angle between the current and target
    zone axes. Two consequences:

    - The error **scales with the length of the hop**. A 5 degree calibration
      error costs 0.44 degrees over a 5 degree hop but 5 degrees over a 90 degree
      one, so several short hops with re-indexing beat one long open-loop move.
    - At ``dphi = 180 deg`` the miss is ``2 theta``: both tilt angles come out
      negated and the specimen goes exactly the wrong way, while the calculation
      still reports a clean zero residual. That is the failure this whole
      calibration exists to prevent.

    Parameters
    ----------
    rotation_error_deg : float
        The error in the assumed diffraction rotation.
    hop_angle_deg : float
        Angle between the current and target zone axes.

    Returns
    -------
    float
        Angle between the achieved direction and the beam, in degrees.

    Examples
    --------
    A five-degree calibration error over a short hop is negligible, and over a
    long one is not::

        >>> round(residual_from_rotation_error_deg(5.0, 5.0), 3)
        0.436
        >>> round(residual_from_rotation_error_deg(5.0, 90.0), 3)
        5.0
    """

    half = math.radians(rotation_error_deg) / 2.0
    argument = math.sin(half) * math.sin(math.radians(hop_angle_deg))
    return float(math.degrees(2.0 * math.asin(max(-1.0, min(1.0, argument)))))


def predicted_excursion_azimuth_deg(
    axis: str, diffraction_rotation_deg: float, *, mirrored: bool = False
) -> float:
    """Azimuth at which a Kikuchi feature moves under a positive tilt excursion.

    Purpose
    -------
    The prediction the two-excursion calibration is compared against, and the
    quantity that lets a report tell an operator *what they should see* rather
    than merely asking what they saw.

    Parameters
    ----------
    axis : {"alpha", "beta"}
    diffraction_rotation_deg : float
    mirrored : bool, default False
        Whether the stored pattern is a mirrored rendering.

    Returns
    -------
    float
        Azimuth in the stored pattern, degrees, in ``[0, 360)``.
    """

    if axis not in ("alpha", "beta"):
        raise ValueError('predicted_excursion_azimuth_deg axis must be "alpha" or "beta".')
    # Laboratory-frame motion: -y for +alpha (azimuth 270), +x for +beta (azimuth 0).
    laboratory_azimuth = 270.0 if axis == "alpha" else 0.0
    azimuth = laboratory_azimuth - diffraction_rotation_deg
    if mirrored:
        azimuth = -azimuth
    return float(azimuth % 360.0)


@dataclass(frozen=True, slots=True)
class TiltExcursionObservation:
    """A measured pattern displacement under a known single-axis tilt.

    Purpose
    -------
    The raw datum of the two-excursion calibration: *I tilted this axis by this
    much, and the tracked feature moved in this direction in the stored image.*

    Attributes
    ----------
    axis : {"alpha", "beta"}
        Which axis was moved.
    applied_deg : float
        The tilt applied. Sign matters: running the procedure with both signs is
        what determines the vendor readout convention.
    azimuth_deg : float
        Direction the tracked feature moved, measured in the **stored pattern**
        with zero along ``+x`` and increasing counter-clockwise.
    displacement_deg : float or None
        Angular magnitude of the motion, when measurable. A tracked pole must
        move by exactly ``|applied_deg|``, so this is a free check on the
        angular scale of the pattern.
    feature : str
        What was tracked, carried into the report.
    """

    axis: str
    applied_deg: float
    azimuth_deg: float
    displacement_deg: float | None = None
    feature: str = "Kikuchi pole"

    def __post_init__(self) -> None:
        if self.axis not in ("alpha", "beta"):
            raise ValueError('TiltExcursionObservation.axis must be "alpha" or "beta".')
        if self.applied_deg == 0.0:
            raise ValueError(
                "TiltExcursionObservation.applied_deg must be non-zero; a zero excursion "
                "measures nothing."
            )


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """What the calibration determined, and how well it hangs together.

    Attributes
    ----------
    calibration : StageCalibration
        The result, ready to attach to a stage model.
    orthogonality_residual_deg : float
        Departure of the measured alpha and beta azimuths from the 90 degrees
        they must differ by. A large value indicts the tracked feature or the
        axis geometry, and is the procedure's built-in self-check.
    scale_residual_deg : float or None
        Difference between the measured feature displacement and the applied
        tilt, when displacement was measured.
    is_consistent : bool
    notes : tuple of str
    """

    calibration: StageCalibration
    orthogonality_residual_deg: float
    scale_residual_deg: float | None
    is_consistent: bool
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "notes", tuple(self.notes))

    def describe(self) -> str:
        """Prose account of what was measured and whether it is trustworthy."""

        head = (
            f"Two-excursion stage calibration: diffraction rotation "
            f"{self.calibration.diffraction_rotation_deg:+.2f} deg, stored pattern "
            f"{'mirrored' if self.calibration.pattern_is_mirrored else 'not mirrored'}."
        )
        check = (
            f" Self-check: the alpha and beta motion azimuths differ from the required "
            f"90 deg by {self.orthogonality_residual_deg:.2f} deg."
        )
        if not self.is_consistent:
            check += (
                " That is large enough to indict the measurement: either the tracked "
                "feature was misidentified between exposures, or the tilt axes are not "
                "orthogonal. Re-run before relying on this."
            )
        scale = ""
        if self.scale_residual_deg is not None:
            scale = (
                f" The tracked feature moved within {self.scale_residual_deg:.2f} deg of "
                "the applied tilt, which also confirms the angular scale of the pattern "
                "and hence the camera length."
            )
        conditions = ""
        if self.calibration.camera_length_mm is not None:
            conditions = (
                f" Valid at {self.calibration.camera_length_mm:g} mm camera length"
                + (
                    f" and {self.calibration.accelerating_voltage_kv:g} kV"
                    if self.calibration.accelerating_voltage_kv is not None
                    else ""
                )
                + " only: the diffraction rotation is hysteretic in the lens settings, "
                "so this value must not be carried to another camera length."
            )
        trailing = (" " + " ".join(self.notes)) if self.notes else ""
        return head + check + scale + conditions + trailing

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "schema": "pytex.stage_calibration/1",
            "calibration": self.calibration.to_json_dict(),
            "orthogonality_residual_deg": self.orthogonality_residual_deg,
            "scale_residual_deg": self.scale_residual_deg,
            "is_consistent": self.is_consistent,
            "notes": list(self.notes),
        }


def calibrate_from_tilt_excursions(
    alpha_excursion: TiltExcursionObservation,
    beta_excursion: TiltExcursionObservation,
    *,
    base: StageCalibration | None = None,
    camera_length_mm: float | None = None,
    accelerating_voltage_kv: float | None = None,
    consistency_tolerance_deg: float = 5.0,
    provenance: ProvenanceRecord | None = None,
) -> CalibrationResult:
    """Determine the diffraction rotation and parity from two known tilts.

    Purpose
    -------
    Turns two exposures into the calibration that single-pattern navigation
    needs, and — more importantly — into the discriminating observation that
    resolves the instrumental-rotation and parity ambiguities. Roughly one minute
    at the microscope buys immunity to the failure mode that otherwise wastes a
    session.

    When to use
    -----------
    Once per camera length, at the start of a session, before any long excursion
    planned from a single indexed pattern. Not needed at all if the two-zone
    reconstruction path is used.

    Algorithm
    ---------
    A rigid crystal rotation moves a Kikuchi pole at the pattern centre along
    ``-y_lab`` for positive alpha and ``+x_lab`` for positive beta. Comparing the
    measured azimuths with those predictions gives ``phi_D = -psi_beta``
    directly, and the sign of ``psi_alpha - psi_beta`` gives the parity: ``-90``
    degrees for an unmirrored rendering, ``+90`` for a mirrored one.

    Parameters
    ----------
    alpha_excursion, beta_excursion : TiltExcursionObservation
        Must be for different axes.
    base : StageCalibration, optional
        Existing calibration to update; its axis geometry, signs and backlash are
        preserved.
    camera_length_mm, accelerating_voltage_kv : float, optional
        Recorded with the result so it cannot later be applied at conditions it
        was not measured at.
    consistency_tolerance_deg : float, default 5
        How far the alpha-versus-beta azimuth difference may stray from 90
        degrees before the result is flagged inconsistent.
    provenance : ProvenanceRecord, optional

    Returns
    -------
    CalibrationResult
        Inspect :attr:`CalibrationResult.is_consistent` before adopting it.

    Raises
    ------
    ValueError
        If both observations name the same axis, in which case the pair
        determines the rotation but not the parity.
    """

    if alpha_excursion.axis == beta_excursion.axis:
        raise ValueError(
            "The two excursions must be on different axes: one alpha and one beta. "
            "A single axis fixes the diffraction rotation but leaves the parity "
            "undetermined, which is half the point of the procedure."
        )
    if alpha_excursion.axis != "alpha":
        alpha_excursion, beta_excursion = beta_excursion, alpha_excursion

    # Normalize to what a *positive* excursion would have shown.
    alpha_azimuth = alpha_excursion.azimuth_deg + (
        180.0 if alpha_excursion.applied_deg < 0.0 else 0.0
    )
    beta_azimuth = beta_excursion.azimuth_deg + (
        180.0 if beta_excursion.applied_deg < 0.0 else 0.0
    )

    difference = (alpha_azimuth - beta_azimuth) % 360.0
    # Unmirrored: alpha azimuth is 90 deg clockwise of beta's, i.e. -90 (or 270).
    mirrored = bool(abs(_signed_difference(difference, 90.0)) < abs(
        _signed_difference(difference, 270.0)
    ))
    expected = 90.0 if mirrored else 270.0
    orthogonality_residual = abs(_signed_difference(difference, expected))

    rotation = (-beta_azimuth if not mirrored else beta_azimuth) % 360.0
    if rotation > 180.0:
        rotation -= 360.0

    scale_residual: float | None = None
    displacements = [
        (observation.displacement_deg, abs(observation.applied_deg))
        for observation in (alpha_excursion, beta_excursion)
        if observation.displacement_deg is not None
    ]
    if displacements:
        scale_residual = float(
            max(abs(measured - applied) for measured, applied in displacements)
        )

    template = base or StageCalibration()
    calibration = StageCalibration(
        axes=template.axes,
        alpha_sign=template.alpha_sign,
        beta_sign=template.beta_sign,
        alpha_zero_deg=template.alpha_zero_deg,
        beta_zero_deg=template.beta_zero_deg,
        diffraction_rotation_deg=float(rotation),
        pattern_is_mirrored=mirrored,
        camera_length_mm=(
            camera_length_mm if camera_length_mm is not None else template.camera_length_mm
        ),
        accelerating_voltage_kv=(
            accelerating_voltage_kv
            if accelerating_voltage_kv is not None
            else template.accelerating_voltage_kv
        ),
        backlash_deg=template.backlash_deg,
        angular_uncertainty_deg=template.angular_uncertainty_deg,
        notes=(
            *template.notes,
            f"Diffraction rotation measured by two tilt excursions tracking "
            f"'{alpha_excursion.feature}'.",
        ),
        provenance=provenance or template.provenance,
    )
    is_consistent = orthogonality_residual <= consistency_tolerance_deg
    notes: list[str] = []
    if not is_consistent:
        notes.append(
            "The measured azimuths are not 90 degrees apart, so this calibration should "
            "not be adopted as measured."
        )
    return CalibrationResult(
        calibration=calibration,
        orthogonality_residual_deg=float(orthogonality_residual),
        scale_residual_deg=scale_residual,
        is_consistent=is_consistent,
        notes=tuple(notes),
    )


def _signed_difference(value_deg: float, reference_deg: float) -> float:
    """Smallest signed difference between two azimuths, in ``(-180, 180]``."""

    difference = (value_deg - reference_deg) % 360.0
    return difference - 360.0 if difference > 180.0 else difference


def fit_stage_and_orientation(
    observations: tuple[ZoneAxisObservation, ...],
    stage: StageModel,
    *,
    fit_axes: bool = False,
    max_iterations: int = 60,
) -> tuple[CurrentState, StageCalibration]:
    """Jointly fit the orientation and, optionally, the stage axis geometry.

    Purpose
    -------
    With several indexed zones the problem becomes over-determined, so the axis
    geometry can in principle be measured rather than assumed. The word
    *optionally* is load-bearing: a spuriously fitted non-orthogonality is worse
    than an assumed-orthogonal stage, because it looks like knowledge and absorbs
    indexing error into a parameter that then corrupts every later solution.

    ``fit_axes`` is therefore refused below five observations, and the caller is
    expected to adopt a fitted deviation only when it is clearly significant.

    Parameters
    ----------
    observations : tuple of ZoneAxisObservation
        Three or more for the orientation; five or more to fit the axes.
    stage : StageModel
    fit_axes : bool, default False
        Fit the axis geometry as well as the orientation.
    max_iterations : int, default 60

    Returns
    -------
    tuple of (CurrentState, StageCalibration)
        The fitted orientation, and the calibration used or determined.

    Raises
    ------
    ValueError
        If ``fit_axes`` is requested with fewer than five observations.
    """

    if len(observations) < 3:
        raise ValueError("fit_stage_and_orientation needs at least three observations.")
    if fit_axes and len(observations) < 5:
        raise ValueError(
            "Fitting the stage axis geometry needs at least five indexed zone axes. "
            "Below that the problem is under-determined and the fit will absorb "
            "indexing error into an axis deviation that looks like a measurement but "
            "is not. Run with fit_axes=False, or index more zones."
        )

    if not fit_axes:
        state = CurrentState.from_zone_axes(observations, stage)
        return state, stage.calibration

    from scipy.optimize import least_squares  # local import: heavy, rarely needed

    template = stage.calibration

    def build(parameters: np.ndarray) -> StageCalibration:
        tilt_alpha, tilt_beta, coupling = (float(value) for value in parameters)
        alpha_axis = np.array([math.cos(tilt_alpha), math.sin(tilt_alpha), 0.0])
        beta_axis = np.array([math.sin(tilt_beta), math.cos(tilt_beta), 0.0])
        return StageCalibration(
            axes=GeneralStageAxes(
                alpha_axis=alpha_axis, beta_axis=beta_axis, coupling=coupling
            ),
            alpha_sign=template.alpha_sign,
            beta_sign=template.beta_sign,
            alpha_zero_deg=template.alpha_zero_deg,
            beta_zero_deg=template.beta_zero_deg,
            diffraction_rotation_deg=template.diffraction_rotation_deg,
            pattern_is_mirrored=template.pattern_is_mirrored,
            camera_length_mm=template.camera_length_mm,
            accelerating_voltage_kv=template.accelerating_voltage_kv,
            backlash_deg=template.backlash_deg,
            angular_uncertainty_deg=template.angular_uncertainty_deg,
            notes=template.notes,
        )

    def residuals(parameters: np.ndarray) -> np.ndarray:
        candidate = _with_calibration(stage, build(parameters))
        state = CurrentState.from_zone_axes(observations, candidate)
        matrix = state.matrix
        errors = []
        for observation in observations:
            beam = candidate.beam_direction(
                observation.position.alpha_deg, observation.position.beta_deg
            )
            achieved = matrix @ observation.unit_vector
            errors.append(
                math.degrees(
                    math.acos(
                        max(-1.0, min(1.0, abs(float(np.dot(achieved, beam)))))
                    )
                )
            )
        return np.asarray(errors, dtype=np.float64)

    solution = least_squares(
        residuals,
        x0=np.zeros(3),
        max_nfev=max_iterations,
        method="lm",
    )
    calibration = build(solution.x)
    fitted_stage = _with_calibration(stage, calibration)
    state = CurrentState.from_zone_axes(observations, fitted_stage)
    return state, calibration


def _with_calibration(stage: StageModel, calibration: StageCalibration) -> StageModel:
    """A copy of ``stage`` carrying a different calibration."""

    from dataclasses import replace

    return replace(stage, calibration=calibration)  # type: ignore[type-var]
