"""Getting there: the route from the current stage position to the target one.

The endpoint is not the answer. A straight line in ``(alpha, beta)`` space is not
a straight line on the crystal sphere — the beta contribution is scaled by
``cos(alpha)`` and the two rotations do not commute — so different interpolations
trace genuinely different trajectories through orientation space, and they differ
in whether the operator can follow along.

The geodesic path is the Kikuchi band
-------------------------------------
The shortest arc between two zone axes on the crystal sphere lies in the plane
they span. That plane's normal rationalizes to a low-index ``(hkl)``, and the
Kikuchi band of that reflection is exactly the band joining the two poles. So the
mathematically optimal path is the one experienced operators already follow by
eye, and the planner can tell the user which band to follow:

    "Follow the (1-10) Kikuchi band from [001] toward [111]; total travel 54.7 deg."

That is why :data:`PathStrategy.GEODESIC` is the default. The sequential
strategies remain available because a stage that cannot drive both axes at once
will execute one leg at a time whatever the plan says, and a plan that does not
match what the operator will actually do is worse than a longer one.

Multi-hop routing
-----------------
An error ``dphi`` in the diffraction rotation costs a residual of about
``dphi sin(theta)`` where ``theta`` is the angular length of the hop. Long
excursions are therefore fragile, and the fix is structural rather than
numerical: route through intermediate low-index zone axes and re-solve the
orientation at each. Each hop is short, so each is robust, and re-indexing at
every waypoint converts an open-loop calculation into a self-correcting
procedure.

See ``docs/architecture/tem_tilt_navigation_foundation.md`` section 11.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import as_float_array, normalize_vector
from pytex.core.lattice import CrystalDirection, MillerIndex, Phase
from pytex.core.notation import format_plane_indices
from pytex.tem.stage import GIMBAL_TOLERANCE, StageModel, StagePosition

__all__ = [
    "DEFAULT_MARGIN_FLOOR_DEG",
    "DEFAULT_PATH_SAMPLES",
    "DEFAULT_RECENTRE_INTERVAL_DEG",
    "PathStrategy",
    "TiltPath",
    "TiltPathSample",
    "connecting_band",
    "plan_path",
    "suggest_waypoints",
]

#: Number of samples along a planned path, including both endpoints.
#:
#: Chosen so that the stereographic trajectory reads as a dotted arc rather than
#: a dashed line at publication size, and so that envelope violations narrower
#: than about a degree cannot slip between samples on a typical excursion.
DEFAULT_PATH_SAMPLES = 25

#: Minimum envelope margin a path sample must retain, in degrees.
#:
#: A path that grazes a mechanical limit will hit it once backlash is included,
#: so the planner keeps a working clearance rather than merely staying inside.
DEFAULT_MARGIN_FLOOR_DEG = 2.0

#: Travel between prompts to re-centre the area of interest, in degrees.
#:
#: Being off eucentric height introduces no orientation error, but it does
#: translate the region out of the selected-area aperture — so the operator loses
#: the crystal rather than mis-orienting it. Segmenting the path is the remedy.
DEFAULT_RECENTRE_INTERVAL_DEG = 15.0


def _int_list(values: ArrayLike) -> list[int]:
    """Round an index array to the plain integer list the notation helpers take."""

    return [round(float(value)) for value in np.asarray(values).ravel()]


class PathStrategy(StrEnum):
    """How to interpolate between the current and target stage positions."""

    #: Shortest arc of the beam direction on the crystal sphere, inverted to
    #: stage angles at each sample. Minimum crystal travel, and the path that
    #: follows the connecting Kikuchi band. The default.
    GEODESIC = "geodesic"
    #: Straight line in ``(alpha, beta)``. Simple; no special virtue.
    LINEAR = "linear"
    #: Beta to completion, then alpha. Easiest to execute by hand.
    BETA_THEN_ALPHA = "beta_then_alpha"
    #: Alpha to completion, then beta.
    ALPHA_THEN_BETA = "alpha_then_beta"


@dataclass(frozen=True, slots=True)
class TiltPathSample:
    """One point along a planned tilt path, with everything a consumer needs.

    Purpose
    -------
    Each sample carries the stage angles, the resulting beam direction in crystal
    coordinates, the nearest low-index zone axis and the angle to it, the
    envelope margin, and the cumulative travel. Different consumers need
    different fields — the operator reads the table, the plotter draws the
    dots, the validator checks the margins — and storing all of them prevents
    any consumer from re-deriving a quantity and drifting from the engine.

    Attributes
    ----------
    fraction : float
        Position along the path, 0 at the start and 1 at the target.
    position : StagePosition
    beam_direction_crystal : np.ndarray
        Unit vector: what the operator is looking down, crystallographically.
    envelope_margin_deg : float
    cumulative_travel_deg : float
        Arc length swept by the beam direction in the crystal frame, which is
        the angle a Kikuchi pattern actually moves through.
    """

    fraction: float
    position: StagePosition
    beam_direction_crystal: np.ndarray
    envelope_margin_deg: float
    cumulative_travel_deg: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "beam_direction_crystal",
            as_float_array(self.beam_direction_crystal, shape=(3,)),
        )


@dataclass(frozen=True, slots=True)
class TiltPath:
    """A validated route from the current stage position to a target position.

    Purpose
    -------
    The plan the operator works from and the array the stereographic figure
    draws. It is deliberately the *same* array: the figure plots what the engine
    produced rather than recomputing the trajectory, so a drawing that agrees
    with the engine is evidence about the engine and not about a second copy of
    the same formula.

    Attributes
    ----------
    samples : tuple of TiltPathSample
    strategy : PathStrategy
    is_valid : bool
        Whether every sample stayed inside the envelope with the required margin.
    violation_reason : str
        Why not, when ``is_valid`` is false.
    total_travel_deg : float
        Arc swept by the beam direction in the crystal frame.
    stage_rotation_deg : float
        Angle of the actual specimen rotation between the endpoints — the
        rotation angle of ``R_target R_current^T``. This is a genuine angle of a
        genuine rotation, unlike the naive ``hypot(d_alpha, d_beta)``.
    connecting_band : MillerIndex or None
        The Kikuchi band joining the start and end zone axes, when it
        rationalizes to a sensible low-index plane.
    recentre_fractions : tuple of float
        Path fractions at which the operator should re-centre the area of
        interest.
    approach_note : str
        The backlash-aware approach instruction, when one applies.
    """

    samples: tuple[TiltPathSample, ...]
    strategy: PathStrategy
    is_valid: bool
    violation_reason: str
    total_travel_deg: float
    stage_rotation_deg: float
    connecting_band: MillerIndex | None = None
    recentre_fractions: tuple[float, ...] = ()
    approach_note: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "samples", tuple(self.samples))
        object.__setattr__(self, "recentre_fractions", tuple(self.recentre_fractions))
        if not self.samples:
            raise ValueError("A TiltPath must contain at least one sample.")

    @property
    def start(self) -> StagePosition:
        """The stage position the path begins at."""

        return self.samples[0].position

    @property
    def end(self) -> StagePosition:
        """The stage position the path ends at."""

        return self.samples[-1].position

    @property
    def minimum_margin_deg(self) -> float:
        """The tightest envelope clearance anywhere along the path."""

        return float(min(sample.envelope_margin_deg for sample in self.samples))

    def beam_directions_crystal(self) -> np.ndarray:
        """All sampled beam directions as an ``(n, 3)`` array.

        The array the stereographic projection plots as its trajectory dots.
        """

        return np.stack([sample.beam_direction_crystal for sample in self.samples])

    def stage_positions(self) -> np.ndarray:
        """All sampled stage angles as an ``(n, 2)`` array in degrees."""

        return np.asarray(
            [[s.position.alpha_deg, s.position.beta_deg] for s in self.samples],
            dtype=np.float64,
        )

    def describe(self) -> str:
        """Prose route description an operator can follow at the microscope."""

        head = (
            f"Path ({self.strategy.value}) from {self.start} to {self.end}: "
            f"{self.total_travel_deg:.1f} deg of crystal travel, a "
            f"{self.stage_rotation_deg:.1f} deg specimen rotation, in "
            f"{len(self.samples)} steps."
        )
        band = ""
        if self.connecting_band is not None:
            band_text = format_plane_indices(
                _int_list(self.connecting_band.indices), style="plain"
            )
            band = (
                f" Follow the {band_text} "
                "Kikuchi band, which is the band joining the two zone axes: the "
                "geodesic path and the band are the same great circle."
            )
        margin = f" Tightest envelope clearance along the way: {self.minimum_margin_deg:.1f} deg."
        recentre = ""
        if self.recentre_fractions:
            recentre = (
                f" Re-centre the area of interest at "
                f"{len(self.recentre_fractions)} point(s) along the path; off-eucentric "
                "tilting does not mis-orient the crystal but it does translate the "
                "region out of the aperture."
            )
        approach = f" {self.approach_note}" if self.approach_note else ""
        validity = (
            ""
            if self.is_valid
            else f" PATH REJECTED: {self.violation_reason}"
        )
        return head + band + margin + recentre + approach + validity

    def to_json_dict(self) -> dict[str, Any]:
        """Serializable payload, kept in lockstep with :meth:`describe`."""

        return {
            "strategy": self.strategy.value,
            "is_valid": self.is_valid,
            "violation_reason": self.violation_reason,
            "total_travel_deg": self.total_travel_deg,
            "stage_rotation_deg": self.stage_rotation_deg,
            "minimum_margin_deg": self.minimum_margin_deg,
            "connecting_band": (
                None
                if self.connecting_band is None
                else [int(value) for value in self.connecting_band.indices]
            ),
            "recentre_fractions": list(self.recentre_fractions),
            "approach_note": self.approach_note,
            "samples": [
                {
                    "fraction": sample.fraction,
                    "alpha_deg": sample.position.alpha_deg,
                    "beta_deg": sample.position.beta_deg,
                    "beam_direction_crystal": sample.beam_direction_crystal.tolist(),
                    "envelope_margin_deg": sample.envelope_margin_deg,
                    "cumulative_travel_deg": sample.cumulative_travel_deg,
                }
                for sample in self.samples
            ],
        }


def _stage_angles_for_holder_direction(direction: ArrayLike) -> tuple[float, float]:
    """Invert the ideal beam-direction formula: holder direction to stage angles.

    The closed form of equation (S). Used to sample a geodesic, which is defined
    in direction space and must be pulled back to stage angles at every point.
    """

    w = normalize_vector(direction)
    rho = math.hypot(float(w[0]), float(w[2]))
    if rho < GIMBAL_TOLERANCE:
        # Along the beta axis: beta is indeterminate, alpha is +/-90 degrees.
        return (90.0 if w[1] > 0.0 else -90.0, 0.0)
    return (
        math.degrees(math.atan2(float(w[1]), rho)),
        math.degrees(math.atan2(-float(w[0]), float(w[2]))),
    )


def _slerp(start: np.ndarray, end: np.ndarray, fractions: np.ndarray) -> np.ndarray:
    """Great-circle interpolation between two unit vectors.

    Falls back to linear interpolation with renormalization when the endpoints
    are nearly parallel, where the spherical formula is ill-conditioned and the
    two agree to machine precision anyway.
    """

    dot = float(np.clip(np.dot(start, end), -1.0, 1.0))
    omega = math.acos(dot)
    if omega < 1e-9:
        return np.tile(start, (fractions.size, 1))
    sin_omega = math.sin(omega)
    a = np.sin((1.0 - fractions) * omega) / sin_omega
    b = np.sin(fractions * omega) / sin_omega
    return np.asarray(
        a[:, None] * start[None, :] + b[:, None] * end[None, :], dtype=np.float64
    )


def connecting_band(
    phase: Phase,
    first_direction: ArrayLike,
    second_direction: ArrayLike,
    *,
    max_index: int = 6,
    tolerance_deg: float = 3.0,
) -> MillerIndex | None:
    """The Kikuchi band joining two zone axes, as a low-index reflection.

    Purpose
    -------
    Turns the geodesic into an instruction an operator can act on. Two zone axes
    span a plane; the reflection whose normal is that plane's normal is the
    Kikuchi band running between the two poles, and following it is exactly what
    a skilled operator does by eye. Reporting it makes the calculated path
    legible at the microscope instead of being a pair of numbers.

    Parameters
    ----------
    phase : Phase
    first_direction, second_direction : array_like
        Cartesian crystal-frame directions; need not be normalized.
    max_index : int, default 6
        Largest absolute index considered when rationalizing the normal.
    tolerance_deg : float, default 3
        Largest deviation accepted between the true normal and the rationalized
        one. Beyond this the band would be mislabelled, and ``None`` — no
        recommendation — is the honest answer.

    Returns
    -------
    MillerIndex or None
        ``None`` when the directions are parallel, so no plane is defined, or
        when no sufficiently low-index plane matches.

    Examples
    --------
    In a cubic crystal the band joining ``[001]`` and ``[111]`` is ``(1-10)``:
    the two poles both lie in that plane, which is why tilting between them
    tracks a single band.
    """

    first = np.asarray(first_direction, dtype=np.float64)
    second = np.asarray(second_direction, dtype=np.float64)
    normal = np.cross(first, second)
    if float(np.linalg.norm(normal)) < 1e-9:
        return None
    normal = np.asarray(normal / float(np.linalg.norm(normal)), dtype=np.float64)

    reciprocal = phase.lattice.reciprocal_basis().matrix
    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1)
    hkl = grid.reshape(-1, 3)
    hkl = hkl[np.any(hkl != 0, axis=1)]
    cartesian = hkl.astype(np.float64) @ reciprocal.T
    norms = np.linalg.norm(cartesian, axis=1)
    cosines = np.abs(cartesian @ normal) / norms
    deviations = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
    admissible = deviations <= tolerance_deg
    if not np.any(admissible):
        return None
    candidates = hkl[admissible]
    # Prefer the lowest-index member of the admissible set: a band is named by
    # its simplest indices, and the higher orders are the same physical band.
    order = np.lexsort(
        (
            np.max(np.abs(candidates), axis=1),
            np.sum(np.abs(candidates), axis=1),
        )
    )
    best = candidates[order[0]]
    if int(np.sum(best)) < 0 or (int(np.sum(best)) == 0 and int(best[0]) < 0):
        best = -best
    return MillerIndex(best, phase=phase)


def plan_path(
    start: StagePosition,
    target: StagePosition,
    stage: StageModel,
    crystal_to_holder: ArrayLike,
    *,
    strategy: PathStrategy = PathStrategy.GEODESIC,
    samples: int = DEFAULT_PATH_SAMPLES,
    margin_floor_deg: float = DEFAULT_MARGIN_FLOOR_DEG,
    recentre_interval_deg: float = DEFAULT_RECENTRE_INTERVAL_DEG,
    phase: Phase | None = None,
    backlash_deg: float | None = None,
) -> TiltPath:
    """Plan and validate a route between two stage positions.

    Purpose
    -------
    Produces the trajectory the operator follows and the figure draws, with every
    sample checked against the holder envelope. Validation is part of planning,
    not a later step: a path that leaves the envelope halfway is not a path, and
    returning it with a warning would invite it to be used.

    When to use
    -----------
    Called automatically by :func:`pytex.tem.navigation.plan_tilt_to_zone_axis`
    for every candidate solution. Call it directly to compare strategies, or to
    replan after an operator has moved partway.

    Parameters
    ----------
    start, target : StagePosition
    stage : StageModel
    crystal_to_holder : array_like
        The 3x3 orientation ``U``. Needed to express each sample's beam direction
        in crystal coordinates, which is what makes the path interpretable.
    strategy : PathStrategy, default GEODESIC
        The geodesic minimizes crystal travel and coincides with the connecting
        Kikuchi band.
    samples : int, default 25
    margin_floor_deg : float, default 2
        Required clearance from the envelope boundary along the whole path.
    recentre_interval_deg : float, default 15
        Travel between re-centring prompts.
    phase : Phase, optional
        Enables the connecting-band recommendation.
    backlash_deg : float, optional
        Overrides the stage calibration's measured backlash for the approach
        note.

    Returns
    -------
    TiltPath
        Check :attr:`TiltPath.is_valid` before acting on it.

    Notes
    -----
    The geodesic strategy interpolates in *direction* space and inverts to stage
    angles at each sample, so the trajectory is the true shortest arc rather than
    an approximation of it in angle space.
    """

    if samples < 2:
        raise ValueError("A path needs at least two samples.")
    matrix = as_float_array(crystal_to_holder, shape=(3, 3))
    fractions = np.linspace(0.0, 1.0, samples)

    start_beam = stage.beam_direction(start.alpha_deg, start.beta_deg)
    target_beam = stage.beam_direction(target.alpha_deg, target.beta_deg)

    if strategy is PathStrategy.GEODESIC:
        directions = _slerp(start_beam, target_beam, fractions)
        angle_pairs = [_stage_angles_for_holder_direction(row) for row in directions]
    elif strategy is PathStrategy.LINEAR:
        angle_pairs = [
            (
                start.alpha_deg + f * (target.alpha_deg - start.alpha_deg),
                start.beta_deg + f * (target.beta_deg - start.beta_deg),
            )
            for f in fractions
        ]
    elif strategy is PathStrategy.BETA_THEN_ALPHA:
        angle_pairs = _two_leg_angles(start, target, fractions, beta_first=True)
    else:
        angle_pairs = _two_leg_angles(start, target, fractions, beta_first=False)

    path_samples: list[TiltPathSample] = []
    cumulative = 0.0
    previous_crystal: np.ndarray | None = None
    for fraction, (alpha, beta) in zip(fractions, angle_pairs, strict=True):
        beam_holder = stage.beam_direction(alpha, beta)
        beam_crystal = matrix.T @ beam_holder
        if previous_crystal is not None:
            cumulative += math.degrees(
                math.acos(
                    max(-1.0, min(1.0, float(np.dot(previous_crystal, beam_crystal))))
                )
            )
        previous_crystal = beam_crystal
        path_samples.append(
            TiltPathSample(
                fraction=float(fraction),
                position=StagePosition(float(alpha), float(beta)),
                beam_direction_crystal=beam_crystal,
                envelope_margin_deg=float(stage.envelope.margin_deg(alpha, beta)),
                cumulative_travel_deg=float(cumulative),
            )
        )

    violations = [
        sample for sample in path_samples if sample.envelope_margin_deg < margin_floor_deg
    ]
    outside = [sample for sample in path_samples if sample.envelope_margin_deg < 0.0]
    if outside:
        worst = min(outside, key=lambda s: s.envelope_margin_deg)
        reason = (
            f"the path leaves the holder envelope at {worst.position} "
            f"(margin {worst.envelope_margin_deg:+.1f} deg). Try a different "
            "symmetry-equivalent target, or a sequential strategy whose legs stay "
            "inside."
        )
    elif violations:
        worst = min(violations, key=lambda s: s.envelope_margin_deg)
        reason = (
            f"the path approaches within {worst.envelope_margin_deg:.1f} deg of a "
            f"mechanical limit at {worst.position}, below the {margin_floor_deg:.1f} deg "
            "working clearance. Backlash alone can carry the stage into the stop from "
            "there."
        )
    else:
        reason = ""

    stage_rotation = float(
        _rotation_angle_between(
            stage.rotation_matrix(start.alpha_deg, start.beta_deg),
            stage.rotation_matrix(target.alpha_deg, target.beta_deg),
        )
    )
    band = (
        connecting_band(
            phase,
            matrix.T @ start_beam,
            matrix.T @ target_beam,
        )
        if phase is not None
        else None
    )
    total_travel = path_samples[-1].cumulative_travel_deg
    recentre = tuple(
        float(f)
        for f in np.arange(recentre_interval_deg, total_travel, recentre_interval_deg)
        / max(total_travel, 1e-12)
    )
    backlash = (
        stage.calibration.backlash_deg if backlash_deg is None else float(backlash_deg)
    )
    approach = ""
    if backlash > 0.0:
        overshoot = max(2.0 * backlash, 1.0)
        approach = (
            f"Approach the final position from the negative side on both axes: "
            f"overshoot by {overshoot:.1f} deg and return. Measured backlash is "
            f"{backlash:.2f} deg, which is not correctable in open loop, only avoided "
            "by a consistent approach direction."
        )

    return TiltPath(
        samples=tuple(path_samples),
        strategy=strategy,
        is_valid=not violations,
        violation_reason=reason,
        total_travel_deg=float(total_travel),
        stage_rotation_deg=stage_rotation,
        connecting_band=band,
        recentre_fractions=recentre,
        approach_note=approach,
    )


def _two_leg_angles(
    start: StagePosition,
    target: StagePosition,
    fractions: np.ndarray,
    *,
    beta_first: bool,
) -> list[tuple[float, float]]:
    """Sample a two-leg sequential path, first axis complete before the second."""

    pairs: list[tuple[float, float]] = []
    for fraction in fractions:
        if fraction <= 0.5:
            leg = fraction * 2.0
            if beta_first:
                pairs.append(
                    (start.alpha_deg, start.beta_deg + leg * (target.beta_deg - start.beta_deg))
                )
            else:
                pairs.append(
                    (
                        start.alpha_deg + leg * (target.alpha_deg - start.alpha_deg),
                        start.beta_deg,
                    )
                )
        else:
            leg = (fraction - 0.5) * 2.0
            if beta_first:
                pairs.append(
                    (
                        start.alpha_deg + leg * (target.alpha_deg - start.alpha_deg),
                        target.beta_deg,
                    )
                )
            else:
                pairs.append(
                    (
                        target.alpha_deg,
                        start.beta_deg + leg * (target.beta_deg - start.beta_deg),
                    )
                )
    return pairs


def _rotation_angle_between(first: np.ndarray, second: np.ndarray) -> float:
    """Angle of the rotation carrying ``first`` onto ``second``, in degrees."""

    relative = second @ first.T
    trace = float(np.trace(relative))
    return float(math.degrees(math.acos(max(-1.0, min(1.0, (trace - 1.0) / 2.0)))))


def suggest_waypoints(
    phase: Phase,
    start_direction: ArrayLike,
    target_direction: ArrayLike,
    *,
    max_index: int = 3,
    max_waypoints: int = 3,
    corridor_deg: float = 12.0,
    minimum_separation_deg: float = 12.0,
) -> tuple[CrystalDirection, ...]:
    """Low-index zone axes lying near the geodesic, for multi-hop routing.

    Purpose
    -------
    Long excursions are fragile: an error ``dphi`` in the diffraction rotation
    costs about ``dphi sin(theta)`` for a hop of angular length ``theta``. Routing
    through intermediate zones and re-indexing at each keeps every hop short, so
    accumulated calibration error is reset rather than compounded. That converts
    an open-loop calculation into a closed-loop procedure.

    Waypoints are scored by how close they lie to the geodesic and how low their
    indices are — a low-index zone is easier to recognize and to index, which is
    the whole point of stopping there.

    Parameters
    ----------
    phase : Phase
    start_direction, target_direction : array_like
        Cartesian crystal-frame directions.
    max_index : int, default 3
        Largest index considered. Higher-index zones make poor waypoints because
        they are hard to recognize.
    max_waypoints : int, default 3
    corridor_deg : float, default 12
        How far off the geodesic a zone may lie and still be worth a stop.
    minimum_separation_deg : float, default 12
        Waypoints closer than this to an endpoint or to each other are dropped;
        a stop that saves two degrees is not worth an exposure.

    Returns
    -------
    tuple of CrystalDirection
        Ordered from the start toward the target. Empty when the hop is short
        enough not to need them.
    """

    start = normalize_vector(start_direction)
    target = normalize_vector(target_direction)
    normal = np.cross(start, target)
    if float(np.linalg.norm(normal)) < 1e-9:
        return ()
    normal = normal / float(np.linalg.norm(normal))
    span = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(start, target))))))

    direct = phase.lattice.direct_basis().matrix
    values = np.arange(-max_index, max_index + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(values, values, values, indexing="ij"), axis=-1)
    uvw = grid.reshape(-1, 3)
    uvw = uvw[np.any(uvw != 0, axis=1)]
    cartesian = uvw.astype(np.float64) @ direct.T
    cartesian = cartesian / np.linalg.norm(cartesian, axis=1)[:, None]

    # Distance from the great circle, and position along it.
    off_plane_deg = np.degrees(np.arcsin(np.clip(np.abs(cartesian @ normal), -1.0, 1.0)))
    along_deg = np.degrees(
        np.arccos(np.clip(cartesian @ start, -1.0, 1.0))
    )
    toward_target = cartesian @ target
    admissible = (
        (off_plane_deg <= corridor_deg)
        & (along_deg > minimum_separation_deg)
        & (along_deg < span - minimum_separation_deg)
        & (toward_target > 0.0)
    )
    if not np.any(admissible):
        return ()

    candidates = uvw[admissible]
    scores = off_plane_deg[admissible] + 2.0 * np.max(np.abs(candidates), axis=1)
    order = np.argsort(scores)

    chosen: list[np.ndarray] = []
    chosen_along: list[float] = []
    admissible_along = along_deg[admissible]
    for index in order:
        position = float(admissible_along[index])
        if any(abs(position - other) < minimum_separation_deg for other in chosen_along):
            continue
        chosen.append(candidates[index])
        chosen_along.append(position)
        if len(chosen) >= max_waypoints:
            break

    order_along = sorted(range(len(chosen)), key=lambda i: chosen_along[i])
    return tuple(
        CrystalDirection(chosen[i].astype(np.float64), phase=phase) for i in order_along
    )
