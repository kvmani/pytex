"""Scoring an indexed SAED solution: what disagrees, by how much, and how much it matters.

Why a score at all
------------------
:class:`pytex.diffraction.solving.PatternSolution` already ranks candidates by a
sort key — matched fraction, then mean residual — and that key is deliberately
*not* offered as a quality: it orders solutions and says nothing about whether
the best one is any good. A user comparing two candidates needs more than an
order. They need to know **what** disagrees (a length? an angle?), **by how
much**, in units they measure in, and **how the disagreements combine**.

This module supplies that, and keeps the two halves apart:

**Deviations are measurements.** For every indexed spot, the measured d-spacing
against the calculated one; for every pair of indexed spots, the measured angle
between them against the calculated one. Both are reported per spot and per pair,
with no weighting and no opinion. They are the evidence.

**The score is a policy.** Combining a 1.2% length error with a 0.4 degree angle
error into one number requires saying how much each matters, and there is no
universal answer — it depends on the instrument, the calibration, and what the
answer is for. So the policy lives in :class:`ScoringWeights`, every term is
documented with the reasoning behind its default, and the report carries the
weights that produced it. A score whose policy is not visible is not a
measurement, it is an assertion.

What the terms mean
-------------------
Three things distinguish a right answer from a plausible one, and they fail
independently:

``length``
    Relative d-spacing agreement. Sensitive to the camera constant, which is the
    one calibration that can be wrong while everything stays self-consistent.
``angle``
    Interspot angle agreement. **Calibration-free** — angles do not depend on the
    camera constant at all — which makes this the term that catches a wrong
    *phase* rather than a wrong scale.
``coverage``
    The fraction of picked spots the solution indexes. A solution explaining half
    the spots perfectly is usually a coincidence on a sub-lattice.

Each is mapped to ``[0, 1]`` by a soft, scale-set-by-tolerance curve rather than
a hard pass/fail, because a solution 1.01 times over tolerance is not
categorically worse than one at 0.99, and a cliff there would make the ranking
jump around under a change nobody can see.

Limits
------
Kinematic geometry only. Intensities are not scored: relative intensities in a
real pattern are dynamical, vary with thickness and tilt, and a solution scored
on them would be scored on the specimen rather than the crystallography. And the
score cannot resolve what one pattern cannot: a zone axis and its reverse give
identical spot positions, so both score identically and always will.

See :doc:`/workflows/tem_pattern_indexing` for the workflow and
:doc:`/theory/saed_ratio_angle_indexing` for the indexing this scores.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

__all__ = [
    "AngleDeviation",
    "ScoringWeights",
    "SolutionScore",
    "SpotDeviation",
    "score_solution",
]


@dataclass(frozen=True, slots=True)
class ScoringWeights:
    """How much each kind of disagreement counts, and what counts as agreement.

    Purpose
    -------
    Makes the scoring policy an explicit, inspectable, overridable object rather
    than a set of constants buried in a function. Two laboratories can disagree
    about the weights and still compare deviations, because the deviations do not
    depend on them.

    Attributes
    ----------
    length : float
        Weight on relative d-spacing agreement. Defaults to 1.0, the reference
        against which the others are set.
    angle : float
        Weight on interspot-angle agreement. Defaults to 1.5, **above** length,
        because angles are calibration-free: a wrong camera constant scales every
        length and leaves every angle alone, so an angle disagreement is evidence
        about the crystallography rather than about the instrument.
    coverage : float
        Weight on the fraction of picked spots indexed. Defaults to 2.0, the
        largest, because an unindexed spot is unexplained evidence and no amount
        of precision on the remaining spots makes it go away.
    length_tolerance : float
        Relative d-spacing deviation scoring 0.5 on the length term. Defaults to
        0.02 — two percent, about what a well-calibrated instrument achieves.
    angle_tolerance_deg : float
        Angular deviation scoring 0.5 on the angle term. Defaults to 1.0 degree,
        roughly the precision of picking two spots and measuring between them.
    sharpness : float
        Exponent of the soft scoring curve, ``1 / (1 + (x/t)^s)``. Defaults to
        2.0. Raising it makes the score behave more like a pass/fail at the
        tolerance; lowering it makes it more forgiving of a single bad spot.

    Examples
    --------
    A house policy that trusts its calibration and cares mostly about geometry::

        >>> weights = ScoringWeights(length=0.5, angle=2.0, length_tolerance=0.05)
        >>> weights.total
        4.5
    """

    length: float = 1.0
    angle: float = 1.5
    coverage: float = 2.0
    length_tolerance: float = 0.02
    angle_tolerance_deg: float = 1.0
    sharpness: float = 2.0

    def __post_init__(self) -> None:
        for name in ("length", "angle", "coverage"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"ScoringWeights.{name} must be finite and non-negative.")
        if self.length + self.angle + self.coverage <= 0.0:
            raise ValueError("At least one ScoringWeights term must carry weight.")
        if not self.length_tolerance > 0.0:
            raise ValueError("ScoringWeights.length_tolerance must be strictly positive.")
        if not self.angle_tolerance_deg > 0.0:
            raise ValueError("ScoringWeights.angle_tolerance_deg must be strictly positive.")
        if not self.sharpness > 0.0:
            raise ValueError("ScoringWeights.sharpness must be strictly positive.")

    @property
    def total(self) -> float:
        """Sum of the three term weights, the normalizer of the fused score."""

        return float(self.length + self.angle + self.coverage)

    def agreement(self, deviation: float, tolerance: float) -> float:
        """Map a deviation to ``[0, 1]``, scoring ``0.5`` at the tolerance.

        ``1 / (1 + (x/t)^s)``. Chosen over a hard threshold because a threshold
        makes the ranking discontinuous exactly where candidates are hardest to
        tell apart, and over a Gaussian because this decays polynomially, so a
        badly wrong solution keeps a small non-zero score and stays comparable
        with another badly wrong one instead of both underflowing to zero.
        """

        if not math.isfinite(deviation):
            return 0.0
        ratio = abs(float(deviation)) / float(tolerance)
        return float(1.0 / (1.0 + ratio ** float(self.sharpness)))

    def to_json(self) -> dict[str, float]:
        """The policy as JSON-ready data, so a score travels with its weights."""

        return {
            "length": float(self.length),
            "angle": float(self.angle),
            "coverage": float(self.coverage),
            "length_tolerance": float(self.length_tolerance),
            "angle_tolerance_deg": float(self.angle_tolerance_deg),
            "sharpness": float(self.sharpness),
        }

    def replace(self, **changes: float) -> ScoringWeights:
        """A copy with some terms changed, validated like any other instance."""

        return replace(self, **changes)


@dataclass(frozen=True, slots=True)
class SpotDeviation:
    """Measured against calculated, for one indexed spot."""

    spot: int
    hkl: tuple[int, int, int]
    label: str
    d_measured_angstrom: float
    d_calculated_angstrom: float

    @property
    def absolute_deviation_angstrom(self) -> float:
        """Signed difference, measured minus calculated."""

        return float(self.d_measured_angstrom - self.d_calculated_angstrom)

    @property
    def relative_deviation(self) -> float:
        """Signed fractional difference. The quantity a camera constant biases."""

        if self.d_calculated_angstrom == 0.0:
            return float("inf")
        return float(self.absolute_deviation_angstrom / self.d_calculated_angstrom)

    def to_json(self) -> dict[str, Any]:
        """The row as JSON-ready data."""

        return {
            "spot": int(self.spot),
            "hkl": [int(value) for value in self.hkl],
            "label": self.label,
            "d_measured": float(self.d_measured_angstrom),
            "d_calculated": float(self.d_calculated_angstrom),
            "deviation_angstrom": float(self.absolute_deviation_angstrom),
            "relative_deviation": float(self.relative_deviation),
        }


@dataclass(frozen=True, slots=True)
class AngleDeviation:
    """Measured against calculated, for the angle between two indexed spots."""

    first_spot: int
    second_spot: int
    first_label: str
    second_label: str
    measured_deg: float
    calculated_deg: float

    @property
    def deviation_deg(self) -> float:
        """Signed difference, measured minus calculated."""

        return float(self.measured_deg - self.calculated_deg)

    def to_json(self) -> dict[str, Any]:
        """The row as JSON-ready data."""

        return {
            "first_spot": int(self.first_spot),
            "second_spot": int(self.second_spot),
            "pair": f"{self.first_label} to {self.second_label}",
            "measured_deg": float(self.measured_deg),
            "calculated_deg": float(self.calculated_deg),
            "deviation_deg": float(self.deviation_deg),
        }


@dataclass(frozen=True, slots=True)
class SolutionScore:
    """Deviations, the terms they produce, and the one number they fuse to.

    Attributes
    ----------
    spot_deviations : tuple of SpotDeviation
    angle_deviations : tuple of AngleDeviation
    matched_fraction : float
    weights : ScoringWeights
        The policy that produced ``score``, carried so a number can be traced to
        the opinion behind it.
    """

    spot_deviations: tuple[SpotDeviation, ...]
    angle_deviations: tuple[AngleDeviation, ...]
    matched_fraction: float
    weights: ScoringWeights

    def __post_init__(self) -> None:
        object.__setattr__(self, "spot_deviations", tuple(self.spot_deviations))
        object.__setattr__(self, "angle_deviations", tuple(self.angle_deviations))
        if not 0.0 <= self.matched_fraction <= 1.0:
            raise ValueError("SolutionScore.matched_fraction must lie in [0, 1].")

    @property
    def rms_relative_length_deviation(self) -> float:
        """Root-mean-square fractional d-spacing disagreement."""

        if not self.spot_deviations:
            return float("inf")
        values = np.asarray(
            [spot.relative_deviation for spot in self.spot_deviations], dtype=float
        )
        return float(np.sqrt(np.mean(values**2)))

    @property
    def max_relative_length_deviation(self) -> float:
        """The worst single d-spacing disagreement, which an r.m.s. can hide."""

        if not self.spot_deviations:
            return float("inf")
        return float(max(abs(spot.relative_deviation) for spot in self.spot_deviations))

    @property
    def rms_angle_deviation_deg(self) -> float:
        """Root-mean-square interspot-angle disagreement, in degrees."""

        if not self.angle_deviations:
            return float("inf")
        values = np.asarray([pair.deviation_deg for pair in self.angle_deviations], dtype=float)
        return float(np.sqrt(np.mean(values**2)))

    @property
    def max_angle_deviation_deg(self) -> float:
        """The worst single angular disagreement."""

        if not self.angle_deviations:
            return float("inf")
        return float(max(abs(pair.deviation_deg) for pair in self.angle_deviations))

    @property
    def length_agreement(self) -> float:
        """The length term, in ``[0, 1]``."""

        return self.weights.agreement(
            self.rms_relative_length_deviation, self.weights.length_tolerance
        )

    @property
    def angle_agreement(self) -> float:
        """The angle term, in ``[0, 1]``.

        A pattern with a single indexed spot has no pair to measure an angle
        between. That is missing evidence, not disagreement, so the term is
        neutral at ``0.5`` rather than 0 — scoring it as total disagreement would
        punish a solution for a spot the *user* did not pick.
        """

        if not self.angle_deviations:
            return 0.5
        return self.weights.agreement(
            self.rms_angle_deviation_deg, self.weights.angle_tolerance_deg
        )

    @property
    def coverage_agreement(self) -> float:
        """The coverage term, in ``[0, 1]``: simply the matched fraction."""

        return float(self.matched_fraction)

    @property
    def score(self) -> float:
        """The fused accuracy score in ``[0, 1]``; higher is better.

        A weighted mean of the three terms, normalized by the total weight so
        that the number keeps its meaning when the policy changes: 1 is perfect
        agreement on everything picked, 0.5 is disagreement at the stated
        tolerances, and 0 is nothing explained.
        """

        weights = self.weights
        return float(
            (
                weights.length * self.length_agreement
                + weights.angle * self.angle_agreement
                + weights.coverage * self.coverage_agreement
            )
            / weights.total
        )

    def describe(self) -> str:
        """Prose stating the score, the evidence, and the policy behind it."""

        weights = self.weights
        head = (
            f"Accuracy score {self.score:.3f} of 1: "
            f"lengths agree to {100.0 * self.rms_relative_length_deviation:.2f} percent r.m.s. "
            f"(worst {100.0 * self.max_relative_length_deviation:.2f} percent), "
        )
        angles = (
            "no interspot angle could be measured, so that term is held neutral, "
            if not self.angle_deviations
            else (
                f"angles to {self.rms_angle_deviation_deg:.2f} degrees r.m.s. "
                f"(worst {self.max_angle_deviation_deg:.2f} degrees), "
            )
        )
        coverage = f"and {100.0 * self.matched_fraction:.0f} percent of picked spots are indexed. "
        policy = (
            f"Scored with weights length {weights.length:g}, angle {weights.angle:g}, coverage "
            f"{weights.coverage:g}, where a deviation of {100.0 * weights.length_tolerance:g} "
            f"percent in length or {weights.angle_tolerance_deg:g} degrees in angle scores one "
            "half on its term. Angles are weighted above lengths because they do not depend on "
            "the camera constant, so an angular disagreement is evidence about the "
            "crystallography rather than about the calibration."
        )
        return head + angles + coverage + policy

    def to_json(self) -> dict[str, Any]:
        """The score, its components, its evidence, and its policy."""

        return {
            "score": float(self.score),
            "length_agreement": float(self.length_agreement),
            "angle_agreement": float(self.angle_agreement),
            "coverage_agreement": float(self.coverage_agreement),
            "rms_relative_length_deviation": float(self.rms_relative_length_deviation),
            "max_relative_length_deviation": float(self.max_relative_length_deviation),
            "rms_angle_deviation_deg": float(self.rms_angle_deviation_deg),
            "max_angle_deviation_deg": float(self.max_angle_deviation_deg),
            "matched_fraction": float(self.matched_fraction),
            "weights": self.weights.to_json(),
            "spot_deviations": [spot.to_json() for spot in self.spot_deviations],
            "angle_deviations": [pair.to_json() for pair in self.angle_deviations],
            "describe": self.describe(),
        }


def _angle_between(first: np.ndarray, second: np.ndarray) -> float:
    """Angle in degrees between two planar vectors."""

    norms = float(np.linalg.norm(first)) * float(np.linalg.norm(second))
    if norms <= 0.0:
        return float("nan")
    cosine = float(np.dot(first, second)) / norms
    return float(math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0)))))


def score_solution(
    solution: Any,
    measured_g: Any,
    *,
    weights: ScoringWeights | None = None,
    max_pairs: int = 24,
) -> SolutionScore:
    """Measure a solution's disagreements and fuse them into one score.

    Purpose
    -------
    Turn an indexing into something a user can weigh: what each spot's d-spacing
    should have been, what each pair's angle should have been, how far each is
    out, and a single number that respects a stated policy about which of those
    matters more.

    When and where to use it
    ------------------------
    On every candidate a solver returns, to rank them by something more
    defensible than a sort key, and on the accepted one, to record how well it
    actually fits. The deviations are also what a user reads when a score is
    mediocre and they want to know *why* — a large length deviation with small
    angles points at the camera constant; the reverse points at the phase.

    Parameters
    ----------
    solution : PatternSolution
        As returned by :func:`pytex.diffraction.solving.solve_saed_pattern`. Its
        ``solved_spots`` carry the assigned indices and the predicted in-plane
        ``g`` vector for each.
    measured_g : array_like
        ``(n, 2)`` measured in-plane ``g`` vectors in inverse angstroms, indexed
        as the *measured* spots are — that is, what
        :meth:`pytex.diffraction.solving.MeasuredSAEDPattern.g_vectors_inv_angstrom`
        returns. Each solved spot's ``measured_index`` selects its row.
    weights : ScoringWeights, optional
        The scoring policy. Defaults to :class:`ScoringWeights`.
    max_pairs : int
        Cap on the number of spot pairs whose angles are measured. Pairs grow as
        the square of the spot count, and the first two dozen — taken from the
        shortest, best-determined vectors — carry essentially all the
        information, since a long vector's direction is measured no better than a
        short one's but its angle to another long vector is highly correlated
        with pairs already counted.

    Returns
    -------
    SolutionScore

    Raises
    ------
    ValueError
        If ``measured_g`` is not an ``(n, 2)`` array, if a solved spot points
        outside it, or if ``max_pairs`` is not positive.

    Notes
    -----
    Angles are measured between *measured* vectors and compared with the angle
    between the corresponding *calculated* ones, so the comparison never passes
    through the camera constant and is unaffected by it. That is what makes the
    angle term evidence about the phase rather than about the instrument.
    """

    policy = weights if weights is not None else ScoringWeights()
    if max_pairs <= 0:
        raise ValueError("max_pairs must be strictly positive.")
    vectors = np.asarray(measured_g, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != 2:
        raise ValueError("measured_g must be an (n, 2) array of in-plane g vectors.")

    solved = list(solution.solved_spots)
    for spot in solved:
        if not 0 <= int(spot.measured_index) < vectors.shape[0]:
            raise ValueError(
                f"Solved spot {spot.measured_index} lies outside the {vectors.shape[0]} measured "
                "spots supplied."
            )

    spot_deviations: list[SpotDeviation] = []
    measured_vectors: list[np.ndarray] = []
    calculated_vectors: list[np.ndarray] = []
    for spot in solved:
        index = int(spot.measured_index)
        measured = vectors[index]
        calculated = np.asarray(spot.predicted_g_inv_angstrom, dtype=float).reshape(2)
        indices = tuple(int(value) for value in np.asarray(spot.hkl, dtype=int).reshape(3))
        measured_length = float(np.linalg.norm(measured))
        calculated_length = float(np.linalg.norm(calculated))
        if measured_length <= 0.0 or calculated_length <= 0.0:
            continue
        measured_vectors.append(measured)
        calculated_vectors.append(calculated)
        spot_deviations.append(
            SpotDeviation(
                spot=index + 1,
                hkl=(indices[0], indices[1], indices[2]),
                label=str(spot.label),
                d_measured_angstrom=1.0 / measured_length,
                d_calculated_angstrom=1.0 / calculated_length,
            )
        )

    angle_deviations: list[AngleDeviation] = []
    if len(measured_vectors) >= 2:
        # Shortest vectors first: their directions are the best determined, and a
        # pair of long vectors adds an angle that is largely implied by pairs
        # already counted.
        order = sorted(
            range(len(measured_vectors)), key=lambda index: np.linalg.norm(measured_vectors[index])
        )
        for position, first in enumerate(order):
            for second in order[position + 1 :]:
                if len(angle_deviations) >= max_pairs:
                    break
                measured_angle = _angle_between(
                    measured_vectors[first], measured_vectors[second]
                )
                calculated_angle = _angle_between(
                    calculated_vectors[first], calculated_vectors[second]
                )
                if not (math.isfinite(measured_angle) and math.isfinite(calculated_angle)):
                    continue
                angle_deviations.append(
                    AngleDeviation(
                        first_spot=spot_deviations[first].spot,
                        second_spot=spot_deviations[second].spot,
                        first_label=spot_deviations[first].label,
                        second_label=spot_deviations[second].label,
                        measured_deg=measured_angle,
                        calculated_deg=calculated_angle,
                    )
                )
            if len(angle_deviations) >= max_pairs:
                break

    return SolutionScore(
        spot_deviations=tuple(spot_deviations),
        angle_deviations=tuple(angle_deviations),
        matched_fraction=float(solution.matched_fraction),
        weights=policy,
    )
