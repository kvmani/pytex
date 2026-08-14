"""Scoring an indexed solution: the evidence, and the policy that fuses it.

The tests split the same way the module does. The deviations are measurements
and are checked against constructions — a spot placed five percent off must
report five percent. The score is a policy and is checked for the properties a
policy must have: bounded, monotone in each disagreement, unchanged by anything
it should not depend on, and carrying the weights that produced it.

The sharpest test here is the calibration one. A wrong camera constant scales
every length and leaves every angle alone, so the angle term must not move at
all — that invariance is the whole reason angles are weighted above lengths.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import pairwise

import numpy as np
import pytest

from pytex.diffraction.solution_scoring import (
    AngleDeviation,
    ScoringWeights,
    SolutionScore,
    SpotDeviation,
    score_solution,
)


@dataclass(frozen=True)
class FakeSpot:
    """The three fields `score_solution` reads from a solved spot."""

    measured_index: int
    hkl: tuple[int, int, int]
    label: str
    predicted_g_inv_angstrom: tuple[float, float]


@dataclass(frozen=True)
class FakeSolution:
    solved_spots: tuple[FakeSpot, ...]
    matched_fraction: float = 1.0


def square_solution(scale: float = 1.0, *, rotate_deg: float = 0.0):
    """Four spots on a square reciprocal lattice, measurable exactly.

    ``scale`` stretches the *measured* vectors relative to the calculated ones,
    which is precisely what a mis-set camera constant does.
    """

    calculated = np.array([[0.5, 0.0], [0.0, 0.5], [0.5, 0.5], [1.0, 0.0]])
    angle = math.radians(rotate_deg)
    rotation = np.array(
        [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
    )
    measured = scale * (calculated @ rotation.T)
    solution = FakeSolution(
        solved_spots=tuple(
            FakeSpot(
                measured_index=index,
                hkl=(2, 0, 0),
                label=f"g{index + 1}",
                predicted_g_inv_angstrom=tuple(calculated[index]),
            )
            for index in range(len(calculated))
        )
    )
    return solution, measured


# ------------------------------------------------------------- the policy


def test_the_agreement_curve_scores_a_half_at_the_tolerance() -> None:
    weights = ScoringWeights()
    assert weights.agreement(0.0, 0.02) == pytest.approx(1.0)
    assert weights.agreement(0.02, 0.02) == pytest.approx(0.5)
    assert weights.agreement(0.04, 0.02) < 0.5
    assert weights.agreement(1.0, 0.02) > 0.0


def test_the_agreement_curve_is_monotone_and_bounded() -> None:
    weights = ScoringWeights()
    values = [weights.agreement(deviation, 0.02) for deviation in np.linspace(0.0, 0.5, 40)]
    assert all(0.0 <= value <= 1.0 for value in values)
    assert all(later <= earlier for earlier, later in pairwise(values))


def test_a_non_finite_deviation_scores_zero_rather_than_raising() -> None:
    assert ScoringWeights().agreement(float("inf"), 0.02) == 0.0


def test_the_default_policy_weights_angles_above_lengths() -> None:
    """Angles are calibration-free, so they are the better evidence about phase."""

    weights = ScoringWeights()
    assert weights.angle > weights.length
    assert weights.coverage >= weights.angle
    assert weights.total == pytest.approx(weights.length + weights.angle + weights.coverage)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"length": -1.0},
        {"angle": float("nan")},
        {"length": 0.0, "angle": 0.0, "coverage": 0.0},
        {"length_tolerance": 0.0},
        {"angle_tolerance_deg": -1.0},
        {"sharpness": 0.0},
    ],
)
def test_an_impossible_policy_is_refused(kwargs) -> None:
    with pytest.raises(ValueError):
        ScoringWeights(**kwargs)


def test_a_policy_can_be_adjusted_without_losing_its_validation() -> None:
    weights = ScoringWeights().replace(angle=3.0)
    assert weights.angle == 3.0
    with pytest.raises(ValueError):
        ScoringWeights().replace(length=-2.0)


# ---------------------------------------------------------- the deviations


def test_a_perfect_indexing_scores_essentially_one() -> None:
    solution, measured = square_solution()
    score = score_solution(solution, measured)
    assert score.rms_relative_length_deviation == pytest.approx(0.0, abs=1e-12)
    assert score.rms_angle_deviation_deg == pytest.approx(0.0, abs=1e-9)
    assert score.score == pytest.approx(1.0, abs=1e-9)


def test_a_stretched_pattern_reports_the_stretch_as_a_length_deviation() -> None:
    """A five percent camera-constant error is a five percent length deviation."""

    solution, measured = square_solution(scale=1.05)
    score = score_solution(solution, measured)
    # d = 1/|g|, so measured g larger by 1.05 makes measured d smaller by 1/1.05.
    expected = 1.0 / 1.05 - 1.0
    for deviation in score.spot_deviations:
        assert deviation.relative_deviation == pytest.approx(expected, rel=1e-9)
    assert score.rms_relative_length_deviation == pytest.approx(abs(expected), rel=1e-9)


def test_the_angle_term_does_not_move_when_the_calibration_does() -> None:
    """The invariance the whole weighting rests on."""

    baseline = score_solution(*square_solution())
    stretched = score_solution(*square_solution(scale=1.05))
    squashed = score_solution(*square_solution(scale=0.8))
    assert stretched.rms_angle_deviation_deg == pytest.approx(
        baseline.rms_angle_deviation_deg, abs=1e-9
    )
    assert squashed.rms_angle_deviation_deg == pytest.approx(
        baseline.rms_angle_deviation_deg, abs=1e-9
    )
    assert stretched.angle_agreement == pytest.approx(baseline.angle_agreement, abs=1e-9)
    # And the fused score does fall, because lengths are part of it.
    assert stretched.score < baseline.score


def test_a_rotated_pattern_is_not_penalised() -> None:
    """One pattern cannot fix the roll about the beam; scoring must not pretend."""

    baseline = score_solution(*square_solution())
    rolled = score_solution(*square_solution(rotate_deg=37.0))
    assert rolled.rms_angle_deviation_deg == pytest.approx(
        baseline.rms_angle_deviation_deg, abs=1e-9
    )
    assert rolled.score == pytest.approx(baseline.score, abs=1e-9)


def test_a_sheared_pattern_shows_up_in_the_angles_not_the_lengths() -> None:
    solution, measured = square_solution()
    # Move one spot along its perpendicular: same length, different angle.
    measured = measured.copy()
    measured[1] = np.array([0.08, 0.5])
    score = score_solution(solution, measured)
    assert score.rms_angle_deviation_deg > 5.0
    assert score.max_angle_deviation_deg >= score.rms_angle_deviation_deg


def test_coverage_is_the_matched_fraction() -> None:
    solution, measured = square_solution()
    partial = FakeSolution(solved_spots=solution.solved_spots, matched_fraction=0.5)
    score = score_solution(partial, measured)
    assert score.coverage_agreement == pytest.approx(0.5)
    assert score.score < score_solution(solution, measured).score


def test_a_single_spot_leaves_the_angle_term_neutral() -> None:
    """Missing evidence is not disagreement, and must not be scored as it."""

    solution, measured = square_solution()
    single = FakeSolution(solved_spots=solution.solved_spots[:1])
    score = score_solution(single, measured)
    assert score.angle_deviations == ()
    assert score.angle_agreement == pytest.approx(0.5)
    assert 0.0 < score.score < 1.0


def test_pairs_are_capped_and_taken_from_the_shortest_vectors_first() -> None:
    solution, measured = square_solution()
    capped = score_solution(solution, measured, max_pairs=2)
    assert len(capped.angle_deviations) == 2
    full = score_solution(solution, measured)
    assert len(full.angle_deviations) == 6  # four spots, all pairs
    first = full.angle_deviations[0]
    assert {first.first_spot, first.second_spot} == {1, 2}


# --------------------------------------------------------- the presentation


def test_the_score_carries_the_policy_that_produced_it() -> None:
    solution, measured = square_solution()
    strict = ScoringWeights(length_tolerance=0.001)
    score = score_solution(solution, measured, weights=strict)
    assert score.weights == strict
    assert score.to_json()["weights"]["length_tolerance"] == pytest.approx(0.001)


def test_a_stricter_policy_scores_the_same_evidence_lower() -> None:
    solution, measured = square_solution(scale=1.02)
    lenient = score_solution(solution, measured, weights=ScoringWeights(length_tolerance=0.10))
    strict = score_solution(solution, measured, weights=ScoringWeights(length_tolerance=0.002))
    assert strict.score < lenient.score
    # The evidence itself is identical; only the reading of it changed.
    assert strict.rms_relative_length_deviation == pytest.approx(
        lenient.rms_relative_length_deviation
    )


def test_describe_states_the_score_the_evidence_and_the_policy() -> None:
    text = score_solution(*square_solution(scale=1.03)).describe()
    assert "Accuracy score" in text
    assert "percent r.m.s." in text
    assert "camera constant" in text


def test_json_carries_components_evidence_and_weights() -> None:
    payload = score_solution(*square_solution()).to_json()
    assert set(payload) >= {
        "score",
        "length_agreement",
        "angle_agreement",
        "coverage_agreement",
        "rms_relative_length_deviation",
        "max_relative_length_deviation",
        "rms_angle_deviation_deg",
        "max_angle_deviation_deg",
        "matched_fraction",
        "weights",
        "spot_deviations",
        "angle_deviations",
        "describe",
    }
    assert set(payload["spot_deviations"][0]) == {
        "spot",
        "hkl",
        "label",
        "d_measured",
        "d_calculated",
        "deviation_angstrom",
        "relative_deviation",
    }
    assert set(payload["angle_deviations"][0]) == {
        "first_spot",
        "second_spot",
        "pair",
        "measured_deg",
        "calculated_deg",
        "deviation_deg",
    }


def test_the_score_stays_inside_its_stated_range() -> None:
    for scale in (0.2, 0.9, 1.0, 1.1, 5.0):
        score = score_solution(*square_solution(scale=scale))
        assert 0.0 <= score.score <= 1.0


# ----------------------------------------------------------- the guardrails


def test_a_malformed_measurement_array_is_refused() -> None:
    solution, _ = square_solution()
    with pytest.raises(ValueError, match=r"\(n, 2\)"):
        score_solution(solution, np.zeros((4, 3)))


def test_a_solved_spot_outside_the_measurements_is_refused() -> None:
    solution, measured = square_solution()
    with pytest.raises(ValueError, match="outside"):
        score_solution(solution, measured[:2])


def test_max_pairs_must_be_positive() -> None:
    solution, measured = square_solution()
    with pytest.raises(ValueError, match="strictly positive"):
        score_solution(solution, measured, max_pairs=0)


def test_the_row_types_validate_their_own_contents() -> None:
    with pytest.raises(ValueError, match="matched_fraction"):
        SolutionScore(
            spot_deviations=(), angle_deviations=(), matched_fraction=1.5, weights=ScoringWeights()
        )
    spot = SpotDeviation(1, (2, 0, 0), "200", 2.0, 2.02)
    assert spot.relative_deviation == pytest.approx(-0.02 / 2.02, rel=1e-12)
    pair = AngleDeviation(1, 2, "200", "020", 90.5, 90.0)
    assert pair.deviation_deg == pytest.approx(0.5)


# ------------------------------------------------- against a simulated plate


@pytest.mark.parametrize(
    ("entry", "phase_key"),
    [("fcc_al_001", "al_fcc"), ("bcc_fe_110", "fe_bcc"), ("hcp_zr_2-1-10", "zr_hcp")],
)
def test_a_practice_plate_scores_near_one(entry: str, phase_key: str) -> None:
    """End to end: a correct indexing of an exact pattern must score highly."""

    pytest.importorskip("matplotlib", reason="the diffraction stack pulls in the plotting layer")
    from pytex.app import REGISTRY
    from pytex.app.phases import builtin_phase
    from pytex.diffraction.solving import (
        MeasuredSAEDPattern,
        MeasuredSpot,
        PatternCalibration,
        solve_saed_pattern,
    )

    opened = REGISTRY.call("tem.gallery_pattern", {"pattern": entry, "realistic_scatter": False})
    calibration = opened["data"]["calibration"]
    picks = opened["data"]["suggested_picks"]
    pattern = MeasuredSAEDPattern(
        name=entry,
        spots=tuple(MeasuredSpot(position=(spot["x"], spot["y"])) for spot in picks["spots"]),
        calibration=PatternCalibration(
            units="px",
            centre=tuple(picks["centre"]),
            camera_constant_mm_angstrom=calibration["camera_constant_mm_angstrom"],
            pixel_size_mm=calibration["pixel_size_mm"],
        ),
    )
    report = solve_saed_pattern(pattern, [builtin_phase(phase_key).to_phase()], max_index=4)
    score = score_solution(
        report.best(), np.asarray(pattern.g_vectors_inv_angstrom(), dtype=float)
    )
    assert score.score > 0.98
    assert score.rms_relative_length_deviation < 0.01
    assert score.rms_angle_deviation_deg < 0.5
