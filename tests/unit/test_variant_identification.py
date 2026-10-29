"""Variant identification from one measured pair, checked against constructed truth.

Every child below is *built* to have a known answer: ``C = P V_kᵀ S_c`` for a
chosen variant ``k`` and a chosen child point-group operator ``S_c``, which is
the same grain written in a different but equivalent crystal setting. The
identification must return ``k`` at zero misfit whatever ``S_c`` was, and a
child turned away by a known angle must report about that angle. Agreement with
:func:`pytex.core.transformation.or_deviation`, which assigns variants without
keeping the intermediates, is checked too, so the two paths cannot drift.

The spacing figures are published: the smallest Kurdjumov-Sachs intervariant
disorientation is 10.53 degrees (Morito et al., Acta Mater. 51 (2003) 1789), so
an exact KS child has a runner-up at no less than that.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.app.phases import phase_from_request
from pytex.core.frame_catalog import specimen_frame
from pytex.core.orientation import Orientation, OrientationSet, Rotation
from pytex.core.transformation import OrientationRelationship, or_deviation
from pytex.core.variant_identification import (
    DEFAULT_ASSIGNMENT_TOLERANCE_DEG,
    identify_variant,
    operator_axis_angle,
)

MORITO_SMALLEST_KS_DEG = 10.53
PARENT_EULER = (30.0, 40.0, 10.0)


def _phases(parent: str, child: str):
    _, parent_phase = phase_from_request({"builtin": parent})
    _, child_phase = phase_from_request({"builtin": child})
    return parent_phase, child_phase


@pytest.fixture(scope="module")
def ks():
    parent_phase, child_phase = _phases("austenite_fcc", "fe_bcc")
    return OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )


@pytest.fixture(scope="module")
def burgers():
    parent_phase, child_phase = _phases("zr_bcc_beta", "zr_hcp")
    return OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )


def _orientation(matrix: np.ndarray, phase) -> Orientation:
    return Orientation(
        rotation=Rotation.from_matrix(matrix),
        crystal_frame=phase.crystal_frame,
        specimen_frame=specimen_frame(),
        symmetry=phase.symmetry,
        phase=phase,
    )


def _pair(relationship, variant_position: int, child_operator: int = 0, turn_deg: float = 0.0):
    parent_phase = relationship.parent_phase
    child_phase = relationship.child_phase
    parent_matrix = Rotation.from_bunge_euler(*PARENT_EULER).as_matrix()
    variant = relationship.generate_variants()[variant_position]
    operator = child_phase.symmetry.operators[child_operator]
    child_matrix = parent_matrix @ variant.parent_to_child_rotation.as_matrix().T @ operator
    if turn_deg:
        child_matrix = (
            Rotation.from_axis_angle([1.0, 2.0, 3.0], np.radians(turn_deg)).as_matrix()
            @ child_matrix
        )
    return _orientation(parent_matrix, parent_phase), _orientation(child_matrix, child_phase)


class TestExactChildren:
    @pytest.mark.parametrize("position", [0, 3, 11, 16, 23])
    @pytest.mark.parametrize("child_operator", [0, 5, 17])
    def test_a_constructed_ks_child_is_its_own_variant(self, ks, position, child_operator):
        parent, child = _pair(ks, position, child_operator)
        report = identify_variant(parent, child, ks)
        assert report.variant_index == position + 1
        assert report.misfit_deg < 1e-6
        assert report.verdict == "assigned"

    def test_the_runner_up_is_at_least_the_smallest_intervariant_angle(self, ks):
        parent, child = _pair(ks, 4)
        report = identify_variant(parent, child, ks)
        assert report.runner_up_misfit_deg is not None
        assert report.runner_up_misfit_deg == pytest.approx(MORITO_SMALLEST_KS_DEG, abs=0.01)

    @pytest.mark.parametrize("position", [0, 5, 9])
    def test_a_hexagonal_child_is_identified_under_burgers(self, burgers, position):
        parent, child = _pair(burgers, position, child_operator=7)
        report = identify_variant(parent, child, burgers)
        assert report.variant_index == position + 1
        assert report.misfit_deg < 1e-6

    def test_the_matching_operator_undoes_the_one_the_child_was_written_in(self, ks):
        """``S_c*`` is ``S_c⁻¹``: it returns the child to the variant's own setting."""

        parent, child = _pair(ks, 2, child_operator=9)
        report = identify_variant(parent, child, ks)
        expected = ks.child_phase.symmetry.operators[9].T
        assert np.allclose(report.child_operator, expected, atol=1e-9)
        assert np.allclose(report.residual_matrix, np.eye(3), atol=1e-9)


class TestIntermediates:
    def test_the_measured_rotation_is_ct_p(self, ks):
        parent, child = _pair(ks, 6)
        report = identify_variant(parent, child, ks)
        assert np.allclose(
            report.measured_rotation, child.as_matrix().T @ parent.as_matrix(), atol=1e-12
        )

    def test_the_predicted_child_is_p_vt(self, ks):
        parent, child = _pair(ks, 6)
        report = identify_variant(parent, child, ks)
        assert np.allclose(
            report.predicted_child_matrix,
            report.parent_matrix @ report.variant_matrix.T,
            atol=1e-12,
        )
        assert np.allclose(report.predicted_child_matrix, child.as_matrix(), atol=1e-9)

    def test_the_variant_is_r_times_the_transposed_parent_operator(self, ks):
        parent, child = _pair(ks, 13)
        report = identify_variant(parent, child, ks)
        base = ks.parent_to_child_rotation.as_matrix()
        assert np.allclose(report.variant_matrix, base @ report.parent_operator.T, atol=1e-9)

    def test_every_misfit_is_the_minimum_over_the_child_group(self, ks):
        parent, child = _pair(ks, 1, turn_deg=3.0)
        report = identify_variant(parent, child, ks)
        for position, variant in enumerate(report.variants):
            matrix = variant.parent_to_child_rotation.as_matrix()
            angles = [
                np.degrees(
                    np.arccos(
                        np.clip(
                            (np.trace(matrix @ report.measured_rotation.T @ operator) - 1) / 2,
                            -1,
                            1,
                        )
                    )
                )
                for operator in report.child_operators
            ]
            assert report.misfits_deg[position] == pytest.approx(min(angles), abs=1e-6)


class TestMeasuredDeviation:
    def test_a_child_turned_by_three_degrees_reports_three_degrees(self, ks):
        parent, child = _pair(ks, 8, child_operator=4, turn_deg=3.0)
        report = identify_variant(parent, child, ks)
        assert report.variant_index == 9
        assert report.misfit_deg == pytest.approx(3.0, abs=1e-6)

    def test_the_parallelism_angles_are_bounded_by_the_misfit(self, ks):
        parent, child = _pair(ks, 8, turn_deg=3.0)
        report = identify_variant(parent, child, ks)
        assert {check.kind for check in report.parallelism_checks} == {"plane", "direction"}
        for check in report.parallelism_checks:
            assert check.angle_deg <= report.misfit_deg + 1e-9

    def test_an_exact_child_realizes_its_parallelism_exactly(self, burgers):
        parent, child = _pair(burgers, 3, child_operator=2)
        report = identify_variant(parent, child, burgers)
        for check in report.parallelism_checks:
            assert check.angle_deg < 1e-6

    def test_a_large_turn_is_no_fit(self, ks):
        parent, child = _pair(ks, 0, turn_deg=7.0)
        report = identify_variant(parent, child, ks)
        assert report.verdict == "no_fit"

    def test_a_loose_tolerance_makes_the_answer_ambiguous(self, ks):
        parent, child = _pair(ks, 0)
        report = identify_variant(parent, child, ks, tolerance_deg=12.0)
        assert report.verdict == "ambiguous"


class TestAgreement:
    def test_it_agrees_with_or_deviation(self, ks):
        generator = np.random.default_rng(20261029)
        for position in range(0, 24, 5):
            parent, child = _pair(
                ks,
                position,
                child_operator=int(generator.integers(24)),
                turn_deg=float(generator.uniform(0.0, 2.0)),
            )
            report = identify_variant(parent, child, ks)
            deviation = or_deviation(
                OrientationSet.from_orientations([parent]),
                OrientationSet.from_orientations([child]),
                ks,
            )
            assert report.variant_index == int(deviation.best_variant_indices[0])
            assert report.misfit_deg == pytest.approx(float(deviation.deviations_deg[0]), abs=1e-6)


class TestContract:
    def test_mismatched_phases_are_refused(self, ks, burgers):
        parent, child = _pair(burgers, 0)
        with pytest.raises(ValueError, match="phase"):
            identify_variant(parent, child, ks)

    def test_a_non_positive_tolerance_is_refused(self, ks):
        parent, child = _pair(ks, 0)
        with pytest.raises(ValueError, match="tolerance"):
            identify_variant(parent, child, ks, tolerance_deg=0.0)

    def test_describe_and_json_agree(self, ks):
        parent, child = _pair(ks, 10, turn_deg=1.0)
        report = identify_variant(parent, child, ks)
        record = report.to_json_dict()
        assert record["variant_index"] == 11
        assert record["verdict"] == "assigned"
        assert "variant 11" in report.describe()
        assert record["describe"] == report.describe()
        assert report.tolerance_deg == DEFAULT_ASSIGNMENT_TOLERANCE_DEG

    def test_operator_axis_angle_reads_a_fourfold(self):
        axis, angle = operator_axis_angle(
            Rotation.from_axis_angle([0, 0, 1], np.pi / 2).as_matrix()
        )
        assert angle == pytest.approx(90.0)
        assert np.allclose(axis, [0, 0, 1])
        axis, angle = operator_axis_angle(np.eye(3))
        assert angle == 0.0
