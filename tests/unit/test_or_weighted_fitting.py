"""Independent circular-mean identities and measurement-evidence contracts."""

from dataclasses import replace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    OrientationRelationship,
    OrientationSet,
    Rotation,
    characterize_orientation_relationship,
    fit_orientation_relationship,
    orientation_relationship_from_euler,
    specimen_frame,
)
from tests.unit.test_or_characterization import _planted_pairs
from tests.unit.test_transformation import make_phases


def measured_pairs(angles_deg: tuple[float, ...] = (0.0, 2.0)) -> tuple:
    _, _, parent, child = make_phases()
    nominal = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    measured = np.stack([
        Rotation.from_axis_angle([1, 0, 0], np.radians(angle)).as_matrix()
        @ nominal.parent_to_child_rotation.as_matrix() for angle in angles_deg
    ])
    parents = OrientationSet.from_matrices(
        np.repeat(np.eye(3)[None], len(measured), axis=0),
        specimen_frame=specimen_frame(), phase=parent,
    )
    children = OrientationSet.from_matrices(
        measured.transpose(0, 2, 1), specimen_frame=specimen_frame(), phase=child,
    )
    return parents, children, nominal


def test_weighted_fit_matches_analytic_circular_mean() -> None:
    parents, children, nominal = measured_pairs()
    report = fit_orientation_relationship(parents, children, nominal, pair_weights=[3, 1])
    # For a common axis Markley's objective reduces to the weighted circular mean.
    theta = np.arctan2(np.sin(np.radians(2)), 3 + np.cos(np.radians(2)))
    expected = Rotation.from_axis_angle([1, 0, 0], theta).as_matrix()
    expected = expected @ nominal.parent_to_child_rotation.as_matrix()
    assert_allclose(report.relationship.parent_to_child_rotation.as_matrix(), expected, atol=1e-12)
    residuals = [np.degrees(theta), 2 - np.degrees(theta)]
    assert_allclose(report.residuals_deg, residuals, atol=1e-10)
    assert_allclose(report.weighted_mean_residual_deg, np.average(residuals, weights=[3, 1]))
    assert report.effective_pair_count == pytest.approx(1.6)
    assert "Markley" in report.describe()


@pytest.mark.parametrize("weights", [[0, 0], [-1, 1], [np.nan, 1], [np.inf, 1], [1], [[1, 1]]])
def test_invalid_weights_rejected(weights: list) -> None:
    parents, children, nominal = measured_pairs()
    for function, args in (
        (fit_orientation_relationship, (parents, children, nominal)),
        (characterize_orientation_relationship, (parents, children)),
    ):
        with pytest.raises(ValueError, match="pair_weights"):
            function(*args, pair_weights=weights)


def test_weights_are_scale_invariant_owned_and_read_only() -> None:
    parents, children, nominal = measured_pairs()
    weights = np.array([1e308, 1e308])
    report = fit_orientation_relationship(parents, children, nominal, pair_weights=weights)
    ordinary = fit_orientation_relationship(parents, children, nominal)
    weights[:] = 0
    assert_allclose(report.pair_weights, [0.5, 0.5])
    assert_allclose(report.residuals_deg, ordinary.residuals_deg, atol=1e-12)
    assert not report.pair_weights.flags.writeable


def test_zero_weight_first_pair_does_not_seed_or_bias_characterization() -> None:
    parents, children, nominal = measured_pairs((12.0, 0.0, 0.0))
    report = characterize_orientation_relationship(
        parents, children, pair_weights=[0, 1, 1], max_index=1, max_statements=1
    )
    assert report.best_catalog_name == nominal.name
    assert report.best_catalog_deviation_deg < 1e-10
    assert report.residuals_deg[0] > 1
    assert report.weighted_mean_residual_deg < 1e-10
    assert report.mean_residual_deg > 0.1
    assert report.is_conclusive
    assert report.to_json_dict()["pair_weights"] == [0, 0.5, 0.5]
    assert report.to_json_dict()["residuals_deg"] == report.residuals_deg.tolist()


@pytest.mark.parametrize("iterations", [0, -1, 1.5, True])
def test_invalid_iteration_budget(iterations: int) -> None:
    parents, children, nominal = measured_pairs()
    with pytest.raises(ValueError, match="max_iterations"):
        fit_orientation_relationship(parents, children, nominal, max_iterations=iterations)


@pytest.mark.parametrize("tolerance", [-1, np.nan, np.inf])
def test_invalid_convergence_tolerance(tolerance: float) -> None:
    parents, children, nominal = measured_pairs()
    with pytest.raises(ValueError, match="convergence_tol_deg"):
        fit_orientation_relationship(parents, children, nominal, convergence_tol_deg=tolerance)


def test_chunk_size_does_not_change_multi_variant_fit(monkeypatch: pytest.MonkeyPatch) -> None:
    import pytex.core.transformation as module

    _, _, nominal = measured_pairs()
    parents, children = _planted_pairs(
        nominal, variant_indices=tuple(range(1, 25)), scatter_deg=0.3, seed=3
    )
    full = fit_orientation_relationship(parents, children, nominal)
    monkeypatch.setattr(module, "_OR_ALIGNMENT_CHUNK_SIZE", 3)
    chunked = fit_orientation_relationship(parents, children, nominal)
    assert_allclose(chunked.residuals_deg, full.residuals_deg, atol=1e-10)
    assert_allclose(chunked.relationship.parent_to_child_rotation.as_matrix(),
                    full.relationship.parent_to_child_rotation.as_matrix(), atol=1e-12)


def test_unconverged_or_scattered_single_catalog_fit_cannot_be_conclusive() -> None:
    parents, children, nominal = measured_pairs()
    report = characterize_orientation_relationship(
        parents, children, nominal=nominal, catalog=(nominal,), max_index=1
    )
    assert report.is_conclusive
    assert not replace(report, converged=False).is_conclusive
    assert "NOT converge" in replace(report, converged=False).describe()
    assert not replace(report, residuals_deg=np.array([10.0, 10.0])).is_conclusive


def test_nominal_and_catalog_must_match_measured_phases() -> None:
    parents, children, nominal = measured_pairs()
    wrong = OrientationRelationship(
        name="wrong phases", parent_phase=nominal.child_phase, child_phase=nominal.parent_phase,
        parent_to_child_rotation=Rotation.identity(),
    )
    for kwargs in ({"nominal": wrong}, {"catalog": (wrong,)}):
        with pytest.raises(ValueError, match="match the measured phases"):
            characterize_orientation_relationship(parents, children, **kwargs)


def test_exhausted_budget_reports_final_fit_residuals_and_no_conclusion() -> None:
    parents, children, nominal = measured_pairs()
    report = characterize_orientation_relationship(
        parents, children, nominal=nominal, pair_weights=[3, 1], max_iterations=1,
        catalog=(nominal,), max_index=1,
    )
    theta = np.degrees(np.arctan2(np.sin(np.radians(2)), 3 + np.cos(np.radians(2))))
    assert not report.converged
    assert not report.is_conclusive
    assert_allclose(report.residuals_deg, [theta, 2-theta], atol=1e-10)


def test_euler_entry_point_preserves_weights_and_exclusions() -> None:
    parents, children, nominal = measured_pairs()
    report = orientation_relationship_from_euler(
        parents.as_bunge_euler(), children.as_bunge_euler(),
        parent_phase=nominal.parent_phase, child_phase=nominal.child_phase,
        pair_weights=[1, 0], max_index=1,
    )
    assert_allclose(report.pair_weights, [1, 0])
    assert_allclose(report.residuals_deg, [0, 2], atol=1e-10)


@pytest.mark.parametrize("name", ["catalog_tolerance_deg", "parallelism_tolerance_deg"])
@pytest.mark.parametrize("value", [-1, np.inf, np.nan])
def test_characterization_tolerances_are_validated(name: str, value: float) -> None:
    parents, children, _ = measured_pairs()
    with pytest.raises(ValueError, match=name):
        characterize_orientation_relationship(parents, children, **{name: value})
