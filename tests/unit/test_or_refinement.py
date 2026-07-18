from __future__ import annotations

import numpy as np
import pytest

from pytex.core import OrientationRelationship
from pytex.core.transformation import _symmetry_reduced_angle_between_deg
from pytex.experimental import refine_orientation_relationship_from_boundaries
from tests.unit.test_or_identification import _microstructure
from tests.unit.test_parent_grain_reconstruction import _phases


def test_refinement_recovers_gt_rotation_from_ks_nominal() -> None:
    parent_phase, child_phase, _ = _phases()
    children, edges = _microstructure(
        OrientationRelationship.from_greninger_troiano_correspondence
    )
    nominal = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    report = refine_orientation_relationship_from_boundaries(children, edges, nominal)
    assert report.converged
    # Exact GT boundaries: the refined rotation reproduces them exactly (to
    # the ~1e-6 deg matrix<->quaternion round-trip noise floor).
    assert report.refined_mean_distance_deg == pytest.approx(0.0, abs=1e-5)
    assert float(np.max(report.edge_distances_deg)) == pytest.approx(0.0, abs=1e-5)
    # The refined rotation IS the GT rotation (symmetry-reduced), recovered
    # from boundaries alone starting 2.4 deg away at KS.
    true_gt = OrientationRelationship.from_greninger_troiano_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    distance_to_truth = _symmetry_reduced_angle_between_deg(
        report.relationship.parent_to_child_rotation.as_matrix(),
        true_gt.parent_to_child_rotation.as_matrix(),
        child_operators=child_phase.symmetry.operators,
        parent_operators=parent_phase.symmetry.operators,
    )
    assert distance_to_truth == pytest.approx(0.0, abs=1e-5)
    # The reported update matches the documented KS-GT separation.
    assert report.rotation_update_deg == pytest.approx(2.404, abs=0.01)
    assert report.initial_mean_distance_deg > 1.0
    assert report.relationship.name == "kurdjumov_sachs_refined"


def test_refinement_with_noise_stays_near_truth_and_describes() -> None:
    parent_phase, child_phase, _ = _phases()
    children, edges = _microstructure(
        OrientationRelationship.from_kurdjumov_sachs_correspondence, noise_deg=0.3
    )
    nominal = OrientationRelationship.from_nishiyama_wassermann_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    report = refine_orientation_relationship_from_boundaries(children, edges, nominal)
    true_ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    distance_to_truth = _symmetry_reduced_angle_between_deg(
        report.relationship.parent_to_child_rotation.as_matrix(),
        true_ks.parent_to_child_rotation.as_matrix(),
        child_operators=child_phase.symmetry.operators,
        parent_operators=parent_phase.symmetry.operators,
    )
    # Boundary noise limits the fit; the refined rotation must still be far
    # closer to the true KS than the NW nominal (2.86 deg away) is.
    assert distance_to_truth < 0.5
    assert report.refined_mean_distance_deg < 1.0
    assert report.refined_mean_distance_deg < report.initial_mean_distance_deg
    text = report.describe()
    for expected in (
        "Boundary-based refinement",
        "no parent orientations",
        "symmetry-reduced",
        "Experimental surface",
    ):
        assert expected in text


def test_refinement_input_validation() -> None:
    parent_phase, child_phase, _ = _phases()
    children, edges = _microstructure(
        OrientationRelationship.from_kurdjumov_sachs_correspondence
    )
    nominal = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    with pytest.raises(ValueError, match="at least one edge"):
        refine_orientation_relationship_from_boundaries(
            children, np.empty((0, 2), dtype=np.int64), nominal
        )
    with pytest.raises(ValueError, match="self-edges"):
        refine_orientation_relationship_from_boundaries(
            children, np.array([[0, 0]]), nominal
        )
    with pytest.raises(ValueError, match="must reference"):
        refine_orientation_relationship_from_boundaries(
            children, np.array([[0, 99]]), nominal
        )
