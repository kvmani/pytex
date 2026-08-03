"""TX1: determining the orientation relationship from measured orientations.

Every expected value here has independent provenance — a definitional
parallelism of a named relationship, a published angular separation between two
named relationships, or an analytic identity — never a copied program output.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    OrientationRelationship,
    OrientationSet,
    Phase,
    Rotation,
    characterize_orientation_relationship,
    default_relationship_catalog,
    describe_orientation_relationship,
    orientation_relationship_from_euler,
    specimen_frame,
)
from tests.unit.test_composite_saed import make_bcc_hcp_phases
from tests.unit.test_transformation import make_phases

# The rotation used as the parent orientation throughout. Nothing depends on
# its value — the relationship is a crystal-to-crystal object — but pinning it
# keeps the tests deterministic.
PARENT_ROTATION = Rotation.from_axis_angle([1.0, 2.0, 3.0], 0.7)

#: Published angular separations between named fcc->bcc relationships, used to
#: check that identification discriminates rather than guesses. Kurdjumov-Sachs
#: and Nishiyama-Wassermann differ by 5.26 deg; Greninger-Troiano sits between
#: them, 2.40 deg from KS.
KS_TO_NW_DEG = 5.26
KS_TO_GT_DEG = 2.40


def _fcc_bcc() -> tuple[Phase, Phase]:
    _, _, parent, child = make_phases()
    return parent, child


def _planted_pairs(
    relationship: OrientationRelationship,
    *,
    variant_indices: tuple[int, ...],
    scatter_deg: float = 0.0,
    seed: int = 0,
) -> tuple[OrientationSet, OrientationSet]:
    """Parent/child orientation pairs built through known variants of an OR.

    The parent orientation is shared; each child is the prediction of one
    variant, optionally perturbed by a random rotation of the given magnitude
    so noise robustness can be exercised against known ground truth.
    """

    frame = specimen_frame()
    variants = relationship.generate_variants()
    available = {variant.variant_index: variant for variant in variants}
    parent_matrix = PARENT_ROTATION.as_matrix()
    rng = np.random.default_rng(seed)
    child_matrices = []
    for index in variant_indices:
        # Canonical convention: C = P V^T.
        matrix = parent_matrix @ available[index].parent_to_child_rotation.as_matrix().T
        if scatter_deg > 0.0:
            axis = rng.normal(size=3)
            axis /= np.linalg.norm(axis)
            angle = np.radians(rng.normal(0.0, scatter_deg))
            matrix = Rotation.from_axis_angle(axis, angle).as_matrix() @ matrix
        child_matrices.append(matrix)
    parents = OrientationSet.from_matrices(
        np.stack([parent_matrix] * len(variant_indices)),
        specimen_frame=frame,
        phase=relationship.parent_phase,
    )
    children = OrientationSet.from_matrices(
        np.stack(child_matrices), specimen_frame=frame, phase=relationship.child_phase
    )
    return parents, children


class TestParallelismExtraction:
    """describe_orientation_relationship recovers the defining statement."""

    def test_kurdjumov_sachs_statement_is_close_packed_planes_and_directions(self) -> None:
        """KS is defined by {111}_fcc || {011}_bcc and <110>_fcc || <111>_bcc.

        Definitional (Kurdjumov and Sachs 1930), so the recovered statement must
        name those families at zero deviation.
        """

        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        planes, directions = describe_orientation_relationship(relationship)
        assert planes and directions
        plane = planes[0]
        direction = directions[0]
        assert plane.deviation_deg < 1e-4
        assert direction.deviation_deg < 1e-4
        # {111} parent: all three indices of unit magnitude. {011} child: one zero.
        assert sorted(np.abs(plane.parent_indices).tolist()) == [1, 1, 1]
        assert sorted(np.abs(plane.child_indices).tolist()) == [0, 1, 1]
        # <110> parent, <111> child.
        assert sorted(np.abs(direction.parent_indices).tolist()) == [0, 1, 1]
        assert sorted(np.abs(direction.child_indices).tolist()) == [1, 1, 1]

    def test_nishiyama_wassermann_statement_pairs_close_packed_plane_with_cube_direction(
        self,
    ) -> None:
        """NW shares KS's {111} || {011} plane but pairs <110>_fcc with <100>_bcc."""

        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_nishiyama_wassermann_correspondence(
            parent_phase=parent, child_phase=child
        )
        planes, directions = describe_orientation_relationship(relationship)
        assert sorted(np.abs(planes[0].parent_indices).tolist()) == [1, 1, 1]
        assert sorted(np.abs(planes[0].child_indices).tolist()) == [0, 1, 1]
        assert sorted(np.abs(directions[0].parent_indices).tolist()) == [0, 1, 1]
        assert sorted(np.abs(directions[0].child_indices).tolist()) == [0, 0, 1]

    def test_burgers_statement_uses_four_index_hexagonal_labels(self) -> None:
        """Burgers is (011)_bcc || (0001)_hcp with <111>_bcc || <11-20>_hcp.

        The hexagonal side must be labeled in four-index Miller-Bravais form,
        which is how the hcp literature states it.
        """

        parent, child = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=parent, child_phase=child
        )
        planes, directions = describe_orientation_relationship(relationship)
        plane = planes[0]
        direction = directions[0]
        assert plane.deviation_deg < 1e-4
        assert direction.deviation_deg < 1e-4
        assert sorted(np.abs(plane.parent_indices).tolist()) == [0, 1, 1]
        # (0001) basal plane of the hexagonal child.
        assert plane.child_label.replace(" ", "") == "(0001)"
        assert sorted(np.abs(direction.parent_indices).tolist()) == [1, 1, 1]
        # A <11-20> direction: the four-index label has three nonzero entries
        # summing (over the first three) to zero and a zero c component.
        assert direction.child_label.count(" ") == 3

    def test_statements_are_deduplicated_by_index_family(self) -> None:
        """The same statement written with a different family member is not new."""

        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        planes, _ = describe_orientation_relationship(relationship, max_statements=6)
        keys = {
            (
                tuple(sorted(np.abs(statement.parent_indices).tolist())),
                tuple(sorted(np.abs(statement.child_indices).tolist())),
            )
            for statement in planes
        }
        assert len(keys) == len(planes)

    def test_tight_tolerance_admits_only_exact_parallelisms(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        planes, directions = describe_orientation_relationship(
            relationship, tolerance_deg=1e-3, max_statements=8
        )
        for statement in planes + directions:
            assert statement.deviation_deg <= 1e-3

    def test_rejects_non_positive_statement_count(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_bain_correspondence(
            parent_phase=parent, child_phase=child
        )
        with pytest.raises(ValueError, match="max_statements must be at least 1"):
            describe_orientation_relationship(relationship, max_statements=0)


class TestCatalogDispatch:
    def test_cubic_cubic_resolves_to_the_fcc_bcc_family(self) -> None:
        parent, child = _fcc_bcc()
        catalog = default_relationship_catalog(parent_phase=parent, child_phase=child)
        assert catalog is not None
        assert set(catalog.names()) == {
            "bain",
            "kurdjumov_sachs",
            "nishiyama_wassermann",
            "greninger_troiano",
            "pitsch",
        }

    def test_cubic_hexagonal_resolves_to_burgers_and_shoji_nishiyama(self) -> None:
        parent, child = make_bcc_hcp_phases()
        catalog = default_relationship_catalog(parent_phase=parent, child_phase=child)
        assert catalog is not None
        assert set(catalog.names()) == {"burgers", "shoji_nishiyama"}

    def test_unsupported_pair_returns_none_rather_than_a_wrong_catalog(self) -> None:
        """A hexagonal->hexagonal pair has no standard catalog; say so."""

        _, hexagonal = make_bcc_hcp_phases()
        assert (
            default_relationship_catalog(parent_phase=hexagonal, child_phase=hexagonal)
            is None
        )


class TestCharacterization:
    def test_exact_kurdjumov_sachs_pairs_are_identified_conclusively(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(
            relationship, variant_indices=(1, 5, 9, 14, 18, 23)
        )
        report = characterize_orientation_relationship(parents, children)
        assert report.best_catalog_name == "kurdjumov_sachs"
        assert report.best_catalog_deviation_deg < 1e-4
        assert report.mean_residual_deg < 1e-4
        assert report.matches_catalog
        assert report.is_conclusive
        assert report.converged

    def test_the_fit_is_seeded_from_the_data_not_from_a_nominal(self) -> None:
        """No nominal relationship is supplied anywhere in these tests.

        Children drawn from six *different* variants still reduce to one
        rotation, which is the property the double-coset seed relies on.
        """

        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_nishiyama_wassermann_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 2, 3, 4, 5, 6))
        report = characterize_orientation_relationship(parents, children)
        assert report.best_catalog_name == "nishiyama_wassermann"
        assert report.best_catalog_deviation_deg < 1e-4

    def test_named_relationships_are_recovered_from_their_own_variants(self) -> None:
        parent, child = _fcc_bcc()
        builders = {
            "bain": OrientationRelationship.from_bain_correspondence,
            "kurdjumov_sachs": (
                OrientationRelationship.from_kurdjumov_sachs_correspondence
            ),
            "nishiyama_wassermann": (
                OrientationRelationship.from_nishiyama_wassermann_correspondence
            ),
            "greninger_troiano": (
                OrientationRelationship.from_greninger_troiano_correspondence
            ),
            "pitsch": OrientationRelationship.from_pitsch_correspondence,
        }
        for name, builder in builders.items():
            relationship = builder(parent_phase=parent, child_phase=child)
            count = min(4, len(relationship.generate_variants()))
            parents, children = _planted_pairs(
                relationship, variant_indices=tuple(range(1, count + 1))
            )
            report = characterize_orientation_relationship(parents, children)
            assert report.best_catalog_name == name, name
            assert report.best_catalog_deviation_deg < 1e-4, name

    def test_bain_survives_the_double_coset_tie(self) -> None:
        """Regression: a symmetric relationship must not be averaged into nonsense.

        Bain is 45 deg about <100> and has only three variants. Its maximum-trace
        double-coset representative is *not* unique, so reducing every pair
        independently and averaging the results lands on a rotation none of them
        shows (a measured 26.9 deg, misread as Kurdjumov-Sachs). Seeding from one
        pair and aligning the rest against it is what fixes this, and this test
        fails if that changes back.
        """

        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_bain_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 2, 3))
        report = characterize_orientation_relationship(parents, children)
        assert report.best_catalog_name == "bain"
        assert report.best_catalog_deviation_deg < 1e-4
        # The Bain correspondence is a 45 deg rotation; the fit must land there.
        assert_allclose(report.relationship.misorientation().angle_deg, 45.0, atol=1e-3)

    def test_catalog_separation_matches_the_published_angles(self) -> None:
        """KS-to-NW is 5.26 deg and KS-to-GT is 2.40 deg in the literature.

        Ranking a planted KS dataset must reproduce those separations, which is
        what makes the margin meaningful.
        """

        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 2, 3))
        report = characterize_orientation_relationship(parents, children)
        deviations = dict(
            zip(report.catalog_names, report.catalog_deviations_deg, strict=True)
        )
        assert_allclose(deviations["nishiyama_wassermann"], KS_TO_NW_DEG, atol=0.01)
        assert_allclose(deviations["greninger_troiano"], KS_TO_GT_DEG, atol=0.01)
        assert_allclose(report.margin_deg, KS_TO_GT_DEG, atol=0.01)

    def test_burgers_is_identified_and_stated_in_four_index_notation(self) -> None:
        parent, child = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 4, 7, 10))
        report = characterize_orientation_relationship(parents, children)
        assert report.best_catalog_name == "burgers"
        assert report.best_catalog_deviation_deg < 1e-4
        assert report.is_conclusive
        assert "(0001)" in report.statement_text().replace(" ", "")

    def test_moderate_scatter_still_identifies_the_relationship(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(
            relationship,
            variant_indices=tuple(range(1, 13)),
            scatter_deg=0.5,
            seed=11,
        )
        report = characterize_orientation_relationship(parents, children)
        assert report.best_catalog_name == "kurdjumov_sachs"
        # The fit averages the noise: the recovered rotation must sit far
        # closer to KS than the per-pair scatter does.
        assert report.best_catalog_deviation_deg < report.mean_residual_deg
        assert report.is_conclusive

    def test_scatter_comparable_to_the_catalog_spacing_is_reported_inconclusive(
        self,
    ) -> None:
        """The failure mode must be an admitted 'don't know', not a wrong name.

        At 5 deg scatter the data cannot separate relationships that are
        themselves only 2.4 deg apart, and the report must say so.
        """

        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(
            relationship,
            variant_indices=tuple(range(1, 13)),
            scatter_deg=5.0,
            seed=11,
        )
        report = characterize_orientation_relationship(parents, children)
        assert not report.is_conclusive
        assert "NOT conclusively identified" in report.describe()

    def test_a_random_rotation_matches_no_catalog_relationship(self) -> None:
        parent, child = _fcc_bcc()
        frame = specimen_frame()
        parent_matrix = PARENT_ROTATION.as_matrix()
        arbitrary = Rotation.from_axis_angle([0.3, 0.9, 0.31], 1.31).as_matrix()
        parents = OrientationSet.from_matrices(
            parent_matrix[None], specimen_frame=frame, phase=parent
        )
        children = OrientationSet.from_matrices(
            (parent_matrix @ arbitrary.T)[None], specimen_frame=frame, phase=child
        )
        report = characterize_orientation_relationship(parents, children)
        assert not report.matches_catalog
        assert not report.is_conclusive
        assert "matches no catalog relationship" in report.describe()

    def test_single_pair_report_states_that_scatter_is_uninformative(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1,))
        report = characterize_orientation_relationship(parents, children)
        assert report.pair_count == 1
        assert "says nothing about measurement quality" in report.describe()

    def test_explicit_catalog_overrides_the_default(self) -> None:
        parent, child = _fcc_bcc()
        bain = OrientationRelationship.from_bain_correspondence(
            parent_phase=parent, child_phase=child
        )
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 2))
        report = characterize_orientation_relationship(parents, children, catalog=(bain,))
        assert report.catalog_names == ("bain",)
        # A single candidate cannot be out-ranked, so the margin is infinite and
        # the verdict rests on the fit alone; KS is 11 deg from Bain, so it fails.
        assert not report.matches_catalog
        assert not report.is_conclusive

    def test_supplying_a_nominal_reproduces_the_seedless_answer(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 6, 11))
        seedless = characterize_orientation_relationship(parents, children)
        seeded = characterize_orientation_relationship(
            parents,
            children,
            nominal=OrientationRelationship.from_pitsch_correspondence(
                parent_phase=parent, child_phase=child
            ),
        )
        assert seeded.best_catalog_name == seedless.best_catalog_name
        assert_allclose(
            seeded.best_catalog_deviation_deg,
            seedless.best_catalog_deviation_deg,
            atol=1e-6,
        )


class TestEulerEntryPoint:
    def test_euler_angles_give_the_same_answer_as_orientation_sets(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 7, 13, 19))
        reference = characterize_orientation_relationship(parents, children)
        report = orientation_relationship_from_euler(
            parents.as_bunge_euler(degrees=True),
            children.as_bunge_euler(degrees=True),
            parent_phase=parent,
            child_phase=child,
        )
        assert report.best_catalog_name == reference.best_catalog_name
        assert_allclose(
            report.best_catalog_deviation_deg,
            reference.best_catalog_deviation_deg,
            atol=1e-6,
        )


class TestReportContract:
    def test_describe_and_json_agree(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 2, 3))
        report = characterize_orientation_relationship(parents, children)
        payload = report.to_json_dict()
        assert payload["schema"] == "pytex.or_characterization_report/1"
        assert payload["best_catalog_name"] == report.best_catalog_name
        assert payload["is_conclusive"] is report.is_conclusive
        assert payload["statement_text"] == report.statement_text()
        assert len(payload["plane_statements"]) == len(report.plane_statements)
        assert payload["best_catalog_name"] in report.describe()

    def test_arrays_are_read_only(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 2))
        report = characterize_orientation_relationship(parents, children)
        with pytest.raises(ValueError):
            report.residuals_deg[0] = 1.0
        with pytest.raises(ValueError):
            report.plane_statements[0].parent_indices[0] = 9


class TestValidation:
    def test_rejects_mismatched_pair_counts(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 2, 3))
        with pytest.raises(ValueError, match="must be paired"):
            characterize_orientation_relationship(parents, children[:2])

    def test_rejects_orientation_sets_in_different_specimen_frames(self) -> None:
        parent, child = _fcc_bcc()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        parents, children = _planted_pairs(relationship, variant_indices=(1, 2))
        from pytex.core import sample_frame

        relabelled = OrientationSet.from_matrices(
            children.as_matrices(), specimen_frame=sample_frame(), phase=child
        )
        with pytest.raises(ValueError, match="share a specimen frame"):
            characterize_orientation_relationship(parents, relabelled)
