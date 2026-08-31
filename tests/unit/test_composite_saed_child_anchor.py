"""TX4: composite patterns anchored on a product-variant zone axis.

The central expectation is an identity, not a measurement: anchoring on variant
``k``'s image of a parent zone must reproduce the parent-anchored pattern for
that zone exactly, because both routes build the same detector basis about the
same parent direction. Anything else is a defect in one of the two paths.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import CrystalDirection, MillerIndex, OrientationRelationship, ZoneAxis
from pytex.diffraction import (
    KinematicSimulationConfig,
    composite_reflection_table,
    simulate_composite_saed,
    simulate_composite_saed_from_child_zone,
    simulate_zone_axis_spots,
)
from tests.unit.test_composite_saed import make_bcc_hcp_phases, make_fcc_bcc_phases


@pytest.fixture(scope="module")
def burgers() -> OrientationRelationship:
    parent, child = make_bcc_hcp_phases()
    return OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent, child_phase=child
    )


class TestConsistencyIdentity:
    """Both anchoring routes must produce the same detector geometry."""

    @pytest.mark.parametrize("anchor", [1, 2, 3, 4])
    def test_anchoring_on_a_variants_own_child_zone_reproduces_the_parent_pattern(
        self, burgers: OrientationRelationship, anchor: int
    ) -> None:
        """Definitional: R_k^T (R_k z_p) = z_p, so the geometry is unchanged.

        Both entry points build the shared basis through the same
        `zone_basis_from_axis` call about the same parent direction, so every
        sub-pattern must land on identical detector coordinates. The residual
        is the floating-point round trip through the rotation, not a tolerance.
        """

        parent_zone = ZoneAxis(np.array([1, 1, 0]), phase=burgers.parent_phase)
        selection = (1, 2, 3, 4)
        reference = simulate_composite_saed(
            burgers, parent_zone, variant_indices=selection
        )
        recovered = simulate_composite_saed_from_child_zone(
            burgers,
            reference.variant_pattern(anchor).zone_axis_child,
            anchor_variant_index=anchor,
            variant_indices=selection,
        )
        assert reference.parent_spots is not None
        assert recovered.parent_spots is not None
        assert np.array_equal(reference.parent_spots.hkl, recovered.parent_spots.hkl)
        assert_allclose(
            recovered.parent_spots.detector_mm,
            reference.parent_spots.detector_mm,
            atol=1e-10,
        )
        for index in selection:
            expected = reference.variant_pattern(index).spots
            actual = recovered.variant_pattern(index).spots
            assert np.array_equal(expected.hkl, actual.hkl)
            assert_allclose(actual.detector_mm, expected.detector_mm, atol=1e-10)

    def test_the_recovered_parent_direction_is_the_original_zone_axis(
        self, burgers: OrientationRelationship
    ) -> None:
        parent_zone = ZoneAxis(np.array([1, 1, 0]), phase=burgers.parent_phase)
        reference = simulate_composite_saed(burgers, parent_zone, variant_indices=(1,))
        recovered = simulate_composite_saed_from_child_zone(
            burgers,
            reference.variant_pattern(1).zone_axis_child,
            anchor_variant_index=1,
            variant_indices=(1,),
        )
        assert_allclose(
            recovered.parent_zone_axis.unit_vector, parent_zone.unit_vector, atol=1e-12
        )
        assert recovered.nearest_parent_zone_axis is not None
        assert recovered.nearest_parent_zone_axis.deviation_deg < 1e-9


class TestSpotOrderStability:
    """The deterministic sort must survive floating-point ties.

    Symmetry-equivalent reflections have mathematically equal intensity and
    radius that differ in the last few ULPs depending on how the detector basis
    was built. Before the sort keys were quantized, that noise decided the row
    order, so the two anchoring routes produced correctly-positioned but
    *permuted* spot tables — and any exported table or pinned figure inherited
    the permutation.
    """

    def test_equivalent_bases_give_identical_row_order(
        self, burgers: OrientationRelationship
    ) -> None:
        parent = burgers.parent_phase
        exact = ZoneAxis(np.array([1, 1, 0]), phase=parent)
        # The same direction reached with a perturbation far below any
        # crystallographic meaning but far above the sort quantum's target.
        perturbed = CrystalDirection(
            np.array([1.0, 1.0, 0.0]) + np.array([1e-15, -1e-15, 1e-15]), phase=parent
        )
        config = KinematicSimulationConfig()
        first = simulate_zone_axis_spots(parent, exact, config=config)
        second = simulate_zone_axis_spots(parent, perturbed, config=config)
        assert np.array_equal(first.hkl, second.hkl)

    def test_repeated_simulation_is_bitwise_stable(
        self, burgers: OrientationRelationship
    ) -> None:
        parent_zone = ZoneAxis(np.array([1, 1, 1]), phase=burgers.parent_phase)
        first = simulate_composite_saed(burgers, parent_zone, variant_indices=(1, 2))
        second = simulate_composite_saed(burgers, parent_zone, variant_indices=(1, 2))
        assert first.parent_spots is not None
        assert second.parent_spots is not None
        assert np.array_equal(first.parent_spots.hkl, second.parent_spots.hkl)


class TestChildAnchoredGeometry:
    def test_anchoring_on_the_basal_zone_gives_the_expected_parent_direction(
        self, burgers: OrientationRelationship
    ) -> None:
        """Burgers pairs (0001) alpha with a {110} beta plane.

        Anchoring on the alpha basal zone must therefore recover a parent
        direction parallel to a <110> beta direction, exactly — that pairing is
        the definition of the relationship.
        """

        child_zone = ZoneAxis(np.array([0, 0, 1]), phase=burgers.child_phase)
        pattern = simulate_composite_saed_from_child_zone(
            burgers, child_zone, anchor_variant_index=1, variant_indices=(1,)
        )
        assert pattern.nearest_parent_zone_axis is not None
        assert pattern.nearest_parent_zone_axis.deviation_deg < 1e-9
        indices = pattern.nearest_parent_zone_axis.indices
        assert sorted(np.abs(indices).tolist()) == [0, 1, 1]

    def test_the_anchor_variants_own_zone_axis_is_exactly_rational(
        self, burgers: OrientationRelationship
    ) -> None:
        """The variant the geometry was anchored on sees exactly what was asked for."""

        child_zone = ZoneAxis(np.array([0, 0, 1]), phase=burgers.child_phase)
        pattern = simulate_composite_saed_from_child_zone(
            burgers, child_zone, anchor_variant_index=2, variant_indices=(2,)
        )
        anchored = pattern.variant_pattern(2)
        assert anchored.nearest_zone_axis.deviation_deg < 1e-9
        assert_allclose(
            anchored.zone_axis_child.unit_vector, child_zone.unit_vector, atol=1e-12
        )

    def test_the_anchor_is_recorded_and_described(
        self, burgers: OrientationRelationship
    ) -> None:
        pattern = simulate_composite_saed_from_child_zone(
            burgers,
            ZoneAxis(np.array([0, 0, 1]), phase=burgers.child_phase),
            anchor_variant_index=3,
            variant_indices=(1, 3),
        )
        assert pattern.anchor_variant_index == 3
        assert "variant 3 of child" in pattern.anchor_description()
        text = pattern.describe()
        assert "anchored on the zone axis of variant 3" in text
        assert "nearest" in pattern.parent_zone_axis_label()

    def test_a_parent_anchored_pattern_records_no_anchor_variant(
        self, burgers: OrientationRelationship
    ) -> None:
        pattern = simulate_composite_saed(
            burgers,
            ZoneAxis(np.array([1, 1, 0]), phase=burgers.parent_phase),
            variant_indices=(1,),
        )
        assert pattern.anchor_variant_index is None
        assert pattern.nearest_parent_zone_axis is None
        assert pattern.parent_zone_axis_label() == "[110]"


class TestInPlaneAlignment:
    def test_aligning_a_child_reflection_places_it_along_plus_u(self) -> None:
        """`align_child_g` must act in the child's own indices.

        The requested child reflection is mapped to the parent frame before the
        basis is built, so the caller never has to convert it by hand. Placing
        it along +u means its detector coordinate has a positive first
        component and a vanishing second.
        """

        parent, child = make_fcc_bcc_phases()
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent, child_phase=child
        )
        target = MillerIndex(np.array([1, 1, 0]), phase=child)
        pattern = simulate_composite_saed_from_child_zone(
            relationship,
            ZoneAxis(np.array([0, 0, 1]), phase=child),
            anchor_variant_index=1,
            variant_indices=(1,),
            align_child_g=target,
        )
        spots = pattern.variant_pattern(1).spots
        matches = [
            position
            for position in range(len(spots))
            if tuple(int(value) for value in spots.hkl[position]) == (1, 1, 0)
        ]
        assert matches, "the aligned reflection must be present in the pattern"
        coordinates = spots.detector_mm[matches[0]]
        assert coordinates[0] > 0.0
        assert abs(coordinates[1]) < 1e-9


class TestValidation:
    def test_rejects_a_child_zone_axis_of_the_wrong_phase(
        self, burgers: OrientationRelationship
    ) -> None:
        with pytest.raises(ValueError, match=r"child_zone_axis\.phase must match"):
            simulate_composite_saed_from_child_zone(
                burgers, ZoneAxis(np.array([1, 1, 0]), phase=burgers.parent_phase)
            )

    def test_rejects_an_unknown_anchor_variant(
        self, burgers: OrientationRelationship
    ) -> None:
        with pytest.raises(ValueError, match="Unknown anchor_variant_index"):
            simulate_composite_saed_from_child_zone(
                burgers,
                ZoneAxis(np.array([0, 0, 1]), phase=burgers.child_phase),
                anchor_variant_index=99,
            )

    def test_rejects_two_conflicting_alignment_references(
        self, burgers: OrientationRelationship
    ) -> None:
        with pytest.raises(ValueError, match="at most one of align_parent_g"):
            simulate_composite_saed(
                burgers,
                ZoneAxis(np.array([1, 1, 0]), phase=burgers.parent_phase),
                align_parent_g=MillerIndex(np.array([0, 0, 2]), phase=burgers.parent_phase),
                align_g_cartesian=np.array([1.0, 0.0, 0.0]),
            )

    def test_rejects_a_non_positive_anchor_index_on_the_pattern(
        self, burgers: OrientationRelationship
    ) -> None:
        pattern = simulate_composite_saed(
            burgers,
            ZoneAxis(np.array([1, 1, 0]), phase=burgers.parent_phase),
            variant_indices=(1,),
        )
        from dataclasses import replace

        with pytest.raises(ValueError, match="anchor_variant_index must be a positive"):
            replace(pattern, anchor_variant_index=0)


class TestExportCarriesTheAnchor:
    def test_the_reflection_table_label_reports_the_irrational_parent_zone(
        self, burgers: OrientationRelationship
    ) -> None:
        pattern = simulate_composite_saed_from_child_zone(
            burgers,
            ZoneAxis(np.array([1, 1, 0]), phase=burgers.child_phase),
            anchor_variant_index=1,
            variant_indices=(1,),
        )
        table = composite_reflection_table(pattern)
        assert "nearest" in table.parent_zone_axis_label

    def test_the_manifest_records_the_anchor_and_the_nearest_parent_zone(
        self, burgers: OrientationRelationship, tmp_path
    ) -> None:
        import json

        import jsonschema

        from pytex.adapters import composite_saed_manifest_schema_path
        from pytex.diffraction import export_composite_saed
        pattern = simulate_composite_saed_from_child_zone(
            burgers,
            ZoneAxis(np.array([0, 0, 1]), phase=burgers.child_phase),
            anchor_variant_index=2,
            variant_indices=(1, 2),
        )
        export = export_composite_saed(pattern, tmp_path, figure_formats=())
        payload = json.loads(export.manifest_path.read_text(encoding="utf-8"))
        jsonschema.validate(
            payload,
            json.loads(composite_saed_manifest_schema_path().read_text(encoding="utf-8")),
        )
        assert payload["anchor_variant_index"] == 2
        assert payload["parent_zone_axis_nearest"] is not None
        assert payload["parent_zone_axis_nearest"]["deviation_deg"] < 1e-6
