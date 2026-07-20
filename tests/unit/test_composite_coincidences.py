"""Regression tests for spot-coincidence analysis and zone sweeps (CD5).

Pinned reference: for the KS composite along parent [0 1 -1] (austenite
a = 3.6 A, martensite a = 2.87 A, camera constant 180 mm*angstrom) the
closest parent/child superpositions are the close-packed-plane pairs
{111}_fcc || {011}_bcc: |g(111)| = sqrt(3)/3.6 = 0.48113 1/A and
|g(011)| = sqrt(2)/2.87 = 0.49276 1/A are parallel for the exactly-oriented
variants, separated by (0.49276 - 0.48113) * 180 = 2.0938 mm on the detector.
Each exactly-oriented variant therefore contributes exactly its two
antipodal close-packed pairs at tolerance 2.5 mm.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.lattice import ZoneAxis
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import (
    CompositeSAEDPattern,
    SpotCoincidence,
    SpotCoincidenceReport,
    find_spot_coincidences,
    simulate_composite_saed,
    sweep_parent_zone_axes,
)
from tests.unit.test_composite_saed import make_bcc_hcp_phases, make_fcc_bcc_phases


@pytest.fixture(scope="module")
def ks() -> OrientationRelationship:
    parent, child = make_fcc_bcc_phases()
    return OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )


@pytest.fixture(scope="module")
def ks_composite(ks: OrientationRelationship) -> CompositeSAEDPattern:
    zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
    return simulate_composite_saed(ks, zone)


class TestFindSpotCoincidences:
    def test_pinned_total_at_2p5_mm(self, ks_composite: CompositeSAEDPattern) -> None:
        report = find_spot_coincidences(ks_composite, tolerance_mm=2.5)
        assert len(report.coincidences) == 24

    def test_exact_variants_have_two_close_packed_pairs(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        exact = [
            pattern.variant_index
            for pattern in ks_composite.variant_patterns
            if pattern.nearest_zone_axis.deviation_deg < 1e-6
        ]
        assert exact
        report = find_spot_coincidences(ks_composite, tolerance_mm=2.5)
        for variant_index in exact:
            pairs = report.coincidences_for_variant(variant_index)
            assert len(pairs) == 2
            for pair in pairs:
                # {111}_p || {011}_c close-packed superposition at 2.0938 mm
                assert tuple(sorted(np.abs(pair.parent_hkl))) == (1, 1, 1)
                assert tuple(sorted(np.abs(pair.child_hkl))) == (0, 1, 1)
                assert pair.separation_mm == pytest.approx(2.09378, abs=2e-4)

    def test_separation_derivation_close_packed_pair(self) -> None:
        g_parent = np.sqrt(3.0) / 3.6
        g_child = np.sqrt(2.0) / 2.87
        assert (g_child - g_parent) * 180.0 == pytest.approx(2.09378, abs=2e-4)

    def test_brute_force_parity(self, ks_composite: CompositeSAEDPattern) -> None:
        tolerance = 3.7
        report = find_spot_coincidences(ks_composite, tolerance_mm=tolerance)
        found = {
            (
                item.variant_index,
                tuple(int(v) for v in item.parent_hkl),
                tuple(int(v) for v in item.child_hkl),
            )
            for item in report.coincidences
        }
        assert ks_composite.parent_spots is not None
        parent_table = ks_composite.parent_spots
        expected = set()
        for pattern in ks_composite.variant_patterns:
            child_table = pattern.spots
            if not len(child_table):
                continue
            deltas = (
                parent_table.detector_mm[:, None, :] - child_table.detector_mm[None, :, :]
            )
            distances = np.linalg.norm(deltas, axis=2)
            rows, cols = np.nonzero(distances <= tolerance)
            for parent_row, child_row in zip(rows, cols, strict=True):
                expected.add(
                    (
                        pattern.variant_index,
                        tuple(int(v) for v in parent_table.hkl[parent_row]),
                        tuple(int(v) for v in child_table.hkl[child_row]),
                    )
                )
        assert found == expected

    def test_tolerance_monotonicity(self, ks_composite: CompositeSAEDPattern) -> None:
        counts = [
            len(find_spot_coincidences(ks_composite, tolerance_mm=tol).coincidences)
            for tol in (1.0, 2.5, 5.0)
        ]
        assert counts[0] == 0
        assert counts[0] <= counts[1] <= counts[2]

    def test_sorted_by_separation(self, ks_composite: CompositeSAEDPattern) -> None:
        report = find_spot_coincidences(ks_composite, tolerance_mm=5.0)
        separations = [item.separation_mm for item in report.coincidences]
        assert separations == sorted(separations)

    def test_all_pairs_respect_tolerance(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        report = find_spot_coincidences(ks_composite, tolerance_mm=2.5)
        assert all(item.separation_mm <= 2.5 for item in report.coincidences)

    def test_requires_parent_spots(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(
            ks, zone, include_parent=False, variant_indices=(1,)
        )
        with pytest.raises(ValueError, match="include_parent"):
            find_spot_coincidences(composite)

    def test_invalid_tolerance_raises(self, ks_composite: CompositeSAEDPattern) -> None:
        with pytest.raises(ValueError, match="tolerance_mm"):
            find_spot_coincidences(ks_composite, tolerance_mm=0.0)

    def test_describe_content(self, ks_composite: CompositeSAEDPattern) -> None:
        report = find_spot_coincidences(ks_composite, tolerance_mm=2.5)
        text = report.describe()
        assert "kurdjumov_sachs" in text
        assert "[0 1 -1]" in text
        assert "2.5 mm" in text
        assert "24" in text
        assert "V2:" in text
        assert "signature" in text

    def test_variant_counts_cover_all_variants(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        report = find_spot_coincidences(ks_composite, tolerance_mm=2.5)
        assert report.variant_indices == ks_composite.variant_indices
        counts = dict(report.variant_spot_counts)
        for pattern in ks_composite.variant_patterns:
            assert counts[pattern.variant_index] == len(pattern.spots)


class TestBurgersCoincidences:
    """Burgers beta->alpha: the canonical hexagonal coincidence case.

    The Burgers plane parallelism {110}_bcc || (0001)_hcp shows up in
    reciprocal space as a near-exact superposition of the {110}_bcc and
    (0002)_hcp reflections, because d(110)_bcc = a/sqrt(2) = 2.3381 A is
    almost exactly d(0002)_hcp = c/2 = 2.3428 A for beta/alpha titanium.
    The residual detector separation is therefore
    (sqrt(2)/a_bcc - 2/c_hcp) * camera_constant = 0.1545 mm, i.e. the two
    reflections sit well inside a typical spot diameter — which is why the
    Burgers OR reads as a decorated single pattern in the basal view.
    """

    @pytest.fixture(scope="class")
    def burgers_composite(self) -> CompositeSAEDPattern:
        beta, alpha = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=beta, child_phase=alpha
        )
        zone = ZoneAxis(np.array([1, 1, 0]), phase=beta)
        return simulate_composite_saed(relationship, zone)

    def test_closest_coincidence_is_110_bcc_on_0002_hcp(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        report = find_spot_coincidences(burgers_composite, tolerance_mm=1.0)
        assert report.coincidences
        closest = report.coincidences[0]
        assert tuple(sorted(np.abs(closest.parent_hkl))) == (0, 1, 1)
        assert tuple(np.abs(closest.child_hkl)) == (0, 0, 2)

    def test_pinned_analytic_separation(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        expected_mm = (np.sqrt(2.0) / 3.3065 - 2.0 / 4.6855) * 180.0
        assert expected_mm == pytest.approx(0.15450, abs=1e-4)
        report = find_spot_coincidences(burgers_composite, tolerance_mm=1.0)
        assert report.coincidences[0].separation_mm == pytest.approx(expected_mm, abs=1e-6)

    def test_hexagonal_child_labels_are_four_index(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        report = find_spot_coincidences(burgers_composite, tolerance_mm=1.0)
        label = report.coincidences[0].label()
        assert "(0 0 0 2)" in label or "(0 0 0 -2)" in label
        # the cubic parent keeps three-index notation in the same label
        parent_token = label.split("_p")[0]
        assert parent_token.count(" ") == 2

    def test_tolerance_monotonicity(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        counts = [
            len(find_spot_coincidences(burgers_composite, tolerance_mm=tol).coincidences)
            for tol in (1.0, 2.0, 4.0)
        ]
        assert counts[0] < counts[1] < counts[2]

    def test_describe_mentions_burgers(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        report = find_spot_coincidences(burgers_composite, tolerance_mm=1.0)
        text = report.describe()
        assert "burgers" in text
        assert "[1 1 0]" in text


class TestReportValidation:
    def test_coincidence_validation(self) -> None:
        with pytest.raises(ValueError, match="variant_index"):
            SpotCoincidence(
                variant_index=0,
                parent_hkl=np.array([1, 1, 1]),
                child_hkl=np.array([0, 1, 1]),
                parent_detector_mm=np.zeros(2),
                child_detector_mm=np.zeros(2),
                separation_mm=0.5,
            )
        with pytest.raises(ValueError, match="separation_mm"):
            SpotCoincidence(
                variant_index=1,
                parent_hkl=np.array([1, 1, 1]),
                child_hkl=np.array([0, 1, 1]),
                parent_detector_mm=np.zeros(2),
                child_detector_mm=np.zeros(2),
                separation_mm=-1.0,
            )

    def test_report_rejects_pair_beyond_tolerance(self) -> None:
        pair = SpotCoincidence(
            variant_index=1,
            parent_hkl=np.array([1, 1, 1]),
            child_hkl=np.array([0, 1, 1]),
            parent_detector_mm=np.zeros(2),
            child_detector_mm=np.array([3.0, 0.0]),
            separation_mm=3.0,
        )
        with pytest.raises(ValueError, match="tolerance"):
            SpotCoincidenceReport(
                relationship_name="ks",
                parent_zone_label="[0 1 -1]",
                tolerance_mm=1.0,
                parent_spot_count=10,
                variant_spot_counts=((1, 5),),
                coincidences=(pair,),
            )

    def test_report_rejects_unknown_variant_reference(self) -> None:
        pair = SpotCoincidence(
            variant_index=2,
            parent_hkl=np.array([1, 1, 1]),
            child_hkl=np.array([0, 1, 1]),
            parent_detector_mm=np.zeros(2),
            child_detector_mm=np.zeros(2),
            separation_mm=0.0,
        )
        with pytest.raises(ValueError, match="missing"):
            SpotCoincidenceReport(
                relationship_name="ks",
                parent_zone_label="[0 1 -1]",
                tolerance_mm=1.0,
                parent_spot_count=10,
                variant_spot_counts=((1, 5),),
                coincidences=(pair,),
            )


class TestZoneSweep:
    def test_sweep_yields_composites_in_order(self, ks: OrientationRelationship) -> None:
        zones = [
            ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase),
            ZoneAxis(np.array([0, 0, 1]), phase=ks.parent_phase),
        ]
        composites = list(
            sweep_parent_zone_axes(ks, zones, variant_indices=(1, 2))
        )
        assert len(composites) == 2
        for composite, zone in zip(composites, zones, strict=True):
            assert np.array_equal(composite.parent_zone_axis.indices, zone.indices)
            assert composite.variant_indices == (1, 2)

    def test_sweep_is_lazy(self, ks: OrientationRelationship) -> None:
        zones = iter(
            [
                ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase),
                ZoneAxis(np.array([1, 1, 1]), phase=ks.parent_phase),
            ]
        )
        iterator = sweep_parent_zone_axes(ks, zones, variant_indices=())
        first = next(iterator)
        assert np.array_equal(first.parent_zone_axis.indices, [0, 1, -1])
