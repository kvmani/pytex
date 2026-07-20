"""Regression tests for composite OR SAED assembly (CD2).

Pinned references:

- Kurdjumov-Sachs parallelism <-101>_fcc || <-1-11>_bcc: a parent [0 1 -1]
  zone axis must map exactly onto a <111>-type child zone for at least one
  variant (Kurdjumov & Sachs 1930; Morito et al. 2003 variant convention).
- Nishiyama-Wassermann parallelism [1 -1 0]_fcc || [1 0 0]_bcc: a parent
  [1 -1 0] zone maps exactly onto a <100>-type child zone.
- Bain correspondence (001)_fcc || (001)_bcc with [110]_fcc || [100]_bcc:
  along the common [001] zone the child (200) reflection is collinear with
  the parent (220) reflection (the classic 45-degree in-plane relation).
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import (
    AtomicSite,
    CrystalDirection,
    Lattice,
    MillerIndex,
    Phase,
    UnitCell,
    ZoneAxis,
)
from pytex.core.symmetry import SymmetrySpec
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import (
    CompositeSAEDPattern,
    rationalize_zone_axis,
    simulate_composite_saed,
)
from pytex.diffraction.kinematic import KinematicSimulationConfig


def make_fcc_bcc_phases() -> tuple[Phase, Phase]:
    parent_frame = ReferenceFrame(
        "parent_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
    )
    child_frame = ReferenceFrame(
        "child_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
    )
    parent_lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame)
    child_lattice = Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame)
    parent_sites = tuple(
        AtomicSite(label=f"Fe{i}", species="Fe", fractional_coordinates=np.array(coords))
        for i, coords in enumerate(
            [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
        )
    )
    child_sites = (
        AtomicSite(label="Fe1", species="Fe", fractional_coordinates=np.array([0.0, 0.0, 0.0])),
        AtomicSite(label="Fe2", species="Fe", fractional_coordinates=np.array([0.5, 0.5, 0.5])),
    )
    parent = Phase(
        "austenite",
        lattice=parent_lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
        crystal_frame=parent_frame,
        unit_cell=UnitCell(lattice=parent_lattice, sites=parent_sites),
        space_group_symbol="Fm-3m",
    )
    child = Phase(
        "martensite",
        lattice=child_lattice,
        symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
        crystal_frame=child_frame,
        unit_cell=UnitCell(lattice=child_lattice, sites=child_sites),
        space_group_symbol="Im-3m",
    )
    return parent, child


@pytest.fixture(scope="module")
def fcc_bcc() -> tuple[Phase, Phase]:
    return make_fcc_bcc_phases()


@pytest.fixture(scope="module")
def ks(fcc_bcc: tuple[Phase, Phase]) -> OrientationRelationship:
    parent, child = fcc_bcc
    return OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )


class TestRationalizeZoneAxis:
    def test_exact_rational_direction_recovered(self, fcc_bcc: tuple[Phase, Phase]) -> None:
        parent, _ = fcc_bcc
        direction = CrystalDirection(np.array([1.0, 1.0, -1.0]), phase=parent)
        result = rationalize_zone_axis(direction)
        assert tuple(result.indices) == (1, 1, -1)
        assert result.deviation_deg == pytest.approx(0.0, abs=1e-9)

    def test_perturbed_direction_snaps_with_small_deviation(
        self, fcc_bcc: tuple[Phase, Phase]
    ) -> None:
        parent, _ = fcc_bcc
        direction = CrystalDirection(np.array([0.0, 1.0, 1.02]), phase=parent)
        result = rationalize_zone_axis(direction, max_index=3)
        assert tuple(result.indices) == (0, 1, 1)
        assert 0.0 < result.deviation_deg < 1.0

    def test_sign_sensitive(self, fcc_bcc: tuple[Phase, Phase]) -> None:
        parent, _ = fcc_bcc
        direction = CrystalDirection(np.array([0.0, 0.0, -1.0]), phase=parent)
        result = rationalize_zone_axis(direction)
        assert tuple(result.indices) == (0, 0, -1)

    def test_invalid_max_index_raises(self, fcc_bcc: tuple[Phase, Phase]) -> None:
        parent, _ = fcc_bcc
        direction = CrystalDirection(np.array([0.0, 0.0, 1.0]), phase=parent)
        with pytest.raises(ValueError, match="max_index"):
            rationalize_zone_axis(direction, max_index=0)


class TestKurdjumovSachsComposite:
    def test_parent_011_zone_maps_exactly_to_child_111_for_some_variant(
        self, ks: OrientationRelationship
    ) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(ks, zone)
        assert len(composite.variant_patterns) == 24
        exact = [
            pattern
            for pattern in composite.variant_patterns
            if pattern.nearest_zone_axis.deviation_deg < 1e-6
        ]
        assert exact, "KS <-101>_p || <-1-11>_c demands an exact <111> child zone"
        for pattern in exact:
            assert tuple(sorted(np.abs(pattern.nearest_zone_axis.indices))) == (1, 1, 1)

    def test_child_patterns_share_parent_detector_basis(
        self, ks: OrientationRelationship
    ) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(ks, zone, variant_indices=(1, 2))
        for pattern in composite.variant_patterns:
            rotation = pattern.variant.parent_to_child_rotation.as_matrix()
            expected_basis = rotation @ composite.zone_basis_parent
            assert np.allclose(pattern.spots.basis, expected_basis, atol=1e-12)

    def test_variant_subset_and_order(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(ks, zone, variant_indices=(5, 1, 24))
        assert composite.variant_indices == (5, 1, 24)

    def test_unknown_variant_index_raises(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        with pytest.raises(ValueError, match="Unknown variant"):
            simulate_composite_saed(ks, zone, variant_indices=(25,))

    def test_include_parent_false(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(
            ks, zone, include_parent=False, variant_indices=(1,)
        )
        assert composite.parent_spots is None
        assert composite.spot_count() == len(composite.variant_patterns[0].spots)

    def test_irrational_child_zones_still_simulate(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 0, 1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(ks, zone, variant_indices=(1, 2, 3))
        for pattern in composite.variant_patterns:
            table = pattern.spots
            assert np.all(
                np.abs(table.excitation_error_inv_angstrom)
                <= table.config.max_excitation_error_inv_angstrom + 1e-12
            )

    def test_describe_content(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(ks, zone, variant_indices=(1, 2))
        text = composite.describe()
        assert "kurdjumov_sachs" in text
        assert "austenite" in text
        assert "martensite" in text
        assert "[0 1 -1]" in text
        assert "Variant 1" in text
        assert "excitation-error" in text
        assert "parent-anchored detector basis" in text


class TestNishiyamaWassermannComposite:
    def test_parent_110_zone_maps_exactly_to_child_100(
        self, fcc_bcc: tuple[Phase, Phase]
    ) -> None:
        parent, child = fcc_bcc
        nw = OrientationRelationship.from_nishiyama_wassermann_correspondence(
            parent_phase=parent, child_phase=child
        )
        zone = ZoneAxis(np.array([1, -1, 0]), phase=parent)
        composite = simulate_composite_saed(nw, zone)
        assert len(composite.variant_patterns) == 12
        exact = [
            pattern
            for pattern in composite.variant_patterns
            if pattern.nearest_zone_axis.deviation_deg < 1e-6
        ]
        assert exact, "NW [1 -1 0]_p || [1 0 0]_c demands an exact <100> child zone"
        for pattern in exact:
            assert tuple(sorted(np.abs(pattern.nearest_zone_axis.indices))) == (0, 0, 1)


class TestBainComposite:
    def test_bain_001_zone_45_degree_in_plane_relation(
        self, fcc_bcc: tuple[Phase, Phase]
    ) -> None:
        parent, child = fcc_bcc
        bain = OrientationRelationship.from_bain_correspondence(
            parent_phase=parent, child_phase=child
        )
        zone = ZoneAxis(np.array([0, 0, 1]), phase=parent)
        composite = simulate_composite_saed(bain, zone)
        assert len(composite.variant_patterns) == 3
        exact = [
            pattern
            for pattern in composite.variant_patterns
            if pattern.nearest_zone_axis.deviation_deg < 1e-6
            and tuple(sorted(np.abs(pattern.nearest_zone_axis.indices))) == (0, 0, 1)
        ]
        assert exact, "Bain (001)_p || (001)_c demands an exact [001] child zone"
        pattern = exact[0]
        assert composite.parent_spots is not None
        parent_rows = {
            tuple(int(v) for v in row): i
            for i, row in enumerate(composite.parent_spots.hkl)
        }
        child_rows = {
            tuple(int(v) for v in row): i for i, row in enumerate(pattern.spots.hkl)
        }
        assert (2, 2, 0) in parent_rows
        parent_vec = composite.parent_spots.detector_mm[parent_rows[(2, 2, 0)]]
        candidates = [
            key for key in child_rows if tuple(sorted(np.abs(key))) == (0, 0, 2)
        ]
        assert candidates
        best_cosine = max(
            float(
                np.dot(parent_vec, pattern.spots.detector_mm[child_rows[key]])
                / (
                    np.linalg.norm(parent_vec)
                    * np.linalg.norm(pattern.spots.detector_mm[child_rows[key]])
                )
            )
            for key in candidates
        )
        assert best_cosine == pytest.approx(1.0, abs=1e-9)

    def test_bain_45_degree_between_parent_and_child_200(
        self, fcc_bcc: tuple[Phase, Phase]
    ) -> None:
        parent, child = fcc_bcc
        bain = OrientationRelationship.from_bain_correspondence(
            parent_phase=parent, child_phase=child
        )
        zone = ZoneAxis(np.array([0, 0, 1]), phase=parent)
        composite = simulate_composite_saed(bain, zone)
        pattern = next(
            p
            for p in composite.variant_patterns
            if p.nearest_zone_axis.deviation_deg < 1e-6
            and tuple(sorted(np.abs(p.nearest_zone_axis.indices))) == (0, 0, 1)
        )
        assert composite.parent_spots is not None
        parent_rows = {
            tuple(int(v) for v in row): i
            for i, row in enumerate(composite.parent_spots.hkl)
        }
        child_rows = {
            tuple(int(v) for v in row): i for i, row in enumerate(pattern.spots.hkl)
        }
        parent_vec = composite.parent_spots.detector_mm[parent_rows[(2, 0, 0)]]
        angles = []
        for key in child_rows:
            if tuple(sorted(np.abs(key))) != (0, 0, 2):
                continue
            child_vec = pattern.spots.detector_mm[child_rows[key]]
            cosine = float(
                np.dot(parent_vec, child_vec)
                / (np.linalg.norm(parent_vec) * np.linalg.norm(child_vec))
            )
            angles.append(np.degrees(np.arccos(np.clip(abs(cosine), 0.0, 1.0))))
        assert min(angles) == pytest.approx(45.0, abs=1e-6)


class TestSharedGeometryControls:
    def test_in_plane_rotation_rotates_whole_composite(
        self, ks: OrientationRelationship
    ) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        base = simulate_composite_saed(ks, zone, variant_indices=(1,))
        rotated = simulate_composite_saed(
            ks, zone, variant_indices=(1,), in_plane_rotation_deg=30.0
        )
        angle = np.deg2rad(30.0)
        rotation = np.array(
            [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        )
        for (_, base_table), (_, rotated_table) in zip(
            base.iter_spot_tables(), rotated.iter_spot_tables(), strict=True
        ):
            base_map = dict(zip(base_table.hkl_labels(), base_table.detector_mm, strict=True))
            for label, coords in zip(
                rotated_table.hkl_labels(), rotated_table.detector_mm, strict=True
            ):
                assert np.allclose(coords, rotation @ base_map[label], atol=1e-9)

    def test_align_parent_g_places_reflection_on_positive_u(
        self, ks: OrientationRelationship
    ) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        g_align = MillerIndex(np.array([2, 0, 0]), phase=ks.parent_phase)
        composite = simulate_composite_saed(
            ks, zone, variant_indices=(), align_parent_g=g_align
        )
        assert composite.parent_spots is not None
        rows = {
            tuple(int(v) for v in row): i
            for i, row in enumerate(composite.parent_spots.hkl)
        }
        coords = composite.parent_spots.detector_mm[rows[(2, 0, 0)]]
        assert coords[0] > 0.0
        assert coords[1] == pytest.approx(0.0, abs=1e-9)

    def test_child_config_override(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        child_config = KinematicSimulationConfig(g_max_inv_angstrom=0.6)
        composite = simulate_composite_saed(
            ks, zone, variant_indices=(1,), child_config=child_config
        )
        child_g = np.linalg.norm(composite.variant_patterns[0].spots.g_crystal, axis=1)
        assert np.all(child_g <= 0.6 + 1e-12)
        assert composite.parent_spots is not None
        parent_g = np.linalg.norm(composite.parent_spots.g_crystal, axis=1)
        assert np.max(parent_g) > 0.6

    def test_select_variants_returns_subset(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(ks, zone, variant_indices=(1, 2, 3))
        subset = composite.select_variants([3, 1])
        assert subset.variant_indices == (3, 1)
        assert subset.parent_spots is composite.parent_spots
        with pytest.raises(KeyError):
            composite.select_variants([9])

    def test_all_detector_coordinates_stacks_every_subpattern(
        self, ks: OrientationRelationship
    ) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(ks, zone, variant_indices=(1, 2))
        coords = composite.all_detector_coordinates()
        assert coords.shape == (composite.spot_count(), 2)

    def test_zone_phase_mismatch_raises(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.child_phase)
        with pytest.raises(ValueError, match="parent_zone_axis.phase"):
            simulate_composite_saed(ks, zone)

    def test_align_g_phase_mismatch_raises(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        bad = MillerIndex(np.array([1, 1, 0]), phase=ks.child_phase)
        with pytest.raises(ValueError, match="align_parent_g"):
            simulate_composite_saed(ks, zone, align_parent_g=bad)

    def test_duplicate_variant_patterns_rejected(self, ks: OrientationRelationship) -> None:
        zone = ZoneAxis(np.array([0, 1, -1]), phase=ks.parent_phase)
        composite = simulate_composite_saed(ks, zone, variant_indices=(1,))
        with pytest.raises(ValueError, match="unique"):
            CompositeSAEDPattern(
                relationship=composite.relationship,
                parent_zone_axis=composite.parent_zone_axis,
                parent_spots=composite.parent_spots,
                variant_patterns=composite.variant_patterns * 2,
                zone_basis_parent=composite.zone_basis_parent,
                config=composite.config,
                child_config=composite.child_config,
            )
