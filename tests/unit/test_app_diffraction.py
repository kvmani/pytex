"""The composite-SAED service, checked against what the relationship asserts.

The load-bearing assertions are structural rather than numerical: a
Kurdjumov-Sachs composite pattern must contain spots from twenty-four distinct
variants and no others, every spot must carry the row that its hover card and
its CSV row are both built from, and the detector radius of a spot must equal
the camera constant divided by its d-spacing — which is the relation the whole
measurement rests on.
"""

from __future__ import annotations

import math

import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError

AUSTENITE = {"builtin": "austenite_fcc"}
FERRITE = {"builtin": "fe_bcc"}
BETA_ZIRCONIUM = {"builtin": "zr_bcc_beta"}


def simulate(**request: object) -> dict:
    payload = {"phase": AUSTENITE, "child_phase": FERRITE, **request}
    return REGISTRY.call("diffraction.composite_saed", payload)


class TestCompositePattern:
    def test_all_kurdjumov_sachs_variants_appear(self) -> None:
        result = simulate(zone_axis=[0, 0, 1])
        variants = {
            spot["variant"] for spot in result["data"]["spots"] if spot["source"] == "variant"
        }
        assert variants == set(range(1, 25))

    def test_nishiyama_wassermann_gives_twelve(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], relationship="nishiyama_wassermann")
        variants = {
            spot["variant"] for spot in result["data"]["spots"] if spot["source"] == "variant"
        }
        assert variants == set(range(1, 13))

    def test_a_variant_selection_is_honoured(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="2 5")
        variants = {
            spot["variant"] for spot in result["data"]["spots"] if spot["source"] == "variant"
        }
        assert variants == {2, 5}

    def test_a_variant_range_is_honoured(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1-3")
        variants = {
            spot["variant"] for spot in result["data"]["spots"] if spot["source"] == "variant"
        }
        assert variants == {1, 2, 3}

    def test_a_variant_outside_the_range_is_refused(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            simulate(zone_axis=[0, 0, 1], variants="99")
        assert excinfo.value.details["field"] == "variants"

    def test_nonsense_in_the_variant_field_is_refused(self) -> None:
        with pytest.raises(InvalidInputError, match="variant number"):
            simulate(zone_axis=[0, 0, 1], variants="third")

    def test_the_parent_can_be_excluded(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1", include_parent=False)
        assert all(spot["source"] == "variant" for spot in result["data"]["spots"])

    def test_the_parent_pattern_is_the_fcc_zone(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1")
        parent = [spot for spot in result["data"]["spots"] if spot["source"] == "parent"]
        for spot in parent:
            # Down [001] every visible reflection has l = 0 and unmixed parity.
            assert spot["l"] == 0
            parities = {spot["h"] % 2, spot["k"] % 2, spot["l"] % 2}
            assert len(parities) == 1


class TestSpotGeometry:
    def test_detector_radius_follows_the_camera_constant(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1")
        camera = result["data"]["camera_constant_mm_angstrom"]
        for spot in result["data"]["spots"]:
            radius = math.hypot(spot["detector_x_mm"], spot["detector_y_mm"])
            assert radius == pytest.approx(spot["detector_r_mm"], rel=1e-9)
            # r = (camera constant) / d only for a reflection exactly on the
            # Ewald sphere: the radius is set by the in-plane component of g, so
            # a non-zero excitation error puts the spot marginally inside it.
            assert radius <= camera / spot["d_angstrom"] + 1e-9

    def test_g_is_the_reciprocal_of_d(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1")
        for spot in result["data"]["spots"]:
            assert spot["g_inv_angstrom"] == pytest.approx(1.0 / spot["d_angstrom"], rel=1e-9)

    def test_a_larger_camera_constant_scales_the_pattern(self) -> None:
        near = simulate(zone_axis=[0, 0, 1], variants="1")
        far = simulate(zone_axis=[0, 0, 1], variants="1", camera_constant_mm_angstrom=360.0)
        assert far["data"]["detector_radius_mm"] == pytest.approx(
            2.0 * near["data"]["detector_radius_mm"], rel=1e-9
        )

    def test_the_reported_radius_bounds_every_spot(self) -> None:
        result = simulate(zone_axis=[0, 0, 1])
        radius = result["data"]["detector_radius_mm"]
        for spot in result["data"]["spots"]:
            assert math.hypot(spot["detector_x_mm"], spot["detector_y_mm"]) <= radius + 1e-6


class TestSpotRows:
    def test_every_spot_carries_every_declared_column(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1")
        keys = {column["key"] for column in result["data"]["columns"]}
        for spot in result["data"]["spots"]:
            assert keys <= set(spot)

    def test_the_hover_columns_are_the_table_columns(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1")
        assert result["data"]["columns"] == result["table"]["columns"]

    def test_a_variant_spot_names_its_variant(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="5")
        variant_spots = [spot for spot in result["data"]["spots"] if spot["source"] == "variant"]
        assert variant_spots
        assert all(spot["origin"] == "Variant 5" for spot in variant_spots)

    def test_a_variant_spot_reports_a_rationalized_zone_axis(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1")
        variant_spot = next(spot for spot in result["data"]["spots"] if spot["source"] == "variant")
        # An exactly rational parent axis maps to an irrational child axis, so
        # the label must carry the deviation rather than pretend it is exact.
        assert "off)" in variant_spot["zone_axis_indexed"]

    def test_sources_are_grouped_per_variant(self) -> None:
        result = simulate(zone_axis=[0, 0, 1])
        labels = [source["label"] for source in result["data"]["sources"]]
        assert labels[0] == "Parent"
        assert len(labels) == 25
        assert sum(len(source["spots"]) for source in result["data"]["sources"]) == len(
            result["data"]["spots"]
        )


class TestLimits:
    def test_an_impossible_threshold_is_reported(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            simulate(zone_axis=[0, 0, 1], min_intensity=1.5)
        assert excinfo.value.details["field"] == "min_intensity"

    def test_raising_the_threshold_removes_weak_spots(self) -> None:
        # Rock salt, not iron: the intensity model uses the atomic number with no
        # angular dependence, so a monatomic phase gives every allowed
        # reflection the same intensity and the threshold has nothing to cut.
        base = {
            "phase": {"builtin": "nacl"},
            "child_phase": {"builtin": "nacl", "name": "Halite (product)"},
        }
        many = REGISTRY.call(
            "diffraction.composite_saed",
            {**base, "relationship": "bain", "zone_axis": [0, 1, 1], "min_intensity": 0.0},
        )
        few = REGISTRY.call(
            "diffraction.composite_saed",
            {**base, "relationship": "bain", "zone_axis": [0, 1, 1], "min_intensity": 0.5},
        )
        assert len(few["data"]["spots"]) < len(many["data"]["spots"])

    def test_a_monatomic_phase_says_its_intensities_are_flat(self) -> None:
        result = simulate(zone_axis=[0, 0, 1], variants="1")
        assert any("monatomic" in note for note in result["notes"])
        intensities = {spot["relative_intensity"] for spot in result["data"]["spots"]}
        assert intensities == {1.0}

    def test_an_inapplicable_relationship_explains_itself(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            REGISTRY.call(
                "diffraction.composite_saed",
                {
                    "phase": {"builtin": "zr_hcp"},
                    "child_phase": {"builtin": "zr_hcp"},
                    "relationship": "kurdjumov_sachs",
                    "zone_axis": [0, 0, 1],
                },
            )
        assert excinfo.value.details["field"] == "relationship"


class TestBurgers:
    def test_a_hexagonal_child_gives_twelve_variants(self) -> None:
        result = REGISTRY.call(
            "diffraction.composite_saed",
            {
                "phase": BETA_ZIRCONIUM,
                "child_phase": {"builtin": "zr_hcp"},
                "relationship": "burgers",
                "zone_axis": [1, 1, 1],
            },
        )
        variants = {
            spot["variant"] for spot in result["data"]["spots"] if spot["source"] == "variant"
        }
        assert variants == set(range(1, 13))
