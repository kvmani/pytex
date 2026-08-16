"""Scientific and wire-contract tests for the EBSD workbench module.

Each gallery dataset is a construction whose answer is known before the
calculation runs, so these tests check the pipeline against the construction
rather than against a stored output: 40 degrees across the bicrystal boundary,
60 degrees across every twin, twelve grains in the polycrystal, and a KAM inside
a linear gradient that is exactly half the per-step rotation.
"""

from __future__ import annotations

import base64

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.ebsd_gallery import (
    BICRYSTAL_BOUNDARY_DEG,
    GALLERY,
    GRADIENT_DEG_PER_STEP,
    SIGMA3_ANGLE_DEG,
    build_map,
    entry_ids,
    get_entry,
)
from pytex.app.errors import InvalidInputError

BASE_REQUEST: dict[str, object] = {
    "dataset": "equiaxed_polycrystal",
    "colouring": "ipf",
    "ipf_direction": "Z",
    "modulate_by": "none",
    "modulation_floor": 0.25,
    "colour_map": "viridis",
    "show_boundaries": True,
    "grain_threshold_deg": 5.0,
    "high_angle_threshold_deg": 15.0,
    "kam_threshold_deg": 5.0,
    "kam_order": 1,
    # Small, because these tests exercise contracts rather than performance and
    # every one of them pays for a symmetry reduction over the whole map.
    "grid_points": 40,
}


def draw(**overrides: object) -> dict:
    return REGISTRY.call("ebsd.map", {**BASE_REQUEST, **overrides})


def decode_rgb(image: dict) -> np.ndarray:
    raw = np.frombuffer(base64.b64decode(image["data"]), dtype=np.uint8)
    return raw.reshape((image["height"], image["width"], 3))


class TestGalleryConstructions:
    """The datasets are built to have answers; these are the answers."""

    def test_every_entry_builds_a_square_map_with_the_three_channels(self) -> None:
        for entry_id in entry_ids():
            crystal_map = build_map(entry_id, grid=24)
            assert crystal_map.grid_shape == (24, 24)
            assert crystal_map.property_names == ("confidence_index", "fit", "image_quality")

    def test_the_bicrystal_boundary_is_exactly_forty_degrees(self) -> None:
        crystal_map = build_map("bicrystal_gradient", grid=40)
        segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0)
        network = segmentation.boundary_network()
        angles = np.array([segment.misorientation_deg for segment in network.segments])
        assert len(segmentation.grains) == 2
        np.testing.assert_allclose(angles, BICRYSTAL_BOUNDARY_DEG, atol=1e-6)

    def test_every_twin_boundary_is_exactly_sixty_degrees(self) -> None:
        """The coherent twin of a cubic metal, quoted in every textbook."""

        crystal_map = build_map("sigma3_twin", grid=40)
        segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0)
        network = segmentation.boundary_network()
        angles = np.array([segment.misorientation_deg for segment in network.segments])
        assert angles.size > 0
        np.testing.assert_allclose(angles, SIGMA3_ANGLE_DEG, atol=1e-6)

    def test_the_polycrystal_has_the_twelve_grains_it_was_built_from(self) -> None:
        crystal_map = build_map("equiaxed_polycrystal", grid=60)
        segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0)
        assert len(segmentation.grains) == 12

    def test_kam_in_a_linear_gradient_is_half_the_per_step_rotation(self) -> None:
        """The claim the bicrystal entry exists to make.

        Half, not all of it: the four-neighbour kernel averages two neighbours a
        full step along the gradient with two neighbours across it that are
        identical to the centre point. A KAM is an average over a kernel, not a
        gradient magnitude.
        """

        crystal_map = build_map("bicrystal_gradient", grid=40)
        kam = crystal_map.kernel_average_misorientation_deg(threshold_deg=5.0)
        # Well inside the graded grain and away from the map edge, where a point
        # has all four neighbours.
        interior = kam[5:-5, 25:-5]
        np.testing.assert_allclose(interior, GRADIENT_DEG_PER_STEP / 2.0, atol=1e-6)

    def test_grod_is_a_deviation_rather_than_a_ramp(self) -> None:
        """It falls to zero at the grain's reference point and rises either side."""

        crystal_map = build_map("bicrystal_gradient", grid=40)
        segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0)
        grod = segmentation.grod_map_deg()
        graded = grod[20, 20:]
        assert graded.min() == pytest.approx(0.0, abs=1e-6)
        assert graded.max() > 1.0

    def test_the_quality_channels_fall_at_the_boundaries(self) -> None:
        """The modulation feature is pointless unless the channels behave."""

        crystal_map = build_map("equiaxed_polycrystal", grid=60)
        segmentation = crystal_map.segment_grains(max_misorientation_deg=5.0)
        labels = np.asarray(segmentation.labels).reshape(60, 60)
        edge = np.zeros((60, 60), dtype=bool)
        edge[:-1, :] |= labels[:-1, :] != labels[1:, :]
        edge[:, :-1] |= labels[:, :-1] != labels[:, 1:]
        confidence = np.asarray(crystal_map.get_property("confidence_index")).reshape(60, 60)
        fit = np.asarray(crystal_map.get_property("fit")).reshape(60, 60)
        assert confidence[edge].mean() < confidence[~edge].mean()
        assert fit[edge].mean() > fit[~edge].mean()

    def test_an_unknown_entry_names_the_ones_that_exist(self) -> None:
        with pytest.raises(KeyError, match="bicrystal_gradient"):
            get_entry("no_such_dataset")

    def test_every_entry_states_a_known_answer(self) -> None:
        for entry in GALLERY:
            assert entry.known_answer.strip()
            assert entry.teaches.strip()


class TestColourings:
    """The four independent choices, each of which must actually change the map."""

    def test_an_ipf_map_returns_no_scalar_scale(self) -> None:
        result = draw(colouring="ipf")
        assert result["data"]["colour_scale"] is None
        assert result["data"]["colouring"] == "ipf"

    def test_the_three_ipf_directions_give_three_different_maps(self) -> None:
        """One direction does not fix an orientation; three together do."""

        images = {
            direction: decode_rgb(draw(colouring="ipf", ipf_direction=direction)["data"]["image"])
            for direction in ("X", "Y", "Z")
        }
        assert not np.array_equal(images["X"], images["Y"])
        assert not np.array_equal(images["Y"], images["Z"])

    def test_a_scalar_colouring_reports_the_range_its_colours_stand_for(self) -> None:
        """A colour bar without numbers on it is decoration."""

        scale = draw(colouring="grod")["data"]["colour_scale"]
        assert scale["label"] == "GROD"
        assert scale["units"] == "°"
        assert scale["maximum"] >= scale["minimum"]
        assert scale["colour_map"] == "viridis"

    def test_kam_and_grod_are_different_fields(self) -> None:
        kam = decode_rgb(draw(dataset="bicrystal_gradient", colouring="kam")["data"]["image"])
        grod = decode_rgb(draw(dataset="bicrystal_gradient", colouring="grod")["data"]["image"])
        assert not np.array_equal(kam, grod)

    def test_a_grain_map_gives_each_grain_its_own_colour(self) -> None:
        result = draw(colouring="grain")
        image = decode_rgb(result["data"]["image"])
        colours = {tuple(pixel) for pixel in image.reshape(-1, 3)}
        assert len(colours) == result["data"]["grain_count"]

    def test_each_measured_channel_can_be_the_colouring(self) -> None:
        for channel, label in (
            ("confidence_index", "Confidence index"),
            ("fit", "Fit"),
            ("image_quality", "Image quality"),
        ):
            scale = draw(colouring=channel)["data"]["colour_scale"]
            assert scale["label"] == label


class TestModulation:
    """Greyscaling a coloured map by a scalar field."""

    def test_modulation_darkens_without_changing_the_map_shape(self) -> None:
        plain = decode_rgb(draw(modulate_by="none")["data"]["image"]).astype(float)
        greyed = decode_rgb(
            draw(modulate_by="confidence_index", modulation_floor=0.1)["data"]["image"]
        ).astype(float)
        assert greyed.shape == plain.shape
        assert greyed.mean() < plain.mean()
        # Darkening only: no pixel may become brighter than it was.
        assert np.all(greyed <= plain + 1.0)

    def test_the_floor_bounds_how_dark_the_worst_pixel_becomes(self) -> None:
        dark = decode_rgb(
            draw(modulate_by="confidence_index", modulation_floor=0.0)["data"]["image"]
        ).astype(float)
        light = decode_rgb(
            draw(modulate_by="confidence_index", modulation_floor=0.6)["data"]["image"]
        ).astype(float)
        assert light.mean() > dark.mean()

    def test_fit_modulation_is_inverted_because_fit_is_an_error(self) -> None:
        """A large fit is a worse measurement, so it must darken rather than brighten.

        Checked where it matters: the boundary pixels, which the dataset builds
        with the worst fit. Modulating without the inversion would darken
        exactly the pixels that were indexed best.
        """

        crystal_map = build_map("equiaxed_polycrystal", grid=40)
        fit = np.asarray(crystal_map.get_property("fit")).reshape(40, 40)
        worst = fit >= np.quantile(fit, 0.9)
        best = fit <= np.quantile(fit, 0.1)

        plain = decode_rgb(draw(modulate_by="none")["data"]["image"]).astype(float).sum(axis=2)
        greyed = (
            decode_rgb(draw(modulate_by="fit", modulation_floor=0.0)["data"]["image"])
            .astype(float)
            .sum(axis=2)
        )
        # Compare the retained fraction of brightness rather than raw values, so
        # the underlying IPF colour of each region cancels out.
        retained_worst = greyed[worst].sum() / max(plain[worst].sum(), 1.0)
        retained_best = greyed[best].sum() / max(plain[best].sum(), 1.0)
        assert retained_worst < retained_best

    def test_modulation_applies_to_every_colouring(self) -> None:
        """Orthogonality is the whole design; a combination must not be special-cased."""

        for colouring in ("ipf", "grain", "grod", "kam"):
            plain = decode_rgb(draw(colouring=colouring)["data"]["image"]).astype(float)
            greyed = decode_rgb(
                draw(colouring=colouring, modulate_by="confidence_index")["data"]["image"]
            ).astype(float)
            assert greyed.mean() < plain.mean(), colouring


class TestBoundaries:
    """The network drawn over whatever is underneath."""

    def test_boundaries_can_be_superimposed_on_any_colouring(self) -> None:
        for colouring in ("ipf", "grain", "grod", "kam", "confidence_index"):
            result = draw(colouring=colouring, show_boundaries=True)
            assert result["data"]["boundaries"], colouring

    def test_turning_them_off_returns_none_without_changing_the_map(self) -> None:
        with_lines = draw(show_boundaries=True)
        without = draw(show_boundaries=False)
        assert without["data"]["boundaries"] == []
        # The pixels are the map; the boundaries are an overlay on top of it.
        assert with_lines["data"]["image"]["data"] == without["data"]["image"]["data"]

    def test_each_segment_is_one_step_long_and_carries_its_misorientation(self) -> None:
        result = draw()
        step = result["data"]["step_um"]
        for line in result["data"]["boundaries"][:40]:
            length = np.hypot(line["x2"] - line["x1"], line["y2"] - line["y1"])
            assert length == pytest.approx(step, rel=1e-9)
            assert line["misorientation_deg"] >= 0.0

    def test_every_twin_segment_classifies_as_high_angle(self) -> None:
        result = draw(dataset="sigma3_twin", colouring="grain")
        summary = {row["character"]: row for row in result["data"]["boundary_summary"]}
        high = next(row for key, row in summary.items() if key.startswith("High-angle"))
        assert high["fraction"] == pytest.approx(1.0)
        assert high["mean_misorientation_deg"] == pytest.approx(SIGMA3_ANGLE_DEG, abs=1e-6)

    def test_boundary_fractions_are_measured_by_length_and_sum_to_one(self) -> None:
        rows = draw()["data"]["boundary_summary"]
        assert sum(row["fraction"] for row in rows) == pytest.approx(1.0)

    def test_raising_the_grain_threshold_cannot_create_grains(self) -> None:
        """A coarser criterion merges; it never splits."""

        fine = draw(grain_threshold_deg=2.0)["data"]["grain_count"]
        coarse = draw(grain_threshold_deg=20.0)["data"]["grain_count"]
        assert coarse <= fine


class TestWireContract:
    """What the panel is handed."""

    def test_the_image_is_rgb_at_the_native_grid_resolution(self) -> None:
        """One pixel is one measurement, so nothing is interpolated anywhere."""

        result = draw(grid_points=32)
        image = result["data"]["image"]
        assert image["encoding"] == "base64-rgb8"
        assert (image["height"], image["width"]) == tuple(result["data"]["grid_shape"])
        assert decode_rgb(image).shape == (32, 32, 3)

    def test_hover_columns_are_the_export_table_columns(self) -> None:
        result = draw()
        assert result["data"]["columns"] == result["table"]["columns"]

    def test_grain_rows_are_largest_first_and_carry_physical_sizes(self) -> None:
        rows = draw()["table"]["rows"]
        assert [row["size"] for row in rows] == sorted(
            (row["size"] for row in rows), reverse=True
        )
        for row in rows:
            expected = 2.0 * np.sqrt(row["area_um2"] / np.pi)
            assert row["equivalent_diameter_um"] == pytest.approx(expected)

    def test_the_extent_is_reported_in_specimen_units(self) -> None:
        result = draw(grid_points=32)
        step = result["data"]["step_um"]
        assert result["data"]["extent_um"] == [0.0, 0.0, 31 * step, 31 * step]

    def test_the_result_carries_the_dataset_known_answer(self) -> None:
        """A practice map is only a test if its answer travels with it."""

        result = draw(dataset="sigma3_twin")
        assert "60 degrees" in result["data"]["dataset"]["known_answer"]
        assert any("Known answer" in note for note in result["notes"])

    def test_an_unknown_dataset_names_the_field_and_the_choices(self) -> None:
        with pytest.raises(InvalidInputError) as raised:
            draw(dataset="not_a_dataset")
        assert raised.value.details["field"] == "dataset"


class TestNarration:
    """The EBSD module reports into the centralized log like every other."""

    def test_the_map_reports_its_grain_and_boundary_counts(self) -> None:
        from pytex.app.contracts import execute

        envelope, status = execute("ebsd.map", {**BASE_REQUEST, "grid_points": 32})
        assert status == 200
        messages = [record["message"] for record in envelope["log"]]
        assert any("grains found at" in message for message in messages)
        assert any("boundary segments" in message for message in messages)
        assert any("fundamental sector" in message for message in messages)
