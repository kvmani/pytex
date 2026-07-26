"""Structural regression tests for composite SAED rendering (CD3).

Per repo policy these are semantic assertions on the matplotlib artists
(collections, gids, colors, sizes, labels, limits) — no image baselines.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import PathCollection
from matplotlib.figure import Figure

from pytex.core.lattice import ZoneAxis
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import CompositeSAEDPattern, simulate_composite_saed
from pytex.plotting.composite_saed import (
    GID_PREFIX,
    VARIANT_COLOR_PALETTE,
    CompositeSAEDPlotConfig,
    SpotStyle,
    render_composite_saed,
)
from tests.unit.test_composite_saed import make_fcc_bcc_phases


@pytest.fixture(scope="module")
def composite() -> CompositeSAEDPattern:
    parent, child = make_fcc_bcc_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    zone = ZoneAxis(np.array([0, 1, -1]), phase=parent)
    return simulate_composite_saed(ks, zone, variant_indices=(1, 2, 3))


def _collections_by_gid(fig: Figure) -> dict[str, PathCollection]:
    result: dict[str, PathCollection] = {}
    for axes in fig.axes:
        for collection in axes.collections:
            gid = collection.get_gid()
            if gid and gid.startswith(GID_PREFIX):
                result[gid] = collection
    return result


class TestRenderStructure:
    def test_one_collection_per_subpattern_plus_beam(
        self, composite: CompositeSAEDPattern
    ) -> None:
        fig = render_composite_saed(composite)
        by_gid = _collections_by_gid(fig)
        assert f"{GID_PREFIX}:parent" in by_gid
        assert f"{GID_PREFIX}:transmitted-beam" in by_gid
        for index in (1, 2, 3):
            assert f"{GID_PREFIX}:variant:{index}" in by_gid
        assert len(by_gid) == 5

    def test_offsets_match_spot_tables(self, composite: CompositeSAEDPattern) -> None:
        fig = render_composite_saed(composite)
        by_gid = _collections_by_gid(fig)
        assert composite.parent_spots is not None
        parent_offsets = np.asarray(by_gid[f"{GID_PREFIX}:parent"].get_offsets())
        assert np.allclose(parent_offsets, composite.parent_spots.detector_mm)
        for pattern in composite.variant_patterns:
            offsets = np.asarray(
                by_gid[f"{GID_PREFIX}:variant:{pattern.variant_index}"].get_offsets()
            )
            assert np.allclose(offsets, pattern.spots.detector_mm)

    def test_variant_subset_config(self, composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(variant_indices=(2,))
        fig = render_composite_saed(composite, config=config)
        by_gid = _collections_by_gid(fig)
        assert f"{GID_PREFIX}:variant:2" in by_gid
        assert f"{GID_PREFIX}:variant:1" not in by_gid
        assert f"{GID_PREFIX}:variant:3" not in by_gid

    def test_show_parent_false(self, composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(show_parent=False)
        fig = render_composite_saed(composite, config=config)
        assert f"{GID_PREFIX}:parent" not in _collections_by_gid(fig)

    def test_show_transmitted_beam_false(self, composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(show_transmitted_beam=False)
        fig = render_composite_saed(composite, config=config)
        assert f"{GID_PREFIX}:transmitted-beam" not in _collections_by_gid(fig)

    def test_render_on_existing_axes(self, composite: CompositeSAEDPattern) -> None:
        fig, axes = plt.subplots()
        result = render_composite_saed(composite, ax=axes)
        assert result is fig
        assert _collections_by_gid(fig)

    def test_equal_aspect_and_symmetric_limits(
        self, composite: CompositeSAEDPattern
    ) -> None:
        fig = render_composite_saed(composite)
        axes = fig.axes[0]
        assert axes.get_aspect() == 1.0
        x_low, x_high = axes.get_xlim()
        assert x_low == pytest.approx(-x_high)
        coords = composite.all_detector_coordinates()
        assert x_high >= float(np.max(np.abs(coords)))


class TestStyles:
    def test_explicit_variant_style_override(self, composite: CompositeSAEDPattern) -> None:
        custom = SpotStyle(marker="X", color="#00ff00", size_mode="constant", size_scale=50.0)
        config = CompositeSAEDPlotConfig(variant_styles={2: custom})
        fig = render_composite_saed(composite, config=config)
        by_gid = _collections_by_gid(fig)
        collection = by_gid[f"{GID_PREFIX}:variant:2"]
        assert np.allclose(collection.get_sizes(), 50.0)
        face = collection.get_facecolor()
        assert np.allclose(face[0][:3], (0.0, 1.0, 0.0), atol=1e-6)

    def test_palette_cycling_by_render_position(
        self, composite: CompositeSAEDPattern
    ) -> None:
        fig = render_composite_saed(composite)
        by_gid = _collections_by_gid(fig)
        from matplotlib.colors import to_rgb

        for position, index in enumerate((1, 2, 3)):
            face = by_gid[f"{GID_PREFIX}:variant:{index}"].get_facecolor()
            expected = to_rgb(VARIANT_COLOR_PALETTE[position])
            assert np.allclose(face[0][:3], expected, atol=1e-6)

    def test_parent_hollow_markers(self, composite: CompositeSAEDPattern) -> None:
        fig = render_composite_saed(composite)
        by_gid = _collections_by_gid(fig)
        parent = by_gid[f"{GID_PREFIX}:parent"]
        assert parent.get_facecolor().size == 0

    def test_intensity_area_size_mode(self, composite: CompositeSAEDPattern) -> None:
        style = SpotStyle(size_mode="intensity_area", size_scale=100.0, min_size_pt2=1.0)
        config = CompositeSAEDPlotConfig(variant_styles={1: style})
        fig = render_composite_saed(composite, config=config)
        by_gid = _collections_by_gid(fig)
        sizes = by_gid[f"{GID_PREFIX}:variant:1"].get_sizes()
        intensity = composite.variant_patterns[0].spots.intensity
        assert np.allclose(sizes, np.maximum(100.0 * intensity, 1.0))

    def test_intensity_radius_size_mode(self, composite: CompositeSAEDPattern) -> None:
        style = SpotStyle(size_mode="intensity_radius", size_scale=100.0, min_size_pt2=1.0)
        sizes = style.marker_sizes_pt2(np.array([1.0, 0.5]))
        assert sizes[0] == pytest.approx(100.0)
        assert sizes[1] == pytest.approx(25.0)

    def test_invalid_style_raises(self) -> None:
        with pytest.raises(ValueError):
            SpotStyle(size_scale=0.0)
        with pytest.raises(ValueError):
            SpotStyle(size_mode="huge")  # type: ignore[arg-type]
        with pytest.raises(ValueError):
            SpotStyle(alpha=0.0)


class TestAxesUnitsAndLegend:
    def test_inv_angstrom_units_rescale_offsets(
        self, composite: CompositeSAEDPattern
    ) -> None:
        config = CompositeSAEDPlotConfig(axes_units="inv_angstrom")
        fig = render_composite_saed(composite, config=config)
        by_gid = _collections_by_gid(fig)
        assert composite.parent_spots is not None
        offsets = np.asarray(by_gid[f"{GID_PREFIX}:parent"].get_offsets())
        assert np.allclose(offsets, composite.parent_spots.g_detector_inv_angstrom)
        assert "mm" not in fig.axes[0].get_xlabel()

    def test_mm_axis_labels(self, composite: CompositeSAEDPattern) -> None:
        fig = render_composite_saed(composite)
        assert fig.axes[0].get_xlabel() == "detector u (mm)"
        assert fig.axes[0].get_ylabel() == "detector v (mm)"

    def test_legend_contains_phase_and_variant_labels(
        self, composite: CompositeSAEDPattern
    ) -> None:
        fig = render_composite_saed(composite)
        legend = fig.axes[0].get_legend()
        assert legend is not None
        labels = [entry.get_text() for entry in legend.get_texts()]
        assert any("austenite" in label for label in labels)
        assert any("martensite V1" in label for label in labels)
        assert any("transmitted beam" in label for label in labels)

    def test_legend_disabled(self, composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(show_legend=False)
        fig = render_composite_saed(composite, config=config)
        assert fig.axes[0].get_legend() is None

    def test_legend_entry_cap(self, composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(legend_max_entries=2, legend_outside=False)
        fig = render_composite_saed(composite, config=config)
        legend = fig.axes[0].get_legend()
        assert legend is not None
        assert len(legend.get_texts()) == 2

    def test_default_title_mentions_relationship_and_zone(
        self, composite: CompositeSAEDPattern
    ) -> None:
        fig = render_composite_saed(composite)
        title = fig.axes[0].get_title()
        assert "kurdjumov_sachs" in title
        assert "[0 1 -1]" in title

    def test_custom_title_and_disabled_title(self, composite: CompositeSAEDPattern) -> None:
        fig = render_composite_saed(
            composite, config=CompositeSAEDPlotConfig(title="My pattern")
        )
        assert fig.axes[0].get_title() == "My pattern"
        fig2 = render_composite_saed(
            composite, config=CompositeSAEDPlotConfig(show_title=False)
        )
        assert fig2.axes[0].get_title() == ""


class TestConfigValidation:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"variant_color_palette": ()},
            {"variant_marker_cycle": ()},
            {"axes_units": "pixels"},
            {"legend_max_entries": 0},
            {"limit_padding_fraction": 1.5},
            {"transmitted_beam_size_pt2": 0.0},
            {"dpi": 0},
        ],
    )
    def test_invalid_config_raises(self, kwargs: dict[str, object]) -> None:
        with pytest.raises(ValueError):
            CompositeSAEDPlotConfig(**kwargs)  # type: ignore[arg-type]

    def test_style_for_variant_override_and_cycle(self) -> None:
        custom = SpotStyle(color="#123456")
        config = CompositeSAEDPlotConfig(variant_styles={7: custom})
        assert config.style_for_variant(7, 0) is custom
        cycled = config.style_for_variant(3, 1)
        assert cycled.color == VARIANT_COLOR_PALETTE[1]


def test_composite_frame_indicator_is_opt_in(composite: CompositeSAEDPattern) -> None:
    default_figure = render_composite_saed(composite)
    try:
        assert list(default_figure.axes[0].child_axes) == []
    finally:
        plt.close(default_figure)


def test_composite_frame_indicator_shows_parent_crystal_axes_on_the_detector(
    composite: CompositeSAEDPattern,
) -> None:
    """The gizmo must report the parent crystal axes, not the detector axes.

    A composite pattern's detector axes are trivially the page axes; what a
    reader actually needs is where the parent crystal's a/b/c point on this
    detector, which the pattern's parent-anchored zone basis supplies.
    """

    figure = render_composite_saed(
        composite, config=CompositeSAEDPlotConfig(show_frame_indicator=True)
    )
    try:
        insets = list(figure.axes[0].child_axes)
        assert len(insets) == 1
        labels = {text.get_text() for text in insets[0].texts if text.get_text()}
        parent_frame = composite.relationship.parent_phase.crystal_frame
        assert set(parent_frame.axes) <= labels
        assert parent_frame.name in labels
    finally:
        plt.close(figure)
