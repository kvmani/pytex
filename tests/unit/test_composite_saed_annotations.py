"""Regression tests for composite SAED spot annotation (CD4).

Structural assertions per repo policy: label texts, gids, pairwise-disjoint
rendered label boxes, merged coincident labels, budget/floor behavior,
determinism — no image baselines.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from dataclasses import replace as dataclass_replace

import numpy as np
import pytest
from matplotlib.figure import Figure
from matplotlib.text import Annotation

from pytex.core.lattice import Phase, UnitCell, ZoneAxis
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import CompositeSAEDPattern, simulate_composite_saed
from pytex.plotting.composite_saed import (
    GID_PREFIX,
    AnnotationResult,
    CompositeSAEDPlotConfig,
    SpotAnnotationConfig,
    format_hkl,
    render_composite_saed,
)
from tests.unit.test_composite_saed import make_bcc_hcp_phases, make_fcc_bcc_phases


def _with_b_iso(phase: Phase, b_iso: float) -> Phase:
    assert phase.unit_cell is not None
    sites = tuple(
        dataclass_replace(site, b_iso=b_iso) for site in phase.unit_cell.sites
    )
    return dataclass_replace(
        phase, unit_cell=UnitCell(lattice=phase.lattice, sites=sites)
    )


@pytest.fixture(scope="module")
def ks_composite() -> CompositeSAEDPattern:
    parent, child = make_fcc_bcc_phases()
    ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
        parent_phase=parent, child_phase=child
    )
    zone = ZoneAxis(np.array([0, 1, -1]), phase=parent)
    return simulate_composite_saed(ks, zone, variant_indices=(1, 2, 3, 4))


def _annotation_artists(fig: Figure) -> list[Annotation]:
    artists = []
    for axes in fig.axes:
        for child in axes.get_children():
            gid = child.get_gid()
            if (
                isinstance(child, Annotation)
                and gid
                and gid.startswith(f"{GID_PREFIX}:annotation:")
            ):
                artists.append(child)
    return artists


def _leader_artists(fig: Figure) -> list[Annotation]:
    artists = []
    for axes in fig.axes:
        for child in axes.get_children():
            gid = child.get_gid()
            if gid and gid.startswith(f"{GID_PREFIX}:leader:"):
                artists.append(child)
    return artists


class TestFormatHkl:
    def test_plain_format(self) -> None:
        assert format_hkl([1, 1, -1], index_format="plain") == "(1 1 -1)"
        assert format_hkl([12, -1, 1], index_format="plain") == "(12 -1 1)"

    def test_overline_compact_single_digit(self) -> None:
        assert format_hkl([1, 1, -1]) == r"$(11\bar{1})$"
        assert format_hkl([0, 2, 0]) == "$(020)$"
        assert format_hkl([-1, -1, -1]) == r"$(\bar{1}\bar{1}\bar{1})$"

    def test_overline_spaced_multi_digit(self) -> None:
        assert format_hkl([12, -1, 1]) == r"$(12\;\bar{1}\;1)$"

    def test_invalid_format_raises(self) -> None:
        with pytest.raises(ValueError, match="index_format"):
            format_hkl([1, 0, 0], index_format="fancy")  # type: ignore[arg-type]

    def test_bravais_expands_to_four_indices(self) -> None:
        # (h k l) -> (h k i l) with i = -(h + k): the hexagonal convention.
        assert format_hkl([1, 0, 0], index_format="plain", bravais=True) == "(1 0 -1 0)"
        assert format_hkl([0, 0, 2], index_format="plain", bravais=True) == "(0 0 0 2)"
        assert format_hkl([1, 1, 0], index_format="plain", bravais=True) == "(1 1 -2 0)"

    def test_bravais_overline_form(self) -> None:
        assert format_hkl([1, 0, 0], bravais=True) == r"$(10\bar{1}0)$"


class TestHexagonalAnnotations:
    """Burgers beta->alpha: hexagonal child labels must use four indices."""

    @pytest.fixture(scope="class")
    def burgers_composite(self) -> CompositeSAEDPattern:
        beta, alpha = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=beta, child_phase=alpha
        )
        zone = ZoneAxis(np.array([1, 1, 0]), phase=beta)
        return simulate_composite_saed(relationship, zone, variant_indices=(1, 2))

    def test_child_labels_four_index_parent_labels_three_index(
        self, burgers_composite: CompositeSAEDPattern
    ) -> None:
        config = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(index_format="plain", max_labels=40)
        )
        _, result = render_composite_saed(
            burgers_composite, config=config, return_annotations=True
        )
        child_lines = [
            line
            for text in result.texts
            for line in text.split("\n")
            if line.endswith(("V1", "V2"))
        ]
        parent_lines = [
            line for text in result.texts for line in text.split("\n") if line.endswith(" p")
        ]
        assert child_lines and parent_lines
        for line in child_lines:
            indices = line.split(")")[0].lstrip("(").split()
            assert len(indices) == 4, f"hexagonal child label not four-index: {line}"
        for line in parent_lines:
            indices = line.split(")")[0].lstrip("(").split()
            assert len(indices) == 3, f"cubic parent label not three-index: {line}"


class TestAnnotationRendering:
    def test_labels_placed_with_gids_and_result(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        fig, result = render_composite_saed(ks_composite, return_annotations=True)
        artists = _annotation_artists(fig)
        assert isinstance(result, AnnotationResult)
        assert result.placed_count == len(artists) == len(result.texts)
        assert result.placed_count > 0
        assert result.placed_count + result.skipped_count == result.cluster_count
        assert result.cluster_count <= ks_composite.spot_count()

    def test_merged_coincident_labels_multiline_phase_tagged(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        config = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(
                coincidence_tolerance_mm=3.0, max_labels=40, min_intensity=0.0
            )
        )
        fig, result = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        assert result.merged_cluster_count > 0
        merged = [text for text in result.texts if "\n" in text]
        assert merged
        assert any(
            " p" in text and " V" in text for text in merged
        ), "expected at least one merged parent/variant label"

    def test_label_boxes_pairwise_disjoint(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        config = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(max_labels=40, min_intensity=0.0)
        )
        fig, _ = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        extents = [
            artist.get_window_extent(renderer=renderer)
            for artist in _annotation_artists(fig)
        ]
        assert len(extents) > 5
        for i in range(len(extents)):
            for j in range(i + 1, len(extents)):
                left, right = extents[i], extents[j]
                overlap = not (
                    left.x1 <= right.x0
                    or right.x1 <= left.x0
                    or left.y1 <= right.y0
                    or right.y1 <= left.y0
                )
                assert not overlap, f"label boxes {i} and {j} overlap"

    def test_max_labels_budget(self, ks_composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(max_labels=5, min_intensity=0.0)
        )
        fig, result = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        assert result.cluster_count <= 5
        assert len(_annotation_artists(fig)) <= 5

    def test_min_intensity_floor(self) -> None:
        # The Z-only scattering proxy with b_iso = 0 gives flat unit
        # intensities, so the floor needs Debye-Waller-damped phases to bite;
        # the deliberately huge B value spreads normalized intensities well
        # below the 0.9 floor within the default g range.
        parent, child = make_fcc_bcc_phases()
        damped_parent = _with_b_iso(parent, 25.0)
        damped_child = _with_b_iso(child, 25.0)
        ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=damped_parent, child_phase=damped_child
        )
        zone = ZoneAxis(np.array([0, 1, -1]), phase=damped_parent)
        composite = simulate_composite_saed(ks, zone, variant_indices=(2,))
        low = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(min_intensity=0.0, max_labels=500)
        )
        high = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(min_intensity=0.9, max_labels=500)
        )
        _, low_result = render_composite_saed(composite, config=low, return_annotations=True)
        _, high_result = render_composite_saed(
            composite, config=high, return_annotations=True
        )
        assert high_result.cluster_count < low_result.cluster_count

    def test_disabled_annotations(self, ks_composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(annotation=SpotAnnotationConfig(enabled=False))
        fig, result = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        assert not _annotation_artists(fig)
        assert result.cluster_count == 0
        assert result.placed_count == 0

    def test_deterministic_output(self, ks_composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(max_labels=30, min_intensity=0.0)
        )
        _, first = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        _, second = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        assert first.texts == second.texts
        assert np.array_equal(first.positions_data, second.positions_data)

    def test_plain_format_flows_into_labels(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        config = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(index_format="plain", max_labels=10)
        )
        _, result = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        assert result.texts
        assert all("$" not in text for text in result.texts)

    def test_overline_format_flows_into_labels(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        _, result = render_composite_saed(ks_composite, return_annotations=True)
        assert any("\\bar" in text for text in result.texts)

    def test_leader_lines_emitted_for_outer_ring_placements(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        config = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(
                leader_lines=True,
                max_labels=60,
                min_intensity=0.0,
                coincidence_tolerance_mm=0.5,
                font_size=11.0,
            )
        )
        fig, result = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        leaders = _leader_artists(fig)
        assert len(leaders) <= result.placed_count
        # Dense composite with large fonts must push some labels off ring 1
        # (or drop them); accept either but require the pass to have engaged.
        assert result.skipped_count + len(leaders) > 0

    def test_describe_content(self, ks_composite: CompositeSAEDPattern) -> None:
        config = CompositeSAEDPlotConfig(
            annotation=SpotAnnotationConfig(coincidence_tolerance_mm=3.0, max_labels=40)
        )
        _, result = render_composite_saed(
            ks_composite, config=config, return_annotations=True
        )
        text = result.describe()
        assert "label" in text
        assert str(result.placed_count) in text
        assert "coincidence tolerance" in text

    def test_default_render_returns_figure_only(
        self, ks_composite: CompositeSAEDPattern
    ) -> None:
        fig = render_composite_saed(ks_composite)
        assert isinstance(fig, Figure)


class TestAnnotationConfigValidation:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"max_labels": 0},
            {"min_intensity": -0.1},
            {"min_intensity": 1.5},
            {"index_format": "fancy"},
            {"coincidence_tolerance_mm": -1.0},
            {"offset_pt": 0.0},
            {"font_size": 0.0},
            {"bbox_alpha": 1.2},
        ],
    )
    def test_invalid_config_raises(self, kwargs: dict[str, object]) -> None:
        with pytest.raises(ValueError):
            SpotAnnotationConfig(**kwargs)  # type: ignore[arg-type]

    def test_annotation_result_consistency_checks(self) -> None:
        with pytest.raises(ValueError, match="cluster_count"):
            AnnotationResult(
                texts=("a",),
                positions_data=np.zeros((1, 2)),
                cluster_count=3,
                placed_count=1,
                skipped_count=1,
                merged_cluster_count=0,
            )
        with pytest.raises(ValueError, match="placed_count"):
            AnnotationResult(
                texts=("a", "b"),
                positions_data=np.zeros((2, 2)),
                cluster_count=2,
                placed_count=1,
                skipped_count=1,
                merged_cluster_count=0,
            )
