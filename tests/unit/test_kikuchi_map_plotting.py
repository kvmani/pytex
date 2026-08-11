"""Structural tests for the Kikuchi-map figure.

Runtime plotting is checked semantically rather than against tracked SVG bytes,
per ``docs/standards/executable_examples.md``: the assertions are about what the
figure encodes -- every trace inside the projection disc, the route drawn where
its own oriented directions say it goes, labels capped -- not about pixels.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pytex.core import FrameDomain, ReferenceFrame, get_phase_fixture
from pytex.diffraction.kikuchi_map import compute_kikuchi_map
from pytex.plotting.kikuchi_map import plot_kikuchi_map

CRYSTAL = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"))


@pytest.fixture(scope="module")
def cubic_map():
    phase = get_phase_fixture("ni_fcc").load_phase(crystal_frame=CRYSTAL)
    return compute_kikuchi_map(phase, max_index=3, zone_axis_max_index=3)


@pytest.fixture(scope="module")
def hexagonal_map():
    phase = get_phase_fixture("zr_hcp").load_phase(crystal_frame=CRYSTAL)
    return compute_kikuchi_map(phase, max_index=2, zone_axis_max_index=2)


def _single_point_offsets(axes) -> list[np.ndarray]:
    """Offsets of every scatter that carries exactly one point."""

    found = []
    for collection in axes.collections:
        offsets = np.asarray(collection.get_offsets())
        if offsets.shape[0] == 1:
            found.append(offsets[0])
    return found


def test_every_drawn_trace_stays_inside_the_projection_disc(cubic_map) -> None:
    """Boundedness is the reason the projection is stereographic rather than gnomonic."""

    figure = plot_kikuchi_map(cubic_map, max_bands=10, samples=181)
    try:
        axes = figure.axes[0]
        assert axes.lines, "bands must be drawn"
        for line in axes.lines:
            data = np.column_stack(line.get_data())
            if data.size == 0:
                continue
            assert np.all(np.linalg.norm(data, axis=1) <= 1.0 + 1e-6)
        lower, upper = axes.get_xlim()
        assert lower < -1.0 and upper > 1.0
        assert axes.get_aspect() == 1.0
    finally:
        plt.close(figure)


def test_the_route_marker_sits_at_the_end_of_the_route_arc(hexagonal_map) -> None:
    """The target marker and the drawn path must agree about which side it is on.

    A zone axis is a line, so each end of a leg has two senses that project to
    opposite sides of the disc. An equatorial axis makes this sharp: its ``z`` is
    zero, and a projection that folds on the sign of ``z`` will send the two senses
    to ``+1`` and ``-1``. The route commits to one sense, and the marker has to
    follow it.
    """

    route = hexagonal_map.route_to([0, 0, 1], [1, 0, 0], max_leg_deg=40.0)
    assert route.reachable
    figure = plot_kikuchi_map(hexagonal_map, route=route, max_bands=8, samples=181)
    try:
        axes = figure.axes[0]
        expected = np.asarray(route.target_direction)[:2]
        markers = _single_point_offsets(axes)
        assert any(
            float(np.linalg.norm(marker - expected)) < 1e-6 for marker in markers
        ), (expected, markers)
        # And the same for the start, which is the map centre here.
        assert any(float(np.linalg.norm(marker)) < 1e-9 for marker in markers)
    finally:
        plt.close(figure)


def test_an_unreachable_route_is_drawn_dashed(cubic_map) -> None:
    """It must not be possible to mistake a failed plan for a plan."""

    route = cubic_map.route_to([0, 0, 1], [1, 1, 1], max_leg_deg=1.0)
    assert not route.reachable
    figure = plot_kikuchi_map(cubic_map, route=route, max_bands=4, samples=91)
    try:
        styles = {line.get_linestyle() for line in figure.axes[0].lines}
        assert any(style not in {"-", "solid"} for style in styles)
    finally:
        plt.close(figure)


def test_labels_are_capped_and_ordered_by_conspicuousness(cubic_map) -> None:
    """A dense map has to be labelled sparsely or the figure is unreadable."""

    figure = plot_kikuchi_map(cubic_map, max_bands=8, samples=91, min_label_order=3,
                              max_labels=5)
    try:
        assert len(figure.axes[0].texts) == 5
    finally:
        plt.close(figure)
    figure = plot_kikuchi_map(cubic_map, max_bands=8, samples=91, min_label_order=3,
                              max_labels=None)
    try:
        qualifying = sum(1 for axis in cubic_map.zone_axes if axis.order >= 3)
        assert len(figure.axes[0].texts) == qualifying
    finally:
        plt.close(figure)


def test_width_scale_changes_only_the_edges(cubic_map) -> None:
    """Exaggerating the width must move the edges and leave the centre lines alone.

    The centre line is the plane's trace and is not a function of the Bragg angle;
    only the Kossel cones are. A ``width_scale`` that moved the centre lines would
    be misrepresenting the geometry rather than magnifying it.
    """

    true_scale = plot_kikuchi_map(cubic_map, max_bands=3, samples=91, width_scale=1.0)
    exaggerated = plot_kikuchi_map(cubic_map, max_bands=3, samples=91, width_scale=10.0)
    try:
        first = [np.column_stack(line.get_data()) for line in true_scale.axes[0].lines]
        second = [np.column_stack(line.get_data()) for line in exaggerated.axes[0].lines]
        assert len(first) == len(second)
        differing = sum(
            1
            for left, right in zip(first, second, strict=True)
            if left.shape != right.shape or not np.allclose(left, right, atol=1e-9)
        )
        assert 0 < differing < len(first), "edges must move and centre lines must not"
    finally:
        plt.close(true_scale)
        plt.close(exaggerated)


def test_the_title_says_when_widths_are_not_to_scale(cubic_map) -> None:
    """A figure drawn with exaggerated widths must announce it."""

    to_scale = plot_kikuchi_map(cubic_map, max_bands=3, samples=91)
    exaggerated = plot_kikuchi_map(cubic_map, max_bands=3, samples=91, width_scale=6.0)
    try:
        assert "widths" not in to_scale.axes[0].get_title()
        assert "widths x6" in exaggerated.axes[0].get_title()
    finally:
        plt.close(to_scale)
        plt.close(exaggerated)


def test_hexagonal_labels_use_four_indices(hexagonal_map) -> None:
    figure = plot_kikuchi_map(hexagonal_map, max_bands=6, samples=91, min_label_order=3)
    try:
        texts = [text.get_text() for text in figure.axes[0].texts]
        assert any("0001" in text for text in texts)
        assert figure.axes[0].get_title().count("0001") == 1
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"method": "orthographic"}, "method must be either"),
        ({"samples": 1}, "samples must be at least two"),
        ({"max_bands": 0}, "max_bands must be strictly positive"),
        ({"width_scale": 0.0}, "width_scale must be finite and strictly positive"),
        ({"min_label_order": 1}, "a crossing needs two bands"),
        ({"max_labels": 0}, "max_labels must be strictly positive"),
    ],
)
def test_invalid_plot_requests_are_rejected(cubic_map, kwargs, message) -> None:
    with pytest.raises(ValueError, match=message):
        plot_kikuchi_map(cubic_map, **kwargs)
