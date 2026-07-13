from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pytex.plotting import resolve_style
from pytex.plotting.figure import (
    PanelGrid,
    add_scale_bar,
    export_figure,
    label_panels,
    publication_style,
    rc_params_from_style,
)


def test_rc_params_reflect_theme_values() -> None:
    style = resolve_style(theme="journal")
    rc = rc_params_from_style(style)
    common = style["common"]
    assert rc["figure.dpi"] == float(common["figure"]["dpi"])
    assert rc["font.family"] == common["font"]["family"]
    assert rc["font.size"] == float(common["font"]["size"])
    # publication export requirements: editable SVG text, TrueType PDF fonts
    assert rc["svg.fonttype"] == "none"
    assert rc["pdf.fonttype"] == 42
    assert rc["savefig.bbox"] == "tight"


def test_publication_style_applies_and_restores_rcparams() -> None:
    baseline = matplotlib.rcParams["font.size"]
    with publication_style(theme="journal") as style:
        assert matplotlib.rcParams["font.size"] == float(style["common"]["font"]["size"])
        assert matplotlib.rcParams["svg.fonttype"] == "none"
    assert matplotlib.rcParams["font.size"] == baseline


def test_panel_grid_shape_labels_and_shared_colorbar() -> None:
    grid = PanelGrid(2, 3, panel_size=(2.0, 2.0))
    try:
        assert grid.axes.shape == (2, 3)
        assert len(grid.axes_flat) == 6
        image = grid.axes_flat[0].imshow(np.arange(9.0).reshape(3, 3))
        labels = grid.label()
        assert [artist.get_text() for artist in labels] == [
            "(a)",
            "(b)",
            "(c)",
            "(d)",
            "(e)",
            "(f)",
        ]
        colorbar = grid.shared_colorbar(image, label="intensity (m.r.d.)")
        assert colorbar.ax.get_ylabel() == "intensity (m.r.d.)"
        grid.hide_unused(4)
        assert not grid.axes_flat[5].get_visible()
    finally:
        plt.close(grid.figure)


def test_panel_grid_rejects_empty_grid() -> None:
    with pytest.raises(ValueError, match="at least one"):
        PanelGrid(0, 2)


def test_label_panels_custom_format_and_limit() -> None:
    fig, axes = plt.subplots(1, 2)
    try:
        artists = label_panels(list(axes), fmt="{})")
        assert [artist.get_text() for artist in artists] == ["a)", "b)"]
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="26"):
        label_panels([object()] * 27)


def test_add_scale_bar_attaches_artist() -> None:
    fig, axis = plt.subplots()
    try:
        axis.imshow(np.zeros((10, 10)), extent=(0.0, 25.0, 0.0, 25.0))
        artist = add_scale_bar(axis, 5.0, units="µm")
        assert artist in axis.artists
        with pytest.raises(ValueError, match="positive"):
            add_scale_bar(axis, 0.0)
    finally:
        plt.close(fig)


def test_export_figure_writes_requested_formats(tmp_path: Path) -> None:
    fig, axis = plt.subplots()
    try:
        axis.plot([0.0, 1.0], [0.0, 1.0])
        written = export_figure(fig, tmp_path / "figures" / "demo", formats=("svg", "png"))
        assert tuple(path.suffix for path in written) == (".svg", ".png")
        for path in written:
            assert path.exists() and path.stat().st_size > 0
        # explicit suffix exports exactly that format
        (single,) = export_figure(fig, tmp_path / "demo.pdf")
        assert single.suffix == ".pdf" and single.exists()
        with pytest.raises(ValueError, match="Unsupported"):
            export_figure(fig, tmp_path / "demo", formats=("tiff",))
    finally:
        plt.close(fig)
