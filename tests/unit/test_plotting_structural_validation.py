from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from pytex.plotting._plotting_validation_cases import build_plotting_validation_figures


def _legend_labels(ax: object) -> list[str]:
    legend = ax.get_legend()
    if legend is None:
        return []
    return [text.get_text() for text in legend.get_texts()]


def test_structural_plotting_validation_cases_match_expected_surface_properties() -> None:
    figures = build_plotting_validation_figures()

    xrd = figures["xrd_ni_fcc_journal"]
    xrd_ax = xrd.axes[0]
    assert len(xrd.axes) == 1
    assert xrd_ax.get_xlabel() == r"2$\theta$ (deg)"
    assert xrd_ax.get_ylabel() == "normalized intensity"
    assert len(xrd_ax.lines) >= 10
    assert len(xrd_ax.collections) == 1
    assert len(xrd_ax.texts) >= 10

    saed = figures["saed_ni_fcc_dark"]
    saed_ax = saed.axes[0]
    assert len(saed.axes) == 1
    assert saed_ax.get_xlabel() == "detector u (mm)"
    assert saed_ax.get_ylabel() == "detector v (mm)"
    assert len(saed_ax.collections) == 1
    assert len(saed_ax.texts) >= 8

    crystal = figures["crystal_zr_hcp_journal"]
    crystal_ax = crystal.axes[0]
    assert len(crystal.axes) == 1
    assert crystal_ax.name == "3d"
    assert len(crystal_ax.lines) >= 50
    assert len(crystal_ax.collections) >= 8
    assert len(crystal_ax.texts) >= 3

    ipf = figures["ipf_ni_fcc_journal"]
    assert len(ipf.axes) == 2
    ipf_ax = ipf.axes[0]
    colorbar_ax = ipf.axes[1]
    assert ipf_ax.get_title() == "Inverse Pole Figure"
    assert _legend_labels(ipf_ax) == ["directions", "fundamental sector"]
    assert len(ipf_ax.collections) == 1
    assert colorbar_ax.get_ylabel() == "intensity"

    direction = figures["stereonet_directions_ni_fcc_journal"]
    direction_ax = direction.axes[0]
    assert direction_ax.get_title() == "FCC Direction Stereonet"
    assert len(direction_ax.collections) == 1
    assert len(direction_ax.lines) >= 100
    assert len(direction_ax.texts) == 3
    assert _legend_labels(direction_ax) == ["directions"]

    plane = figures["stereonet_planes_ni_fcc_journal"]
    plane_ax = plane.axes[0]
    assert plane_ax.get_title() == "FCC Plane Traces"
    assert len(plane_ax.collections) == 1
    assert len(plane_ax.lines) >= 100
    assert len(plane_ax.texts) == 2
    assert _legend_labels(plane_ax) == ["plane poles"]

    symmetry = figures["symmetry_elements_ni_fcc_journal"]
    symmetry_ax = symmetry.axes[0]
    assert symmetry_ax.get_title() == "FCC Symmetry Elements"
    assert len(symmetry_ax.collections) == 3
    assert len(symmetry_ax.lines) >= 100
    assert len(symmetry_ax.texts) >= 10
    assert _legend_labels(symmetry_ax) == ["2-fold axes", "3-fold axes", "4-fold axes"]

    for figure in figures.values():
        plt.close(figure)
