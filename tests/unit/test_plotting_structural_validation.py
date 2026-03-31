from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pytex import (
    CrystalMap,
    FrameDomain,
    Handedness,
    Lattice,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    plot_ipf_map,
    plot_kam_map,
)
from pytex.plotting._plotting_validation_cases import build_plotting_validation_figures


def _legend_labels(ax: object) -> list[str]:
    legend = ax.get_legend()
    if legend is None:
        return []
    return [text.get_text() for text in legend.get_texts()]


def test_structural_plotting_validation_cases_match_expected_surface_properties() -> None:
    pytest.importorskip("pymatgen.core")
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


def test_ebsd_runtime_plotting_surfaces_have_expected_structure() -> None:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    orientations = OrientationSet.from_euler_angles(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [30.0, 0.0, 0.0],
                [0.0, 30.0, 0.0],
                [30.0, 30.0, 0.0],
            ]
        ),
        crystal_frame=crystal,
        specimen_frame=specimen,
        phase=phase,
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )

    ipf = plot_ipf_map(crystal_map, direction="z")
    assert len(ipf.axes) == 1
    assert len(ipf.axes[0].images) == 1
    assert ipf.axes[0].get_title() == "IPF Map (z)"

    kam = plot_kam_map(crystal_map)
    assert len(kam.axes) == 2
    assert len(kam.axes[0].images) == 1
    assert kam.axes[1].get_ylabel() == "KAM (deg)"

    plt.close(ipf)
    plt.close(kam)
