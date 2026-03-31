from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    generate_stereonet_grid,
    plot_crystal_directions,
    plot_crystal_planes,
    plot_symmetry_elements,
    plot_wulff_net,
    project_great_circle_trace,
)
from pytex.core.lattice import MillerIndex
from pytex.diffraction.stereonets import projection_boundary_radius
from pytex.plotting.spherical import build_symmetry_elements_figure_spec


def make_phase() -> Phase:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    return Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)


def test_project_great_circle_trace_for_basal_plane_hits_projection_boundary() -> None:
    trace = project_great_circle_trace(np.array([0.0, 0.0, 1.0]), method="stereographic")
    radii = np.linalg.norm(trace, axis=1)
    assert np.allclose(radii, 1.0, atol=1e-8)


def test_generate_stereonet_grid_stays_inside_projection_boundary() -> None:
    grid = generate_stereonet_grid(method="stereographic", major_step_deg=15.0, minor_step_deg=5.0)
    boundary = projection_boundary_radius("stereographic") + 1e-9
    for line in grid.major_lines + grid.minor_lines:
        assert np.all(np.linalg.norm(line, axis=1) <= boundary)


def test_generate_equal_area_stereonet_grid_stays_inside_projection_boundary() -> None:
    grid = generate_stereonet_grid(method="equal_area", major_step_deg=15.0, minor_step_deg=5.0)
    boundary = projection_boundary_radius("equal_area") + 1e-9
    for line in grid.major_lines + grid.minor_lines:
        assert np.all(np.linalg.norm(line, axis=1) <= boundary)


def test_symmetry_elements_builder_groups_orders_into_publication_symbols() -> None:
    phase = make_phase()
    spec = build_symmetry_elements_figure_spec(phase.symmetry)
    markers = {layer.marker for layer in spec.marker_layers}
    labels = {layer.label for layer in spec.marker_layers}
    assert {"^", "s", "D"} <= markers
    assert {"2-fold axes", "3-fold axes", "4-fold axes"} <= labels


def test_stereographic_plotters_return_matplotlib_figures() -> None:
    phase = make_phase()
    direction_figure = plot_crystal_directions(
        (
            CrystalDirection(np.array([1.0, 0.0, 0.0]), phase=phase),
            CrystalDirection(np.array([1.0, 1.0, 1.0]), phase=phase),
        ),
        labels=((1, 0, 0), (1, 1, 1)),
    )
    plane_figure = plot_crystal_planes(
        (
            CrystalPlane(MillerIndex([1, 1, 1], phase=phase), phase=phase),
            CrystalPlane(MillerIndex([1, 0, 0], phase=phase), phase=phase),
        ),
        labels=((1, 1, 1), (1, 0, 0)),
        render="both",
    )
    symmetry_figure = plot_symmetry_elements(phase.symmetry, annotate_axes=True)
    net_figure = plot_wulff_net()
    for figure in (direction_figure, plane_figure, symmetry_figure, net_figure):
        assert isinstance(figure, Figure)
        plt.close(figure)


def test_equal_area_plotting_and_hcp_labels_render_structurally() -> None:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(3.23, 3.23, 5.15, 90.0, 90.0, 120.0, crystal_frame=crystal)
    symmetry = SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal)
    phase = Phase("hcp_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)

    plane_figure = plot_crystal_planes(
        (
            CrystalPlane.from_miller_bravais((0, 0, 0, 1), phase=phase),
            CrystalPlane.from_miller_bravais((1, 0, -1, 0), phase=phase),
            CrystalPlane.from_miller_bravais((1, 0, -1, 1), phase=phase),
        ),
        labels=((0, 0, 0, 1), (1, 0, -1, 0), (1, 0, -1, 1)),
        method="equal_area",
        render="both",
        title="HCP Equal-Area Plane Plot",
    )
    symmetry_figure = plot_symmetry_elements(
        phase.symmetry,
        method="equal_area",
        annotate_axes=True,
        title="HCP Equal-Area Symmetry",
    )

    plane_ax = plane_figure.axes[0]
    texts = [text.get_text() for text in plane_ax.texts]
    assert any("\\bar{1}" in text for text in texts)
    assert plane_ax.get_title() == "HCP Equal-Area Plane Plot"
    assert len(plane_ax.collections) == 1
    assert len(plane_ax.lines) > 20

    symmetry_ax = symmetry_figure.axes[0]
    legend = symmetry_ax.get_legend()
    assert legend is not None
    assert {text.get_text() for text in legend.get_texts()} == {"2-fold axes", "6-fold axes"}

    plt.close(plane_figure)
    plt.close(symmetry_figure)
