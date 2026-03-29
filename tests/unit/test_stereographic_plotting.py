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
