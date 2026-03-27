from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.figure import Figure

from pytex import (
    ODF,
    FrameDomain,
    Handedness,
    InversePoleFigure,
    KernelSpec,
    Lattice,
    Orientation,
    OrientationSet,
    Phase,
    PoleFigure,
    QuaternionSet,
    ReferenceFrame,
    Rotation,
    RotationSet,
    SymmetrySpec,
    VectorSet,
    plot_euler_set,
    plot_inverse_pole_figure,
    plot_odf,
    plot_orientations,
    plot_pole_figure,
    plot_quaternion_set,
    plot_rotations,
    plot_symmetry_elements,
    plot_symmetry_orbit,
    plot_vector_set,
    save_documentation_figure_svg,
)
from pytex.core.lattice import CrystalPlane, MillerIndex
from pytex.plotting._render import MultiFigureSpec2D
from pytex.plotting.builders import (
    build_odf_figure_spec,
    build_pole_figure_spec,
    build_symmetry_orbit_figure_spec,
    build_vector_figure_spec,
)


def make_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    specimen = ReferenceFrame(
        "specimen",
        FrameDomain.SPECIMEN,
        ("x", "y", "z"),
        Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def make_orientation_set() -> tuple[OrientationSet, Phase]:
    crystal, specimen, phase = make_context()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            ),
            Orientation(
                rotation=Rotation.from_bunge_euler(35.0, 25.0, 10.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            ),
            Orientation(
                rotation=Rotation.from_bunge_euler(90.0, 0.0, 0.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            ),
        ]
    )
    return orientations, phase


def test_build_vector_figure_spec_preserves_frame_labels() -> None:
    crystal, _, _ = make_context()
    vectors = VectorSet(values=np.eye(3), reference_frame=crystal)
    spec = build_vector_figure_spec(vectors, normalize=True)
    assert spec.xlabel == "a"
    assert spec.ylabel == "b"
    assert spec.zlabel == "c"
    assert spec.unit_sphere_radius == 1.0


def test_build_symmetry_orbit_figure_spec_contains_sector_and_reduced_point() -> None:
    crystal, _, phase = make_context()
    spec = build_symmetry_orbit_figure_spec(
        phase.symmetry,
        VectorSet(values=np.array([[1.0, 0.0, 0.0]]), reference_frame=crystal),
    )
    assert len(spec.scatter_layers) == 2
    assert len(spec.line_layers) == 1
    assert spec.boundary_circle_radius is not None


def test_runtime_plotters_return_matplotlib_figures() -> None:
    orientations, phase = make_orientation_set()
    crystal, specimen, _ = make_context()
    vectors = VectorSet(values=np.eye(3), reference_frame=crystal)
    pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase)
    pole_figure = PoleFigure.from_orientations(
        orientations,
        pole,
        include_symmetry_family=False,
        antipodal=False,
    )
    inverse_pole_figure = InversePoleFigure.from_orientations(
        orientations,
        np.array([0.0, 0.0, 1.0]),
        reduce_by_symmetry=True,
    )
    odf = ODF.from_orientations(
        orientations,
        weights=[4.0, 2.0, 1.0],
        kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=10.0),
    )
    figures = [
        plot_vector_set(vectors, normalize=True),
        plot_symmetry_orbit(phase.symmetry, np.array([1.0, 0.0, 0.0])),
        plot_symmetry_elements(phase.symmetry),
        plot_euler_set(orientations.as_euler_set()),
        plot_quaternion_set(QuaternionSet(orientations.quaternions)),
        plot_rotations(RotationSet.from_rotations([Rotation.identity()])),
        plot_orientations(orientations),
        plot_pole_figure(pole_figure),
        plot_inverse_pole_figure(inverse_pole_figure),
        plot_odf(odf),
        plot_vector_set(np.eye(3), reference_frame=crystal),
    ]
    for figure in figures:
        assert isinstance(figure, Figure)
        plt.close(figure)
    del specimen


def test_plot_symmetry_orbit_rejects_frame_mismatch() -> None:
    crystal, specimen, phase = make_context()
    del crystal
    wrong_seed = VectorSet(values=np.array([[1.0, 0.0, 0.0]]), reference_frame=specimen)
    with pytest.raises(ValueError):
        plot_symmetry_orbit(phase.symmetry, wrong_seed)


def test_plot_odf_builder_and_runtime_png_export(tmp_path: Path) -> None:
    orientations, _ = make_orientation_set()
    odf = ODF.from_orientations(orientations, weights=[3.0, 2.0, 1.0])
    spec = build_odf_figure_spec(odf)
    assert spec.scatter_layers[0].values is not None
    figure = plot_odf(odf)
    png_path = tmp_path / "odf.png"
    figure.savefig(png_path, dpi=120)
    assert png_path.exists()
    plt.close(figure)


def test_pole_figure_histogram_builder_uses_raster_layer() -> None:
    orientations, phase = make_orientation_set()
    pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
    pole_figure = PoleFigure.from_orientations(
        orientations,
        pole,
        include_symmetry_family=False,
        antipodal=False,
    )
    spec = build_pole_figure_spec(pole_figure, kind="histogram")
    assert len(spec.image_layers) == 1


def test_pole_figure_contour_builder_uses_contour_layer() -> None:
    orientations, phase = make_orientation_set()
    pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
    pole_figure = PoleFigure.from_orientations(
        orientations,
        pole,
        include_symmetry_family=False,
        antipodal=False,
    )
    spec = build_pole_figure_spec(pole_figure, kind="contour")
    assert len(spec.contour_layers) == 1
    assert np.isnan(spec.contour_layers[0].values).any()


def test_plot_contour_pole_figure_and_odf_return_figures() -> None:
    orientations, phase = make_orientation_set()
    pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase)
    pole_figure = PoleFigure.from_orientations(
        orientations,
        pole,
        include_symmetry_family=False,
        antipodal=False,
    )
    odf = ODF.from_orientations(orientations, weights=[3.0, 2.0, 1.0])
    figures = [
        plot_pole_figure(pole_figure, kind="contour"),
        plot_odf(odf, kind="contour"),
    ]
    for figure in figures:
        assert isinstance(figure, Figure)
        plt.close(figure)


def test_build_odf_sections_returns_multi_panel_spec() -> None:
    orientations, _ = make_orientation_set()
    odf = ODF.from_orientations(orientations, weights=[3.0, 2.0, 1.0])
    spec = build_odf_figure_spec(odf, kind="sections", section_phi2_deg=(0.0, 45.0))
    assert isinstance(spec, MultiFigureSpec2D)
    assert len(spec.panels) == 2
    assert all(len(panel.contour_layers) == 1 for panel in spec.panels)


def test_plot_odf_sections_returns_matplotlib_figure() -> None:
    orientations, _ = make_orientation_set()
    odf = ODF.from_orientations(orientations, weights=[3.0, 2.0, 1.0])
    figure = plot_odf(odf, kind="sections", section_phi2_deg=(0.0, 45.0, 65.0))
    assert isinstance(figure, Figure)
    plt.close(figure)


def test_save_documentation_figure_svg_enforces_svg_suffix(tmp_path: Path) -> None:
    crystal, _, _ = make_context()
    vectors = VectorSet(values=np.eye(3), reference_frame=crystal)
    figure = plot_vector_set(vectors)
    svg_path = save_documentation_figure_svg(figure, tmp_path / "vectors.svg")
    assert svg_path.exists()
    with pytest.raises(ValueError):
        save_documentation_figure_svg(figure, tmp_path / "vectors.png")
    plt.close(figure)
