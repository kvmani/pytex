from __future__ import annotations

from typing import Any

import numpy as np

from pytex import (
    CrystalCellOverlay,
    CrystalDirection,
    CrystalDirectionOverlay,
    CrystalPlane,
    CrystalPlaneOverlay,
    DirectionAnnotationStyle,
    FrameDomain,
    Handedness,
    InversePoleFigure,
    MillerIndex,
    Orientation,
    OrientationSet,
    PlaneAnnotationStyle,
    ReferenceFrame,
    Rotation,
    ZoneAxis,
    build_crystal_scene,
    generate_saed_pattern,
    generate_xrd_pattern,
    get_phase_fixture,
    plot_crystal_directions,
    plot_crystal_planes,
    plot_crystal_structure_3d,
    plot_inverse_pole_figure,
    plot_saed_pattern,
    plot_symmetry_elements,
    plot_xrd_pattern,
)


def _make_crystal_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


def _make_specimen_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )


def _fixture_phase(fixture_id: str) -> Any:
    return get_phase_fixture(fixture_id).load_phase(crystal_frame=_make_crystal_frame())


def build_plotting_validation_figures() -> dict[str, Any]:
    ni_fcc = _fixture_phase("ni_fcc")
    zr_hcp = _fixture_phase("zr_hcp")
    specimen = _make_specimen_frame()

    xrd_pattern = generate_xrd_pattern(
        ni_fcc,
        two_theta_range_deg=(20.0, 120.0),
        resolution_deg=0.02,
        max_index=6,
        broadening_fwhm_deg=0.16,
    )
    xrd_figure = plot_xrd_pattern(xrd_pattern, theme="journal")

    saed_pattern = generate_saed_pattern(
        ni_fcc,
        ZoneAxis(np.array([0, 0, 1]), phase=ni_fcc),
        camera_constant_mm_angstrom=180.0,
        max_index=4,
        max_g_inv_angstrom=1.0,
    )
    saed_figure = plot_saed_pattern(saed_pattern, theme="dark")

    zr_scene = build_crystal_scene(
        zr_hcp,
        repeats=(2, 2, 1),
        show_unit_cells=True,
        cell_overlays=(
            CrystalCellOverlay(
                kind="hexagonal_prism",
                anchor_fractional=np.array([1.0, 1.0, 0.0]),
                color="#0f766e",
                alpha=0.95,
                linewidth=1.2,
                show_faces=True,
                face_alpha=0.08,
            ),
        ),
        plane_overlays=(
            CrystalPlaneOverlay(
                plane=CrystalPlane.from_miller_bravais((1, 1, -2, 1), phase=zr_hcp),
                label_indices=(1, 1, -2, 1),
                color="#ef4444",
                alpha=0.20,
                annotation_style=PlaneAnnotationStyle(fontsize=10.0),
            ),
            CrystalPlaneOverlay(
                plane=CrystalPlane.from_miller_bravais((0, 0, 0, 1), phase=zr_hcp),
                label_indices=(0, 0, 0, 1),
                color="#0ea5e9",
                alpha=0.16,
                annotation_style=PlaneAnnotationStyle(fontsize=10.0),
            ),
        ),
        direction_overlays=(
            CrystalDirectionOverlay(
                direction=CrystalDirection.from_miller_bravais((2, -1, -1, 0), phase=zr_hcp),
                anchor_fractional=np.array([0.0, 0.0, 0.0]),
                label_indices=(2, -1, -1, 0),
                color="#f59e0b",
                annotation_style=DirectionAnnotationStyle(fontsize=10.0),
            ),
        ),
        theme="journal",
    )
    crystal_figure = plot_crystal_structure_3d(
        zr_scene,
        projection="persp",
        elev_deg=22.0,
        azim_deg=34.0,
    )

    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=Rotation.identity(),
                crystal_frame=ni_fcc.crystal_frame,
                specimen_frame=specimen,
                symmetry=ni_fcc.symmetry,
                phase=ni_fcc,
            ),
            Orientation(
                rotation=Rotation.from_bunge_euler(35.0, 25.0, 10.0),
                crystal_frame=ni_fcc.crystal_frame,
                specimen_frame=specimen,
                symmetry=ni_fcc.symmetry,
                phase=ni_fcc,
            ),
            Orientation(
                rotation=Rotation.from_bunge_euler(70.0, 35.0, 40.0),
                crystal_frame=ni_fcc.crystal_frame,
                specimen_frame=specimen,
                symmetry=ni_fcc.symmetry,
                phase=ni_fcc,
            ),
        ]
    )
    ipf = InversePoleFigure.from_orientations(
        orientations,
        np.array([0.0, 0.0, 1.0]),
        reduce_by_symmetry=True,
        antipodal=True,
    )
    ipf_figure = plot_inverse_pole_figure(ipf)

    direction_figure = plot_crystal_directions(
        (
            CrystalDirection(np.array([1.0, 0.0, 0.0]), phase=ni_fcc),
            CrystalDirection(np.array([1.0, 1.0, 1.0]), phase=ni_fcc),
            CrystalDirection(np.array([1.0, 1.0, 0.0]), phase=ni_fcc),
        ),
        labels=((1, 0, 0), (1, 1, 1), (1, 1, 0)),
        theme="journal",
        title="FCC Direction Stereonet",
    )
    plane_figure = plot_crystal_planes(
        (
            CrystalPlane(MillerIndex([1, 1, 1], phase=ni_fcc), phase=ni_fcc),
            CrystalPlane(MillerIndex([1, 0, 0], phase=ni_fcc), phase=ni_fcc),
        ),
        labels=((1, 1, 1), (1, 0, 0)),
        render="both",
        theme="journal",
        title="FCC Plane Traces",
    )
    symmetry_figure = plot_symmetry_elements(
        ni_fcc.symmetry,
        annotate_axes=True,
        theme="journal",
        title="FCC Symmetry Elements",
    )

    return {
        "xrd_ni_fcc_journal": xrd_figure,
        "saed_ni_fcc_dark": saed_figure,
        "crystal_zr_hcp_journal": crystal_figure,
        "ipf_ni_fcc_journal": ipf_figure,
        "stereonet_directions_ni_fcc_journal": direction_figure,
        "stereonet_planes_ni_fcc_journal": plane_figure,
        "symmetry_elements_ni_fcc_journal": symmetry_figure,
    }
