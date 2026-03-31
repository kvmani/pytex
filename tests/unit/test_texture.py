from __future__ import annotations

import pytest
from numpy.testing import assert_allclose

from pytex.core import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Orientation,
    OrientationSet,
    Phase,
    ProvenanceRecord,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.texture import ODF, InversePoleFigure, KernelSpec, ODFInversionReport, PoleFigure


def make_orientation_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name="fcc-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def make_orientation_set() -> tuple[OrientationSet, Phase]:
    crystal, specimen, phase = make_orientation_context()
    provenance = ProvenanceRecord.minimal("demo-orientations")
    orientations = [
        Orientation(
            Rotation.identity(),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
            provenance=provenance,
        ),
        Orientation(
            Rotation.from_bunge_euler(90.0, 0.0, 0.0),
            crystal_frame=crystal,
            specimen_frame=specimen,
            symmetry=phase.symmetry,
            phase=phase,
            provenance=provenance,
        ),
    ]
    return OrientationSet.from_orientations(orientations), phase


def test_pole_figure_from_orientations_maps_plane_normal_into_specimen() -> None:
    orientations, phase = make_orientation_set()
    pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase)
    pole_figure = PoleFigure.from_orientations(
        orientations,
        pole,
        include_symmetry_family=False,
        antipodal=False,
    )
    assert pole_figure.sample_directions.shape == (2, 3)
    assert_allclose(pole_figure.sample_directions[0], [1.0, 0.0, 0.0], atol=1e-8)
    assert_allclose(pole_figure.sample_directions[1], [0.0, 1.0, 0.0], atol=1e-8)


def test_inverse_pole_figure_without_symmetry_reduction_preserves_crystal_vectors() -> None:
    orientations, _ = make_orientation_set()
    inverse_pole_figure = InversePoleFigure.from_orientations(
        orientations,
        [1.0, 0.0, 0.0],
        reduce_by_symmetry=False,
        antipodal=False,
    )
    assert inverse_pole_figure.crystal_directions.shape == (2, 3)
    assert_allclose(inverse_pole_figure.crystal_directions[0], [1.0, 0.0, 0.0], atol=1e-8)
    assert_allclose(inverse_pole_figure.crystal_directions[1], [0.0, -1.0, 0.0], atol=1e-8)


def test_inverse_pole_figure_with_symmetry_reduction_uses_class_specific_sector() -> None:
    orientations, _ = make_orientation_set()
    inverse_pole_figure = InversePoleFigure.from_orientations(
        orientations,
        [0.0, 0.0, 1.0],
        reduce_by_symmetry=True,
        antipodal=True,
    )
    assert inverse_pole_figure.sector_vertices is not None
    assert inverse_pole_figure.project_sector_vertices() is not None
    for direction in inverse_pole_figure.crystal_directions:
        assert orientations.symmetry is not None
        assert orientations.symmetry.vector_in_fundamental_sector(direction, antipodal=True)


def test_odf_evaluation_is_larger_near_support_orientation() -> None:
    orientations, _ = make_orientation_set()
    odf = ODF.from_orientations(
        orientations,
        kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=15.0),
    )
    near = Orientation(
        Rotation.identity(),
        crystal_frame=orientations.crystal_frame,
        specimen_frame=orientations.specimen_frame,
        symmetry=orientations.symmetry,
        phase=orientations.phase,
    )
    far = Orientation(
        Rotation.from_bunge_euler(45.0, 45.0, 45.0),
        crystal_frame=orientations.crystal_frame,
        specimen_frame=orientations.specimen_frame,
        symmetry=orientations.symmetry,
        phase=orientations.phase,
    )
    assert odf.evaluate(near) > odf.evaluate(far)


def test_odf_volume_fraction_tracks_weighted_neighborhood() -> None:
    orientations, _ = make_orientation_set()
    odf = ODF.from_orientations(
        orientations,
        weights=[3.0, 1.0],
        kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=20.0),
    )
    center = Orientation(
        Rotation.identity(),
        crystal_frame=orientations.crystal_frame,
        specimen_frame=orientations.specimen_frame,
        symmetry=orientations.symmetry,
        phase=orientations.phase,
    )
    assert_allclose(
        odf.volume_fraction(center, max_angle_deg=10.0, symmetry_aware=False),
        0.75,
        atol=1e-8,
    )


def test_pole_figure_rejects_mismatched_phase() -> None:
    orientations, phase = make_orientation_set()
    crystal = phase.crystal_frame
    lattice = Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    other_phase = Phase(
        name="other",
        lattice=lattice,
        symmetry=phase.symmetry,
        crystal_frame=crystal,
    )
    pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=other_phase), phase=other_phase)
    with pytest.raises(ValueError):
        PoleFigure.from_orientations(orientations, pole)


def test_odf_rejects_specimen_symmetry_with_wrong_reference_frame() -> None:
    orientations, phase = make_orientation_set()
    with pytest.raises(ValueError):
        ODF.from_orientations(
            orientations,
            specimen_symmetry=SymmetrySpec.from_point_group(
                "mmm",
                reference_frame=phase.crystal_frame,
            ),
        )


def test_derived_texture_models_preserve_orientation_set_provenance_by_default() -> None:
    orientations, phase = make_orientation_set()
    pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase)
    pole_figure = PoleFigure.from_orientations(orientations, pole)
    inverse_pole_figure = InversePoleFigure.from_orientations(orientations, [1.0, 0.0, 0.0])
    odf = ODF.from_orientations(orientations)
    assert pole_figure.provenance == orientations.provenance
    assert inverse_pole_figure.provenance == orientations.provenance
    assert odf.provenance == orientations.provenance


def test_odf_inversion_recovers_dictionary_weights_from_synthetic_pole_figures() -> None:
    orientations, phase = make_orientation_set()
    true_odf = ODF.from_orientations(
        orientations,
        weights=[3.0, 1.0],
        kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=8.0),
    )
    poles = (
        CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase),
        CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase),
    )
    pole_figures = true_odf.reconstruct_pole_figures(
        poles,
        include_symmetry_family=False,
        antipodal=False,
    )
    report = ODF.invert_pole_figures(
        pole_figures,
        orientation_dictionary=orientations,
        kernel=true_odf.kernel,
        regularization=1e-8,
        include_symmetry_family=False,
        max_iterations=1000,
        tolerance=1e-10,
    )
    assert isinstance(report, ODFInversionReport)
    assert report.observation_count == sum(
        len(pole_figure.intensities) for pole_figure in pole_figures
    )
    assert report.predicted_intensities.shape == (report.observation_count,)
    assert report.relative_residual_norm >= 0.0
    assert report.mean_absolute_error >= 0.0
    assert report.max_absolute_error >= 0.0
    assert report.dictionary_coverage_ratio == report.observation_count / report.dictionary_size
    assert_allclose(report.odf.normalized_weights, true_odf.normalized_weights, atol=5e-2)


def test_odf_inversion_rejects_mismatched_specimen_frames() -> None:
    orientations, phase = make_orientation_set()
    pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase)
    pole_figure = PoleFigure.from_orientations(
        orientations,
        pole,
        include_symmetry_family=False,
        antipodal=False,
    )
    other_specimen = ReferenceFrame(
        name="other_specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    mismatched_dictionary = OrientationSet.from_orientations(
        [
            Orientation(
                rotation=orientation.rotation,
                crystal_frame=orientation.crystal_frame,
                specimen_frame=other_specimen,
                symmetry=orientation.symmetry,
                phase=orientation.phase,
                provenance=orientation.provenance,
            )
            for orientation in [orientations[0], orientations[1]]
        ]
    )
    with pytest.raises(ValueError):
        ODF.invert_pole_figures(
            [pole_figure],
            orientation_dictionary=mismatched_dictionary,
        )


def test_orientation_set_from_bunge_grid_builds_expected_support_size() -> None:
    crystal, specimen, phase = make_orientation_context()
    dictionary = OrientationSet.from_bunge_grid(
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
        phi1_step_deg=90.0,
        big_phi_step_deg=45.0,
        phi2_step_deg=90.0,
    )
    assert len(dictionary) == 12
    assert dictionary.phase == phase
