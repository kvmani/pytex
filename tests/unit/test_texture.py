from __future__ import annotations

import numpy as np
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


def test_kernel_spec_exposes_dvp_so3_kernel_and_bandwidth() -> None:
    spec = KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=10.0)
    so3 = spec.as_so3_kernel()
    assert so3.halfwidth_deg == 10.0
    # sharper kernels need a higher harmonic bandwidth
    sharp = KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=5.0)
    assert sharp.bandwidth() > spec.bandwidth()
    # the SO(3) kernel is unavailable for the von Mises-Fisher spec
    with pytest.raises(ValueError, match="de la Vallee Poussin"):
        KernelSpec(name="von_mises_fisher", halfwidth_deg=10.0).as_so3_kernel()


def test_normalized_dvp_evaluation_rescales_by_kernel_constant() -> None:
    spec = KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=12.0)
    angles = np.deg2rad([0.0, 5.0, 20.0])
    raw = spec.evaluate(angles)
    normed = spec.evaluate(angles, normalized=True)
    # normalization is a positive multiplicative constant equal to the
    # dedicated SO(3) kernel's normalization factor
    factor = spec.as_so3_kernel().normalization
    assert_allclose(normed, raw * factor)


def test_odf_normalized_evaluation_matches_unit_scaled_density() -> None:
    orientations, _ = make_orientation_set()
    kernel = KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=15.0)
    odf = ODF.from_orientations(orientations, kernel=kernel)
    near = Orientation(
        Rotation.identity(),
        crystal_frame=orientations.crystal_frame,
        specimen_frame=orientations.specimen_frame,
        symmetry=orientations.symmetry,
        phase=orientations.phase,
    )
    factor = kernel.as_so3_kernel().normalization
    assert_allclose(
        odf.evaluate(near, normalized=True),
        odf.evaluate(near) * factor,
    )


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


def test_odf_phi2_sections_peak_at_planted_component() -> None:
    crystal, specimen, phase = make_orientation_context()
    component = OrientationSet.from_quaternions(
        Rotation.from_bunge_euler(35.0, 45.0, 0.0).quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    odf = ODF.from_orientations(
        component, kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=8.0)
    )
    sections = odf.phi2_sections(phi2_deg=[0.0, 45.0], resolution_deg=5.0)
    assert sections.densities.shape == (2, 19, 19)
    assert sections.section_count == 2
    # the phi2 = 0 section must peak at the planted (phi1=35, Phi=45) component
    peak = np.unravel_index(int(np.argmax(sections.densities[0])), sections.densities[0].shape)
    assert sections.phi1_deg[peak[1]] == pytest.approx(35.0)
    assert sections.big_phi_deg[peak[0]] == pytest.approx(45.0)


def test_odf_phi2_sections_are_nonnegative_and_broaden_with_kernel() -> None:
    crystal, specimen, phase = make_orientation_context()
    component = OrientationSet.from_quaternions(
        Rotation.from_bunge_euler(0.0, 45.0, 0.0).quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )

    def peak_to_mean(halfwidth: float) -> float:
        odf = ODF.from_orientations(
            component, kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=halfwidth)
        )
        density = odf.phi2_sections(phi2_deg=[0.0], resolution_deg=10.0).densities
        assert np.all(np.isfinite(density)) and np.all(density >= 0.0)
        return float(np.max(density) / np.mean(density))

    # a broader kernel spreads the texture out, lowering the peak-to-mean ratio
    assert peak_to_mean(30.0) < peak_to_mean(8.0)


def test_plot_odf_phi2_sections_returns_panel() -> None:
    import matplotlib

    matplotlib.use("Agg")
    from pytex import plot_odf_phi2_sections

    crystal, specimen, phase = make_orientation_context()
    component = OrientationSet.from_quaternions(
        Rotation.from_bunge_euler(0.0, 45.0, 0.0).quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    odf = ODF.from_orientations(component)
    sections = odf.phi2_sections(phi2_deg=[0.0, 30.0, 60.0], resolution_deg=15.0)
    figure = plot_odf_phi2_sections(sections)
    # three section panels + one shared colorbar axis
    assert len(figure.axes) == 4


def test_plot_odf_phi2_sections_publication_defaults() -> None:
    import matplotlib

    matplotlib.use("Agg")
    from pytex import plot_odf_phi2_sections

    crystal, specimen, phase = make_orientation_context()
    component = OrientationSet.from_quaternions(
        Rotation.from_bunge_euler(0.0, 45.0, 0.0).quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    odf = ODF.from_orientations(component)
    sections = odf.phi2_sections(phi2_deg=[0.0, 30.0, 60.0, 90.0], resolution_deg=15.0)
    figure = plot_odf_phi2_sections(sections, max_cols=2)
    # 2x2 section grid + one shared colorbar axis
    assert len(figure.axes) == 5
    assert figure.axes[0].get_title() == r"$\varphi_2 = 0^\circ$"
    assert figure.axes[0].get_ylabel() == r"$\Phi$ (deg)"
    assert figure.axes[2].get_xlabel() == r"$\varphi_1$ (deg)"
    # automatic panel labels
    labels = [text.get_text() for axis in figure.axes[:4] for text in axis.texts]
    assert "(a)" in labels and "(d)" in labels
    # shared colorbar in multiples-of-random-density units
    assert figure.axes[-1].get_ylabel() == "ODF density (m.r.d.)"


def test_sigma_sections_peak_at_planted_component() -> None:
    crystal, specimen, phase = make_orientation_context()
    component = OrientationSet.from_quaternions(
        Rotation.from_bunge_euler(20.0, 45.0, 30.0).quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    odf = ODF.from_orientations(component)
    sections = odf.sigma_sections(sigma_deg=[20.0, 50.0], resolution_deg=5.0)
    assert sections.section_kind == "sigma"
    assert sections.section_values_deg.tolist() == [20.0, 50.0]
    # the component at (phi1, Phi, phi2) = (20, 45, 30) lives on the
    # sigma = phi1 + phi2 = 50 section, peaked at (phi1=20, Phi=45)
    density = sections.densities[1]
    peak = np.unravel_index(int(np.argmax(density)), density.shape)
    assert sections.big_phi_deg[peak[0]] == pytest.approx(45.0)
    assert sections.phi1_deg[peak[1]] == pytest.approx(20.0)
    assert density.max() > sections.densities[0].max()


def test_plot_sigma_sections_uses_sigma_labels() -> None:
    import matplotlib

    matplotlib.use("Agg")
    from pytex import plot_odf_phi2_sections

    crystal, specimen, phase = make_orientation_context()
    component = OrientationSet.from_quaternions(
        Rotation.from_bunge_euler(20.0, 45.0, 30.0).quaternion[None, :],
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    sections = ODF.from_orientations(component).sigma_sections(
        sigma_deg=[0.0, 45.0], resolution_deg=15.0
    )
    figure = plot_odf_phi2_sections(sections)
    assert figure.axes[0].get_title() == r"$\sigma = 0^\circ$"
