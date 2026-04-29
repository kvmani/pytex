from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose

from pytex.core import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    Orientation,
    Phase,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)
from pytex.diffraction.stereonets import spherical_angles_to_directions
from pytex.plotting.builders import build_odf_figure_spec
from pytex.texture import HarmonicODF, KernelSpec, residual_reports_for_pole_figures
from pytex.texture.harmonics import (
    _bunge_quadrature,
    _enumerate_terms,
    _orthonormalize_weighted_basis,
    _symmetry_projected_raw_basis,
)


def make_harmonic_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase, SymmetrySpec]:
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
    crystal_symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    specimen_symmetry = SymmetrySpec.from_point_group("mmm", reference_frame=specimen)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(
        name="fcc-demo",
        lattice=lattice,
        symmetry=crystal_symmetry,
        crystal_frame=crystal,
    )
    return crystal, specimen, phase, specimen_symmetry


def make_synthetic_harmonic_odf() -> HarmonicODF:
    crystal, specimen, phase, specimen_symmetry = make_harmonic_context()
    quadrature_orientations, quadrature_weights = _bunge_quadrature(
        crystal_frame=crystal,
        specimen_frame=specimen,
        crystal_symmetry=phase.symmetry,
        phase=phase,
        phi1_step_deg=60.0,
        big_phi_step_deg=60.0,
        phi2_step_deg=60.0,
        provenance=None,
    )
    basis_terms = _enumerate_terms(degree_bandlimit=2, even_degrees_only=True)
    raw_basis = _symmetry_projected_raw_basis(
        quadrature_orientations,
        terms=basis_terms,
        crystal_symmetry=phase.symmetry,
        specimen_symmetry=specimen_symmetry,
    )
    quadrature_basis_values, basis_transform = _orthonormalize_weighted_basis(
        raw_basis,
        quadrature_weights,
        tolerance=1e-10,
    )
    coefficients = np.zeros(quadrature_basis_values.shape[1], dtype=np.float64)
    coefficients[: min(4, coefficients.size)] = np.array([1.0, 0.35, -0.2, 0.1], dtype=np.float64)[
        : min(4, coefficients.size)
    ]
    return HarmonicODF(
        coefficients=coefficients,
        basis_terms=basis_terms,
        basis_transform=basis_transform,
        quadrature_orientations=quadrature_orientations,
        quadrature_weights=quadrature_weights,
        quadrature_basis_values=quadrature_basis_values,
        degree_bandlimit=2,
        crystal_symmetry=phase.symmetry,
        specimen_symmetry=specimen_symmetry,
        phase=phase,
        pole_kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=7.5),
        even_degrees_only=True,
        provenance=None,
    )


def make_measurement_grid() -> np.ndarray:
    phi_values = np.arange(0.0, 360.0, 60.0, dtype=np.float64)
    psi_values = np.arange(0.0, 91.0, 30.0, dtype=np.float64)
    psi_grid, phi_grid = np.meshgrid(psi_values, phi_values, indexing="ij")
    return spherical_angles_to_directions(psi_grid, phi_grid).reshape(-1, 3)


def test_harmonic_odf_is_invariant_under_crystal_and_sample_symmetry_actions() -> None:
    harmonic_odf = make_synthetic_harmonic_odf()
    orientation = Orientation(
        rotation=Rotation.from_bunge_euler(35.0, 40.0, 15.0),
        crystal_frame=harmonic_odf.crystal_frame,
        specimen_frame=harmonic_odf.specimen_frame,
        symmetry=harmonic_odf.crystal_symmetry,
        phase=harmonic_odf.phase,
    )
    base_value = float(harmonic_odf.evaluate(orientation))
    matrix = orientation.as_matrix()
    sample_operator = harmonic_odf.specimen_symmetry.operators[1]
    crystal_operator = harmonic_odf.crystal_symmetry.operators[1]
    equivalent = Orientation(
        rotation=Rotation.from_matrix(sample_operator @ matrix @ crystal_operator),
        crystal_frame=harmonic_odf.crystal_frame,
        specimen_frame=harmonic_odf.specimen_frame,
        symmetry=harmonic_odf.crystal_symmetry,
        phase=harmonic_odf.phase,
    )
    assert_allclose(harmonic_odf.evaluate(equivalent), base_value, atol=1e-10)


def test_harmonic_odf_inversion_recovers_synthetic_pole_density_response() -> None:
    harmonic_odf = make_synthetic_harmonic_odf()
    sample_directions = make_measurement_grid()
    phase = harmonic_odf.phase
    assert phase is not None
    poles = (
        CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase),
        CrystalPlane(miller=MillerIndex([1, 1, 0], phase=phase), phase=phase),
        CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase),
    )
    pole_figures = tuple(
        harmonic_odf.reconstruct_pole_figure(
            pole,
            sample_directions=sample_directions,
            include_symmetry_family=True,
            antipodal=True,
        )
        for pole in poles
    )
    report = HarmonicODF.invert_pole_figures(
        pole_figures,
        degree_bandlimit=2,
        regularization=1e-8,
        include_symmetry_family=True,
        pole_kernel=harmonic_odf.pole_kernel,
        specimen_symmetry=harmonic_odf.specimen_symmetry,
        phi1_step_deg=60.0,
        big_phi_step_deg=60.0,
        phi2_step_deg=60.0,
        basis_tolerance=1e-10,
    )
    observations = np.concatenate([pole_figure.intensities for pole_figure in pole_figures])
    assert report.even_degrees_only is True
    assert all(term.degree % 2 == 0 for term in report.odf.basis_terms)
    assert report.relative_residual_norm < 1e-6
    assert_allclose(report.predicted_intensities, observations, atol=5e-8)
    query_orientations = report.odf.quadrature_orientations
    assert_allclose(
        report.odf.evaluate(query_orientations),
        harmonic_odf.evaluate(query_orientations),
        atol=1e-5,
    )


def test_harmonic_odf_section_plot_builder_returns_panel_grid() -> None:
    harmonic_odf = make_synthetic_harmonic_odf()
    spec = build_odf_figure_spec(
        harmonic_odf,
        kind="sections",
        section_phi2_deg=(0.0, 45.0),
        section_phi1_steps=25,
        section_big_phi_steps=13,
    )
    assert len(spec.panels) == 2


def test_harmonic_reconstruction_report_exposes_basis_and_density_diagnostics() -> None:
    harmonic_odf = make_synthetic_harmonic_odf()
    sample_directions = make_measurement_grid()
    phase = harmonic_odf.phase
    assert phase is not None
    pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase)
    report = HarmonicODF.invert_pole_figures(
        (
            harmonic_odf.reconstruct_pole_figure(
                pole,
                sample_directions=sample_directions,
                include_symmetry_family=True,
                antipodal=True,
            ),
        ),
        degree_bandlimit=2,
        regularization=1e-8,
        include_symmetry_family=True,
        pole_kernel=harmonic_odf.pole_kernel,
        specimen_symmetry=harmonic_odf.specimen_symmetry,
        phi1_step_deg=60.0,
        big_phi_step_deg=60.0,
        phi2_step_deg=60.0,
        basis_tolerance=1e-10,
    )
    assert report.basis_size > 0
    assert report.raw_basis_size >= report.basis_size
    assert report.matrix_rank > 0
    assert report.quadrature_size == len(report.odf.quadrature_orientations)
    assert report.crystal_symmetry_order == harmonic_odf.crystal_symmetry.order
    assert report.specimen_symmetry_order == harmonic_odf.specimen_symmetry.order
    assert report.coefficient_l2_norm >= 0.0
    assert report.coefficient_max_abs >= 0.0
    assert 0.0 <= report.negative_density_fraction <= 1.0


def test_harmonic_residual_reports_match_reconstruction_surface() -> None:
    harmonic_odf = make_synthetic_harmonic_odf()
    sample_directions = make_measurement_grid()
    phase = harmonic_odf.phase
    assert phase is not None
    poles = (
        CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase),
        CrystalPlane(miller=MillerIndex([1, 1, 0], phase=phase), phase=phase),
    )
    pole_figures = tuple(
        harmonic_odf.reconstruct_pole_figure(
            pole,
            sample_directions=sample_directions,
            include_symmetry_family=True,
            antipodal=True,
        )
        for pole in poles
    )
    reports = residual_reports_for_pole_figures(harmonic_odf, pole_figures)
    assert len(reports) == 2
    for report, pole_figure in zip(reports, pole_figures, strict=True):
        assert report.observation_count == pole_figure.intensities.size
        assert_allclose(report.predicted_intensities, pole_figure.intensities, atol=1e-8)
        assert_allclose(report.residuals, np.zeros_like(report.residuals), atol=1e-8)
        assert report.relative_residual_norm <= 1e-8
