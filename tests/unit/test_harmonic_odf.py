from __future__ import annotations

import math

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
    Phase,
    ReferenceFrame,
    Rotation,
    S2Grid,
    SymmetrySpec,
)
from pytex.core.orientation import OrientationSet
from pytex.diffraction.stereonets import spherical_angles_to_directions
from pytex.plotting.builders import build_odf_figure_spec
from pytex.texture import (
    HarmonicODF,
    KernelSpec,
    PoleFigure,
    residual_reports_for_pole_figures,
)
from pytex.texture.harmonics import (
    _bunge_quadrature,
    _enumerate_terms,
    _orthonormalize_weighted_basis,
    _symmetry_projected_raw_basis,
    _weighted_mean,
    _wigner_small_d,
)


def _build_harmonic_odf(
    perturbations: dict[int, float],
    *,
    step_deg: float = 30.0,
    halfwidth_deg: float = 10.0,
) -> HarmonicODF:
    """Build a valid (mean-density-1) HarmonicODF with identity symmetry.

    ``perturbations`` maps non-constant orthonormal basis columns to their
    coefficient; the constant column is set to 1 so the ODF integrates to 1.
    ``step_deg`` sets the Bunge quadrature spacing and ``halfwidth_deg`` the pole
    kernel; the defaults are deliberately coarse and fast, and a test that needs
    the pole-density quadrature to be accurate must refine both.
    """

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
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name="demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    quadrature_orientations, quadrature_weights = _bunge_quadrature(
        crystal_frame=crystal,
        specimen_frame=specimen,
        crystal_symmetry=symmetry,
        phase=phase,
        phi1_step_deg=step_deg,
        big_phi_step_deg=step_deg,
        phi2_step_deg=step_deg,
        provenance=None,
    )
    basis_terms = _enumerate_terms(degree_bandlimit=2, even_degrees_only=False)
    raw_basis = _symmetry_projected_raw_basis(
        quadrature_orientations,
        terms=basis_terms,
        crystal_symmetry=symmetry,
        specimen_symmetry=None,
    )
    quadrature_basis_values, basis_transform = _orthonormalize_weighted_basis(
        raw_basis, quadrature_weights, tolerance=1e-10
    )
    means = np.array(
        [
            _weighted_mean(quadrature_basis_values[:, k], quadrature_weights)
            for k in range(quadrature_basis_values.shape[1])
        ]
    )
    constant_column = int(np.argmax(np.abs(means)))
    coefficients = np.zeros(quadrature_basis_values.shape[1], dtype=np.float64)
    coefficients[constant_column] = 1.0
    for column, value in perturbations.items():
        coefficients[column] = value
    return HarmonicODF(
        coefficients=coefficients,
        basis_terms=basis_terms,
        basis_transform=basis_transform,
        quadrature_orientations=quadrature_orientations,
        quadrature_weights=quadrature_weights,
        quadrature_basis_values=quadrature_basis_values,
        degree_bandlimit=2,
        crystal_symmetry=symmetry,
        specimen_symmetry=None,
        phase=phase,
        pole_kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=halfwidth_deg),
        even_degrees_only=False,
        provenance=None,
    )


def test_uniform_harmonic_odf_has_unit_texture_index_and_zero_entropy() -> None:
    odf = _build_harmonic_odf({})
    # small deviations from exactly 1 are coarse-quadrature artefacts
    assert odf.mean_density == pytest.approx(1.0, abs=1e-3)
    assert odf.texture_index == pytest.approx(1.0, abs=1e-3)
    assert odf.entropy() == pytest.approx(0.0, abs=1e-6)


def test_textured_harmonic_odf_index_and_entropy_grow_with_sharpness() -> None:
    # Perturb two non-constant orthonormal columns; the texture index equals
    # 1 + sum of squared perturbations (orthonormal basis).
    mild = _build_harmonic_odf({0: 0.3, 1: -0.2})
    sharp = _build_harmonic_odf({0: 0.6, 1: -0.4})
    assert mild.mean_density == pytest.approx(1.0, abs=1e-2)
    assert mild.texture_index == pytest.approx(1.0 + 0.3**2 + 0.2**2, abs=1e-3)
    assert mild.texture_index > 1.0
    assert mild.entropy() > 0.0
    # a sharper texture has a larger index and entropy
    assert sharp.texture_index > mild.texture_index
    assert sharp.entropy() > mild.entropy()


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


def _reference_wigner_small_d(
    degree: int, sample_order: int, crystal_order: int, beta_rad: np.ndarray
) -> np.ndarray:
    """Wigner small-d straight from the factorial formula, in exact integers.

    This is the textbook expression evaluated with Python's arbitrary-precision
    integers and converted to float only at the last moment, so it is exact to
    rounding. It is unusable in the library — the integers grow past what NumPy
    can hold — but it is the right thing to check the log-gamma form against.
    """

    prefactor = math.sqrt(
        float(
            math.factorial(degree + sample_order)
            * math.factorial(degree - sample_order)
            * math.factorial(degree + crystal_order)
            * math.factorial(degree - crystal_order)
        )
    )
    values = np.zeros_like(beta_rad, dtype=np.float64)
    k_min = max(0, crystal_order - sample_order)
    k_max = min(degree - sample_order, degree + crystal_order)
    for k in range(k_min, k_max + 1):
        denominator = float(
            math.factorial(degree + crystal_order - k)
            * math.factorial(k)
            * math.factorial(sample_order - crystal_order + k)
            * math.factorial(degree - sample_order - k)
        )
        values += (
            ((-1) ** (k - sample_order + crystal_order))
            * prefactor
            / denominator
            * np.cos(beta_rad / 2.0) ** (2 * degree + crystal_order - sample_order - 2 * k)
            * np.sin(beta_rad / 2.0) ** (sample_order - crystal_order + 2 * k)
        )
    return values


def test_wigner_small_d_matches_the_exact_factorial_form() -> None:
    """The log-gamma evaluation must agree with exact integer arithmetic.

    Degrees 0-6 are the range where the factorial products still fit in the
    integer types NumPy can hold, so the exact form is available to check
    against for every (degree, sample order, crystal order) triple.
    """

    beta = np.linspace(0.05, np.pi - 0.05, 17)
    worst = 0.0
    for degree in range(7):
        for sample_order in range(-degree, degree + 1):
            for crystal_order in range(-degree, degree + 1):
                computed = _wigner_small_d(degree, sample_order, crystal_order, beta)
                expected = _reference_wigner_small_d(degree, sample_order, crystal_order, beta)
                worst = max(worst, float(np.max(np.abs(computed - expected))))
    assert worst < 1e-12


def test_wigner_small_d_survives_degrees_that_overflow_integer_factorials() -> None:
    """Above degree 6 the factorial products exceed int64 and used to crash.

    ``factorial(2*degree)`` alone reaches 8.7e18 at degree 7, and the formula
    multiplies four such terms; the result became a Python big integer that
    NumPy could hold only as an object, and ``np.sqrt`` on an object array
    raises ``TypeError``. The bandlimit default is 6, so this surfaced only when
    a user raised the bandwidth for a sharp texture — exactly when they need it.
    """

    beta = np.linspace(0.05, np.pi - 0.05, 9)
    for degree in (7, 10, 16, 24):
        values = _wigner_small_d(degree, degree, -degree, beta)
        assert np.all(np.isfinite(values))
    # d^l_{00}(pi/2) has the closed form 0 for odd l and, for even l = 2m,
    # (-1)^m * C(2m, m) / 4^m. Checking it pins the normalization, not just
    # finiteness. The tolerance loosens with degree because the sum alternates
    # in sign over l + 1 terms and loses a few digits to cancellation; at degree
    # 16 the agreement is 7e-12 against a value of 0.196.
    half_pi = np.array([np.pi / 2.0])
    assert _wigner_small_d(7, 0, 0, half_pi)[0] == pytest.approx(0.0, abs=1e-12)
    for degree in (8, 10, 16):
        order = degree // 2
        expected = ((-1) ** order) * math.comb(degree, order) / (2.0**degree)
        assert _wigner_small_d(degree, 0, 0, half_pi)[0] == pytest.approx(expected, abs=1e-10)


def test_harmonic_inversion_accepts_a_bandlimit_above_six() -> None:
    """The end-to-end path that the overflow made unreachable."""

    crystal, specimen, phase, _ = make_harmonic_context()
    orientations = OrientationSet.from_euler_angles(
        np.array([[0.0, 0.0, 0.0], [30.0, 20.0, 10.0]]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
    grid = S2Grid.equispaced(
        20.0, reference_frame=specimen, hemisphere="upper", antipodal=True
    )
    measured = PoleFigure.from_orientations(orientations, pole).on_grid(
        grid, halfwidth_deg=20.0
    )
    report = HarmonicODF.invert_pole_figures(
        (measured,), degree_bandlimit=8, regularization=1e-3
    )
    assert np.all(np.isfinite(report.odf.coefficients))


def test_uniform_harmonic_odf_pole_density_is_exactly_one_mrd() -> None:
    """A uniform ODF must give unit pole density, on the m.r.d. scale.

    This is an exact identity rather than a tolerance against a prior run: if
    every orientation is equally likely then every specimen direction is equally
    likely to carry a given pole, so the pole figure is flat at 1.0 multiples of
    random. It is the check that pins the scale of
    :meth:`HarmonicODF.evaluate_pole_density`, whose kernel-weighted quadrature
    sum would otherwise return the kernel's spherical mean — of order 0.006 —
    for exactly this case.
    """

    # The quadrature has to resolve the kernel for the identity to hold
    # numerically, so both are refined relative to the module default.
    odf = _build_harmonic_odf({}, step_deg=12.0, halfwidth_deg=25.0)
    phase = odf.phase
    assert phase is not None
    pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
    densities = odf.evaluate_pole_density(pole, make_measurement_grid())
    assert_allclose(densities, np.ones_like(densities), atol=5e-3)


def test_harmonic_inversion_returns_an_odf_on_the_mrd_scale() -> None:
    """Inversion of m.r.d. pole figures must return an m.r.d. ODF.

    The forward operator and the returned density have to share one scale. They
    did not: the operator carried the raw kernel response, so the fitted
    coefficients absorbed its reciprocal and every density the ODF reported was
    about 163 times too large at the default halfwidth. The pole-figure round
    trip still closed, which is why the error was invisible from the residual
    alone — only the ODF's own mean density exposes it.
    """

    crystal, specimen, phase, _ = make_harmonic_context()
    orientations = OrientationSet.from_euler_angles(
        np.array([[0.0, 0.0, 0.0], [35.0, 45.0, 0.0], [55.0, 30.0, 65.0]]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
    )
    grid = S2Grid.equispaced(15.0, reference_frame=specimen, hemisphere="upper")
    pole_figures = tuple(
        PoleFigure.from_orientations(
            orientations, CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase)
        ).on_grid(grid, halfwidth_deg=15.0)
        for indices in ([1, 1, 1], [2, 0, 0], [2, 2, 0])
    )
    # Every input figure is normalized to multiples of random, so its spherical
    # mean is 1 and the recovered ODF's mean density must be too.
    for pole_figure in pole_figures:
        assert pole_figure.spherical_mean() == pytest.approx(1.0, abs=0.1)
    report = HarmonicODF.invert_pole_figures(pole_figures, degree_bandlimit=6)
    assert report.odf.mean_density == pytest.approx(1.0, abs=0.1)
    assert report.odf.texture_index > 1.0
