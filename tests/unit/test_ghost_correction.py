"""Ghost correction: the odd part a pole figure cannot measure.

The premise of every test here is the same physical statement. A diffraction
pole figure obeys Friedel's law, so it cannot distinguish a plane normal from
its opposite; the forward operator therefore annihilates every odd-degree
harmonic, and an inversion returns the even part alone. That even part is not
the texture: it carries false maxima where the specimen is empty and depressed
maxima where it is not. Positivity is the one piece of information the
experiment did not supply and physics does, and it is what the correction
spends.

The headline test is quantitative and uses a known answer: a texture is built,
its pole figures are computed, the even-only inversion is compared with the
truth, and the correction must move measurably *towards* it.
"""

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
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)
from pytex.core.orientation import OrientationSet
from pytex.diffraction.stereonets import spherical_angles_to_directions
from pytex.texture import (
    ODF,
    GhostCorrectionSpec,
    HarmonicODF,
    HarmonicODFReconstructionReport,
    KernelSpec,
    PoleFigure,
    correct_ghosts,
    random_pole_density,
)
from pytex.texture.harmonics import _symmetry_projected_raw_basis
from pytex.texture.models import _pole_density_response_matrix

#: Truth halfwidth, expansion bandlimit and quadrature step for the worked case.
#: The truth is deliberately broad: a texture that a degree-4 expansion can
#: represent isolates the ghost problem from the truncation problem, which is a
#: different defect with a different remedy.
TRUTH_HALFWIDTH_DEG = 50.0
BANDLIMIT = 4
QUADRATURE_STEP_DEG = 15.0
#: The response kernel used when modelling the measurement. It is much sharper
#: than the texture, so the forward operator does not itself destroy the detail
#: the expansion is meant to recover.
RESPONSE_HALFWIDTH_DEG = 10.0


def _orthorhombic_context() -> tuple[ReferenceFrame, ReferenceFrame, Phase]:
    """An orthorhombic phase, whose proper point group admits odd harmonics.

    Symmetry choice is not incidental. The rotation group 432 of a cubic
    material admits no odd-degree invariant below degree 9, so a cubic material
    expanded to a modest bandlimit has no ghost problem to correct. 222 admits
    odd terms from degree 3, which is what makes the effect visible at all.
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
    symmetry = SymmetrySpec.from_point_group("222", reference_frame=crystal)
    lattice = Lattice(3.0, 4.0, 5.0, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name="ortho-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return crystal, specimen, phase


def _measurement_directions() -> np.ndarray:
    azimuth = np.arange(0.0, 360.0, 15.0, dtype=np.float64)
    polar = np.arange(0.0, 91.0, 15.0, dtype=np.float64)
    polar_grid, azimuth_grid = np.meshgrid(polar, azimuth, indexing="ij")
    return spherical_angles_to_directions(polar_grid, azimuth_grid).reshape(-1, 3)


@pytest.fixture(scope="module")
def worked_case() -> tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport]:
    """A known texture, its pole figures, and the even-only inversion of them."""

    crystal, specimen, phase = _orthorhombic_context()
    orientations = OrientationSet.from_euler_angles(
        np.array([[35.0, 45.0, 20.0]]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
        convention="bunge",
        degrees=True,
    )
    truth = ODF(
        orientations=orientations,
        weights=np.array([1.0]),
        kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=TRUTH_HALFWIDTH_DEG),
    )
    directions = _measurement_directions()
    # A measurement is in multiples of a random distribution; the discrete ODF
    # returns a raw kernel response, so the caller divides by the response a
    # random texture would give under the same antipodal convention.
    scale = random_pole_density(truth.kernel, antipodal=True)
    pole_figures = tuple(
        PoleFigure(
            pole=CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase),
            sample_directions=directions,
            intensities=truth.evaluate_pole_density(
                CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase),
                directions,
                antipodal=True,
            )
            / scale,
            specimen_frame=specimen,
            antipodal=True,
            includes_symmetry_family=True,
            sampling="sampled_density",
        )
        for indices in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 1, 1])
    )
    report = HarmonicODF.invert_pole_figures(
        pole_figures,
        degree_bandlimit=BANDLIMIT,
        regularization=1e-6,
        pole_kernel=KernelSpec(
            name="de_la_vallee_poussin", halfwidth_deg=RESPONSE_HALFWIDTH_DEG
        ),
        phi1_step_deg=QUADRATURE_STEP_DEG,
        big_phi_step_deg=QUADRATURE_STEP_DEG,
        phi2_step_deg=QUADRATURE_STEP_DEG,
    )
    return truth, pole_figures, report


def _truth_on_quadrature(truth: ODF, odf: HarmonicODF) -> np.ndarray:
    """The true density on the ODF's quadrature, normalized to unit mean.

    The discrete and harmonic representations normalize over different domains,
    so the comparison is made after rescaling the truth to the mean density of 1
    that the harmonic convention fixes. Only the *shape* is under test.
    """

    density = np.asarray(truth.evaluate(odf.quadrature_orientations, normalized=True))
    return density / float(np.sum(odf.quadrature_weights * density))


def _weighted_distance(left: np.ndarray, right: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(np.sum(weights * (left - right) ** 2)))


def test_odd_harmonics_are_invisible_to_a_friedel_symmetric_pole_figure() -> None:
    """The premise of the whole module, checked directly on the operator.

    An odd-degree basis function must produce no pole density at all once
    opposite normals are identified — that is why pole-figure data cannot
    determine it, and why a correction is entitled to add one without touching
    the fit. Without the folding the same function produces a signal two orders
    of magnitude larger, which is the defect this check would have caught: a
    forward model that ignores the antipodal convention its own pole figure
    declares appears to measure something no diffraction experiment can.
    """

    crystal, specimen, phase = _orthorhombic_context()
    quadrature = OrientationSet.from_bunge_grid(
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=phase.symmetry,
        phase=phase,
        phi1_range_deg=(0.0, 360.0),
        big_phi_range_deg=(0.0, 180.0),
        phi2_range_deg=(0.0, 360.0),
        phi1_step_deg=15.0,
        big_phi_step_deg=15.0,
        phi2_step_deg=15.0,
    )
    weights = np.full(len(quadrature), 1.0 / len(quadrature))
    weights = weights * np.sin(
        np.deg2rad(quadrature.as_euler_set(convention="bunge", degrees=True).angles[:, 1])
    )
    weights = weights / np.sum(weights)
    from pytex.texture.ghosts import _odd_terms

    odd_values = _symmetry_projected_raw_basis(
        quadrature,
        terms=_odd_terms(3),
        crystal_symmetry=phase.symmetry,
        specimen_symmetry=None,
    )
    # The single odd raw column with the largest spread over the quadrature is
    # the strongest candidate for a visible signal.
    column = odd_values[:, int(np.argmax(np.std(odd_values, axis=0)))]
    column = column / float(np.max(np.abs(column)))
    pole = CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase)
    directions = _measurement_directions()
    kernel = KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=RESPONSE_HALFWIDTH_DEG)
    folded = _pole_density_response_matrix(
        quadrature,
        pole=pole,
        sample_directions=directions,
        kernel=kernel,
        include_symmetry_family=True,
        antipodal=True,
    ) @ (weights * column)
    unfolded = _pole_density_response_matrix(
        quadrature,
        pole=pole,
        sample_directions=directions,
        kernel=kernel,
        include_symmetry_family=True,
        antipodal=False,
    ) @ (weights * column)
    folded_level = float(np.max(np.abs(folded))) / random_pole_density(kernel, antipodal=True)
    unfolded_level = float(np.max(np.abs(unfolded))) / random_pole_density(kernel)
    assert folded_level < 0.01
    assert unfolded_level > 20.0 * folded_level


def test_ghost_correction_moves_the_solution_towards_the_true_texture(
    worked_case: tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport],
) -> None:
    """The headline claim, against a texture whose answer is known by construction.

    The even-only inversion of a texture built from a single broad component
    reports a density that is negative over ~9% of orientation space and a
    maximum well below the true one. Both are the classical ghost signature.
    Correcting must reduce the distance to the true distribution substantially,
    not merely make the picture prettier.
    """

    truth, pole_figures, report = worked_case
    odf = report.odf
    weights = odf.quadrature_weights
    true_density = _truth_on_quadrature(truth, odf)

    # The even-only solution shows the ghost problem in the first place.
    assert report.minimum_density < -0.3
    assert report.negative_density_fraction > 0.05
    assert report.maximum_density < float(np.max(true_density))

    corrected = correct_ghosts(odf, pole_figures=pole_figures)
    before = _weighted_distance(odf.quadrature_densities, true_density, weights)
    after = _weighted_distance(corrected.odf.quadrature_densities, true_density, weights)
    assert after < 0.6 * before
    # The depressed maximum is restored towards the true one.
    assert corrected.maximum_density_after > report.maximum_density
    assert corrected.maximum_density_after == pytest.approx(float(np.max(true_density)), rel=0.15)
    # And the false negative lobes are gone.
    assert corrected.negative_density_fraction_after < 0.01
    assert corrected.minimum_density_after > -0.01


def test_ghost_correction_does_not_disturb_the_measured_fit(
    worked_case: tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport],
) -> None:
    """A correction that changed the fit would be spending data agreement.

    The odd part is invisible to the measurement, so adding one must leave every
    predicted pole density where it was. The tolerance is the quadrature error,
    not zero: the odd basis is orthonormalized on a discrete grid, so its
    orthogonality to the even part is exact only in the continuum.
    """

    _, pole_figures, report = worked_case
    corrected = correct_ghosts(report.odf, pole_figures=pole_figures)
    observed = np.concatenate([figure.intensities for figure in pole_figures])
    assert corrected.pole_figure_max_change is not None
    assert corrected.pole_figure_max_change < 0.01 * float(np.max(observed))


def test_ghost_correction_preserves_the_mean_density(
    worked_case: tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport],
) -> None:
    """Odd-degree harmonics integrate to zero, so the normalization cannot move."""

    _, _, report = worked_case
    corrected = correct_ghosts(report.odf)
    assert corrected.mean_density_after == pytest.approx(corrected.mean_density_before, abs=1e-12)
    assert corrected.odf.mean_density == pytest.approx(report.odf.mean_density, abs=1e-12)


def test_ghost_correction_reduces_the_violation_it_minimizes(
    worked_case: tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport],
) -> None:
    """Positivity is the objective, so it must improve and cannot get worse."""

    _, _, report = worked_case
    corrected = correct_ghosts(report.odf)
    assert corrected.infeasibility_before > 0.0
    assert corrected.infeasibility_after < 1e-3 * corrected.infeasibility_before
    assert corrected.converged


def test_zero_range_holds_the_declared_empty_range_nearer_zero(
    worked_case: tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport],
) -> None:
    """The zero-range method spends the odd part on emptying the empty range.

    Positivity alone is satisfied by any non-negative density, including one
    that leaves a broad low plateau where the measurement says there is nothing.
    Declaring a zero range asks for more, and the density inside the declared
    range must end up lower than positivity alone leaves it.
    """

    _, _, report = worked_case
    odf = report.odf
    threshold = 0.15
    inside = np.asarray(odf.quadrature_densities) < threshold
    assert np.any(inside)
    positivity = correct_ghosts(odf, spec=GhostCorrectionSpec(method="positivity"))
    zero_range = correct_ghosts(
        odf,
        spec=GhostCorrectionSpec(method="zero_range", zero_range_threshold=threshold),
    )
    weights = odf.quadrature_weights[inside]
    positivity_mass = float(np.sum(weights * np.abs(positivity.odf.quadrature_densities[inside])))
    zero_range_mass = float(np.sum(weights * np.abs(zero_range.odf.quadrature_densities[inside])))
    assert zero_range_mass < positivity_mass
    assert zero_range.zero_range_fraction == pytest.approx(float(np.mean(inside)))


def test_cubic_symmetry_admits_no_odd_terms_at_a_modest_bandlimit() -> None:
    """For a cubic material below degree 9 there is no ghost part to correct.

    The rotation group 432 has no odd-degree invariant until degree 9, so a
    cubic ODF expanded to degree 6 is already as complete as the symmetry
    allows. The correction must say so rather than pretend to have done work.
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
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name="fcc-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    orientations = OrientationSet.from_euler_angles(
        np.array([[0.0, 0.0, 0.0], [35.0, 45.0, 0.0]]),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    )
    directions = _measurement_directions()
    truth = ODF(
        orientations=orientations,
        weights=np.array([1.0, 1.0]),
        kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=30.0),
    )
    scale = random_pole_density(truth.kernel, antipodal=True)
    pole_figures = tuple(
        PoleFigure(
            pole=CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase),
            sample_directions=directions,
            intensities=truth.evaluate_pole_density(
                CrystalPlane(miller=MillerIndex(indices, phase=phase), phase=phase),
                directions,
                antipodal=True,
            )
            / scale,
            specimen_frame=specimen,
            antipodal=True,
            includes_symmetry_family=True,
            sampling="sampled_density",
        )
        for indices in ([1, 1, 1], [2, 0, 0])
    )
    report = HarmonicODF.invert_pole_figures(
        pole_figures,
        degree_bandlimit=6,
        regularization=1e-4,
        pole_kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=15.0),
        phi1_step_deg=20.0,
        big_phi_step_deg=20.0,
        phi2_step_deg=20.0,
        ghost_correction=True,
    )
    correction = report.ghost_correction
    assert correction is not None
    assert correction.odd_basis_size == 0
    assert correction.ghost_amplitude_ratio == 0.0
    assert_allclose(
        correction.odf.quadrature_densities,
        report.odf.quadrature_densities,
        atol=1e-12,
    )
    assert "no odd-degree harmonic terms" in correction.describe()


def test_inversion_attaches_the_correction_and_offers_the_corrected_odf(
    worked_case: tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport],
) -> None:
    """``ghost_correction`` must change what downstream code reads.

    A correction that only produced a report nobody consults would leave every
    quoted number — texture index, maximum density, Kearns parameter — computed
    from the uncorrected ODF. ``final_odf`` is the surface that makes the
    request effective, while ``odf`` keeps describing the data alone.
    """

    _, pole_figures, _ = worked_case
    corrected_report = HarmonicODF.invert_pole_figures(
        pole_figures,
        degree_bandlimit=BANDLIMIT,
        regularization=1e-6,
        pole_kernel=KernelSpec(
            name="de_la_vallee_poussin", halfwidth_deg=RESPONSE_HALFWIDTH_DEG
        ),
        phi1_step_deg=QUADRATURE_STEP_DEG,
        big_phi_step_deg=QUADRATURE_STEP_DEG,
        phi2_step_deg=QUADRATURE_STEP_DEG,
        ghost_correction=True,
    )
    assert corrected_report.ghost_correction is not None
    assert corrected_report.odf.even_degrees_only is True
    assert corrected_report.final_odf.even_degrees_only is False
    assert corrected_report.final_odf is corrected_report.ghost_correction.odf
    assert float(np.min(corrected_report.final_odf.quadrature_densities)) > -0.01

    plain_report = HarmonicODF.invert_pole_figures(
        pole_figures,
        degree_bandlimit=BANDLIMIT,
        regularization=1e-6,
        pole_kernel=KernelSpec(
            name="de_la_vallee_poussin", halfwidth_deg=RESPONSE_HALFWIDTH_DEG
        ),
        phi1_step_deg=QUADRATURE_STEP_DEG,
        big_phi_step_deg=QUADRATURE_STEP_DEG,
        phi2_step_deg=QUADRATURE_STEP_DEG,
    )
    assert plain_report.ghost_correction is None
    assert plain_report.final_odf is plain_report.odf


def test_ghost_correction_refuses_an_odf_that_already_carries_odd_degrees(
    worked_case: tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport],
) -> None:
    """Correcting a corrected ODF would discard what determined its odd part."""

    _, _, report = worked_case
    corrected = correct_ghosts(report.odf)
    with pytest.raises(ValueError, match="already carries odd degrees"):
        correct_ghosts(corrected.odf)


def test_describe_states_the_assumption_and_the_cost(
    worked_case: tuple[ODF, tuple[PoleFigure, ...], HarmonicODFReconstructionReport],
) -> None:
    """The prose must name the inference, not just the numbers."""

    _, pole_figures, report = worked_case
    prose = correct_ghosts(report.odf, pole_figures=pole_figures).describe()
    assert "positivity alone" in prose
    assert "ghost amplitude ratio" in prose
    assert "inference" in prose
    assert "no pole-figure experiment can confirm or refute it" in prose


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"method": "wishful"}, "must be one of"),
        ({"zero_range_threshold": -0.1}, "non-negative"),
        ({"max_iterations": 0}, "strictly positive"),
        ({"tolerance": 0.0}, "strictly positive"),
        ({"odd_regularization": -1.0}, "non-negative"),
        ({"degree_bandlimit": -1}, "non-negative"),
        ({"basis_tolerance": 0.0}, "strictly positive"),
    ],
)
def test_spec_rejects_settings_that_cannot_mean_anything(
    kwargs: dict[str, object], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        GhostCorrectionSpec(**kwargs)  # type: ignore[arg-type]
