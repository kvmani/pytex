"""Precise lattice-parameter determination, checked against injected truth.

Every accuracy assertion compares a determined cell with the cell that
generated the pattern, and every systematic-error assertion compares a refined
aberration with the aberration that was deliberately injected. Nothing here is
compared against a stored output of this code.

The central claim the module exists to support is tested directly in
:func:`test_extrapolation_beats_averaging_against_an_injected_displacement`:
with a 100 micrometre specimen displacement, the naive per-reflection average
is wrong by 4 parts in 1e4 and a Cohen determination against the matching
extrapolation function is wrong by 1 part in 1e7.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_corrections import specimen_displacement_shift_deg
from pytex.diffraction.xrd_indexing import index_peaks
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_lattice_parameter import (
    EXTRAPOLATION_FUNCTIONS,
    LATTICE_PARAMETER_SCHEMA,
    crystal_system_of,
    determine_lattice_parameters,
    determine_lattice_parameters_from_pattern,
    determine_lattice_parameters_le_bail,
    extrapolation_values,
    nelson_riley_extrapolation,
)
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import detect_and_fit_peaks

FWHM = 0.12
RADIUS_MM = 240.0


def _scan(
    identifier: str,
    *,
    displacement_mm: float = 0.0,
    zero_deg: float = 0.0,
    peak_counts: float = 30000.0,
    seed: int = 5,
    two_theta_range_deg: tuple[float, float] = (25.0, 150.0),
) -> MeasuredPowderPattern:
    """Return a synthetic scan with a *known* aberration written into its axis."""

    phase = builtin_phase(identifier).to_phase()
    radiation = RadiationSpec.cu_ka_doublet()
    pattern = generate_xrd_pattern(
        phase,
        radiation=radiation,
        two_theta_range_deg=two_theta_range_deg,
        resolution_deg=0.01,
        broadening_fwhm_deg=FWHM,
        profile="pseudo_voigt",
        max_index=6,
    )
    axis = np.asarray(pattern.two_theta_grid_deg, dtype=float)
    profile = np.asarray(pattern.intensity_grid, dtype=float)
    expected = profile / profile.max() * peak_counts + 150.0
    counts = np.random.default_rng(seed).poisson(expected).astype(float)
    shifted = axis + zero_deg
    if displacement_mm != 0.0:
        shifted = shifted + specimen_displacement_shift_deg(
            axis, displacement_mm=displacement_mm, goniometer_radius_mm=RADIUS_MM
        )
    return MeasuredPowderPattern(
        name=f"{identifier} synthetic",
        two_theta_deg=shifted,
        intensity=counts,
        radiation=radiation,
        synthetic=True,
    )


def _determine(identifier: str, measured: MeasuredPowderPattern, **kwargs: object):
    phase = builtin_phase(identifier).to_phase()
    return determine_lattice_parameters_from_pattern(
        measured,
        phase,
        instrument=InstrumentBroadening.ideal(FWHM),
        phase_name=identifier,
        **kwargs,  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# The cell parameterization
# ---------------------------------------------------------------------------


def test_crystal_system_comes_from_symmetry_not_from_the_cell_numbers() -> None:
    assert crystal_system_of(builtin_phase("ni_fcc").to_phase()) == "cubic"
    assert crystal_system_of(builtin_phase("ti_hcp").to_phase()) == "hexagonal"
    assert crystal_system_of(builtin_phase("w_bcc").to_phase()) == "cubic"


@pytest.mark.parametrize(
    ("identifier", "free"),
    [("ni_fcc", 1), ("w_bcc", 1), ("ti_hcp", 2), ("mg_hcp", 2)],
)
def test_the_system_fixes_how_many_cell_parameters_may_vary(
    identifier: str, free: int
) -> None:
    result, _ = _determine(identifier, _scan(identifier))
    assert len(result.free_parameter_names) == free


def test_hexagonal_constraint_reproduces_the_textbook_quadratic_form() -> None:
    """A (h^2 + hk + k^2) + C l^2 must fall out of G*12 = G*11 / 2, not be written in."""

    phase = builtin_phase("ti_hcp").to_phase()
    reciprocal = phase.lattice.reciprocal_metric_tensor()
    assert reciprocal[0, 0] == pytest.approx(reciprocal[1, 1])
    assert reciprocal[0, 1] == pytest.approx(0.5 * reciprocal[0, 0])
    assert reciprocal[0, 2] == pytest.approx(0.0, abs=1e-12)
    assert reciprocal[1, 2] == pytest.approx(0.0, abs=1e-12)
    for indices in [(1, 0, 0), (1, 1, 0), (1, 0, 2), (2, 1, 3)]:
        h, k, l_index = indices
        vector = np.array(indices, dtype=float)
        quadratic = float(vector @ reciprocal @ vector)
        textbook = reciprocal[0, 0] * (h * h + h * k + k * k) + reciprocal[2, 2] * l_index**2
        assert quadratic == pytest.approx(textbook)


# ---------------------------------------------------------------------------
# The extrapolation functions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("function", ["nelson_riley", "bradley_jay", "cos_squared_over_sin"])
def test_every_extrapolation_function_vanishes_at_back_reflection(function: str) -> None:
    """Vanishing at theta = 90 degrees is the entire reason extrapolation works."""

    values = extrapolation_values(
        np.array([30.0, 90.0, 150.0, 179.9]),
        function=function,  # type: ignore[arg-type]
    )
    assert np.all(np.diff(values) < 0.0)
    assert values[-1] == pytest.approx(0.0, abs=1e-5)


def test_bradley_jay_drift_column_is_cohens_classical_sin_squared_two_theta() -> None:
    """sin^2(theta) cos^2(theta) = sin^2(2 theta) / 4, which is Cohen's column."""

    angles = np.linspace(20.0, 160.0, 25)
    theta = np.deg2rad(0.5 * angles)
    column = np.square(np.sin(theta)) * extrapolation_values(angles, function="bradley_jay")
    assert np.allclose(column, np.square(np.sin(2.0 * theta)) / 4.0)


def test_none_extrapolation_is_identically_zero() -> None:
    assert np.all(extrapolation_values(np.linspace(20.0, 160.0, 9), function="none") == 0.0)


def test_extrapolation_rejects_unknown_functions_and_divergent_angles() -> None:
    with pytest.raises(ValueError, match="function in"):
        extrapolation_values([45.0], function="taylor_sinclair")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="diverge"):
        extrapolation_values([0.0], function="nelson_riley")


# ---------------------------------------------------------------------------
# Accuracy on a clean pattern
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("identifier", ["ni_fcc", "w_bcc", "ti_hcp", "mg_hcp"])
def test_a_clean_pattern_recovers_the_generating_cell(identifier: str) -> None:
    phase = builtin_phase(identifier).to_phase()
    result, _ = _determine(identifier, _scan(identifier))
    assert result.a == pytest.approx(phase.lattice.a, rel=2.0e-5)
    assert result.c == pytest.approx(phase.lattice.c, rel=2.0e-5)
    assert result.relative_uncertainty < 1.0e-4


def test_uncertainties_are_real_numbers_not_placeholders() -> None:
    result, _ = _determine("ni_fcc", _scan("ni_fcc"))
    phase = builtin_phase("ni_fcc").to_phase()
    assert result.a_standard_uncertainty > 0.0
    # The determined value must sit within a few of its own stated sigmas of
    # the truth; an uncertainty that fails this is decoration, not a result.
    deviation = abs(result.a - phase.lattice.a) / result.a_standard_uncertainty
    assert deviation < 10.0


def test_more_counts_give_a_smaller_uncertainty() -> None:
    weak, _ = _determine("ni_fcc", _scan("ni_fcc", peak_counts=1_000.0, seed=11))
    strong, _ = _determine("ni_fcc", _scan("ni_fcc", peak_counts=1_000_000.0, seed=11))
    assert strong.a_standard_uncertainty < weak.a_standard_uncertainty


# ---------------------------------------------------------------------------
# The central claim: systematic error is not removed by averaging
# ---------------------------------------------------------------------------


def test_extrapolation_beats_averaging_against_an_injected_displacement() -> None:
    """The whole argument of the module, on data whose true answer is known."""

    identifier = "ni_fcc"
    phase = builtin_phase(identifier).to_phase()
    truth = phase.lattice.a
    measured = _scan(identifier, displacement_mm=0.10)

    def relative_error(**kwargs: object) -> float:
        result, _ = _determine(identifier, measured, **kwargs)
        return abs(result.a - truth) / truth

    naive = relative_error(method="average", extrapolation="none")
    uncorrected = relative_error(method="cohen", extrapolation="none")
    matched = relative_error(method="cohen", extrapolation="cos_squared_over_sin")
    nelson_riley = relative_error(method="cohen", extrapolation="nelson_riley")

    # Averaging cannot remove a theta-dependent error, so it is the worst.
    assert naive > 2.0e-4
    # A joint solution without a drift term is better but still systematically
    # wrong: it fits the same wrong positions, only more carefully.
    assert uncorrected < naive
    assert uncorrected > 1.0e-5
    # The function that matches the aberration's own angular form removes it
    # essentially completely.
    assert matched < 1.0e-6
    # Nelson-Riley approximates a combination of aberrations, so it is close
    # but not exact when only displacement is present.
    assert nelson_riley < 1.0e-5
    assert matched < nelson_riley


def test_restricting_to_high_angles_helps_but_does_not_replace_extrapolation() -> None:
    """cot(theta) shrinks towards back-reflection, so high angles are better."""

    identifier = "ni_fcc"
    truth = builtin_phase(identifier).to_phase().lattice.a
    measured = _scan(identifier, displacement_mm=0.10)
    everything, _ = _determine(identifier, measured, method="average", extrapolation="none")
    high_only, _ = _determine(
        identifier,
        measured,
        method="average",
        extrapolation="none",
        minimum_two_theta_deg=100.0,
    )
    matched, _ = _determine(
        identifier, measured, method="cohen", extrapolation="cos_squared_over_sin"
    )
    assert abs(high_only.a - truth) < abs(everything.a - truth)
    assert abs(matched.a - truth) < abs(high_only.a - truth)


def test_the_drift_coefficient_reports_the_shift_it_removed() -> None:
    identifier = "ni_fcc"
    measured = _scan(identifier, displacement_mm=0.10)
    result, _ = _determine(
        identifier, measured, method="cohen", extrapolation="cos_squared_over_sin"
    )
    shift = result.systematic_shift_deg
    injected = specimen_displacement_shift_deg(
        result.two_theta_deg, displacement_mm=0.10, goniometer_radius_mm=RADIUS_MM
    )
    # The removed shift must reproduce the injected aberration, up to the
    # constant offset that the cell itself absorbs.
    assert np.corrcoef(shift, injected)[0, 1] > 0.999
    assert float(np.max(np.abs(shift))) > 0.02
    assert abs(result.drift_coefficient) > 2.0 * result.drift_standard_uncertainty
    assert "significantly different from zero" in result.describe()


def test_a_well_aligned_specimen_gives_an_insignificant_drift_term() -> None:
    result, _ = _determine("ni_fcc", _scan("ni_fcc"), extrapolation="nelson_riley")
    assert abs(result.drift_coefficient) < 5.0 * result.drift_standard_uncertainty + 1e-6
    assert float(np.max(np.abs(result.systematic_shift_deg))) < 0.01


# ---------------------------------------------------------------------------
# The naive average, kept as a teaching comparison
# ---------------------------------------------------------------------------


def test_the_average_method_refuses_a_non_cubic_cell_and_says_why() -> None:
    with pytest.raises(ValueError, match="defined only for a cubic cell"):
        _determine("ti_hcp", _scan("ti_hcp"), method="average")


def test_the_average_of_a_clean_pattern_is_still_correct() -> None:
    """Averaging is not wrong; it is only defenceless against systematics."""

    phase = builtin_phase("ni_fcc").to_phase()
    result, _ = _determine("ni_fcc", _scan("ni_fcc"), method="average")
    assert result.a == pytest.approx(phase.lattice.a, rel=1.0e-5)
    assert result.extrapolation == "none"
    assert "cannot remove a systematic error" in result.describe()


def test_nelson_riley_plot_intercept_agrees_with_the_least_squares_answer() -> None:
    identifier = "ni_fcc"
    measured = _scan(identifier, displacement_mm=0.10)
    phase = builtin_phase(identifier).to_phase()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))
    indexing = index_peaks(table, phase, phase_name=identifier)
    plot = nelson_riley_extrapolation(indexing, phase)
    assert plot["extrapolation_function"].size == len(indexing)
    # The graphical intercept and the weighted joint solution are the same
    # method; they must agree far better than the uncorrected value does.
    joint, _ = _determine(identifier, measured, extrapolation="nelson_riley")
    assert float(plot["intercept"]) == pytest.approx(joint.a, rel=2.0e-5)
    # The per-reflection values must genuinely trend, otherwise the plot is
    # not showing what it claims to show.
    assert abs(float(plot["slope"])) > 1.0e-4


def test_nelson_riley_plot_refuses_a_non_cubic_phase() -> None:
    measured = _scan("ti_hcp")
    phase = builtin_phase("ti_hcp").to_phase()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))
    indexing = index_peaks(table, phase, phase_name="Ti")
    with pytest.raises(ValueError, match="only for a cubic cell"):
        nelson_riley_extrapolation(indexing, phase)


# ---------------------------------------------------------------------------
# Le Bail
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("identifier", ["ni_fcc", "ti_hcp"])
def test_le_bail_recovers_the_generating_cell_on_a_clean_pattern(identifier: str) -> None:
    phase = builtin_phase(identifier).to_phase()
    result = determine_lattice_parameters_le_bail(
        _scan(identifier), phase, systematic="none", cycles=10
    )
    assert result.a == pytest.approx(phase.lattice.a, rel=5.0e-6)
    assert result.c == pytest.approx(phase.lattice.c, rel=5.0e-6)
    # A correct whole-pattern model on Poisson data sits near unity.
    assert 0.5 < result.reduced_chi_squared < 3.0


@pytest.mark.parametrize("identifier", ["ni_fcc", "ti_hcp"])
def test_le_bail_refines_the_injected_displacement_itself(identifier: str) -> None:
    """Not merely 'a is right': the aberration comes back in millimetres."""

    phase = builtin_phase(identifier).to_phase()
    measured = _scan(identifier, displacement_mm=0.10)
    result = determine_lattice_parameters_le_bail(
        measured,
        phase,
        systematic="displacement",
        goniometer_radius_mm=RADIUS_MM,
        cycles=10,
    )
    assert result.drift_coefficient == pytest.approx(0.10, rel=0.02)
    assert result.a == pytest.approx(phase.lattice.a, rel=2.0e-6)
    assert result.c == pytest.approx(phase.lattice.c, rel=2.0e-6)
    assert 0.5 < result.reduced_chi_squared < 3.0


def test_le_bail_chi_squared_exposes_an_unmodelled_aberration() -> None:
    """The goodness of fit must fail loudly when the model is wrong."""

    phase = builtin_phase("ni_fcc").to_phase()
    measured = _scan("ni_fcc", displacement_mm=0.10)
    ignored = determine_lattice_parameters_le_bail(
        measured, phase, systematic="none", cycles=10
    )
    modelled = determine_lattice_parameters_le_bail(
        measured, phase, systematic="displacement", cycles=10
    )
    assert ignored.reduced_chi_squared > 5.0 * modelled.reduced_chi_squared
    assert abs(ignored.a - phase.lattice.a) > 20.0 * abs(modelled.a - phase.lattice.a)


def test_le_bail_and_cohen_agree_on_the_same_pattern() -> None:
    """Two methods with different failure modes must give the same answer."""

    identifier = "ti_hcp"
    phase = builtin_phase(identifier).to_phase()
    measured = _scan(identifier, displacement_mm=0.10)
    whole_pattern = determine_lattice_parameters_le_bail(
        measured, phase, systematic="displacement", cycles=10
    )
    positions, _ = _determine(
        identifier, measured, method="cohen", extrapolation="cos_squared_over_sin"
    )
    assert whole_pattern.a == pytest.approx(positions.a, rel=5.0e-5)
    assert whole_pattern.c == pytest.approx(positions.c, rel=5.0e-5)


def test_le_bail_validates_its_inputs() -> None:
    phase = builtin_phase("ni_fcc").to_phase()
    measured = _scan("ni_fcc")
    with pytest.raises(ValueError, match="systematic in"):
        determine_lattice_parameters_le_bail(measured, phase, systematic="both")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="at least one cycle"):
        determine_lattice_parameters_le_bail(measured, phase, cycles=0)
    bare = MeasuredPowderPattern(
        name="no radiation",
        two_theta_deg=measured.two_theta_deg,
        intensity=measured.intensity,
        synthetic=True,
    )
    with pytest.raises(ValueError, match="needs a radiation"):
        determine_lattice_parameters_le_bail(bare, phase)


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------


def test_strain_is_reported_against_the_reference_cell_and_labelled_a_strain() -> None:
    identifier = "ni_fcc"
    phase = builtin_phase(identifier).to_phase()
    result, _ = _determine(identifier, _scan(identifier))
    strain = result.strain_relative_to_reference
    assert strain is not None
    assert strain == pytest.approx((result.a - phase.lattice.a) / phase.lattice.a)
    prose = result.describe()
    assert "lattice strain along a" in prose
    # Refusing to call a strain a stress is a scientific requirement, not a
    # stylistic one.
    assert "X-ray elastic constants" in prose
    assert "several specimen tilts" in prose


def test_the_determined_cell_round_trips_into_a_lattice() -> None:
    result, _ = _determine("ti_hcp", _scan("ti_hcp"))
    lattice = result.to_lattice()
    assert lattice.a == pytest.approx(result.a)
    assert lattice.c == pytest.approx(result.c)
    assert lattice.gamma_deg == pytest.approx(120.0, abs=1e-6)
    assert lattice.crystal_frame is builtin_phase("ti_hcp").to_phase().lattice.crystal_frame


def test_json_contract_carries_the_schema_and_the_per_reflection_detail() -> None:
    result, _ = _determine("ni_fcc", _scan("ni_fcc", displacement_mm=0.10))
    payload = result.to_json()
    assert payload["schema"] == LATTICE_PARAMETER_SCHEMA
    assert payload["cell"]["a"] == pytest.approx(result.a)
    assert payload["standard_uncertainty"]["a"] == pytest.approx(result.a_standard_uncertainty)
    assert len(payload["reflections"]) == result.reflection_count
    assert payload["reflections"][0]["systematic_shift_deg"] != 0.0
    assert payload["extrapolation"] in EXTRAPOLATION_FUNCTIONS


def test_describe_names_the_method_the_system_and_the_precision_reached() -> None:
    result, _ = _determine("ti_hcp", _scan("ti_hcp"))
    prose = result.describe()
    assert "reciprocal metric tensor" in prose
    assert "hexagonal system leaves 2 free" in prose
    assert "c/a" in prose
    assert "Cullity" in prose
    assert "Cohen" in prose


def test_a_result_without_a_drift_term_says_the_error_is_still_inside_the_cell() -> None:
    result, _ = _determine("ni_fcc", _scan("ni_fcc"), extrapolation="none")
    assert "inside the quoted cell rather than removed from it" in result.describe()


def test_determination_validates_its_inputs() -> None:
    measured = _scan("ni_fcc")
    phase = builtin_phase("ni_fcc").to_phase()
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))
    indexing = index_peaks(table, phase, phase_name="Ni")
    with pytest.raises(ValueError, match="determine_lattice_parameters_le_bail"):
        determine_lattice_parameters(indexing, phase, method="le_bail")
    with pytest.raises(ValueError, match="extrapolation in"):
        determine_lattice_parameters(indexing, phase, extrapolation="taylor")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="survives the angular restriction"):
        determine_lattice_parameters(indexing, phase, minimum_two_theta_deg=179.0)


# ---------------------------------------------------------------------------
# Lower symmetry: three free parameters, and crowded patterns
# ---------------------------------------------------------------------------


def test_an_orthorhombic_cell_determines_three_independent_edges() -> None:
    """Nothing above exercises a system with more than two free parameters.

    Alpha-uranium is orthorhombic, so `a`, `b` and `c` are independent and the
    constraint matrix has three columns. A determination that silently tied any
    two of them together would still look plausible on a cubic or hexagonal
    test.
    """

    phase = builtin_phase("alpha_u").to_phase()
    assert crystal_system_of(phase) == "orthorhombic"
    result, indexing = _determine("alpha_u", _scan("alpha_u", peak_counts=40_000.0, seed=3))
    assert len(result.free_parameter_names) == 3
    assert indexing is not None
    assert indexing.indexed_count > 20
    for determined, truth in (
        (result.a, phase.lattice.a),
        (result.b, phase.lattice.b),
        (result.c, phase.lattice.c),
    ):
        assert determined == pytest.approx(truth, rel=5.0e-5)
    # The three edges are genuinely different, so a tied parameterization could
    # not have produced this.
    assert len({round(value, 3) for value in (result.a, result.b, result.c)}) == 3


def test_le_bail_beats_peak_fitting_on_a_crowded_low_symmetry_pattern() -> None:
    """The reason the whole-pattern method exists, measured.

    On a heavily overlapped orthorhombic pattern the single-peak route runs out
    of resolvable lines and its goodness of fit degrades badly, while the
    whole-pattern fit stays near unity. Both must still land on the same cell:
    two methods with different failure modes agreeing is the strongest evidence
    available without an external standard.
    """

    phase = builtin_phase("alpha_u").to_phase()
    measured = _scan("alpha_u", peak_counts=40_000.0, seed=3)
    positions, _ = _determine("alpha_u", measured)
    whole_pattern = determine_lattice_parameters_le_bail(
        measured, phase, systematic="none", cycles=8
    )
    assert whole_pattern.reduced_chi_squared < 0.2 * positions.reduced_chi_squared
    assert whole_pattern.reduced_chi_squared < 3.0
    for from_profile, from_positions in (
        (whole_pattern.a, positions.a),
        (whole_pattern.b, positions.b),
        (whole_pattern.c, positions.c),
    ):
        assert from_profile == pytest.approx(from_positions, rel=1.0e-4)


def test_a_trigonal_cell_uses_the_hexagonal_constraint() -> None:
    phase = builtin_phase("quartz_alpha").to_phase()
    assert crystal_system_of(phase) == "trigonal"
    result = determine_lattice_parameters_le_bail(
        _scan("quartz_alpha", peak_counts=40_000.0, seed=3),
        phase,
        systematic="none",
        cycles=8,
    )
    assert len(result.free_parameter_names) == 2
    assert result.a == pytest.approx(phase.lattice.a, rel=2.0e-5)
    assert result.c == pytest.approx(phase.lattice.c, rel=2.0e-5)
    assert result.gamma_deg == pytest.approx(120.0, abs=1.0e-6)


def test_the_le_bail_cell_bounds_survive_a_negative_metric_component() -> None:
    """A relative bound built by multiplication inverts on a negative parameter.

    Off-diagonal components of G* are negative whenever the corresponding
    reciprocal angle is obtuse, which a triclinic cell routinely is. Building
    the bounds as `start * 0.9` and `start * 1.1` would then put the lower bound
    above the upper one and the optimizer would refuse the problem outright.
    The bounds are built by magnitude instead; this checks the arithmetic
    directly, since no built-in phase is triclinic.
    """

    start = np.array([0.12, -0.03, 0.0])
    span = np.maximum(0.1 * np.abs(start), 1.0e-4)
    lower = start - span
    upper = start + span
    assert np.all(lower < upper)
    assert np.all(lower <= start)
    assert np.all(start <= upper)
    # A component starting at exactly zero still gets a usable window.
    assert upper[2] - lower[2] == pytest.approx(2.0e-4)


def test_le_bail_describes_the_aberration_it_refined_in_physical_units() -> None:
    """A whole-pattern fit refines the aberration itself, not a drift coefficient.

    Its `extrapolation` field is "none" because no extrapolation function was
    used -- which is not the same thing as no correction having been made.
    Reporting the first as the second would tell a reader the displacement was
    left in the cell when it had just been taken out of it.
    """

    phase = builtin_phase("ni_fcc").to_phase()
    measured = _scan("ni_fcc", displacement_mm=0.10)

    refined = determine_lattice_parameters_le_bail(
        measured, phase, systematic="displacement", cycles=10
    )
    prose = refined.describe()
    assert "specimen displacement of" in prose
    assert "mm was" in prose
    assert "rather than absorbed into the cell" in prose
    assert "No systematic-error term was refined" not in prose

    zero = determine_lattice_parameters_le_bail(
        measured, phase, systematic="zero", cycles=6
    )
    assert "detector zero of" in zero.describe()
    assert "degrees 2*theta" in zero.describe()

    ignored = determine_lattice_parameters_le_bail(
        measured, phase, systematic="none", cycles=6
    )
    assert "No instrumental aberration was refined" in ignored.describe()
    assert "inside the quoted cell" in ignored.describe()
