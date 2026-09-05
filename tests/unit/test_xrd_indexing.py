"""Indexing measured peaks against a known phase, and scoring the assignment.

The correctness assertions here are structural rather than numerical wherever
possible: a cubic F-lattice may only produce unmixed-parity indices, ``Q`` must
be proportional to ``h^2 + k^2 + l^2``, and a global assignment must be
one-to-one. Those hold whatever the fitted angles happen to be.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_indexing import (
    INDEXED_REFLECTION_SCHEMA,
    PEAK_INDEXING_SCHEMA,
    PeakIndexing,
    index_peaks,
)
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import PeakTable, detect_and_fit_peaks

FWHM = 0.12


def _scan(
    identifier: str,
    *,
    two_theta_range_deg: tuple[float, float] = (25.0, 150.0),
    peak_counts: float = 30000.0,
    seed: int = 5,
) -> tuple[MeasuredPowderPattern, tuple[float, ...]]:
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
    measured = MeasuredPowderPattern(
        name=f"{identifier} synthetic",
        two_theta_deg=axis,
        intensity=np.random.default_rng(seed).poisson(expected).astype(float),
        radiation=radiation,
        synthetic=True,
    )
    return measured, tuple(float(item.two_theta_deg) for item in pattern.reflections)


def _indexed(identifier: str, **kwargs: object) -> PeakIndexing:
    measured, _ = _scan(identifier)
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))
    return index_peaks(
        table,
        builtin_phase(identifier).to_phase(),
        phase_name=identifier,
        **kwargs,  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# The assignment itself
# ---------------------------------------------------------------------------


def test_a_clean_cubic_pattern_indexes_completely() -> None:
    indexing = _indexed("ni_fcc")
    assert indexing.indexed_count == 7
    assert indexing.unindexed_peaks == ()
    assert indexing.indexed_fraction == pytest.approx(1.0)
    assert indexing.mean_absolute_delta_two_theta_deg < 1.0e-3


def test_face_centred_cubic_indices_are_unmixed_in_parity() -> None:
    """The centring is a structural fact, so it is checked structurally."""

    indexing = _indexed("ni_fcc")
    for reflection in indexing:
        parities = {abs(value) % 2 for value in reflection.miller_indices}
        assert len(parities) == 1, f"{reflection.miller_indices} mixes odd and even"


def test_body_centred_cubic_indices_have_an_even_index_sum() -> None:
    indexing = _indexed("w_bcc")
    assert indexing.indexed_count >= 5
    for reflection in indexing:
        assert sum(reflection.miller_indices) % 2 == 0


def test_cubic_q_is_proportional_to_the_index_sum_of_squares() -> None:
    """1/d^2 = (h^2 + k^2 + l^2) / a^2 is the whole of cubic indexing."""

    indexing = _indexed("ni_fcc")
    sums = np.array(
        [sum(value**2 for value in item.miller_indices) for item in indexing]
    )
    q_values = np.array([1.0 / item.d_observed_angstrom**2 for item in indexing])
    ratios = q_values / sums
    # Every reflection must give the same a, to a few parts in 1e5.
    assert float(np.std(ratios) / np.mean(ratios)) < 1.0e-4
    lattice_parameter = float(np.mean(1.0 / np.sqrt(ratios)))
    assert lattice_parameter == pytest.approx(
        builtin_phase("ni_fcc").to_phase().lattice.a, rel=1.0e-4
    )


def test_the_assignment_is_one_to_one() -> None:
    indexing = _indexed("ti_hcp")
    indices = [item.miller_indices for item in indexing]
    assert len(indices) == len(set(indices))
    positions = [item.peak.two_theta_deg for item in indexing]
    assert len(positions) == len(set(positions))


def test_global_assignment_beats_greedy_on_a_crowded_pattern() -> None:
    """A greedy pass would give two peaks the same index; this must not."""

    indexing = _indexed("mg_hcp")
    indices = [item.miller_indices for item in indexing]
    assert len(indices) == len(set(indices))
    assert indexing.indexed_count >= 24


def test_indexing_is_monotonic_in_angle_and_spacing() -> None:
    indexing = _indexed("ni_fcc")
    angles = [item.peak.two_theta_deg for item in indexing]
    spacings = [item.d_observed_angstrom for item in indexing]
    assert angles == sorted(angles)
    assert spacings == sorted(spacings, reverse=True)


def test_overlapped_hexagonal_reflections_are_reported_not_hidden() -> None:
    """Magnesium 015 and 122 sit 0.23 degrees apart; single-peak fitting strains.

    The point of the test is not that the strain is absent, but that the report
    names it: a badly fitted residue appears as an unindexed peak and the line
    it failed to resolve appears as unobserved. Recovering both is what the
    whole-pattern method exists for.
    """

    indexing = _indexed("mg_hcp")
    assert indexing.unindexed_peaks or indexing.unobserved_indices
    for peak in indexing.unindexed_peaks:
        # Anything left over here is a failed fit, not a good peak thrown away.
        assert peak.reduced_chi_squared > 5.0
    prose = indexing.describe()
    if indexing.unindexed_peaks:
        assert "not indexed" in prose
    if indexing.unobserved_indices:
        assert "not observed" in prose


# ---------------------------------------------------------------------------
# Residuals carry the diagnosis
# ---------------------------------------------------------------------------


def test_a_wrong_zero_shows_as_same_signed_residuals_and_is_named() -> None:
    measured, _ = _scan("ni_fcc")
    shifted = MeasuredPowderPattern(
        name="zero-shifted",
        two_theta_deg=np.asarray(measured.two_theta_deg) + 0.05,
        intensity=measured.intensity,
        radiation=measured.radiation,
        synthetic=True,
    )
    table = detect_and_fit_peaks(shifted, instrument=InstrumentBroadening.ideal(FWHM))
    indexing = index_peaks(table, builtin_phase("ni_fcc").to_phase(), phase_name="Ni")
    residuals = indexing.delta_two_theta_deg
    assert np.all(residuals > 0.0)
    assert np.allclose(residuals, 0.05, atol=2.0e-3)
    assert "same sign" in indexing.describe()
    assert "zero-point or specimen-displacement" in indexing.describe()


def test_normalized_residual_flags_what_counting_statistics_cannot_explain() -> None:
    indexing = _indexed("ni_fcc")
    # With the generating cell, every residual is within a few sigma.
    assert max(abs(item.normalized_residual) for item in indexing) < 25.0
    measured, _ = _scan("ni_fcc")
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))
    # Against a deliberately wrong cell, they are not.
    wrong = index_peaks(
        table, builtin_phase("cu_fcc").to_phase(), phase_name="Cu", tolerance_deg=1.0
    )
    assert max(abs(item.normalized_residual) for item in wrong) > 100.0


def test_delta_q_is_the_discrepancy_in_reciprocal_square_spacing() -> None:
    indexing = _indexed("ni_fcc")
    for item in indexing:
        expected = 1.0 / item.d_observed_angstrom**2 - 1.0 / item.d_calculated_angstrom**2
        assert item.delta_q == pytest.approx(expected)


# ---------------------------------------------------------------------------
# Figures of merit
# ---------------------------------------------------------------------------


def test_figures_of_merit_reward_a_correct_cell_over_a_wrong_one() -> None:
    measured, _ = _scan("ni_fcc")
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))
    right = index_peaks(table, builtin_phase("ni_fcc").to_phase(), phase_name="Ni")
    wrong = index_peaks(
        table, builtin_phase("cu_fcc").to_phase(), phase_name="Cu", tolerance_deg=1.5
    )
    assert right.figure_of_merit_m()[0] > 100.0 * wrong.figure_of_merit_m()[0]
    assert right.figure_of_merit_f()[0] > 100.0 * wrong.figure_of_merit_f()[0]


def test_figures_of_merit_report_the_n_they_were_computed_over() -> None:
    """M_7 and M_20 are not comparable, so the N must travel with the value."""

    indexing = _indexed("ni_fcc")
    value, count = indexing.figure_of_merit_m(count=20)
    assert count == 7  # only seven lines exist to index
    assert np.isfinite(value)
    _, capped = indexing.figure_of_merit_f(count=30)
    assert capped == 7
    _, requested = indexing.figure_of_merit_m(count=3)
    assert requested == 3


def test_figure_of_merit_penalises_a_cell_that_predicts_too_many_lines() -> None:
    """N_poss in the denominator is the part that does the work."""

    indexing = _indexed("ni_fcc")
    baseline, _ = indexing.figure_of_merit_m()
    inflated = PeakIndexing(
        name=indexing.name,
        phase_name=indexing.phase_name,
        reflections=indexing.reflections,
        tolerance_deg=indexing.tolerance_deg,
        radiation=indexing.radiation,
        settings={**dict(indexing.settings), "possible_lines": 4.0 * float(
            indexing.settings["possible_lines"]
        )},
    )
    quadrupled, _ = inflated.figure_of_merit_m()
    assert quadrupled == pytest.approx(baseline / 4.0)


def test_empty_indexing_describes_itself_without_pretending() -> None:
    measured, _ = _scan("ni_fcc")
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))
    impossible = index_peaks(
        table,
        builtin_phase("nacl").to_phase(),
        phase_name="NaCl",
        tolerance_deg=1.0e-4,
    )
    assert impossible.indexed_count == 0
    assert np.isnan(impossible.figure_of_merit_m()[0])
    assert np.isnan(impossible.mean_absolute_delta_two_theta_deg)
    prose = impossible.describe()
    assert "assigned no reflection" in prose
    assert "tolerance is" in prose


# ---------------------------------------------------------------------------
# Contracts and validation
# ---------------------------------------------------------------------------


def test_json_contract_carries_the_schemas_and_the_residuals() -> None:
    indexing = _indexed("ni_fcc")
    payload = indexing.to_json()
    assert payload["schema"] == PEAK_INDEXING_SCHEMA
    assert payload["indexed_count"] == indexing.indexed_count
    assert len(payload["reflections"]) == len(indexing)
    first = payload["reflections"][0]
    assert first["schema"] == INDEXED_REFLECTION_SCHEMA
    assert first["miller_indices"] == list(indexing.reflections[0].miller_indices)
    assert first["delta_two_theta_deg"] == pytest.approx(
        indexing.reflections[0].delta_two_theta_deg
    )
    assert payload["figure_of_merit_m"]["count"] == 7


def test_describe_states_the_method_and_cites_the_figures_of_merit() -> None:
    prose = _indexed("ni_fcc").describe()
    assert "global one-to-one assignment" in prose
    assert "de Wolff" in prose
    assert "Smith & Snyder" in prose or "Smith and Snyder" in prose
    assert "M_7" in prose


def test_indexing_validates_its_inputs() -> None:
    measured, _ = _scan("ni_fcc")
    table = detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))
    phase = builtin_phase("ni_fcc").to_phase()
    with pytest.raises(ValueError, match="empty peak table"):
        index_peaks(PeakTable(name="empty", peaks=()), phase)
    with pytest.raises(ValueError, match="positive tolerance_deg"):
        index_peaks(table, phase, tolerance_deg=0.0)
    bare = PeakTable(name="no radiation", peaks=table.peaks)
    with pytest.raises(ValueError, match="needs a radiation"):
        index_peaks(bare, phase)
    with pytest.raises(ValueError, match="below minimum_relative_intensity"):
        index_peaks(table, phase, minimum_relative_intensity=2.0)


def test_indexed_reflection_validates_its_invariants() -> None:
    indexing = _indexed("ni_fcc")
    from dataclasses import replace

    reflection = indexing.reflections[0]
    with pytest.raises(ValueError, match="three integers"):
        replace(reflection, miller_indices=(1, 1))
    with pytest.raises(ValueError, match="multiplicity"):
        replace(reflection, multiplicity=0)
    with pytest.raises(ValueError, match="spacings"):
        replace(reflection, d_observed_angstrom=0.0)
    with pytest.raises(ValueError, match="name must be non-empty"):
        PeakIndexing(name="  ", phase_name="Ni", reflections=())
