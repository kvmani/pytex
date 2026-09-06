"""Ranking several candidate phases against one measured powder pattern.

The assertions here are chosen so that they follow from crystallography rather
than from the particular numbers a fit happens to produce. A pattern generated
from a phase must rank that phase first. A face-centred candidate offered a
body-centred pattern must fail on *completeness*, because the two differ in
which lines are absent, and that is a structural fact about the centrings, not
a property of this scan. Two cells that differ by less than the matching
tolerance must not be reported as distinguished, because they are not.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from pytex.app.phases import builtin_phase
from pytex.core.lattice import Phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern
from pytex.diffraction.xrd_instrument import InstrumentBroadening
from pytex.diffraction.xrd_measurement import MeasuredPowderPattern
from pytex.diffraction.xrd_peaks import PeakTable, detect_and_fit_peaks
from pytex.diffraction.xrd_phase_identification import (
    CRITERION_NAMES,
    DEFAULT_CRITERION_WEIGHTS,
    PHASE_CANDIDATE_SCORE_SCHEMA,
    PHASE_IDENTIFICATION_SCHEMA,
    PhaseCandidateScore,
    PhaseIdentification,
    identify_phase,
    identify_phase_from_pattern,
)

FWHM = 0.12
RANGE = (25.0, 140.0)


def _scan(
    identifier: str,
    *,
    peak_counts: float = 30000.0,
    seed: int = 11,
) -> MeasuredPowderPattern:
    """Return a Poisson-noised synthetic diffractogram of a built-in phase."""

    phase = builtin_phase(identifier).to_phase()
    radiation = RadiationSpec.cu_ka()
    pattern = generate_xrd_pattern(
        phase,
        radiation=radiation,
        two_theta_range_deg=RANGE,
        resolution_deg=0.01,
        broadening_fwhm_deg=FWHM,
        profile="pseudo_voigt",
        max_index=6,
    )
    profile = np.asarray(pattern.intensity_grid, dtype=float)
    expected = profile / profile.max() * peak_counts + 150.0
    return MeasuredPowderPattern(
        name=f"{identifier} synthetic",
        two_theta_deg=np.asarray(pattern.two_theta_grid_deg, dtype=float),
        intensity=np.random.default_rng(seed).poisson(expected).astype(float),
        radiation=radiation,
        synthetic=True,
    )


def _table(identifier: str, **kwargs: object) -> PeakTable:
    """Return the fitted peaks of a synthetic scan of a built-in phase."""

    measured = _scan(identifier, **kwargs)  # type: ignore[arg-type]
    return detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))


def _dilated(phase: Phase, scale: float, *, name: str) -> Phase:
    """Return the same structure on a uniformly dilated cell.

    Used to build a candidate that is genuinely indistinguishable from another
    on a laboratory scan, which is the case an identification has to be able to
    decline rather than resolve.
    """

    lattice = replace(
        phase.lattice,
        a=phase.lattice.a * scale,
        b=phase.lattice.b * scale,
        c=phase.lattice.c * scale,
    )
    return replace(
        phase,
        lattice=lattice,
        unit_cell=replace(phase.unit_cell, lattice=lattice),
        name=name,
    )


def _candidates(*identifiers: str) -> dict[str, Phase]:
    """Return built-in phases keyed by their catalogue names."""

    return {
        builtin_phase(identifier).name: builtin_phase(identifier).to_phase()
        for identifier in identifiers
    }


# ---------------------------------------------------------------------------
# The ranking
# ---------------------------------------------------------------------------


def test_the_generating_phase_is_ranked_first() -> None:
    report = identify_phase(
        _table("ni_fcc"), _candidates("cu_fcc", "fe_bcc", "ni_fcc", "al_fcc")
    )
    assert report.best.phase_name == "Nickel (fcc)"
    assert report.is_conclusive
    assert report.is_decisive


def test_the_ranking_is_independent_of_the_order_candidates_were_offered() -> None:
    """A ranking that depended on input order would not be a measurement."""

    table = _table("ni_fcc")
    forward = identify_phase(table, _candidates("ni_fcc", "cu_fcc", "al_fcc"))
    reverse = identify_phase(table, _candidates("al_fcc", "cu_fcc", "ni_fcc"))
    assert [item.phase_name for item in forward] == [item.phase_name for item in reverse]
    assert forward.best.score == pytest.approx(reverse.best.score)


def test_every_criterion_is_bounded_to_the_unit_interval() -> None:
    report = identify_phase(
        _table("ti_hcp"), _candidates("ti_hcp", "zr_hcp", "mg_hcp", "fe_bcc")
    )
    for candidate in report:
        for name, value in candidate.criteria.items():
            assert name in CRITERION_NAMES
            assert np.isnan(value) or 0.0 <= value <= 1.0, f"{name} = {value}"
        assert 0.0 <= candidate.score <= 1.0


def test_the_correct_phase_scores_near_one_on_its_own_pattern() -> None:
    """A noise-free-in-position synthetic pattern leaves no room for excuses."""

    report = identify_phase(_table("ni_fcc"), _candidates("ni_fcc"))
    best = report.best
    assert best.score > 0.95
    assert best.explained_intensity_fraction == pytest.approx(1.0)
    assert best.completeness == pytest.approx(1.0)
    assert best.strongest_unexplained_fraction == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# What separates the candidates, criterion by criterion
# ---------------------------------------------------------------------------


def test_a_centring_is_caught_by_completeness_not_by_position() -> None:
    """The fcc/bcc distinction is about absent lines, so completeness must carry it.

    A body-centred candidate offered a face-centred pattern predicts a line at
    ``h+k+l`` even where the face-centred structure is extinct. Some of its
    lines will coincide with observed peaks by arithmetic accident, so a
    position-only score is a poor discriminator; the reflections it insists on
    and the scan does not show are the evidence.
    """

    report = identify_phase(_table("ni_fcc"), _candidates("ni_fcc", "fe_bcc"))
    correct = next(item for item in report if item.phase_name == "Nickel (fcc)")
    wrong = next(item for item in report if item.phase_name == "Ferrite (bcc Fe)")
    assert correct.completeness > wrong.completeness
    assert wrong.strongest_unobserved_relative_intensity > 0.5


def test_an_unexplained_strong_peak_is_reported_as_such() -> None:
    """A candidate whose cell is far wrong leaves the strongest peak unindexed."""

    report = identify_phase(_table("ni_fcc"), _candidates("ni_fcc", "nacl"))
    wrong = next(item for item in report if item.phase_name == "Halite (NaCl)")
    assert wrong.strongest_unexplained_fraction > 0.2
    assert wrong.explained_intensity_fraction < 0.9


def test_a_candidate_predicting_nothing_in_range_is_recorded_not_raised() -> None:
    """One impossible candidate must not abort the comparison of the others.

    A phase whose lines all fall outside the measured range cannot be indexed.
    That is a decisive result about *that* candidate, so it is scored zero with
    a stated reason, and the remaining candidates stay rankable -- which is the
    entire purpose of offering several at once.
    """

    nickel = builtin_phase("ni_fcc").to_phase()
    # A cell this small puts every d spacing below lambda/2, so Bragg's law has
    # no solution at any angle and the candidate predicts nothing at all.
    impossible = _dilated(nickel, 0.30, name="Unphysically small cell")
    report = identify_phase(
        _table("ni_fcc"), {"Nickel (fcc)": nickel, impossible.name: impossible}
    )
    rejected = next(
        item for item in report if item.phase_name == "Unphysically small cell"
    )
    assert rejected.indexing is None
    assert rejected.score == 0.0
    assert "predicts no reflection" in rejected.rejection
    assert "could not be indexed" in rejected.describe()
    assert report.best.phase_name == "Nickel (fcc)"


def test_a_candidate_matching_no_peak_scores_zero_without_being_rejected() -> None:
    """Predicting lines that miss every peak is a different failure from predicting none.

    Copper's cell is 2.6 per cent larger than nickel's, so on a nickel pattern
    it puts lines everywhere *except* where the peaks are. Indexing succeeds and
    assigns nothing, which is a measurement about copper rather than an error,
    and it must score zero on that evidence rather than on a missing result.

    The cell-scale search is switched off here, because with it on the question
    being asked changes: a candidate free to dilate is no longer the same
    candidate, and copper stretched onto nickel's lines is a different finding
    from copper missing them.
    """

    report = identify_phase(
        _table("ni_fcc"),
        _candidates("ni_fcc", "cu_fcc"),
        tolerance_deg=0.05,
        cell_scale_range=0.0,
    )
    copper = next(item for item in report if item.phase_name == "Copper (fcc)")
    assert copper.indexing is not None
    assert copper.indexed_count == 0
    assert copper.score == 0.0
    assert np.isnan(copper.position_score)
    assert report.best.phase_name == "Nickel (fcc)"


def test_the_position_score_measures_the_discrepancy_against_the_tolerance() -> None:
    """The documented definition, checked against the reported residual."""

    candidate = identify_phase(
        _table("ni_fcc"), _candidates("ni_fcc"), tolerance_deg=0.2
    ).best
    assert candidate.indexing is not None
    expected = 1.0 - candidate.indexing.mean_absolute_delta_two_theta_deg / 0.2
    assert candidate.position_score == pytest.approx(expected)


def test_widening_the_tolerance_does_not_promote_a_wrong_candidate() -> None:
    """A user who widens the window until their preferred phase indexes gains nothing.

    A laxer window admits more matches for every candidate, including sloppy
    ones, and each discrepancy is then judged against the laxer standard it was
    admitted under. The ranking must therefore survive a tolerance far wider
    than any real instrument error.
    """

    table = _table("ni_fcc")
    candidates = _candidates("ni_fcc", "fe_bcc", "al_fcc", "nacl")
    for tolerance in (0.1, 0.3, 1.0):
        report = identify_phase(table, candidates, tolerance_deg=tolerance)
        assert report.best.phase_name == "Nickel (fcc)", f"at {tolerance} degrees"
        assert report.is_conclusive and report.is_decisive


# ---------------------------------------------------------------------------
# The two qualifications on the winner
# ---------------------------------------------------------------------------


def test_a_pattern_none_of_the_candidates_explains_is_called_inconclusive() -> None:
    """A ranking always has a winner; an identification must be able to decline."""

    report = identify_phase(_table("ni_fcc"), _candidates("nacl", "quartz_alpha"))
    assert not report.is_conclusive
    assert "none of the candidates offered" in report.describe()


def test_two_indistinguishable_candidates_are_reported_as_not_decisive() -> None:
    """Cells closer than the matching tolerance are not separated, and it says so."""

    nickel = builtin_phase("ni_fcc").to_phase()
    twin = _dilated(nickel, 1.0002, name="Nickel, dilated")
    report = identify_phase(
        _table("ni_fcc"),
        {"Nickel (fcc)": nickel, "Nickel, dilated": twin},
        decisive_margin=0.05,
    )
    assert report.is_conclusive
    assert not report.is_decisive
    assert "not distinguished by this scan" in report.describe()


def test_a_single_candidate_has_no_margin_and_is_decisive_by_convention() -> None:
    report = identify_phase(_table("ni_fcc"), _candidates("ni_fcc"))
    assert report.runner_up is None
    assert np.isnan(report.margin)
    assert report.is_decisive
    assert "a check on one phase" in report.describe()


# ---------------------------------------------------------------------------
# Scoring mechanics
# ---------------------------------------------------------------------------


def test_an_undefined_criterion_is_omitted_rather_than_scored_zero() -> None:
    """"Not measurable here" and "measured and bad" are different findings.

    Intensity agreement is undefined on a single indexed line, since a
    similarity between two one-element distributions is identically one. A
    candidate perfect on the three criteria that *are* defined must therefore
    score one, not 0.85 -- the missing weight is renormalized away rather than
    counted against it.
    """

    indexing = identify_phase(_table("ni_fcc"), _candidates("ni_fcc")).best.indexing
    assert indexing is not None
    perfect = PhaseCandidateScore(
        phase_name="stub",
        indexing=indexing,
        explained_intensity_fraction=1.0,
        completeness=1.0,
        position_score=1.0,
        intensity_agreement=float("nan"),
    )
    assert perfect.score == pytest.approx(1.0)

    # A candidate that could not be indexed at all scores zero outright, which
    # is a different statement from every criterion being undefined.
    rejected = PhaseCandidateScore(
        phase_name="stub",
        indexing=None,
        rejection="predicts nothing in range",
    )
    assert rejected.score == 0.0


def test_criterion_weights_are_validated() -> None:
    table = _table("ni_fcc")
    with pytest.raises(ValueError, match="Unknown scoring criterion"):
        identify_phase(table, _candidates("ni_fcc"), weights={"nonsense": 1.0})
    with pytest.raises(ValueError, match="finite and non-negative"):
        identify_phase(table, _candidates("ni_fcc"), weights={"completeness": -1.0})
    with pytest.raises(ValueError, match="strictly positive"):
        identify_phase(table, _candidates("ni_fcc"), weights={"completeness": 0.0})


def test_the_default_weights_favour_positions_over_intensities() -> None:
    """Preferred orientation moves intensities without moving positions."""

    assert (
        DEFAULT_CRITERION_WEIGHTS["intensity_agreement"]
        < DEFAULT_CRITERION_WEIGHTS["explained_intensity_fraction"]
    )
    assert sum(DEFAULT_CRITERION_WEIGHTS.values()) == pytest.approx(1.0)


def test_reweighting_changes_the_score_it_is_supposed_to_change() -> None:
    table = _table("ni_fcc")
    candidates = _candidates("ni_fcc", "fe_bcc")
    balanced = identify_phase(table, candidates)
    positions_only = identify_phase(
        table, candidates, weights={"position_score": 1.0}
    )
    wrong_balanced = next(
        item for item in balanced if item.phase_name == "Ferrite (bcc Fe)"
    )
    wrong_positions = next(
        item for item in positions_only if item.phase_name == "Ferrite (bcc Fe)"
    )
    assert wrong_positions.score != pytest.approx(wrong_balanced.score)
    assert wrong_positions.score == pytest.approx(wrong_positions.position_score)


def test_duplicate_candidate_names_are_numbered_rather_than_merged() -> None:
    """Two uploaded CIFs can legitimately carry the same formula as their name."""

    nickel = builtin_phase("ni_fcc").to_phase()
    report = identify_phase(
        _table("ni_fcc"), [("Ni", nickel), ("Ni", builtin_phase("cu_fcc").to_phase())]
    )
    assert sorted(item.phase_name for item in report) == ["Ni", "Ni (2)"]


def test_candidates_may_be_given_as_a_bare_sequence_of_phases() -> None:
    report = identify_phase(
        _table("ni_fcc"),
        [builtin_phase("ni_fcc").to_phase(), builtin_phase("fe_bcc").to_phase()],
    )
    assert report.best.phase_name == "Nickel (fcc)"


def test_sources_are_carried_into_the_report_for_provenance() -> None:
    report = identify_phase(
        _table("ni_fcc"),
        _candidates("ni_fcc"),
        sources={"Nickel (fcc)": "nickel.cif"},
    )
    assert report.best.source == "nickel.cif"
    assert "nickel.cif" in report.best.describe()


# ---------------------------------------------------------------------------
# Contracts and prose
# ---------------------------------------------------------------------------


def test_the_json_contract_carries_the_schema_and_every_criterion() -> None:
    report = identify_phase(_table("ni_fcc"), _candidates("ni_fcc", "fe_bcc"))
    payload = report.to_json()
    assert payload["schema"] == PHASE_IDENTIFICATION_SCHEMA
    assert payload["best_phase_name"] == "Nickel (fcc)"
    assert payload["candidate_count"] == 2
    assert payload["is_conclusive"] is True
    first = payload["candidates"][0]
    assert first["schema"] == PHASE_CANDIDATE_SCORE_SCHEMA
    assert set(first["criteria"]) == set(CRITERION_NAMES)
    assert first["figure_of_merit_m"]["value"] > 0.0
    assert first["indexing"]["phase_name"] == "Nickel (fcc)"


def test_describe_names_the_winner_the_margin_and_its_citations() -> None:
    report = identify_phase(_table("ni_fcc"), _candidates("ni_fcc", "fe_bcc"))
    prose = report.describe()
    assert "Nickel (fcc)" in prose
    assert "doi:10.1021/ac50125a001" in prose  # Hanawalt, Rinn & Frevel
    assert "doi:10.1107/S0021889886089458" in prose  # Dollase, on why intensities
    assert "preferred orientation" in prose


def test_candidate_describe_explains_a_loss_rather_than_only_stating_it() -> None:
    report = identify_phase(_table("ni_fcc"), _candidates("ni_fcc", "fe_bcc"))
    prose = next(
        item for item in report if item.phase_name == "Ferrite (bcc Fe)"
    ).describe()
    assert "predicts at" in prose and "not observed" in prose


# ---------------------------------------------------------------------------
# The whole pipeline, and input validation
# ---------------------------------------------------------------------------


def test_the_pipeline_detects_fits_and_ranks_in_one_call() -> None:
    report, table = identify_phase_from_pattern(
        _scan("ni_fcc"), _candidates("ni_fcc", "cu_fcc", "fe_bcc")
    )
    assert len(table) >= 5
    assert report.peak_count == len(table)
    assert report.best.phase_name == "Nickel (fcc)"
    # The table is returned rather than consumed: peak detection is where an
    # identification most often goes wrong, and it must remain inspectable.
    assert all(peak.converged for peak in table)


def test_an_empty_peak_table_is_refused() -> None:
    with pytest.raises(ValueError, match="empty peak table"):
        identify_phase(
            PeakTable(name="nothing", peaks=(), radiation=RadiationSpec.cu_ka()),
            _candidates("ni_fcc"),
        )


def test_a_table_without_radiation_is_refused() -> None:
    table = _table("ni_fcc")
    stripped = PeakTable(name=table.name, peaks=table.peaks, radiation=None)
    with pytest.raises(ValueError, match="needs a radiation"):
        identify_phase(stripped, _candidates("ni_fcc"))


def test_no_candidates_is_refused() -> None:
    with pytest.raises(ValueError, match="at least one candidate"):
        identify_phase(_table("ni_fcc"), {})


def test_a_non_phase_candidate_is_refused_with_a_pointer_to_the_cif_reader() -> None:
    with pytest.raises(TypeError, match=r"Phase.from_cif"):
        identify_phase(_table("ni_fcc"), {"not a phase": "Ni"})  # type: ignore[dict-item]


def test_an_impossible_tolerance_is_refused() -> None:
    with pytest.raises(ValueError, match="positive tolerance_deg"):
        identify_phase(_table("ni_fcc"), _candidates("ni_fcc"), tolerance_deg=0.0)


def test_an_impossible_strong_line_threshold_is_refused() -> None:
    with pytest.raises(ValueError, match="strong_line_threshold"):
        identify_phase(
            _table("ni_fcc"), _candidates("ni_fcc"), strong_line_threshold=1.5
        )


def test_an_identification_needs_at_least_one_candidate() -> None:
    with pytest.raises(ValueError, match="at least one candidate"):
        PhaseIdentification(name="empty", candidates=())


def test_a_candidate_without_an_indexing_must_state_why() -> None:
    with pytest.raises(ValueError, match="rejection reason"):
        PhaseCandidateScore(phase_name="stub", indexing=None)


# ---------------------------------------------------------------------------
# The cell-scale refinement: matching a tabulated CIF to a real specimen
# ---------------------------------------------------------------------------


def _dilated_scan(identifier: str, scale: float, *, seed: int = 11) -> PeakTable:
    """Return fitted peaks of a specimen whose cell differs from the tabulated one."""

    phase = _dilated(
        builtin_phase(identifier).to_phase(), scale, name=f"{identifier} specimen"
    )
    radiation = RadiationSpec.cu_ka()
    pattern = generate_xrd_pattern(
        phase,
        radiation=radiation,
        two_theta_range_deg=RANGE,
        resolution_deg=0.01,
        broadening_fwhm_deg=FWHM,
        profile="pseudo_voigt",
        max_index=6,
    )
    profile = np.asarray(pattern.intensity_grid, dtype=float)
    measured = MeasuredPowderPattern(
        name=f"{identifier} dilated by {scale}",
        two_theta_deg=np.asarray(pattern.two_theta_grid_deg, dtype=float),
        intensity=np.random.default_rng(seed)
        .poisson(profile / profile.max() * 30000.0 + 150.0)
        .astype(float),
        radiation=radiation,
        synthetic=True,
    )
    return detect_and_fit_peaks(measured, instrument=InstrumentBroadening.ideal(FWHM))


def test_a_specimen_whose_cell_differs_from_the_cif_is_still_identified() -> None:
    """The case every real identification meets: a solid solution against a tabulated cell.

    ``Delta(2 theta) = 2 e tan(theta)`` makes a three-parts-in-a-thousand cell
    difference -- ordinary for an alloy against a pure-element CIF -- displace a
    back-reflection line by more than half a degree. Without the scale
    refinement the true phase loses exactly the high-angle lines that would
    have confirmed it.
    """

    table = _dilated_scan("ni_fcc", 1.003)
    unrefined = identify_phase(
        table, _candidates("ni_fcc"), tolerance_deg=0.3, cell_scale_range=0.0
    ).best
    refined = identify_phase(table, _candidates("ni_fcc"), tolerance_deg=0.3).best

    assert refined.indexed_count > unrefined.indexed_count
    assert refined.explained_intensity_fraction > unrefined.explained_intensity_fraction
    assert refined.score > unrefined.score


def test_the_refined_scale_recovers_the_dilation_that_was_applied() -> None:
    """The reported factor is a measurement, so it is checked against the truth."""

    refined = identify_phase(
        _dilated_scan("ni_fcc", 1.003), _candidates("ni_fcc")
    ).best
    assert refined.cell_scale == pytest.approx(1.003, abs=5.0e-4)
    assert f"{refined.cell_scale:.5f}" in refined.describe()


def test_a_uniform_dilation_cannot_make_a_wrong_structure_fit() -> None:
    """Scaling preserves every ratio of d spacings, and ratios are what indexing tests.

    This is the property that makes the refinement safe. If a scale factor could
    rescue an arbitrary candidate the criterion would be worthless, so the fcc
    pattern is offered a body-centred and a rock-salt candidate with the full
    search range available to both.
    """

    report = identify_phase(
        _table("ni_fcc"), _candidates("ni_fcc", "fe_bcc", "nacl"), cell_scale_range=0.02
    )
    assert report.best.phase_name == "Nickel (fcc)"
    assert report.is_decisive
    for loser in report.candidates[1:]:
        assert loser.completeness < report.best.completeness


def test_the_refinement_is_off_when_the_range_is_zero() -> None:
    candidate = identify_phase(
        _table("ni_fcc"), _candidates("ni_fcc"), cell_scale_range=0.0
    ).best
    assert candidate.cell_scale == 1.0
    assert "dilated by a factor" not in candidate.describe()


def test_an_impossible_cell_scale_range_is_refused() -> None:
    with pytest.raises(ValueError, match="cell_scale_range"):
        identify_phase(_table("ni_fcc"), _candidates("ni_fcc"), cell_scale_range=1.5)


def test_the_refined_scale_travels_in_the_json_contract() -> None:
    payload = identify_phase(
        _dilated_scan("ni_fcc", 1.003), _candidates("ni_fcc")
    ).to_json()
    assert payload["candidates"][0]["cell_scale"] == pytest.approx(1.003, abs=5.0e-4)
    assert payload["settings"]["cell_scale_range"] == pytest.approx(0.02)
