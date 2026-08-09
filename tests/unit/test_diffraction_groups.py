"""Tests for `pytex.diffraction.diffraction_groups`.

The 31 diffraction groups and the table relating them to the 32 point groups are
*derived* in PyTex rather than transcribed, so the tests check the derivation
against what the literature reports rather than against a copy of it:

- The construction must produce exactly **31** groups from PyTex's own operator
  tables. Buxton, Eades, Steeds and Rackham (1976) derived 31; arriving at 30 or
  32 would mean the construction or the operators are wrong.
- Individual entries are checked against the ones every textbook treatment
  quotes: ``m-3m`` down a four-fold gives ``4mm1_R``, ``-43m`` down the same
  direction gives ``4_Rmm_R``, ``-6m2`` down its three-fold gives ``3m1_R`` with
  a *six*-fold bright-field disc, and so on.
- The centrosymmetry theorem is checked exhaustively: ``2_R`` appears at every
  beam direction of every centrosymmetric point group and at none of the others.
  That statement, over all 32 groups, is the capability the module exists for.
- The observation that decides centrosymmetry must be **only** the ``+-g``
  relation: knowing it alone must leave exactly the 21 non-centrosymmetric point
  groups, and not knowing it must leave the verdict open however much disc and
  pattern symmetry is supplied.
"""

from __future__ import annotations

import numpy as np
import pytest

from pytex.core.frame_catalog import crystal_frame
from pytex.core.lattice import AtomicSite, Lattice, Phase, SpaceGroupSpec, UnitCell, ZoneAxis
from pytex.core.point_groups import PointGroup, all_point_group_symbols
from pytex.core.symmetry import SymmetrySpec
from pytex.diffraction.diffraction_groups import (
    POINT_GROUP_DETERMINATION_SCHEMA,
    SymmetryObservations,
    determine_point_group,
    diffraction_group_for,
    diffraction_group_for_zone_axis,
    diffraction_group_symbols,
    diffraction_group_table,
)

#: The 31 diffraction groups of Buxton, Eades, Steeds and Rackham (1976),
#: *Phil. Trans. R. Soc. Lond. A* **281**, 171. Listed here as the count and
#: membership the derivation must reproduce, not as a table it reads from.
BUXTON_DIFFRACTION_GROUPS = frozenset(
    {
        "1",
        "1_R",
        "2",
        "2_R",
        "21_R",
        "m_R",
        "m",
        "m1_R",
        "2m_Rm_R",
        "2mm",
        "2_Rmm_R",
        "2mm1_R",
        "4",
        "4_R",
        "41_R",
        "4m_Rm_R",
        "4mm",
        "4_Rmm_R",
        "4mm1_R",
        "3",
        "31_R",
        "3m_R",
        "3m",
        "3m1_R",
        "6",
        "6_R",
        "61_R",
        "6m_Rm_R",
        "6mm",
        "6_Rmm_R",
        "6mm1_R",
    }
)

#: Entries every textbook treatment of CBED symmetry quotes, as
#: ``(point group, beam direction) -> (diffraction group, BF, WP)``.
CANONICAL_ENTRIES = (
    ("m-3m", (0, 0, 1), "4mm1_R", "4mm", "4mm"),
    ("m-3m", (1, 1, 1), "6_Rmm_R", "3m", "3m"),
    ("m-3m", (1, 1, 0), "2mm1_R", "2mm", "2mm"),
    ("432", (0, 0, 1), "4m_Rm_R", "4mm", "4"),
    ("432", (1, 1, 1), "3m_R", "3m", "3"),
    ("-43m", (0, 0, 1), "4_Rmm_R", "4mm", "2mm"),
    ("-43m", (1, 1, 1), "3m", "3m", "3m"),
    ("6/mmm", (0, 0, 1), "6mm1_R", "6mm", "6mm"),
    ("6mm", (0, 0, 1), "6mm", "6mm", "6mm"),
    ("-6m2", (0, 0, 1), "3m1_R", "6mm", "3m"),
    ("4/mmm", (0, 0, 1), "4mm1_R", "4mm", "4mm"),
    ("-42m", (0, 0, 1), "4_Rmm_R", "4mm", "2mm"),
    ("4mm", (0, 0, 1), "4mm", "4mm", "4mm"),
    ("-1", (0, 0, 1), "2_R", "1", "1"),
    ("1", (0, 0, 1), "1", "1", "1"),
)


# --------------------------------------------------------------------------- #
# The derivation reproduces the literature
# --------------------------------------------------------------------------- #


def test_the_construction_produces_exactly_the_thirty_one_diffraction_groups() -> None:
    """31 groups, derived from PyTex's operators, matching Buxton's set exactly.

    The count is the check that makes transcription unnecessary. Subgroups of a
    plane point group crossed with ``Z_2`` fall into three families — no tagged
    element (10), the full direct product (10), and the graph of a surjection
    onto ``Z_2`` (11) — and the scan over all 32 crystallographic point groups
    must realize every one of them.
    """

    derived = frozenset(diffraction_group_symbols())
    assert len(derived) == 31
    assert derived == BUXTON_DIFFRACTION_GROUPS


@pytest.mark.parametrize(
    ("point_group", "direction", "expected", "bright_field", "whole_pattern"),
    CANONICAL_ENTRIES,
    ids=[f"{entry[0]}-{entry[1]}" for entry in CANONICAL_ENTRIES],
)
def test_canonical_zone_axis_entries_match_the_literature(
    point_group: str,
    direction: tuple[int, int, int],
    expected: str,
    bright_field: str,
    whole_pattern: str,
) -> None:
    """The named cases every treatment of CBED symmetry quotes.

    ``-6m2`` down its three-fold is the sharp one: the bright-field disc shows a
    *six*-fold that the whole pattern does not, because the horizontal mirror
    contributes ``1_R`` and reciprocity supplies the extra inversion inside the
    direct disc. A construction that mixed up the bright-field homomorphism
    would report ``3m`` there and pass every other case in this list.
    """

    group = diffraction_group_for(point_group, direction)
    assert group.symbol == expected
    assert group.bright_field_symbol == bright_field
    assert group.whole_pattern_symbol == whole_pattern


def test_the_bright_field_symmetry_always_contains_the_whole_pattern_symmetry() -> None:
    """WP is the untagged subgroup; BF is the image of the whole group.

    Their orders must therefore satisfy ``|WP| <= |BF|``, with equality exactly
    when no tagged element contributes anything new. This holds for every point
    group at every characteristic direction, and it is the structural reason a
    CBED pattern carries more symmetry information than its spot positions do.
    """

    orders = {
        "1": 1, "2": 2, "3": 3, "4": 4, "6": 6,
        "m": 2, "2mm": 4, "3m": 6, "4mm": 8, "6mm": 12,
    }
    for symbol in all_point_group_symbols():
        for entry in diffraction_group_table(symbol):
            group = entry.diffraction_group
            assert orders[group.whole_pattern_symbol] <= orders[group.bright_field_symbol], (
                f"{symbol} at {entry.beam_direction}: WP {group.whole_pattern_symbol} is not "
                f"contained in BF {group.bright_field_symbol}"
            )


# --------------------------------------------------------------------------- #
# The centrosymmetry theorem
# --------------------------------------------------------------------------- #


def test_two_r_is_present_exactly_for_the_centrosymmetric_point_groups() -> None:
    """``2_R`` at every beam direction of a centric crystal, and at none of an acentric one.

    ``2_R`` demands an operator acting as ``-1`` on the beam direction *and* as
    ``-1`` on the transverse plane, and the inversion is the only operator that
    does both. So this is not a statistical tendency: it is an exact
    correspondence, and checking it over all 32 point groups at every
    characteristic direction is checking the capability the module exists for.
    """

    for symbol in all_point_group_symbols():
        centrosymmetric = PointGroup.from_symbol(symbol).is_centrosymmetric
        entries = diffraction_group_table(symbol)
        assert entries
        for entry in entries:
            assert entry.diffraction_group.has_friedel_symmetry is centrosymmetric, (
                f"{symbol} at {entry.beam_direction} gave "
                f"{entry.diffraction_group.symbol}"
            )


def test_two_r_is_invisible_in_bright_field_and_whole_pattern() -> None:
    """Why the ``+-g`` observation is indispensable rather than merely helpful.

    ``phi(2, tagged) = -2 = 1``, so the inversion contributes nothing to either
    disc symmetry. A crystal and its centrosymmetric relative can therefore show
    identical bright-field and whole-pattern symmetry, and be separated only by
    the two-fold relation between the ``+g`` and ``-g`` discs.
    """

    centric = diffraction_group_for("-1", (0, 0, 1))
    acentric = diffraction_group_for("1", (0, 0, 1))
    assert centric.bright_field_symbol == acentric.bright_field_symbol == "1"
    assert centric.whole_pattern_symbol == acentric.whole_pattern_symbol == "1"
    assert centric.has_friedel_symmetry is True
    assert acentric.has_friedel_symmetry is False


def test_one_r_comes_from_a_mirror_perpendicular_to_the_beam() -> None:
    """A horizontal mirror gives a two-fold in the direct disc alone.

    ``mm2`` down its polar axis has the two vertical mirrors and no horizontal
    one, so it shows ``2mm`` with no reciprocity element; ``m`` viewed along its
    mirror *normal* has the horizontal mirror and nothing else, and gives
    ``1_R``: a bright-field two-fold with no whole-pattern symmetry at all.
    """

    polar = diffraction_group_for("mm2", (0, 0, 1))
    assert polar.symbol == "2mm"
    assert polar.has_projection_reciprocity is False

    horizontal = diffraction_group_for("m", (0, 0, 1))
    assert horizontal.symbol == "1_R"
    assert horizontal.has_projection_reciprocity is True
    assert horizontal.bright_field_symbol == "2"
    assert horizontal.whole_pattern_symbol == "1"


def test_m_r_comes_from_a_two_fold_axis_perpendicular_to_the_beam() -> None:
    """A two-fold across the beam gives a bright-field mirror the pattern lacks."""

    group = diffraction_group_for("2", (1, 0, 0))
    assert group.symbol == "m_R"
    assert group.bright_field_symbol == "m"
    assert group.whole_pattern_symbol == "1"


# --------------------------------------------------------------------------- #
# The inverse problem
# --------------------------------------------------------------------------- #


def test_the_friedel_observation_alone_splits_the_thirty_two_point_groups() -> None:
    """21 acentric groups and 11 centric ones, from one observation.

    This is the arithmetic of the whole technique: 32 point groups fall into 11
    Laue classes, which is where selected-area diffraction stops, and the
    ``+-g`` observation recovers the distinction that Friedel's law destroyed.
    """

    acentric = determine_point_group(SymmetryObservations(friedel_pair_two_fold=False))
    centric = determine_point_group(SymmetryObservations(friedel_pair_two_fold=True))
    assert len(acentric.point_groups) == 21
    assert len(centric.point_groups) == 11
    assert acentric.is_centrosymmetric is False
    assert centric.is_centrosymmetric is True
    assert set(acentric.point_groups) | set(centric.point_groups) == set(all_point_group_symbols())
    assert not set(acentric.point_groups) & set(centric.point_groups)


def test_disc_symmetry_without_the_friedel_observation_leaves_the_centre_open() -> None:
    """``4mm`` bright field, ``2mm`` whole pattern: acentric, but which one?

    Both ``-42m`` and ``-43m`` produce ``4_Rmm_R``. The determination reports
    both, says the crystal is acentric because they agree on that, and says what
    would separate them rather than picking one.
    """

    result = determine_point_group(
        SymmetryObservations(bright_field="4mm", whole_pattern="2mm", friedel_pair_two_fold=False)
    )
    assert result.diffraction_groups == ("4_Rmm_R",)
    assert set(result.point_groups) == {"-42m", "-43m"}
    assert result.is_centrosymmetric is False
    assert result.is_unique is False
    assert "second zone axis" in result.describe()


def test_a_unique_determination_is_reported_as_such() -> None:
    """``6mm`` bright field with ``3m`` whole pattern can only be ``-6m2``."""

    result = determine_point_group(
        SymmetryObservations(bright_field="6mm", whole_pattern="3m")
    )
    assert result.diffraction_groups == ("3m1_R",)
    assert result.point_groups == ("-6m2",)
    assert result.is_unique is True
    assert result.is_centrosymmetric is False


def test_prior_knowledge_narrows_the_candidate_set() -> None:
    """Restricting to a known crystal system is usually what makes it unique."""

    observations = SymmetryObservations(
        bright_field="4mm", whole_pattern="2mm", friedel_pair_two_fold=False
    )
    tetragonal = determine_point_group(
        observations, candidate_point_groups=("4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm")
    )
    assert tetragonal.point_groups == ("-42m",)
    assert tetragonal.is_unique is True


def test_impossible_observations_are_reported_as_impossible() -> None:
    """A six-fold whole pattern inside a four-fold disc cannot happen."""

    result = determine_point_group(
        SymmetryObservations(bright_field="4mm", whole_pattern="6mm")
    )
    assert result.point_groups == ()
    assert result.is_centrosymmetric is None
    assert "No crystallographic point group" in result.describe()


def test_empty_observations_are_refused() -> None:
    """Every point group would be consistent, so the answer would have no content."""

    with pytest.raises(ValueError, match="At least one symmetry observation"):
        determine_point_group(SymmetryObservations())


def test_a_symmetry_that_is_not_a_plane_point_group_is_refused() -> None:
    """The observable is one of ten groups; anything else is a misreading."""

    with pytest.raises(ValueError, match="two-dimensional crystallographic point"):
        SymmetryObservations(bright_field="5mm")


# --------------------------------------------------------------------------- #
# Phases, zone axes, and reporting
# --------------------------------------------------------------------------- #


def _cubic_phase(name: str, point_group: str, space_group: tuple[str, int]) -> Phase:
    frame = crystal_frame()
    lattice = Lattice(5.6535, 5.6535, 5.6535, 90.0, 90.0, 90.0, crystal_frame=frame)
    sites = (
        AtomicSite(label="Ga", species="Ga", fractional_coordinates=np.zeros(3)),
        AtomicSite(
            label="As", species="As", fractional_coordinates=np.full(3, 0.25, dtype=np.float64)
        ),
    )
    return Phase(
        name,
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
        crystal_frame=frame,
        unit_cell=UnitCell(lattice=lattice, sites=sites),
        space_group=SpaceGroupSpec(
            symbol=space_group[0], number=space_group[1], reference_frame=frame
        ),
    )


def test_a_zone_axis_of_a_phase_resolves_to_its_diffraction_group() -> None:
    """The practical entry point keeps the beam direction in the crystal frame."""

    phase = _cubic_phase("gaas", "-43m", ("F-43m", 216))
    down_001 = diffraction_group_for_zone_axis(
        phase, ZoneAxis(indices=(0, 0, 1), phase=phase)
    )
    assert down_001.symbol == "4_Rmm_R"
    assert down_001.has_friedel_symmetry is False

    down_111 = diffraction_group_for_zone_axis(
        phase, ZoneAxis(indices=(1, 1, 1), phase=phase)
    )
    assert down_111.symbol == "3m"


def test_a_hexagonal_zone_axis_uses_the_lattice_cartesian_direction() -> None:
    """Hexagonal settings are where a fixed index list would go wrong.

    The in-plane two-folds of a hexagonal group sit at 0, 60 and 120 degrees in
    the Cartesian setting and are not integer triples, so the beam direction has
    to come from the lattice basis rather than from the Miller indices read as
    Cartesian components.
    """

    frame = crystal_frame()
    lattice = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=frame)
    phase = Phase(
        "zirconium-hcp",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=frame),
        crystal_frame=frame,
        unit_cell=UnitCell(
            lattice=lattice,
            sites=(
                AtomicSite(label="Zr1", species="Zr", fractional_coordinates=np.zeros(3)),
                AtomicSite(
                    label="Zr2",
                    species="Zr",
                    fractional_coordinates=np.array([1 / 3, 2 / 3, 0.5]),
                ),
            ),
        ),
        space_group=SpaceGroupSpec(symbol="P6_3/mmc", number=194, reference_frame=frame),
    )
    basal = diffraction_group_for_zone_axis(phase, ZoneAxis(indices=(0, 0, 1), phase=phase))
    assert basal.symbol == "6mm1_R"
    assert basal.bright_field_symbol == "6mm"
    assert basal.has_friedel_symmetry is True


def test_zone_axis_from_another_phase_is_refused() -> None:
    """Frame discipline: a zone axis carries the phase it belongs to."""

    first = _cubic_phase("gaas", "-43m", ("F-43m", 216))
    second = _cubic_phase("other", "m-3m", ("Fm-3m", 225))
    with pytest.raises(ValueError, match=r"zone_axis\.phase must match phase"):
        diffraction_group_for_zone_axis(first, ZoneAxis(indices=(0, 0, 1), phase=second))


def test_the_table_offers_a_beam_direction_for_every_reachable_group() -> None:
    """The forward table is what an experiment plan is built from."""

    entries = diffraction_group_table("-43m")
    symbols = [entry.diffraction_group.symbol for entry in entries]
    assert symbols[0] == "4_Rmm_R", "the most informative zone axis should come first"
    assert "3m" in symbols
    assert "1" in symbols
    for entry in entries:
        assert float(np.linalg.norm(entry.beam_direction)) == pytest.approx(1.0, abs=1e-12)
        recomputed = diffraction_group_for("-43m", entry.beam_direction)
        assert recomputed.symbol == entry.diffraction_group.symbol
    assert entries[0].is_special is True
    assert entries[-1].is_special is False


def test_describe_and_json_stay_in_lockstep() -> None:
    """The explainable-results contract, including the unimplemented observations."""

    group = diffraction_group_for("-43m", (0, 0, 1))
    prose = group.describe()
    assert "4_Rmm_R" in prose
    assert "reciprocity" in prose
    assert "no centre of symmetry" in prose

    result = determine_point_group(
        SymmetryObservations(bright_field="4mm", whole_pattern="2mm", friedel_pair_two_fold=False)
    )
    payload = result.to_json_dict()
    assert payload["schema"] == POINT_GROUP_DETERMINATION_SCHEMA
    assert payload["is_centrosymmetric"] is False
    assert payload["diffraction_groups"] == list(result.diffraction_groups)
    assert payload["point_groups"] == list(result.point_groups)
    assert payload["observations"]["friedel_pair_two_fold"] is False
    assert "not implemented here" in result.describe()

    unknown = determine_point_group(SymmetryObservations(bright_field="1"))
    assert unknown.is_centrosymmetric is None
    assert "must be read" in unknown.describe()
