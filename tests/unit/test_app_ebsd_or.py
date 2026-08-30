"""OR determination from measured grains, checked against constructed truth.

Nothing here compares against a previously recorded output of this code. Every
expectation is either a published value — the Burgers relationship's 12 variants
and its 45.29 degree disorientation — or a property of a measurement *built* to
have a known answer: pairs synthesised through a known variant must be named as
that relationship at essentially zero deviation, and the same pairs perturbed by
a known angle must move by about that angle and no more.

The construction is worth stating once. With a parent orientation ``P`` and a
variant rotation ``R``, the child orientation is ``C = P Rᵀ``: the
characterization forms each measured rotation as ``Cᵀ P``, so this is the
composition that makes ``Cᵀ P = R`` exactly. Getting it backwards produces a
plausible-looking set of angles that no relationship matches, which is how the
convention was pinned down in the first place.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.core.frame_catalog import specimen_frame
from pytex.core.orientation import Orientation, Rotation
from pytex.core.transformation import OrientationRelationship

BETA_ZIRCONIUM = {"builtin": "zr_bcc_beta"}
ZIRCONIUM = {"builtin": "zr_hcp"}
AUSTENITE = {"builtin": "austenite_fcc"}
FERRITE = {"builtin": "fe_bcc"}

#: The Burgers disorientation, in degrees. Burgers, Physica 1 (1934) 561.
BURGERS_DISORIENTATION_DEG = 45.29

#: The parent grain every constructed pair below descends from.
PARENT_EULER = (30.0, 40.0, 10.0)


def call(operation: str, **request: object) -> dict:
    return REGISTRY.call(operation, request)


def _index_count(label: str) -> int:
    """How many indices a label carries, whichever way it is written.

    The formatter writes ``(0001)`` when every index is a single digit and
    ``(1 -2 1 0)`` when one of them is not, so counting the tokens alone would
    call a four-index basal plane a one-index one.
    """

    body = label.strip("()[]<>{}")
    if " " in body:
        return len(body.split())
    return len(body.lstrip("-").replace("-", ""))


def child_euler(variant_index: int, parent: tuple[float, float, float] = PARENT_EULER):
    """The child Euler angles a given Burgers variant produces from ``parent``.

    One-based ``variant_index``, matching ``generate_variants()`` order.
    """

    _, parent_phase = phase_from_request(BETA_ZIRCONIUM)
    _, child_phase = phase_from_request(ZIRCONIUM)
    relationship = OrientationRelationship.from_burgers_correspondence(
        parent_phase=parent_phase, child_phase=child_phase
    )
    variant = relationship.generate_variants()[variant_index - 1]
    frame = specimen_frame()
    orientation = Orientation.from_euler(
        *parent,
        specimen_frame=frame,
        symmetry=parent_phase.symmetry,
        phase=parent_phase,
        convention="bunge",
        degrees=True,
    )
    matrix = np.asarray(orientation.rotation.as_matrix(), dtype=float)
    rotation = np.asarray(variant.parent_to_child_rotation.as_matrix(), dtype=float)
    return Rotation.from_matrix(matrix @ rotation.T).to_bunge_euler(degrees=True)


def pairs_text(variant_indices, *, noise_deg: float = 0.0) -> str:
    """A grain-pair table built through the named variants, optionally perturbed."""

    generator = np.random.default_rng(20260830)
    lines = []
    for index in variant_indices:
        angles = list(child_euler(index))
        if noise_deg:
            angles = [
                value + float(generator.normal(0.0, noise_deg)) for value in angles
            ]
        lines.append(
            " ".join(f"{value:.6f}" for value in (*PARENT_EULER, *angles))
        )
    return "\n".join(lines)


def run(pairs: str, **overrides: object) -> dict:
    request: dict[str, object] = {
        "phase": BETA_ZIRCONIUM,
        "child_phase": ZIRCONIUM,
        "pairs": pairs,
        "euler_convention": "bunge",
        "catalog_tolerance_deg": 3.0,
        "max_index": 3,
        "max_statements": 5,
    }
    request.update(overrides)
    return call("ebsd.or_from_grains", **request)


class TestDefaults:
    """The panel opens on an answer that can be checked."""

    def test_the_defaults_are_an_exact_burgers_measurement(self) -> None:
        """Three exact pairs, so the opening screen states a known result.

        A panel whose defaults produce a plausible but unverifiable number
        teaches the reader to trust an unverified answer; one that opens on a
        constructed case teaches them what a clean answer looks like.
        """

        operation = REGISTRY.get("ebsd.or_from_grains")
        data = REGISTRY.call(
            "ebsd.or_from_grains",
            {parameter.name: parameter.default for parameter in operation.parameters},
        )["data"]
        assert data["naming"]["best"] == "burgers"
        assert data["naming"]["is_conclusive"]
        assert data["naming"]["best_deviation_deg"] < 0.01
        assert data["fit"]["pair_count"] == 3

    def test_the_default_fit_is_the_published_burgers_disorientation(self) -> None:
        operation = REGISTRY.get("ebsd.or_from_grains")
        data = REGISTRY.call(
            "ebsd.or_from_grains",
            {parameter.name: parameter.default for parameter in operation.parameters},
        )["data"]
        assert data["fit"]["angle_deg"] == pytest.approx(BURGERS_DISORIENTATION_DEG, abs=0.01)


class TestNaming:
    """What the fit is called, and when the panel refuses to call it anything."""

    def test_pairs_from_one_variant_are_named_burgers(self) -> None:
        data = run(pairs_text([1]))["data"]
        assert data["naming"]["best"] == "burgers"
        assert data["naming"]["best_deviation_deg"] == pytest.approx(0.0, abs=1e-6)

    def test_pairs_from_different_variants_average_rather_than_fight(self) -> None:
        """The point of the symmetry reduction, and the reason to allow them.

        Variants differ by a parent symmetry operator, which the double-coset
        reduction absorbs before the pairs are compared. Three pairs on three
        different variants therefore fit *one* rotation with no scatter — if
        they did not, the panel would be unusable on any real map, where the
        product grains inside one parent are on different variants by
        definition.
        """

        data = run(pairs_text([1, 5, 9]))["data"]
        assert data["fit"]["pair_count"] == 3
        assert data["fit"]["mean_residual_deg"] == pytest.approx(0.0, abs=1e-6)
        assert data["naming"]["best"] == "burgers"

    def test_a_perturbed_measurement_still_names_burgers_and_gains_a_scatter(self) -> None:
        """Half a degree of orientation noise: the name survives, the scatter appears.

        Both halves matter. A tenth of a degree of noise must not unname a
        relationship whose runner-up is forty degrees away; and the scatter must
        stop being zero, because a zero scatter on noisy data would mean the
        panel was reporting the arithmetic rather than the measurement.
        """

        data = run(pairs_text([1, 5, 9], noise_deg=0.5))["data"]
        assert data["naming"]["best"] == "burgers"
        assert data["naming"]["is_conclusive"]
        assert 0.05 < data["fit"]["mean_residual_deg"] < 3.0

    def test_one_pair_has_no_scatter_whatever_the_grains_were(self) -> None:
        """The reason the panel takes a table and not six boxes.

        A single pair fits one rotation exactly, so its residual is zero by
        construction and carries no information about whether the relationship
        is real. The service must report that zero rather than concealing it,
        and the panel says beside it what it means.
        """

        data = run(pairs_text([7]))["data"]
        assert data["fit"]["pair_count"] == 1
        assert data["fit"]["mean_residual_deg"] == pytest.approx(0.0, abs=1e-9)

    def test_the_catalogue_is_the_one_for_these_two_crystal_systems(self) -> None:
        """bcc to hcp is ranked against Burgers and Shoji-Nishiyama, not martensite.

        A relationship defined between two cubic phases cannot be compared with
        a measurement of a hexagonal product, and offering it would invite a
        reader to conclude a bcc-to-hcp transformation was Kurdjumov-Sachs.
        """

        rows = run(pairs_text([1]))["table"]["rows"]
        names = {row["relationship"] for row in rows}
        assert names == {"Burgers", "Shoji-Nishiyama"}


class TestStatement:
    """The relationship written the way a paper writes it."""

    def test_the_statement_recovers_the_defining_parallelism(self) -> None:
        """Burgers is {110}_bcc || (0001)_hcp with <-111>_bcc || <11-20>_hcp.

        Recovered from the measurement rather than looked up, so this is a test
        of the rationalization and not of the catalogue: the fit knows nothing
        about which pair defined it.
        """

        statement = run(pairs_text([1, 5, 9]))["data"]["statement"]
        assert statement is not None
        assert statement["plane"]["child"] == "(0001)"
        assert statement["plane"]["deviation_deg"] < 0.01
        assert statement["cost_deg"] < 0.01

    def test_the_hexagonal_side_is_indexed_four_ways_throughout(self) -> None:
        """One crystal, one notation.

        The core labels both sides with three indices. A hexagonal child shown
        as (001) in one column and (0001) in the next is two notations for one
        plane, and a reader is entitled to conclude they are two planes.
        """

        data = run(pairs_text([1, 5, 9]))["data"]
        labels = [row["child"] for row in data["coincidences"]["planes"]]
        assert labels, "expected at least one candidate plane pair"
        for label in labels:
            assert _index_count(label) == 4, label

    def test_the_runners_up_are_reported_with_the_winner(self) -> None:
        """A winner alone cannot say whether it won clearly or out of a tie.

        The ranking is the evidence for the chosen pair, and it is what makes
        "best coincident direction" a finding rather than a choice.
        """

        data = run(pairs_text([1, 5, 9]), max_statements=5)["data"]
        directions = data["coincidences"]["directions"]
        assert len(directions) > 1
        deviations = [row["deviation_deg"] for row in directions]
        assert deviations == sorted(deviations)


class TestVariantAssignment:
    """Which admissible child each measured grain actually is."""

    def test_each_pair_is_assigned_a_variant_it_sits_on_exactly(self) -> None:
        data = run(pairs_text([1, 5, 9]))["data"]
        assert len(data["pairs"]) == 3
        for row in data["pairs"]:
            assert row["distance_deg"] < 1e-3
            assert 1 <= row["variant"] <= row["variant_count"]

    def test_distinct_variants_are_assigned_distinct_indices(self) -> None:
        """Three grains on three variants must not collapse to one index.

        The failure this guards against is an assignment that always returns
        the nearest-by-construction first variant, which would look right on a
        single-pair run and be useless on a map.
        """

        data = run(pairs_text([1, 5, 9]))["data"]
        assert len({row["variant"] for row in data["pairs"]}) == 3


class TestInput:
    """The table format, and what it says when it cannot read a line."""

    def test_comments_and_blank_lines_are_ignored(self) -> None:
        text = "# a header\n\n" + pairs_text([1]) + "\n\n"
        assert run(text)["data"]["fit"]["pair_count"] == 1

    def test_commas_count_as_separators(self) -> None:
        """A CSV pasted whole should work, because that is what people have."""

        text = ",".join(f"{value}" for value in (*PARENT_EULER, *child_euler(1)))
        assert run(text)["data"]["fit"]["pair_count"] == 1

    def test_a_short_row_is_refused_by_line_number(self) -> None:
        """"Line 2 has five numbers" is actionable; "invalid input" is not."""

        with pytest.raises(InvalidInputError) as error:
            run(pairs_text([1]) + "\n30 40 10 1 2")
        assert "Line 2" in str(error.value)
        assert error.value.details["field"] == "pairs"

    def test_a_non_numeric_row_is_refused_by_line_number(self) -> None:
        with pytest.raises(InvalidInputError) as error:
            run("30 40 10 nonsense 2 3")
        assert "Line 1" in str(error.value)

    def test_an_empty_table_is_refused(self) -> None:
        with pytest.raises(InvalidInputError) as error:
            run("# nothing but a comment")
        assert error.value.details["field"] == "pairs"

    def test_two_identical_phases_are_refused_with_the_reason(self) -> None:
        with pytest.raises(InvalidInputError) as error:
            run(pairs_text([1]), child_phase=BETA_ZIRCONIUM)
        assert "distinct phases" in (error.value.hint or "")


class TestOtherSystems:
    """The panel is not Burgers-only; the canonical case is only the default."""

    def test_an_fcc_to_bcc_pair_is_named_from_the_martensite_catalogue(self) -> None:
        _, parent_phase = phase_from_request(AUSTENITE)
        _, child_phase = phase_from_request(FERRITE)
        relationship = OrientationRelationship.from_kurdjumov_sachs_correspondence(
            parent_phase=parent_phase, child_phase=child_phase
        )
        variant = relationship.generate_variants()[0]
        frame = specimen_frame()
        parent = Orientation.from_euler(
            *PARENT_EULER,
            specimen_frame=frame,
            symmetry=parent_phase.symmetry,
            phase=parent_phase,
            convention="bunge",
            degrees=True,
        )
        matrix = np.asarray(parent.rotation.as_matrix(), dtype=float)
        rotation = np.asarray(variant.parent_to_child_rotation.as_matrix(), dtype=float)
        angles = Rotation.from_matrix(matrix @ rotation.T).to_bunge_euler(degrees=True)
        text = " ".join(f"{value:.6f}" for value in (*PARENT_EULER, *angles))
        data = run(text, phase=AUSTENITE, child_phase=FERRITE)["data"]
        assert data["naming"]["best"] == "kurdjumov_sachs"
        assert data["fit"]["angle_deg"] == pytest.approx(42.85, abs=0.01)


class TestAxisNaming:
    """The fitted axis, named in both crystal bases."""

    def test_the_axis_is_named_against_both_bases_with_its_residual(self) -> None:
        """One physical axis, two indexings, and neither claimed to be exact.

        The axis has the same Cartesian components in both crystal frames — it
        is the fixed vector of the map between them — so what differs is only
        how it is indexed. An OR axis is not in general rational in either
        basis, so a label without its residual would assert something false.
        """

        fit = run(pairs_text([1, 5, 9]))["data"]["fit"]
        for entry in (fit["axis_parent"], fit["axis_child"]):
            assert entry["label"]
            assert math.isfinite(entry["deviation_deg"])
            assert entry["deviation_deg"] >= 0.0
        assert fit["axis_parent"]["cartesian"] == fit["axis_child"]["cartesian"]
        # Three indices against the cubic basis, four against the hexagonal one.
        assert _index_count(fit["axis_parent"]["label"]) == 3
        assert _index_count(fit["axis_child"]["label"]) == 4
