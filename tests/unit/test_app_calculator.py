"""The calculator is checked against answers that are known without it.

Every assertion here is an angle, a multiplicity, or a spacing that a textbook
or an exact geometric argument fixes independently: 45° between (100) and (110)
in any cubic crystal, 90° between basal and prism planes in any hexagonal one,
the fcc powder sequence 111/200/220/311 with multiplicities 8/6/12/24. Nothing
is compared against a previously recorded output of this code.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.phases import BUILTIN_PHASES, PhaseSpec, builtin_phase
from pytex.app.services.calculator import KS_VARIANT_ANGLE_DEG, KS_VARIANT_AXIS

#: An explicit null rotation, for tests that mean "no relative orientation".
#: The operation's own default is the Kurdjumov-Sachs variant, not the identity,
#: so a test that wants the identity has to say so.
NO_ROTATION: dict[str, object] = {
    "rotation_convention": "bunge",
    "rotation_1": 0.0,
    "rotation_2": 0.0,
    "rotation_3": 0.0,
    "rotation_angle_deg": 0.0,
}

NICKEL = {"builtin": "ni_fcc"}
ZIRCONIUM = {"builtin": "zr_hcp"}
BETA_ZIRCONIUM = {"builtin": "zr_bcc_beta"}


def call(operation: str, **request: object) -> dict:
    return REGISTRY.call(operation, request)


def defaults_for(operation: str) -> dict[str, object]:
    """Return the request an untouched control panel would send.

    This is what the first press of the button submits, which is the state most
    users see and the one least often tested.
    """

    spec = REGISTRY.get(operation)
    return {
        parameter.name: parameter.default
        for parameter in spec.parameters
        if parameter.default is not None
    }


def angle_between(rows: list[dict], left: str, right: str) -> float:
    for row in rows:
        if row["left"] == left and row["right"] == right:
            return float(row["angle_deg"])
    raise AssertionError(f"no row for {left} against {right} in {rows}")


class TestPhaseSpec:
    """The specification is where bad materials are stopped."""

    def test_builtin_phases_all_build(self) -> None:
        for identifier in BUILTIN_PHASES:
            phase = builtin_phase(identifier).to_phase()
            assert phase.lattice.a > 0.0

    def test_cubic_phase_rejects_unequal_edges(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            PhaseSpec(
                name="not cubic",
                a=3.0,
                b=4.0,
                c=3.0,
                alpha=90.0,
                beta=90.0,
                gamma=90.0,
                point_group="m-3m",
            )
        assert excinfo.value.details["field"] == "b"

    def test_hexagonal_phase_requires_gamma_120(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            PhaseSpec(
                name="not hexagonal",
                a=3.0,
                b=3.0,
                c=5.0,
                alpha=90.0,
                beta=90.0,
                gamma=90.0,
                point_group="6/mmm",
            )
        assert excinfo.value.details["field"] == "gamma"

    def test_unknown_point_group_is_rejected(self) -> None:
        with pytest.raises(InvalidInputError, match="32 crystallographic point groups"):
            PhaseSpec(
                name="nonsense",
                a=3.0,
                b=3.0,
                c=3.0,
                alpha=90.0,
                beta=90.0,
                gamma=90.0,
                point_group="m-3n",
            )

    def test_a_space_group_from_another_crystal_system_is_refused(self) -> None:
        """The pairing a user reaches by editing a catalogue phase.

        Start from nickel, set the point group to 4/mmm and the cell to
        4 x 4 x 6, and the space group inherited from the entry is still cubic
        Fm-3m. Accepted, its F centring deletes every mixed-parity family from a
        tetragonal phase that has no F centring at all — a wrong pattern with
        nothing on screen to suggest it.
        """

        with pytest.raises(InvalidInputError) as excinfo:
            PhaseSpec(
                name="tetragonal with a cubic space group",
                a=4.0,
                b=4.0,
                c=6.0,
                alpha=90.0,
                beta=90.0,
                gamma=90.0,
                point_group="4/mmm",
                space_group_symbol="Fm-3m",
                space_group_number=225,
            )
        assert excinfo.value.details["field"] == "space_group_number"
        assert "cubic" in str(excinfo.value) and "tetragonal" in str(excinfo.value)

    def test_every_builtin_phase_agrees_with_its_own_space_group(self) -> None:
        for identifier in BUILTIN_PHASES:
            spec = builtin_phase(identifier)
            assert spec.point_group  # constructed, therefore already validated

    def test_a_space_group_of_the_stated_system_is_accepted(self) -> None:
        spec = PhaseSpec(
            name="rutile-like",
            a=4.594,
            b=4.594,
            c=2.959,
            alpha=90.0,
            beta=90.0,
            gamma=90.0,
            point_group="4/mmm",
            space_group_symbol="P4_2/mnm",
            space_group_number=136,
        )
        assert spec.crystal_system == "tetragonal"

    def test_space_group_symbol_and_number_travel_together(self) -> None:
        with pytest.raises(InvalidInputError, match="together"):
            PhaseSpec(
                name="half a space group",
                a=3.0,
                b=3.0,
                c=3.0,
                alpha=90.0,
                beta=90.0,
                gamma=90.0,
                point_group="m-3m",
                space_group_symbol="Fm-3m",
            )

    def test_builtin_can_be_overridden_without_retyping_the_structure(self) -> None:
        spec = PhaseSpec.from_json({"builtin": "ni_fcc", "a": 3.6, "b": 3.6, "c": 3.6})
        assert spec.a == pytest.approx(3.6)
        assert len(spec.sites) == 4
        assert spec.point_group == "m-3m"

    def test_round_trips_through_json(self) -> None:
        original = builtin_phase("zr_hcp")
        restored = PhaseSpec.from_json(original.to_json())
        assert restored == original

    def test_negative_cell_edge_names_the_field(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            PhaseSpec(
                name="impossible",
                a=-1.0,
                b=1.0,
                c=1.0,
                alpha=90.0,
                beta=90.0,
                gamma=90.0,
                point_group="1",
            )
        assert excinfo.value.details["field"] == "a"


class TestInterplanarAngles:
    """Cubic and hexagonal angles that hold for every lattice parameter."""

    def test_cubic_100_110_is_exactly_45_degrees(self) -> None:
        result = call("calc.plane_angles", phase=NICKEL, planes=[[1, 0, 0], [1, 1, 0]])
        assert angle_between(result["table"]["rows"], "(100)", "(110)") == pytest.approx(45.0)

    def test_cubic_100_111_is_the_body_diagonal_angle(self) -> None:
        result = call("calc.plane_angles", phase=NICKEL, planes=[[1, 0, 0], [1, 1, 1]])
        expected = math.degrees(math.acos(1.0 / math.sqrt(3.0)))
        angle = angle_between(result["table"]["rows"], "(100)", "(111)")
        assert angle == pytest.approx(expected, abs=1e-9)

    def test_cubic_angles_do_not_depend_on_the_lattice_parameter(self) -> None:
        stretched = {"builtin": "ni_fcc", "a": 7.0, "b": 7.0, "c": 7.0}
        first = call("calc.plane_angles", phase=NICKEL, planes=[[1, 0, 0], [1, 1, 1]])
        second = call("calc.plane_angles", phase=stretched, planes=[[1, 0, 0], [1, 1, 1]])
        assert angle_between(first["table"]["rows"], "(100)", "(111)") == pytest.approx(
            angle_between(second["table"]["rows"], "(100)", "(111)")
        )

    def test_hexagonal_basal_and_prism_are_perpendicular(self) -> None:
        result = call("calc.plane_angles", phase=ZIRCONIUM, planes=[[0, 0, 1], [1, 0, 0]])
        row = result["table"]["rows"][0]
        assert float(row["angle_deg"]) == pytest.approx(90.0)
        assert row["left"] == "(0001)"

    def test_spacing_columns_match_the_bragg_planes(self) -> None:
        result = call("calc.plane_angles", phase=NICKEL, planes=[[1, 0, 0], [2, 0, 0]])
        rows = result["table"]["rows"]
        assert float(rows[0]["left_d_angstrom"]) == pytest.approx(
            2.0 * float(rows[0]["right_d_angstrom"])
        )

    def test_a_single_plane_is_refused_with_an_explanation(self) -> None:
        with pytest.raises(InvalidInputError, match="at least two planes"):
            call("calc.plane_angles", phase=NICKEL, planes=[[1, 1, 1]])

    def test_rectangular_comparison_keeps_every_pair(self) -> None:
        result = call(
            "calc.plane_angles",
            phase=NICKEL,
            planes=[[1, 0, 0], [1, 1, 0]],
            against=[[1, 1, 1]],
        )
        assert len(result["table"]["rows"]) == 2


class TestDirectionAngles:
    def test_cubic_directions_match_the_plane_result(self) -> None:
        planes = call("calc.plane_angles", phase=NICKEL, planes=[[1, 0, 0], [1, 1, 1]])
        directions = call("calc.direction_angles", phase=NICKEL, directions=[[1, 0, 0], [1, 1, 1]])
        assert float(directions["table"]["rows"][0]["angle_deg"]) == pytest.approx(
            float(planes["table"]["rows"][0]["angle_deg"])
        )

    def test_hexagonal_direction_angle_differs_from_the_plane_angle(self) -> None:
        planes = call("calc.plane_angles", phase=ZIRCONIUM, planes=[[1, 0, 0], [1, 0, 1]])
        directions = call(
            "calc.direction_angles", phase=ZIRCONIUM, directions=[[1, 0, 0], [1, 0, 1]]
        )
        assert float(planes["table"]["rows"][0]["angle_deg"]) != pytest.approx(
            float(directions["table"]["rows"][0]["angle_deg"]), abs=1e-6
        )

    def test_antipodal_off_allows_obtuse_angles(self) -> None:
        result = call(
            "calc.direction_angles",
            phase=NICKEL,
            directions=[[1, 0, 0], [-1, 0, 0]],
            antipodal=False,
        )
        assert float(result["table"]["rows"][0]["angle_deg"]) == pytest.approx(180.0)


class TestPlaneDirectionGeometry:
    def test_zone_law_zero_means_the_direction_lies_in_the_plane(self) -> None:
        result = call(
            "calc.plane_direction_angles",
            phase=NICKEL,
            planes=[[1, 1, 1]],
            directions=[[1, -1, 0], [1, 1, 1]],
        )
        rows = {row["direction"]: row for row in result["table"]["rows"]}
        in_plane = rows["[1 -1 0]"]
        assert in_plane["zone_law"] == 0
        assert in_plane["in_zone"] is True
        assert float(in_plane["normal_angle_deg"]) == pytest.approx(90.0)
        assert float(in_plane["inclination_deg"]) == pytest.approx(0.0, abs=1e-9)

    def test_normal_and_inclination_are_complementary(self) -> None:
        result = call(
            "calc.plane_direction_angles",
            phase=ZIRCONIUM,
            planes=[[1, 0, 1]],
            directions=[[2, 1, 3]],
        )
        row = result["table"]["rows"][0]
        assert float(row["normal_angle_deg"]) + float(row["inclination_deg"]) == pytest.approx(90.0)


class TestSymmetryFamilies:
    def test_cubic_100_family_has_three_antipodal_members(self) -> None:
        result = call("calc.symmetry_family", phase=NICKEL, indices=[1, 0, 0])
        assert result["data"]["multiplicity"] == 3
        assert result["data"]["point_group_order"] == 48
        assert result["data"]["self_mapping_operations"] == 16

    def test_cubic_111_family_has_four_antipodal_members(self) -> None:
        result = call("calc.symmetry_family", phase=NICKEL, indices=[1, 1, 1])
        assert result["data"]["multiplicity"] == 4

    def test_family_members_share_one_spacing(self) -> None:
        result = call("calc.symmetry_family", phase=ZIRCONIUM, indices=[1, 0, 0])
        assert result["data"]["multiplicity"] == 3
        assert result["data"]["d_spacing_angstrom"] == pytest.approx(
            builtin_phase("zr_hcp").a * math.sqrt(3.0) / 2.0
        )

    def test_hexagonal_families_are_labelled_with_four_indices(self) -> None:
        result = call("calc.symmetry_family", phase=ZIRCONIUM, indices=[1, 0, 0])
        assert result["data"]["family_label"] == "{1 0 -1 0}"

    def test_direction_family_without_antipodal_equivalence_doubles(self) -> None:
        folded = call("calc.symmetry_family", phase=NICKEL, family="direction", indices=[1, 0, 0])
        unfolded = call(
            "calc.symmetry_family",
            phase=NICKEL,
            family="direction",
            indices=[1, 0, 0],
            antipodal=False,
        )
        assert unfolded["data"]["multiplicity"] == 2 * folded["data"]["multiplicity"]


class TestZoneAxis:
    def test_zone_axis_is_the_cross_product_reduced(self) -> None:
        result = call("calc.zone_axis", phase=NICKEL, first=[1, 1, 1], second=[1, -1, 0])
        assert result["data"]["zone_axis"] in ([1, 1, -2], [-1, -1, 2])

    def test_every_listed_reflection_satisfies_the_zone_law(self) -> None:
        result = call("calc.zone_axis", phase=NICKEL, first=[1, 1, 1], second=[1, -1, 0])
        axis = result["data"]["zone_axis"]
        for row in result["table"]["rows"]:
            assert row["h"] * axis[0] + row["k"] * axis[1] + row["l"] * axis[2] == 0

    def test_parallel_planes_are_refused(self) -> None:
        with pytest.raises(InvalidInputError, match="parallel"):
            call("calc.zone_axis", phase=NICKEL, first=[1, 1, 1], second=[2, 2, 2])

    def test_reflections_are_sorted_by_decreasing_spacing(self) -> None:
        result = call("calc.zone_axis", phase=NICKEL, first=[1, 0, 0], second=[0, 1, 0])
        spacings = [float(row["d_angstrom"]) for row in result["table"]["rows"]]
        assert spacings == sorted(spacings, reverse=True)


class TestDSpacings:
    def test_fcc_gives_the_textbook_powder_sequence(self) -> None:
        result = call("calc.d_spacings", phase=NICKEL, max_index=4, min_d_angstrom=0.6)
        families = [row["family"] for row in result["table"]["rows"]]
        assert families[:5] == ["{111}", "{200}", "{220}", "{311}", "{222}"]

    def test_fcc_multiplicities_are_the_published_ones(self) -> None:
        result = call("calc.d_spacings", phase=NICKEL, max_index=4, min_d_angstrom=0.6)
        multiplicity = {row["family"]: row["multiplicity"] for row in result["table"]["rows"]}
        assert multiplicity["{111}"] == 8
        assert multiplicity["{200}"] == 6
        assert multiplicity["{220}"] == 12
        assert multiplicity["{311}"] == 24

    def test_fcc_spacings_follow_the_cubic_law(self) -> None:
        result = call("calc.d_spacings", phase=NICKEL, max_index=2, min_d_angstrom=1.0)
        a = builtin_phase("ni_fcc").a
        for row in result["table"]["rows"]:
            squared = row["h"] ** 2 + row["k"] ** 2 + row["l"] ** 2
            assert float(row["d_angstrom"]) == pytest.approx(a / math.sqrt(squared))

    def test_bcc_removes_the_odd_sum_reflections(self) -> None:
        result = call(
            "calc.d_spacings", phase={"builtin": "fe_bcc"}, max_index=3, min_d_angstrom=0.8
        )
        for row in result["table"]["rows"]:
            assert (row["h"] + row["k"] + row["l"]) % 2 == 0

    def test_fcc_keeps_only_unmixed_parity(self) -> None:
        result = call("calc.d_spacings", phase=NICKEL, max_index=3, min_d_angstrom=0.8)
        for row in result["table"]["rows"]:
            parities = {row["h"] % 2, row["k"] % 2, row["l"] % 2}
            assert len(parities) == 1

    def test_a_phase_without_a_space_group_says_so(self) -> None:
        phase = {
            "name": "bare cubic",
            "a": 4.0,
            "b": 4.0,
            "c": 4.0,
            "alpha": 90.0,
            "beta": 90.0,
            "gamma": 90.0,
            "point_group": "m-3m",
        }
        result = call("calc.d_spacings", phase=phase, max_index=2, min_d_angstrom=1.0)
        assert any("No space group" in note for note in result["notes"])

    def test_an_impossible_spacing_filter_is_reported(self) -> None:
        with pytest.raises(InvalidInputError, match="No reflection"):
            call("calc.d_spacings", phase=NICKEL, max_index=1, min_d_angstrom=99.0)


class TestInterphaseAngles:
    """Angles across a stated orientation relationship."""

    def test_identity_rotation_reproduces_the_single_phase_answer(self) -> None:
        result = call(
            "calc.interphase_angles",
            phase=NICKEL,
            other_phase=NICKEL,
            first_indices=[[1, 0, 0]],
            second_indices=[[1, 1, 0]],
            **NO_ROTATION,
        )
        assert float(result["table"]["rows"][0]["angle_deg"]) == pytest.approx(45.0)

    def test_a_rotation_about_the_shared_normal_leaves_it_parallel(self) -> None:
        result = call(
            "calc.interphase_angles",
            phase=NICKEL,
            other_phase=NICKEL,
            first_indices=[[0, 0, 1]],
            second_indices=[[0, 0, 1]],
            rotation_convention="axis_angle",
            rotation_1=0.0,
            rotation_2=0.0,
            rotation_3=1.0,
            rotation_angle_deg=37.0,
        )
        assert float(result["table"]["rows"][0]["angle_deg"]) == pytest.approx(0.0, abs=1e-9)

    def test_reported_rotation_angle_matches_the_request(self) -> None:
        result = call(
            "calc.interphase_angles",
            phase={"builtin": "austenite_fcc"},
            other_phase={"builtin": "fe_bcc"},
            first_indices=[[1, 1, 1]],
            second_indices=[[1, 1, 0]],
            rotation_convention="axis_angle",
            rotation_1=1.0,
            rotation_2=1.0,
            rotation_3=1.0,
            rotation_angle_deg=42.85,
        )
        assert float(result["data"]["rotation_angle_deg"]) == pytest.approx(42.85, abs=1e-6)

    def test_kurdjumov_sachs_default_matches_the_computed_relationship(self) -> None:
        """The opening rotation is the KS misorientation, recomputed here.

        The parameter defaults are literals because a manifest is static, so this
        test is what keeps them honest: it builds the relationship from
        :mod:`pytex.core.transformation` and fails if the literals drift from
        what the library computes.
        """

        from pytex.core.transformation import OrientationRelationship

        rotation = (
            OrientationRelationship.from_kurdjumov_sachs_correspondence(
                parent_phase=builtin_phase("austenite_fcc").to_phase(),
                child_phase=builtin_phase("fe_bcc").to_phase(),
            )
            .generate_variants()[0]
            .parent_to_child_rotation
        )
        assert float(rotation.angle_deg) == pytest.approx(KS_VARIANT_ANGLE_DEG, abs=1e-5)
        # The default is the inverse rotation, written by negating the axis.
        inverse_axis = -np.asarray(rotation.axis, dtype=float)
        assert inverse_axis == pytest.approx(np.asarray(KS_VARIANT_AXIS), abs=1e-5)

    def test_the_opening_press_compares_two_different_phases(self) -> None:
        """Running from the declared defaults must not compare a phase with itself.

        The defaults exist so that the first press answers the question the help
        text poses: does Kurdjumov-Sachs really put a {111} of austenite on a
        {110} of ferrite? It does — the residual here is 1e-5 degrees, which is
        the six-decimal rounding of the stored axis and nothing else.
        """

        result = call("calc.interphase_angles", **defaults_for("calc.interphase_angles"))
        assert result["title"] == "Austenite (fcc Fe) against Ferrite (bcc Fe)"
        top = result["table"]["rows"][0]
        assert (top["first"], top["second"]) == ("(111)", "(011)")
        assert float(top["angle_deg"]) == pytest.approx(0.0, abs=1e-4)

    def test_a_zero_axis_is_refused(self) -> None:
        with pytest.raises(InvalidInputError, match="zero vector"):
            call(
                "calc.interphase_angles",
                phase=NICKEL,
                first_indices=[[1, 0, 0]],
                second_indices=[[1, 0, 0]],
                rotation_convention="axis_angle",
                rotation_1=0.0,
                rotation_2=0.0,
                rotation_3=0.0,
                rotation_angle_deg=10.0,
            )


class TestResultShape:
    """Every result can be read as prose and exported as rows."""

    @pytest.mark.parametrize(
        ("operation", "payload"),
        [
            ("calc.catalog", {}),
            ("calc.phase_summary", {"phase": NICKEL}),
            ("calc.plane_angles", {"phase": NICKEL}),
            ("calc.direction_angles", {"phase": NICKEL}),
            ("calc.plane_direction_angles", {"phase": NICKEL}),
            ("calc.symmetry_family", {"phase": NICKEL}),
            ("calc.zone_axis", {"phase": NICKEL}),
            ("calc.d_spacings", {"phase": NICKEL}),
        ],
    )
    def test_result_carries_a_summary_and_a_table(self, operation: str, payload: dict) -> None:
        result = REGISTRY.call(operation, payload)
        assert len(result["summary"]) > 40
        assert result["table"]["columns"]
        for row in result["table"]["rows"]:
            assert set(column["key"] for column in result["table"]["columns"]) <= set(row)


class TestOrientationRelationships:
    """Variant counts and misorientations that the literature fixes."""

    def test_kurdjumov_sachs_gives_24_variants_at_42_85_degrees(self) -> None:
        result = call(
            "calc.orientation_relationship",
            phase={"builtin": "austenite_fcc"},
            child_phase={"builtin": "fe_bcc"},
            relationship="kurdjumov_sachs",
        )
        assert result["data"]["variant_count"] == 24
        assert float(result["data"]["misorientation_angle_deg"]) == pytest.approx(42.85, abs=0.01)

    def test_nishiyama_wassermann_gives_12_variants(self) -> None:
        result = call(
            "calc.orientation_relationship",
            phase={"builtin": "austenite_fcc"},
            child_phase={"builtin": "fe_bcc"},
            relationship="nishiyama_wassermann",
        )
        assert result["data"]["variant_count"] == 12

    def test_bain_gives_three_variants(self) -> None:
        result = call(
            "calc.orientation_relationship",
            phase={"builtin": "austenite_fcc"},
            child_phase={"builtin": "fe_bcc"},
            relationship="bain",
        )
        assert result["data"]["variant_count"] == 3

    def test_burgers_gives_12_variants_into_a_hexagonal_child(self) -> None:
        result = call(
            "calc.orientation_relationship",
            phase=BETA_ZIRCONIUM,
            child_phase={"builtin": "zr_hcp"},
            relationship="burgers",
        )
        assert result["data"]["variant_count"] == 12
        assert result["data"]["parallel_planes"][0]["child"] == "(0001)"

    def test_the_defining_parallelism_is_reported(self) -> None:
        result = call(
            "calc.orientation_relationship",
            phase={"builtin": "austenite_fcc"},
            child_phase={"builtin": "fe_bcc"},
        )
        planes = result["data"]["parallel_planes"][0]
        assert planes["parent"] == "(111)"
        assert planes["child"] == "(011)"

    def test_an_inapplicable_relationship_explains_itself(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            call(
                "calc.orientation_relationship",
                phase={"builtin": "zr_hcp"},
                child_phase={"builtin": "zr_hcp"},
                relationship="kurdjumov_sachs",
            )
        assert excinfo.value.details["field"] == "relationship"
