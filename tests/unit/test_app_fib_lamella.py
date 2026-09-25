"""The workbench operation `fib.lamella_plan`.

What is pinned here is the contract between the library and the page: the result
states the answer in the operator's terms (grain, residual with its uncertainty,
azimuths in three frames, feasibility), carries every figure of section 9, the
work order, and the report JSON that validates against its schema, and refuses
bad input with a message that names the field.
"""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest

import pytex.app.services  # noqa: F401  (registers the operations)
from pytex.app.contracts import dumps
from pytex.app.errors import InvalidInputError
from pytex.app.registry import REGISTRY

SCHEMA = json.loads(
    (Path(__file__).resolve().parents[2] / "schemas" / "fib_lamella_plan.schema.json").read_text(
        encoding="utf-8"
    )
)


@pytest.fixture(scope="module")
def default_result() -> dict:
    return json.loads(dumps(REGISTRY.call("fib.lamella_plan", {})))


class TestDefaultRun:
    def test_the_result_names_the_answer(self, default_result) -> None:
        assert default_result["title"].startswith("FIB lamella for <011>")
        labels = [metric["label"] for metric in default_result["highlights"]]
        assert labels[0] == "Recommended grain"
        assert any(label.startswith("Residual tilt") for label in labels)
        assert any(label.startswith("FIB pattern rotation") for label in labels)
        assert "Feasibility" in labels

    def test_an_uncalibrated_chamber_is_a_warning_on_the_result(self, default_result) -> None:
        assert any("uncalibrated" in warning for warning in default_result["warnings"])

    def test_every_figure_of_section_nine_is_present(self, default_result) -> None:
        keys = {
            figure["key"]
            for stage in default_result["stages"]
            for figure in stage.get("figures", [])
        }
        assert keys == {
            "fib_plan_view",
            "fib_stereogram",
            "fib_phi_feasibility",
            "fib_section",
            "fib_preparability",
            "fib_saed",
        }

    def test_the_report_validates_against_its_schema(self, default_result) -> None:
        jsonschema.validate(default_result["data"]["report"], SCHEMA)

    def test_the_work_order_is_a_printable_page(self, default_result) -> None:
        page = default_result["data"]["work_order_html"]
        assert page.startswith("<!DOCTYPE html>")
        assert "@page" in page and "UNCALIBRATED" in page
        assert "<svg" in page  # the plan view and section travel with it

    def test_the_table_ranks_feasible_before_infeasible(self, default_result) -> None:
        order = {"guaranteed": 0, "probabilistic": 1, "unreachable": 2}
        classes = [order[row["feasibility"]] for row in default_result["table"]["rows"]]
        assert classes == sorted(classes)
        assert [row["rank"] for row in default_result["table"]["rows"]] == list(
            range(1, len(classes) + 1)
        )

    def test_the_stages_follow_the_report_order(self, default_result) -> None:
        sections = [stage.get("section") for stage in default_result["stages"]]
        order = ["result", "evidence", "diagnostics", "method", "audit"]
        assert [order.index(section) for section in sections] == sorted(
            order.index(section) for section in sections
        )


class TestSettings:
    def test_a_calibrated_chamber_moves_the_fib_rotation_and_clears_the_warning(self) -> None:
        request = {
            "zone_axis": [0, 0, 1],
            "grain_id": 0,
            "rotation_sense": "-1",
            "rotation_offset_deg": 90.0,
            "chamber_calibrated": True,
        }
        result = REGISTRY.call("fib.lamella_plan", request)
        plan = result["data"]["report"]["plans"][0]
        expected = (-(plan["geometry"]["theta_sample_deg"]) + 90.0) % 360.0
        assert plan["theta_ion_deg"] == pytest.approx(expected)
        assert not any("uncalibrated" in warning for warning in result.get("warnings", []))
        assert [row["grain_id"] for row in result["table"]["rows"]] == [0]

    def test_an_affine_registration_widens_the_budget(self) -> None:
        points = "0 0 20 10\n60 0 139 25\n0 60 5 130\n60 60 125 146\n30 30 72 79"
        result = REGISTRY.call(
            "fib.lamella_plan", {"registration": "affine", "control_points": points}
        )
        plan = result["data"]["report"]["plans"][0]
        assert plan["registration"]["method"] == "affine"
        terms = {item["name"]: item["value_deg"] for item in plan["uncertainty"]["components"]}
        assert terms["registration residual"] > 0.0

    def test_a_similarity_registration_rotates_the_image_azimuth(self) -> None:
        result = REGISTRY.call(
            "fib.lamella_plan",
            {"registration": "similarity", "image_rotation_deg": 20.0, "grain_id": 3},
        )
        plan = result["data"]["report"]["plans"][0]
        assert plan["theta_image_deg"] == pytest.approx(
            (plan["geometry"]["theta_sample_deg"] + 20.0) % 360.0
        )


class TestRefusals:
    def test_too_few_control_points_name_the_field(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            REGISTRY.call(
                "fib.lamella_plan", {"registration": "affine", "control_points": "0 0 1 1"}
            )
        assert excinfo.value.details["field"] == "control_points"

    def test_an_unplanned_grain_is_refused_with_a_hint(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            REGISTRY.call("fib.lamella_plan", {"grain_id": 999})
        assert excinfo.value.details["field"] == "grain_id"

    def test_an_unknown_phase_name_lists_the_known_ones(self) -> None:
        with pytest.raises(InvalidInputError) as excinfo:
            REGISTRY.call("fib.lamella_plan", {"phase_name": "unobtainium"})
        assert excinfo.value.details["field"] == "phase_name"
        assert "Nickel" in (excinfo.value.hint or "")

    def test_an_impossible_lamella_is_refused(self) -> None:
        with pytest.raises(InvalidInputError):
            REGISTRY.call("fib.lamella_plan", {"length_um": 1.0, "width_um": 5.0})
