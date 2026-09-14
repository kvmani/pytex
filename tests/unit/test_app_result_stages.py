"""Intermediate stages: the evidence behind a result, reported and exported.

A result that shows only its final number hides the step that went wrong. These
tests pin the stage contract itself, then check that the lattice-parameter
operation reports each step of its computation with numbers that agree with the
final answer — the stages must describe the run that produced the result, not a
separate one.
"""

from __future__ import annotations

import io
import zipfile

import pytest

from pytex.app import REGISTRY
from pytex.app.export import result_to_markdown, result_to_xlsx
from pytex.app.results import AppResult, ResultMetric, ResultStage

_LATTICE = "xrd.lattice_parameters"


def _lattice(**overrides: object) -> dict:
    request: dict[str, object] = {"phase": {"builtin": "ni_fcc"}}
    request.update(overrides)
    return REGISTRY.call(_LATTICE, request)


@pytest.fixture(scope="module")
def cohen() -> dict:
    return _lattice(specimen_displacement_mm=0.05)


class TestStageContract:
    def test_a_result_without_stages_keeps_its_wire_form(self) -> None:
        payload = AppResult(title="t", summary="s").to_json()
        assert "stages" not in payload

    def test_stages_serialize_in_order_with_metrics(self) -> None:
        result = AppResult(
            title="t",
            summary="s",
            stages=(
                ResultStage(key="one", title="1", summary="first", status="info"),
                ResultStage(
                    key="two",
                    title="2",
                    summary="second",
                    metrics=(ResultMetric("x", 1.5, "Å", "help"),),
                ),
            ),
        )
        stages = result.to_json()["stages"]
        assert [stage["key"] for stage in stages] == ["one", "two"]
        assert stages[1]["metrics"] == [{"label": "x", "value": 1.5, "units": "Å", "help": "help"}]
        assert "- 2: second" in result.describe()

    def test_duplicate_keys_and_unknown_statuses_are_refused(self) -> None:
        with pytest.raises(ValueError, match="unique"):
            AppResult(
                title="t",
                summary="s",
                stages=(
                    ResultStage(key="a", title="1", summary=""),
                    ResultStage(key="a", title="2", summary=""),
                ),
            )
        with pytest.raises(ValueError, match="status"):
            ResultStage(key="a", title="1", summary="", status="great")


class TestLatticeStages:
    def test_every_step_of_the_computation_is_reported(self, cohen: dict) -> None:
        keys = [stage["key"] for stage in cohen["stages"]]
        assert keys == [
            "scan",
            "peaks",
            "passes",
            "assignment",
            "least_squares",
            "cell",
            "cross_check",
        ]
        assert all(stage["explanation"] for stage in cohen["stages"])

    def test_the_assignment_stage_lists_the_reflections_the_cell_was_fitted_to(
        self, cohen: dict
    ) -> None:
        stage = next(item for item in cohen["stages"] if item["key"] == "assignment")
        assert len(stage["table"]["rows"]) == cohen["data"]["reflection_count"]

    def test_the_passes_end_at_the_final_assignment(self, cohen: dict) -> None:
        passes = next(item for item in cohen["stages"] if item["key"] == "passes")
        taken = [row for row in passes["table"]["rows"] if row["outcome"] == "taken"]
        assert taken[-1]["indexed_count"] == cohen["data"]["reflection_count"]
        counts = [row["indexed_count"] for row in taken]
        assert counts == sorted(counts)

    def test_the_correlation_matrix_is_symmetric_with_a_unit_diagonal(self, cohen: dict) -> None:
        stage = next(item for item in cohen["stages"] if item["key"] == "least_squares")
        rows = stage["table"]["rows"]
        keys = [column["key"] for column in stage["table"]["columns"][1:]]
        for i, row in enumerate(rows):
            assert row[keys[i]] == pytest.approx(1.0)
            for j, key in enumerate(keys):
                assert row[key] == pytest.approx(rows[j][keys[i]])
                assert -1.0 <= row[key] <= 1.0

    def test_the_cell_stage_reports_the_answer_in_the_summary(self, cohen: dict) -> None:
        stage = next(item for item in cohen["stages"] if item["key"] == "cell")
        metric = next(item for item in stage["metrics"] if item["label"] == "a")
        assert metric["value"] == pytest.approx(cohen["data"]["a"])

    def test_the_cross_check_reference_row_is_the_reported_value(self, cohen: dict) -> None:
        stage = next(item for item in cohen["stages"] if item["key"] == "cross_check")
        reported = stage["table"]["rows"][0]
        assert reported["a_angstrom"] == pytest.approx(cohen["data"]["a"])
        assert reported["difference_ppm"] == 0.0

    def test_le_bail_reports_its_whole_pattern_stage_and_no_peak_list(self) -> None:
        keys = [stage["key"] for stage in _lattice(method="le_bail")["stages"]]
        assert keys == ["scan", "whole_pattern", "cell"]

    def test_a_hexagonal_cell_reports_its_two_parameters_in_the_cell_stage(self) -> None:
        result = _lattice(phase={"builtin": "ti_hcp"})
        stage = next(item for item in result["stages"] if item["key"] == "cell")
        labels = {metric["label"] for metric in stage["metrics"]}
        assert {"a", "c", "c/a"} <= labels

    def test_stages_reach_the_markdown_and_workbook_exports(self, cohen: dict) -> None:
        markdown = result_to_markdown(cohen).decode("utf-8")
        assert "## How the result was reached" in markdown
        assert "### 4. Reflection assignment" in markdown
        with zipfile.ZipFile(io.BytesIO(result_to_xlsx(cohen))) as archive:
            workbook = archive.read("xl/workbook.xml").decode("utf-8")
        assert 'name="Stages"' in workbook
        assert 'name="Stage 4 assignment"' in workbook


_IDENTIFY = "xrd.phase_identification"


@pytest.fixture(scope="module")
def identification() -> dict:
    candidates = {
        "phases": [
            {"phase": {"builtin": identifier}} for identifier in ("ni_fcc", "cu_fcc", "fe_bcc")
        ]
    }
    return REGISTRY.call(_IDENTIFY, {"phase": {"builtin": "ni_fcc"}, "candidates": candidates})


class TestIdentificationStages:
    def test_every_step_of_the_identification_is_reported(self, identification: dict) -> None:
        keys = [stage["key"] for stage in identification["stages"]]
        assert keys == [
            "scan",
            "peaks",
            "cell_search",
            "best_match",
            "runner_up_match",
            "scoring",
            "decision",
        ]

    def test_the_weighted_contributions_sum_to_each_candidate_score(
        self, identification: dict
    ) -> None:
        stage = next(item for item in identification["stages"] if item["key"] == "scoring")
        for row in stage["table"]["rows"]:
            parts = [value for key, value in row.items() if key.startswith("from_") and value]
            if parts:
                assert sum(parts) == pytest.approx(row["score"])

    def test_the_peak_stage_matches_the_peaks_carried_for_the_plot(
        self, identification: dict
    ) -> None:
        stage = next(item for item in identification["stages"] if item["key"] == "peaks")
        assert len(stage["table"]["rows"]) == len(identification["data"]["peaks"])

    def test_the_decision_stage_agrees_with_the_verdict(self, identification: dict) -> None:
        stage = next(item for item in identification["stages"] if item["key"] == "decision")
        metrics = {metric["label"]: metric["value"] for metric in stage["metrics"]}
        assert metrics["Conclusive"] is identification["data"]["is_conclusive"]
        assert metrics["Decisive"] is identification["data"]["is_decisive"]
        assert stage["status"] == "ok"

    def test_the_leading_assignment_is_the_best_candidate(self, identification: dict) -> None:
        stage = next(item for item in identification["stages"] if item["key"] == "best_match")
        assert stage["summary"].startswith(identification["data"]["best_phase_name"])
