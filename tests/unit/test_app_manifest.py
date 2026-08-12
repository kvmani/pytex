"""The manifest is the user interface, so it is tested like one.

Every claim the frontend relies on is checked here: that operations are
reachable, that every operation and every parameter carries help a user can
read, that declared defaults are values the operation itself accepts, and that
unknown or malformed input is rejected with a message rather than absorbed.
"""

from __future__ import annotations

import json

import pytest

from pytex.app import REGISTRY
from pytex.app.contracts import dumps, execute, to_jsonable
from pytex.app.errors import InvalidInputError, ServiceError
from pytex.app.registry import (
    ChoiceParameter,
    IndicesListParameter,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    ServiceRegistry,
)
from pytex.app.results import APP_RESULT_SCHEMA


def test_manifest_lists_every_registered_operation() -> None:
    manifest = REGISTRY.manifest()
    assert manifest["schema"] == "pytex.app_manifest/1"
    listed = {entry["id"] for entry in manifest["operations"]}
    assert listed == set(REGISTRY.ids())
    assert listed, "the application must expose at least one operation"


@pytest.mark.parametrize("spec", REGISTRY.operations(), ids=lambda spec: spec.id)
def test_every_operation_is_documented(spec) -> None:  # type: ignore[no-untyped-def]
    assert spec.title.strip(), f"{spec.id} has no title"
    assert len(spec.summary.strip()) > 10, f"{spec.id} has no usable summary"
    assert len(spec.help_text.strip()) > 40, f"{spec.id} has no usable help text"
    assert spec.returns.strip(), f"{spec.id} does not say what it returns"
    assert spec.panel.strip(), f"{spec.id} is not assigned to a panel"


@pytest.mark.parametrize("spec", REGISTRY.operations(), ids=lambda spec: spec.id)
def test_every_parameter_is_documented(spec) -> None:  # type: ignore[no-untyped-def]
    for parameter in spec.parameters:
        assert parameter.label.strip(), f"{spec.id}.{parameter.name} has no label"
        assert len(parameter.help_text.strip()) > 10, (
            f"{spec.id}.{parameter.name} has no usable help text"
        )


@pytest.mark.parametrize("spec", REGISTRY.operations(), ids=lambda spec: spec.id)
def test_declared_defaults_are_accepted_by_their_own_validator(spec) -> None:  # type: ignore[no-untyped-def]
    for parameter in spec.parameters:
        if parameter.default is not None:
            parameter.coerce(parameter.default)


@pytest.mark.parametrize("spec", REGISTRY.operations(), ids=lambda spec: spec.id)
def test_manifest_entry_is_json_serialisable(spec) -> None:  # type: ignore[no-untyped-def]
    json.loads(dumps(spec.describe()))


def test_unknown_operation_is_reported_with_the_known_list() -> None:
    envelope, status = execute("calc.does_not_exist", {})
    assert status == 404
    assert envelope["ok"] is False
    assert envelope["error"]["code"] == "operation.unknown"
    assert "calc.plane_angles" in envelope["error"]["details"]["known"]


def test_unknown_parameter_is_rejected_rather_than_ignored() -> None:
    envelope, status = execute(
        "calc.plane_angles", {"phase": {"builtin": "ni_fcc"}, "planez": [[1, 1, 1]]}
    )
    assert status == 400
    assert envelope["error"]["code"] == "input.invalid"
    assert "planez" in envelope["error"]["message"]


def test_successful_call_returns_a_result_envelope() -> None:
    envelope, status = execute("calc.catalog", {})
    assert status == 200
    assert envelope["ok"] is True
    assert envelope["result"]["schema"] == APP_RESULT_SCHEMA


def test_registry_rejects_duplicate_identifiers() -> None:
    registry = ServiceRegistry()

    @registry.operation(
        "demo.echo",
        title="Echo",
        summary="Return the input unchanged.",
        help_text="Exists only to exercise the registry in tests, and echoes its input.",
    )
    def _echo(request: dict[str, object]) -> dict[str, object]:
        return dict(request)

    with pytest.raises(ValueError, match="already registered"):
        registry.register(registry.get("demo.echo"))


class TestParameterValidation:
    """Each parameter kind rejects what it should, with the field named."""

    def test_number_respects_bounds(self) -> None:
        parameter = NumberParameter("x", label="X", help_text="A bounded number.", minimum=1.0)
        assert parameter.coerce("2.5") == 2.5
        with pytest.raises(InvalidInputError) as excinfo:
            parameter.coerce(0.5)
        assert excinfo.value.details["field"] == "x"

    def test_number_rejects_non_finite(self) -> None:
        parameter = NumberParameter("x", label="X", help_text="A finite number.")
        with pytest.raises(InvalidInputError):
            parameter.coerce(float("inf"))

    def test_integer_rejects_fractional_values(self) -> None:
        parameter = IntegerParameter("n", label="N", help_text="A whole number.")
        assert parameter.coerce(3.0) == 3
        with pytest.raises(InvalidInputError):
            parameter.coerce(3.5)

    def test_choice_rejects_unlisted_values(self) -> None:
        parameter = ChoiceParameter(
            "kind",
            label="Kind",
            help_text="One of two things.",
            options=(("a", "A", "The first."), ("b", "B", "The second.")),
        )
        assert parameter.coerce("b") == "b"
        with pytest.raises(InvalidInputError):
            parameter.coerce("c")

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("1 1 1", (1, 1, 1)),
            ("1,1,0", (1, 1, 0)),
            ("(1 -1 0)", (1, -1, 0)),
            ("[112]", (1, 1, 2)),
            ("111", (1, 1, 1)),
        ],
    )
    def test_indices_accept_the_ways_people_type_them(
        self, text: str, expected: tuple[int, ...]
    ) -> None:
        parameter = IndicesParameter("hkl", label="Plane", help_text="Three Miller indices.")
        assert parameter.coerce(text) == expected

    def test_indices_reject_the_zero_row(self) -> None:
        parameter = IndicesParameter("hkl", label="Plane", help_text="Three Miller indices.")
        with pytest.raises(InvalidInputError, match="all-zero"):
            parameter.coerce("0 0 0")

    def test_indices_report_the_wrong_count_with_a_hint(self) -> None:
        parameter = IndicesParameter("hkl", label="Plane", help_text="Three Miller indices.")
        with pytest.raises(InvalidInputError) as excinfo:
            parameter.coerce("1 1")
        assert excinfo.value.hint is not None

    def test_indices_list_accepts_newline_separated_text(self) -> None:
        parameter = IndicesListParameter("planes", label="Planes", help_text="Any number of rows.")
        assert parameter.coerce("1 1 1\n1 -1 0") == ((1, 1, 1), (1, -1, 0))

    def test_indices_list_accepts_the_empty_input(self) -> None:
        parameter = IndicesListParameter(
            "planes", label="Planes", help_text="Any number of rows.", required=False
        )
        assert parameter.coerce("") == ()


class TestJsonConversion:
    """NumPy leaves the process as JSON, and non-finite floats leave as null."""

    def test_numpy_scalars_and_arrays_convert(self) -> None:
        import numpy as np

        payload = to_jsonable({"a": np.float64(1.5), "b": np.arange(3)})
        assert payload == {"a": 1.5, "b": [0, 1, 2]}

    def test_non_finite_floats_become_null(self) -> None:
        assert to_jsonable(float("nan")) is None
        assert json.loads(dumps({"x": float("inf")})) == {"x": None}

    def test_objects_with_to_json_are_used(self) -> None:
        class Thing:
            def to_json(self) -> dict[str, int]:
                return {"n": 1}

        assert to_jsonable(Thing()) == {"n": 1}


def test_service_errors_carry_a_message_and_a_hint() -> None:
    error = ServiceError("Something specific went wrong.", hint="Try the other thing.")
    payload = error.to_json()
    assert payload["message"].endswith(".")
    assert payload["hint"] == "Try the other thing."


class TestCanonicalExamples:
    """Every panel ships runnable examples, and they are run here.

    An example that no longer executes is worse than no example: it is the
    first thing a new user clicks. These tests execute each one against the
    real operation, so an example cannot drift out of step with the service it
    demonstrates.
    """

    def test_every_panel_has_at_least_three_examples(self) -> None:
        for panel in REGISTRY.panels():
            examples = REGISTRY.examples(panel=panel)
            assert len(examples) >= 3, f"panel {panel!r} has only {len(examples)} example(s)"

    @pytest.mark.parametrize("example", REGISTRY.examples(), ids=lambda example: example.id)
    def test_every_example_runs(self, example) -> None:  # type: ignore[no-untyped-def]
        result = REGISTRY.call(example.operation, example.request)
        assert result["summary"]

    @pytest.mark.parametrize("example", REGISTRY.examples(), ids=lambda example: example.id)
    def test_every_example_explains_what_it_teaches(self, example) -> None:  # type: ignore[no-untyped-def]
        assert example.title.strip()
        assert len(example.summary.strip()) > 10
        assert len(example.teaches.strip()) > 60, (
            f"{example.id} does not say what to notice once it has run"
        )
        assert example.operation in REGISTRY.ids()

    def test_the_canonical_materials_are_all_reachable(self) -> None:
        from pytex.app.phases import BUILTIN_PHASES

        for identifier in ("nacl", "austenite_fcc", "fe_bcc", "zr_hcp"):
            assert identifier in BUILTIN_PHASES

    def test_examples_appear_in_the_manifest(self) -> None:
        manifest = REGISTRY.manifest()
        listed = {entry["id"] for entry in manifest["examples"]}
        assert listed == {example.id for example in REGISTRY.examples()}
        assert manifest["panels"] == list(REGISTRY.panels())
