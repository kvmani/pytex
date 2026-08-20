"""Self-describing operations: the single source of truth for the user interface.

An *operation* is one thing a user can ask the application to do. It is declared
once, here, together with everything anyone needs to know about it: what it is
called, what it does, what it takes, what it gives back, and where the science
comes from. The frontend does not contain a list of operations, a form layout, or
a help string; it fetches :func:`ServiceRegistry.manifest` and builds itself.

Why this is worth the ceremony
------------------------------
Hand-written forms and hand-written help drift away from the code they describe,
and the drift is silent — a stale tooltip misleads a researcher with no error
anywhere. Declaring the parameter *once*, next to the function that consumes it,
makes the label, the units, the default, the validation, and the help text the
same object. `tests/unit/test_app_manifest.py` then enforces that every operation
and every parameter is actually documented, so an undocumented surface fails the
suite rather than shipping.

Example
-------
>>> registry = ServiceRegistry()
>>> @registry.operation(
...     "demo.double",
...     title="Double a number",
...     summary="Return twice the input.",
...     help_text="Exists to demonstrate the registry; multiplies by two.",
...     parameters=(NumberParameter("value", label="Value", help_text="Any real number."),),
...     returns="An object with the doubled value.",
... )
... def _double(request: dict[str, object]) -> dict[str, object]:
...     return {"value": float(request["value"]) * 2.0}  # type: ignore[arg-type]
>>> registry.call("demo.double", {"value": 3})
{'value': 6.0}
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from pytex.app.errors import InvalidInputError, UnknownOperationError

__all__ = [
    "REGISTRY",
    "BooleanParameter",
    "ChoiceParameter",
    "ExampleScenario",
    "IndicesListParameter",
    "IndicesParameter",
    "IntegerParameter",
    "NumberParameter",
    "ObjectParameter",
    "OperationSpec",
    "Parameter",
    "ServiceRegistry",
    "TextParameter",
]

#: A handler takes a validated parameter mapping and returns a JSON-ready result.
Handler = Callable[[dict[str, Any]], dict[str, Any]]

_SPHINX_SOURCE_BASE = "https://github.com/kvmani/pytex/blob/main/docs/site/"


@dataclass(frozen=True)
class DocumentationLink:
    """A canonical Sphinx page attached to operation help.

    ``path`` is the document name below ``docs/site`` without ``.md``. Keeping
    the source-relative name in the manifest lets tests prove the page exists;
    ``url`` points at GitHub's rendered view of that canonical MyST source until
    the project publishes a dedicated Sphinx host.
    """

    title: str
    path: str

    def describe(self) -> dict[str, str]:
        """Return the browser-facing documentation target."""

        return {
            "title": self.title,
            "path": self.path,
            "url": f"{_SPHINX_SOURCE_BASE}{self.path}.md",
        }


_PANEL_DOCUMENTATION = {
    "calculator": DocumentationLink(
        "Miller planes and directions", "concepts/miller_planes_directions"
    ),
    "crystal": DocumentationLink(
        "Crystal visualization workflow", "workflows/crystal_visualization"
    ),
    "diffraction": DocumentationLink(
        "Composite orientation-relationship diffraction", "workflows/composite_or_diffraction"
    ),
    "xrd": DocumentationLink("Powder XRD generation", "workflows/xrd_generation"),
    "tem": DocumentationLink("TEM specimen tilt navigation", "theory/tem_specimen_tilt_navigation"),
    "tem_simulator": DocumentationLink("SAED pattern generation", "workflows/saed_generation"),
    "cbed": DocumentationLink(
        "Convergent-beam electron diffraction", "theory/convergent_beam_electron_diffraction"
    ),
    "ebsd": DocumentationLink(
        "EBSD grain segmentation and GROD", "theory/ebsd_grain_segmentation_and_grod"
    ),
    "texture": DocumentationLink("Orientation and texture", "concepts/orientation_texture"),
    "variants": DocumentationLink(
        "Orientation relationships", "concepts/orientation_relationships"
    ),
}


@dataclass(frozen=True)
class Parameter:
    """One input of an operation, described well enough to render a control.

    Purpose
    -------
    Base class of the parameter kinds below. Subclasses add type-specific
    validation and the extra descriptors a control needs (bounds, options,
    index width). Instances are frozen and shared; they hold no per-call state.

    Attributes
    ----------
    name : str
        Key in the request mapping.
    label : str
        Short human label for the control. Sentence case, no trailing colon.
    help_text : str
        One or two sentences explaining what the value means and how to choose
        it. Shown in the inline help popover. Required — an empty help string
        fails the manifest test.
    required : bool
        When false, ``default`` is substituted for a missing value.
    default : object, optional
        Pre-filled value. Must itself pass this parameter's validation.
    units : str, optional
        Physical unit rendered as a suffix inside the control, e.g. ``"Å"``.
    group : str, optional
        Control-panel section this parameter belongs to, so a long form can be
        split without the frontend inventing a grouping.
    advanced : bool
        Hidden behind the panel's "Advanced" disclosure. Use for parameters
        with a good default that most users should not touch.
    """

    name: str
    label: str
    help_text: str
    required: bool = True
    default: Any = None
    units: str | None = None
    group: str | None = None
    advanced: bool = False

    #: Discriminator the frontend switches on to choose a control widget.
    kind: str = field(default="text", init=False)

    def describe(self) -> dict[str, Any]:
        """Return the manifest entry for this parameter."""

        payload: dict[str, Any] = {
            "name": self.name,
            "kind": self.kind,
            "label": self.label,
            "help": self.help_text,
            "required": self.required,
            "advanced": self.advanced,
        }
        if self.default is not None:
            payload["default"] = self.default
        if self.units is not None:
            payload["units"] = self.units
        if self.group is not None:
            payload["group"] = self.group
        payload.update(self._describe_extra())
        return payload

    def _describe_extra(self) -> dict[str, Any]:
        return {}

    def coerce(self, value: Any) -> Any:
        """Validate and convert one supplied value.

        Raises
        ------
        InvalidInputError
            If the value cannot be interpreted as this parameter's type or
            violates a declared constraint.
        """

        return value

    def _invalid(self, message: str, *, hint: str | None = None) -> InvalidInputError:
        return InvalidInputError(
            f"{self.label}: {message}",
            field=self.name,
            hint=hint or self.help_text,
        )


@dataclass(frozen=True)
class NumberParameter(Parameter):
    """A real-valued input, optionally bounded."""

    minimum: float | None = None
    maximum: float | None = None
    step: float | None = None
    exclusive_minimum: bool = False

    kind: str = field(default="number", init=False)

    def _describe_extra(self) -> dict[str, Any]:
        extra: dict[str, Any] = {}
        if self.minimum is not None:
            extra["minimum"] = self.minimum
            extra["exclusiveMinimum"] = self.exclusive_minimum
        if self.maximum is not None:
            extra["maximum"] = self.maximum
        if self.step is not None:
            extra["step"] = self.step
        return extra

    def coerce(self, value: Any) -> float:
        if isinstance(value, bool) or not isinstance(value, int | float | str):
            raise self._invalid(f"expected a number, got {type(value).__name__}.")
        try:
            number = float(value)
        except ValueError as exc:
            raise self._invalid(f"{value!r} is not a number.") from exc
        if not math.isfinite(number):
            raise self._invalid("must be finite.")
        if self.minimum is not None:
            too_small = number <= self.minimum if self.exclusive_minimum else number < self.minimum
            if too_small:
                relation = "greater than" if self.exclusive_minimum else "at least"
                raise self._invalid(f"must be {relation} {self.minimum:g}; got {number:g}.")
        if self.maximum is not None and number > self.maximum:
            raise self._invalid(f"must be at most {self.maximum:g}; got {number:g}.")
        return number


@dataclass(frozen=True)
class IntegerParameter(Parameter):
    """An integer-valued input, optionally bounded."""

    minimum: int | None = None
    maximum: int | None = None

    kind: str = field(default="integer", init=False)

    def _describe_extra(self) -> dict[str, Any]:
        extra: dict[str, Any] = {}
        if self.minimum is not None:
            extra["minimum"] = self.minimum
        if self.maximum is not None:
            extra["maximum"] = self.maximum
        return extra

    def coerce(self, value: Any) -> int:
        if isinstance(value, bool):
            raise self._invalid("expected an integer, got a boolean.")
        if isinstance(value, float) and not value.is_integer():
            raise self._invalid(f"must be a whole number; got {value!r}.")
        try:
            number = int(value)
        except (TypeError, ValueError) as exc:
            raise self._invalid(f"{value!r} is not an integer.") from exc
        if self.minimum is not None and number < self.minimum:
            raise self._invalid(f"must be at least {self.minimum}; got {number}.")
        if self.maximum is not None and number > self.maximum:
            raise self._invalid(f"must be at most {self.maximum}; got {number}.")
        return number


@dataclass(frozen=True)
class TextParameter(Parameter):
    """A free-text input."""

    max_length: int | None = None
    multiline: bool = False
    placeholder: str | None = None

    kind: str = field(default="text", init=False)

    def _describe_extra(self) -> dict[str, Any]:
        extra: dict[str, Any] = {"multiline": self.multiline}
        if self.max_length is not None:
            extra["maxLength"] = self.max_length
        if self.placeholder is not None:
            extra["placeholder"] = self.placeholder
        return extra

    def coerce(self, value: Any) -> str:
        if not isinstance(value, str):
            raise self._invalid(f"expected text, got {type(value).__name__}.")
        text = value.strip()
        if self.required and not text:
            raise self._invalid("must not be empty.")
        if self.max_length is not None and len(text) > self.max_length:
            raise self._invalid(f"must be at most {self.max_length} characters.")
        return text


@dataclass(frozen=True)
class BooleanParameter(Parameter):
    """A toggle."""

    kind: str = field(default="boolean", init=False)

    def coerce(self, value: Any) -> bool:
        if isinstance(value, bool):
            return value
        if isinstance(value, str) and value.lower() in {"true", "false"}:
            return value.lower() == "true"
        raise self._invalid(f"expected true or false, got {value!r}.")


@dataclass(frozen=True)
class ChoiceParameter(Parameter):
    """A selection from a fixed set of options.

    ``options`` carries a label and a help string per choice, so a select
    control can explain each entry rather than showing bare identifiers.
    """

    options: tuple[tuple[str, str, str], ...] = ()

    kind: str = field(default="choice", init=False)

    def _describe_extra(self) -> dict[str, Any]:
        return {
            "options": [
                {"value": value, "label": label, "help": help_text}
                for value, label, help_text in self.options
            ]
        }

    @property
    def values(self) -> tuple[str, ...]:
        """The accepted option identifiers."""

        return tuple(value for value, _, _ in self.options)

    def coerce(self, value: Any) -> str:
        if not isinstance(value, str) or value not in self.values:
            allowed = ", ".join(self.values)
            raise self._invalid(f"must be one of {allowed}; got {value!r}.")
        return value


@dataclass(frozen=True)
class IndicesParameter(Parameter):
    """One set of Miller indices, entered as text or as a list.

    Accepts ``"1 1 1"``, ``"1,1,1"``, ``"(1 -1 0)"``, ``"1 -1 0"`` or
    ``[1, 1, 1]``. Overbars are not accepted on input — a minus sign is what a
    keyboard has — but output is formatted through :mod:`pytex.core.notation`,
    which does use overbars.
    """

    width: int = 3
    allow_zero: bool = False

    kind: str = field(default="indices", init=False)

    def _describe_extra(self) -> dict[str, Any]:
        return {"width": self.width, "allowZero": self.allow_zero}

    def coerce(self, value: Any) -> tuple[int, ...]:
        return _coerce_index_row(value, parameter=self)


@dataclass(frozen=True)
class IndicesListParameter(Parameter):
    """Any number of index rows, one per line or as a list of lists.

    This is the parameter behind "superimpose an arbitrary number of planes":
    the control is a repeatable row editor, and the empty list is legal.
    """

    width: int = 3
    allow_zero: bool = False
    max_rows: int | None = None

    kind: str = field(default="indices-list", init=False)

    def _describe_extra(self) -> dict[str, Any]:
        extra: dict[str, Any] = {"width": self.width, "allowZero": self.allow_zero}
        if self.max_rows is not None:
            extra["maxRows"] = self.max_rows
        return extra

    def coerce(self, value: Any) -> tuple[tuple[int, ...], ...]:
        rows: list[Any]
        if value is None:
            rows = []
        elif isinstance(value, str):
            rows = [line for line in value.replace(";", "\n").splitlines() if line.strip()]
        elif isinstance(value, Sequence):
            rows = list(value)
        else:
            raise self._invalid(f"expected a list of index rows, got {type(value).__name__}.")
        if self.max_rows is not None and len(rows) > self.max_rows:
            raise self._invalid(f"accepts at most {self.max_rows} rows; got {len(rows)}.")
        return tuple(_coerce_index_row(row, parameter=self) for row in rows)


@dataclass(frozen=True)
class ObjectParameter(Parameter):
    """A nested JSON object validated by the service itself.

    Used where the payload is a structured document — a phase specification, a
    picked-spot list — that has its own constructor and its own error messages.
    The frontend renders these with a bespoke control named by ``editor``.
    """

    editor: str = "json"

    kind: str = field(default="object", init=False)

    def _describe_extra(self) -> dict[str, Any]:
        return {"editor": self.editor}

    def coerce(self, value: Any) -> dict[str, Any]:
        if not isinstance(value, Mapping):
            raise self._invalid(f"expected an object, got {type(value).__name__}.")
        return dict(value)


def _coerce_index_row(value: Any, *, parameter: Parameter) -> tuple[int, ...]:
    """Parse one row of Miller indices from text or a sequence."""

    width = getattr(parameter, "width", 3)
    allow_zero = getattr(parameter, "allow_zero", False)
    items: list[Any]
    if isinstance(value, str):
        cleaned = value.strip().strip("()[]{}<>")
        cleaned = cleaned.replace(",", " ")
        items = cleaned.split()
        if len(items) == 1 and len(items[0]) == width and _is_compact_row(items[0]):
            items = list(items[0])
    elif isinstance(value, Sequence) and not isinstance(value, bytes | bytearray):
        items = list(value)
    else:
        raise parameter._invalid(f"expected {width} indices, got {type(value).__name__}.")
    if len(items) != width:
        raise parameter._invalid(
            f"expected {width} indices; got {len(items)} ({value!r}).",
            hint=f"Enter {width} whole numbers separated by spaces, e.g. "
            + " ".join(["1"] * width)
            + ".",
        )
    indices: list[int] = []
    for item in items:
        if isinstance(item, bool):
            raise parameter._invalid("indices must be whole numbers.")
        if isinstance(item, float) and not item.is_integer():
            raise parameter._invalid(f"indices must be whole numbers; got {item!r}.")
        try:
            indices.append(int(item))
        except (TypeError, ValueError) as exc:
            raise parameter._invalid(f"{item!r} is not a whole number.") from exc
    if not allow_zero and not any(indices):
        raise parameter._invalid("the all-zero index row does not describe a plane or direction.")
    return tuple(indices)


def _is_compact_row(text: str) -> bool:
    """True for compact forms such as ``"111"`` that carry no sign."""

    return text.isdigit()


@dataclass(frozen=True)
class OperationSpec:
    """Everything the application knows about one operation.

    Purpose
    -------
    The manifest entry the frontend renders and the contract the server
    validates against, in one object beside the handler it describes.

    Attributes
    ----------
    id : str
        Dotted, stable, and the key the client calls: ``"calc.plane_angles"``.
    title : str
        Menu- and heading-level name.
    summary : str
        One line, shown in the command palette and as the panel subtitle.
    help_text : str
        Long form, Markdown, shown in the help drawer. Explain *when to use it*
        and how to read the result, not only what it computes.
    parameters : tuple of Parameter
    returns : str
        Prose description of the response object.
    panel : str
        Tab this operation belongs to, e.g. ``"calculator"``. Drives grouping in
        the capability index.
    citations : tuple of str
        Normative sources for the science, per `AGENTS.md`.
    tags : tuple of str
        Extra search keywords for feature discovery.
    documentation : DocumentationLink or None
        Closest canonical Sphinx page for the inline help drawer.
    """

    id: str
    title: str
    summary: str
    help_text: str
    handler: Handler
    parameters: tuple[Parameter, ...] = ()
    returns: str = ""
    panel: str = "general"
    citations: tuple[str, ...] = ()
    tags: tuple[str, ...] = ()
    documentation: DocumentationLink | None = None

    def describe(self) -> dict[str, Any]:
        """Return the manifest entry for this operation."""

        return {
            "id": self.id,
            "title": self.title,
            "summary": self.summary,
            "help": self.help_text,
            "panel": self.panel,
            "returns": self.returns,
            "citations": list(self.citations),
            "tags": list(self.tags),
            "documentation": self.documentation.describe() if self.documentation else None,
            "parameters": [parameter.describe() for parameter in self.parameters],
        }

    def bind(self, request: Mapping[str, Any]) -> dict[str, Any]:
        """Validate a raw request against this operation's parameters.

        Unknown keys are rejected rather than ignored: a typo in a parameter
        name is a mistake worth reporting, and silently dropping it produces a
        result computed from a default the user did not intend.
        """

        declared = {parameter.name for parameter in self.parameters}
        unexpected = sorted(set(request) - declared)
        if unexpected:
            raise InvalidInputError(
                f"{self.title}: unknown parameter(s) {', '.join(unexpected)}.",
                hint=f"Accepted parameters: {', '.join(sorted(declared)) or 'none'}.",
                details={"unexpected": unexpected, "accepted": sorted(declared)},
            )
        bound: dict[str, Any] = {}
        for parameter in self.parameters:
            if parameter.name in request and request[parameter.name] is not None:
                bound[parameter.name] = parameter.coerce(request[parameter.name])
            elif parameter.required and parameter.default is None:
                raise InvalidInputError(
                    f"{self.title}: {parameter.label} is required.",
                    field=parameter.name,
                    hint=parameter.help_text,
                )
            else:
                bound[parameter.name] = (
                    parameter.coerce(parameter.default) if parameter.default is not None else None
                )
        return bound


@dataclass(frozen=True)
class ExampleScenario:
    """A canonical, runnable demonstration attached to a panel.

    Purpose
    -------
    Someone opening a panel for the first time may have no data of their own,
    and an empty form teaches nothing. An example is a complete set of inputs
    that runs the *real* operation and leaves the user in a working state they
    can then modify — not a screenshot, and not a separate demo code path that
    can rot while the operation moves on.

    Every panel carries at least three, built on the shared canonical materials
    (NaCl, Fe-fcc, Fe-bcc, Zr-hcp) so familiarity accumulates across tabs
    instead of restarting in each one.

    Attributes
    ----------
    id : str
        Stable identifier, dotted like an operation.
    title : str
        Name shown on the "Try an example" menu.
    panel : str
        Tab this example belongs to.
    summary : str
        One line: what this example shows.
    teaches : str
        The thing to notice once it has run. This is the sentence that turns a
        demonstration into a lesson, so it states a specific observation rather
        than repeating the title.
    operation : str
        Operation to run.
    request : mapping
        Parameters to run it with. Validated by the operation like any request,
        and exercised by the test suite, so an example cannot drift out of step
        with the operation it demonstrates.
    """

    id: str
    title: str
    panel: str
    summary: str
    teaches: str
    operation: str
    request: Mapping[str, Any] = field(default_factory=dict)

    def describe(self) -> dict[str, Any]:
        """Return the manifest entry for this example."""

        return {
            "id": self.id,
            "title": self.title,
            "panel": self.panel,
            "summary": self.summary,
            "teaches": self.teaches,
            "operation": self.operation,
            "request": dict(self.request),
        }


class ServiceRegistry:
    """The set of operations an application instance exposes.

    Purpose
    -------
    Holds operations and their examples, produces the manifest, and dispatches
    calls. Deliberately the only lookup path: there is no way to invoke a
    service that is not in the manifest, so "documented" and "reachable" cannot
    diverge.
    """

    def __init__(self) -> None:
        self._operations: dict[str, OperationSpec] = {}
        self._examples: dict[str, ExampleScenario] = {}

    def operation(
        self,
        operation_id: str,
        *,
        title: str,
        summary: str,
        help_text: str,
        parameters: Iterable[Parameter] = (),
        returns: str = "",
        panel: str = "general",
        citations: Iterable[str] = (),
        tags: Iterable[str] = (),
        documentation: DocumentationLink | None = None,
    ) -> Callable[[Handler], Handler]:
        """Register a handler, returning it unchanged.

        Used as a decorator directly above the handler, so the declaration and
        the implementation are read together.
        """

        def decorate(handler: Handler) -> Handler:
            spec = OperationSpec(
                id=operation_id,
                title=title,
                summary=summary,
                help_text=help_text,
                handler=handler,
                parameters=tuple(parameters),
                returns=returns,
                panel=panel,
                citations=tuple(citations),
                tags=tuple(tags),
                documentation=documentation or _PANEL_DOCUMENTATION.get(panel),
            )
            self.register(spec)
            return handler

        return decorate

    def register(self, spec: OperationSpec) -> None:
        """Add an operation, rejecting duplicate identifiers."""

        if spec.id in self._operations:
            raise ValueError(f"Operation {spec.id!r} is already registered.")
        self._operations[spec.id] = spec

    def get(self, operation_id: str) -> OperationSpec:
        """Return one operation or raise :class:`UnknownOperationError`."""

        try:
            return self._operations[operation_id]
        except KeyError:
            raise UnknownOperationError(operation_id, known=self.ids()) from None

    def ids(self) -> tuple[str, ...]:
        """Registered operation identifiers, sorted."""

        return tuple(sorted(self._operations))

    def operations(self) -> tuple[OperationSpec, ...]:
        """Registered operations, sorted by identifier."""

        return tuple(self._operations[key] for key in self.ids())

    def add_example(self, example: ExampleScenario) -> None:
        """Register a runnable example, rejecting duplicate identifiers."""

        if example.id in self._examples:
            raise ValueError(f"Example {example.id!r} is already registered.")
        self._examples[example.id] = example

    def add_examples(self, examples: Iterable[ExampleScenario]) -> None:
        """Register several examples in declaration order."""

        for example in examples:
            self.add_example(example)

    def examples(self, *, panel: str | None = None) -> tuple[ExampleScenario, ...]:
        """Registered examples in declaration order, optionally filtered to a panel.

        Declaration order, not sorted order: the first example a panel offers is
        the one a newcomer runs, so which one comes first is a teaching decision
        that belongs to whoever wrote them.
        """

        found = tuple(self._examples.values())
        if panel is None:
            return found
        return tuple(example for example in found if example.panel == panel)

    def panels(self) -> tuple[str, ...]:
        """Panels that have at least one operation."""

        return tuple(sorted({spec.panel for spec in self.operations()}))

    def manifest(self) -> dict[str, Any]:
        """Return the document the frontend builds itself from."""

        from pytex import __version__
        from pytex.app.about import about_document
        from pytex.app.export import EXPORT_FORMATS

        return {
            "schema": "pytex.app_manifest/1",
            "version": __version__,
            # Identity, authorship and licence travel with the manifest the
            # frontend already fetches at startup, so the About panel needs no
            # request of its own and cannot show a version the running code was
            # not built from.
            "about": about_document(),
            "panels": list(self.panels()),
            "operations": [spec.describe() for spec in self.operations()],
            "examples": [example.describe() for example in self.examples()],
            # Published rather than hard-coded in the browser: a format added in
            # Python then appears on every result in every panel at once, which
            # is the same rule the operations themselves follow.
            "export_formats": [
                {
                    "id": identifier,
                    "label": spec["label"],
                    "description": spec["description"],
                    "extension": spec["extension"],
                }
                for identifier, spec in EXPORT_FORMATS.items()
            ],
        }

    def call(self, operation_id: str, request: Mapping[str, Any] | None = None) -> dict[str, Any]:
        """Validate and execute one operation."""

        spec = self.get(operation_id)
        return spec.handler(spec.bind(request or {}))


#: The application-wide registry. Importing :mod:`pytex.app.services` fills it.
REGISTRY = ServiceRegistry()
