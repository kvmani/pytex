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
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from pytex.app.errors import InvalidInputError, UnknownOperationError
from pytex.core.symbols import has_symbol, symbol

__all__ = [
    "FIELD_WIDTHS",
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

#: How wide a control needs to be, named by the content rather than by pixels.
#:
#: A parameter knows what it holds; it does not know how wide the rail is, what
#: the font is, or whether the window is a laptop or a phone. So it declares the
#: former and the stylesheet owns the latter. The tokens, narrowest first:
#:
#: ``index``
#:     A single Miller index or a similarly tiny integer. Three characters,
#:     because that is the widest an index ever is in practice and a box sized
#:     for more invites a run of digits (Decision 9 of the application platform
#:     document).
#: ``tiny``
#:     A small bounded count: a repeat number, a harmonic degree, an index
#:     limit. Up to four characters.
#: ``short``
#:     A physical quantity with a known scale and a decimal or two: an angle in
#:     degrees, a lattice edge in angstrom, a tilt. Up to six characters.
#: ``medium``
#:     The default, and the width for anything whose magnitude is not bounded
#:     by its meaning.
#: ``long``
#:     A value that is genuinely text-like and needs room to be read: a file
#:     name, a formula, a comma-separated list.
#: ``full``
#:     The whole rail. For controls that are not a single value at all --- a
#:     text area, a JSON editor, a repeatable row stack.
FIELD_WIDTHS: tuple[str, ...] = ("index", "tiny", "short", "medium", "long", "full")


@dataclass(frozen=True)
class DocumentationLink:
    """A canonical Sphinx page attached to operation help.

    ``path`` is the document name below ``docs/site`` without ``.md``. Keeping
    the source-relative name in the manifest lets tests prove the page exists;
    ``url`` points at the local in-app HTML document served directly by the PyTex
    workbench, ensuring publication-quality offline MathJax mathematics
    rendering on airgapped or office network environments.
    """

    title: str
    path: str

    def describe(self) -> dict[str, str]:
        """Return the browser-facing documentation target."""

        return {
            "title": self.title,
            "path": self.path,
            "url": f"/docs/{self.path}.html",
            "source_url": f"{_SPHINX_SOURCE_BASE}{self.path}.md",
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
    "ebsd_kikuchi": DocumentationLink("Kikuchi band geometry", "workflows/kikuchi_geometry"),
    "ebsd_or": DocumentationLink(
        "Orientation relationships", "concepts/orientation_relationships"
    ),
    "ecci": DocumentationLink("Kikuchi band geometry", "workflows/kikuchi_geometry"),
    "kearns": DocumentationLink(
        "The Kearns parameter and basal-pole texture",
        "theory/kearns_parameter_and_basal_pole_texture",
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
    group_collapsed : bool
        Start this parameter's group closed. For a section that matters only
        under one setting of another control -- the parallelisms of a custom
        orientation relationship, which are ignored unless the relationship is
        *Custom* -- so it costs one line of rail rather than five when it does
        not apply. Declared on the parameter rather than on the group because
        there is no group object: a group is whatever set of parameters names
        it, and one member asking for it is enough.
    advanced : bool
        Hidden behind the panel's "Advanced" disclosure. Use for parameters
        with a good default that most users should not touch.
    symbol : str, optional
        Name of a registered symbol in :mod:`pytex.core.symbols`. The control
        *shows* the symbol --- the two characters of a first Bunge angle rather
        than the six of ``phi1`` --- while ``label`` remains what the control is
        *called*, and is carried into the accessible name, the tooltip and every
        error message. This is what lets three Euler angles share one line: the
        label was most of the width, not the value. Registration is checked at
        construction, so a misspelled name fails the import rather than
        shipping a control labelled with a raw registry key.
    row : str, optional
        Name of a layout row. Parameters that declare the same ``row`` render
        side by side on one line instead of one under the other. Use it for a
        set of quantities that are read and entered *together* and are short
        enough to fit --- the three Bunge angles, a tilt pair, the cell edges.
        It is a presentation grouping only: it changes no name, no default and
        no validation, and a client that ignores it renders a correct form.
    field_width : str, optional
        How much width the control needs, as one of :data:`FIELD_WIDTHS`.
        Declared in characters-of-content terms rather than pixels, because
        that is the thing the parameter actually knows: a Miller index is never
        more than a few characters however wide the rail is. When omitted, a
        width is derived from the declared bounds --- see
        :meth:`_derived_field_width` --- so most parameters need not say
        anything and still stop short of the full rail.

    See Also
    --------
    pytex.core.symbols : The registry ``symbol`` names.
    """

    name: str
    label: str
    help_text: str
    required: bool = True
    default: Any = None
    units: str | None = None
    group: str | None = None
    group_collapsed: bool = False
    advanced: bool = False
    symbol: str | None = None
    row: str | None = None
    field_width: str | None = None

    #: Discriminator the frontend switches on to choose a control widget.
    kind: str = field(default="text", init=False)

    def __post_init__(self) -> None:
        # Construction-time, per AGENTS.md: an invariant checked where the
        # mistake is made names the parameter that made it. Checked downstream
        # it would surface as a control labelled with a raw registry key, which
        # nothing fails on.
        if self.symbol is not None and not has_symbol(self.symbol):
            try:
                symbol(self.symbol)
            except KeyError as exc:
                raise ValueError(f"Parameter {self.name!r}: {exc.args[0]}") from None
        if self.field_width is not None and self.field_width not in FIELD_WIDTHS:
            allowed = ", ".join(FIELD_WIDTHS)
            raise ValueError(
                f"Parameter {self.name!r}: field_width must be one of {allowed}; "
                f"got {self.field_width!r}."
            )

    def describe(self) -> dict[str, Any]:
        """Return the manifest entry for this parameter."""

        payload: dict[str, Any] = {
            "name": self.name,
            "kind": self.kind,
            "label": self.label,
            "help": self.help_text,
            "required": self.required,
            "advanced": self.advanced,
            "fieldWidth": self.resolved_field_width(),
        }
        if self.default is not None:
            payload["default"] = self.default
        if self.units is not None:
            payload["units"] = self.units
        if self.group is not None:
            payload["group"] = self.group
            payload["groupCollapsed"] = self.group_collapsed
        if self.symbol is not None:
            registered = symbol(self.symbol)
            payload["symbol"] = registered.text
            payload["symbolLatex"] = registered.latex
        if self.row is not None:
            payload["row"] = self.row
        payload.update(self._describe_extra())
        return payload

    def resolved_field_width(self) -> str:
        """The width token this control will be rendered at.

        The declared ``field_width`` if there is one, otherwise whatever the
        parameter's own constraints imply. Published in the manifest so the
        frontend never has to re-derive it and the two cannot disagree.
        """

        return self.field_width or self._derived_field_width()

    def _derived_field_width(self) -> str:
        """The width implied by this parameter's constraints.

        The base answer is ``"medium"``: a kind with no declared bounds could
        hold anything, and guessing narrow there truncates a value the user can
        no longer read. Subclasses that *do* know their content narrow it.
        """

        return "medium"

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

    def _derived_field_width(self) -> str:
        """A bounded real is a short box; an unbounded one keeps the default.

        Almost every real quantity in this application is an angle, a length in
        angstrom, a voltage in kilovolts or a fraction, and every one of those
        is at most six characters with the decimals anyone types. Declaring
        bounds is therefore enough to earn the narrow box, which is why so few
        parameters need to say ``field_width`` at all.
        """

        if self.minimum is None or self.maximum is None:
            return "medium"
        span = max(abs(self.minimum), abs(self.maximum))
        return "short" if span < 1.0e4 else "medium"

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

    def _derived_field_width(self) -> str:
        """A bounded integer is exactly as wide as its widest legal value."""

        if self.minimum is None or self.maximum is None:
            return "medium"
        digits = max(len(str(self.minimum)), len(str(self.maximum)))
        if digits <= 3:
            return "index"
        if digits <= 5:
            return "tiny"
        return "medium"

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

    def _derived_field_width(self) -> str:
        """A text area takes the rail; a one-line box takes what it can use.

        ``max_length`` is the only thing free text declares about its size, so
        a short bounded string --- a two-letter element symbol, a Laue class
        name --- gets a box its content fits rather than one the rail fits.
        """

        if self.multiline:
            return "full"
        if self.max_length is not None and self.max_length <= 6:
            return "short"
        return "long"

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

    def _derived_field_width(self) -> str:
        """A checkbox carries its label beside the box and owns the line."""

        return "full"

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

    A run of digits with no separator, such as ``"110"``, is **refused**: see
    :func:`_separator_hint`. The workbench renders this parameter as one box per
    index, so the refusal is a backstop for the API rather than something a user
    of the page can reach.
    """

    width: int = 3
    allow_zero: bool = False

    kind: str = field(default="indices", init=False)

    def _describe_extra(self) -> dict[str, Any]:
        return {
            "width": self.width,
            "allowZero": self.allow_zero,
            "symbols": list(_index_symbols(self.label, self.width)),
        }

    def _derived_field_width(self) -> str:
        """One row of index boxes, each three characters wide.

        The row is the narrowest control in the application and is why the
        cardinal layout rule exists: a set of indices sized to the rail wasted
        most of a line to hold ``1 1 0``, and three of those wasted three.
        """

        return "index"

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
        extra: dict[str, Any] = {
            "width": self.width,
            "allowZero": self.allow_zero,
            "symbols": list(_index_symbols(self.label, self.width)),
        }
        if self.max_rows is not None:
            extra["maxRows"] = self.max_rows
        return extra

    def _derived_field_width(self) -> str:
        """A stack that grows downwards owns its line.

        Each row is still a row of three-character boxes; it is the *stack*,
        with its ordinals and its add and remove buttons, that needs the rail.
        """

        return "full"

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

    def _derived_field_width(self) -> str:
        """A structured document is not a value on a line."""

        return "full"

    def coerce(self, value: Any) -> dict[str, Any]:
        if not isinstance(value, Mapping):
            raise self._invalid(f"expected an object, got {type(value).__name__}.")
        return dict(value)


#: The letters a bracketed index group in a label may be written with, e.g. the
#: ``hkl`` of "Planes (hkl)" or the ``uvw`` of "Current zone axis [uvw]".
_INDEX_GROUP = re.compile(r"[(\[{<]\s*([A-Za-z]{2,4})\s*[)\]}>]")

#: What to call each index when the label does not say. Three is ``hkl``; four
#: is the Bravais-Miller ``hkil`` of a hexagonal plane, which is the only
#: four-index form the application declares.
_DEFAULT_SYMBOLS: dict[int, tuple[str, ...]] = {
    1: ("h",),
    2: ("h", "k"),
    3: ("h", "k", "l"),
    4: ("h", "k", "i", "l"),
}


def _index_symbols(label: str, width: int) -> tuple[str, ...]:
    """Name each index of a row, for the one-box-per-index control.

    Purpose
    -------
    The workbench gives every index its own box, and a box needs a name: an
    unlabelled row of three identical fields is exactly the ambiguity this
    change exists to remove. The name comes from the parameter's own label
    wherever the label already states it — "Current zone axis [uvw]" names its
    indices ``u``, ``v``, ``w`` — and otherwise from the width.

    Parameters
    ----------
    label : str
        The parameter's label, as shown above the control.
    width : int
        How many indices the row carries.

    Returns
    -------
    tuple of str
        Exactly ``width`` single-letter names.

    Examples
    --------
    >>> _index_symbols("Current zone axis [uvw]", 3)
    ('u', 'v', 'w')
    >>> _index_symbols("Planes (hkl)", 3)
    ('h', 'k', 'l')
    >>> _index_symbols("Child plane to plot", 3)
    ('h', 'k', 'l')
    >>> _index_symbols("Direction", 4)
    ('h', 'k', 'i', 'l')
    """

    match = _INDEX_GROUP.search(label)
    if match is not None and len(match.group(1)) == width:
        return tuple(match.group(1).lower())
    default = _DEFAULT_SYMBOLS.get(width)
    if default is not None:
        return default
    return tuple(f"i{position + 1}" for position in range(width))


def _coerce_index_row(value: Any, *, parameter: Parameter) -> tuple[int, ...]:
    """Parse one row of Miller indices from text or a sequence."""

    width = getattr(parameter, "width", 3)
    allow_zero = getattr(parameter, "allow_zero", False)
    items: list[Any]
    if isinstance(value, str):
        cleaned = value.strip().strip("()[]{}<>")
        cleaned = cleaned.replace(",", " ")
        items = cleaned.split()
        if len(items) == 1 and _is_compact_row(items[0]):
            raise parameter._invalid(
                f"{items[0]!r} runs the indices together, so it does not say which digits "
                "belong to which index.",
                hint=_separator_hint(items[0], width),
            )
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
    """True for a run of digits such as ``"110"`` offered as a whole index row.

    A single digit is not one of these: ``"1"`` for a width-3 row is a row with
    two indices missing, and saying so is more useful than saying it is run
    together. Anything longer is refused, including the cases that used to be
    guessed correctly — see :func:`_separator_hint` for why guessing was wrong.
    """

    return len(text) > 1 and text.isdigit()


def _separator_hint(text: str, width: int) -> str:
    """Explain why a run of digits cannot be read, using the user's own input.

    Why this is a refusal rather than a reading
    -------------------------------------------
    ``"110"`` for a width-3 row has exactly one reading, and for years it got
    it. ``"10100"`` does not: it is ``(10, 10, 0)``, ``(1, 0, 100)``,
    ``(101, 0, 0)`` and several more, and no rule distinguishes them. A parser
    that accepts the easy case teaches the habit that produces the ambiguous
    one, and the ambiguous one fails *silently* — an index nobody typed, in a
    calculation that returns a plausible number. So the run is refused whatever
    its length, and the workbench does not offer a control that can produce one.
    """

    example = " ".join(["1"] * width)
    if len(text) == width:
        spaced = " ".join(text)
        return (
            f"Separate them with spaces: {spaced}. Every index is a whole number of its own, "
            f"and a two-digit index such as 10 makes a run like {text} ambiguous, so PyTex asks "
            "rather than guesses."
        )
    return (
        f"Enter {width} whole numbers separated by spaces, for example {example}. "
        "Negative indices take a minus sign, as in 1 -1 0."
    )


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
        from pytex.core.symbols import registered_symbols

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
            # The symbol table, for the controls the frontend builds itself
            # rather than from a parameter declaration -- the cell editor of the
            # phase control is the one that matters, whose six boxes are the
            # lattice parameters and were labelled with the words "alpha",
            # "beta", "gamma". Published here for the same reason the export
            # formats are: so there is one registry, in Python, and no Greek
            # letter is ever typed into a JavaScript file where nothing can
            # check it against the standards document.
            "symbols": {
                name: {
                    "text": entry.text,
                    "latex": entry.latex,
                    "meaning": entry.meaning,
                }
                for name, entry in registered_symbols().items()
            },
        }

    def call(self, operation_id: str, request: Mapping[str, Any] | None = None) -> dict[str, Any]:
        """Validate and execute one operation."""

        spec = self.get(operation_id)
        return spec.handler(spec.bind(request or {}))


#: The application-wide registry. Importing :mod:`pytex.app.services` fills it.
REGISTRY = ServiceRegistry()
