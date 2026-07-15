"""Executable-worked-example framework for PyTex documentation.

A :class:`WorkedExample` is the atomic unit of PyTex's *documentation-as-test*
system. Each example bundles, in one place:

* a scientific *scenario* that states when and why a reader would run the code;
* a single Python *source* (``setup`` + ``code``) that is both rendered
  verbatim in the documentation and executed to produce a number;
* an *expected* value drawn from an independently-known result (an analytic
  identity, an International-Tables convention, or a cited textbook value);
* a numerical *tolerance*, provenance note, and citation;
* cross-links to the theory note, API surface, and concept page.

The same object drives three surfaces:

#. **Tests** iterate the registry and assert ``computed`` matches ``expected``
   within ``tolerance``. Documentation therefore cannot silently drift from the
   implementation.
#. **Docs** render each example into a Sphinx gallery that shows the computed
   value and the expected value side by side.
#. **Notation** is kept consistent because every example cites symbols from the
   central terminology-and-symbol registry.

Because the *displayed* code is exactly the code that is *executed*, a reader
can copy any snippet and reproduce the tabulated output verbatim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

__all__ = [
    "ExampleGroup",
    "SeeAlso",
    "SymbolUse",
    "WorkedExample",
    "WorkedExampleResult",
    "format_value",
    "validate_unique_ids",
]


def format_value(value: Any, fmt: str = "{:.4f}") -> str:
    """Render a scalar or small array as a compact, deterministic string.

    Scalars use ``fmt``; 1-D integer arrays render as bracketed integer lists;
    other arrays render element-wise with ``fmt``. The output is stable across
    platforms so it can be embedded in generated documentation and diffed.
    """

    array = np.asarray(value)
    if array.ndim == 0:
        if np.issubdtype(array.dtype, np.integer):
            return str(int(array))
        return fmt.format(float(array))
    flat = array.reshape(-1)
    if np.issubdtype(array.dtype, np.integer):
        body = ", ".join(str(int(item)) for item in flat)
    else:
        body = ", ".join(fmt.format(float(item)) for item in flat)
    return f"[{body}]"


@dataclass(frozen=True)
class SymbolUse:
    """A registry symbol referenced by an example.

    Attributes
    ----------
    latex:
        The LaTeX form exactly as fixed in the terminology-and-symbol registry.
    meaning:
        A short gloss for the reader; must not redefine the registry meaning.
    """

    latex: str
    meaning: str


@dataclass(frozen=True)
class SeeAlso:
    """A navigable cross-link emitted into the rendered example."""

    label: str
    target: str


@dataclass(frozen=True)
class WorkedExampleResult:
    """Outcome of executing and checking a :class:`WorkedExample`."""

    computed: Any
    expected: Any
    max_abs_deviation: float
    within_tolerance: bool

    @property
    def status_text(self) -> str:
        return "PASS" if self.within_tolerance else "FAIL"


@dataclass(frozen=True)
class WorkedExample:
    """A single executable, self-validating documentation example.

    Parameters
    ----------
    id:
        Stable, kebab-or-snake identifier, unique across the whole registry.
    title:
        Human-readable one-line title.
    domain:
        Registry domain slug, e.g. ``"core"``, ``"orientation"``, ``"diffraction"``.
    scenario:
        Prose that states when and why a user would compute this quantity.
    setup:
        Python source executed before ``code`` (imports and object construction).
        Rendered as a collapsible preamble so the highlighted snippet stays short.
    code:
        The highlighted Python snippet. It must bind a variable named ``result``
        holding the value that is compared against ``expected``.
    expected:
        The independently-known reference value (scalar or array-like).
    unit:
        Physical unit string for display (``""`` for dimensionless quantities).
    tolerance:
        Absolute tolerance for the ``computed`` vs ``expected`` comparison.
    reference:
        The provenance of ``expected`` (the analytic identity or standard value).
    citation:
        A normative or informative citation supporting ``reference``.
    symbols:
        Registry symbols this example uses, kept consistent with the registry.
    see_also:
        Navigable cross-links to theory, API, and concept surfaces.
    result_format:
        ``str.format`` spec used when rendering scalar/array values.
    """

    id: str
    title: str
    domain: str
    scenario: str
    setup: str
    code: str
    expected: Any
    unit: str
    tolerance: float
    reference: str
    citation: str
    symbols: tuple[SymbolUse, ...] = ()
    see_also: tuple[SeeAlso, ...] = ()
    result_format: str = "{:.4f}"

    def run(self, namespace: Mapping[str, Any] | None = None) -> Any:
        """Execute ``setup`` then ``code`` and return the bound ``result``."""

        scope: dict[str, Any] = dict(namespace) if namespace is not None else {}
        source = f"{self.setup.strip()}\n{self.code.strip()}\n"
        exec(compile(source, f"<worked-example:{self.id}>", "exec"), scope, scope)
        if "result" not in scope:
            raise AssertionError(
                f"Worked example '{self.id}' did not bind a `result` variable."
            )
        return scope["result"]

    def check(self, namespace: Mapping[str, Any] | None = None) -> WorkedExampleResult:
        """Run the example and compare the computed value against ``expected``."""

        computed = self.run(namespace)
        computed_array = np.asarray(computed, dtype=float)
        expected_array = np.asarray(self.expected, dtype=float)
        if computed_array.shape != expected_array.shape:
            return WorkedExampleResult(
                computed=computed,
                expected=self.expected,
                max_abs_deviation=float("inf"),
                within_tolerance=False,
            )
        if computed_array.size:
            deviation = float(np.max(np.abs(computed_array - expected_array)))
        else:
            deviation = 0.0
        within = bool(np.allclose(computed_array, expected_array, atol=self.tolerance, rtol=0.0))
        return WorkedExampleResult(
            computed=computed,
            expected=self.expected,
            max_abs_deviation=deviation,
            within_tolerance=within,
        )

    def computed_display(self, namespace: Mapping[str, Any] | None = None) -> str:
        return format_value(self.run(namespace), self.result_format)

    def expected_display(self) -> str:
        return format_value(self.expected, self.result_format)


@dataclass(frozen=True)
class ExampleGroup:
    """A titled collection of worked examples that render as one page."""

    slug: str
    title: str
    summary: str
    examples: tuple[WorkedExample, ...] = field(default_factory=tuple)

    def ids(self) -> list[str]:
        return [example.id for example in self.examples]


def validate_unique_ids(groups: Sequence[ExampleGroup]) -> None:
    """Raise if any example id is duplicated across the supplied groups."""

    seen: dict[str, str] = {}
    for group in groups:
        for example in group.examples:
            if example.id in seen:
                raise AssertionError(
                    f"Duplicate worked-example id '{example.id}' in groups "
                    f"'{seen[example.id]}' and '{group.slug}'."
                )
            seen[example.id] = group.slug
