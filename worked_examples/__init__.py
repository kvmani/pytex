"""PyTex executable worked examples: documentation that is also a live test.

This package holds the single source of truth for PyTex's computed reference
examples. Each :class:`~worked_examples.framework.WorkedExample` stores one
Python source that is both rendered verbatim in the Sphinx documentation and
executed to produce a value, which is then checked against an independently
known reference. The same objects drive:

* the regression test ``tests/unit/test_worked_examples.py``;
* the generated Sphinx gallery under ``docs/site/examples/`` produced by
  ``scripts/generate_worked_examples.py``.

Because the displayed code is the executed code, every tabulated output is
exactly reproducible by copying the snippet.
"""

from __future__ import annotations

from .framework import (
    ExampleGroup,
    SeeAlso,
    SymbolUse,
    WorkedExample,
    WorkedExampleResult,
    format_residue_scale,
    format_value,
    validate_unique_ids,
)
from .registry import all_examples, all_groups, example_by_id, iter_examples

__all__ = [
    "ExampleGroup",
    "SeeAlso",
    "SymbolUse",
    "WorkedExample",
    "WorkedExampleResult",
    "all_examples",
    "all_groups",
    "example_by_id",
    "format_residue_scale",
    "format_value",
    "iter_examples",
    "validate_unique_ids",
]
