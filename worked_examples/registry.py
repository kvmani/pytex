"""Aggregation and validation of the worked-example registry."""

from __future__ import annotations

from collections.abc import Iterator

from .examples import GROUPS
from .framework import ExampleGroup, WorkedExample, validate_unique_ids

__all__ = [
    "all_examples",
    "all_groups",
    "example_by_id",
    "iter_examples",
]


def all_groups() -> tuple[ExampleGroup, ...]:
    """Return every registered example group, after checking id uniqueness."""

    validate_unique_ids(GROUPS)
    return GROUPS


def iter_examples() -> Iterator[WorkedExample]:
    """Yield every registered worked example across all groups."""

    for group in all_groups():
        yield from group.examples


def all_examples() -> tuple[WorkedExample, ...]:
    """Return every registered worked example as a tuple."""

    return tuple(iter_examples())


def example_by_id(example_id: str) -> WorkedExample:
    """Return the worked example with ``example_id`` or raise ``KeyError``."""

    for example in iter_examples():
        if example.id == example_id:
            return example
    raise KeyError(f"No worked example registered with id '{example_id}'.")
