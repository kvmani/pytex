"""Executable-worked-example regression tests.

These tests make the documentation gallery self-validating:

* every registered example must compute a value that agrees with its
  independently known reference within its declared tolerance;
* the generated Sphinx pages must be byte-identical to a fresh render, so the
  published numbers can never drift from the live code.

Because the examples only touch the pure-Python ``pytex`` core plus NumPy, this
module runs without any optional scientific dependency.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.generate_worked_examples import (
    EXAMPLES_DIR,
    GENERATED_DIR,
    render_group,
    render_index,
)
from worked_examples import all_examples, all_groups

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_registry_is_non_empty_and_ids_unique() -> None:
    groups = all_groups()
    assert groups, "The worked-example registry must not be empty."
    examples = all_examples()
    assert len(examples) >= 8
    ids = [example.id for example in examples]
    assert len(ids) == len(set(ids)), "Worked-example ids must be unique."


@pytest.mark.parametrize("example", all_examples(), ids=lambda ex: ex.id)
def test_worked_example_matches_expected(example) -> None:  # type: ignore[no-untyped-def]
    result = example.check()
    assert result.within_tolerance, (
        f"Worked example '{example.id}' diverged from its reference value: "
        f"computed={example.computed_display()} expected={example.expected_display()} "
        f"max_abs_deviation={result.max_abs_deviation:.3e} tolerance={example.tolerance:.3e}"
    )


def test_every_example_binds_a_result() -> None:
    for example in all_examples():
        value = example.run()
        assert value is not None, f"Example '{example.id}' produced a None result."


def test_generated_group_pages_are_current() -> None:
    for group in all_groups():
        path = GENERATED_DIR / f"{group.slug}.md"
        assert path.exists(), f"Missing generated page for group '{group.slug}'."
        expected = render_group(group)
        actual = path.read_text(encoding="utf-8")
        assert actual == expected, (
            f"Generated page for '{group.slug}' is stale. "
            "Run `python scripts/generate_worked_examples.py`."
        )


def test_generated_index_is_current() -> None:
    index_path = EXAMPLES_DIR / "index.md"
    assert index_path.exists()
    assert index_path.read_text(encoding="utf-8") == render_index(), (
        "Worked-examples index is stale. Run `python scripts/generate_worked_examples.py`."
    )


def test_worked_examples_are_wired_into_site_index() -> None:
    site_index = (REPO_ROOT / "docs" / "site" / "index.md").read_text(encoding="utf-8")
    assert "examples/index" in site_index


def test_examples_declare_provenance_and_scenario() -> None:
    for example in all_examples():
        assert example.scenario.strip(), f"Example '{example.id}' needs a scenario."
        assert example.reference.strip(), f"Example '{example.id}' needs a reference derivation."
        assert example.citation.strip(), f"Example '{example.id}' needs a citation."
