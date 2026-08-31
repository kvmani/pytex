"""Render the worked-example registry into the Sphinx documentation gallery.

This script is the *documentation leg* of the executable-worked-examples system.
It executes every registered example against the live PyTex API and writes the
computed value alongside the independently-known expected value into MyST
Markdown pages under ``docs/site/examples/``.

The output is deterministic: it contains no timestamps and formats all numbers
through :func:`worked_examples.framework.format_value`, so regenerating on an
unchanged registry produces byte-identical files. ``tests/unit/test_worked_examples.py``
relies on that property to detect stale generated documentation.

Run manually with::

    python scripts/generate_worked_examples.py

or import :func:`render_group` / :func:`render_index` for testing.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from worked_examples import (  # noqa: E402
    ExampleGroup,
    WorkedExample,
    all_groups,
    format_residue_scale,
)

EXAMPLES_DIR = REPO_ROOT / "docs" / "site" / "examples"
GENERATED_DIR = EXAMPLES_DIR / "generated"

_BANNER = (
    "<!-- GENERATED FILE. Do not edit by hand.\n"
    "     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).\n"
    "     Run `python scripts/generate_worked_examples.py` to regenerate. -->\n"
)

_INTRO = (
    "```{note}\n"
    "Every number on this page is computed live from the public PyTex API when the documentation "
    "is regenerated, then checked against an independently known reference value by "
    "`tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the "
    "computed value, so you can copy any snippet and reproduce the tabulated output.\n"
    "```\n"
)


def _render_symbols(example: WorkedExample) -> str:
    if not example.symbols:
        return ""
    lines = ["**Symbols**", ""]
    for symbol in example.symbols:
        lines.append(f"- ${symbol.latex}$ &mdash; {symbol.meaning}")
    lines.append("")
    return "\n".join(lines) + "\n"


def _render_see_also(example: WorkedExample) -> str:
    if not example.see_also:
        return ""
    links = ", ".join(f"{{doc}}`{link.label} <{link.target}>`" for link in example.see_also)
    return f"**See also**: {links}\n"


def _render_result_table(example: WorkedExample) -> str:
    result = example.check()
    unit = example.unit if example.unit else "&mdash;"
    computed = example.computed_display()
    expected = example.expected_display()
    if example.tolerance == 0.0:
        deviation = "exact"
        tolerance = "exact"
    else:
        # The deviation of an exact identity is rounding, and its digits differ
        # between platforms; see worked_examples.framework.format_residue_scale.
        deviation = format_residue_scale(result.max_abs_deviation) or (
            f"{result.max_abs_deviation:.2e}"
        )
        tolerance = f"{example.tolerance:.0e}"
    status = "✅ pass" if result.within_tolerance else "❌ FAIL"
    header = "| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |"
    divider = "| --- | --- | --- | --- | --- | --- | --- |"
    row = f"| `{example.id}` | {computed} | {expected} | {unit} | {deviation} | {tolerance} | {status} |"
    return "\n".join([header, divider, row]) + "\n"


def render_example(example: WorkedExample) -> str:
    parts: list[str] = [f"## {example.title}\n"]
    parts.append(example.scenario + "\n")
    symbols = _render_symbols(example)
    if symbols:
        parts.append(symbols)
    parts.append(":::{dropdown} Setup (imports and object construction)\n")
    parts.append("```python\n" + example.setup.strip() + "\n```\n")
    parts.append(":::\n")
    parts.append("**Compute**\n")
    parts.append("```python\n" + example.code.strip() + "\n```\n")
    parts.append("**Result**\n")
    parts.append(_render_result_table(example))
    parts.append(f"**Why this value**: {example.reference}\n")
    parts.append(f"**Citation**: {example.citation}\n")
    see_also = _render_see_also(example)
    if see_also:
        parts.append(see_also)
    return "\n".join(parts)


def render_group(group: ExampleGroup) -> str:
    parts: list[str] = [_BANNER, f"# {group.title}\n", group.summary + "\n", _INTRO]
    for example in group.examples:
        parts.append(render_example(example))
    return "\n".join(parts).rstrip() + "\n"


def _ledger_rows() -> list[str]:
    rows: list[str] = []
    for group in all_groups():
        for example in group.examples:
            result = example.check()
            status = "✅" if result.within_tolerance else "❌"
            unit = example.unit if example.unit else "&mdash;"
            rows.append(
                f"| `{example.id}` | {example.computed_display()} | "
                f"{example.expected_display()} | {unit} | {status} |"
            )
    return rows


def render_index() -> str:
    groups = all_groups()
    toctree_entries = "\n".join(f"generated/{group.slug}" for group in groups)
    ledger = "\n".join(
        [
            "| Example | Computed (live) | Expected (reference) | Unit | Status |",
            "| --- | --- | --- | --- | --- |",
            *_ledger_rows(),
        ]
    )
    section_list = "\n".join(
        f"- {{doc}}`{group.title} <generated/{group.slug}>` &mdash; {group.summary}"
        for group in groups
    )
    return (
        f"{_BANNER}"
        "# Worked Examples\n\n"
        "This section is PyTex's *documentation-as-test* surface. Each worked example bundles a "
        "scientific scenario, a runnable snippet, the value that snippet computes from the live "
        "code, and an independently known reference value. The examples are the single source of "
        "truth for both this gallery and the regression test "
        "`tests/unit/test_worked_examples.py`.\n\n"
        f"{_INTRO}\n"
        "## Reference-value ledger\n\n"
        "The complete set of computed-versus-expected values at a glance:\n\n"
        f"{ledger}\n\n"
        "## Example groups\n\n"
        f"{section_list}\n\n"
        "```{toctree}\n"
        ":maxdepth: 1\n:hidden:\n\n"
        f"{toctree_entries}\n"
        "```\n"
    )


def write_docs() -> list[Path]:
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for group in all_groups():
        path = GENERATED_DIR / f"{group.slug}.md"
        path.write_text(render_group(group), encoding="utf-8")
        written.append(path)
    index_path = EXAMPLES_DIR / "index.md"
    index_path.write_text(render_index(), encoding="utf-8")
    written.append(index_path)
    return written


def main() -> None:
    written = write_docs()
    for path in written:
        print(f"wrote {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
