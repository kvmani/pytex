from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_foundational_docs_agree_on_hybrid_documentation_policy() -> None:
    foundational_docs = [
        "README.md",
        "mission.md",
        "specifications.md",
        "AGENTS.md",
        "docs/standards/documentation_architecture.md",
    ]
    for path in foundational_docs:
        content = _read(path).lower()
        assert "sphinx" in content
        assert "latex" in content
        assert "svg" in content


def test_docs_site_placeholder_encodes_public_entry_point() -> None:
    content = _read("docs/site/README.md").lower()
    assert "public documentation entry point" in content
    assert "concept" in content
    assert "tutorial" in content
    assert "api" in content


def test_latex_standard_points_back_to_documentation_architecture() -> None:
    content = _read("docs/standards/latex_and_figures.md").lower()
    assert "documentation_architecture.md" in content
    assert "sphinx" in content
