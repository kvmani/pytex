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
        "docs/standards/terminology_and_symbol_registry.md",
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
    assert "notebook" in content


def test_latex_standard_points_back_to_documentation_architecture() -> None:
    content = _read("docs/standards/latex_and_figures.md").lower()
    assert "documentation_architecture.md" in content
    assert "sphinx" in content


def test_foundational_docs_encode_notebook_policy() -> None:
    foundational_docs = [
        "mission.md",
        "specifications.md",
        "docs/standards/documentation_architecture.md",
        "docs/standards/development_principles.md",
    ]
    for path in foundational_docs:
        content = _read(path).lower()
        assert "notebook" in content, path


def test_foundational_docs_encode_dual_plotting_output_policy() -> None:
    policy_docs = [
        "specifications.md",
        "docs/standards/documentation_architecture.md",
        "docs/standards/latex_and_figures.md",
    ]
    for path in policy_docs:
        content = _read(path).lower()
        assert "matplotlib" in content, path
        assert "svg" in content, path


def test_installation_guide_covers_platforms_and_docs_builds() -> None:
    content = _read("docs/site/tutorials/installation_and_build.md").lower()
    for token in ("windows", "macos", "linux", "sphinx", "latex", "jupyter"):
        assert token in content


def test_site_index_references_glossary_and_installation_pages() -> None:
    content = _read("docs/site/index.md")
    assert "concepts/technical_glossary_and_symbols" in content
    assert "tutorials/installation_and_build" in content


def test_site_index_and_readme_expose_active_roadmap() -> None:
    site_index = _read("docs/site/index.md").lower()
    site_readme = _read("docs/site/README.md").lower()
    assert "roadmap/implementation_roadmap" in site_index
    assert "implementation roadmap" in site_readme


def test_interop_docs_state_validation_boundaries_explicitly() -> None:
    interop = _read("docs/site/workflows/orix_kikuchipy_interop.md").lower()
    ebsd = _read("docs/site/workflows/ebsd_import_normalization.md").lower()
    assert "what is executable today" in interop
    assert "what is not being claimed" in interop
    assert "current limits" in interop
    assert "what is verified today" in ebsd
    assert "interpretation rule" in ebsd
