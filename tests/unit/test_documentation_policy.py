from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_foundational_docs_agree_on_layered_documentation_policy() -> None:
    """Every foundational doc must name the same three documentation pillars.

    Sphinx is the browsable surface, the scientific notes are canonical prose
    and mathematics, and SVG is canonical for figures. The notes were LaTeX
    until 2026-08-11; they are now MyST under docs/site/theory/ so that the
    derivations render on the site instead of downloading as .tex sources.
    """
    foundational_docs = [
        "README.md",
        "mission.md",
        "specifications.md",
        "AGENTS.md",
        "docs/standards/documentation_architecture.md",
        "docs/standards/terminology_and_symbol_registry.md",
        "docs/standards/visualization_style_guide.md",
    ]
    for path in foundational_docs:
        content = _read(path).lower()
        assert "sphinx" in content, path
        assert "notes" in content, path
        assert "svg" in content, path


def test_no_foundational_doc_still_calls_latex_canonical() -> None:
    """The retired policy must not survive anywhere in the governing docs."""
    for path in (
        "README.md",
        "mission.md",
        "specifications.md",
        "AGENTS.md",
        "docs/README.md",
        "docs/standards/scientific_notes_and_figures.md",
        "docs/standards/documentation_architecture.md",
    ):
        content = _read(path).lower()
        assert "docs/tex/" not in content or "until 2026-08-11" in content, path
        assert "latex is the canonical" not in content, path
        assert "canonically in latex" not in content, path


def test_docs_site_placeholder_encodes_public_entry_point() -> None:
    content = _read("docs/site/README.md").lower()
    assert "public documentation entry point" in content
    assert "concept" in content
    assert "tutorial" in content
    assert "api" in content
    assert "notebook" in content


def test_notes_standard_points_back_to_documentation_architecture() -> None:
    content = _read("docs/standards/scientific_notes_and_figures.md").lower()
    assert "documentation_architecture.md" in content
    assert "visualization_style_guide.md" in content
    assert "sphinx" in content


def test_visualization_style_guide_defines_canonical_svg_policy() -> None:
    content = _read("docs/standards/visualization_style_guide.md").lower()
    for token in ("#07122f", "#2563eb", "#7c3aed", "mermaid", "title", "desc"):
        assert token in content


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
        "docs/standards/scientific_notes_and_figures.md",
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


def test_executable_examples_standard_is_encoded_and_cross_linked() -> None:
    standard = _read("docs/standards/executable_examples.md").lower()
    for token in ("worked example", "computed", "reference value", "tolerance", "citation"):
        assert token in standard, token

    agents = _read("AGENTS.md")
    assert "docs/standards/executable_examples.md" in agents
    assert "worked_examples/" in agents

    doc_arch = _read("docs/standards/documentation_architecture.md").lower()
    assert "executable worked example" in doc_arch
    assert "docstring contract" in doc_arch

    site_index = _read("docs/site/index.md")
    assert "examples/index" in site_index

    standards_index = _read("docs/site/standards/index.md")
    assert "executable_examples" in standards_index


def test_foundational_docs_encode_executable_example_policy() -> None:
    for path in ("mission.md", "specifications.md", "AGENTS.md"):
        content = _read(path).lower()
        assert "worked example" in content, path


def test_interop_docs_state_validation_boundaries_explicitly() -> None:
    interop = _read("docs/site/workflows/orix_kikuchipy_interop.md").lower()
    ebsd = _read("docs/site/workflows/ebsd_import_normalization.md").lower()
    assert "what is executable today" in interop
    assert "what is not being claimed" in interop
    assert "current limits" in interop
    assert "what is verified today" in ebsd
    assert "interpretation rule" in ebsd
