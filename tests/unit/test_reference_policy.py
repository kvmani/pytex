from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_foundational_docs_declare_normative_and_informative_references() -> None:
    foundational_docs = [
        "mission.md",
        "specifications.md",
        "docs/architecture/overview.md",
        "docs/architecture/canonical_data_model.md",
        "docs/architecture/orientation_and_texture_foundation.md",
        "docs/architecture/ebsd_foundation.md",
        "docs/architecture/diffraction_foundation.md",
        "docs/architecture/multimodal_characterization_foundation.md",
        "docs/architecture/phase_transformation_foundation.md",
        "docs/architecture/repo_review_2026_foundation_audit.md",
        "docs/standards/engineering_governance.md",
        "docs/standards/notation_and_conventions.md",
        "docs/standards/documentation_architecture.md",
        "docs/standards/latex_and_figures.md",
        "docs/standards/scientific_citation_policy.md",
        "docs/standards/development_principles.md",
        "docs/standards/data_contracts_and_manifests.md",
        "docs/standards/reference_canon.md",
        "docs/testing/strategy.md",
        "docs/testing/mtex_parity_matrix.md",
        "docs/testing/diffraction_validation_matrix.md",
        "docs/roadmap/implementation_roadmap.md",
    ]
    for path in foundational_docs:
        content = _read(path)
        assert "### Normative" in content, path
        assert "### Informative" in content, path
