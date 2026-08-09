from __future__ import annotations

from pathlib import Path

import pytest

from scripts import check_repo_integrity
from scripts.check_repo_integrity import _check_repository_content, main

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_repo_integrity_script_passes() -> None:
    assert main() == 0


#: One tracked path per excluded category, and the substring the report must name.
VIOLATIONS = [
    ("docs/site/_build/html/index.html", "build output"),
    ("docs/_build/doctrees/index.doctree", "build output"),
    ("src/pytex/__pycache__/lattice.cpython-311.pyc", "bytecode"),
    ("docs/site/tutorials/notebooks/.ipynb_checkpoints/01-checkpoint.ipynb", "checkpoint"),
    ("inspection_outputs/atlas_preview.html", "inspection"),
    ("outputs/run_2026_08_09.json", "run output"),
    ("benchmarks/run.log", "run log"),
    ("scratch.tmp", "scratch"),
    ("src/pytex/core/lattice.py.bak", "backup"),
    ("references/a_new_paper.pdf", "reference_index"),
]


@pytest.mark.parametrize(("tracked", "expected"), VIOLATIONS)
def test_regenerable_and_scratch_artifacts_are_rejected(
    tracked: str, expected: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The cardinal repository-content rule in AGENTS.md, enforced rather than remembered.

    An artifact committed once is carried by every clone forever, even after it
    is deleted, so the only effective moment to catch it is the commit that
    introduces it.
    """

    monkeypatch.setattr(check_repo_integrity, "_tracked_files", lambda _root: (tracked,))
    issues = _check_repository_content(REPO_ROOT)
    assert issues, f"{tracked} should have been rejected"
    assert expected in issues[0]


@pytest.mark.parametrize(
    "tracked",
    [
        "src/pytex/core/lattice.py",
        "docs/figures/class_model_core.svg",
        "docs/site/architecture/class_model_atlas.md",
        "docs/site/examples/index.md",
        "fixtures/mtex_parity/cases/cubic.json",
        "tests/unit/test_repo_integrity.py",
    ],
)
def test_sources_and_canonical_assets_are_accepted(
    tracked: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Generated files that documentation or tests name by hand stay committed.

    The rule excludes what is *regenerable and unreferenced*, not everything a
    script writes: the class-model and reference-frame SVGs are compared
    byte-for-byte by the suite and embedded by name in the docs.
    """

    monkeypatch.setattr(check_repo_integrity, "_tracked_files", lambda _root: (tracked,))
    assert _check_repository_content(REPO_ROOT) == []


def test_the_grandfathered_pdf_list_does_not_grow() -> None:
    """Eight PDFs predate the reference-corpus rule; no ninth may be added.

    Untracking them is a separate decision — the bytes are in history either way
    — but the list is frozen so the exception cannot quietly become the rule.
    """

    assert len(check_repo_integrity._GRANDFATHERED_PDFS) == 8
    assert all(
        path.startswith("references/") and path.endswith(".pdf")
        for path in check_repo_integrity._GRANDFATHERED_PDFS
    )
