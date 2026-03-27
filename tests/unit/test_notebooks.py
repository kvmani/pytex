from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK_ROOT = REPO_ROOT / "docs" / "site" / "tutorials" / "notebooks"


def test_tutorial_notebook_catalog_exists_and_is_populated() -> None:
    notebook_paths = sorted(NOTEBOOK_ROOT.glob("*.ipynb"))
    assert len(notebook_paths) >= 9


def test_each_tutorial_notebook_contains_markdown_and_code_cells() -> None:
    for notebook_path in sorted(NOTEBOOK_ROOT.glob("*.ipynb")):
        payload = json.loads(notebook_path.read_text(encoding="utf-8"))
        assert payload["nbformat"] == 4
        cell_types = [cell["cell_type"] for cell in payload["cells"]]
        assert "markdown" in cell_types, notebook_path.name
        assert "code" in cell_types, notebook_path.name


def test_notebook_index_references_all_tutorial_notebooks() -> None:
    content = (REPO_ROOT / "docs" / "site" / "tutorials" / "notebooks.md").read_text(
        encoding="utf-8"
    )
    for notebook_path in sorted(NOTEBOOK_ROOT.glob("*.ipynb")):
        stem = notebook_path.stem
        assert stem in content
