from __future__ import annotations

import contextlib
import json
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK_ROOT = REPO_ROOT / "docs" / "site" / "tutorials" / "notebooks"
PRIORITY_NOTEBOOKS = (
    "04_phases_lattices_space_groups_and_cif.ipynb",
    "11_powder_xrd_workflows.ipynb",
    "12_saed_workflows.ipynb",
    "13_crystal_visualization_workflows.ipynb",
    "15_structure_diffraction_visualization_pipeline.ipynb",
)


def _notebook_payload(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _cell_source(cell: dict[str, object]) -> str:
    source = cell["source"]
    if isinstance(source, list):
        return "".join(str(fragment) for fragment in source)
    return str(source)


def _execute_notebook_code_cells(path: Path) -> None:
    payload = _notebook_payload(path)
    namespace: dict[str, object] = {"__name__": "__notebook__"}
    with contextlib.chdir(REPO_ROOT):
        for cell in payload["cells"]:
            if cell["cell_type"] != "code":
                continue
            source = _cell_source(cell)
            if not source.strip():
                continue
            exec(compile(source, f"{path.name}", "exec"), namespace, namespace)
    plt.close("all")


def test_tutorial_notebook_catalog_exists_and_is_populated() -> None:
    notebook_paths = sorted(NOTEBOOK_ROOT.glob("*.ipynb"))
    assert len(notebook_paths) >= 9


def test_each_tutorial_notebook_contains_markdown_and_code_cells() -> None:
    for notebook_path in sorted(NOTEBOOK_ROOT.glob("*.ipynb")):
        payload = _notebook_payload(notebook_path)
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


def test_priority_teaching_notebooks_use_fixture_corpus_and_manifest_trail() -> None:
    for notebook_name in PRIORITY_NOTEBOOKS:
        content = (NOTEBOOK_ROOT / notebook_name).read_text(encoding="utf-8")
        assert "get_phase_fixture" in content or "list_phase_fixtures" in content

    pipeline_content = (
        NOTEBOOK_ROOT / "15_structure_diffraction_visualization_pipeline.ipynb"
    ).read_text(encoding="utf-8")
    assert "read_workflow_result_manifest" in pipeline_content
    assert "read_validation_manifest" in pipeline_content


def test_priority_teaching_notebooks_smoke_execute() -> None:
    pytest.importorskip("pymatgen.core")
    for notebook_name in PRIORITY_NOTEBOOKS:
        _execute_notebook_code_cells(NOTEBOOK_ROOT / notebook_name)


def test_no_notebook_is_committed_with_outputs() -> None:
    """A committed notebook must carry its source only, never its outputs.

    Outputs are the record of one particular run on one particular machine.
    Committing them makes every rerun a diff, makes review impossible (a
    one-line code change arrives as a thousand lines of changed base64), and
    inflates the repository without adding a single reviewable fact — the
    stored figures were 13 MB against 0.45 MB of actual notebook source.

    The site does not need them: ``docs/site/conf.py`` sets
    ``nb_execution_mode = "cache"``, so myst-nb runs each notebook at build
    time and renders what it produces. A notebook that no longer runs fails the
    docs build, which is a stronger guarantee than a stored output could give —
    a stored output proves only that the notebook ran once, against whatever
    the library looked like then.
    """

    offenders: list[str] = []
    for notebook_path in sorted(NOTEBOOK_ROOT.glob("*.ipynb")):
        payload = _notebook_payload(notebook_path)
        for index, cell in enumerate(payload["cells"]):
            if cell["cell_type"] != "code":
                continue
            if cell.get("outputs"):
                offenders.append(f"{notebook_path.name} cell {index}: stored outputs")
            if cell.get("execution_count") is not None:
                offenders.append(f"{notebook_path.name} cell {index}: execution count")
    assert not offenders, (
        "These notebook cells were committed with run artifacts. Clear all outputs before "
        "committing (in VS Code: 'Notebook: Clear All Outputs'; in Jupyter: Kernel > Restart "
        "and Clear Output). Offenders: " + ", ".join(offenders)
    )


def test_no_notebook_carries_run_specific_metadata() -> None:
    """Per-run metadata is as unreviewable as an output, and just as transient.

    ``metadata.execution`` records wall-clock timestamps of one run, and
    notebook-level ``widgets`` state records the state of interactive widgets in
    one session. Both change on every execution while meaning nothing to a
    reader, so they are stripped alongside outputs.
    """

    offenders: list[str] = []
    for notebook_path in sorted(NOTEBOOK_ROOT.glob("*.ipynb")):
        payload = _notebook_payload(notebook_path)
        notebook_metadata = payload.get("metadata", {})
        assert isinstance(notebook_metadata, dict)
        if "widgets" in notebook_metadata:
            offenders.append(f"{notebook_path.name}: notebook widget state")
        for index, cell in enumerate(payload["cells"]):
            cell_metadata = cell.get("metadata", {})
            assert isinstance(cell_metadata, dict)
            if "execution" in cell_metadata:
                offenders.append(f"{notebook_path.name} cell {index}: execution timings")
    assert not offenders, (
        "These notebooks carry run-specific metadata: " + ", ".join(offenders)
    )
