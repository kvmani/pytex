"""Execute tutorial notebooks in place so committed outputs render on the site.

The ``.ipynb`` files under ``docs/site/tutorials/notebooks`` are the source of
truth: they are hand-authored and edited directly, like any other source file.
There is no generator step.

The Sphinx site builds with ``nb_execution_mode = "off"`` (myst-nb renders the
stored outputs), so a notebook's figures and printed results appear on the site
only if it is committed *executed*. After editing a notebook, run::

    python scripts/execute_notebooks.py --only 21

``--only`` matches filename prefixes; omit it to execute every notebook.
``tests/unit/test_notebooks.py`` enforces that committed notebooks are executed
and error-free, so a forgotten run fails CI rather than silently shipping a
tutorial page with no outputs.

``--only`` matches filename prefixes; with no argument every notebook under
``docs/site/tutorials/notebooks`` is executed. Execution runs with the
repository root as working directory and a headless Matplotlib backend.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_ROOT = REPO_ROOT / "docs" / "site" / "tutorials" / "notebooks"


def execute_notebook(path: Path, timeout: int) -> None:
    notebook = nbformat.read(path, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(REPO_ROOT)}},
    )
    client.execute()
    nbformat.write(notebook, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only",
        nargs="*",
        default=None,
        help="filename prefixes to execute (e.g. 18 19); default: all notebooks",
    )
    parser.add_argument("--timeout", type=int, default=600)
    arguments = parser.parse_args()

    # Let ipykernel select its inline backend so figures are captured as
    # notebook outputs (forcing Agg would silently drop them).
    os.environ.pop("MPLBACKEND", None)
    notebook_paths = sorted(NOTEBOOK_ROOT.glob("*.ipynb"))
    if arguments.only:
        notebook_paths = [
            path
            for path in notebook_paths
            if any(path.name.startswith(prefix) for prefix in arguments.only)
        ]
    if not notebook_paths:
        print("No notebooks matched.", file=sys.stderr)
        return 1
    for path in notebook_paths:
        print(f"Executing {path.relative_to(REPO_ROOT)} ...")
        execute_notebook(path, arguments.timeout)
    print(f"Executed {len(notebook_paths)} notebook(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
