# Installation And Build

This page explains how to install PyTex, run the test suite, build the Sphinx documentation, and produce a PDF of it on Windows, macOS, and Linux.

## Supported Python

PyTex currently targets Python `3.11+`.

## Common Repository Setup

Clone the repository, then move into it:

```bash
git clone <your-pytex-repo-url>
cd pytex
```

## Install On macOS Or Linux

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

Install the base contributor lane:

```bash
python -m pip install -e ".[dev,docs]"
```
This base lane covers repository integrity, docs builds, linting, type checking, and the default
lightweight test suite.

## Install On Windows

Use PowerShell from the repo root:

```powershell
py -3.11 -m venv .venv
.venv\\Scripts\\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev,docs]"
```
This base lane covers repository integrity, docs builds, linting, type checking, and the default
lightweight test suite.

If PowerShell script execution is blocked, activate from `cmd.exe` with:

```bat
.venv\\Scripts\\activate.bat
```

## Contributor Lanes

### Base Lane

Use the base lane for normal contributor work:

```bash
python scripts/check_repo_integrity.py
python -m ruff check .
python -m mypy src
python -m pytest -q
python -m sphinx -b html docs/site docs/_build/html
```

### Full Scientific Lane

Use the full scientific lane when you need optional structure-import or interoperability coverage:

```bash
python -m pip install -e ".[dev,docs,adapters]"
python -m pytest -q -rs
```

The full scientific lane now runs without the previous `pymatgen`-gated skips and is the
controlling environment for CIF-backed phase construction, pinned diffraction external baselines,
and the heavier notebook smoke path.

## Run The Core Quality Gates

These commands should work the same way on Windows, macOS, and Linux once the environment is active:

```bash
python scripts/check_repo_integrity.py
python -m ruff check .
python -m mypy src
python -m pytest -q
```

## Build The Sphinx HTML Docs

```bash
sphinx-build -b html docs/site docs/_build/html
```

Open the built site at:

- `docs/_build/html/index.html`

See also: {doc}`../concepts/library_structure`, {doc}`../concepts/technical_glossary_and_symbols`.

## Work On Tutorial Notebooks

The `.ipynb` files under `docs/site/tutorials/notebooks/` are the source of truth. Edit them
directly — in Jupyter, in your editor, or by hand — like any other source file. There is no
generator step.

Notebooks are hand-authored and committed **without outputs**. The Sphinx site builds with
`nb_execution_mode = "cache"`, so myst-nb executes each notebook itself and renders what it
produces; `nb_execution_raise_on_error = True` means a tutorial that no longer runs fails the
build rather than publishing a traceback as a result.

So after editing a notebook there is nothing to regenerate — clear its outputs and commit the
source. `tests/unit/test_notebooks.py` fails if any output, execution count, or run-specific
metadata is committed.

In VS Code, clear outputs with the command palette entry **Notebook: Clear All Outputs**; in
Jupyter, **Kernel > Restart Kernel and Clear All Outputs**.

## Build A PDF

The canonical scientific notes live under `docs/site/theory/` as ordinary MyST pages, so they
need no separate toolchain to read: they render with the rest of the site. A typeset PDF of the
whole documentation set comes from Sphinx's own LaTeX builder:

```bash
python -m sphinx -b latexpdf docs/site docs/_build/latex
```

This is the only supported PDF path. Because it renders the same sources as the HTML site, the
print and web forms cannot drift apart.

### LaTeX distribution

`latexpdf` shells out to a real LaTeX installation, so that step needs one:

- TeX Live on Linux
- MacTeX on macOS
- MiKTeX or TeX Live on Windows

If your TeX setup lacks the required packages, install them through your distribution package
manager or TeX package manager first. To produce the intermediate `.tex` without running LaTeX,
use `-b latex` instead and build it yourself.

## Optional Jupyter Use

To run the checked-in notebooks locally:

```bash
python -m jupyter lab docs/site/tutorials/notebooks
```

## Troubleshooting

### `ImportError` from broader optional adapters

The base lane intentionally excludes the heavier optional scientific adapters. Install the full
scientific lane when you need CIF-backed phase loading, ORIX bridges, or other adapter-heavy
validation:

```bash
python -m pip install -e ".[dev,docs,adapters]"
```

### Sphinx builds but a notebook page shows code with no outputs

The site executes notebooks itself, so this means execution was skipped or served from a stale
cache. Clear the cache and rebuild:

```bash
rm -rf docs/site/_build && python -m sphinx -b html docs/site docs/site/_build/html
```

### PDF build fails

Check that your TeX distribution includes `latexmk`, which Sphinx's `latexpdf` builder drives, and
the packages Sphinx's LaTeX output requires. To isolate the failure, run `-b latex` first: if that
succeeds, the problem is in the TeX toolchain rather than in the documentation sources.

## Related Material

- {doc}`notebooks`
- {doc}`../concepts/technical_glossary_and_symbols`
- {doc}`../theory/index`

## References

### Normative

- {doc}`../standards/documentation_architecture`
- {doc}`../standards/scientific_notes_and_figures`

### Informative

- <a href="../../development/local_development.md">Local Development</a>
- <a href="../README.md">Sphinx Site README</a>
