# Installation And Build

This page explains how to install PyTex, run the test suite, build the Sphinx documentation, and generate the LaTeX or PDF scientific notes on Windows, macOS, and Linux.

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

## Regenerate Tutorial Notebooks

The checked-in notebooks are generated from the canonical notebook generator:

```bash
python scripts/generate_tutorial_notebooks.py
```

After regenerating notebooks, rebuild the Sphinx site so the rendered notebook pages stay synchronized with the runtime API.

## Build LaTeX Or PDF Notes

The canonical scientific notes live under `docs/tex/`.

### Basic LaTeX expectation

You need a modern LaTeX distribution such as:

- TeX Live on Linux
- MacTeX on macOS
- MiKTeX or TeX Live on Windows

### Typical PDF build

From `docs/tex/`, use a modern toolchain such as:

```bash
latexmk -pdf -shell-escape <document>.tex
```

If your TeX setup lacks the required packages, install them through your distribution package manager or TeX package manager first.

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

### Sphinx builds but notebook content looks stale

Regenerate the notebooks first:

```bash
python scripts/generate_tutorial_notebooks.py
```

### LaTeX PDF build fails

Check that your TeX distribution includes `latexmk` and any packages needed by the note you are building.

## Related Material

- {doc}`notebooks`
- {doc}`../concepts/technical_glossary_and_symbols`
- {doc}`../theory/index`

## References

### Normative

- {doc}`../standards/documentation_architecture`
- {doc}`../standards/latex_and_figures`

### Informative

- <a href="../../development/local_development.md">Local Development</a>
- <a href="../README.md">Sphinx Site README</a>
