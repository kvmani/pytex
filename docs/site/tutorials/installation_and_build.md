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

Install the full local development toolchain:

```bash
python -m pip install -e ".[dev,docs]"
```
This standard development install includes the CIF-backed structure-import stack used by the
default test suite. The broader `adapters` extra remains optional for heavier interoperability
work beyond the normal contributor workflow.

## Install On Windows

Use PowerShell from the repo root:

```powershell
py -3.11 -m venv .venv
.venv\\Scripts\\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev,docs]"
```
This standard development install includes the CIF-backed structure-import stack used by the
default test suite. The broader `adapters` extra remains optional for heavier interoperability
work beyond the normal contributor workflow.

If PowerShell script execution is blocked, activate from `cmd.exe` with:

```bat
.venv\\Scripts\\activate.bat
```

## Run The Core Quality Gates

These commands should work the same way on Windows, macOS, and Linux once the environment is active:

```bash
python scripts/check_repo_integrity.py
ruff check .
mypy src
pytest -q
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

The default development install already includes the CIF-backed structure-import support used by
the normal test suite. For heavier optional interoperability adapters, install:

```bash
python -m pip install -e ".[adapters]"
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

- `../../standards/documentation_architecture.md`
- `../../standards/latex_and_figures.md`

### Informative

- `../../development/local_development.md`
- `../README.md`
