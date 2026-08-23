# Local Development

PyTex targets Python `3.11+` and should remain friendly to normal local development on macOS, Linux, and Windows.

## Bootstrap

```bash
python -m pip install -e '.[dev,docs]'
python scripts/check_repo_integrity.py
python -m ruff check .
python -m mypy src
python -m pytest -q
python scripts/check_sphinx_warnings.py --max-warnings 0
```

This is the base lane. It is the default contributor environment for integrity checks, docs builds,
type checking, linting, and the lightweight test suite.

The workbench's critical real-browser journeys use a test-only Node dependency; they do not add a
JavaScript build or a browser runtime dependency:

```bash
npm ci
npx playwright install chromium
npm run test:browser
```

Generated `node_modules/`, Playwright reports, screenshots, and traces stay outside repository
history. CI installs Chromium in a dedicated Ubuntu/Python 3.11/Node 22 lane.

Install the full scientific lane when you need the optional CIF-backed structure-import path and
the heavier interoperability or external-baseline tests:

```bash
python -m pip install -e '.[dev,docs,adapters]'
python -m pytest -q -rs
```

For a fuller user-facing setup guide, including Windows activation details, notebook use, and PDF build notes, see [Installation And Build](../site/tutorials/installation_and_build.md).

## Working Expectations

- keep the core package importable without heavy optional scientific dependencies
- use optional extras for adapters
- keep generated artifacts out of source directories
- treat documentation and figures as repo-tracked assets, not temporary outputs
