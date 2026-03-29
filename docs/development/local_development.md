# Local Development

PyTex targets Python `3.11+` and should remain friendly to normal local development on macOS, Linux, and Windows.

## Bootstrap

```bash
python -m pip install -e '.[dev,docs]'
python scripts/check_repo_integrity.py
pytest
sphinx-build -b html docs/site docs/_build/html
```

This default local-development install includes the CIF-backed structure-import path used by the
normal unit suite, so a clean `pytest` run does not depend on a second undocumented install step.

For a fuller user-facing setup guide, including Windows activation details, notebook use, and PDF build notes, see `../site/tutorials/installation_and_build.md`.

## Working Expectations

- keep the core package importable without heavy optional scientific dependencies
- use optional extras for adapters
- keep generated artifacts out of source directories
- treat documentation and figures as repo-tracked assets, not temporary outputs
