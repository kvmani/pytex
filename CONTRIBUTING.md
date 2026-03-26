# Contributing

PyTex is being built as a scientific software program, not as an ungoverned prototype.

Before making substantial changes, read:

- `mission.md`
- `specifications.md`
- `AGENTS.md`
- `docs/README.md`
- `docs/standards/engineering_governance.md`
- `docs/standards/notation_and_conventions.md`
- `docs/testing/strategy.md`

## Working Rules

- Stable public APIs must use canonical PyTex domain primitives where frame or symmetry meaning matters.
- Significant scientific behavior changes require documentation updates in the same change.
- Stable features are incomplete until theory, implementation notes, validation notes, and figures exist.
- New numerical behavior requires tests; bug fixes require regression tests.
- LaTeX documents and SVG figures are product assets and must be treated with the same care as source code.
- Do not introduce a new frame or symmetry convention locally in a subsystem. Extend the shared core model instead.

## Local Development

```bash
python -m pip install -e '.[dev]'
python scripts/check_repo_integrity.py
pytest
ruff check .
mypy src
```

## Pull Request Expectations

- explain the scientific or architectural motivation
- list documentation updates
- list validation or tests added
- call out any remaining non-applicability items in the MTEX parity matrix when relevant
