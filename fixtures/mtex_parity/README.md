# MTEX Parity Fixtures

This tree contains pinned parity fixtures derived from MTEX-documented conventions and workflows.

Each fixture bundle records:

- the intended MTEX behavior category
- the pinned MTEX version target
- the source script name reserved for future regeneration
- numerical tolerances used by the parity tests

These fixtures are intentionally small and human-auditable. They are the deterministic oracle for the parity tests under `tests/parity/`.
