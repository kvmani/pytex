# MTEX Parity Fixtures

This tree contains pinned parity fixtures derived from MTEX-documented conventions and workflows.

Each fixture bundle records:

- the intended MTEX behavior category
- the pinned MTEX version target
- the source script name reserved for future regeneration
- numerical tolerances used by the parity tests

These fixtures are intentionally small and human-auditable. They are the deterministic oracle for the parity tests under `tests/parity/`.

## Campaign-Based Validation

The newer campaign files under `campaigns/` define shared input JSON that both MTEX and PyTex can
read. Each ecosystem writes `pytex.parity_result` JSON with identical field names so outputs can be
compared field by field.

The active campaigns currently cover:

- orientation construction and operations
- IPF color calculations
- Miller geometry and symmetry families
- sampled discrete ODF density evaluation

The XRDML pole-figure and XRDML ODF campaigns are present but pending until cubic and hexagonal
fixture files are provided.
