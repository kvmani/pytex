# Benchmark And Tolerance Governance

## Purpose

Scientific software backtracks expensively when tolerances, fixtures, and comparison policies are inconsistent. PyTex therefore fixes these rules centrally.

## Benchmark Rules

- Every benchmark input must have provenance.
- Every benchmark must record whether it is:
  - parity against external software
  - literature-backed numerical validation
  - PyTex-only regression coverage
- Every benchmark result note must state:
  - expected behavior
  - numerical tolerances
  - reasons for any accepted discrepancy

## Tolerance Rules

- Tolerances must be named and justified in tests when they are not obvious machine-precision checks.
- Rotations, transforms, and orthonormality checks should default to tight floating-point tolerances.
- Geometry and indexing interoperability tests may require looser tolerances, but those must be documented.
- Visual regression tests should document whether they are pixel, vector-asset, or semantic checks.

## Determinism Rules

- Randomized workflows must expose a deterministic seed path.
- Regression tests should prefer fixed fixtures over ad hoc generated inputs unless the generation process is itself under test.
