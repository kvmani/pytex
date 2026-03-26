# Testing Strategy

PyTex uses a layered validation program.

## Layers

1. Unit tests
   Verify invariants of canonical primitives, transforms, normalizations, and container validation.

2. Numerical regression tests
   Verify stable outputs for representative scientific cases.

3. MTEX parity tests
   Mirror relevant MTEX public test categories and examples where functionality overlaps.

4. Interoperability tests
   Verify that optional adapters reconcile external-tool semantics into canonical PyTex primitives correctly.

5. Documentation and asset tests
   Verify that required docs, LaTeX sources, and SVG figures exist and remain synchronized with repository structure.

6. Convention tests
   Verify notation conversions, canonicalization rules, and special-system indexing behavior such as HCP 3-index and 4-index mappings.

## Stable Feature Exit Criteria

A stable feature is incomplete until all of the following exist:

- theory documentation
- implementation documentation
- validation note
- benchmark fixture or explicit placeholder
- deterministic automated tests
- example workflow

## Scientific Validation Policy

- MTEX is the floor for relevant validation.
- Parity is recorded in `mtex_parity_matrix.md`.
- The parity ledger may use `foundational` when PyTex has a correct base implementation but not yet full behavioral parity.
- Non-applicability requires explanation, not omission.
- PyTex must also cover cases MTEX does not, especially provenance, interoperability, and explicit convention handling.
- Numerical tolerances and benchmark governance are centralized in `../standards/benchmark_and_tolerance_governance.md`.
