# Testing Strategy

PyTex treats automated testing, validation ledgers, benchmark placeholders, and documentation integrity as one scientific quality system.

## Layers

- unit tests for invariants, conversions, and algorithmic behavior
- parity tests for MTEX-backed texture and EBSD categories
- integration tests for optional adapters and command-line entry points
- documentation policy tests for references and documentation architecture rules
- benchmark fixtures and manifests for comparison-oriented workflows

## Stable Feature Exit Criteria

A stable feature is not complete without:

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

## External-Baseline Policy Beyond MTEX

MTEX is not the only validation authority PyTex needs.

- Texture and EBSD:
  MTEX remains the validation floor where categories overlap.
- Structure import and CIF semantics:
  Validation must be anchored to IUCr and International Tables semantics, with parser behavior checked against documented structure-library expectations.
- Diffraction geometry and kinematic workflows:
  Validation must be tracked in `diffraction_validation_matrix.md` through literature-backed checks, geometry invariants, and future adapter or external-tool comparisons.
- Future phase-transformation workflows:
  Validation must be anchored to literature-backed orientation-relationship and variant-generation references, not only to tool parity.

## Current Review Note

As of the current hardening pass, the repository has passing tests and docs builds, but validation breadth is still uneven:

- texture and EBSD have explicit MTEX-backed ledgers
- diffraction has foundational implementation and internal tests, but its external-baseline program is still emerging
- structure import is implemented, but it needs a more explicit validation story than parser-success tests alone

## References

### Normative

- `mtex_parity_matrix.md`
- `diffraction_validation_matrix.md`
- `../standards/reference_canon.md`

### Informative

- `../tex/validation/validation_program.tex`
