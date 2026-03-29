# Testing Strategy

PyTex treats automated testing, validation ledgers, stable benchmark and validation manifests, and documentation integrity as one scientific quality system.

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
- Powder XRD, SAED, and crystal-visualization workflows:
  Validation must cover geometry and unit invariants, deterministic style handling, and clear statement of which intensity and rendering choices are pedagogical approximations rather than full physical models.
- Stable plotting and texture-inspection workflows:
  Validation must be tracked in `plotting_validation_matrix.md` through semantic input checks, deterministic builder coverage, and future parity-oriented visual baselines.
- Future phase-transformation workflows:
  Validation must be anchored to literature-backed orientation-relationship and variant-generation references, not only to tool parity.

## Current Review Note

As of the current hardening pass, the repository has passing integrity checks, tests, and docs
builds for the immediate roadmap surface, but validation breadth is still uneven:

- texture and EBSD have explicit MTEX-backed ledgers
- structure import now has a hash-pinned fixture corpus, manifest-backed audit workflow, and
  default test-suite coverage, but broader IUCr-style external baselines are still ahead
- diffraction now has pinned open-source external baselines for one powder XRD case and one SAED
  case, but broader material and orientation coverage remains ahead
- plotting now has deterministic SVG visual-regression coverage for the highest-value runtime
  surfaces, but external visual-parity work remains ahead
- the priority teaching notebooks now smoke-execute in the default suite, but the full notebook
  atlas still follows a lighter validation path than the primary roadmap sequence

## References

### Normative

- `mtex_parity_matrix.md`
- `diffraction_validation_matrix.md`
- `structure_validation_matrix.md`
- `plotting_validation_matrix.md`
- `../standards/reference_canon.md`

### Informative

- `../tex/validation/validation_program.tex`
