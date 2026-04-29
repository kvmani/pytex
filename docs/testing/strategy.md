# Testing Strategy

PyTex treats automated testing, validation ledgers, stable benchmark and validation manifests, and
documentation integrity as one scientific quality system.

## Layers

- unit tests for invariants, conversions, and algorithmic behavior
- parity tests for MTEX-backed texture and EBSD categories
- integration tests for optional adapters and command-line entry points
- documentation policy tests for references and documentation architecture rules
- benchmark fixtures and manifests for comparison-oriented workflows

## Execution Lanes

PyTex now treats the test and quality surface as two explicit lanes:

- `base lane`: `.[dev,docs]` for integrity, linting, typing, docs builds, and the default lightweight test suite
- `full scientific lane`: `.[dev,docs,adapters]` for CIF-backed structure import, pinned diffraction external baselines, and adapter-heavy interoperability coverage

The base lane is the default contributor environment. The full scientific lane is the controlling
environment for claims that depend on optional scientific packages such as `pymatgen`, ORIX,
KikuchiPy, PyEBSDIndex, or diffsims.

## Stable Feature Exit Criteria

A stable feature is not complete without:

- theory documentation
- implementation documentation
- validation note
- benchmark fixture or explicit placeholder
- deterministic automated tests
- human-auditable automated test documentation for major formulas, conventions, or workflows
- example workflow

For pinned phase fixtures, deterministic repository hygiene also includes the catalog digests in
`fixtures/phases/catalog.json`. When those hashes must be refreshed, use
`python scripts/regenerate_phase_fixture_catalog_hashes.py`; the script rewrites only the digest
fields in catalog order with stable JSON formatting so `python scripts/check_repo_integrity.py`
remains reproducible across contributors and CI.

## Scientific Validation Policy

- MTEX is the floor for relevant validation.
- Parity is recorded in `mtex_parity_matrix.md`.
- The parity ledger may use `foundational` when PyTex has a correct base implementation but not yet full behavioral parity.
- Non-applicability requires explanation, not omission.
- PyTex must also cover cases MTEX does not, especially provenance, interoperability, and explicit convention handling.
- Numerical tolerances and benchmark governance are centralized in `../standards/benchmark_and_tolerance_governance.md`.

## Human-Auditable Automated Test Documentation

Important tests should be readable by domain experts who do not want to reverse-engineer the
implementation first.

- The Sphinx validation surface must include documented test cases for major formulas, conversions, and conventions.
- Each documented test case should identify the code surface, the authoritative reference, the governing formula, at least one worked example, the expected automated assertion, and the last verified code output.
- For formula-heavy pathways, the documented case should include a domain-audit table with at least: quantity, source-derived expected value, current code output, and interpretation.
- When reference notation and PyTex notation differ, the documented case should state the mismatch explicitly instead of silently rewriting the source into PyTex terminology.
- When current code output differs from the most obvious textbook expression only by mapping direction, transpose, basis choice, or another convention issue, the documentation should explain that distinction directly.
- When a future task changes a cited numerical pathway, the corresponding test documentation should be updated in the same change set.
- The working reference corpus for such documentation is indexed in `references/reference_index.md`, and future tasks should consult `references/formulation_summary.md` before introducing new formulations or examples.

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
  Validation must be tracked in `plotting_validation_matrix.md` through semantic input checks, deterministic builder coverage, structural figure-property assertions, and future parity-oriented comparisons when scientifically justified.
- XRDML texture import:
  Validation must combine at least one pinned real vendor-style XML fixture with deterministic synthetic inversion cases so import semantics and reconstruction semantics can be audited separately.
- Open legacy pole-figure formats:
  Validation may use public multi-material `PPF` or `EPF` exercise datasets through external acquisition scripts when redistribution terms are unclear, but the provenance, format semantics, and non-bundled status must be documented explicitly.
- Future phase-transformation workflows:
  Validation must be anchored to literature-backed orientation-relationship and variant-generation references, not only to tool parity.

## Current Review Note

The repository now has passing integrity checks, tests, docs builds, and a passing full
scientific lane for the immediate roadmap surface, but validation breadth is still uneven:

- texture and EBSD have explicit MTEX-backed ledgers
- texture now includes both explicit dictionary inversion and a first harmonic PF-to-ODF reconstruction path, but higher-order external parity and correction workflows remain ahead
- structure import now has a hash-pinned fixture corpus, manifest-backed audit workflow, and default test-suite coverage, but broader IUCr-style external baselines are still ahead
- diffraction now has pinned open-source external baselines for one powder XRD workflow family and one SAED workflow family across `ni_fcc` and `fe_bcc`, but broader material and orientation coverage remains ahead
- plotting now has structural validation coverage for the highest-value runtime surfaces, but external visual-parity work remains ahead
- XRDML import now has a pinned open-source regression fixture plus synthetic ODF-inversion coverage, but broader vendor and correction-workflow coverage remains ahead
- phase transformation now has a dedicated manifest schema, benchmark identity, validation ledger, and experimental parent-candidate scoring surface, but literature breadth and full reconstruction workflows remain ahead
- the priority teaching notebooks now smoke-execute in the default suite, but the full notebook atlas still follows a lighter validation path than the primary roadmap sequence

The main quality-gap question is therefore no longer whether the repo can build and test at all, or
whether the optional `pymatgen`-gated paths execute end to end. The remaining question is how far
the current evidence justifies broader scientific claims for structure import, diffraction,
transformation, and interoperability.

## References

### Normative

- <a href="mtex_parity_matrix.html">MTEX Parity Matrix</a>
- <a href="diffraction_validation_matrix.html">Diffraction Validation Matrix</a>
- <a href="structure_validation_matrix.html">Structure Validation Matrix</a>
- <a href="plotting_validation_matrix.html">Plotting Validation Matrix</a>
- <a href="phase_transformation_validation_matrix.html">Phase Transformation Validation Matrix</a>
- <a href="automated_test_cases.html">Automated Test Cases</a>
- <a href="../standards/reference_canon.html">Reference Canon</a>

### Informative

- <a href="../tex/validation/validation_program.tex">Validation Program</a>
