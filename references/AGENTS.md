# References Folder Guide

This file defines local instructions for work that relies on the PDF corpus under `references/`.

## Purpose

The PDFs in this folder are the repository's working corpus for:

- crystallographic formulas and derivations
- canonical worked examples
- orientation, texture, EBSD, and diffraction notation
- figure redraws for orientation, IPF, pole-figure, reciprocal-space, and frame-geometry topics
- human-auditable automated test documentation

This corpus does not override the repository-wide source hierarchy in `../docs/standards/reference_canon.md`.
If a higher-ranked IUCr or International Tables source conflicts with a textbook or paper here, stop and reconcile the conflict explicitly.

## Mandatory Reading Order

Before using any PDF in this folder for implementation, docs, tests, or figures:

1. Read `formulation_summary.md`.
2. Read `reference_index.md`.
3. Read the relevant PDF pages listed there.
4. Compare the source against the current PyTex code, docs, figures, and tests before changing anything.

Do not jump straight into a PDF search and start coding from memory.

## Topic Routing

Use the following sources first for each topic:

- rotations, Euler/quaternion/Rodrigues conventions:
  `MathsOfrotations_RolletDegraef.pdf`
- texture descriptors, pole figures, inverse pole figures, Euler space, Rodrigues space, EBSD/orientation mapping:
  `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`
- hexagonal 4-index and zone-law mathematics:
  `hexagnoal 4index mathematics.pdf`
- reciprocal lattice, metric tensor, interplanar spacing, direct/reciprocal transforms:
  `crystallographY_calcualtions.pdf`, `Kelly & Groves.pdf`, `Bhadesia_crystallography.pdf`
- Kikuchi geometry and map-based indexing references:
  `kikuchi maps of cubic and hexgonal crystals.pdf`, `williamsandcarter.pdf`

## Extraction Rules

- Preserve the source notation as closely as practical when recording formulas.
- Always state the mapping direction, frame meaning, and basis meaning.
- When a source uses a different symbol convention from PyTex, record both:
  source notation first, then the PyTex crosswalk.
- Record exact page anchors in `reference_index.md`.
- Prefer printed page numbers when the PDF exposes the book pagination clearly.
- If only PDF page numbers are reliable, say so explicitly.
- Do not copy long prose from the PDFs; summarize in PyTex's own words.

## Implementation Rules

- New code should cite the exact source pages in the nearest theory, validation, or test documentation.
- If a formula already exists in `formulation_summary.md`, future tasks must consult that document before introducing a new implementation or changing notation.
- If current code disagrees with the cited source, document the discrepancy before changing behavior.
- If a formula is implemented, add or update deterministic tests with source-backed examples.
- If a figure encodes a convention or geometry, redraw the SVG to match the cited reference more closely rather than leaving an approximate sketch in place.

## Documentation Rules

- Major scientific docs should keep symbols, axis names, and index ordering close to the cited source material in this folder.
- The Sphinx validation surface must include human-auditable test documentation:
  formula, source, worked example, expected output, and last verified code output.
- For major formulas and convention-heavy pathways, the validation page should also include
  a domain-audit table with:
  source-derived expected value, current code output, and interpretation.
- If a current result differs from a textbook value because of transpose, active/passive
  choice, frame mapping direction, basis choice, or another convention issue, state that
  distinction explicitly instead of letting readers infer it.
- When adding a new formula summary or feature opportunity, update `reference_index.md` if the page anchor is useful for future discovery.

## Expected Outputs In This Folder

- `reference_index.md`
  Searchable map of important topics to PDFs and pages.
- `formulation_summary.md`
  Ready-to-use formulas and notation crosswalks for implementation and docs.
- `feature_opportunities.md`
  PyTex gap analysis and feature ideas grounded in the reference corpus.

## References

### Normative

- [Reference Canon](../docs/standards/reference_canon.md)
- [Notation And Conventions](../docs/standards/notation_and_conventions.md)
- [Terminology And Symbol Registry](../docs/standards/terminology_and_symbol_registry.md)
- [Testing Strategy](../docs/testing/strategy.md)

### Informative

- [Automated Test Cases](../docs/testing/automated_test_cases.md)
