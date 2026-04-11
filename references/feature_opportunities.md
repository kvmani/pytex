# PyTex Feature Opportunities From The Reference Corpus

This document compares the current PyTex surface against the authoritative material in `references/` and records feature ideas that are justified by those sources.

## Current Strengths

PyTex already has meaningful coverage in the following areas:

- Reciprocal-basis, plane-normal, and d-spacing semantics via `src/pytex/core/lattice.py` and `src/pytex/core/miller.py`.
- Hexagonal 3-index/4-index conversions in `src/pytex/core/hexagonal.py` and vectorized variants in `src/pytex/core/miller.py`.
- Orientation representations and conversions in `src/pytex/core/orientation.py`.
- Texture-facing plotting and IPF color key support in `src/pytex/plotting/ipf.py`.
- Pedagogical powder-XRD and SAED generation in `src/pytex/diffraction/xrd.py` and `src/pytex/diffraction/saed.py`.
- EBSD crystal-map, grain, and KAM workflows in `src/pytex/ebsd/`.

The main gap is no longer the total absence of core crystallographic logic. The gap is that many features are not yet tied tightly enough to the reference corpus through examples, figure geometry, and human-auditable validation docs.

## High-Value Opportunities

### 1. First-Class Miller-Bravais Semantics

Reference basis:

- `hexagnoal 4index mathematics.pdf`, pp. 1-7
- `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, Appendix 1

Current PyTex status:

- Scalar and vectorized conversion helpers exist.
- Public semantic types are still primarily 3-index first.

Opportunity:

- Add explicit `MillerBravaisPlane`, `MillerBravaisDirection`, and 4-index zone-law helpers.
- Preserve both 3-index and 4-index notation in docs and formatting utilities.
- Expose round-trip-safe formatters for HCP teaching and validation surfaces.

Why it matters:

- Hexagonal notation is already a stable repo concern.
- The reference corpus provides exact formulas, constraints, and worked examples.

### 2. Metric-Tensor And Direct/Reciprocal Transform APIs

Reference basis:

- `crystallographY_calcualtions.pdf`, book pp. 10-20
- `Kelly & Groves.pdf`, Appendix 1
- `Bhadesia_crystallography.pdf`, pp. 2-6

Current PyTex status:

- Reciprocal basis and d-spacing exist.
- Explicit metric tensors and named component transforms are not first-class public tools.

Opportunity:

- Add public helpers for direct metric tensor, reciprocal metric tensor, and direct<->reciprocal component transforms.
- Surface these in docs and notebooks with cited worked examples.

Why it matters:

- It removes hidden matrix algebra from downstream algorithms.
- It makes reciprocal-space reasoning auditable by domain experts.

### 3. Gnomonic Projection And Kikuchi Geometry Utilities

Reference basis:

- `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 150-159
- `kikuchi maps of cubic and hexgonal crystals.pdf`, Appendix 2
- `williamsandcarter.pdf`, pp. 317-318

Current PyTex status:

- SAED spots exist.
- No explicit Kikuchi-band generator, gnomonic-projection API, or band-index overlay utility exists.

Opportunity:

- Add a `KikuchiBand`, `KikuchiPattern`, or related geometry surface.
- Implement gnomonic projection helpers and pattern-center-aware line geometry.
- Add overlay builders for cubic and hexagonal Kikuchi maps.

Why it matters:

- This would connect EBSD and TEM pedagogy to explicit crystallographic geometry.
- The repo already values publication-quality geometry figures.

### 4. Reference-Backed EBSD Calibration Helpers

Reference basis:

- `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 165-185
- `williamsandcarter.pdf`, pp. 192-195

Current PyTex status:

- Acquisition and calibration metadata objects exist.
- No explicit EBSD pattern-center or camera-calibration workflow API exists.

Opportunity:

- Add calibration records and helper routines for EBSD detector geometry, pattern center, and scan/frame registration.
- Tie them to acquisition manifests and human-readable validation docs.

Why it matters:

- This is a natural extension of the explicit frame model already present in PyTex.

### 5. IPF And Pole-Figure Geometry Beyond The Current Basic Sector View

Reference basis:

- `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 43-50
- `Kelly & Groves.pdf`, Appendix 2

Current PyTex status:

- IPF plotting exists.
- Canonical SVGs and docs still under-explain how the spherical sector, stereographic projection, and symmetry reduction fit together.

Opportunity:

- Add better PF/IPF concept pages and notebook sections using source-backed geometry.
- Extend IPF plotting docs to low-symmetry and HCP cases.
- Add source-backed SVG redraws for spherical-sector reduction and inverse-pole interpretation.

Why it matters:

- The current code is stronger than the documentation surface suggests.

### 6. Explicit Misorientation Spaces And Fundamental Regions

Reference basis:

- `MathsOfrotations_RolletDegraef.pdf`
- `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`, pp. 49-50 and Part 3 Chapter 9

Current PyTex status:

- Misorientation and symmetry-aware orientation objects exist.
- Rodrigues-space and angle-axis-space docs/examples are still thin.

Opportunity:

- Add public conversions and plots for Rodrigues, homochoric, and angle-axis spaces where justified.
- Add symmetry-reduced misorientation-region utilities and cited worked examples.

Why it matters:

- The reference corpus treats these spaces as central texture-analysis tools.

### 7. Reference-Backed Orientation-Relationship Solvers

Reference basis:

- `Bhadesia_crystallography.pdf`, pp. 12-24

Current PyTex status:

- `Orientation.from_plane_direction` exists.
- Variant-family and orientation-relationship workflows are still early.

Opportunity:

- Add named orientation-relationship builders with cited plane-direction inputs.
- Support variant enumeration and interface-plane reasoning for transformation workflows.

Why it matters:

- This would connect core orientation math to transformation science, not just texture display.

### 8. Diffraction Indexing And Zone-Axis Workflows

Reference basis:

- `williamsandcarter.pdf`, pp. 78-79, 317-318
- `crystallographY_calcualtions.pdf`, book pp. 10-20

Current PyTex status:

- XRD and SAED generation exist.
- General indexing helpers and zone-axis inference workflows are limited.

Opportunity:

- Add explicit spot-pattern indexing helpers, zone-axis estimators, and indexed-report objects.
- Use reciprocal-space formulas and Bragg/Laue reasoning already documented in the references.

Why it matters:

- It would move the diffraction subsystem from generation-only toward interpretation.

### 9. Human-Auditable Automated Test Documentation

Reference basis:

- All of the above, especially Donnay, De Graef, Rowenhorst, and Randle/Engler

Current PyTex status:

- Strong automated tests exist.
- Domain experts still have to read code or infer the mathematical basis indirectly.

Opportunity:

- Maintain Sphinx pages that pair each important tested method with formula, source, example, expected output, and last verified code output.

Why it matters:

- This closes the gap between domain review and code review.
- It matches the repo's teaching-and-research dual mandate.

### 10. Figure Regeneration Tooling For Canonical SVGs

Reference basis:

- `Kelly & Groves.pdf`, Appendix 2
- `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf`
- `kikuchi maps of cubic and hexgonal crystals.pdf`

Current PyTex status:

- Canonical SVGs exist, but several are hand-maintained and only loosely tied to explicit page-level references.

Opportunity:

- Add figure source metadata or lightweight generation scripts for canonical geometry figures.
- Keep a source note in the adjacent docs and use exact symbols from the reference corpus.

Why it matters:

- It reduces drift between code, figures, and theory notes.

## Priority Order

1. Human-auditable automated test documentation.
2. First-class Miller-Bravais semantics.
3. Metric-tensor and direct/reciprocal transform APIs.
4. Gnomonic/Kikuchi geometry utilities.
5. Reference-backed EBSD calibration helpers.
6. Orientation-relationship and variant solvers.
7. Diffraction indexing workflows.
8. Figure regeneration tooling.

## References

### Normative

- [Reference Canon](../docs/standards/reference_canon.md)
- [Development Principles](../docs/standards/development_principles.md)

### Informative

- [Reference Index](reference_index.md)
- [Formulation Summary](formulation_summary.md)
