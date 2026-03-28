# Implementation Roadmap

The roadmap is now expressed as capability ladders rather than only chronological phases.

## Capability Ladders

### Semantic Foundation

| Subsystem | Status | Notes |
| --- | --- | --- |
| Documentation and standards | implemented | Mission, standards, Sphinx, LaTeX, SVG, and validation doctrine are in place. |
| Canonical core model | implemented | Frames, transforms, symmetry, lattice, provenance, orientations, and reciprocal primitives exist. |
| Multimodal characterization doctrine | implemented | Shared cross-modality semantics are documented and the core acquisition/calibration/quality primitive family now exists. |
| Phase-transformation doctrine | implemented | Foundation docs exist and the first stable transformation primitive family now exists in the core model. |

### Validated Foundational Implementation

| Subsystem | Status | Notes |
| --- | --- | --- |
| Orientation and texture | implemented | Rotations, misorientation, PF/IPF, symmetry reduction, discrete ODF foundations, and runtime contour or section plotting surfaces exist with tests and MTEX-backed ledgers. |
| EBSD regular-grid workflows | implemented | KAM, segmentation, GROD, boundaries, cleanup, graph aggregation, and manifest-backed normalization exist. |
| Diffraction foundations | foundational | Geometry, powder XRD generation, SAED spot workflows, and diffraction plotting exist with unit coverage, but the external-baseline validation program is still growing. |
| Scientific visualization | foundational | Shared YAML style themes, runtime diffraction plotting, and VESTA-like crystal visualization now exist, but broader publication presets and external visual baselines are still ahead. |
| CIF and structure import | foundational | Core-model construction, space-group semantics, and a dedicated validation ledger exist, but broader external baselines remain ahead. |
| Multimodal acquisition core | foundational | Acquisition geometry, calibration, quality, scattering, experiment manifests, and workflow entry points now exist, but broader modality coverage remains ahead. |
| Phase-transformation foundation | foundational | Core transformation primitives now exist, but literature-backed validation and richer algorithms remain ahead. |

### Research-Grade Algorithmic Expansion

| Subsystem | Status | Notes |
| --- | --- | --- |
| Exact orientation-space boundary catalogs | planned | Required for broader class-by-class parity claims. |
| Harmonic ODF inversion and richer reconstruction | planned | Discrete/kernel foundations exist; harmonic inversion remains ahead. |
| Rich diffraction refinement and intensity models | planned | Current implementation is geometric and kinematic, not full physical modeling. |
| Phase transformation and parent reconstruction | planned | Stable transformation semantics now exist; algorithmic breadth and validation remain ahead. |

### Teaching-Grade Explanatory Surface

| Subsystem | Status | Notes |
| --- | --- | --- |
| Sphinx concepts and workflows | implemented | Public entry point is live and buildable. |
| Canonical LaTeX theory notes | implemented | Major foundation notes exist and are cross-linked from the site. |
| SVG geometry figures | implemented | Core orientation, diffraction, and EBSD figures exist. |
| Multimodal and transformation teaching notes | foundational | Architectural prose is now defined; broader workflow coverage remains ahead. |

## Immediate Execution Roadmap

This section is the active near-term development program for the next few days and weeks.
Unless an explicit decision supersedes it, immediate implementation work should align to this sequence rather than expanding breadth opportunistically.

The governing rule for near-term execution is:

1. strengthen scientific defensibility before adding major new surface area
2. strengthen reproducibility before making stronger public claims
3. strengthen publication and teaching quality alongside runtime implementation

### Phase 1: Validation Corpus And Reproducibility Backbone

Priority: highest

Goals:

- formalize the bundled phase fixtures as the default structure-validation corpus
- make benchmark and validation manifests first-class in real workflows
- tighten repo-integrity checks around scientific assets

Expected implementation focus:

- expand the bundled CIF fixture corpus and fixture metadata contracts
- enrich structure-import benchmark and validation manifests with real fixture references and expected semantics
- ensure every bundled phase fixture is consumed by tests, benchmark manifests, and at least one documentation or workflow surface
- keep repo-integrity checks strict for provenance, metadata completeness, and manifest references

Definition of done for this phase:

- every bundled phase fixture has provenance, citation, redistribution status, and expected symmetry metadata
- structure-import benchmarks run from pinned in-repo assets
- fixture and manifest integrity failures are enforced automatically

### Phase 2: Diffraction External-Baseline Program

Priority: very high

Goals:

- move diffraction from internally coherent to externally benchmarked
- establish pinned validation cases for powder XRD and SAED

Expected implementation focus:

- expand diffraction benchmark and validation manifests
- add canonical benchmark cases derived from the pinned fixture corpus where appropriate
- record expected peak-position and spot-geometry reference data from literature or trusted open datasets
- add comparison tests that clearly separate geometric agreement from physical-intensity limitations

Definition of done for this phase:

- at least one powder XRD case has pinned external baseline data
- at least one SAED case has pinned external baseline data
- diffraction validation status can advance beyond purely foundational internal checks

### Phase 3: Publication-Grade Visualization Program

Priority: high

Goals:

- make plotting and crystal visualization consistently publication-ready
- tighten alignment between scientific semantics and visual presentation

Expected implementation focus:

- add pinned visual regression cases for core plots and crystal scenes
- grow stronger publication and lecture house styles on top of the shared YAML theme system
- extend figure-quality validation around Miller annotation formatting, overlay clipping, and scene composition
- continue auditing canonical SVG figures for scientific and layout correctness

Definition of done for this phase:

- core plots and crystal scenes support stable publication-style output
- visual regressions catch broken labels, overlays, and style drift
- docs, notebooks, and runtime examples use the same higher-quality style presets

### Phase 4: Teaching Surface Consolidation

Priority: high

Goals:

- make notebooks and Sphinx workflows behave as one coherent teaching system
- demonstrate the richer visualization and validation surface through canonical examples

Expected implementation focus:

- upgrade the highest-value notebooks for symmetry, phases/CIF, diffraction, and crystal visualization
- add explicit callouts when a figure or workflow is pedagogical rather than exact
- improve cross-linking between concept pages, workflows, notebooks, theory notes, and validation pages
- standardize publication-quality figure export patterns in tutorial material

Definition of done for this phase:

- the main learning path from canonical structure to visualization to diffraction to manifest-backed workflow is coherent
- major notebooks and pages no longer depend on outdated visuals or terminology

### Phase 5: Algorithmic Expansion Preparation

Priority: after Phases 1-4 are materially advanced

Goals:

- prepare the next research-grade algorithm tier without rushing implementation ahead of validation

Expected implementation focus:

- extend theory and algorithm notes for harmonic ODF inversion, richer diffraction intensity/refinement, parent reconstruction, and exact orientation-space boundary catalogs
- define benchmark placeholders and validation plans before major new algorithmic releases
- choose the next algorithmic feature based on scientific leverage and validation feasibility

Definition of done for this phase:

- the next major algorithmic feature has theory, benchmark shape, and validation plan defined before implementation

## Immediate Default Development Order

Unless explicitly redirected, immediate development should follow this order:

1. engineering hygiene: keep tests, docs build, and repository integrity green
2. structure-import benchmark hardening and fixture-corpus strengthening
3. diffraction external-baseline validation
4. manifest-backed reproducibility expansion
5. publication-grade plotting and visual validation
6. teaching-surface consolidation
7. next algorithmic expansion only after the above are materially stronger

## Near-Term Non-Goals

The current roadmap intentionally deprioritizes the following until the hardening program above is further advanced:

- broad new modality expansion
- major stable adapter-surface expansion
- full physical diffraction intensity modeling before geometric and benchmark validation are stronger
- large new algorithmic breadth without theory, benchmark, and validation scaffolding

## References

### Normative

- `../standards/engineering_governance.md`
- `../standards/reference_canon.md`

### Informative

- `../architecture/repo_review_2026_foundation_audit.md`
