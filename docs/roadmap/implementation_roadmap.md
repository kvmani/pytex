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
| Diffraction foundations | foundational | Geometry and kinematic spot workflows exist with unit coverage, but the external-baseline validation program is still growing. |
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

## Immediate Next Steps

1. Keep engineering hygiene green in CI: lint, types, tests, and docs build.
2. Grow structure-import validation into a literature-backed and benchmark-backed program.
3. Grow the new manifest family into richer experiment, benchmark, validation, and workflow-result use cases.
4. Grow diffraction validation into a literature-backed and benchmark-backed program.
5. Extend the transformation primitive family into literature-backed variant-generation and parent-reconstruction workflows.

## References

### Normative

- `../standards/engineering_governance.md`
- `../standards/reference_canon.md`

### Informative

- `../architecture/repo_review_2026_foundation_audit.md`
