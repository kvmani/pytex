# Implementation Roadmap

The roadmap is now expressed as capability ladders rather than only chronological phases.

## Capability Ladders

### Semantic Foundation

| Subsystem | Status | Notes |
| --- | --- | --- |
| Documentation and standards | implemented | Mission, standards, Sphinx, LaTeX, SVG, and validation doctrine are in place. |
| Canonical core model | implemented | Frames, transforms, symmetry, lattice, provenance, orientations, and reciprocal primitives exist. |
| Multimodal characterization doctrine | foundational | Shared cross-modality semantics are now documented, but the corresponding stable API family is not fully implemented. |
| Phase-transformation doctrine | planned | Foundation docs exist, but stable transformation primitives are still ahead. |

### Validated Foundational Implementation

| Subsystem | Status | Notes |
| --- | --- | --- |
| Orientation and texture | implemented | Rotations, misorientation, PF/IPF, symmetry reduction, and discrete ODF foundations exist with tests and MTEX-backed ledgers. |
| EBSD regular-grid workflows | implemented | KAM, segmentation, GROD, boundaries, cleanup, graph aggregation, and manifest-backed normalization exist. |
| Diffraction foundations | foundational | Geometry and kinematic spot workflows exist with unit coverage, but the external-baseline validation program is still growing. |
| CIF and structure import | foundational | Core-model construction exists, but broader validation and space-group semantics remain ahead. |

### Research-Grade Algorithmic Expansion

| Subsystem | Status | Notes |
| --- | --- | --- |
| Exact orientation-space boundary catalogs | planned | Required for broader class-by-class parity claims. |
| Harmonic ODF inversion and richer reconstruction | planned | Discrete/kernel foundations exist; harmonic inversion remains ahead. |
| Rich diffraction refinement and intensity models | planned | Current implementation is geometric and kinematic, not full physical modeling. |
| Phase transformation and parent reconstruction | planned | Requires stable transformation semantics first. |

### Teaching-Grade Explanatory Surface

| Subsystem | Status | Notes |
| --- | --- | --- |
| Sphinx concepts and workflows | implemented | Public entry point is live and buildable. |
| Canonical LaTeX theory notes | implemented | Major foundation notes exist and are cross-linked from the site. |
| SVG geometry figures | implemented | Core orientation, diffraction, and EBSD figures exist. |
| Multimodal and transformation teaching notes | foundational | Architectural prose is now defined; broader workflow coverage remains ahead. |

## Immediate Next Steps

1. Keep engineering hygiene green in CI: lint, types, tests, and docs build.
2. Harden space-group and structure semantics before broadening stable structure-facing APIs.
3. Add acquisition, calibration, and uncertainty foundations for multimodal workflows.
4. Grow diffraction validation into a literature-backed and benchmark-backed program.
5. Define and then implement transformation primitives and validation posture.

## References

### Normative

- `../standards/engineering_governance.md`
- `../standards/reference_canon.md`

### Informative

- `../architecture/repo_review_2026_foundation_audit.md`
