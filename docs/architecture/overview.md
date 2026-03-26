# Architecture Overview

PyTex is organized as a modular scientific library built around a canonical internal data model.

## Architectural Priorities

1. one shared model for frames, symmetry, provenance, and structure semantics
2. stable APIs that communicate scientific meaning explicitly
3. optional adapter boundaries for heavy external tools
4. documentation and validation artifacts treated as first-class outputs
5. foundational semantics explained through mathematics and graphics, not only code and prose

## Module Map

```text
src/pytex/
  core/
    conventions.py
    provenance.py
    frames.py
    symmetry.py
    lattice.py
    orientation.py
  texture/
    models.py
  ebsd/
    models.py
  diffraction/
    models.py
  adapters/
  plotting/
  experimental/
```

## Current Posture

The repository already contains real foundational implementations in:

- `core/` for frames, symmetry, lattice, provenance, and orientation primitives
- `core/` semantic batch primitives for vectors, Euler angles, quaternions, rotations, and orientations
- `core/` multimodal acquisition primitives for geometry, calibration, measurement quality, and scattering context
- `core/` transformation primitives for orientation relationships, variants, and phase-transformation records
- `texture/` for PF, IPF, ODF, and color-key foundations
- `ebsd/` for regular-grid neighborhood, KAM, segmentation, GROD, boundaries, cleanup, and graph aggregation
- `diffraction/` for geometry, reciprocal-space primitives, kinematic spots, family grouping, and indexing scaffolding
- `adapters/` for stable manifest contracts spanning import, experiment, benchmark, validation, and workflow-result interchange

The current architectural risk is not lack of code. It is adding new stable surfaces before the multimodal semantics are written down centrally.

## Layering Rules

- `core/` owns canonical primitives and low-level math semantics.
- `core/` also owns the canonical batch-semantic layer for vectorized scientific operations.
- `texture/`, `ebsd/`, and `diffraction/` build on `core/`.
- `adapters/` may depend on external libraries, but stable core types must remain usable without them.
- `experimental/` can depend on stable types but must not weaken or bypass stable semantics.

## Data Flow Principle

Imports from vendors or external scientific libraries should be normalized into PyTex canonical primitives at the boundary. Internal algorithms should then work on PyTex types, not on source-specific representations.

## Adjacent Architecture Notes

- `multimodal_characterization_foundation.md` defines the shared semantic backbone across EBSD, XRD, neutron, and TEM.
- `phase_transformation_foundation.md` defines the parent-child and variant contracts that now anchor the stable transformation primitive family.
- `repo_review_2026_foundation_audit.md` records the current state of the repository and the immediate hardening priorities.

## References

### Normative

- `canonical_data_model.md`
- `../standards/reference_canon.md`
- `../standards/engineering_governance.md`

### Informative

- `orientation_and_texture_foundation.md`
- `ebsd_foundation.md`
- `diffraction_foundation.md`
