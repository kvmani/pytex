# Phase Transformation Validation Matrix

This document is the authoritative validation ledger for PyTex phase-transformation workflows.

## Status Keys

- `implemented`: automated coverage and validation notes exist for the current category
- `foundational`: the implementation exists and is scientifically structured, but the external or
  literature-backed validation surface is not yet complete
- `planned`: the category is accepted but not yet validated adequately
- `n/a`: not applicable to current PyTex scope, with explanation

## Matrix

| Area | Baseline | Status | Notes |
| --- | --- | --- | --- |
| Orientation-relationship semantic contracts | Core invariant tests and canonical data-model doctrine | implemented | Parent and child phases, frames, and provenance are checked through automated tests. |
| Orientation-relationship construction from explicit plane-direction correspondence | Deterministic right-handed basis construction plus unit tests | implemented | `OrientationRelationship.from_parallel_plane_direction(...)` is covered for matrix recovery and phase mismatch rejection. |
| Named Bain correspondence helper | Deterministic correspondence construction plus cubic-phase guards | implemented | `OrientationRelationship.from_bain_correspondence(...)` is covered for the stated `(001)_p || (001)_c`, `[110]_p || [100]_c` mapping and rejects non-cubic parents. |
| Named Nishiyama-Wassermann helper | Deterministic correspondence construction plus cubic-phase guards | implemented | `OrientationRelationship.from_nishiyama_wassermann_correspondence(...)` is covered for the stated `(111)_p || (011)_c`, `[1-10]_p || [100]_c` mapping and rejects non-cubic children. |
| Variant generation and uniqueness | Internal deterministic tests plus symmetry-backed invariants | implemented | `TransformationVariant` generation is covered by unit tests and benchmark manifests. |
| Variant-indexed predicted child orientations | Internal deterministic tests and manifest-backed workflow identity | implemented | `PhaseTransformationRecord.predicted_child_orientations()` now respects explicit variant assignments. |
| Transformation manifest schema | Stable JSON schema plus round-trip tests | implemented | `TransformationManifest` now records dedicated transformation workflow context. |
| Experimental parent-candidate scoring | Internal deterministic tests and benchmark identity | implemented | Candidate-parent ranking is staged under `pytex.experimental` with explicit non-stable status. |
| Literature-tracked starter transformation families | Textbook Bain and Nishiyama-Wassermann correspondences | foundational | PyTex now tracks named literature-facing helpers and tests, but broader family breadth and curated datasets remain ahead. |
| Full parent reconstruction workflows | Literature-backed reconstruction studies and curated datasets | planned | PyTex intentionally keeps full reconstruction outside the stable API until benchmark breadth is stronger. |

## Current Posture

PyTex now has explicit transformation primitives, a dedicated manifest contract, named literature
starter correspondences, and a bounded experimental reconstruction-scoring surface. The semantics
are no longer ad hoc, but the validation posture is still foundational because broad
literature-backed families and curated reconstruction datasets are not yet in place.

## References

### Normative

- `strategy.md`
- `../architecture/phase_transformation_foundation.md`

### Informative

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*
- `../../benchmarks/transformation/variant_prediction_benchmark_manifest.json`
