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
| Variant generation and uniqueness | Internal deterministic tests plus symmetry-backed invariants | implemented | `TransformationVariant` generation is covered by unit tests and benchmark manifests. |
| Variant-indexed predicted child orientations | Internal deterministic tests and manifest-backed workflow identity | implemented | `PhaseTransformationRecord.predicted_child_orientations()` now respects explicit variant assignments. |
| Transformation manifest schema | Stable JSON schema plus round-trip tests | implemented | `TransformationManifest` now records dedicated transformation workflow context. |
| Experimental parent-candidate scoring | Internal deterministic tests and benchmark identity | implemented | Candidate-parent ranking is staged under `pytex.experimental` with explicit non-stable status. |
| Literature-backed transformation-family parity | Textbook and peer-reviewed orientation-relationship references | foundational | The first validation ledger exists, but broad literature-backed case coverage remains ahead. |
| Full parent reconstruction workflows | Literature-backed reconstruction studies and curated datasets | planned | PyTex intentionally keeps full reconstruction outside the stable API until benchmark breadth is stronger. |

## Current Posture

PyTex now has explicit transformation primitives, a dedicated manifest contract, and a bounded
experimental reconstruction-scoring surface. The semantics are no longer ad hoc, but the
validation posture is still foundational because broad literature-backed families and curated
reconstruction datasets are not yet in place.

## References

### Normative

- `strategy.md`
- `../architecture/phase_transformation_foundation.md`

### Informative

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*
- `../../benchmarks/transformation/variant_prediction_benchmark_manifest.json`
