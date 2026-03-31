# Phase Transformation Manifests And Experimental Scoring

PyTex now has a dedicated transformation manifest contract plus an explicitly experimental
parent-candidate scoring surface.

## Stable Transformation Surfaces

- `OrientationRelationship`
- `TransformationVariant`
- `PhaseTransformationRecord`
- `TransformationManifest`

## Experimental Surface

- `pytex.experimental.score_parent_orientations(...)`
- `pytex.experimental.ParentReconstructionResult`

## Manifest Example

```python
from pytex import TransformationManifest

manifest = TransformationManifest.from_phase_transformation_record(
    record,
    referenced_files=("indexed_map.h5",),
    notes=("Variant assignments came from a bounded indexing workflow.",),
)
```

The dedicated transformation manifest keeps the transformation identifier, parent and child phase
declarations, variant-count metadata, and referenced files in one schema-backed artifact.

## Experimental Candidate Scoring

```python
from pytex.experimental import score_parent_orientations

result = score_parent_orientations(
    record,
    candidate_parent_orientations,
    symmetry_aware=True,
    reduction="mean",
)

best_parent = result.best_parent_orientation()
predicted_children = result.predicted_child_orientations()
```

This is intentionally not a full parent-reconstruction claim. It is a bounded ranking primitive for
research workflows that already have a transformation record and a set of candidate parents.

## Why This Split Matters

- transformation semantics and provenance are stable
- broad reconstruction algorithms and literature coverage are still ahead
- researchers can still run bounded candidate-ranking studies without pretending the whole workflow
  is already stable

## Related Material

- {doc}`../architecture/phase_transformation_foundation`
- {doc}`../validation/phase_transformation_validation_matrix`
- [../../tex/algorithms/experimental_parent_candidate_scoring.tex](../../tex/algorithms/experimental_parent_candidate_scoring.tex)
