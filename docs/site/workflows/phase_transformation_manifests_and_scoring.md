# Phase Transformation Manifests And Experimental Scoring

PyTex now has a dedicated transformation manifest contract plus an explicitly experimental
parent-candidate scoring surface.

## Stable Transformation Surfaces

- `OrientationRelationship`
- `TransformationVariant`
- `PhaseTransformationRecord`
- `TransformationManifest`

## Constructing A Relationship From Crystallographic Correspondence

When one parent-plane / child-plane correspondence plus one parent-direction / child-direction
correspondence is known, PyTex can construct the parent-to-child crystal mapping directly:

```python
from pytex import CrystalDirection, CrystalPlane, MillerIndex, OrientationRelationship

relationship = OrientationRelationship.from_parallel_plane_direction(
    name="demo_or",
    parent_plane=CrystalPlane(MillerIndex([0, 0, 1], phase=parent_phase), phase=parent_phase),
    child_plane=CrystalPlane(MillerIndex([0, 0, 1], phase=child_phase), phase=child_phase),
    parent_direction=CrystalDirection([1.0, 0.0, 0.0], phase=parent_phase),
    child_direction=CrystalDirection([0.0, 1.0, 0.0], phase=child_phase),
)
```

This constructor keeps the crystallographic evidence visible:

- the planes preserve reciprocal-space meaning on the parent and child phases
- the directions preserve direct-space meaning on the parent and child phases
- the resulting rotation is the right-handed parent-crystal to child-crystal mapping implied by
  those correspondences

## Named Bain Helper

PyTex now exposes one named literature-style helper for cubic parent and child phases:

```python
relationship = OrientationRelationship.from_bain_correspondence(
    parent_phase=parent_phase,
    child_phase=child_phase,
)
```

The current builder encodes this explicit correspondence:

- `(001)_parent || (001)_child`
- `[110]_parent || [100]_child`

This is intentionally narrow. It is meant to make one named relationship reproducible and auditable
before a broader catalog of competing literature families is added.

## Named Nishiyama-Wassermann Helper

PyTex now also exposes a Nishiyama-Wassermann helper for cubic parent and child phases:

```python
relationship = OrientationRelationship.from_nishiyama_wassermann_correspondence(
    parent_phase=parent_phase,
    child_phase=child_phase,
)
```

The current builder encodes this explicit correspondence:

- `(111)_parent || (011)_child`
- `[1-10]_parent || [100]_child`

As with Bain, the point is not to imply a complete literature catalog. The point is to make one
named relationship explicit, reproducible, and test-backed.

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

## Current Limits

- The correspondence constructor currently uses one plane plus one in-plane direction pair as the
  defining data. It does not yet resolve overdetermined or conflicting literature correspondences.
- Only a small named-helper set currently exists. Broader literature orientation-relationship
  catalogs are still ahead.
- Full parent reconstruction remains experimental even though relationship semantics are now
  stronger.

## Related Material

- {doc}`../architecture/phase_transformation_foundation`
- {doc}`../validation/phase_transformation_validation_matrix`
- [../../tex/algorithms/phase_transformation_relationship_construction.tex](../../tex/algorithms/phase_transformation_relationship_construction.tex)
- [../../tex/algorithms/experimental_parent_candidate_scoring.tex](../../tex/algorithms/experimental_parent_candidate_scoring.tex)
