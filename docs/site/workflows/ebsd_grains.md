# EBSD: Grain Segmentation And GROD

PyTex includes a minimal but scientifically structured Phase 3 grain workflow for regular-grid EBSD maps: thresholded grain segmentation, grain reference orientation deviation (GROD), grain-boundary extraction, and deterministic small-grain cleanup.

## Conceptual Flow

The workflow deliberately reuses the same core semantics as the rest of the library:

1. start from a `CrystalMap` with explicit map-frame and orientation semantics
2. define neighborhood adjacency on the regular grid
3. segment grains by a misorientation threshold
4. choose a representative orientation within each grain
5. compute GROD relative to that representative orientation
6. extract boundary segments and optionally clean small grains

![EBSD Grain Workflow](../../figures/ebsd_grain_workflow.svg)

## Scope

- thresholded connected-component segmentation on 2D regular grids
- 4- and 8-connectivity support
- representative grain orientation selection from measured orientations
- GROD maps relative to that representative orientation
- grain-boundary extraction from cross-grain neighbor pairs
- per-segment boundary length from `step_sizes`
- low-angle or high-angle classification through a user threshold
- adjacency-based small-grain merging with deterministic tie-breaking
- grain-graph aggregation from the extracted boundary network
- majority-vote label smoothing for regular-grid denoising

## GROD In This Foundation Layer

PyTex currently defines GROD relative to a representative measured orientation for each grain. That choice is conservative: it avoids exposing a public “grain mean orientation” surface until the averaging rule and its symmetry behavior are documented tightly enough for stable parity testing.

```{note}
This is a deliberate API boundary. PyTex prefers a narrower but explicit scientific contract over exposing a broader workflow surface with underspecified averaging semantics.
```

## Example

```python
segmentation = crystal_map.segment_grains(
    max_misorientation_deg=5.0,
    symmetry_aware=True,
    connectivity=4,
)

labels = segmentation.label_grid
grod = segmentation.grod_map_deg()
boundaries = segmentation.boundary_network(high_angle_threshold_deg=15.0)
cleaned = segmentation.merge_small_grains(
    min_size=4,
    until_stable=True,
)
```

## How To Read The Outputs

### Grain Labels

The label grid is the most direct segmentation surface. It is useful for checking threshold sensitivity, connectivity behavior, and cleanup effects.

### GROD

GROD is a within-grain misorientation metric. It is useful for detecting internal grain deformation structure, but its meaning depends on the chosen grain reference orientation.

### Boundary Network

Boundary segments represent cross-grain adjacency on the map. Length is derived from regular-grid spacing, and classification can be split into low-angle versus high-angle boundaries through a threshold.

### Grain Graph

PyTex can now aggregate the boundary network into a grain graph. Each graph edge records:

- the adjacent grain pair
- the contributing boundary-segment indices
- total shared boundary length
- mean boundary misorientation
- high-angle fraction under the chosen threshold

### Cleanup

Small-grain cleanup is intended as a controlled post-processing step, not an excuse to ignore segmentation sensitivity. PyTex keeps the merge rule deterministic so that parity tests and user workflows stay reproducible.

### Label Smoothing

PyTex now also exposes majority-vote label smoothing through `GrainSegmentation.majority_smoothed()`. This is intentionally a label-space denoising operation:

- it acts on an existing segmentation rather than on raw orientations
- it uses neighbor label votes on the regular grid
- it is deterministic and iteration-controlled

This makes it useful for removing isolated segmentation noise while keeping the denoising policy explicit.

## Current Limits

- no alternate grain-reference definitions beyond the current representative measured orientation
- vendor-backed object bridges now exist for normalization, but detector/pattern metadata normalization is still ahead

## Related Material

- `docs/architecture/ebsd_foundation.md`
- [../../tex/algorithms/ebsd_grain_segmentation_and_grod.tex](../../tex/algorithms/ebsd_grain_segmentation_and_grod.tex)
- [../../tex/algorithms/ebsd_boundaries_and_cleanup.tex](../../tex/algorithms/ebsd_boundaries_and_cleanup.tex)
- [../../figures/ebsd_grain_workflow.svg](../../figures/ebsd_grain_workflow.svg)

## References

### Informative

- MTEX documentation: [Grain Reference Orientation Deviation (GROD)](https://mtex-toolbox.github.io/EBSDGROD.html)
- MTEX documentation: [Grain Orientation Parameters](https://mtex-toolbox.github.io/GrainOrientationParameters.html)
- MTEX documentation: [Selecting Grains](https://mtex-toolbox.github.io/SelectingGrains.html)
