# EBSD: Kernel-Average Misorientation

Kernel-average misorientation (KAM) is a local orientation-gradient metric on an EBSD map. PyTex
now computes it on a shared neighbor-graph substrate rather than only on a fixed regular-grid
implementation path.

## Definition

For a measurement site $i$ with neighbor set $\mathcal{N}(i)$, PyTex uses

```{math}
\mathrm{KAM}(i) =
\frac{1}{|\mathcal{N}(i)|}
\sum_{j \in \mathcal{N}(i)} \omega\!\left(g_i, g_j\right),
```

where $\omega$ is the misorientation angle or the symmetry-reduced disorientation angle,
depending on `symmetry_aware`.

## What PyTex Exposes

- regular-grid adjacency for 4- and 8-connectivity
- graph-backed adjacency for irregular coordinate sets
- explicit neighbor order
- optional thresholding
- symmetry-aware or raw-angle evaluation
- mean-style and max-style aggregation
- automatic exclusion of cross-phase pairs on multiphase maps
- optional restriction to an existing grain segmentation

## Example

```python
kam = crystal_map.kernel_average_misorientation_deg(
    order=1,
    threshold_deg=5.0,
    symmetry_aware=True,
    statistic="mean",
)
```

If `grid_shape` is not available, PyTex falls back to graph mode and infers a neighborhood radius
from the scan geometry.

## Phase Semantics

For multiphase maps, PyTex preserves the full coordinate graph but only evaluates KAM on
same-phase pairs. This keeps topology visible without mixing physically incompatible crystal
symmetries into one local-angular metric.

## Related Material

- {doc}`ebsd_import_normalization`
- {doc}`ebsd_grains`
- [../../tex/algorithms/ebsd_local_misorientation.tex](../../tex/algorithms/ebsd_local_misorientation.tex)
- [../../tex/algorithms/ebsd_kam_parameterization.tex](../../tex/algorithms/ebsd_kam_parameterization.tex)
- [../../tex/algorithms/multiphase_ebsd_graph_workflows.tex](../../tex/algorithms/multiphase_ebsd_graph_workflows.tex)
