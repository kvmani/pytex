# EBSD: Kernel-Average Misorientation

Kernel-average misorientation (KAM) is a local orientation-gradient metric on an EBSD map. PyTex
now computes it on a shared neighbor-graph substrate rather than only on a fixed regular-grid
implementation path.

## Definition and Physical Meaning

{ref}`Kernel-average misorientation (KAM) <term-kam>` quantifies local intragranular lattice distortion and orientation gradients from 2D EBSD maps. For a measurement pixel $i$ with valid neighbor set $\mathcal{N}(i)$, PyTex computes:

```{math}
\mathrm{KAM}(i) =
\frac{1}{|\mathcal{N}(i)|}
\sum_{j \in \mathcal{N}(i)} \omega\!\left(g_i, g_j\right),
```

where $\omega(g_i, g_j)$ is the {ref}`disorientation angle <term-disorientation>` (when `symmetry_aware=True`) or the unsymmetrized misorientation angle.

### Relation to Lattice Curvature and Dislocation Density

In deformed crystals, local misorientations across neighboring pixels reflect underlying lattice curvature $\boldsymbol{\kappa} = \nabla\boldsymbol{\omega}$. To leading order, the scalar misorientation angle per neighbor step $\Delta x$ approximates local curvature:

$$
\kappa \approx \frac{\mathrm{KAM}}{\Delta x}
$$

Through the Nye dislocation tensor $\boldsymbol{\alpha} = \boldsymbol{\kappa}^\mathsf{T} - \operatorname{tr}(\boldsymbol{\kappa})\mathbf{I}$, this curvature provides an experimental proxy for {ref}`Geometrically Necessary Dislocation (GND) <term-gnd>` density:

$$
\rho_{\mathrm{GND}} \approx \frac{\alpha_{\mathrm{geom}} \, \mathrm{KAM}}{b \, \Delta x}
$$

where $b$ is the dislocation Burgers vector magnitude (in $\text{m}$), $\Delta x$ is the grid spacing (in $\text{m}$), and $\alpha_{\mathrm{geom}}$ is a dimensionless geometric factor ($\sim 2-3$ depending on boundary type and slip geometry).

High-angle grain boundaries represent discrete structural interfaces rather than continuous lattice curvature. Setting an upper threshold angle $\theta_{\mathrm{threshold}}$ (typically $2^\circ$ to $5^\circ$) excludes intergranular boundaries so that the computed metric selectively reflects intragranular deformation.

## What PyTex Exposes

- regular-grid adjacency for 4- and 8-connectivity
- staggered hexagonal-grid adjacency with the natural 6-connectivity
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

For `grid_kind="hexagonal"`, the default follows the exact six-neighbour logical-row graph and the
result remains a flat value per measured point. If neither rectangular nor hexagonal topology is
available, PyTex falls back to graph mode and infers a neighborhood radius from the coordinates.

## Phase Semantics

For multiphase maps, PyTex preserves the full coordinate graph but only evaluates KAM on
same-phase pairs. This keeps topology visible without mixing physically incompatible crystal
symmetries into one local-angular metric.

## Related Material

- {doc}`ebsd_import_normalization`
- {doc}`ebsd_grains`
- {doc}`/theory/ebsd_local_misorientation`
- {doc}`/theory/ebsd_kam_parameterization`
- {doc}`/theory/multiphase_ebsd_graph_workflows`
