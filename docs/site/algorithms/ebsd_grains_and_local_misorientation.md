# Grains, Local Misorientation, And Dislocation Density From An EBSD Map

**Surface:** `CrystalMap.segment_grains`, `GrainSegmentation`, `Grain`,
`GrainBoundaryNetwork`, `CrystalMap.kernel_average_misorientation_deg`,
the GROD/GOS/GAM family on `GrainSegmentation`, and
`pytex.ebsd.gnd.lattice_curvature_tensor`,
`nye_dislocation_density_tensor`,
`geometrically_necessary_dislocation_density`; exposed as the workbench
operations `ebsd.map`, `ebsd.grod`, `ebsd.kam` and `ebsd.scan_summary`.

An EBSD scan is a grid of orientations. Everything a materials scientist reads
off it — grain size, subgrain structure, stored deformation, dislocation
content — is *derived*, and every derived number carries the parameters that
produced it. This page states each derivation and the parameter that decides it,
because a grain size quoted without its threshold, or a KAM map without its
kernel, is not reproducible.

```{figure} ../../figures/ebsd_grain_metrics_algorithm.svg
:alt: Four-lane flow sheet. Lane 1 turns the orientation grid into a neighbour
  graph and computes symmetry-reduced pair disorientations. Lane 2 thresholds
  the edges and takes connected components. Lane 3 derives KAM, GROD, GOS and
  GAM, distinguished by what each compares a point with. Lane 4 forms the
  curvature and Nye tensors and reports a lower-bound GND density.
:width: 100%

Segmentation, the local metrics, and the route to dislocation density.
```

## 1. Segmentation: from points to grains

### 1.1 The algorithm

```text
input : orientations on a grid, threshold theta_c, connectivity, symmetry_aware

1  build the first-shell neighbour graph (4/8 on a square grid, 6 on hexagonal)
2  drop every pair that crosses a phase boundary
3  for each remaining pair, compute the misorientation angle
       symmetry_aware=True  -> disorientation (minimum over the symmetry orbit)
       symmetry_aware=False -> raw rotation angle
4  keep the pairs with angle <= theta_c as graph edges
5  grains are the connected components of that edge set
6  number the grains by each component's lowest member index
```

Step 5 is a flood fill, computed as connected components of a sparse adjacency
matrix rather than a Python union-find — the traversal is the same, the compiled
version is what makes a full-size map tractable. Step 6 is done explicitly
rather than inherited from the traversal order, so grain ids are a stable
function of the data and not of a library's internals.

### 1.2 The threshold decides the answer

$\theta_c$ is conventionally $5\text{–}15^\circ$, and **it is not a detail**:
it decides whether subgrains are resolved as separate grains, so it moves grain
size, grain count, and every distribution derived from them. `GrainSegmentation`
stores it rather than merely applying it, so a downstream metric cannot be
reported without the criterion that produced it.

### 1.3 The failure mode that is inherent, not a bug

Flood fill merges points connected by a *chain* of small steps. **A grain with a
continuous orientation gradient can therefore exceed the threshold end to end**
while never exceeding it between neighbours — a heavily deformed grain is one
grain by this criterion however far it has rotated across its width. This is a
property of the definition, shared by every flood-fill segmentation in the
field, and it is the reason the local-misorientation family of section 2 exists:
the gradient the segmentation absorbs is exactly what GROD then measures.

### 1.4 Symmetry awareness

With `symmetry_aware=True` (the default) the pair angle is the **disorientation**
— the minimum rotation angle over the symmetry orbit of the misorientation. For
a cubic phase the raw angle can be up to $180^\circ$ where the disorientation is
a few degrees, so segmenting on raw angles fragments grains at random. Turn it
off only to reproduce a tool that does the same.

## 2. Local misorientation: GROD, KAM, GOS, GAM

Four numbers, often confused, differing in *what each point is compared with*.

| Metric | Compares a point with | Answers |
| --- | --- | --- |
| **KAM** | its immediate neighbours | how sharp is the local gradient |
| **GROD** | its own grain's reference orientation | how far has this point rotated within its grain |
| **GOS** | grain mean, averaged over the grain | how deformed is this grain overall |
| **GAM** | neighbours, averaged over the grain | how much local gradient does this grain hold |

KAM is a *short-wavelength* measure and GROD a *long-wavelength* one. A grain
with a smooth, large total rotation has high GROD and low KAM; a grain with a
sharp subgrain wall has the reverse. Reporting one and calling it "deformation"
loses the distinction that matters.

### 2.1 KAM

```text
input : neighbour shell `order`, connectivity, threshold_deg, statistic, segmentation

1  build the neighbour graph at the requested shell
2  drop pairs crossing a phase boundary            -- always
3  drop pairs above threshold_deg                  -- if given
   or drop pairs whose ends are in different grains -- if a segmentation is given
4  per point: mean (or max) of the surviving pair misorientations
5  points with no admissible neighbour report zero
```

**Excluding the boundary is not optional in practice.** Without step 3 a pixel
on a grain boundary reports the *boundary misorientation* — tens of degrees —
rather than the local gradient, and the KAM map becomes a boundary map with a
deformation map faintly visible underneath. Two ways to exclude are offered, and
they are not equivalent: `threshold_deg` is the conventional one; passing a
`segmentation` is stricter and more physical, because it excludes by grain
membership rather than by angle, and it is *required* for the GAM definition.

`order` sets the kernel radius. A larger kernel smooths and lowers KAM, so **KAM
values are comparable only at equal step size and equal order** — the quantity
depends on the measurement grid, not only on the material.

## 3. GND density: from curvature to dislocations

### 3.1 The Nye route (`method="curvature"`, the default)

```text
1  form the lattice curvature tensor from the orientation gradient
       kappa_ij = d(omega_i) / d(x_j),  omega the rotation vector field
2  convert to the Nye dislocation density tensor
       alpha = kappa^T - trace(kappa) I
3  sum the absolute values of the measurable components
4  divide by the Burgers vector magnitude
```

A two-dimensional surface scan measures gradients in two directions only, so
**not every component of $\boldsymbol{\alpha}$ is measurable**: the derivative
normal to the surface is unavailable. The result is therefore a **lower bound**
on the GND density, and it is labelled as one. It is also, by construction, only
the *geometrically necessary* content — statistically stored dislocations
produce no net lattice curvature and are invisible to the method entirely.

### 3.2 The KAM route (`method="kam"`)

$$
\rho \;\approx\; \frac{2\theta}{b\,u}
$$

with $\theta$ the KAM in radians, $b$ the Burgers vector and $u$ the step size.
Cruder — it discards the *direction* of the gradient, which is the whole content
of the Nye tensor — but it is what a large part of the EBSD literature reports,
so it is provided for comparability rather than because it is better. When
comparing with a published number, check which route that number used.

### 3.3 Units and the things that silently corrupt them

| Parameter | Meaning | Failure if wrong |
| --- | --- | --- |
| `burgers_vector_nm` | Cu 0.2556, α-Fe 0.2483, Al 0.2863 | density scales as $1/b$ |
| `step_scale_m` | metres per map coordinate unit; default treats them as µm | density scales as $1/u$; a factor of $10^{6}$ if the map is in metres |
| `kam_threshold_deg` | boundary exclusion for the KAM route | boundary pixels return meaningless densities |

`NaN` is returned where the gradient could not be measured — across a phase
boundary, for instance — rather than zero, because zero is a physical claim and
"not measurable here" is not.

## 4. Boundaries

`GrainBoundaryNetwork` carries the segments between grains with their
misorientation, so boundary character can be classified: low-angle against
high-angle at the segmentation threshold, and CSL relationships through
`pytex.ebsd.csl`. Boundary *length* inherits the grid, so a boundary traced on a
square grid is quantised to the step and its length is a step-dependent number.

## 5. Reading the parameters back

Every number on this page depends on at least one choice. For a result to be
reproducible, all of these must travel with it:

| Metric | Must be reported with |
| --- | --- |
| grain size, count | $\theta_c$, connectivity, symmetry_aware |
| KAM | step size, `order`, `threshold_deg` or the segmentation, `statistic` |
| GROD, GOS, GAM | the segmentation, hence $\theta_c$ |
| GND | method, $b$, step size, and (for the KAM route) the threshold |

`GrainSegmentation` stores its settings for exactly this reason, and the
workbench result objects carry them into `describe()` so the prose beside a
figure states the conventions the figure was computed under.

## 6. Cost

| Stage | Cost |
| --- | --- |
| Neighbour graph | $O(n)$ on a regular grid, built vectorised rather than per point |
| Pair misorientations | $O(n_{\text{pairs}} \cdot |G|)$ for the symmetry orbit, vectorised |
| Connected components | near-linear, in compiled code |
| KAM | one pass over the pair list per shell |
| Curvature/GND | finite differences over the grid, vectorised |

## Verification

- Grain segmentation, KAM and the GROD family against pinned expectations, in
  {doc}`../examples/generated/ebsd`.
- The KAM parameterisation and its dependence on kernel and step, in
  {doc}`../theory/ebsd_kam_parameterization`.

## See also

- {doc}`../theory/ebsd_grain_segmentation_and_grod` — the definitions in full.
- {doc}`../theory/ebsd_local_misorientation` — how the four metrics differ.
- {doc}`../theory/lattice_curvature_and_gnd_density` — the Nye tensor derivation.
- {doc}`../theory/ebsd_boundaries_and_cleanup` — boundary extraction and the
  cleanup operations that precede segmentation.
- {doc}`ipf_coloring` — how the map is coloured once it is segmented.

## References

### Normative

- Nye, J. F. (1953). Some geometrical relations in dislocated crystals. *Acta
  Metallurgica* **1**, 153-162.
  <https://doi.org/10.1016/0001-6160(53)90054-6>
- Pantleon, W. (2008). Resolving the geometrically necessary dislocation content
  by conventional electron backscattering diffraction. *Scripta Materialia*
  **58**, 994-997. <https://doi.org/10.1016/j.scriptamat.2008.01.050>

### Informative

- Wright, S. I., Nowell, M. M. & Field, D. P. (2011). A review of strain
  analysis using electron backscatter diffraction. *Microscopy and
  Microanalysis* **17**, 316-329.
  <https://doi.org/10.1017/S1431927611000055>
- Kamaya, M. (2011). Assessment of local deformation using EBSD. *Ultramicroscopy*
  **111**, 1189-1199. <https://doi.org/10.1016/j.ultramic.2011.02.004>
