# Inverse Pole Figure (IPF) Color Mapping

**Surface:** `pytex.plotting.ipf.IPFColorKey`, `ipf_color`, `ipf_colors`,
`plot_ipf_key`, `SymmetrySpec.fundamental_sector`,
`SymmetrySpec.reduce_vectors_to_fundamental_sector`, and
`pytex.core.sphere.FundamentalSector`, with workbench operations
`ebsd.map` and `texture.inverse_pole_figure`.

Inverse pole figure (IPF) color mapping represents crystallographic orientation data
by assigning a unique RGB color to the crystal direction $\mathbf{h}$ aligned parallel
to a chosen specimen reference axis $\mathbf{y}$ (such as normal direction $\mathrm{ND}$,
rolling direction $\mathrm{RD}$, or transverse direction $\mathrm{TD}$). Widely employed
in electron backscatter diffraction (EBSD) microstructural analysis, IPF maps visualize
spatial distributions of directional texture.

Because an IPF color encodes only the projection $\mathbf{h} = g^{-1}\mathbf{y}$, the
degree of freedom corresponding to crystal rotation about $\mathbf{y}$ is integrated out.
Consequently, grains displaying identical IPF colors may possess substantial relative
misorientations (up to $45^\circ$ in cubic crystals). Accurate interpretation requires
rigorous specification of the reference direction, crystal point-group symmetry, and
fundamental sector reduction conventions.

## 1. Mathematical formulation

### 1.1 Crystal direction projection

Given an orientation matrix $g \in \mathrm{SO}(3)$ transforming crystal coordinates to
the specimen reference frame, the crystal direction $\mathbf{h}$ parallel to a unit
specimen vector $\mathbf{y}$ is obtained via the inverse transformation:

$$
\mathbf{h} = g^{-1}\mathbf{y}.
$$

The vector $\mathbf{h}$ is normalized to unit length on the unit sphere $\mathbb{S}^2$.
IPF coloring parameterizes $\mathbf{h}$ across the fundamental sector of the crystal
point group.

### 1.2 Symmetry reduction to the fundamental sector

The rotational point-group symmetry $\mathcal{G}_{\text{xtal}}$ partitions $\mathbb{S}^2$
into equivalent spherical domains. A canonical **fundamental sector** (or standard
stereographic triangle) contains exactly one representative vector $\tilde{\mathbf{h}}$
from each symmetry orbit:

$$
\tilde{\mathbf{h}} = \operatorname{reduce}(\mathbf{h}, \mathcal{G}_{\text{xtal}}, \text{antipodal}).
$$

For centrosymmetric diffraction phenomena or non-polar direction analysis,
`antipodal=True` applies Friedel symmetry, identifying antipodal vectors
($\mathbf{h} \sim -\mathbf{h}$) and halving the required sector area. For cubic $m\bar{3}m$
symmetry, the standard stereographic triangle is bounded by vertices:

$$
\mathbf{v}_1 = [001], \quad \mathbf{v}_2 = [101], \quad \mathbf{v}_3 = [111].
$$

### 1.3 Barycentric color coordinates

Let $\mathbf{B} = \begin{bmatrix} \mathbf{v}_1 & \mathbf{v}_2 & \mathbf{v}_3 \end{bmatrix} \in \mathbb{R}^{3 \times 3}$
be the basis matrix formed by the normalized corner vectors of the fundamental spherical
triangle. For any reduced direction $\tilde{\mathbf{h}}$, Cartesian coordinate
decomposition yields barycentric weights $\mathbf{c} = [c_1, c_2, c_3]^{\mathsf{T}}$:

$$
\mathbf{B}\,\mathbf{c} = \tilde{\mathbf{h}} \implies \mathbf{c} = \mathbf{B}^{-1}\tilde{\mathbf{h}}.
$$

To guarantee numerical admissibility on the standard 2-simplex:
1. Negative components arising from floating-point rounding along sector boundaries are clipped:
   $$
   c_i \leftarrow \max(c_i, 0).
   $$
2. Weights are normalized to unit sum:
   $$
   c_i \leftarrow \frac{c_i}{\sum_{k=1}^3 c_k}.
   $$

The base RGB triplet is evaluated as a convex combination of assigned vertex colors
$\mathbf{C}_{\text{vertex}} \in \mathbb{R}^{3 \times 3}$:

$$
\mathbf{rgb}_0 = \sum_{i=1}^3 c_i\,\mathbf{C}_{\text{vertex}, i}.
$$

By convention, cubic vertices $[001]$, $[101]$, and $[111]$ are assigned pure red
$[1, 0, 0]$, green $[0, 1, 0]$, and blue $[0, 0, 1]$, respectively.

### 1.4 Gamut saturation and perceptual contrast

Linear barycentric interpolation yields low color saturation near the triangle centroid,
resulting in muted grayish tones that obscure subtle grain misorientations. PyTex applies
a power-law gamma transfer function followed by channel-peak normalization:

$$
\mathbf{rgb}_1 = \mathbf{rgb}_0^{1 / \gamma_{\text{sat}}}, \qquad \mathbf{rgb} = \frac{\mathbf{rgb}_1}{\max(\mathbf{rgb}_1)},
$$

where $\gamma_{\text{sat}}$ denotes the saturation exponent (defaulting to 0.5).

> [!WARNING]
> The transformation $\mathbf{h} \mapsto \mathbf{rgb}$ is non-isometric. Distance in
> RGB color space does not correspond linearly to crystallographic misorientation angle.
> Color similarity must not be interpreted as low misorientation.

## 2. Sector geometry across crystal systems

| Crystal System | Laue Class / Point Group | Fundamental Sector Boundaries | Primary Vertex Colors |
| --- | --- | --- | --- |
| Cubic | $m\bar{3}m$ ($O_h$) | $[001] - [101] - [111]$ | Red $[001]$, Green $[101]$, Blue $[111]$ |
| Hexagonal | $6/mmm$ ($D_{6h}$) | $[0001] - [10\bar{1}0] - [2\bar{1}\bar{1}0]$ | Red $[0001]$, Green $[10\bar{1}0]$, Blue $[2\bar{1}\bar{1}0]$ |
| Tetragonal | $4/mmm$ ($D_{4h}$) | $[001] - [100] - [110]$ | Red $[001]$, Green $[100]$, Blue $[110]$ |
| Orthorhombic | $mmm$ ($D_{2h}$) | $[001] - [100] - [010]$ | Red $[001]$, Green $[100]$, Blue $[010]$ |
| Trigonal / Monoclinic | Various | Extended spherical polygons | Multi-triangle barycentric subdivision or reference octant projection |

For low-symmetry crystal systems where the fundamental domain is bounded by more than
three vertices, the domain is partitioned into contiguous spherical triangles, or
anchored onto canonical Cartesian octants.

## 3. Configuration parameters and diagnostic constraints

| Parameter | Default | Description and Validity Range |
| --- | --- | --- |
| `crystal_symmetry` | (Mandatory) | Crystal point-group symmetry specifying sector topology. |
| `specimen_direction` | `[0, 0, 1]` ($\mathrm{ND}$) | Specimen reference vector $\mathbf{y} \in \mathbb{S}^2$ projected into the crystal frame. |
| `antipodal` | `True` | Enforces Friedel symmetry ($\mathbf{h} \sim -\mathbf{h}$), halving the fundamental sector. |
| `saturation_gamma` | `0.5` | Non-linear gamma parameter ($\gamma > 0$) enhancing color saturation. |
| `resolution_deg` | `1.0` | Angular grid spacing in degrees for rendering the legend mesh ($0 < \Delta \theta \le 15^\circ$). |

### Invariant validations and error conditions

- **Non-crystal coordinate frames:** `IPFColorKey` requires an explicit crystal reference
  frame. Attempting to initialize with a specimen or laboratory frame raises a construction-time error.
- **Out-of-bounds directions:** If a numerical vector fails fundamental sector bounds
  after reduction (barycentric sum $\le 0$), PyTex raises an evaluation exception rather
  than silently substituting uncalibrated colors.

## 4. Downstream applications and visualization

The `IPFColorKey` instance serves as a synchronized color engine across PyTex workflows:

- **EBSD spatial mapping:** `ebsd.map` colors orientation pixel grids, supporting
  band-contrast shading and grain-boundary vector overlays.
- **Inverse pole figure distributions:** `texture.inverse_pole_figure` projects
  discrete orientation scatter points or continuous density functions onto the standard
  stereographic triangle.
- **Standalone color legends:** `plot_ipf_key` generates publication-ready vector legends
  with stereographic boundary arcs and Miller indices labels.

## Verification

- `tests/unit/test_ipf_key.py`: Validates vertex RGB values, rotational symmetry invariance,
  vectorized array mapping, and stereographic projection bounds.
- Executable worked examples:
  - {doc}`../examples/generated/ipf-coloring`
  - {doc}`../examples/generated/crystal_geometry`

## See also

- {doc}`../theory/ipf_color_keys` — Analytical derivations of sector geometries and barycentric coordinates.
- {doc}`../theory/fundamental_region_reduction` — Algorithms for canonical symmetry reduction on $\mathbb{S}^2$.
- {doc}`../concepts/symmetry_and_fundamental_regions` — Point groups, Laue classes, and orientation space topology.

## References

### Normative

- Nolze, G. & Hielscher, R. (2016). Orientations - perfectly colored. *Journal of Applied
  Crystallography* **49**, 1786–1802. <https://doi.org/10.1107/S1600576716012942>

### Informative

- Engler, O. & Randle, V. (2010). *Introduction to Texture Analysis: Macrotexture,
  Microtexture, and Orientation Mapping*, 2nd ed. CRC Press. <https://doi.org/10.1201/9781420063660>
- Schwartz, A. J., Kumar, M., Adams, B. L. & Field, D. P., eds. (2009). *Electron
  Backscatter Diffraction in Materials Science*, 2nd ed. Springer.
  <https://doi.org/10.1007/978-0-387-88136-2>
