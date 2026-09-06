# EBSD Kernel Average Misorientation (KAM) And Dislocation Density

This note fixes the mathematical definitions, topological neighborhood policies, boundary thresholding
constraints, and dislocation density models implemented across `pytex.ebsd`:
`CrystalMap.kernel_average_misorientation_deg`, `CrystalMap.neighbor_graph`, and `pytex.ebsd.gnd`.

Kernel Average Misorientation (KAM) is the primary local misorientation metric used to quantify intragranular
plastic strain, subgrain cell formation, and stored deformation energy from Electron Backscatter
Diffraction (EBSD) orientation maps.

## 1. Mathematical Formulation

### Pairwise Disorientation

Let $x_i$ be an EBSD measurement point with crystal orientation represented by proper rotation matrix
$\mathbf{g}_i \in \mathrm{SO}(3)$ carrying specimen Cartesian coordinates to crystal Cartesian coordinates.
For any neighboring site $x_j$, the orientation difference in the specimen frame is
$\Delta \mathbf{g}_{ij} = \mathbf{g}_j \mathbf{g}_i^\mathsf{T}$.

Under the crystal proper point group symmetry $G = \{ \mathbf{S}_k \}_{k=1}^{N_{\text{sym}}}$, there are
$N_{\text{sym}}$ symmetry-equivalent representations of the orientation. The **disorientation angle**
$\theta(\mathbf{g}_i, \mathbf{g}_j)$ is the minimum rotation angle among all crystallographically
equivalent rotation operations:

$$
\theta(\mathbf{g}_i, \mathbf{g}_j) = \min_{\mathbf{S} \in G} \arccos\left(
    \frac{\operatorname{Tr}(\mathbf{S}\,\mathbf{g}_j\,\mathbf{g}_i^\mathsf{T}) - 1}{2}
\right)
$$ (eq-kam-disorientation)

Equation {eq}`eq-kam-disorientation` satisfies $0 \le \theta \le \theta_{\max}$, where $\theta_{\max}$ is the
maximum disorientation angle of the crystal symmetry's fundamental zone (for example, $62.8^\circ$ for cubic
$m\bar{3}m$ symmetry).

### Kernel Average Misorientation (KAM)

For each point $x_i$, let $\mathcal{N}(x_i)$ denote the set of admissible neighboring measurement sites.
The mean Kernel Average Misorientation is defined as:

$$
\mathrm{KAM}(x_i) = \frac{1}{|\mathcal{N}(x_i)|} \sum_{x_j \in \mathcal{N}(x_i)} \theta(\mathbf{g}_i, \mathbf{g}_j)
$$ (eq-kam-mean)

When the maximum local distortion is sought (e.g., to detect localized micro-shear bands or sub-boundary
nucleation), PyTex also provides the maximum operator:

$$
\mathrm{KAM}_{\max}(x_i) = \max_{x_j \in \mathcal{N}(x_i)} \theta(\mathbf{g}_i, \mathbf{g}_j)
$$ (eq-kam-max)

If $|\mathcal{N}(x_i)| = 0$ (such as an isolated indexed pixel or a pixel where all neighbors exceed the
filtering threshold), $\mathrm{KAM}(x_i) = 0^\circ$ by definition.

## 2. Neighborhood Topologies And Graph Distance

The definition of $\mathcal{N}(x_i)$ depends strictly on the physical grid geometry of the EBSD acquisition:

### Rectangular Grids (Square / Orthogonal)

For scans acquired on a regular Cartesian grid with step sizes $\Delta x, \Delta y$:

1. **4-connected (Manhattan / von Neumann)**:
   Connects the 4 nearest orthogonal neighbors at grid offsets $(0, \pm 1)$ and $(\pm 1, 0)$.
   The topological distance corresponds to the $\ell_1$ norm:
   $$\lVert \Delta \mathbf{r} \rVert_1 = |\Delta u| + |\Delta v| = 1$$

2. **8-connected (Chebyshev / Moore)**:
   Connects both orthogonal and diagonal nearest neighbors (8 points).
   The topological distance corresponds to the $\ell_\infty$ norm:
   $$\lVert \Delta \mathbf{r} \rVert_\infty = \max(|\Delta u|, |\Delta v|) = 1$$

### Hexagonal Grids (Staggered Rows)

Modern EBSD systems frequently scan in hexagonal patterns because each interior site is equidistant from
exactly **six nearest neighbors**, giving isotropic sampling without diagonal metric distortion.

PyTex represents hexagonal grids natively as a staggered-row adjacency graph:
- Even and odd rows are offset by $\frac{1}{2} \Delta x$.
- Each interior point connects to 6 equidistant adjacent neighbors ($k=1$).
- Neighbor shells of order $k \ge 1$ are constructed using cumulative shortest-path graph distance rather
  than Euclidean coordinate approximations, preserving the true topology across ragged row boundaries.

```text
      (u-1, v+1)     (u, v+1)
            \       /
   (u-1, v) — (u, v) — (u+1, v)    [Hexagonal 6-neighbor topology]
            /       \
      (u-1, v-1)     (u, v-1)
```

### Higher-Order Shells ($k > 1$)

For shell order $k > 1$:
- In **cumulative mode** (default), $\mathcal{N}_k(x_i)$ encompasses all nodes with graph distance $1 \le d(x_i, x_j) \le k$.
- Higher orders average over larger physical interaction distances $d = k \cdot \Delta x$, smoothing noise
  at the cost of spatial resolution.

## 3. Boundary Thresholding And Segmentation Constraints

Without filtering, a kernel straddling a high-angle grain boundary (HAGB) includes disorientation angles of
$15^\circ$ to $60^\circ$. This produces an artificial halo of high KAM along all grain boundaries that reflects
grain topology rather than intragranular dislocation structures.

PyTex provides two strict mechanisms to prevent grain boundary contamination:

### Misorientation Threshold ($\theta_{\text{thresh}}$)

An admissible neighbor pair must satisfy:

$$
\theta(\mathbf{g}_i, \mathbf{g}_j) \le \theta_{\text{thresh}}
$$ (eq-kam-threshold)

Pairs with $\theta > \theta_{\text{thresh}}$ are excluded from both the summation and the denominator
$|\mathcal{N}(x_i)|$. In metallographic practice, $\theta_{\text{thresh}}$ is conventionally chosen
between $2^\circ$ and $5^\circ$ to match the low-angle boundary cutoff.

### Grain Segmentation Masking

When a `GrainSegmentation` is provided, neighbors are restricted to points sharing the exact same grain
identifier:

$$
\mathcal{N}_{\text{grain}}(x_i) = \left\{ x_j \in \mathcal{N}(x_i) \;\middle|\; \operatorname{grain}(x_j) = \operatorname{grain}(x_i) \right\}
$$ (eq-kam-segmentation)

This is more physically defensible than an angular threshold alone: it prevents subgrains that happen to differ
by less than $\theta_{\text{thresh}}$ across an actual grain boundary from falsely averaging across the boundary.

Furthermore, pairs across phase boundaries are always excluded: crystallographic misorientation is undefined
between different space groups or lattice metrics.

## 4. Evaluation Of Dislocation Density From KAM

Intragranular misorientation arises from excess dislocations of one sign accommodating lattice curvature:
the **Geometrically Necessary Dislocations (GNDs)**.

### The Read-Shockley Sub-Boundary Relation

A low-angle tilt boundary with misorientation angle $\theta$ and dislocation spacing $h$ has Burgers vector
magnitude $b = \lVert \mathbf{b} \rVert$ related by Frank's formula:

$$
\theta = 2 \arcsin\left(\frac{b}{2h}\right) \approx \frac{b}{h}
$$

For a subgrain cell of diameter / kernel dimension $d$ with dislocation spacing $h$, the line length per
unit volume (dislocation density $\rho$) is:

$$
\rho = \frac{\text{total line length}}{\text{volume}} \approx \frac{\text{area of sub-boundary} \times (1/h)}{\text{volume}}
\approx \frac{1}{h \cdot d} = \frac{\theta}{b \cdot d}
$$

### The Scalar GND Density Model

Accounting for 3D boundary geometry and mixed tilt/twist character, the scalar GND density is estimated
from the KAM via:

$$
\rho_{\mathrm{GND}} \approx \frac{\alpha\,\mathrm{KAM}}{b\,d}
$$ (eq-kam-gnd-density)

where:
- $\mathrm{KAM}$ is the kernel average misorientation expressed in **radians**,
- $b = \lVert \mathbf{b} \rVert$ is the magnitude of the Burgers vector of the primary slip system in meters
  (e.g., $b = \frac{a}{\sqrt{2}} \approx 0.25\,\mathrm{nm}$ for fcc nickel/copper; $b = \frac{a\sqrt{3}}{2} \approx 0.248\,\mathrm{nm}$ for bcc iron),
- $d$ is the effective kernel interaction length in meters ($d = k \cdot \Delta x$),
- $\alpha$ is a dimensionless geometric factor:
  - $\alpha \approx 2$ for pure tilt sub-boundaries,
  - $\alpha \approx 3$ for mixed boundaries,
  - $\alpha \approx 4$ for pure twist networks.

### Resolution Dependence And Lower-Bound Nature

GND density calculated from 2D orientation maps is subject to two fundamental physical limits:

1. **Step Size Invariance Limit**: As step size $\Delta x$ decreases, higher curvature gradients are resolved,
   meaning measured $\rho_{\mathrm{GND}}$ increases. Therefore, GND values are valid for comparison **only at
   identical step size and kernel configuration**.
2. **Statistically Stored Dislocations (SSDs)**: Dislocations that form mutually cancelling multipoles or dipoles
   produce no net lattice curvature. Consequently, $\rho_{\mathrm{GND}}$ is always a **lower bound** to the
   total dislocation density ($\rho_{\text{total}} = \rho_{\mathrm{GND}} + \rho_{\text{SSD}}$).

## 5. Comparison: KAM-Based vs Nye Tensor Curvature Approaches

PyTex provides two independent paths to dislocation content:

| Feature | Kernel Average Misorientation (KAM) | Nye Tensor Lattice Curvature ($\alpha_{ij}$) |
|---|---|---|
| **Module** | `pytex.ebsd.models.CrystalMap.kernel_average_misorientation_deg` | `pytex.ebsd.gnd.geometrically_necessary_dislocation_density` |
| **Input** | Local scalar angle distribution $\theta(g_i, g_j)$ | Spatial gradients of rotation vectors $\partial \omega_i / \partial x_j$ |
| **Symmetry** | **Reduced** to crystallographic fundamental zone | **Unreduced**; preserves continuous vector gradient direction |
| **Output** | Degrees ($^\circ$) or scalar density ($\mathrm{m}^{-2}$) | Full tensor components $\alpha_{ij}$ and scalar sum ($\mathrm{m}^{-2}$) |
| **Spatial Detail** | Directionless isotropic scalar | Resolves in-plane gradient components ($\alpha_{01}, \alpha_{02}, \alpha_{10}, \alpha_{12}$) |
| **Standard Use** | Routine deformation mapping & recovery monitoring | Rigorous plastic strain gradient & stress state modeling |

See [Lattice Curvature and GND Density](lattice_curvature_and_gnd_density.md) for the tensor derivation.

## Normative And Informative References

1. Humphreys, F. J., *Characterisation of fine-scale microstructures by electron backscatter diffraction (EBSD)*,
   Journal of Materials Science 36 (2001) 3833–3854.
2. Brewer, L. N., Field, D. P., Merriman, C. C., *Mapping and analyzing local orientation gradients in EBSD data*,
   in Electron Backscatter Diffraction in Materials Science, Springer (2009) 251–262.
3. Wilkinson, A. J., Randman, D., *Determination of elastic strain fields and geometrically necessary dislocation
   distributions using electron backscatter diffraction*, Philosophical Magazine 90 (2010) 1159–1167.
4. Nye, J. F., *Some geometrical relations in dislocated crystals*, Acta Metallurgica 1 (1953) 153–162.
5. Read, W. T., Shockley, W., *Dislocation models of crystal grain boundaries*, Physical Review 78 (1950) 275–289.
