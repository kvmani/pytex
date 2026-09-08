# Kearns Texture Parameters and the Basal Orientation Tensor

**Surface:** `pytex.texture.kearns.pole_orientation_tensor`,
`kearns_from_orientations`, `kearns_from_pole_figure`, `kearns_from_odf`,
`kernel_axis_shrinkage`, `KearnsReport`, with workbench operations
`kearns.from_odf` and the Kearns analysis interface.

In hexagonal close-packed (HCP) metals such as zirconium and titanium alloys,
the orientation distribution of basal poles $[0001]$ governs fundamental
anisotropic material behavior, including thermal expansion, hydride habit-plane
precipitation, anisotropic plastic flow, and radiation-induced growth (Kearns, 1965;
Holt & Aldridge, 1985). The **Kearns orientation parameter** $f_d$ represents the
effective volume fraction of grains with basal poles resolved along a specified
specimen direction $\mathbf{d}$ (conventionally rolling direction $\mathrm{RD}$,
transverse direction $\mathrm{TD}$, and normal direction $\mathrm{ND}$).

PyTex models Kearns parameters by evaluating the underlying second-moment
**basal orientation tensor** $\mathbf{A} = \langle \mathbf{v}\mathbf{v}^{\mathsf{T}} \rangle$.
This page presents the unified tensor formulation, details the three experimental
evaluation routes (discrete orientations, diffraction pole figures, and reconstructed ODFs),
derives analytical kernel deconvolution, and delineates the numerical validity limits
of each modality.

## 1. The basal orientation tensor

### 1.1 Definition and tensor quadratic form

For a crystal orientation characterized by a unit basal pole $\mathbf{v} \parallel [0001]$
in the specimen reference frame, the angle $\alpha_d$ relative to a specimen direction
$\mathbf{d}$ satisfies $\cos\alpha_d = \mathbf{v} \cdot \mathbf{d}$. The Kearns parameter
$f_d$ is defined as the ensemble-averaged squared direction cosine:

$$
f_d = \bigl\langle \cos^{2}\alpha_{d} \bigr\rangle
    = \bigl\langle (\mathbf{v}\cdot\mathbf{d})^2 \bigr\rangle
    = \mathbf{d}^{\mathsf{T}} \bigl\langle \mathbf{v}\,\mathbf{v}^{\mathsf{T}} \bigr\rangle \mathbf{d}.
$$

Defining the symmetric, second-rank **basal orientation tensor** $\mathbf{A}$:

$$
\mathbf{A} = \bigl\langle \mathbf{v}\,\mathbf{v}^{\mathsf{T}} \bigr\rangle
           = \int_{\mathrm{SO}(3)} \left(g\,\hat{\mathbf{c}}\right)\left(g\,\hat{\mathbf{c}}\right)^{\mathsf{T}} f(g)\,\mathrm{d}g,
$$

where $\hat{\mathbf{c}}$ is the unit $[0001]$ direction in the crystal frame. The Kearns
parameter along any arbitrary specimen unit vector $\mathbf{d}$ is evaluated as the quadratic form:

$$
f_d = \mathbf{d}^{\mathsf{T}}\mathbf{A}\mathbf{d}.
$$

### 1.2 Mathematical properties and physical invariants

1. **Trace identity and unit partition.** Because $\mathbf{v}$ is a unit vector
   ($\lVert\mathbf{v}\rVert^2 = v_x^2 + v_y^2 + v_z^2 = 1$):
   $$
   \operatorname{tr}(\mathbf{A}) = \langle v_x^2 + v_y^2 + v_z^2 \rangle = 1.
   $$
   Consequently, for any mutually orthogonal triad of reference axes $(\mathbf{d}_1, \mathbf{d}_2, \mathbf{d}_3)$:
   $$
   f_{\mathbf{d}_1} + f_{\mathbf{d}_2} + f_{\mathbf{d}_3} = 1.
   $$
   A measured Kearns triple failing this identity indicates experimental sampling or
   quadrature errors.
2. **Spectral decomposition and principal texture axes.** The eigendecomposition
   $\mathbf{A} = \mathbf{V}\mathbf{\Lambda}\mathbf{V}^{\mathsf{T}}$ identifies the
   principal Kearns factors (eigenvalues $\lambda_1, \lambda_2, \lambda_3$) and the
   corresponding principal texture directions (eigenvectors). In rolled or pilgered
   tubing, principal basal axes frequently tilt away from orthogonal processing axes
   by $20^\circ$ to $40^\circ$; tensor diagonalization locates the true peak orientation
   independent of the initial coordinate frame.
3. **Inversion symmetry.** The dyadic product $\mathbf{v}\mathbf{v}^{\mathsf{T}}$ is
   strictly invariant under sign inversion ($\mathbf{v} \to -\mathbf{v}$), naturally
   respecting Friedel symmetry in diffraction measurements.
4. **Isotropic baseline.** For a completely random orientation distribution,
   $\mathbf{A} = \frac{1}{3}\mathbf{I}$, yielding $f_d = \frac{1}{3}$ for all directions.

## 2. Experimental evaluation routes

```
                         Orientations / EBSD
                                 │
                                 ▼
                     Route A: Discrete Average
                                 │
                                 ▼
Diffraction Pole Figure ──► Route B: Spherical Quadrature ──► Basal Orientation Tensor A
                                 ▲                                 │
                                 │                                 ▼
Reconstructed ODF       ──► Route C: Continuous Integral     KearnsReport (f_d, axes)
```

### 2.1 Route A: Discrete orientations (EBSD)

Implemented in `kearns_from_orientations`. When individual orientations $g_i$ and weights
$w_i$ are available from EBSD or diffraction spot indexing, $\mathbf{A}$ is calculated
directly:

$$
\mathbf{A} = \frac{\sum_{i=1}^N w_i\,\mathbf{v}_i\mathbf{v}_i^{\mathsf{T}}}{\sum_{i=1}^N w_i},
$$

where $\mathbf{v}_i = g_i [0001]$. This direct summation requires no spherical quadrature
or inversion steps. However, statistical uncertainty depends directly on grain sampling
density; small EBSD map areas may under-sample macrotexture variations.

### 2.2 Route B: Direct pole figure quadrature

Implemented in `kearns_from_pole_figure`. Following Baron et al. (1990) and the Kern–Bergmann
pseudo-normalization formulation for experimental $(0002)$ pole figures:

$$
f_d = \frac{\int_0^{2\pi}\int_0^{\alpha_{\text{max}}} P(\alpha, \beta) (\mathbf{v}(\alpha,\beta)\cdot\mathbf{d})^2 \sin\alpha\,\mathrm{d}\alpha\,\mathrm{d}\beta}
           {\int_0^{2\pi}\int_0^{\alpha_{\text{max}}} P(\alpha, \beta) \sin\alpha\,\mathrm{d}\alpha\,\mathrm{d}\beta}.
$$

#### Solid-angle weighting versus sampling mode

Experimental pole figures must be integrated using solid-angle weights consistent with
their data structure:
- **`sampled_density` (equiangular grids):** Each grid point requires solid-angle weighting
  $w_{ij} \propto \sin\alpha_i\,\Delta\alpha\,\Delta\beta$. Omitting the $\sin\alpha$ metric
  factor artificially overweights directions near the tilt pole by up to $50\%$.
- **`scattered_poles` (discrete reflections):** Intensities represent integrated counts,
  and are summed directly without spherical coordinate weighting.

#### Limitations of Route B in low-intensity sections

Incomplete pole figures truncated at $\alpha_{\text{max}} < 90^\circ$ rely on the assumption
that unmeasured peripheral densities match the measured mean. Furthermore, in specimen
sections where basal plane normals are nearly absent (such as the $\mathrm{ND}-\mathrm{TD}$
section of strongly textured zirconium tubing), measured $(0002)$ diffraction counts are
near the detector background noise. Normalizing by a vanishing integral amplifies counting
noise and produces severe numerical instability (Mani Krishna et al., 2011). In such
cases, Route C is mathematically preferred.

### 2.3 Route C: Reconstructed ODF integration

Implemented in `kearns_from_odf`. When multiple diffraction pole figures
(e.g., $\{10\bar{1}0\}$, $\{0002\}$, $\{10\bar{1}1\}$) are co-inverted into an ODF,
the complete orientation distribution provides basal pole densities even when direct
$(0002)$ measurements are incomplete or noisy.

#### Kernel deconvolution for discrete ODF models

A discrete ODF models continuous density as a superposition of spherical kernels:
$f(g) = \sum_j w_j K(g; g_j)$. This introduces two distinct tensor quantities:
1. **Support tensor $\mathbf{A}_{\text{support}}$:** The second moment evaluated purely
   at the discrete dictionary nodes $g_j$.
2. **Density tensor $\mathbf{A}_{\text{density}}$:** The second moment evaluated over
   the convolved continuous field $f(g)$.

Convolution with an isotropic kernel of angular dispersion contracts orientation variance
toward isotropy. The exact relationship is governed by the kernel shrinkage coefficient $\rho$:

$$
\mathbf{A}_{\text{density}} = \frac{1}{3}\mathbf{I} + \beta\left(\mathbf{A}_{\text{support}} - \frac{1}{3}\mathbf{I}\right),
\qquad
\beta = \frac{3\rho - 1}{2},
$$

where $\rho = \langle \cos^2\omega \rangle_{\text{kernel}}$ is computed analytically by
`kernel_axis_shrinkage`. The parameter `deconvolve_kernel` governs this selection:
- `deconvolve_kernel=False` (default for fitted ODFs): Reports $\mathbf{A}_{\text{density}}$,
  reflecting the physical continuous distribution that reproduced experimental pole figures.
- `deconvolve_kernel=True` (for ODFs constructed from discrete orientations): Recovers
  $\mathbf{A}_{\text{support}}$, removing artificial kernel broadening.

## 3. Method selection guide

| Experimental Scenario | Recommended Route | Rationale |
| --- | --- | --- |
| Spatially resolved EBSD dataset | **Route A** | Direct summation, independent of diffraction peak overlap or defocus |
| Complete $(0002)$ pole figure with high intensity | **Route B** | Fast direct integration conforming to ASTM/standard industrial methods |
| Negligible $(0002)$ intensity in section plane | **Route C** | Basal orientation recovered via inversion of alternative reflections |
| Several incomplete pole figures | **Route C** | Multi-reflection ODF inversion constrains unmeasured regions |
| Non-standard principal direction analysis | **Any** | Evaluate full tensor $\mathbf{A}$ and diagonalize for principal axes |

## 4. Diagnostic report structure (`KearnsReport`)

The `KearnsReport` encapsulates full tensor diagnostics:

| Field | Interpretation |
| --- | --- |
| `f_triad` | Kearns parameters $[f_x, f_y, f_z]$ along requested axes, verifying $\sum f_i = 1.0$. |
| `tensor` | Full $3\times 3$ symmetric orientation tensor $\mathbf{A}$. |
| `eigenvalues`, `eigenvectors` | Principal Kearns parameters and corresponding unit axes. |
| `route` | Execution modality (`"orientations"`, `"pole_figure"`, or `"odf"`). |
| `deconvolve_kernel` | Status of kernel deconvolution during ODF integration. |

## Verification

- `tests/unit/test_kearns.py`: Verifies unit-trace invariant $\operatorname{tr}(\mathbf{A}) = 1$,
  isotropic limit $f = 1/3$, analytical kernel shrinkage matching numerical quadrature,
  and EBSD-versus-ODF consistency.
- Executable worked examples:
  - {doc}`../examples/generated/kearns-parameter`

## See also

- {doc}`../theory/kearns_parameter_and_basal_pole_texture` — Theoretical background on radiation growth and resolved fractions.
- {doc}`pole_figure_inversion` — Inversion of experimental pole figures into an ODF.
- {doc}`../theory/pole_figure_arithmetic_and_mrd` — Solid-angle integration and normalization conventions.

## References

### Normative

- Kearns, J. J. (1965). *Thermal Expansion and Preferred Orientation in Zircaloy*.
  Report WAPD-TM-472, Bettis Atomic Power Laboratory, Pittsburgh, PA.
- Baron, J. L., Esling, C., Feron, J. L., Gex, D., Glimois, J. L., Guillen, R.,
  Humbert, M., Lemoine, P., Lepape, J., Mardon, J. P., Thil, A. & Uny, G. (1990).
  Interlaboratories tests of textures of Zircaloy-4 tubes. Part 1: pole figure
  measurements and calculation of Kearns coefficients. *Textures and Microstructures*
  **12**, 125–140. <https://doi.org/10.1155/TSM.12.125>

### Informative

- Holt, R. A. & Aldridge, S. A. (1985). Effect of extrusion variables on
  crystallographic texture of Zr-2.5 wt% Nb. *Journal of Nuclear Materials* **135**,
  246–259. <https://doi.org/10.1016/0022-3115(85)90448-3>
- Mani Krishna, K. V., Srivastava, D., Dey, G. K., Hiwarkar, V., Samajdar, I. &
  Saibaba, N. (2011). Comparative study of methods of the determination of Kearns
  parameter in zirconium. *Journal of Nuclear Materials* **414**, 492–497.
  <https://doi.org/10.1016/j.jnucmat.2011.04.065>
