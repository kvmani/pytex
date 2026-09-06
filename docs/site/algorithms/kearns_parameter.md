# The Kearns Parameter: Three Routes To One Tensor

**Surface:** `pytex.texture.kearns.pole_orientation_tensor`,
`kearns_from_orientations`, `kearns_from_pole_figure`, `kearns_from_odf`,
`kernel_axis_shrinkage`, `KearnsReport`, with the workbench operations
`kearns.from_odf` and the Kearns panel.

The Kearns parameter $f$ is the resolved fraction of basal poles along a
specimen direction, and in zirconium alloys it governs irradiation growth,
hydride orientation, and anisotropic creep — which is why a cladding tube's
$f_{\text{RD}}$, $f_{\text{TD}}$, $f_{\text{ND}}$ appear on its specification.
It is also a quantity three different measurements claim to produce, and they
disagree in circumstances that are predictable. This page states the one object
underneath all three, the algorithm of each route, and which route is valid
when.

## 1. One tensor, not three numbers

Kearns' definition resolves the basal-pole distribution onto a direction:

$$
f_{d} \;=\; \bigl\langle \cos^{2}\alpha_{d} \bigr\rangle,
$$

with $\alpha_d$ the angle between a crystal's $[0001]$ and the specimen
direction $d$. Written out, $\cos\alpha_d = \mathbf{v}\cdot\mathbf{d}$ for a unit
basal pole $\mathbf{v}$, so

$$
f_d \;=\; \bigl\langle (\mathbf{v}\cdot\mathbf{d})^2 \bigr\rangle
      \;=\; \mathbf{d}^{\mathsf{T}}
            \underbrace{\bigl\langle \mathbf{v}\,\mathbf{v}^{\mathsf{T}}\bigr\rangle}
                       {}_{\textstyle \mathbf{A}}
            \mathbf{d}.
$$

**Every Kearns route is an estimate of the same second-moment tensor
$\mathbf{A}$**, and $f$ along any direction is its quadratic form. That
reframing is what `pole_orientation_tensor` implements, and it buys three
things immediately:

- **$\operatorname{tr}\mathbf{A} = 1$ identically**, so $f_{\text{RD}} +
  f_{\text{TD}} + f_{\text{ND}} = 1$ for any orthonormal triad. This is a
  *constraint*, not a measurement: three quoted values that fail it are
  arithmetically inconsistent, and PyTex reports the triple together so the sum
  can be seen.
- **The eigenvalues are the extreme Kearns parameters** and the eigenvectors
  their axes. A pilgered tube's basal maxima are not generally on RD/TD/ND, and
  the tensor finds where they are; three numbers on a fixed triad cannot.
- **$\mathbf{A}$ is invariant under $\mathbf{v} \to -\mathbf{v}$**, so antipodal
  data need no folding and a figure measured on one hemisphere is as good as one
  measured on two.

## 2. Route A — from discrete orientations

`kearns_from_orientations`. Available when orientations are known individually,
as from EBSD:

```text
1  resolve each orientation's (0001) pole (and its symmetry family) into the specimen frame
2  A = sum_i w_i v_i v_i^T / sum_i w_i     with v normalised
3  f_d = d^T A d for each requested direction
```

The most direct route, and the one with no quadrature and no inversion between
the measurement and the answer. Its limit is sampling: an EBSD map measures
grains at a surface, and a few thousand grains is a small sample of a texture
component.

## 3. Route B — from a measured pole figure

`kearns_from_pole_figure`, the standard industrial route, implementing Baron
*et al.* Eq. (5):

$$
f_{i} = \frac{\displaystyle\int\!\!\int P(\alpha,\beta)\cos^{2}\alpha_{i}\,\sin\alpha\,\mathrm{d}\alpha\,\mathrm{d}\beta}
             {\displaystyle\int\!\!\int P(\alpha,\beta)\,\sin\alpha\,\mathrm{d}\alpha\,\mathrm{d}\beta}
$$

which is exactly $\mathbf{d}^{\mathsf{T}}\mathbf{A}\mathbf{d}$ with
$\mathbf{A}$ the solid-angle-weighted second moment of the measured directions.

### 3.1 The weighting must follow the figure's sampling

The two kinds of pole figure need different weights, and using the wrong one
**over-counts the pole of a tilt raster by up to 50 per cent**:

| `sampling` | Rows are | Weight |
| --- | --- | --- |
| `scattered_poles` | individual poles carrying their own weight | the intensities *are* the weights; no quadrature |
| `sampled_density` | densities on a raster | intensity $\times$ solid angle of its cell |

PyTex reads this from the figure's own `sampling` attribute rather than
assuming, because the failure is silent: the answer stays plausible and is
systematically wrong toward the tilt axis.

### 3.2 The pseudo-norm, and the section where this route fails

Dividing by the *measured* integral rather than by $2\pi$ is the pseudo-norm of
Kern and Bergmann. It makes an **incomplete** figure usable — real
diffractometers cannot reach high tilt — at the price of assuming the unmeasured
cap has the same mean as the measured region.

This route is **not usable where the $(0002)$ peak is negligible**. In the ND-TD
section of strongly basal-textured zirconium the basal intensity is tiny, the
normalisation divides by noise, and $f$ becomes unstable. Mani Krishna *et al.*
(2011) traced inconsistent $f_{\text{RD}}$ values to exactly this, and the
resolution is route C rather than a better integration.

## 4. Route C — from a reconstructed ODF

`kearns_from_odf`. Invert several pole figures (see
{doc}`pole_figure_inversion`), then resolve the ODF's basal poles. Its advantage
is precisely the case that defeats route B: **alternative reflections plus the
inversion supply the basal information**, so no strong $(0002)$ peak is needed.

### 4.1 The kernel question, which is not a detail

A discrete ODF is a sum of kernels, and it contains **two distinguishable
objects with different Kearns parameters**:

- the **support tensor** $\mathbf{A}_{\text{support}}$ — the second moment of
  the support orientations' poles, which is what route A returns for the same
  weighted set;
- the **density tensor** $\mathbf{A}_{\text{density}}$ — the second moment under
  the continuous density $\sum_j w_j K(g; g_j)$ that the ODF object *is*.

Convolving with the kernel shrinks every departure from isotropy by a
closed-form factor, so they are related **exactly**:

$$
\mathbf{A}_{\text{density}}
  = \tfrac{1}{3}\mathbf{I} + \beta\left(\mathbf{A}_{\text{support}} - \tfrac{1}{3}\mathbf{I}\right),
\qquad
\beta = \frac{3\rho - 1}{2},
$$

with $\rho$ from `kernel_axis_shrinkage`. No numerical deconvolution is needed
in either direction.

**Which to report depends on where the ODF came from**, so it is the
`deconvolve_kernel` parameter rather than a hidden choice:

| ODF origin | Report | Why |
| --- | --- | --- |
| fitted by pole-figure inversion | density (default) | the weights were chosen so the *smoothed* density matches the measurement, so the smoothed density is the fitted quantity |
| built from known discrete orientations | support | the kernel is a display convenience there, not part of the measurement |

Getting this wrong biases $f$ toward $1/3$ — the isotropic value — by an amount
that grows with the kernel halfwidth, which looks like a weaker texture rather
than like an error.

## 5. Choosing a route

| Situation | Route |
| --- | --- |
| EBSD orientations in hand | A |
| A good $(0002)$ pole figure | B |
| $(0002)$ negligible in the measured section | **C** — B divides by noise |
| $f$ wanted along non-triad directions, or principal axes wanted | any; read the tensor, not three numbers |
| Several incomplete figures, no single good one | C |

## 6. What `KearnsReport` carries

The triple, the tensor, the eigenvalues and eigenvectors, the direction labels,
and the route with its parameters. The sum-to-one identity is checkable from the
report itself, which is the point: a Kearns triple that does not sum to one is
detectable without going back to the data.

## Verification

- The sum-to-one identity, the isotropic value $f = 1/3$, and the closed-form
  kernel shrinkage against direct quadrature, in
  {doc}`../examples/generated/kearns-parameter`.

## See also

- {doc}`../theory/kearns_parameter_and_basal_pole_texture` — the derivation,
  the tilt profile, and the volume-fraction reading.
- {doc}`pole_figure_inversion` — route C's first half.
- {doc}`../theory/pole_figure_arithmetic_and_mrd` — the m.r.d. scale and the
  solid-angle weights.

## References

### Normative

- J. J. Kearns, *Thermal Expansion and Preferred Orientation in Zircaloy*,
  WAPD-TM-472, Bettis Atomic Power Laboratory (November 1965). The defining
  report; his Eqs. (1)-(7) and Tables 2 and 3 are the pinned baselines.
- J. L. Baron, C. Esling, J. L. Feron, D. Gex, J. L. Glimois, R. Guillen,
  M. Humbert, P. Lemoine, J. Lepape, J. P. Mardon, A. Thil and G. Uny,
  *Interlaboratories tests of textures of Zircaloy-4 tubes. Part 1: pole figure
  measurements and calculation of Kearns coefficients*, Textures and
  Microstructures **12** (1990) 125-140.
  <https://doi.org/10.1155/TSM.12.125>. The pole-figure route (Eq. 5 above) and
  the incomplete-figure pseudo-norm.

### Informative

- K. V. Mani Krishna, D. Srivastava, G. K. Dey, V. Hiwarkar, I. Samajdar and
  N. Saibaba, *Comparative study of methods of the determination of Kearns
  parameter in zirconium*, Journal of Nuclear Materials **414** (2011) 492-497.
  <https://doi.org/10.1016/j.jnucmat.2011.04.065>. The cross-section dependence
  of the routes, and the ND-TD section where the pole-figure route fails.
- R. A. Holt and S. A. Aldridge, *Effect of extrusion variables on
  crystallographic texture of Zr-2.5 wt% Nb*, Journal of Nuclear Materials
  **135** (1985) 246-259.
  <https://doi.org/10.1016/0022-3115(85)90448-3>. The resolved-basal-pole form
  used throughout the pressure-tube literature.
