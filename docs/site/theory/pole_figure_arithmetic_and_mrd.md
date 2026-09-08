# Pole-Figure Arithmetic And The m.r.d. Scale

A pole figure represents a probability density function defined over the two-dimensional unit
sphere $\mathbb{S}^2$. Because the spherical metric differs fundamentally from Euclidean planar
geometry, treating pole figures as planar images introduces systematic errors in normalization,
averaging, and residual calculation.

This note formalizes the integral definition of the multiples of a random distribution (m.r.d.)
scale, derives the solid-angle weighting required for equiangular diffractometer rasters, and
distinguishes true pole densities from signed residual difference fields.

## m.r.d. Is Defined By An Integral

Densities are reported in **multiples of a random distribution**: the value is 1 wherever the
distribution equals that of a texture-free, macroscopically isotropic aggregate. The definition is
governed by a spherical surface integral,

$$
\frac{1}{4\pi}\oint P_{hkl}(\mathbf{y}) \, \mathrm{d}\Omega = 1 ,
$$ (eq-pf-mrd)

and not by peak height or unweighted sample summation. A figure normalized so that its peak value
equals 1, or so that discrete point samples sum to 1, does not satisfy the m.r.d. metric and cannot
be directly compared with quantitative texture data. A uniform orientation distribution produces
uniformly distributed poles on $\mathbb{S}^2$, yielding a constant value of exactly 1 m.r.d. in
every direction and for every crystallographic plane family.

## Geometric Bias in Equiangular Spherical Rasters

A laboratory diffractometer typically samples a pole figure on an equiangular tilt/rotation grid:
a regular lattice in polar angle $\psi$ (sample tilt) and azimuth $\varphi$. This grid is not
areally uniform on the sphere. A latitude circle at polar angle $\psi$ has circumference
proportional to $\sin\psi$. Consequently, maintaining a constant angular step $\Delta\varphi$
crowds sampling points near the pole while dispersing them near the equator. The differential solid
angle element is

$$
\mathrm{d}\Omega = \sin\psi \, \mathrm{d}\psi \, \mathrm{d}\varphi ,
$$ (eq-pf-solid-angle)

so an unweighted mean over raster points is a mean with respect to $\mathrm{d}\psi\,\mathrm{d}\varphi$,
not $\mathrm{d}\Omega$, and it over-counts the pole.

The size of the error is not small and, crucially, **it does not shrink as the raster is refined**,
because it is a bias in the estimator rather than a discretisation error. Take the field
$f = \cos^{2}\psi$ over a hemisphere, whose two averages are both elementary:

$$
\langle f \rangle_{\mathrm{naive}}
= \frac{\displaystyle\int_{0}^{\pi/2}\!\cos^{2}\psi \, \mathrm{d}\psi}{\pi/2} = \frac{1}{2},
\qquad
\langle f \rangle_{\mathrm{sphere}}
= \frac{\displaystyle\int_{0}^{\pi/2}\!\cos^{2}\psi \, \sin\psi \, \mathrm{d}\psi}
       {\displaystyle\int_{0}^{\pi/2}\!\sin\psi \, \mathrm{d}\psi} = \frac{1}{3} .
$$ (eq-pf-bias)

The unweighted answer is exactly $3/2$ times the correct one — a **+50% error** — and the numbers
confirm it stays there:

| Raster step | Naive mean | Weighted mean | Exact | Naive error | Weighted error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| $5.0^{\circ}$ | 0.50000 | 0.31960 | $1/3$ | $+50.00\%$ | $-4.12\%$ |
| $2.5^{\circ}$ | 0.50000 | 0.32627 | $1/3$ | $+50.00\%$ | $-2.12\%$ |

Halving the step leaves the naive error at exactly 50% and halves the weighted one. This is the
signature of a bias: refinement cannot cure it, only weighting can.

`pytex.core.sphere.raster_solid_angle_weights` supplies the weights. Points are grouped into rings
of equal polar angle and each ring receives the solid angle of the band midway to its neighbours,
$\cos\psi_{\mathrm{lower}} - \cos\psi_{\mathrm{upper}}$, shared equally among its points. On a
$5^{\circ}$ raster the resulting per-point weight at the equator is about **92 times** the weight at
the pole, which is the factor a naive mean is silently applying as 1.

Two consequences worth stating:

- Grids built by `S2Grid` already carry their own weights. Applying raster weights on top would
  double-count; the weights belong to the sampling, not to the figure.
- The bands are clipped to the measured polar range, so a **partial** pole figure — the usual case,
  since defocusing limits the reachable tilt — is averaged over its measured cap. That equals the
  true spherical mean only if the unmeasured cap has the same mean, which is exactly the assumption
  a texture measurement cannot check. It should be stated rather than absorbed.

## Two Readings Of `intensities`, And Two Estimators

The same array can mean two different things, and PyTex records which rather than guessing:

- **`scattered_poles`** — the intensities are per-pole *weights* of a cloud of discrete poles. The
  underlying density is recovered by kernel density estimation: a weighted sum of kernels centred
  on the poles.
- **`sampled_density`** — the intensities are a density already *evaluated* at the given
  directions. Resampling it is interpolation.

These are not interchangeable, and applying the wrong estimator is a category error rather than an
approximation: density-estimating an already-smooth field broadens it by the kernel width a second
time, while interpolating a pole cloud returns a spiky field that depends on where the samples fell.
`PoleFigure.sampling` therefore carries the reading, and the resampling method follows from it.

## Normalisation of Kernel Density Estimates to the m.r.d. Scale

A raw kernel density estimate yields a smoothing response rather than a physical probability density
on the m.r.d. scale. The magnitude of the response depends directly on kernel bandwidth $\psi_{1/2}$.
To convert raw responses to normalized m.r.d., the estimate is divided by the response of a random
distribution evaluated with the identical kernel:

$$
P_{hkl}(\mathbf{y}) = \frac{\hat{P}_{\mathrm{KDE}}(\mathbf{y})}{P_{\mathrm{rand}}},
$$ (eq-pf-mrd-normalise)

where $P_{\mathrm{rand}}$ is evaluated via `random_pole_density(kernel)`. The normalization factor
varies substantially with kernel halfwidth:

| Kernel halfwidth | $P_{\mathrm{rand}}$ |
| ---: | ---: |
| $5^{\circ}$ | 33.80 |
| $10^{\circ}$ | 16.88 |
| $20^{\circ}$ | 8.39 |

Omitting this division scales the recovered field by bandwidth-dependent factors of $10^1$ to $10^2$.
Consequently, unnormalized distributions obtained with different kernel bandwidths cannot be
directly compared. Normalization to the m.r.d. scale ensures cross-dataset consistency and physical
comparability.

## Residual Pole Figures and Differential Densities

Subtracting two pole figures defined over a shared spherical support produces a signed difference field:

$$
\Delta P(\mathbf{y}) = P_{1}(\mathbf{y}) - P_{2}(\mathbf{y}) ,
$$ (eq-pf-difference)

where positive and negative excursions identify directions where $P_1$ exceeds or falls below $P_2$,
respectively. While a physical pole density is strictly non-negative and integrates to $4\pi$ sr
(normalized to 1 m.r.d.), the difference field $\Delta P$ satisfies:

$$
\frac{1}{4\pi}\oint \Delta P(\mathbf{y})\,\mathrm{d}\Omega = 0 .
$$

Consequently, $\Delta P$ represents a signed residual field rather than a valid pole density:
- It requires a divergent color scale centered at zero rather than a sequential m.r.d. colormap.
- It cannot be normalized to unit mean.
- It cannot serve as direct input to standard non-negative ODF inversion algorithms.

A primary application of differential fields is the **residual pole figure**, defined as the difference
between measured experimental intensities and back-calculated densities from a reconstructed ODF.
Residual fields serve as spatial diagnostic metrics: systematic patterns indicate unmodeled texture
components, sample misalignment, or instrument defocusing aberrations rather than random counting noise.

## Assumptions And Limits

- Weights describe the region actually measured. Nothing here reconstructs the unmeasured cap.
- The m.r.d. identity assumes the pole figure covers the sphere, or that the measured portion is
  representative. For a partial figure the normalisation inherits that assumption.
- Kernel density estimation smooths. A sharp texture measured with a broad kernel reports a lower
  peak than it has, and the halfwidth should be reported with any peak intensity.

## References

### Normative

- H.-J. Bunge, *Texture Analysis in Materials Science*, Butterworths (1982). Normalisation of the
  ODF and of pole figures to multiples of a random distribution.

### Informative

- V. Randle and O. Engler, *Introduction to Texture Analysis*, CRC Press. Measured pole-figure
  rasters, defocusing, and partial pole figures.
- R. Hielscher and H. Schaeben, *A novel pole figure inversion method: specification of the MTEX
  algorithm*, Journal of Applied Crystallography **41** (2008) 1024–1037. DOI:
  <https://doi.org/10.1107/S0021889808030112>. Kernel choice and its effect on recovered texture
  strength.

## See Also

- {doc}`discrete_odf_and_pole_figures` — construction and inversion of the figures treated here.
- {doc}`harmonic_odf_reconstruction` — the series alternative, where truncation replaces bandwidth.
- {doc}`/examples/generated/texture` — the uniform-ODF identity computed live.
