# Pole-Figure Arithmetic And The m.r.d. Scale

A pole figure is a density on a sphere, and almost every mistake made with one comes from treating
it as an image instead. Three of those mistakes are silent — they produce a plausible figure with
the wrong numbers on it — and each has an exact closed form showing how wrong.

This note covers the scale on which pole densities are reported, the weights that make an average
over measured data an integral, and what a difference of two pole figures is.

## m.r.d. Is Defined By An Integral

Densities are reported in **multiples of a random distribution**: the value is 1 wherever the
distribution is what a texture-free aggregate would give. The definition is a property of an
integral,

$$
\frac{1}{4\pi}\oint P_{hkl}(\mathbf{y}) \, \mathrm{d}\Omega = 1 ,
$$ (eq-pf-mrd)

and not of a maximum or of a sum. A figure normalised so its *peak* is 1, or so its samples *sum*
to 1, is on neither scale and cannot be compared with a published texture strength. The identity
that fixes it is that a uniform ODF sends poles uniformly over the sphere, so its pole figure is
flat at exactly 1 m.r.d. in every direction and for every plane family — which is the check worth
running whenever a normalisation is in doubt.

## The Raster Trap: A 50% Error That Finer Sampling Does Not Fix

A diffractometer samples a pole figure on a tilt/rotation raster: a regular grid in polar angle
$\psi$ and azimuth. Such a grid is **not** uniform on the sphere. A ring at polar angle $\psi$ has
circumference proportional to $\sin\psi$, so the same number of azimuthal points crowds into a
vanishing ring near the pole and spreads over the full circle at the equator. Solid angle goes as

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

## The Normalisation That Is Two Orders Of Magnitude

A kernel density estimate returns a *response*, not a density on the m.r.d. scale. Its size depends
on the kernel bandwidth, so it is not a physical quantity at all. The conversion divides by the
response a random texture produces,

$$
P_{hkl}(\mathbf{y}) = \frac{\hat{P}_{\mathrm{KDE}}(\mathbf{y})}{P_{\mathrm{rand}}},
$$ (eq-pf-mrd-normalise)

with $P_{\mathrm{rand}}$ from `random_pole_density(kernel)`. The factor is large and
bandwidth-dependent:

| Kernel halfwidth | $P_{\mathrm{rand}}$ |
| ---: | ---: |
| $5^{\circ}$ | 33.80 |
| $10^{\circ}$ | 16.88 |
| $20^{\circ}$ | 8.39 |

Skipping the division inflates every value by these factors — one to two orders of magnitude, not a
rounding matter. Worse, because the factor depends on the halfwidth, **two un-normalised figures
computed with different kernels are not comparable with each other**, and a single un-normalised
figure is not comparable with itself at a different smoothing. Normalisation is what makes the
number mean something outside the run that produced it.

## A Pole-Figure Difference Is Not A Pole Figure

Subtracting two pole figures on a shared support gives the signed field

$$
\Delta P(\mathbf{y}) = P_{1}(\mathbf{y}) - P_{2}(\mathbf{y}) ,
$$ (eq-pf-difference)

and the sign is the whole point: it is where one texture exceeds the other. But a density is
non-negative and integrates to 1 m.r.d., and $\Delta P$ does neither — by {eq}`eq-pf-mrd` its
spherical mean is **zero**, not one. It is therefore not a pole figure and should not be handed to
anything expecting one:

- an m.r.d. colour scale is wrong for it; a diverging scale centred on zero is right;
- normalising it to unit mean is meaningless, since its mean is zero by construction;
- it cannot be inverted to an ODF, because no ODF has negative pole density.

The important special case is the **residual pole figure**: the difference between a measurement
and the figure a fitted ODF recalculates for it. That is a goodness-of-fit check in spatial form,
and reading it as a texture is the mistake the type distinction exists to prevent — structure in a
residual is unmodelled texture or a systematic measurement error, never a physical density.

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
