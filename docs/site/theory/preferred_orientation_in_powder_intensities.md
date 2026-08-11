# Preferred Orientation In Powder Intensities

This note fixes the two preferred-orientation corrections PyTex applies to
powder diffraction intensities, and states exactly what each assumes. It is the
theory behind `pytex.diffraction.preferred_orientation`:
`march_dollase_factors`, `MarchDollaseModel`,
`ODFPreferredOrientationModel`, and the
`preferred_orientation` argument of `generate_xrd_pattern`.

## The Problem

The relative intensities of a powder pattern are derived on the assumption that
the crystallites are randomly oriented, so that the fraction of them presenting
any given plane to the beam depends only on the plane's multiplicity. Real
specimens violate this routinely: platy or acicular powders align on packing,
rolled and drawn material carries deformation texture, and thin films grow with
a fibre axis.

Texture changes intensities without changing peak positions. Uncorrected, that
is easily misread as a wrong structure, a wrong phase fraction, or a wrong site
occupancy — which is why a correction is not optional for quantitative work.

## March–Dollase

The model of March, as adapted to diffraction by Dollase, describes a single
fibre texture with one parameter. For a reflection whose plane normal makes an
angle $\alpha$ with the preferred-orientation axis,

$$
P(\alpha) = \left( r^{2}\cos^{2}\alpha
                     + \frac{\sin^{2}\alpha}{r} \right)^{-3/2}
$$

where $r$ is the March coefficient. Its limits are

$$
P(0) = r^{-3}, \qquad P(\pi/2) = r^{3/2}
$$

so $r < 1$ describes a plate-like habit whose plate normals cluster along the
axis and enhances those reflections, while $r > 1$ describes a needle-like habit
and suppresses them. At $r = 1$ the bracket is $\cos^2\alpha + \sin^2\alpha$ and
the correction is the identity.

**Exact normalization.**

Averaged over a uniform distribution of directions the March function is exactly
1, for every $r$. With $u = \cos\alpha$, $A = r^{2} - r^{-1}$ and $B = r^{-1}$,

$$
\frac{1}{2}\int_{-1}^{1} \bigl(A u^{2} + B\bigr)^{-3/2} \, du
    = \left[ \frac{u}{B\sqrt{A u^{2} + B}} \right]_{0}^{1}
    = \frac{1}{B\sqrt{A + B}}
    = \frac{r}{\sqrt{r^{2}}}
    = 1
$$

This is the statement that preferred orientation *redistributes* diffracted
intensity rather than creating it, and it is what makes a fitted $r$ a
description of texture rather than a free intensity scale. PyTex pins the
identity as a test by equal-area quadrature on the sphere, and as a worked
example by direct integration.

**Family averaging.**

A powder reflection is a symmetry family, not a single plane: every equivalent
$(hkl)$ diffracts at the same Bragg angle, and each sits at its own angle to the
preferred-orientation axis. PyTex therefore averages the March function over the
whole family,

$$
P_{hkl} = \frac{1}{m}\sum_{i=1}^{m} P(\alpha_i)
$$

with $\alpha_i$ the angle between the $i$-th equivalent normal and the axis,
taken to the nearer pole so that a plane and its opposite normal are treated as
one. Without this averaging the correction would depend on which family
representative the reflection enumeration happened to emit, which is an artefact
of the code rather than a property of the specimen.

**What it assumes.**

Rotational symmetry about the specimen axis. A sheet texture with distinct
rolling and transverse behaviour violates that assumption, and March–Dollase
will absorb the discrepancy into a fitted $r$ with no physical meaning. PyTex
does not check the assumption — it cannot, from a powder pattern alone — but
the `describe()` output states it.

## ODF Weighting

The physically direct correction reads the texture instead of parameterizing it.
The intensity of a powder reflection is proportional to the density of $\{hkl\}$
poles lying along the scattering vector, and an orientation distribution
function supplies exactly that quantity:

$$
P_{hkl} = P_{hkl}\bigl(\hat{\mathbf{y}}\bigr)
$$

the $\{hkl\}$ pole density evaluated along the scattering direction
$\hat{\mathbf{y}}$, in multiples of a random distribution. There is no fitted
parameter and no assumption of fibre symmetry: an arbitrarily complex texture is
handled as faithfully as the ODF represents it.

**Normalization.**

`ODF.evaluate_pole_density` returns a kernel-weighted response, not a
value in multiples of random: the smoothing kernel used there peaks at 1 rather
than integrating to 1, so a *uniform* texture yields the kernel's spherical
mean, not unity. PyTex divides the response by that mean,

$$
c = \frac{1}{2}\int_{-1}^{1} k\bigl(\arccos u\bigr)\, du
$$

evaluated by Gauss–Legendre quadrature in $u = \cos\omega$, where the integrand
is smooth. This is what makes the correction reduce to exactly 1 for an
untextured specimen, which is the limiting case that gives the factors their
meaning.

**Geometry.**

In symmetric Bragg–Brentano reflection geometry the scattering vector lies along
the specimen normal at every angle, which is why the default scattering
direction is ND. That is an approximation in asymmetric and transmission
geometries, where the scattering direction moves with $2\theta$; the direction is
an explicit argument so it can be set correctly there.

## Choosing Between Them

- Use **March–Dollase** when preferred orientation is a nuisance to be fitted away, no texture measurement exists, and the specimen plausibly has one fibre axis.
- Use **ODF weighting** when a texture measurement exists, when the texture is not a simple fibre, or when the intensities themselves are the quantity of interest rather than a means to a structure.

Both are exposed through the same protocol, so a pattern simulation accepts
either without knowing which.

## Current Limits

- The corrections scale kinematic intensities; absorption, extinction, and surface roughness are not modelled, and those also alter powder intensities.
- March–Dollase is a single-axis model. Multi-component textures need the ODF route.
- The ODF correction is only as good as the ODF: an under-determined or over-smoothed reconstruction flattens real texture towards 1 and under-corrects. The reconstruction residuals should be checked first.
- Corrected reflection lists are deliberately not renormalized, because the ratios between reflections are the texture information.

## References

March, A., *Z.\ Kristallogr.* **81**, 285–297 (1932).

Dollase, W. A., *J.\ Appl.\ Cryst.* **19**, 267–272 (1986),
`DOI: 10.1107/S0021889886089458`.

Von Dreele, R. B., *J.\ Appl.\ Cryst.* **30**, 517–525 (1997),
`DOI: 10.1107/S0021889897005918`.

Bunge, H.-J., *Texture Analysis in Materials Science*,
`DOI: 10.1016/C2013-0-11769-2`.
