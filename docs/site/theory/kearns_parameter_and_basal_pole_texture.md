# The Kearns Parameter And Basal-Pole Texture

The Kearns parameter $f$ is the most used scalar index of texture in the zirconium industry, and
it is quoted in component specifications, correlated with irradiation growth and creep, and
measured routinely by four different techniques that do not always agree. This note derives it,
shows that all four techniques estimate one and the same tensor, and works out exactly how much
each technique's approximations cost.

## Why One Number Can Replace An Orientation Distribution

An orientation distribution function carries everything about texture, which is precisely why it
is awkward for a structure–property correlation: it is a function on a three-dimensional
non-Euclidean space, and low crystal symmetry makes it worse. Kearns' observation (Kearns 1965)
is that for one important class of properties, almost all of that information is irrelevant.

A hexagonal single crystal is transversely isotropic about $[0001]$. Any second-rank property
tensor of such a crystal is therefore fixed by two numbers, its values parallel and perpendicular
to the $c$ axis, and its value along a direction at angle $\phi$ to $[0001]$ is

$$
P(\phi) = P_{\parallel}\cos^{2}\phi + P_{\perp}\bigl(1 - \cos^{2}\phi\bigr) .
$$ (eq-kearns-single-crystal)

Thermal expansion, irradiation growth, and the second-rank parts of elastic and creep response all
obey this. Now average {eq}`eq-kearns-single-crystal` over a polycrystal, taking the bulk property
in a reference direction to be the volume-weighted sum of the single-crystal contributions. With
$V_i$ the volume fraction of crystals whose $c$ axis lies at $\phi_i$ to that direction and
$\sum_i V_i = 1$,

$$
P_{\mathrm{ref}}
= P_{\parallel}\sum_i V_i \cos^{2}\phi_i
+ P_{\perp}\Bigl(1 - \sum_i V_i \cos^{2}\phi_i\Bigr) .
$$ (eq-kearns-polycrystal)

The whole orientation distribution has collapsed into the single summation, which Kearns names the
**orientation parameter**

$$
f = \sum_i V_i \cos^{2}\phi_i = \bigl\langle \cos^{2}\phi \bigr\rangle ,
\qquad
P_{\mathrm{ref}} = f\,P_{\parallel} + (1 - f)\,P_{\perp} .
$$ (eq-kearns-definition)

So $f$ is not merely *a* texture index: it is the exact and complete texture input to any property
obeying {eq}`eq-kearns-single-crystal`. That is what distinguishes it from the many other scalar
texture measures, and why it survived. Its interpretation — *the effective fraction of basal poles
aligned with the direction of interest* — follows from
{eq}`eq-kearns-polycrystal`: the aggregate behaves exactly as if a fraction $f$ of it were single
crystal with $c$ along the reference direction and the remaining $1 - f$ were single crystal with
$c$ perpendicular to it.

```{warning}
The reduction holds **only** for properties of the form {eq}`eq-kearns-single-crystal`. Yield
strength, fracture toughness, and hydride-orientation susceptibility are not second-rank
properties of this kind, and correlations of those with $f$ are empirical, not derived.
```

## The Tensor Behind It

Write $\mathbf{c}$ for the unit basal-pole direction of a crystal expressed in the specimen frame,
and $\mathbf{d}$ for a unit specimen direction. Then $\cos\phi = \mathbf{d}\cdot\mathbf{c}$ and
{eq}`eq-kearns-definition` becomes a quadratic form:

$$
f(\mathbf{d})
= \bigl\langle (\mathbf{d}\cdot\mathbf{c})^{2} \bigr\rangle
= \mathbf{d}^{\mathsf{T}}
  \underbrace{\bigl\langle \mathbf{c}\,\mathbf{c}^{\mathsf{T}} \bigr\rangle}_{\displaystyle \mathbf{A}}
  \mathbf{d} .
$$ (eq-kearns-tensor)

$\mathbf{A}$ is the second-moment, or **orientation**, tensor of the basal-pole direction
distribution — the same object directional statistics uses to classify a set of axes as clustered,
girdle-like, or uniform (see {doc}`directional_statistics_and_mean_axes`). Everything the
literature states about $f$ as separate empirical facts is a property of $\mathbf{A}$.

**The sum rule is exact.** For any orthonormal specimen triad $\mathbf{d}_1,\mathbf{d}_2,\mathbf{d}_3$,

$$
f_1 + f_2 + f_3 = \operatorname{tr}\mathbf{A}
= \bigl\langle \mathbf{c}\cdot\mathbf{c} \bigr\rangle = 1 ,
$$ (eq-kearns-sum-rule)

because every $\mathbf{c}$ is a unit vector. This holds for *every* texture, with no approximation
whatsoever. It follows that **a measured triad that does not sum to 1 is reporting the systematic
error of the measurement, not a property of the material** — and that is the single most useful
diagnostic available, because it needs no reference specimen.

**The random value is $1/3$.** A random texture sends $\mathbf{c}$ uniformly over the sphere, so
$\mathbf{A} = \mathbf{I}/3$ and $f = 1/3$ in every direction. Combined with the sum rule this fixes
the whole scale: $f = 0$ means no basal poles along the direction, $f = 1/3$ means random, $f = 1$
means all of them.

**Every direction, not only three.** Because {eq}`eq-kearns-tensor` is a quadratic form, one
tensor answers for all directions. Its eigenvalues are the Kearns parameters along the texture's
own principal axes and bound every other direction's value between the smallest and largest. This
matters for pilgered tubing, whose basal maxima sit at $\pm 20$–$40^{\circ}$ from the radial
direction in the $R$–$T$ plane (Baron *et al.* 1990) and therefore *not* on the axes $f_R$, $f_T$,
$f_L$ are quoted along.

In PyTex, {func}`~pytex.texture.kearns.pole_orientation_tensor` computes $\mathbf{A}$ and every
`kearns_from_*` function returns it on the {class}`~pytex.texture.kearns.KearnsReport` when the
route determines it.

## The Tilt Profile: Why Azimuth Does Not Matter

Kearns' second observation is that the azimuthal detail of a pole figure is irrelevant to $f$.
Only $\phi$ enters {eq}`eq-kearns-definition`, so the reference direction may be treated as a
fibre axis and the pole density averaged over the full $360^{\circ}$ of rotation about it. Write
$I(\phi)$ for that azimuthal average. The volume fraction in a band $\mathrm{d}\phi$ is the density
times the band's solid angle, and the band at $\phi$ has circumference proportional to $\sin\phi$:

$$
\mathrm{d}V \propto I(\phi)\,\sin\phi\,\mathrm{d}\phi ,
\qquad
f = \frac{\displaystyle\int_{0}^{\pi/2} I(\phi)\,\sin\phi\,\cos^{2}\phi\,\mathrm{d}\phi}
         {\displaystyle\int_{0}^{\pi/2} I(\phi)\,\sin\phi\,\mathrm{d}\phi} .
$$ (eq-kearns-tilt-integral)

This is Kearns' Equation (5), implemented as
{func}`~pytex.texture.kearns.kearns_from_tilt_profile`. Two consequences of the $\sin\phi$ factor
are worth stating because they are routinely misread off a pole figure:

- **The volume fraction vanishes at the centre** however intense the pole is there. A bright spot
  at $\phi = 0$ contributes nothing, because the band it occupies has no area.
- **High-tilt data matter more than they look.** The weight rises to a maximum near
  $\phi = 55^{\circ}$ and stays large to $90^{\circ}$, exactly where a reflection measurement is
  weakest. Truncating at $75^{\circ}$ discards the region that pulls $f$ down, so a truncated
  figure reports $f$ too high along the section normal.

Because {eq}`eq-kearns-tilt-integral` is a *ratio*, the scale of $I$ cancels: times-random units,
counts per second, or arbitrary units all give the same $f$. What does not cancel is any
$\phi$-dependent distortion — defocusing, absorption, background — which is why those must be
corrected before integration and not after.

## The Four Routes And What Each Assumes

### Discrete orientations (EBSD, simulation)

Evaluate {eq}`eq-kearns-tensor` directly on the measured orientations. No binning, no
interpolation, no truncation, and no normalization assumption enter; the only error is orientation
statistics. Mani Krishna *et al.* (Mani Krishna *et al.* 2011) found this the most consistent route for
recrystallized microstructures, with $f$ identical across all three principal sections — and also
found its weakness, that heavily deformed material indexes poorly, so the orientations that fail
to index are exactly the ones that are differently oriented.

Implemented by {func}`~pytex.texture.kearns.kearns_from_orientations`. Grain-area or volume
weights belong here: an unweighted mean over indexed points is a mean over *area*, which is the
right measure only if grain size does not correlate with orientation.

### The pole-figure route

Baron *et al.* (Baron *et al.* 1990) define the Kearns coefficients as

$$
f_{i} = \frac{\displaystyle\int_{0}^{\alpha_{\max}}\!\!\int_{0}^{2\pi}
              P(\alpha,\beta)\,\cos^{2}\alpha_{i}\,\sin\alpha\,\mathrm{d}\alpha\,\mathrm{d}\beta}
             {\displaystyle\int_{0}^{\alpha_{\max}}\!\!\int_{0}^{2\pi}
              P(\alpha,\beta)\,\sin\alpha\,\mathrm{d}\alpha\,\mathrm{d}\beta} ,
$$ (eq-kearns-pole-figure)

which is {eq}`eq-kearns-tensor` written out with $\alpha,\beta$ the polar and azimuthal
coordinates of the pole figure and $\alpha_i$ the angle from specimen direction $i$. Dividing by
the *measured* integral rather than by $2\pi$ is the pseudo-norm of Kern and Bergmann: it makes an
incomplete figure usable, at the price of assuming the unmeasured cap has the same mean as the
measured one.

Three systematic errors dominate this route.

1. **Truncation.** Reflection geometry defocuses beyond about $75$–$80^{\circ}$ of tilt, so
   $\alpha_{\max} < 90^{\circ}$. The measured cap covers $1 - \cos\alpha_{\max}$ of the
   hemisphere — only $74\%$ at $75^{\circ}$ — and, by the $\sin\phi$ argument above, the missing
   part is where the weight is largest. The mitigation in industrial practice is to measure two
   perpendicular sections and average, which (Mani Krishna *et al.* 2011) shows brings the pole-figure
   route into agreement with EBSD.
2. **Sections with no basal intensity.** In strongly basal-textured zirconium the $(0002)$ peak in
   the ND–TD section is negligible, so {eq}`eq-kearns-pole-figure` divides by noise. This is the
   traced cause of the inconsistent $f_{\mathrm{RD}}$ values from that section
   (Mani Krishna *et al.* 2011); the ODF route is the substitute, not a better pole-figure correction.
3. **Raster quadrature.** A tilt/rotation raster is not a uniform sampling of the sphere, and an
   unweighted mean over it over-counts the centre by up to $50\%$ — see
   {doc}`pole_figure_arithmetic_and_mrd`. The correct weights are
   {func}`~pytex.core.sphere.raster_solid_angle_weights`.

  There is a subtler quadrature point specific to $f$. A hemispherical raster ending exactly at
  $\alpha = 90^{\circ}$ has an outermost ring whose band, if extended outwards by its own half
  step, runs past the equator and claims close to twice the solid angle it owns. Since
  $\cos^{2}\alpha = 0$ there, that ring's excess weight pulls $f$ down: on a $5^{\circ}$ raster the
  spherical mean of $\cos^{2}$ comes out at $0.3196$ instead of $1/3$, a $-4.1\%$ error, which
  bounding the band at $90^{\circ}$ reduces to $-0.06\%$. {func}`~pytex.texture.kearns.kearns_from_pole_figure`
  passes `polar_max_deg=90.0` for antipodal figures for this reason.

Implemented by {func}`~pytex.texture.kearns.kearns_from_pole_figure`.

### The ODF route

Reconstruct an ODF — typically by inverting several incomplete pole figures — and take
$\mathbf{A}$ from it. Its advantage is that it does not need a strong $(0002)$ peak in the measured
section: alternative reflections plus the inversion supply the basal information, which is why it
remains usable where the pole-figure route fails (Mani Krishna *et al.* 2011). Its cost is that
pole-figure inversion is ill-posed, so the result inherits the regularization, the positivity
correction, and the ghost problem (see {doc}`ghost_problem_and_odd_harmonics`).

A kernel-density ODF carries one further bias that is worth deriving, because it is exactly
computable and is otherwise mistaken for material behaviour.

::::{admonition} The kernel shrinkage factor
:class: note

A fitted ODF is a smoothed estimate. Convolving with an isotropic SO(3) kernel spreads each
crystal's basal pole over a cone, pulling $\mathbf{A}$ toward isotropy. Let a smeared pole be
$\mathbf{v} = \mathbf{R}\mathbf{c}$, with $\mathbf{R}$ a rotation of angle $\omega$ about an axis
$\mathbf{a}$ uniform on the sphere. Rodrigues' formula gives

$$
\mathbf{c}\cdot\mathbf{R}\mathbf{c} = \cos\omega + t\,(1 - \cos\omega),
\qquad t = (\mathbf{a}\cdot\mathbf{c})^{2} ,
$$ (eq-kearns-rodrigues)

and for a uniform axis $\mathbb{E}[t] = 1/3$, $\mathbb{E}[t^{2}] = 1/5$. Averaging first over the
axis and then over the kernel's angular density on SO(3) — the kernel value times the Haar factor
$(1 - \cos\omega)/\pi$ — gives

$$
\rho = \bigl\langle \cos^{2}\beta \bigr\rangle
     = \Bigl\langle \cos^{2}\omega + \tfrac{2}{3}\cos\omega\,(1-\cos\omega)
       + \tfrac{1}{5}(1-\cos\omega)^{2} \Bigr\rangle_{K} ,
$$ (eq-kearns-rho)

and, since the smeared distribution is rotationally symmetric about each original pole,

$$
\mathbf{A}_{\text{density}}
= \tfrac{1}{3}\mathbf{I} + \beta\left(\mathbf{A}_{\text{support}} - \tfrac{1}{3}\mathbf{I}\right),
\qquad
\beta = \frac{3\rho - 1}{2} .
$$ (eq-kearns-shrinkage)

Every departure from $1/3$ is scaled by $\beta$, which is $0.98$ at a $5^{\circ}$ halfwidth,
$0.937$ at $10^{\circ}$, $0.775$ at $20^{\circ}$ and $0.577$ at $30^{\circ}$ for the de la Vallée
Poussin kernel. At a $20^{\circ}$ halfwidth an $f$ of $0.70$ reads as $0.62$: the difference is the
kernel, not the material. The result is independent of crystal symmetry, because conjugating the
kernel by a symmetry operation leaves it isotropic, so every symmetry branch shrinks the tensor by
the same factor.
::::

{func}`~pytex.texture.kearns.kernel_axis_shrinkage` returns $\rho$ and
{func}`~pytex.texture.kearns.kearns_from_odf` reports either reading. Which is correct depends on
provenance: for an ODF fitted by pole-figure inversion the weights were chosen so that the
*smoothed* density matches the data, so the smoothed density is the model of the material; for an
ODF built from measured EBSD orientations the support *is* the data and the kernel is estimation
blur. **Report the halfwidth with any $f$ taken from an ODF.**

### The diffractogram (inverse pole figure) route

This is Kearns' own method, and the only one needing no texture goniometer. In a symmetric
$\theta$–$2\theta$ scan only planes parallel to the specimen surface diffract, so the integrated
intensity of the $(hkil)$ peak, compared with the same peak from a random powder, is the pole
density of $(hkil)$ along the section normal in times-random units. Crystals contributing to that
peak have their $(hkil)$ normal along the section normal, and therefore their basal pole at the
fixed angle

$$
\cos\phi_{hkil} =
\frac{\ell/c}{\sqrt{\dfrac{4\left(h^{2} + hk + k^{2}\right)}{3a^{2}} + \dfrac{\ell^{2}}{c^{2}}}}
$$ (eq-kearns-tilt-angle)

to it. Each reflection therefore supplies one point of $I(\phi)$, and the set of them *is* an
inverse pole figure of the section normal — hence the name. In PyTex
{func}`~pytex.texture.kearns.basal_tilt_angle_deg` evaluates {eq}`eq-kearns-tilt-angle` from the
phase's own reciprocal metric rather than from a transcribed table, so it follows the lattice
parameters actually in use. For $\alpha$-Zr at $c/a = 1.593$ it reproduces Kearns' hand-tabulated
values: $(10\bar{1}1)$ at $61.5^{\circ}$ against his $61.4$, $(10\bar{1}2)$ at $42.6$ against
$42.5$, $(20\bar{2}1)$ at $74.8$ against $74.8$.

Four approximations enter, in decreasing order of consequence.

1. **Interpolation across gaps.** The available reflections are not uniformly spread over $\phi$;
   for $\alpha$-Zr there is nothing between $0^{\circ}$ and $20.2^{\circ}$. Assuming a linear
   variation of intensity across that gap — the assumption (Mani Krishna *et al.* 2011) identifies as
   the route's central approximation — over-estimates $I$ at low $\phi$ for a peaked texture and
   therefore over-estimates $f$. This, and not the quadrature, is why the three sections' values
   do not sum to 1 and are normalized before use.
2. **The reference intensities.** Raw peak areas are not pole densities: in a *random* powder
   $(10ar{1}1)$ and $(20ar{2}0)$ differ by a factor of twenty from structure factor and
   multiplicity alone, so every reflection needs an $I_{0}$. It may be measured on a powder
   standard or calculated from the structure — only ratios enter $f$, so arbitrary units serve.
   The classical Harris texture coefficient rescales $I/I_{0}$ to a mean of one over the measured
   reflections; Kearns tested that assumption against his own standard and found the mean ran from
   $1.02$ to $1.56$, averaging $1.23$, so it understates absolute pole densities by about $23\%$.
   It understates $f$ by nothing at all, because $f$ is a ratio of two integrals over the same
   profile and a common scale factor cancels identically. Anyone reaching for the normalization to
   repair a triad that misses 1 is reaching for the wrong tool.
3. **Three sections, three specimens.** A triad needs all three principal sections, which for thin
   product forms such as clad tubing means stacking slices and correcting for the packing. Where
   the material is inhomogeneous, the three values then describe three different regions.
4. **Symmetric geometry.** The derivation assumes only surface-parallel planes diffract. A
   **fixed-$\omega$ detector scan** violates this: the diffraction vector sits at $\theta - \omega$
   from the surface normal and moves through the pattern, so different reflections probe different
   specimen directions. On a four-circle instrument this is easy to do by accident — the data file
   records `scanAxis="2Theta"` with a single `Omega` position, against `scanAxis="Gonio"` with both
   axes ranged for a coupled scan — and at $\omega = 15^{\circ}$ over a $20$–$120^{\circ}$ range
   the diffraction vector swings from $2.5^{\circ}$ to $45^{\circ}$ off the normal.
   {class}`~pytex.texture.kearns.DiffractogramReflection` carries a `specimen_tilt_deg` per
   reflection and {func}`~pytex.texture.kearns.kearns_from_diffractogram` reports the spread, so
   the condition is visible rather than silent.

## Validation

The identities of {eq}`eq-kearns-sum-rule` and the random value are exact, which makes them
calibrations rather than tolerances. `tests/unit/test_kearns_parameter.py` pins:

| Case | Expected | Provenance |
| --- | --- | --- |
| Random orientations | $f = 1/3$ in every direction | Isotropy plus {eq}`eq-kearns-sum-rule` |
| Any texture, any orthonormal triad | $\sum f = 1$ to $10^{-12}$ | {eq}`eq-kearns-sum-rule` |
| Single crystal | $f = (0, 0, 1)$ | {eq}`eq-kearns-definition` |
| Ideal basal girdle in the RD–TD plane | $f = (1/2, 1/2, 0)$ | $\langle\cos^{2}\rangle$ over a great circle |
| Kearns' Table 3, longitudinal section | $f = 0.488$ | Kearns (1965) Table 3 |
| Tilt angles $\phi_{hkil}$ for $\alpha$-Zr | Kearns Table 2, to $0.2^{\circ}$ | (Kearns 1965) Table 2 |
| Kernel shrinkage $\rho$ | Monte-Carlo integral over SO(3) | {eq}`eq-kearns-rho` |
| Diffractogram route on a simulated fibre | true $f$ to $\pm 0.02$ | End-to-end simulation |

```{admonition} An arithmetic slip in the source of record
:class: caution

Kearns' Table 3 reproduces exactly for the longitudinal section: his tabulated intensities give
$f = 0.4879$ against his quoted $0.488$. The transverse-section block does not. Its $70$–$80^{\circ}$
row lists $V_{\Delta\phi}\cos^{2}\bar{\phi} = 0.0214$ where
$0.353 \times \cos^{2}(75^{\circ}) = 0.0237$, and the quoted total $f = 0.0508$ carries the error;
recomputing from his own $I_{\phi}$ column gives $0.0526$. The longitudinal block is therefore what
this repository pins as a regression baseline.
```

## Assumptions And Limits

- $f$ is complete only for properties obeying {eq}`eq-kearns-single-crystal`. It is an index, not a
  derivation, for anything else.
- A single tilt profile determines one direction's value. The tensor, and hence the sum rule,
  needs one profile per section.
- Nothing here reconstructs an unmeasured cap. The pseudo-norm assumes it resembles what was
  measured, and that assumption is recorded in the report's notes rather than hidden.
- The machinery works for any pole, not only $(0001)$; a non-basal result is a legitimate
  resolved-pole fraction but is not *the* Kearns parameter, so the pole travels on the report.

## References

### Normative

- J. J. Kearns, *Thermal Expansion and Preferred Orientation in Zircaloy*, WAPD-TM-472, Bettis
  Atomic Power Laboratory (November 1965). The defining report: Eqs. (1)–(7) here correspond to his
  (1)–(7), and his Tables 2 and 3 are the pinned baselines.
- J. L. Baron, C. Esling, J. L. Feron, D. Gex, J. L. Glimois, R. Guillen, M. Humbert, P. Lemoine,
  J. Lepape, J. P. Mardon, A. Thil and G. Uny, *Interlaboratories tests of textures of Zircaloy-4
  tubes. Part 1: pole figure measurements and calculation of Kearns coefficients*, Textures and
  Microstructures **12** (1990) 125–140. DOI: <https://doi.org/10.1155/TSM.12.125>. The
  pole-figure route and the incomplete-figure pseudo-norm, with an interlaboratory scatter of
  $4$–$10\%$ on $f_R$ and $f_T$ and $40\%$ on the small $f_L$.

### Informative

- R. A. Holt and S. A. Aldridge, *Effect of extrusion variables on crystallographic texture of
  Zr-2.5 wt% Nb*, Journal of Nuclear Materials **135** (1985) 246–259. DOI:
  <https://doi.org/10.1016/0022-3115(85)90448-3>. The resolved-basal-pole form
  $F_d = \sum V(\theta)\cos^{2}\theta$ as used throughout the CANDU pressure-tube literature.
- K. V. Mani Krishna, D. Srivastava, G. K. Dey, V. Hiwarkar, I. Samajdar and N. Saibaba,
  *Comparative study of methods of the determination of Kearns parameter in zirconium*, Journal of
  Nuclear Materials **414** (2011) 492–497. DOI:
  <https://doi.org/10.1016/j.jnucmat.2011.04.065>. Cross-section dependence of the four routes and
  the normalization the diffractogram route needs.
- K. Linga Murty and I. Charit, *Texture development and anisotropic deformation of zircaloys*,
  Progress in Nuclear Energy **48** (2006) 325–359. DOI:
  <https://doi.org/10.1016/j.pnucene.2005.09.011>. Why $f$ is specified for reactor components.

## See Also

- {doc}`pole_figure_arithmetic_and_mrd` — the m.r.d. scale and raster quadrature the pole-figure
  route depends on.
- {doc}`directional_statistics_and_mean_axes` — the orientation tensor as a tool of directional
  statistics.
- {doc}`hexagonal_conventions` — four-index notation and the hexagonal metric behind
  {eq}`eq-kearns-tilt-angle`.
- {doc}`discrete_odf_and_pole_figures` — the ODF and pole-figure objects the routes consume.
- {doc}`/tutorials/notebooks/31_kearns_parameter` — the four routes worked end to end on simulated
  and measured zirconium textures.
- {doc}`/examples/generated/texture` — the identities above computed live.
