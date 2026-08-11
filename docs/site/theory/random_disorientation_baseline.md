# Random Disorientation And The Mackenzie Baseline

Every statement about a misorientation distribution is a statement about a *departure* from
randomness. "Grain boundaries cluster near 60°" means nothing until it is known what a texture-free
aggregate would have produced, because a texture-free aggregate does not produce a flat
distribution — it produces a strongly peaked one, for reasons that are entirely geometric.

This note derives that baseline, gives its closed forms, and shows what
`pytex.core.misorientation_distribution` computes against them. It also fixes a number that is
easy to quote wrongly: the mean and the median of the random cubic disorientation angle differ by
more than a degree and a half, and they are not interchangeable.

## The Unreduced Case Has An Exact Density

Take two orientations drawn independently from the uniform (Haar) distribution on $SO(3)$. Their
misorientation $\Delta g = g_{1}^{-1}g_{2}$ is then itself Haar-distributed, which is why a
baseline can be generated from *one* random rotation rather than two.

The Haar measure written in axis–angle coordinates factorises into a uniform distribution of axes
on the sphere and a non-uniform distribution of angles,

$$
\mathrm{d}\mu = \frac{1}{\pi^{2}}\,(1 - \cos\omega)\,\mathrm{d}\omega\,\mathrm{d}\Omega_{\hat{\mathbf{n}}},
$$ (eq-mdf-haar)

so with no symmetry at all the disorientation-angle density is exactly

$$
p(\omega) = \frac{1 - \cos\omega}{\pi},
\qquad 0 \le \omega \le \pi .
$$ (eq-mdf-triclinic)

The $(1-\cos\omega)$ factor is not a modelling choice; it is the volume of the shell of rotations
at angle $\omega$, and it is why large angles dominate. Integrating gives the mean

$$
\langle \omega \rangle
= \int_{0}^{\pi} \omega \, \frac{1-\cos\omega}{\pi} \, \mathrm{d}\omega
= \frac{\pi}{2} + \frac{2}{\pi}
= 126.4756^{\circ},
$$ (eq-mdf-triclinic-mean)

and the cumulative form $\int_{a}^{b} p = [\omega - \sin\omega]_{a}^{b}/\pi$ makes every band
probability a closed expression. The sampler reproduces it:

| Band | Sampled | Analytic {eq}`eq-mdf-triclinic` | Difference |
| --- | ---: | ---: | ---: |
| 0–30° | 0.00788 | 0.00751 | 0.00037 |
| 30–60° | 0.05055 | 0.05016 | 0.00039 |
| 60–90° | 0.12492 | 0.12402 | 0.00090 |
| 90–120° | 0.20897 | 0.20931 | 0.00034 |
| 120–150° | 0.28245 | 0.28318 | 0.00073 |
| 150–180° | 0.32523 | 0.32582 | 0.00059 |

with $n = 2\times10^{5}$ and a sampled mean of $126.36^{\circ}$ against the exact
$126.4756^{\circ}$. This is the strongest check available on the sampler, because the target is
analytic rather than another simulation.

## Symmetry Reduction, And Why Both Sides

With crystal symmetry the misorientation is no longer a single rotation but an equivalence class.
Both grains carry their own symmetry, so the disorientation is the minimum-angle representative
over the double coset

$$
\omega_{\mathrm{d}}(\Delta g)
= \min_{\mathbf{S}_{a}, \mathbf{S}_{b} \in G}
\ \angle\!\left( \mathbf{S}_{a} \, \Delta g \, \mathbf{S}_{b} \right),
$$ (eq-mdf-disorientation)

which is $|G|^{2}$ candidates per pair — 576 for cubic. PyTex evaluates them as one batched
`einsum` over the operator set rather than a nested loop.

Two points that are easy to get wrong:

- **The operators must be the proper rotations only.** The reduction is over rotations, and an
  improper operation is not one. PyTex uses the 24 rotations of the cubic group and 12 of the
  hexagonal, verified by determinant. Including the improper partners would over-reduce every
  angle and shift the whole baseline downward.
- **Grain-exchange symmetry does not need separate handling.** The disorientation is usually
  defined to include $\Delta g \mapsto \Delta g^{-1}$, but $\Delta g$ and $\Delta g^{-1}$ are
  transposes and a rotation matrix has the same trace as its transpose. The *angle* is therefore
  invariant under exchange, so omitting it changes the axis convention but never the distribution
  treated here.

## The Cubic Maximum Is Exact

The largest possible cubic disorientation is a property of the fundamental zone, not a sampling
outcome. In Rodrigues space $\boldsymbol{\rho} = \hat{\mathbf{n}}\tan(\omega/2)$ the cubic zone is
the intersection of a cube and an octahedron,

$$
|\rho_{i}| \le \sqrt{2} - 1 ,
\qquad
|\rho_{1}| + |\rho_{2}| + |\rho_{3}| \le 1 ,
$$ (eq-mdf-cubic-zone)

and the angle grows monotonically with $|\boldsymbol{\rho}|$, so the maximum sits at the zone
vertex farthest from the origin. That vertex is

$$
\boldsymbol{\rho}_{\max} = \bigl(\sqrt{2}-1,\ \sqrt{2}-1,\ 3-2\sqrt{2}\bigr),
$$ (eq-mdf-cubic-vertex)

which satisfies both constraints with equality — its components sum to exactly 1 and its two
largest equal $\sqrt{2}-1$. Its magnitude is $\sqrt{23 - 16\sqrt{2}}$, so

$$
\omega_{\max} = 2\arctan\sqrt{23 - 16\sqrt{2}} = 62.7994^{\circ},
$$ (eq-mdf-cubic-max)

about the axis $\langle 1,\ 1,\ \sqrt{2}-1 \rangle$, since dividing
{eq}`eq-mdf-cubic-vertex` through by $\sqrt{2}-1$ leaves
$\bigl(1, 1, \sqrt{2}-1\bigr)$. A sample of $3\times10^{5}$ random misorientations reaches
$62.56^{\circ}$; the shortfall is expected, because a maximum converges far more slowly than a
mean, and the exact value is the one to quote.

## Mean, Median, And A Number Worth Not Confusing

For random cubic orientations, over $3\times10^{5}$ samples:

| Quantity | Value |
| --- | ---: |
| mean $\langle \omega_{\mathrm{d}} \rangle$ | $40.71^{\circ}$ |
| median | $42.32^{\circ}$ |
| maximum (exact, {eq}`eq-mdf-cubic-max`) | $62.7994^{\circ}$ |

**The mean and the median differ by 1.6°** because symmetry reduction truncates the long
high-angle tail of {eq}`eq-mdf-triclinic` at $62.8^{\circ}$ and leaves a left-skewed distribution.
A single quoted "average disorientation of a random cubic aggregate" is therefore ambiguous, and
the two values are close enough to be interchanged by accident and far enough apart to matter when
a measured mean is being compared against the baseline.

The cubic mean was checked against an independently written quaternion implementation sharing no
code with PyTex: $40.749^{\circ}$ versus PyTex's $40.731^{\circ}$, a difference of one standard
error of the mean ($0.018^{\circ}$). Agreement between two implementations of the same definition
is weaker evidence than the analytic triclinic check above, but it rules out an error in the
double-coset reduction, which is the part with no closed form to test against.

## What The Baseline Is For

The practical use is as a null hypothesis, and the tail probabilities are more useful than the
mean. For random cubic orientations:

| Threshold | $P(\omega_{\mathrm{d}} < \theta)$ |
| --- | ---: |
| $2^{\circ}$ | $7.3\times10^{-5}$ |
| $5^{\circ}$ | $8.7\times10^{-4}$ |
| $10^{\circ}$ | $6.5\times10^{-3}$ |
| $15^{\circ}$ | $2.2\times10^{-2}$ |
| $> 60^{\circ}$ | $7.7\times10^{-3}$ |

The $15^{\circ}$ row is the one to remember: **about 2.2% of boundaries in a texture-free cubic
aggregate are low-angle by the conventional threshold, purely by chance.** A map reporting 3%
low-angle boundaries has demonstrated essentially nothing; one reporting 30% has. The same table
sets the floor for a Σ3 claim, since the $60^{\circ}/\langle 111 \rangle$ twin sits in the sparse
region above $60^{\circ}$ where the random background is under 1%.

## Why Sample Rather Than Evaluate The Closed Form

Mackenzie's distribution is piecewise analytic for the cubic case, so it could in principle be
evaluated directly. PyTex samples instead, for two reasons. The sampler is one implementation that
serves all 32 point groups, including the low-symmetry groups whose piecewise forms are not
tabulated anywhere convenient; and it extends without modification to a *correlated* baseline over
an explicit neighbour list, where no closed form exists at all. The closed forms in this note are
used as the check on the sampler rather than as the production path — which is the right way round,
since a check that shares code with the thing it checks is not a check.

## Cost And Limits

- **Memory is the binding constraint, not time.** {eq}`eq-mdf-disorientation` materialises
  $n|G|^{2}$ rotation matrices, or $72\,n$ doubles for cubic symmetry. At $n = 2\times10^{5}$ that
  is 7.7 GiB and the allocation fails. Generate a large baseline in chunks of order $2\times10^{4}$
  and concatenate; the samples are independent, so chunking changes nothing statistically.
- The baseline is seeded and therefore reproducible, but it is still Monte Carlo. Quote the exact
  {eq}`eq-mdf-cubic-max` for the maximum rather than a sampled one, and expect roughly
  $\sigma/\sqrt{n} \approx 11.3^{\circ}/\sqrt{n}$ on a cubic mean.
- The baseline assumes *uncorrelated* orientations. A real map has spatial correlation, so an
  uncorrelated MDF built from all unique pairs and a correlated MDF built from neighbour pairs
  answer different questions and should not be compared with each other.
- Specimen symmetry is not applied here; the reduction is by crystal symmetry on both sides.

## References

### Normative

- Th. Hahn (ed.), *International Tables for Crystallography, Volume A: Space-Group Symmetry*,
  IUCr / Springer. DOI: <https://doi.org/10.1107/97809553602060000100>. The point groups whose
  rotation subgroups enter {eq}`eq-mdf-disorientation`.

### Informative

- J. K. Mackenzie, *Second paper on statistics associated with the random disorientation of
  cubes*, Biometrika **45** (1958) 229–240. DOI: <https://doi.org/10.1093/biomet/45.1-2.229>.
  The original derivation of the random cubic disorientation distribution.
- A. Morawiec, *Orientations and Rotations: Computations in Crystallographic Textures*, Springer
  (2004). Rodrigues fundamental zones and the geometry behind {eq}`eq-mdf-cubic-zone`.
- V. Randle and O. Engler, *Introduction to Texture Analysis*, CRC Press. Boundary-character
  statistics and the conventional 15° low-angle threshold.

## See Also

- {doc}`orientation_space_and_disorientation` — the fundamental zone and the reduction itself.
- {doc}`orientation_representations` — the Rodrigues chart {eq}`eq-mdf-cubic-zone` is written in,
  and the same $(1-\cos\omega)$ measure derived from the equal-volume maps.
- {doc}`ebsd_local_misorientation` — where disorientation angles become KAM.
