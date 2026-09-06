# The Ghost Problem: What Pole Figures Cannot Determine

Diffraction pole figures do not determine the orientation distribution function. This is not a
precision limitation to be beaten with better counting statistics — it is an exact structural
degeneracy, it removes roughly **half** the degrees of freedom of the ODF, and it has been known
since the 1970s under the name *ghost problem*.

PyTex's harmonic reconstruction defaults to `even_degrees_only=True` in consequence. This note says
why that default exists, how much information is actually missing, and what does and does not
recover it.

## The Cause: A Pole Figure Is Centrosymmetric

A pole figure records the density of a plane normal $\mathbf{h}$ along specimen directions
$\mathbf{y}$. Two facts collapse it:

- a lattice plane has no sense — $(hkl)$ and $(\bar{h}\bar{k}\bar{l})$ are the same plane, so the
  normal enters as an axis;
- Friedel's law makes the diffracted intensity from $\mathbf{h}$ and $-\mathbf{h}$ identical in the
  absence of anomalous scattering.

So the measured quantity satisfies

$$
P_{\mathbf{h}}(\mathbf{y}) = P_{\mathbf{h}}(-\mathbf{y}),
$$ (eq-ghost-centro)

**whatever the ODF is**. The centrosymmetry is a property of the measurement, not of the specimen.

This is directly visible in the pole sets PyTex builds. Taking a deliberately one-sided orientation
population — 200 orientations with $\varphi_{1}$ confined to $0$–$40^{\circ}$, nothing about it
symmetric — and generating its $\{111\}$ pole figure gives 1600 poles of which **every single one
has its antipode also present** (checked for 300 of them, 300/300). An asymmetric ODF produces a
symmetric pole figure, and the asymmetry is gone.

## The Consequence: Odd Harmonics Are Unobservable

Writing the ODF in the generalized spherical harmonic series
$f(g) = \sum_{\ell} \sum_{m n} c_{\ell}^{mn} T_{\ell}^{mn}(g)$, the pole figure is a linear
projection of $f$ — an integral over the great circle of orientations placing $\mathbf{h}$ along
$\mathbf{y}$. Under {eq}`eq-ghost-centro` that projection annihilates every **odd** degree
$\ell$:

$$
P_{\mathbf{h}}(\mathbf{y}) \ \text{depends only on} \ \left\{ c_{\ell}^{mn} : \ell \ \text{even} \right\} .
$$ (eq-ghost-even)

The odd coefficients are not poorly determined; they are not determined at all. Any two ODFs
differing only in their odd part produce **identical** pole figures, so no amount of data, no
regularization, and no better solver can separate them.

How much is that? Counting basis terms before symmetry projection:

| Bandlimit $L$ | Even-degree terms | All terms | Discarded | Fraction |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 1125 | 1929 | 804 | 41.7% |
| 16 | 7113 | 13073 | 5960 | 45.6% |
| 22 | 17260 | 32407 | 15147 | 46.7% |
| 32 | 50065 | 95777 | 45712 | 47.7% |

Approaching one half, as it must: odd and even degrees contribute comparably as $L$ grows. **A
diffraction texture measurement constrains about half of the ODF and is silent on the rest.**

## Why It Is Called A Ghost

Setting the odd part to zero — the even-only reconstruction — is *a* solution consistent with the
data. It is not *the* ODF, and the difference has a characteristic appearance that gave the problem
its name.

An ODF built from even terms alone inherits their symmetry: every genuine texture component is
accompanied by spurious density elsewhere in orientation space. These are the **ghosts**. They are
not noise, they are not artefacts of the solver, and they do not diminish with more data. A
reconstruction can also go negative, which is impossible for a density and is the clearest
signal that the even-only solution is not physical.

Two consequences for reading a reconstructed ODF:

- **A weak secondary component may be a ghost of a strong primary one.** Before interpreting it,
  check whether it sits where the even-order truncation would place one.
- **Texture index and entropy are affected.** Both are integrals of $f$ or $f\log f$ over
  orientation space, so ghosts and negative lobes contribute to them. They are comparable between
  reconstructions made the same way, and not absolute measures of texture strength.

## What Actually Recovers The Odd Part

Three things, only one of which is a computation:

1. **Positivity.** The true ODF satisfies $f(g) \ge 0$ everywhere, which the even part alone does
   not enforce. This is genuine information about the odd part: it must be whatever makes the total
   non-negative. Positivity-enforcing and maximum-entropy methods exploit exactly this, and are the
   classical route to ghost correction, and it is the one PyTex implements — see
   [The correction PyTex applies](#the-correction-pytex-applies). The regularized least squares
   that produces the even part does **not** enforce positivity; the correction is a separate,
   separately reported step applied to its output.
2. **Measuring orientations instead of pole figures.** EBSD determines each grain's full
   orientation directly rather than projecting it onto a pole figure, so **an ODF estimated from
   EBSD has no ghost problem at all**. The degeneracy belongs to the diffraction *geometry*, not to
   texture analysis. This is a strong practical reason to prefer orientation-resolved data when the
   odd part matters, and it is why the two routes to an ODF are not interchangeable.
3. **Anomalous scattering**, which breaks Friedel's law and is the only way a diffraction
   experiment itself can see the difference. It is rarely exploited for texture.

## What PyTex Does

`HarmonicODF` carries `even_degrees_only` as an explicit, recorded field, defaulting to `True` when
the inputs are antipodal pole figures. The default is the honest one: it declines to invent the
half of the ODF the data cannot support, rather than returning odd coefficients that are artefacts
of the regularizer. The flag is part of the object, so a reconstruction states on its face which
half of orientation space it is speaking about.

The forward operator folds opposite normals together whenever the pole figure declares itself
antipodal, which is Friedel's law written into the model rather than only into the prose. That
folding is what makes the degeneracy exact in the code: an odd-degree basis function produces no
predicted pole density at all, so nothing in the fit can depend on it. Without the folding the
operator would appear to determine part of the odd component — an artefact of the model, not a
measurement — and a correction could not honestly claim to leave the fit alone.

(the-correction-pytex-applies)=
## The Correction PyTex Applies

`pytex.texture.correct_ghosts` recovers an odd part from positivity, and
`HarmonicODF.invert_pole_figures(..., ghost_correction=...)` runs it as part of an inversion.

**The two sets.** Write the corrected density as $f = \tilde f + \hat f$, where $\tilde f$ is the
even part the data determined and $\hat f$ lies in the span of the symmetry-projected odd-degree
harmonics up to the same bandlimit. Two convex sets act on it:

- $C_d$, the densities whose even part equals $\tilde f$. It is affine, and it is exactly the set
  of distributions consistent with the measurement.
- $C_+$, the densities that are physically admissible: non-negative everywhere and, under the
  zero-range method, identically zero over a range the data declare empty.

The correction is the point of $C_d$ closest to admissibility, taken as the minimizer of

$$
\Phi(\hat f) = \tfrac{1}{2}\int_{SO(3)} \big[ v(f(g)) \big]^{2} \,\mathrm{d}g
+ \tfrac{\mu}{2} \lVert \hat f \rVert^{2},
$$

where $v$ is the inadmissible part of the density — $\min(f, 0)$ outside the zero range, and $f$
itself inside it. The functional is smooth and convex in the odd coefficients, so a quasi-Newton
minimizer reaches its solution in tens of iterations: the same solution the classical alternating
projection between $C_d$ and $C_+$ converges to, without its long tail.

**Why the second term is not optional.** Once a corrected density is admissible, every remaining
direction in the odd subspace is free: positivity alone does not determine the odd part, it only
bounds it. An unregularized minimizer stops at whichever admissible point it reaches first, which
is a larger odd part than the data force, reported as though the data forced it. The term
$\mu \lVert \hat f \rVert^{2} / 2$ selects the *smallest* odd part that achieves admissibility.
That is a choice, and it is the defensible one when the alternative is an arbitrary choice; $\mu$
is exposed as `odd_regularization`.

**What it costs.** `GhostCorrectionReport.describe()` states the size of the inference — the odd
coefficient norm against the even one — the negative-density fraction and minimum density before
and after, the change to the texture index and entropy, and the residual violation the correction
could not remove. It also reports the largest change in the pole densities the ODF predicts at the
measured directions, which must be at the level of the quadrature error: a correction that moved
the fit would have bought positivity with data agreement it is not entitled to spend.

**What it recovers.** On an orthorhombic single-component texture whose answer is known by
construction, the even-only inversion is negative over 9% of orientation space with its maximum
depressed to 3.74 m.r.d. against a true 4.06. Positivity correction restores the maximum to
4.24 m.r.d., leaves 0.1% negative, and *halves* the quadrature-weighted distance to the true
distribution, while changing the predicted pole densities by 1e-3 m.r.d. against measured values of
order 1 to 4. Those numbers are computed and pinned by `tests/unit/test_ghost_correction.py`.

## Assumptions And Limits

- The degeneracy assumes Friedel's law. It is exact for the kinematic, non-anomalous case that
  covers essentially all laboratory texture measurement.
- Setting odd coefficients to zero is a choice, not a derivation. It minimizes the norm among
  data-consistent solutions and has no claim to being the physical one.
- **The correction is an inference, not a measurement.** No pole-figure experiment can confirm or
  refute the odd part it supplies. A reported ODF must say whether it was corrected, which is why
  the correction returns its own report rather than silently replacing the distribution.
- Correction addresses the missing odd part, not truncation. A texture too sharp for the bandlimit
  rings, and that ringing is negative for reasons positivity cannot repair with odd terms of the
  same bandlimit; the residual violation the report carries is what exposes such a case.
- A crystal symmetry with no odd-degree invariant below the bandlimit has no ghost part to correct.
  For the rotation group 432 the first one is at degree 9, so a cubic material expanded to degree 6
  or 8 is already as complete as its symmetry allows, and the correction says so rather than
  pretending to have worked.

## References

### Normative

- H.-J. Bunge, *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths (1969).
  The generalized spherical harmonic expansion and the fundamental equation relating an ODF to its
  pole figures.

### Informative

- S. Matthies, *On the reproducibility of the orientation distribution function of texture samples
  from pole figures (ghost phenomena)*, Physica Status Solidi (b) **92** (1979) K135–K138.
  DOI: <https://doi.org/10.1002/pssb.2220920253>. The identification and naming of the problem.
- S. Matthies and G. W. Vinel, *On the reproduction of the orientation distribution function of
  texturized samples from reduced pole figures using the conception of a conditional ghost
  correction*, Physica Status Solidi (b) **112** (1982) K111–K114.
  DOI: <https://doi.org/10.1002/pssb.2221120254>.
- R. Hielscher and H. Schaeben, *A novel pole figure inversion method: specification of the MTEX
  algorithm*, Journal of Applied Crystallography **41** (2008) 1024–1037.
  DOI: <https://doi.org/10.1107/S0021889808030112>. A modern inversion that enforces
  non-negativity.
- M. Dahms and H.-J. Bunge, *The iterative series-expansion method for quantitative texture
  analysis. I. General outline*, Journal of Applied Crystallography **22** (1989) 439-447.
  DOI: <https://doi.org/10.1107/S0021889889005261>. The positivity-driven iterative correction
  whose convex-projection form PyTex implements.
- V. Randle and O. Engler, *Introduction to Texture Analysis*, CRC Press (2000).

## See Also

- {doc}`harmonic_odf_reconstruction` — the reconstruction whose `even_degrees_only` default this
  note explains.
- {doc}`discrete_odf_and_pole_figures` — the discrete route, which shares the degeneracy whenever
  its input is pole figures rather than orientations.
- {doc}`pole_figure_arithmetic_and_mrd` — the other family of silent pole-figure errors.
