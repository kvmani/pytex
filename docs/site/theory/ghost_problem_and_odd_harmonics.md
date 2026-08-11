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
   classical route to ghost correction. PyTex's regularized least squares does **not** enforce
   positivity — hence the note in {doc}`harmonic_odf_reconstruction` that negative lobes can
   appear.
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

## Assumptions And Limits

- The degeneracy assumes Friedel's law. It is exact for the kinematic, non-anomalous case that
  covers essentially all laboratory texture measurement.
- Setting odd coefficients to zero is a choice, not a derivation. It minimizes the norm among
  data-consistent solutions and has no claim to being the physical one.
- No ghost correction is implemented. Reconstructions should be read as the even part of the ODF
  plus whatever the regularizer supplied, and compared only with reconstructions made identically.

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
- V. Randle and O. Engler, *Introduction to Texture Analysis*, CRC Press (2000).

## See Also

- {doc}`harmonic_odf_reconstruction` — the reconstruction whose `even_degrees_only` default this
  note explains.
- {doc}`discrete_odf_and_pole_figures` — the discrete route, which shares the degeneracy whenever
  its input is pole figures rather than orientations.
- {doc}`pole_figure_arithmetic_and_mrd` — the other family of silent pole-figure errors.
