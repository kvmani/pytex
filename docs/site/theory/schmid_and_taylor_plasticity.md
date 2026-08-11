# Schmid Factors And The Taylor Factor

Two numbers describe how hard a crystal is to deform: the Schmid factor, which says how much of an
applied stress reaches a slip system, and the Taylor factor, which says how much slip an imposed
strain costs. The first has a closed-form maximum with a one-line proof. The second is a
constrained optimisation whose classical solution is combinatorial and whose implementation here is
not — and whose answer is unique even though the slip pattern producing it usually is not.

This note covers `pytex.properties.slip` and `pytex.properties.taylor`.

## The Schmid Factor, And Why It Cannot Exceed One Half

A slip system is a plane normal $\hat{\mathbf{n}}$ and a slip direction $\hat{\mathbf{d}}$ lying in
that plane, so $\hat{\mathbf{n}} \cdot \hat{\mathbf{d}} = 0$. Under uniaxial stress $\sigma$ along
$\hat{\mathbf{t}}$ the resolved shear stress on the system is

$$
\tau = \sigma \, m,
\qquad
m = (\hat{\mathbf{t}} \cdot \hat{\mathbf{n}})(\hat{\mathbf{t}} \cdot \hat{\mathbf{d}})
  = \cos\phi \, \cos\lambda ,
$$ (eq-sc-schmid)

with $\phi$ the angle to the plane normal and $\lambda$ the angle to the slip direction. Slip begins
when $\tau$ reaches a critical value, so $m$ ranks systems by how favourably they are oriented.

**The maximum is exactly $1/2$, and the bound is two lines.** Write $a = \hat{\mathbf{t}}\cdot\hat{\mathbf{n}}$
and $b = \hat{\mathbf{t}}\cdot\hat{\mathbf{d}}$. Since $\hat{\mathbf{n}}$ and $\hat{\mathbf{d}}$ are
orthonormal, they are two members of an orthonormal basis, so Bessel's inequality gives
$a^{2} + b^{2} \le \lVert\hat{\mathbf{t}}\rVert^{2} = 1$. Then

$$
m = ab \le \frac{a^{2}+b^{2}}{2} \le \frac{1}{2},
$$ (eq-sc-half)

with equality when $|a| = |b| = 1/\sqrt{2}$ — that is, $\phi = \lambda = 45^{\circ}$ and
$\hat{\mathbf{t}}$ in the plane spanned by $\hat{\mathbf{n}}$ and $\hat{\mathbf{d}}$. Sampling
$2\times10^{5}$ random orientations of fcc octahedral slip reaches $0.4999994$, approaching the
bound as it must.

For the cube orientation under $[001]$ tension, all twelve $\{111\}\langle 110 \rangle$ systems take
$|m| \in \{0,\ 1/\sqrt{6}\}$: eight are equally stressed at $1/\sqrt{6} = 0.408248$ and four carry
nothing. That eightfold degeneracy is the reason a cube-oriented grain has no single preferred slip
system, and it is the first hint of the ambiguity that dominates the Taylor problem.

## The Taylor Problem: Five Constraints, Not Six

Under the Taylor full-constraint hypothesis every grain suffers the *same* strain as the aggregate.
The grain must therefore accommodate an imposed deviatoric strain $\boldsymbol{\varepsilon}$ using
its available slip systems, and among the combinations that can, the operative one is that of least
total slip:

$$
\Gamma = \min \sum_{s} \gamma^{(s)}
\quad\text{subject to}\quad
\sum_{s} \gamma^{(s)} \mathbf{N}^{(s)} = \boldsymbol{\varepsilon},
\quad \gamma^{(s)} \ge 0 ,
$$ (eq-sc-taylor)

where $\mathbf{N}^{(s)} = \tfrac{1}{2}\bigl(\hat{\mathbf{d}}\hat{\mathbf{n}}^{\mathsf{T}} +
\hat{\mathbf{n}}\hat{\mathbf{d}}^{\mathsf{T}}\bigr)$ is the symmetric Schmid tensor in the *sample*
frame. The Taylor factor is that minimum normalised by the von Mises equivalent strain,

$$
M = \frac{\Gamma}{\varepsilon_{\mathrm{eq}}},
\qquad
\varepsilon_{\mathrm{eq}} = \sqrt{\tfrac{2}{3}\,\boldsymbol{\varepsilon} : \boldsymbol{\varepsilon}} .
$$ (eq-sc-taylor-factor)

Two structural points make the implementation what it is.

**The tensor equation is five scalar constraints, not six.** Both sides are symmetric and traceless
— slip is volume-conserving, and the imposed strain is deviatoric by construction — so matching
five independent components forces the sixth. PyTex constrains
$(11), (22), (12), (13), (23)$ and lets $(33)$ follow. Constraining all six would make the system
rank-deficient and the solver's behaviour dependent on how it handles the redundancy.

**Slip is signed, but the variables are not.** A slip system can shear in either sense, yet
{eq}`eq-sc-taylor` requires $\gamma^{(s)} \ge 0$ so that "total slip" is a sum rather than a
cancellation. The implementation appends the negated Schmid tensors as extra columns, turning 12
systems into 24 non-negative variables covering both senses. This is exact, not an approximation:
any signed solution maps to a non-negative one of the same total magnitude.

**It is a linear program.** The classical route enumerates the $\binom{12}{5} = 792$ ways of
choosing five systems from twelve, solves each $5\times5$ system, discards those with negative or
inconsistent slips, and keeps the cheapest — the Bishop–Hill approach. {eq}`eq-sc-taylor` is
exactly the same problem written as an LP, and simplex reaches the same vertex without visiting all
792. PyTex vectorises every part of the setup — sample-frame rotation, Schmid tensors, constraint
matrices — across the whole orientation population and leaves only the LP solve per orientation,
because that part genuinely is per orientation.

## Closed Forms And Known Answers

| Case | Value | Computed |
| --- | ---: | ---: |
| Max Schmid factor, any orientation | $1/2$ | $0.4999994$ (2×10⁵ samples) |
| Schmid factor, fcc $[001]$ tension | $1/\sqrt{6} = 0.408248$ | $0.4082483$ |
| Taylor factor, fcc $[001]$ tension | $\sqrt{6} = 2.449490$ | $2.449490$ |
| Taylor factor, fcc $[111]$ tension | $3\sqrt{6}/2 = 3.674235$ | $3.674235$ |
| Taylor factor, random fcc texture | $\approx 3.06$ (Taylor 1938) | $3.055 \pm 0.009$ |

The two single-orientation values are exact and reproduced to six figures, which tests the LP
formulation, the five-component constraint, and the signed-slip doubling at once. The random-texture
average is the classic result Taylor obtained in 1938 for a randomly oriented fcc aggregate in
tension; $3.055 \pm 0.009$ over 2000 orientations covers $3.06$ within one standard error.

The spread matters as much as the mean. Over a random texture $M$ runs from about $2.29$ to
$3.67$, so **the hardest orientation is roughly 60% harder than the softest**, and $[111]$ is the
hard end while $[001]$ is near the soft end. A textured sheet can therefore differ substantially
in flow stress from a random one at identical composition and grain size, which is the practical
reason for computing $M$ at all.

## The Answer Is Unique; The Slip Pattern Is Not

The minimum $\Gamma$ in {eq}`eq-sc-taylor` is unique, but **the set of slips achieving it usually is
not**. The cube orientation is the clearest case: eight systems share the same Schmid factor, so
many different five-system combinations reach the same total slip. This is the classical Taylor
ambiguity, and it is a property of the model rather than of the solver.

Two consequences follow. The Taylor *factor* is well defined and safe to report. The predicted
**lattice rotation is not**, because different optimal slip combinations rotate the lattice
differently, and texture-evolution predictions built on an arbitrary LP vertex inherit that
arbitrariness. Anyone extending this surface to texture evolution must add an explicit
tie-breaking rule — rate sensitivity, or a maximum-work criterion — rather than taking whichever
vertex the simplex happened to land on.

## Assumptions And Limits

- **Full constraint is an upper bound.** Requiring every grain to match the aggregate strain
  over-constrains real deformation, so Taylor over-predicts flow stress. Relaxed-constraint and
  self-consistent models are less stiff; none is implemented here.
- **Rigid-plastic, rate-independent, no hardening.** All systems share one critical resolved shear
  stress and it does not evolve, so $M$ describes the onset of plasticity and not a flow curve.
- **`inf` is a real answer.** When the supplied family cannot span the imposed strain the LP is
  infeasible and the factor is `inf`. That is correct reporting, not failure: a family with fewer
  than five independent systems genuinely cannot accommodate an arbitrary deviatoric strain, which
  is von Mises's criterion.
- Only $\{111\}\langle 110 \rangle$ fcc and $\{110\}\langle 111 \rangle$ bcc families are
  provided. Pencil glide, non-Schmid effects, and twinning are outside the model.
- Grain interaction, grain shape and neighbourhood are absent by construction.

## References

### Normative

- *International Tables for Crystallography, Volume D: Physical Properties of Crystals*, IUCr.
  DOI: <https://doi.org/10.1107/97809553602060000105>.

### Informative

- G. I. Taylor, *Plastic strain in metals*, Journal of the Institute of Metals **62** (1938)
  307–324. The full-constraint hypothesis and the $M \approx 3.06$ result for random fcc.
- J. F. W. Bishop and R. Hill, *A theory of the plastic distortion of a polycrystalline aggregate
  under combined stresses*, Philosophical Magazine **42** (1951) 414–427.
  DOI: <https://doi.org/10.1080/14786445108561065>. The dual formulation and the enumeration
  {eq}`eq-sc-taylor` replaces.
- U. F. Kocks, C. N. Tomé and H.-R. Wenk, *Texture and Anisotropy*, Cambridge University Press
  (1998). Taylor ambiguity, relaxed-constraint models, and texture evolution.
- R. von Mises, *Mechanik der plastischen Formänderung von Kristallen*, ZAMM **8** (1928) 161–185.
  DOI: <https://doi.org/10.1002/zamm.19280080302>. The five-independent-systems criterion.

## See Also

- {doc}`elastic_anisotropy_and_homogenization` — the elastic counterpart, where the analogous
  bounds are Voigt and Reuss.
- {doc}`discrete_odf_and_pole_figures` — where a weighted orientation population comes from.
