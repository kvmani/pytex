# Inverting Pole Figures Into An ODF

**Surface:** `pytex.texture.models.ODF.invert_pole_figures` (discrete route),
`pytex.texture.harmonics.HarmonicODF.invert_pole_figures` (series route),
`ODFInversionReport`, `HarmonicODFReconstructionReport`,
`pytex.texture.ghosts.correct_ghosts`, with
`pytex.texture.reconstruction` supplying defocus correction and residual
reporting, and the workbench operation `texture.measured_pole_figures`.

A diffractometer measures **pole figures**: the density of one plane normal
$\{hkl\}$ over specimen directions. Every physical model of a polycrystal —
elastic and plastic anisotropy, the Kearns factor, variant selection —
needs the **orientation distribution function** $f(g)$ instead. Recovering
$f$ from a handful of pole figures is the central inverse problem of
quantitative texture analysis, and it is ill-posed in three distinct ways that
this page keeps separate, because they have different cures and only two of them
have any cure at all.

```{figure} ../../figures/pole_figure_inversion_algorithm.svg
:alt: Four-lane flow sheet. Lane 1 takes measured pole figures on the m.r.d.
  scale and applies the defocus correction. Lane 2 builds the kernel response
  operator and rescales it to the observations' scale. Lane 3 solves the
  non-negative, simplex-constrained, regularised least-squares problem by
  projected gradient, with the harmonic series as the alternative unknown.
  Lane 4 reports residuals and recalculates unfitted poles.
:width: 100%

The algorithm, with the constraint governing each stage.
```

## 1. The forward model, which is where the difficulty comes from

The pole density of the plane family $\{hkl\}$ along the specimen direction
$\mathbf{y}$ is the ODF integrated over every orientation that puts some member
of the family along $\mathbf{y}$:

$$
P_{hkl}(\mathbf{y}) \;=\; \frac{1}{m}\sum_{i=1}^{m}\;
\int_{\{g\,:\,g\,\mathbf{h}_i \,\parallel\, \mathbf{y}\}} f(g)\,\mathrm{d}g .
$$

This is a **projection**: a one-parameter family of orientations — the rotations
about $\mathbf{y}$ — is integrated away at every point. Three consequences
follow, and they are the whole subject:

1. **The map is not injective.** Different ODFs give identical pole figures. One
   pole figure never determines $f$; several independent $\{hkl\}$ do better but
   never reach uniqueness.
2. **Friedel's law halves the information.** A diffraction experiment cannot
   distinguish $\mathbf{h}$ from $-\mathbf{h}$, so pole figures are centro
   symmetric and determine only the **even-order** part of $f$. The odd part is
   invisible to the measurement. This is the *ghost problem*, and it is
   structural: no amount of data or regularisation recovers the odd part,
   because no odd information was recorded. See
   {doc}`ghost_correction` for what can be inferred instead.
3. **The data are finite and noisy**, so even the even part is recovered only up
   to a resolution the measurement supports.

PyTex offers two routes, differing in what the unknown *is*.

## 2. Route A — the discrete route

The unknown is a set of non-negative weights on a fixed dictionary of
orientations. `ODF.invert_pole_figures` implements it.

### 2.1 Building the operator

```text
input : pole figures P_1..P_K (intensities in m.r.d.), dictionary G of N orientations,
        kernel psi with halfwidth b

1  for each pole figure k and each measured direction y_i in it:
2      for each dictionary orientation g_j:
3          A[i, j] = sum over the {hkl} family of psi(angle between g_j h and y_i)
4  A <- A / random_pole_density(psi)        -- put the operator on the m.r.d. scale
5  stack the per-figure blocks into one A; stack the intensities into one b
```

Step 4 is not cosmetic and is worth stating because getting it wrong produces a
failure that *reports success*. The observations are pole densities in multiples
of a random distribution, so a random texture reads 1.0 everywhere. The raw
kernel sum for a random texture is not 1 but the kernel's spherical mean — at a
$12^\circ$ halfwidth, a factor of order 64. The weights are constrained to sum
to one (section 2.2), so the model **cannot** absorb that factor into its
amplitude: the solver stalls at a relative residual near 1 and its stationarity
test then sees a step that has stopped moving and declares convergence. Dividing
the operator by the kernel's random level puts both sides on the same scale, and
the system becomes fittable rather than merely mis-scaled.

### 2.2 The constrained least-squares problem

$$
\min_{w}\; \tfrac{1}{2}\lVert A w - b\rVert^{2}
        \;+\; \tfrac{1}{2}\lambda\lVert w\rVert^{2},
\qquad
w \ge 0,\quad \sum_j w_j = 1 .
$$

Both constraints are physics, not numerical convenience. An ODF is a probability
density: it cannot be negative, and it integrates to one. Solving unconstrained
and clipping afterwards gives a different — and worse — answer, because the
negative excursions have already distorted the fitted positive lobes.

The Tikhonov term $\lambda$ trades detail for stability. It is the knob that
decides how much of the noise the solution is allowed to explain.

### 2.3 The solver

Projected gradient on the simplex:

```text
w <- uniform (1/N each)
G <- A^T A ;  r <- A^T b ;  L <- ||G + lambda I||_2
repeat up to max_iterations:
    grad      <- G w - r + lambda w
    candidate <- project_onto_simplex(w - grad / L)
    stationarity <- ||candidate - w|| * L / ||A^T b||
    w <- candidate
    stop when stationarity <= tolerance
```

**Why stationarity is scaled the way it is.** The step length is $1/L$, so the
raw step $\lVert w_{n+1}-w_n\rVert$ is proportional to $1/L$. On a system whose
operator has large entries the very first step is tiny *for that reason alone*,
and testing the raw step against a fixed tolerance declares convergence
immediately and returns the uniform starting guess as the answer — a smooth,
plausible, entirely uninformative ODF. Multiplying by $L$ recovers the projected
-gradient magnitude and dividing by $\lVert A^{\mathsf{T}} b\rVert$ makes it
dimensionless, so one tolerance means the same thing whatever units the pole
densities carry.

### 2.4 What the report carries, and how to read it

`ODFInversionReport` exists so the fit can be judged rather than trusted:

| Field | Read it for |
| --- | --- |
| `relative_residual_norm` | the headline. Near 1 means the model explained nothing |
| `mean_absolute_error`, `max_absolute_error` | in m.r.d., so directly interpretable against the texture strength |
| `converged`, `iterations`, `objective_history` | whether the solver stopped or merely ran out |
| `dictionary_coverage_ratio` | observations per dictionary orientation. Below 1 the system is underdetermined and the regularisation is doing the deciding |
| `predicted_intensities` | for a measured-versus-recalculated pole figure, which is the only honest visual check |

A converged fit with a low residual is **not** evidence of a correct ODF: the
map is non-injective, so a wrong $f$ can reproduce the data exactly. The
recalculated pole figures of *poles that were not fitted* are the check that
carries information.

### 2.5 Cost and limits

| | |
| --- | --- |
| Operator build | $O(K \cdot n_{\text{pts}} \cdot N \cdot m)$, vectorised over the dictionary |
| Per iteration | one Gram product, $O(N^{2})$ |
| Angular detail | bounded by the dictionary resolution **and** by the kernel halfwidth, whichever is coarser |
| Failure mode | one pole figure, or several from nearly parallel poles, leaves the system underdetermined; the answer is then mostly the regularisation |

## 3. Route B — the harmonic route

`HarmonicODF.invert_pole_figures` implements the classical Bunge series
expansion. The unknown is a truncated set of symmetry-projected harmonic
coefficients $C_\ell^{\mu\nu}$ rather than weights on a support:

$$
f(g) = \sum_{\ell=0}^{L}\sum_{\mu,\nu} C_\ell^{\mu\nu}\, \dot{T}_\ell^{\mu\nu}(g).
$$

The response of each basis function at each measured point is built the same
way, and the system is solved as regularised least squares.

Two limits are intrinsic and are *not* worked around:

- **Truncation at `degree_bandlimit`** bounds the recoverable detail. Raising it
  raises cost and noise sensitivity together.
- **Odd degrees are not determined.** `even_degrees_only` defaults to the honest
  choice, because pole figures carry no odd information. Passing
  `ghost_correction` recovers an odd part from positivity and reports what that
  inference cost; without it the odd part is silently zero, which is a specific
  and named error rather than a neutral default.

### Choosing between the routes

| | discrete | harmonic |
| --- | --- | --- |
| Unknown | weights on a dictionary | coefficients to degree $L$ |
| Non-negativity | enforced exactly | not enforced; positivity is a post-check |
| Sharp textures | dictionary resolution limits it | needs high $L$, which rings |
| Smooth textures | needs a large dictionary | compact and natural |
| Ghost correction | not applicable | available |
| Natural output | a discrete ODF, ready for sampling | a series, ready for analytic integration |

The workbench exposes both as the **Inversion route** control on
`texture.measured_pole_figures`, so the same measured data can be inverted both
ways and the answers compared — which is the practical test of whether a feature
of the ODF is real or an artefact of one method.

## 4. Before inverting: what the measurement needs first

An inversion is only as good as the pole figures entering it, and two
corrections belong upstream of everything above.

- **Defocus.** At high tilt the irradiated area leaves the focusing circle and
  the measured intensity falls for geometric reasons that have nothing to do
  with texture. `defocus_from_random_standard` calibrates the fall-off from a
  texture-free standard, and `PoleFigureCorrectionSpec` applies it. Uncorrected,
  the inversion faithfully reproduces an instrumental artefact as a rim of low
  density.
- **Normalisation to m.r.d.** The operator is built on the multiples-of-random
  scale (section 2.1). A scattered pole *cloud* whose intensities are per-pole
  weights is not on that scale; resample it onto a grid with
  `PoleFigure.on_grid` first.

`residual_reports_for_pole_figures` compares measured and recalculated figures
per pole, which is where a bad defocus correction shows up as a systematic
residual in the outer rings rather than as noise.

## 5. How the rest of PyTex uses the result

| Consumer | Uses the ODF for |
| --- | --- |
| `texture.odf_sections` | constant-$\varphi_2$ sections, on one shared contour ladder |
| `fit_odf_components`, `component_volume_fractions` | volume fractions of named components (Cube, Goss, Brass …) |
| `pytex.texture.kearns` | the Kearns factor $f$, by integrating basal-pole density against $\sin\phi\cos^{2}\phi$ |
| `pytex.texture.fibres` | fibre density along a declared axis |
| `ODF.pole_figure` | recalculated pole figures, including of poles never measured — the check of section 2.4 |

## 6. Constraints and failure modes

| Situation | What happens | What to do |
| --- | --- | --- |
| One pole figure | badly underdetermined; result dominated by $\lambda$ | measure at least three independent $\{hkl\}$ |
| Poles nearly parallel | same, less obviously | choose poles spanning the sector |
| Dictionary too coarse | sharp components smear to the dictionary spacing | refine the dictionary, not the kernel |
| Kernel too wide | everything smooth, residual acceptable | reduce halfwidth; watch the residual rise as noise is admitted |
| `converged=False` | ran out of iterations | raise `max_iterations` before trusting anything |
| Mismatched specimen frames | refused at call time | this is a construction-time invariant, not a warning |
| Odd part needed | zero unless ghost-corrected | see {doc}`ghost_correction` |

## Verification

- Pole-figure arithmetic and the m.r.d. scale, in
  {doc}`../examples/generated/pole-figure-arithmetic`.
- A uniform ODF gives pole density 1.000 m.r.d.; a known component mixture is
  recovered exactly by `fit_odf_components`, in
  {doc}`../examples/generated/texture`.
- Ghost correction and what it recovers, in
  {doc}`../examples/generated/ghost-problem`.

## See also

- {doc}`../theory/discrete_odf_and_pole_figures` — the discrete representation
  and its integrals.
- {doc}`../theory/harmonic_odf_reconstruction` — the series expansion in full.
- {doc}`../theory/ghost_problem_and_odd_harmonics` — why the odd part is absent.
- {doc}`../theory/pole_figure_arithmetic_and_mrd` — the m.r.d. convention.
- {doc}`ghost_correction` — the algorithm that infers an odd part.
- {doc}`ipf_coloring` — the other direction: reading an orientation as a colour.

## References

### Normative

- Bunge, H. J. (1982). *Texture Analysis in Materials Science: Mathematical
  Methods*. Butterworths. <https://doi.org/10.1016/C2013-0-11769-2>
- Matthies, S., Vinel, G. W. & Helming, K. (1987). *Standard Distributions in
  Texture Analysis*. Akademie-Verlag.

### Informative

- Hielscher, R. & Schaeben, H. (2008). A novel pole figure inversion method:
  specification of the MTEX algorithm. *Journal of Applied Crystallography*
  **41**, 1024-1037. <https://doi.org/10.1107/S0021889808030112>
- Randle, V. & Engler, O. (2000). *Introduction to Texture Analysis*. CRC Press.
  <https://doi.org/10.1201/9781482287479>
