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
needs the **orientation distribution function** $f(g)$ instead. Recovering $f(g)$ from a finite set of measured pole figures constitutes the
fundamental inverse problem of quantitative texture analysis. This inverse problem
is ill-posed due to projection non-injectivity, diffraction centrosymmetry
(Friedel's law), and experimental noise. This page details the discrete and
harmonic inversion formulations, their regularisation schemes, and their
physical limitations.

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

## 1. Mathematical foundation of the forward projection

The pole density of the plane family $\{hkl\}$ along the specimen direction
$\mathbf{y}$ is the ODF integrated over every orientation that puts some member
of the family along $\mathbf{y}$:

$$
P_{hkl}(\mathbf{y}) \;=\; \frac{1}{m}\sum_{i=1}^{m}\;
\int_{\{g\,:\,g\,\mathbf{h}_i \,\parallel\, \mathbf{y}\}} f(g)\,\mathrm{d}g .
$$

This is a **projection**: a one-parameter family of orientations — the rotations
about $\mathbf{y}$ — is integrated away at every point. Three physical consequences
govern the inversion:

1. **Non-injectivity.** Different ODFs can yield identical pole figures. A single
   pole figure never uniquely determines $f(g)$; multiple independent $\{hkl\}$
   projections constrain the solution space, but cannot restore uniqueness without
   regularizing assumptions.
2. **Friedel's law and the ghost problem.** In conventional X-ray and neutron
   diffraction, Friedel's law ($I_{\mathbf{h}} = I_{-\mathbf{h}}$) renders pole
   figures centrosymmetric. Consequently, measured pole figures contain information
   only about the **even-order** harmonic components of $f(g)$. The odd-order
   components are invisible to kinematic diffraction. See {doc}`ghost_correction`
   for non-negativity and entropy-based recovery methods.
3. **Discrete sampling and noise.** Experimental pole figures are sampled over
   finite angular grids and contain counting noise, limiting the achievable
   angular resolution.

PyTex provides two distinct mathematical paths: a discrete orientation dictionary
formulation and a continuous spherical harmonic series expansion.

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

Step 4 normalizes the kernel projection operator to multiples of a random
distribution (m.r.d.). Measured pole figure intensities are expressed in m.r.d.,
where an isotropic (random) polycrystal yields unit intensity everywhere. Because
the orientation weights $w$ are constrained to the probability simplex
($\sum_j w_j = 1$), the raw sum of spherical kernel densities (which exceeds
unity by the kernel spherical concentration factor) cannot be absorbed by weight
scaling. Normalizing each row block of $\mathbf{A}$ by the kernel's analytical
random pole density ensures metric consistency between the linear operator and
measured intensities, preventing artificial solver stalling.

### 2.2 The constrained least-squares problem

$$
\min_{w}\; \tfrac{1}{2}\lVert A w - b\rVert^{2}
        \;+\; \tfrac{1}{2}\lambda\lVert w\rVert^{2},
\qquad
w \ge 0,\quad \sum_j w_j = 1 .
$$

Both constraints reflect fundamental physical principles: an orientation distribution
function is a probability density, requiring non-negativity ($w_j \ge 0$) and unit
total probability ($\sum_j w_j = 1$). Solving an unconstrained least-squares system
and subsequently clipping negative values introduces severe distortion, as negative
density artifacts alter the amplitude of adjacent positive peaks.

The Tikhonov regularization parameter $\lambda$ balances residual minimization
against solution smoothness, preventing overfitting to experimental counting noise.

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

**Scaling of the stationarity metric.** The projected gradient step length is
$1/L$, meaning the unscaled parameter step $\lVert w_{n+1}-w_n\rVert$ scales
inversely with the Lipschitz constant $L$. For an operator matrix with large
entries, $L$ is large and the unscaled step $\lVert w_{n+1}-w_n\rVert$ is small
regardless of distance from the true stationary point. Evaluating unscaled
parameter increments against a fixed threshold can trigger premature convergence
termination at the initial uniform estimate. Multiplying the increment by $L$
reconstructs the projected gradient magnitude, while dividing by
$\lVert A^{\mathsf{T}} b\rVert$ yields a dimensionless stationarity metric that is
invariant to the absolute intensity scale of the input pole densities.

### 2.4 Diagnostic inversion reporting

The `ODFInversionReport` structure provides diagnostic metrics to evaluate
solution fidelity and numerical condition:

| Field | Interpretation |
| --- | --- |
| `relative_residual_norm` | Normalized residual $\lVert A w - b\rVert / \lVert b\rVert$. Values approaching 1 indicate that the model accounts for negligible measured variation. |
| `mean_absolute_error`, `max_absolute_error` | Residual errors in m.r.d., directly comparable to physical texture intensity peaks. |
| `converged`, `iterations`, `objective_history` | Optimization termination status and convergence trajectory across iterations. |
| `dictionary_coverage_ratio` | Ratio of measured pole directions to dictionary orientations ($M/N$). Values below unity signify an underdetermined system where regularization dominates the solution. |
| `predicted_intensities` | Forward-projected pole figure intensities for direct validation against experimental data. |

Because the spherical Radon projection is non-injective, a low residual norm on
the fitted pole figures is a necessary but insufficient condition for ODF
accuracy. Validating the reconstructed ODF against *unfitted* experimental pole
figures provides an independent verification of solution fidelity.

### 2.5 Computational complexity and resolution bounds

| Attribute | Scaling / Bound |
| --- | --- |
| Operator construction | $O(K \cdot n_{\text{pts}} \cdot N \cdot m)$, vectorized over dictionary orientations |
| Iteration cost | $O(N^{2})$ via precomputed Gram matrix $G = A^{\mathsf{T}} A$ |
| Angular resolution | Bounded jointly by dictionary spacing $\Delta g$ and kernel halfwidth $b$ |
| System conditioning | Severely underdetermined if only one pole figure or coplanar reflection poles are provided |

## 3. Route B — the harmonic route

`HarmonicODF.invert_pole_figures` implements the classical Bunge series
expansion. The unknown is a truncated set of symmetry-projected harmonic
coefficients $C_\ell^{\mu\nu}$ rather than weights on a support:

$$
f(g) = \sum_{\ell=0}^{L}\sum_{\mu,\nu} C_\ell^{\mu\nu}\, \dot{T}_\ell^{\mu\nu}(g).
$$

The projection operator for each generalized spherical harmonic basis function is
evaluated across measurement directions, forming a linear system solved via
regularized least squares.

Two physical and mathematical constraints govern the harmonic approach:

- **Bandwidth truncation at degree $L$ (`degree_bandlimit`).** The maximum series
  degree limits the minimum resolvable angular feature size ($\Delta \omega \approx 360^\circ / L$).
  Increasing $L$ improves angular sharpness but escalates parameter count ($O(L^3)$)
  and amplifies sensitivity to high-frequency experimental noise.
- **Null space of odd harmonic degrees.** In the presence of Friedel symmetry,
  diffraction pole figures provide constraints solely on even harmonic coefficients
  ($\ell = 0, 2, 4, \dots$). The odd coefficients ($\ell = 1, 3, 5, \dots$) reside in the
  null space of the projection operator. By default, `even_degrees_only=True` restricts
  the inversion to even orders ($f_{\text{even}}$), setting odd coefficients to zero.
  Activating `ghost_correction` reconstructs non-zero odd coefficients by enforcing
  physical non-negativity ($f(g) \ge 0$) via iterative entropy or positivity
  regularization; see {doc}`ghost_correction`.

### Method comparison: discrete versus harmonic routes

| Characteristic | Discrete Route (`ODF`) | Harmonic Route (`HarmonicODF`) |
| --- | --- | --- |
| Unknown parameters | Non-negative weights on orientation grid | Generalized harmonic coefficients $C_\ell^{\mu\nu}$ up to degree $L$ |
| Non-negativity | Enforced strictly via simplex projection ($w_j \ge 0$) | Not intrinsically constrained; truncation can yield negative lobes |
| Sharp textures | Limited by grid discretization and kernel width | Requires high $L$, prone to Gibbs ringing phenomena |
| Weak/smooth textures | Requires extensive orientation grid | Highly compact and computationally efficient |
| Ghost artifacts | Suppressed by strict non-negativity constraint | Explicit ghost correction required to resolve odd harmonics |
| Mathematical form | Discrete kernel mixture, ideal for Monte Carlo sampling | Continuous series, admitting exact analytical derivatives and integrals |

The PyTex workbench provides both options via the **Inversion route** selector on
`texture.measured_pole_figures`, allowing direct comparison between discrete and
harmonic reconstructions from identical experimental data.

## 4. Upstream corrections: experimental data preprocessing

The accuracy of ODF inversion depends strictly on experimental pole figure
calibration prior to optimization:

- **Defocus correction.** At elevated specimen tilt angles $\chi$, the beam
  footprint expands beyond the diffractometer receiving slit, causing geometric
  intensity attenuation unrelated to crystallographic texture.
  `defocus_from_random_standard` establishes the attenuation profile from a
  texture-free isotropic standard, and `PoleFigureCorrectionSpec` corrects measured
  intensities. Omitting defocus correction causes artificial intensity deficits at the
  periphery of recalculated pole figures.
- **Normalization to m.r.d.** As established in Section 2.1, the linear operator
  presupposes intensities expressed in multiples of a random distribution (m.r.d.).
  Scattered or unbinned count distributions must be regularized onto a standard
  equispaced spherical grid via `PoleFigure.on_grid` prior to operator assembly.

The function `residual_reports_for_pole_figures` computes point-by-point
discrepancies between measured and recalculated pole figures, identifying
systematic experimental artifacts (such as uncorrected defocus gradients) that
deviate from random Gaussian counting errors.

## 5. Downstream applications in PyTex

| Target Workflow | Application of Reconstructed ODF |
| --- | --- |
| `texture.odf_sections` | Visualization of constant-$\varphi_2$ sections with unified contour scaling |
| `fit_odf_components`, `component_volume_fractions` | Quantitative decomposition into ideal rolling, recrystallization, or shear texture components |
| `pytex.texture.kearns` | Calculation of Kearns orientation parameters $f$ by integrating basal plane densities |
| `pytex.texture.fibres` | Extraction of continuous orientation density along crystallographic fiber axes |
| `ODF.pole_figure` | Forward recalculation of pole figures, including unmeasured reflections for model validation |

## 6. Constraints and diagnostic conditions

| Diagnostic Scenario | Underlying Physical / Numerical Cause | Recommended Action |
| --- | --- | --- |
| Single pole figure | Severe non-uniqueness; solution dominated by Tikhonov prior $\lambda$ | Measure at least three independent crystallographic reflections $\{hkl\}$ |
| Nearly parallel reflection normals | Insufficient angular diversity among projection directions | Select reflection families with divergent reciprocal space angles |
| Discretization smearing | Dictionary grid spacing exceeds the physical orientation spread | Increase dictionary angular sampling density |
| Excess smoothing with low residual | Kernel halfwidth $b$ excessively large, filtering genuine texture peaks | Reduce kernel halfwidth until residual stabilizes near experimental noise |
| `converged=False` | Maximum iterations reached prior to satisfying stationarity tolerance | Increase `max_iterations` or verify numerical conditioning of operator |
| Reference frame mismatch | Specimen coordinate axes of input pole figures are mutually inconsistent | Verify coordinate alignments; PyTex raises a construction-time error |
| Missing odd harmonic components | Friedel's law restricts kinematic diffraction to even harmonic orders | Apply entropy- or positivity-based ghost correction; see {doc}`ghost_correction` |

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
