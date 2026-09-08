# Precise Lattice Parameter Determination

**Surface:** `pytex.diffraction.xrd_peaks.detect_and_fit_peaks`,
`pytex.diffraction.xrd_indexing.index_peaks`,
`pytex.diffraction.xrd_lattice_parameter.determine_lattice_parameters`,
`determine_lattice_parameters_le_bail`,
`determine_lattice_parameters_from_pattern`, with
`pytex.diffraction.xrd_corrections` providing geometric aberration models and
extrapolation functions.

Determining unit-cell parameters with high precision (relative uncertainties below
$10^{-4}$) is fundamental to quantitative crystallography, solid-solution thermodynamics,
thermal expansion measurement, and residual stress analysis. While Bragg's law
($\lambda = 2d\sin\theta$) provides the basic relation between interplanar spacing
and diffraction angle, experimental peak positions are perturbed by instrumental
aberrations, specimen positioning errors, and counting statistics.

Differentiating Bragg's law demonstrates that angular sensitivity diverges toward back-reflection:

$$
\frac{\Delta d}{d} = -\cot\theta\,\Delta\theta \to 0 \quad \text{as } \theta \to 90^\circ.
$$

Consequently, high-angle reflections provide the highest intrinsic metrological precision,
yet they are also the most vulnerable to index misassignment from unrefined starting cells.
PyTex resolves this through two complementary mathematical frameworks:
1. **The Peak-Position Route (Cohen's Method):** Multi-pass iterative indexing coupled
   with generalized linear least-squares extrapolation across all seven crystal systems.
2. **The Whole-Pattern Route (Le Bail Method):** Structure-factor-free profile decomposition
   and non-linear least squares refinement across overlapping reflection envelopes.

## 1. Algorithmic architecture

| Stage | Computational Module | Output Data Structure | Physical Role / Diagnostic Objective |
| --- | --- | --- | --- |
| 1. Detection | `detect_peaks` | Candidate diffraction angles | Wavelet-based peak finding with Anscombe variance stabilization |
| 2. Profile Fitting | `fit_peaks` | Centroids, ESDs, FWHMs | Non-linear least-squares fitting of doublets and local backgrounds |
| 3. Indexing | `index_peaks` | Global $(hkl)$ assignments | Bipartite assignment matching observed peaks to calculated lines |
| 4. Iterative Refinement | `determine_lattice_parameters` | Cell metrics, ESDs, aberration scale | Multi-pass re-indexing eliminating metric divergence at high angles |
| 4'. Whole-Pattern Fit | `determine_lattice_parameters_le_bail` | Cell metrics, profile parameters | Iterative Le Bail intensity partitioning on overlapping reflections |

## 2. Stage 1: Wavelet peak detection with variance stabilization

Diffraction count profiles follow Poisson statistics, where variance equals expected
signal ($\sigma^2 = \mu$). Direct peak searching across decades of intensity variation
leads to over-detection in high-intensity regions and under-detection in low-intensity
tails.

PyTex applies the Anscombe variance-stabilizing transformation:

$$
y_i = 2\sqrt{\max(I_i - B_i, 0) + \tfrac{3}{8}},
$$

where $B_i$ represents the baseline background estimated via the Statistics-sensitive
Non-linear Iterative Peak-clipping (SNIP) algorithm. In the transformed domain, Poisson
noise approximates Gaussian white noise of unit variance ($\sigma \approx 1$).

Peak candidates are identified by continuous wavelet convolution using a normalized
Ricker (Mexican hat) wavelet kernel evaluated across scale parameters matching the
instrumental resolution function:

$$
R(2\theta) = \int y(\tau)\,k_{s(2\theta)}(2\theta - \tau)\,\mathrm{d}\tau.
$$

Local maxima exceeding the significance threshold $R > \text{prominence\_sigma} \cdot \sigma_{\text{noise}}$
(with $\sigma_{\text{noise}} = 1.4826 \cdot \operatorname{MAD}(R)$) are retained and refined
via sub-step parabolic interpolation. Candidate peaks corresponding to secondary $K\alpha_2$
doublets are identified and suppressed during initial peak harvesting.

| Parameter | Default | Physical Interpretation |
| --- | --- | --- |
| `prominence_sigma` | 5.0 | Peak detection cutoff in units of estimated noise standard deviation. |
| `expected_fwhm_deg` | 0.12° | Nominal instrumental FWHM; scaled across $2\theta$ by Caglioti relations. |
| `background_half_window_deg` | 2.0° | Clipping filter half-width for the SNIP baseline estimator. |
| `suppress_kalpha2` | `True` | Suppresses secondary doublet peaks during candidate generation. |

## 3. Stage 2: Constrained doublet profile fitting

Each detected candidate peak is isolated within an adaptive angular window and fitted
using bounded trust-region least squares. The analytical profile models the characteristic
$K\alpha_1 - K\alpha_2$ spectral doublet over a linear background:

$$
y(2\theta) = h\,P(2\theta; c, w, \eta) + h\,r\,P\!\left(2\theta; c_2(c), w, \eta\right) + b_0 + b_1(2\theta - c),
$$

where $P(2\theta)$ is a pseudo-Voigt profile function with Lorentzian mixing fraction $\eta$
and FWHM $w$. The secondary centroid $c_2$ is constrained by the physical wavelength ratio:

$$
\sin\theta_2 = \frac{\lambda_2}{\lambda_1}\sin\theta_1.
$$

The intensity ratio $r = I(K\alpha_2)/I(K\alpha_1)$ is fixed to the theoretical spectral
value (typically 0.5 for copper radiation). Propagated parameter standard deviations
(ESDs) are extracted from the covariance matrix:

$$
\sigma(c) = \sqrt{\left[\left(\mathbf{J}^{\mathsf{T}}\mathbf{W}\mathbf{J}\right)^{-1}\right]_{11}} \cdot \sqrt{\chi_{\text{red}}^2},
$$

where $\mathbf{W} = \operatorname{diag}(1/\sigma_i^2)$ and $\chi_{\text{red}}^2$ is the reduced chi-squared statistic.

## 4. Stage 3: Global bipartite indexing and iterative re-indexing

Observed peak centroids are assigned to theoretical reflection lines by solving a global
bipartite matching problem via the Hungarian algorithm on cost matrix $C_{ij} = |2\theta_{\text{obs}, i} - 2\theta_{\text{calc}, j}|$.

### The necessity of multi-pass re-indexing

An initial cell parameter differing from the true physical value by relative error
$e = \Delta a / a$ induces an angular displacement that diverges with Bragg angle:

$$
\Delta(2\theta) = 2e\tan\theta.
$$

For an ordinary mismatch of $e = 0.3\%$ at $2\theta = 120^\circ$ ($\theta = 60^\circ$),
the displacement is $\Delta(2\theta) \approx 0.6^\circ$, exceeding standard matching tolerances.
A single-pass indexer discards these crucial high-angle reflections.

PyTex implements an iterative re-indexing loop:
1. Low-angle reflections (where $\tan\theta$ is small) are indexed using the initial cell.
2. A preliminary Cohen least-squares fit refines the metric tensor, reducing the residual error $e$.
3. Theoretical line positions are recalculated using the updated cell, successfully capturing
   the high-angle reflections in subsequent passes.
4. Convergence is achieved when the set of indexed reflections stabilizes (typically 2 to 3 iterations).

## 5. Stage 4: Generalized Cohen least-squares refinement

Cohen's method (Cohen, 1935) establishes a unified linear least-squares framework that
simultaneously refines reciprocal metric tensor components and an instrumental drift parameter:

$$
\sin^2\theta_i = \frac{\lambda^2}{4}\mathbf{h}_i^{\mathsf{T}}\mathbf{G}^{*}\mathbf{h}_i + D\,\sin^2\theta_i\,f(\theta_i),
$$

where $\mathbf{G}^{*}$ is the reciprocal metric tensor, $D$ is the aberration coefficient,
and $f(\theta)$ is an analytical aberration extrapolation function.

### Linear system parameterization across crystal systems

By parameterizing $\mathbf{G}^{*}$ via system-specific symmetry constraints, the problem
reduces to the weighted linear least-squares system $\mathbf{X}\mathbf{p} \approx \mathbf{y}$:

| Crystal System | Free Metric Parameters ($n$) | Basis Components in Vector $\mathbf{p}$ |
| --- | --- | --- |
| Cubic | 1 | $a^{*2}$ |
| Tetragonal / Hexagonal | 2 | $a^{*2}, c^{*2}$ |
| Orthorhombic | 3 | $a^{*2}, b^{*2}, c^{*2}$ |
| Monoclinic | 4 | $a^{*2}, b^{*2}, c^{*2}, 2a^{*}c^{*}\cos\beta^{*}$ |
| Triclinic | 6 | $a^{*2}, b^{*2}, c^{*2}, 2b^{*}c^{*}\cos\alpha^{*}, 2a^{*}c^{*}\cos\beta^{*}, 2a^{*}b^{*}\cos\gamma^{*}$ |

The $(n+1)$-th column of $\mathbf{X}$ contains the aberration term $\sin^2\theta_i\,f(\theta_i)$.

### Aberration extrapolation functions

| Extrapolation Function $f(\theta)$ | Physical Aberration Mechanism | Instrumental Geometry |
| --- | --- | --- |
| `cos_squared_over_sin` | Specimen surface displacement ($\Delta h$) | Bragg–Brentano parafocusing |
| `cot_theta` | Detector zero-point angular offset ($\Delta 2\theta_0$) | Universal diffractometer geometry |
| `nelson_riley` | Combined specimen displacement and absorption | Classical standard compromise |
| `bradley_jay` | Specimen absorption in cylindrical cameras | Debye–Scherrer geometry |

Direct matrix inversion of the reciprocal metric tensor yields real-space lattice parameters:
$\mathbf{G} = (\mathbf{G}^{*})^{-1}$, with parameter ESDs evaluated via numerical Jacobian
transformation of the covariance matrix.

## 6. Stage 4': Whole-pattern Le Bail profile refinement

In materials exhibiting low crystal symmetry or severe peak overlap, individual peak
positions cannot be reliably extracted. The Le Bail method (Le Bail et al., 1988) refines
lattice parameters directly from the whole-pattern profile without requiring structural
model coordinates.

The observed background-subtracted profile is iteratively partitioned among candidate reflections:

$$
I_k = \sum_{i} y_{\text{net}, i}\,\frac{I_k\,P_k(2\theta_i)}{\sum_j I_j\,P_j(2\theta_i)},
$$

where $P_k(2\theta)$ is the normalized profile function of reflection $k$. The calculated
profile $y_{\text{calc}}(2\theta) = \sum_k I_k P_k(2\theta) + B(2\theta)$ is optimized
against observed counts using the Trust Region Reflective algorithm over:
- Reciprocal lattice parameters ($n \le 6$).
- Instrumental aberration parameter (specimen displacement or zero shift).
- Caglioti instrumental resolution parameters ($U, V, W$).
- Pseudo-Voigt Lorentzian mixing parameter ($\eta$).
- Residual linear background coefficients ($b_0, b_1$).

All least-squares weights reflect experimental Poisson uncertainties: $w_i = 1 / \sqrt{y_{\text{raw}, i}}$.

## 7. Method selection guide

| Experimental Scenario | Recommended Approach | Preferred Extrapolation / Setting |
| --- | --- | --- |
| Cubic materials with resolved peaks | **Cohen Method** | `f_theta="cos_squared_over_sin"` for Bragg–Brentano geometry |
| Hexagonal, tetragonal, orthorhombic | **Cohen Method** | Multi-pass re-indexing with Nelson–Riley extrapolation |
| Severe peak overlap or low symmetry | **Le Bail Method** | Whole-pattern refinement with refined displacement parameter |
| Calibrated external zero standard | **Cohen Method** | Calibrate zero offset first, then refine displacement only |

## Verification

- `tests/unit/test_xrd_peaks.py`: Validates doublet deconvolution, ESD calculation,
  and Anscombe transformation accuracy.
- `tests/unit/test_xrd_lattice_parameter.py`: Verifies multi-pass re-indexing convergence,
  Cohen metric parameter recovery across all 7 crystal systems, and Le Bail whole-pattern fitting.
- Executable worked examples:
  - {doc}`../examples/generated/lattice-parameters`

## See also

- {doc}`../theory/precise_lattice_parameter_determination` — Comprehensive derivations of aberration functions and metric covariance propagation.
- {doc}`phase_identification` — Upstream identification of candidate crystal structures.
- {doc}`rietveld_refinement` — Quantitative structural refinement incorporating atomic coordinates and site occupancies.

## References

### Normative

- Cohen, M. U. (1935). Calculation of precision lattice constants from X-ray powder photographs.
  *Review of Scientific Instruments* **6**, 68–74. <https://doi.org/10.1063/1.1751937>
- Nelson, J. B. & Riley, D. P. (1945). An experimental investigation of extrapolation methods
  in the derivation of accurate unit-cell dimensions of crystals. *Proceedings of the Physical
  Society* **57**, 160–177. <https://doi.org/10.1088/0959-5309/57/3/302>
- Le Bail, A., Duroy, H. & Fourquet, J. L. (1988). Ab-initio structure determination of
  $\mathrm{LiSbWO}_6$ by X-ray powder diffraction. *Materials Research Bulletin* **23**,
  447–452. <https://doi.org/10.1016/0025-5408(88)90019-0>
- Wilson, A. J. C. (1963). *Mathematical Theory of X-ray Powder Diffractometry*. Philips
  Technical Library.

### Informative

- Cullity, B. D. & Stock, S. R. (2001). *Elements of X-Ray Diffraction*, 3rd ed. Prentice Hall.
- de Wolff, P. M. (1968). A simplified criterion for the reliability of a powder pattern indexing.
  *Journal of Applied Crystallography* **1**, 108–113. <https://doi.org/10.1107/S002188986800508X>
- Smith, G. S. & Snyder, R. L. (1979). $F_N$: a criterion for rating powder diffraction patterns.
  *Journal of Applied Crystallography* **12**, 60–65. <https://doi.org/10.1107/S002188987901178X>
