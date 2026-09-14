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

## 8. Reading the workbench report

The workbench operation **XRD → Determine lattice parameters** returns the answer three ways —
the summary sentence, the residual table and the plot — and, under **How this result was
reached**, one section for every stage of the computation described above. A final cell is only
as believable as the steps behind it, and each of those steps can fail in its own way: a missed
peak, a wrong assignment, a systematic error absorbed into the cell. The stages put that chain of
evidence on the page.

Each stage carries a status mark. **✓** means the stage did its job. **!** means it completed but
left something that must be checked before the stages after it are believed; the first such stage
opens itself. **i** marks a stage that records an input rather than a computation. The same
stages are written into the Markdown report, and into the Excel workbook as a `Stages` sheet of
every stage quantity plus one sheet per stage table. In the library they come from
`lattice_parameter_pipeline`, which returns the peak table and the pass history together with the
`LatticeParameterResult` and the final `PeakIndexing`;
`determine_lattice_parameters_from_pattern` is the same computation returning only the last two.

### Stage 1 — The scan as read

| Quantity | Meaning and what to check |
| --- | --- |
| Points, first 2θ, last 2θ | The profile that was analysed. A range that stops short of back-reflection forfeits the reflections with the smallest $\cot\theta$, which carry most of the precision. |
| Median step | Sampling interval in degrees. A centre can be located to a small fraction of a step only with roughly five or more points across a peak's FWHM. |
| Radiation, wavelength $\lambda$ | Every spacing is computed from $\lambda$; the wrong radiation scales every cell edge by the ratio of the two wavelengths. |
| Injected specimen displacement | Demonstration scans only: the known aberration added so the methods can be judged against it. |

### Stage 2 — Peak detection and profile fitting

The headline quantities are the number of peaks detected, how many fits converged, the detection
threshold (in robust noise standard deviations, `prominence_sigma`), the expected width, and the
median $\sigma(2\theta)$ and FWHM. The table has one row per peak:

| Column | Meaning |
| --- | --- |
| 2θ | Fitted Kα₁ centre, in degrees. |
| σ(2θ) | Standard uncertainty of that centre, in millidegrees, from the fit covariance of Stage 2 above. It becomes the weight of the reflection in the least squares, so an imprecise peak counts for little. |
| Height, Integrated intensity | Fitted peak height and area above the local linear background. |
| FWHM | Fitted full width at half maximum, in degrees. |
| η | Lorentzian fraction of the pseudo-Voigt: 0 is Gaussian, 1 is Lorentzian. |
| χ²ν | Reduced chi-squared of that peak's profile fit. Near 1, the profile describes the peak to within counting noise; far above 1 flags an overlapped, asymmetric or badly backgrounded peak whose position deserves suspicion. |
| Converged | Whether the optimizer converged. A non-converged fit stays in the table but should not be trusted. |

The stage is marked **!** when any fit failed to converge, or when there are no more peaks than
refined parameters, since the least squares then has no redundancy to detect an error with.

### Stage 3 — Index-then-determine passes

One row per pass of the re-indexing loop of Stage 3 above.

| Column | Meaning |
| --- | --- |
| Pass | Pass number; pass 1 indexes against the tabulated cell of the selected phase. |
| a | The cell edge determined from that pass's assignment, in ångströms. For a pass that was not taken, the cell it was indexed against. |
| Indexed, Unindexed | Peaks assigned and left over in that pass. |
| M | de Wolff's figure of merit of that pass's assignment (defined in Stage 4 below). |
| Mean \|Δ2θ\| | Mean absolute position discrepancy of the assignment, in millidegrees. |
| Outcome | *taken*, or *not taken: indexed no more* — the pass that ends the loop. |

The metric *Reflections recovered by re-indexing* is the final indexed count minus the first. A
positive value is the direct evidence that the starting cell misplaced high-angle reflections by
$\Delta(2\theta) = 2e\tan\theta$ and that re-indexing recovered them.

### Stage 4 — Reflection assignment

The final one-to-one assignment the cell is fitted to.

| Column | Meaning |
| --- | --- |
| Reflection | Miller indices of the calculated line (Miller–Bravais for hexagonal phases). |
| 2θ observed, σ(2θ) | The fitted peak centre and its uncertainty. |
| 2θ calculated | Where the cell of the final pass puts the reflection, **before** any systematic correction. |
| Δ2θ | Observed minus calculated, in millidegrees. A smooth trend with angle is expected when a zero or displacement error is present, and Stage 5 removes it; a single reflection far off the trend is a misassignment. |
| d observed, d calculated | Interplanar spacings from $\lambda/(2\sin\theta)$ and from the cell. |
| Multiplicity | Number of symmetry-equivalent planes contributing to the line. |
| I calculated | Calculated relative intensity. It ranks lines for matching and is never used to fit the cell. |

The quantities beside the table are the indexed fraction, de Wolff's
$M_N = Q_N / (2\langle|\Delta Q|\rangle N_{\text{poss}})$ with $Q = 1/d^2$, Smith and Snyder's
$F_N = N / (\langle|\Delta 2\theta|\rangle N_{\text{poss}})$, the number of unindexed peaks, and the
number of calculated lines above the intensity threshold that were not observed. $N$ is 20 for
$M$ and 30 for $F$ unless fewer lines were indexed, and the subscript always shows the $N$ used,
because $M_7$ and $M_{20}$ are not comparable. $M_N > 10$ is a plausible cell and $M_N > 20$ a
convincing one. The stage is marked **!** when $M_N < 10$ or any peak stayed unindexed; an
unindexed strong peak belongs to something the phase does not describe — a second phase, a Kβ or
tungsten line, or the sample holder. When every Δ2θ carries the same sign, the summary says so:
that is the signature of an uncorrected zero or displacement error, not of a wrong cell.

### Stage 5 — Least-squares determination

For Cohen's method, the quantities are:

| Quantity | Meaning |
| --- | --- |
| Method, Extrapolation function | The chosen method and the function $f(\theta)$ of Stage 4 above. |
| Refined parameters | The free reciprocal-metric components the crystal system allows (for example `a*^2`, `c*^2`) and, when a systematic term is refined, `D`. |
| Observations, Degrees of freedom | $N$ reflections and $N - p$ for $p$ refined parameters. |
| Drift coefficient D, σ(D), \|D\|/σ(D) | The refined systematic-error coefficient, its standard uncertainty and their ratio; above about 2 the term is significantly different from zero. |
| Largest systematic shift | The largest angular correction $D$ applied to any reflection, in millidegrees. Compare it with σ(2θ): a correction much larger than the position uncertainties did real work. |
| Reduced χ² | $\chi^2_\nu = \sum_i (r_i/\sigma_i)^2/(N-p)$ in $\sin^2\theta$. About 1: the residuals match the position uncertainties. Much larger: an unmodelled error or a misassignment. Much smaller: the uncertainties are overstated. |
| Angular floor, Reflections discarded by the floor | The `minimum_two_theta_deg` restriction and how many assigned reflections it removed after the passes converged. |

The table is the correlation matrix of the refined parameters,
$r_{jk} = V_{jk}/\sqrt{V_{jj}V_{kk}}$ with
$\mathbf{V} = (\mathbf{X}^{\mathsf{T}}\mathbf{W}\mathbf{X})^{-1}\chi^2_\nu$. Its diagonal is 1 and
every entry lies in $[-1, 1]$. A cell parameter correlated with $D$ beyond about $\pm 0.95$ means
the scan's angular range barely separates a change of cell from the systematic error, which is why
that parameter's uncertainty is larger than the scatter alone suggests. The stage is marked **!**
when $\chi^2_\nu > 3$ or any off-diagonal $|r| > 0.98$.

For the average method there is no joint fit: the stage reports the number of reflections and
the reduced χ² of the per-reflection values, and states that no systematic term could be refined.

### Stage 6 — The determined cell

The cell edges and their standard uncertainties that the crystal system leaves free — $a$ for
cubic; $a$, $c$ and $c/a$ for hexagonal and tetragonal; also $b$ for orthorhombic and lower; and
the angles for monoclinic and triclinic — then the relative uncertainty $\sigma(a)/a$, the
reference $a$ of the selected phase, and the lattice strain $(a - a_{\text{ref}})/a_{\text{ref}}$.
The summary grades the relative uncertainty: below $5\times10^{-5}$ is strain-grade, below
$5\times10^{-4}$ composition-grade, and anything larger identification-grade. The uncertainty is
the precision of this determination on this scan, not the accuracy of the instrument; calibrate
against a certified standard such as NIST SRM 640 (silicon) or SRM 660 (lanthanum hexaboride)
before quoting an absolute value. The strain is along the scattering vector of a symmetric scan,
normal to the surface, and is not a stress.

### Stage 7 — Cross-check against the other methods

The same assignment of Stage 4, and the same angular floor, pushed through the alternative
methods: Cohen with no systematic term, Cohen with the Nelson–Riley term (when that is not already
the reported method), and, for a cubic cell, the average over reflections.

| Column | Meaning |
| --- | --- |
| Method | The first row is the reported determination; the others are alternatives. |
| a, σ(a) | Cell edge and its standard uncertainty from that method. |
| Difference from reported a | $10^6\,(a_{\text{method}} - a_{\text{reported}})/a_{\text{reported}}$, in parts per million. |
| χ²ν | Reduced chi-squared of that method's fit. |

Because the peaks and the assignment are shared, every difference is due to the method alone. A
large gap between the fits with and without a systematic term says the term did real work; a gap
within a few $\sigma(a)$ says the specimen was well aligned. Prefer the method whose $\chi^2_\nu$ is
closest to 1.

### Le Bail runs

A whole-pattern decomposition forms no peak list, so its report has three stages: the scan; the
**whole-pattern decomposition**, reporting the reflections modelled, the reduced χ², $R_{wp}$ and
the refined systematic term (specimen displacement in millimetres or detector zero in degrees 2θ);
and the determined cell. $R_{wp}$ is computed on the background-subtracted profile, so it is
systematically higher than a Rietveld program's $R_{wp}$ on the raw scan and must not be compared
with one. The stage is marked **!** when $\chi^2_\nu > 3$; its diagnostic is the difference curve
in the plot, not a per-reflection residual.

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
