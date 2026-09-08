# Rietveld Whole-Pattern Profile Refinement

**Surface:** `pytex.diffraction.rietveld.refine_rietveld`, `RietveldResult`,
`RefinedParameter`, `REFINABLE_PARAMETERS`, `DEFAULT_REFINEMENT_SET`, with
`pytex.diffraction.xrd_instrument.InstrumentBroadening`,
`xrd_measurement.MeasuredPowderPattern`, and workbench operation `xrd.rietveld`.

Rietveld refinement (Rietveld, 1969) is the foundational whole-pattern profile fitting
technique in powder diffraction crystallography. Rather than attempting to deconvolute
individual integrated peak intensities from overlapping diffraction reflections, the
Rietveld method models the entire experimental step-scan intensity profile point-by-point
using structural crystallographic models, instrumental optics parameters, and empirical
background functions.

Every observation point across the two-theta profile acts as a constraint in a non-linear
least-squares optimization. This page details the physical forward model, profile
broadening functions, bounded parameter refinement strategies, diagnostic goodness-of-fit
statistics, and failure modes.

```{figure} ../../figures/rietveld_refinement_algorithm.svg
:alt: Four-lane flow sheet. Lane 1 takes the raw profile with its background
  intact and enumerates reflections once over a padded window. Lane 2
  recomputes positions, intensities, profiles and background every evaluation.
  Lane 3 minimises the weighted residual by bounded trust-region least squares.
  Lane 4 reports the R factors, Durbin-Watson, and the residual curve.
:width: 100%

Algorithmic execution architecture of whole-pattern Rietveld profile refinement.
```

## 1. The forward physical profile model

At each measured diffraction angle $2\theta_i$, the total calculated intensity
$y_i^{\text{calc}}$ is evaluated as the superposition of crystallographic diffraction
profiles and an underlying background:

$$
y_i^{\text{calc}} = b(2\theta_i) + s \sum_{k} m_k\,\lvert F_k\rvert^{2}\,L(\theta_k)\,P_k(\alpha_k)\,\Omega\!\left(2\theta_i - 2\theta_k\right),
$$

where each term models a distinct physical or instrumental contribution:

| Term | Physical Mechanism | Mathematical Formulation |
| --- | --- | --- |
| $b(2\theta_i)$ | Instrument and sample background | Orthogonal Chebyshev polynomial of degree $N_{\text{bg}} \in [4, 8]$ |
| $s$ | Global phase scale factor | Refined scalar proportional to phase volume fraction |
| $m_k$ | Reflection multiplicity | Number of symmetrically equivalent planes $\{hkl\}$ from point-group symmetry |
| $\lvert F_k\rvert^{2}$ | Kinematic structure factor | Atomic basis summation with isotropic Debye–Waller thermal parameters $B_{\text{iso}}$ |
| $L(\theta_k)$ | Lorentz-polarization factor | $L(\theta) = \frac{1 + \cos^2(2\theta)\cos^2(2\theta_{\text{mon}})}{\sin^2\theta\cos\theta}$ for laboratory geometry |
| $P_k(\alpha_k)$ | Preferred orientation factor | March–Dollase cylindrical fiber texture model |
| $\Omega(\Delta 2\theta)$ | Normalized profile shape | Thompson–Cox–Hastings pseudo-Voigt profile function |
| $2\theta_k$ | Bragg peak centroid | Computed from unit-cell parameters with zero-point error $\Delta 2\theta_0$ |

### 1.1 Instrumental and sample broadening: the Caglioti formulation

The angular dependence of the Gaussian profile component $\Gamma_G$ follows the
Caglioti polynomial relation (Caglioti et al., 1958):

$$
\Gamma_G^2 = U\tan^{2}\theta + V\tan\theta + W,
$$

while the Lorentzian component $\Gamma_L$ accounts for crystallite size ($X$) and
microstrain ($Y$) broadening:

$$
\Gamma_L = \frac{X}{\cos\theta} + Y\tan\theta.
$$

The overall FWHM $\Gamma$ and pseudo-Voigt mixing parameter $\eta$ are evaluated using
the analytical Thompson–Cox–Hastings (TCH) approximations (Thompson et al., 1987).

> [!IMPORTANT]
> Unconstrained simultaneous refinement of $U, V, W, X,$ and $Y$ against a single scan
> without an independent standard leads to severe parameter correlation. Instrumental
> parameters must be calibrated against a certified reference material (e.g., NIST SRM 660
> $\mathrm{LaB}_6$) via `InstrumentBroadening`, holding instrumental components fixed
> while refining sample microstrain and size parameters.

### 1.2 Reflection window padding and derivative continuity

To ensure numerical stability during optimization, reflection families are enumerated
over an angular window padded beyond the experimental limits:
$2\theta_{\text{min}} - \Delta \le 2\theta \le 2\theta_{\text{max}} + \Delta$.
If reflections were strictly truncated at experimental boundaries, small lattice parameter
increments would cause boundary reflections to enter or exit the active set discontinuously,
producing non-differentiable step artifacts in the least-squares objective function.

## 2. Parameter parameterization and optimization strategy

The refinement minimizes the weighted sum of squared residuals:

$$
\min_{\mathbf{p}} \; \Phi(\mathbf{p}) = \sum_{i=1}^N w_i \left[y_i^{\text{obs}} - y_i^{\text{calc}}(\mathbf{p})\right]^2,
\qquad w_i = \frac{1}{\sigma_i^2} = \frac{1}{y_i^{\text{obs}}},
$$

using bounded trust-region least squares (`scipy.optimize.least_squares`).

### 2.1 Parameter bounds and correlation control

Unbounded least-squares refinement can drive strongly correlated parameters into
unphysical parameter domains. PyTex enforces physical parameter bounds (`_PARAMETER_BOUNDS`):
- **Zero-point shift:** $\Delta 2\theta_0 \in [-1.0^\circ, +1.0^\circ]$. Zero errors exceeding
  one degree indicate instrumental misalignment rather than a refinable parameter offset.
- **Lattice parameter dilation:** $s_{\text{cell}} \in [0.9, 1.1]$ ($\pm 10\%$). Dilation
  exceeding this range indicates phase misidentification.
- **Chebyshev background coefficients:** Bounded to prevent negative total background values.

### 2.2 Incremental refinement hierarchy

Because diffraction parameters exhibit strong mutual correlations (e.g., zero shift
versus unit cell; overall thermal parameter $B_{\text{iso}}$ versus scale factor),
parameters must be activated incrementally:
1. **Initial cycle (`DEFAULT_REFINEMENT_SET`):** Refine overall scale factor $s$ and
   Chebyshev background coefficients $b_k$.
2. **Metric cycle:** Activate zero-point shift $\Delta 2\theta_0$ and unit-cell parameters.
3. **Profile cycle:** Activate Gaussian width parameter $W$ or Lorentzian parameter $Y$.
4. **Structural cycle:** Activate atomic displacement parameters $B_{\text{iso}}$ and
   preferred orientation parameters.

### 2.3 Joint background refinement versus pre-subtraction

Background polynomial coefficients must be refined simultaneously with structural
parameters. Subtracting an independently estimated baseline prior to refinement distorts
experimental Poisson counting statistics and eliminates parameter covariance between
background and peak tails, causing the least-squares Hessian to underestimate parameter
estimated standard deviations (ESDs).

## 3. Goodness-of-fit statistics and residual diagnostics

| Figure of Merit | Mathematical Definition | Physical Interpretation |
| --- | --- | --- |
| **Profile Residual** ($R_p$) | $R_p = \frac{\sum |y_i^{\text{obs}} - y_i^{\text{calc}}|}{\sum y_i^{\text{obs}}}$ | Unweighted absolute fractional profile error. |
| **Weighted Profile** ($R_{wp}$) | $R_{wp} = \sqrt{\frac{\sum w_i (y_i^{\text{obs}} - y_i^{\text{calc}})^2}{\sum w_i (y_i^{\text{obs}})^2}}$ | Objective metric minimized during optimization. |
| **Expected Error** ($R_{\text{exp}}$) | $R_{\text{exp}} = \sqrt{\frac{N - P}{\sum w_i (y_i^{\text{obs}})^2}}$ | Statistical lower limit governed purely by Poisson counting noise. |
| **Goodness of Fit** ($\text{GoF}$) | $\text{GoF} = \chi = \frac{R_{wp}}{R_{\text{exp}}} = \sqrt{\chi_{\text{red}}^2}$ | Ratio of actual fit residual to counting statistics floor. Values near $1.0$ indicate model convergence to noise limit. |
| **Bragg Residual** ($R_B$) | $R_B = \frac{\sum |I_k^{\text{obs}} - I_k^{\text{calc}}|}{\sum I_k^{\text{obs}}}$ | Agreement of partitioned integrated reflection intensities. |
| **Durbin–Watson** ($d$) | $d = \frac{\sum_{i=2}^N (e_i - e_{i-1})^2}{\sum_{i=1}^N e_i^2}$ | Serial correlation of normalized residuals $e_i = \sqrt{w_i}(y_i^{\text{obs}} - y_i^{\text{calc}})$. Values near $2.0$ indicate uncorrelated residuals; $d < 1.5$ signifies systematic misfit. |

All residual factors are reported as decimal fractions (e.g., $R_{wp} = 0.085$ corresponds to $8.5\%$).

### Diagnostic interpretation of residual profiles

The spatial distribution of residuals across $2\theta$ reveals specific modeling errors:
- **Derivative-shaped (S-curve) residuals:** Peak position error; indicates unrefined
  zero shift, sample displacement, or inaccurate unit-cell dimensions.
- **Symmetric bell-shaped residuals:** Profile shape mismatch; indicates inappropriate
  Gaussian/Lorentzian mixing $\eta$ or unmodeled size/strain broadening.
- **Monotonic high-angle discrepancy:** Thermal parameter ($B_{\text{iso}}$) inaccuracy or
  uncorrected absorption.
- **Systematic intensity deviations on specific reflection families:** Unmodeled crystallographic
  texture (preferred orientation).

## 4. Preferred orientation modeling

PyTex implements the March–Dollase cylindrical fiber texture model (Dollase, 1986):

$$
P_k(\alpha_k) = \left(r^2\cos^2\alpha_k + \frac{1}{r}\sin^2\alpha_k\right)^{-3/2},
$$

where $\alpha_k$ is the acute angle between the reflection normal $\mathbf{h}_k$ and the
declared preferred orientation direction (e.g., $[001]$ for plate-like crystallites).
The refined March parameter $r$ characterizes texture strength:
- $r = 1.0$: Ideal random isotropic powder.
- $r < 1.0$: Plate-like habit orientation along the designated axis.
- $r > 1.0$: Needle-like (acicular) habit orientation.

> [!CAUTION]
> Refining the March parameter without declaring a crystallographically justified habit plane
> is physically ungrounded and can mask systematic absorption or atomic occupancy errors.

## 5. Constraints and failure modes

| Diagnostic Condition | Root Physical / Numerical Cause | Corrective Action |
| --- | --- | --- |
| $\text{GoF} \gg 2.0$ with low $R_p$ | Over-counting statistics or unmodeled secondary phases | Inspect residual profile for unindexed peaks or asymmetric profiles |
| $\text{GoF} < 1.0$ | Overestimated experimental uncertainties $\sigma_i$ or smoothed data | Verify raw counting data; avoid applying digital smoothing filters |
| Durbin–Watson $d \ll 1.5$ | Strong serial correlation in residuals | Refine peak shape parameters or improve background modeling |
| Refined parameter hits boundary | Parameter divergence due to high correlation or incorrect model | Review model validity; freeze correlated parameter |
| Negative thermal factor ($B_{\text{iso}} < 0$) | Absorption distortion or severe preferred orientation | Verify linear absorption correction and background baseline |

## Verification

- `tests/unit/test_rietveld.py`: Validates forward pattern synthesis, derivative calculation,
  R-factor evaluations, and convergence on synthetic multi-phase standards.
- Executable worked examples:
  - {doc}`../examples/generated/lattice-parameters`

## See also

- {doc}`../theory/powder_xrd_and_saed` — Theoretical formulation of kinematic powder diffraction.
- {doc}`../theory/preferred_orientation_in_powder_intensities` — March–Dollase texture corrections.
- {doc}`precise_lattice_parameter_determination` — Lattice metric refinement without structural intensity modeling.
- {doc}`phase_identification` — Prior identification of constituent crystalline phases.

## References

### Normative

- Rietveld, H. M. (1969). A profile refinement method for nuclear and magnetic structures.
  *Journal of Applied Crystallography* **2**, 65–71. <https://doi.org/10.1107/S0021889869006558>
- Thompson, P., Cox, D. E. & Hastings, J. B. (1987). Rietveld refinement of Debye-Scherrer
  synchrotron X-ray data from $\mathrm{Al}_2\mathrm{O}_3$. *Journal of Applied Crystallography*
  **20**, 79–83. <https://doi.org/10.1107/S0021889887087090>
- Dollase, W. A. (1986). Correction of intensities for preferred orientation in powder
  diffractometry: application of the March model. *Journal of Applied Crystallography*
  **19**, 267–272. <https://doi.org/10.1107/S0021889886089458>

### Informative

- Toby, B. H. (2006). $R$ factors in Rietveld analysis: how good is good enough?
  *Powder Diffraction* **21**, 67–70. <https://doi.org/10.1154/1.2179804>
- Caglioti, G., Paoletti, A. & Ricci, F. P. (1958). Choice of collimators for a crystal
  spectrometer for neutron diffraction. *Nuclear Instruments* **3**, 223–228.
  <https://doi.org/10.1016/0369-643X(58)90029-X>
- McCusker, L. B., Von Dreele, R. B., Cox, D. E., Louër, D. & Scardi, P. (1999). Rietveld
  refinement guidelines. *Journal of Applied Crystallography* **32**, 36–50.
  <https://doi.org/10.1107/S0021889898009856>
