# Powder XRD Phase Identification

**Surface:** `pytex.diffraction.xrd_phase_identification.identify_phase`,
`identify_phase_from_pattern`, `PhaseCandidateScore`, `PhaseIdentification`,
built upon `pytex.diffraction.xrd_peaks.detect_and_fit_peaks` and
`pytex.diffraction.xrd_indexing.index_peaks`. In the workbench: **XRD → Identify the phase**.

Phase identification from powder X-ray diffraction (PXRD) profiles solves the classical
crystallographic search-match problem: determining the constituent crystalline phases
present within an unknown polycrystalline sample by matching observed diffraction peaks
against reference crystal structures. In real polycrystalline specimens, preferred
orientation (texture), solid-solution lattice expansion or contraction, instrumental
peak broadening, and multiphase superposition obscure naive fingerprint matching.

PyTex implements a deterministic, multi-stage search-match algorithm featuring global
bipartite reflection matching, scalar lattice dilation optimization, and a four-parameter
bounded figure of merit. This page details the algorithmic pipeline, specifies operational
parameters, provides diagnostic failure analysis, and establishes verification baselines.

```{figure} ../../figures/phase_identification_algorithm.svg
:alt: Four-lane flow sheet. Lane 1 detects and fits the peaks of the measured
  scan once and takes the candidate structures. Lane 2 runs per candidate,
  refining one cell dilation, enumerating that candidate's reflections and
  assigning peaks to lines by the Hungarian algorithm. Lane 3 scores four
  bounded criteria. Lane 4 ranks by their weighted mean and qualifies the
  winner twice.
:width: 100%

Algorithmic execution pipeline for powder XRD phase identification.
```

## 1. Algorithmic pipeline

| Stage | Computational Module | Output Data Structure | Potential Numerical Anomalies |
| --- | --- | --- | --- |
| 1. Peak detection & fitting | `detect_and_fit_peaks` | Centroid $2\theta_p$, intensity $A_p$, FWHM, ESD | Inappropriate prominence threshold causing peak omissions or fitting background noise |
| 2. Lattice dilation search | Internal grid optimization | Metric scale factor $s \in [1-\delta, 1+\delta]$ | Convergence to search boundary $\pm\delta$ indicating significant lattice mismatch |
| 3. Reflection generation | `generate_powder_reflections` | Kinematic $(hkl)$ lines and relative intensities | Zero reflections within experimental $2\theta$ window |
| 4. Bipartite assignment | `index_peaks` (Hungarian method) | Global one-to-one pairs $(2\theta_p \leftrightarrow hkl)$ | Unindexed peaks exceeding matching tolerance $\varepsilon$ |
| 5. Multi-criteria scoring | `PhaseCandidateScore` | Four normalized metrics in $[0, 1]$ | Incomplete indexing or missing predicted lines |
| 6. Qualification & ranking | `PhaseIdentification` | Sorted candidates, `is_conclusive`, `is_decisive` | Ambiguous score margin between top structural candidates |

Stage 1 executes once on the experimental pattern and is cached. Stages 2–5 execute
independently for each candidate structure, guaranteeing that ranking comparisons are
evaluated against an identical experimental peak list.

## 2. Peak detection and profile fitting

The function `identify_phase_from_pattern` executes automated background subtraction,
local peak detection, and non-linear least-squares pseudo-Voigt profile fitting as
detailed in {doc}`precise_lattice_parameter_determination`. The resulting `PeakTable`
is preserved on the final `PhaseIdentification` report.

Retaining the fitted `PeakTable` ensures analytical traceability: if detection thresholds
are improperly configured, systematic errors affect all candidate structures equally:
- An excessively high prominence threshold omits low-intensity reflections, artificially
  degrading completeness scores for all candidates.
- An excessively low prominence threshold fits diffuse background fluctuations as diffraction
  peaks, depressing the explained intensity fraction across all candidates.

| Parameter | Default | Operational Selection Criteria |
| --- | --- | --- |
| `prominence_sigma` | 5.0 | Peak detection significance threshold above background noise standard deviation. |
| `expected_fwhm_deg` | Estimated | Nominal instrumental FWHM; automatically estimated from scan resolution if omitted. |
| `minimum_two_theta_deg` | None | Low-angle cutoff to exclude direct-beam shadows and air-scattering signals. |
| `max_peaks` | 128 | Maximum peak capacity for indexing and bipartite assignment. |

## 3. Lattice dilation optimization

Solid-solution alloying, thermal expansion, and macroscopic residual stress alter crystal
lattice constants relative to nominal database entries without altering space-group
symmetry. Prior to assignment, PyTex optimizes an isotropic lattice scale factor $s$ by
minimizing the robust truncated distance objective:

$$
\Phi(s) = \sum_{p=1}^P \min\!\left(\min_{j} \left| 2\theta_p - 2\theta_j(s) \right|,\; \varepsilon\right),
\qquad \sin\theta_j(s) = \frac{\lambda}{2 s d_j},
$$

evaluated over an equispaced 401-point grid across $s \in [1-\delta, 1+\delta]$. A grid
search is employed because $\Phi(s)$ is non-convex and piecewise linear, exhibiting local
minima at discrete peak-line coincidences. The threshold parameter $\varepsilon$ prevents
unassociated outliers from biasing the optimum scale factor.

> [!NOTE]
> Isotropic dilation preserves interplanar spacing ratios ($d_i / d_j = \text{const}$),
> ensuring that metric scaling cannot artificially convert an incorrect crystal structure
> into a match. Refined values of $s$ departing significantly from unity indicate either
> chemical substitution or an incompatible structural model.

| Parameter | Default | Interpretation |
| --- | --- | --- |
| `cell_scale_range` ($\delta$) | 0.02 | Permissible isotropic lattice strain search bracket ($\pm 2\%$). Set to 0 to enforce rigid database lattice parameters. |

## 4. Reflection generation and global assignment

### 4.1 Symmetry-constrained reflection generation

Kinematic reflections are generated using candidate point-group symmetry, lattice
centering, and Wyckoff atomic coordinates. Systematic absences (screw axes, glide
planes, non-primitive Bravais centerings) are strictly enforced. Reflections with
theoretical relative intensities below `minimum_relative_intensity` ($0.001$ of peak line)
are excluded from candidate line matching.

### 4.2 Hungarian bipartite matching

Assigning observed peak centroids to theoretical reflection lines requires global
optimization. A greedy nearest-neighbor assignment can yield degenerate assignments,
mapping multiple observed peaks to a single reflection and stranding valid neighbors.
PyTex constructs a cost matrix $C_{pj} = |2\theta_p - 2\theta_j|$ and applies the
**Hungarian algorithm** (Kuhn, 1955) to compute a minimum-cost, strictly one-to-one
bipartite matching.

If a candidate structure predicts zero reflections within the experimental angular range,
it is scored zero and recorded as a rejected candidate with an explicit diagnostic reason
rather than raising an exception, preserving evaluation continuity across user-supplied libraries.

## 5. Composite figure of merit (FOM)

PyTex evaluates candidate structures against four bounded, complementary criteria
defined in $[0, 1]$:

| Metric | Mathematical Definition | Physical Significance |
| --- | --- | --- |
| **Explained Intensity Fraction** ($E$) | $E = \frac{\sum_{p \in \text{indexed}} A_p}{\sum_{p \in \text{all}} A_p}$ | Fraction of total experimental diffraction intensity accounted for by the candidate. |
| **Completeness** ($C$) | $C = \frac{N_{\text{observed}}}{N_{\text{predicted}}}$ | Proportion of theoretically expected reflection lines observed in the pattern. |
| **Position Agreement** ($P$) | $P = \max\!\left(0,\, 1 - \frac{\langle|\Delta 2\theta|\rangle}{\varepsilon}\right)$ | Average angular proximity between observed and calculated diffraction lines. |
| **Intensity Correlation** ($S$) | $S = 1 - \frac{1}{2}\sum_{i} |\hat{o}_i - \hat{c}_i|$ | Normalized $L_1$ correlation between observed and theoretical relative intensity vectors. |

The composite score $\text{FOM}$ represents a normalized weighted sum:

$$
\text{FOM} = w_E E + w_C C + w_P P + w_S S, \qquad \sum w_i = 1.
$$

If a metric is mathematically undefined (e.g., $S$ with fewer than two indexed reflections),
its weight is redistributed proportionally among the remaining active criteria.

### Weighting presets

| Preset Name | $(w_E, w_C, w_P, w_S)$ | Recommended Application Scenario |
| --- | --- | --- |
| **Balanced** | $(0.40, 0.25, 0.20, 0.15)$ | Untextured, well-ground isotropic powder specimens. |
| **Textured Specimen** | $(0.40, 0.30, 0.30, 0.00)$ | Rolled sheet, extruded rods, thin films exhibiting strong preferred orientation. |
| **Positions Only** | $(0.50, 0.00, 0.50, 0.00)$ | Preliminary screening relying solely on peak positions and intensity capture. |

## 6. Diagnostic interpretation of scoring profiles

Decomposing the match into four independent metrics allows unambiguous identification
of structural and microstructural defects:

| Diagnostic Pattern | Physical / Crystallographic Cause |
| --- | --- |
| High $P$, Low $C$ | Incorrect Bravais centering or space-group symmetry (predicting reflections absent in the true structure). |
| High $C$, High $P$, Low $E$ | Multiphase mixture: candidate is present, but unindexed peaks belong to secondary phases. |
| Low $P$, Moderate $E, C$ | Inaccurate lattice parameters or instrumental zero-shift aberration exceeding tolerance $\varepsilon$. |
| High $E, C, P$, Low $S$ | Correct phase exhibiting strong crystallographic preferred orientation (texture). |
| $s$ pinned at $1 \pm \delta$ | Lattice parameter mismatch exceeds permitted search bounds. |

## 7. Statistical qualification boundaries

The `PhaseIdentification` report validates matching decisions via two independent
Boolean criteria:

- **`is_conclusive`:** Evaluates whether the top-ranked candidate explains the pattern
  sufficiently to be credible:
  $$
  \text{FOM}_{\text{top}} \ge \text{minimum\_score} \quad (\text{default: } 0.55).
  $$
- **`is_decisive`:** Evaluates whether the top candidate is statistically distinct from
  competing candidates:
  $$
  \text{FOM}_{\text{top}} - \text{FOM}_{\text{runner-up}} \ge \text{decisive\_margin} \quad (\text{default: } 0.05).
  $$

When `is_decisive=False`, the experimental scan cannot distinguish between candidate
structures with statistical confidence, indicating the need for extended high-angle
counting, wavelength variation, or complementary chemical analysis.

## Verification

- `tests/unit/test_xrd_phase_identification.py`: Verifies ranking correctness on synthetic
  and experimental benchmarks, Bravais lattice discrimination, dilation recovery,
  and candidate order invariance.
- Executable worked examples:
  - {doc}`../examples/generated/phase-identification`

## See also

- {doc}`../theory/phase_identification_from_powder_patterns` — Detailed mathematical derivations of criteria metrics.
- {doc}`precise_lattice_parameter_determination` — Peak fitting and instrumental aberration correction.
- {doc}`rietveld_refinement` — Quantitative multi-phase whole-pattern profile refinement.

## References

### Normative

- Hanawalt, J. D., Rinn, H. W. & Frevel, L. K. (1938). Chemical analysis by X-ray diffraction.
  *Industrial & Engineering Chemistry Analytical Edition* **10**, 457–512. <https://doi.org/10.1021/ac50125a001>
- Smith, G. S. & Snyder, R. L. (1979). $F_N$: a criterion for rating powder diffraction patterns.
  *Journal of Applied Crystallography* **12**, 60–65. <https://doi.org/10.1107/S002188987901178X>
- de Wolff, P. M. (1968). A simplified criterion for the reliability of a powder pattern indexing.
  *Journal of Applied Crystallography* **1**, 108–113. <https://doi.org/10.1107/S002188986800508X>
- Dollase, W. A. (1986). Correction of intensities for preferred orientation in powder diffractometry.
  *Journal of Applied Crystallography* **19**, 267–272. <https://doi.org/10.1107/S0021889886089458>

### Informative

- Cullity, B. D. & Stock, S. R. (2001). *Elements of X-Ray Diffraction*, 3rd ed. Prentice Hall.
- Gates-Rector, S. & Blanton, T. (2019). The Powder Diffraction File. *Powder Diffraction* **34**,
  352–360. <https://doi.org/10.1017/S0885715619000812>
- Kuhn, H. W. (1955). The Hungarian method for the assignment problem. *Naval Research Logistics Quarterly*
  **2**, 83–97. <https://doi.org/10.1002/nav.3800020109>
