# Identifying A Phase Among Candidate Structures

**Surface:** `pytex.diffraction.xrd_phase_identification.identify_phase`,
`identify_phase_from_pattern`, `PhaseCandidateScore`, `PhaseIdentification`, on top of
`pytex.diffraction.xrd_peaks.detect_and_fit_peaks` and
`pytex.diffraction.xrd_indexing.index_peaks`. In the workbench: **XRD → Identify the phase**.

The derivations are in {doc}`../theory/phase_identification_from_powder_patterns`. This page is
the implementation: the stages in order, the settings that matter with the values they take, what
each stage produces on a case with a known answer, and how it fails.

```{figure} ../../figures/phase_identification_algorithm.svg
:alt: Four-lane flow sheet. Lane 1 detects and fits the peaks of the measured
  scan once and takes the candidate structures. Lane 2 runs per candidate,
  refining one cell dilation, enumerating that candidate's reflections and
  assigning peaks to lines by the Hungarian algorithm. Lane 3 scores four
  bounded criteria. Lane 4 ranks by their weighted mean and qualifies the
  winner twice.
:width: 100%

The ranking, with the constraints that keep it honest.
```

## 1. The pipeline, and where each stage can go wrong

| Stage | Surface | Produces | Fails as |
| --- | --- | --- | --- |
| 1. Detect and fit | `detect_and_fit_peaks` | positions, ESDs, integrated intensities | too few peaks (threshold high) or background structure fitted as peaks (threshold low) — penalises every candidate equally and wrongly |
| 2. Refine dilation | internal, per candidate | one scalar $s$ | $s$ pinned at $\pm\delta$: the candidate was stretched as far as allowed and still did not fit |
| 3. Enumerate | `generate_powder_reflections` | that candidate's own $(hkl)$ lines | no line in range at all — recorded as a rejection, never raised |
| 4. Assign | `index_peaks` | $(hkl)$ per peak, $M_N$, $F_N$ | nothing assigned: the cell is wrong by more than the tolerance can absorb |
| 5. Score | `PhaseCandidateScore` | four criteria in $[0,1]$ | see §6 for what each failure means |
| 6. Rank and qualify | `PhaseIdentification` | ranking, `is_conclusive`, `is_decisive` | a winner that is neither — read §7 before believing it |

Stages 2–5 run once per candidate; stage 1 runs once and is shared, so the peaks every candidate is
judged against are literally the same peaks.

## 2. Stage 1 — the peaks, and why they are returned

`identify_phase_from_pattern` runs detection and fitting exactly as
{doc}`precise_lattice_parameter_determination` §2–3 describes, then filters to converged fits, then
ranks. The fitted `PeakTable` is **returned alongside the identification, not consumed inside it**.

That is deliberate. Peak detection is where an identification most often goes wrong, and the
failure is invisible in the ranking: a threshold that swallows the weak lines lowers every
candidate's completeness together, and a threshold that promotes shoulders of the background into
peaks lowers every candidate's explained intensity together. Neither changes the *order*, so the
ranking looks healthy while resting on a peak list nobody looked at. The workbench draws the fitted
positions on the scan for the same reason.

| Setting | Default | Choose it by |
| --- | --- | --- |
| `prominence_sigma` | 5.0 | lower to 3–4 to admit weak lines and expect background structure with them |
| `expected_fwhm_deg` | estimated | only needs to be right within about a factor of two |
| `minimum_two_theta_deg` | none | set it above a beam-stop shadow or air-scatter rise, which are not diffraction |
| `max_peaks` | 128 | rarely binding on a laboratory scan |

## 3. Stage 2 — the cell dilation

For each candidate a single scalar $s$ is refined before anything is indexed, minimising

$$
\Phi(s) = \sum_{p} \min\!\left(\min_{j} \left| 2\theta_p - 2\theta_j(s) \right|,\; \varepsilon\right),
\qquad \sin\theta_j(s) = \frac{\lambda}{2 s d_j},
$$

by a grid of 401 points on $[1-\delta, 1+\delta]$. A grid rather than a gradient method because
$\Phi$ is piecewise linear with a local minimum at every near-coincidence. The clip at $\varepsilon$
keeps a line that is nowhere near any peak from dominating the sum.

Reflections are enumerated once over the full angular range rather than over the measured window,
because a line outside the window at $s = 1$ may be inside it at another $s$. Values of $s$ for
which $\lambda / (2 s d_j) > 1$ are excluded: those reflections are past back-reflection.

**Why this cannot rescue a wrong candidate** is derived in
{doc}`../theory/phase_identification_from_powder_patterns` §4 and checked as an executable worked
example: a uniform dilation preserves every ratio of $d$ spacings exactly, and the ratios are what
indexing tests.

| Setting | Default | Choose it by |
| --- | --- | --- |
| `cell_scale_range` ($\delta$) | 0.02 | covers alloying, thermal expansion and residual stress; set 0 to match the CIF cells exactly |

Reported per candidate as `cell_scale`, and in the workbench as the **Cell** column in per cent.
Read it: a few hundredths of a per cent is an ordinary composition difference; a value at the edge
of the range is a candidate that did not fit.

## 4. Stages 3–4 — enumeration and assignment

Each candidate's reflections come from its *own* symmetry and systematic absences, not from a
generic $(hkl)$ list, which is what makes the completeness criterion meaningful. Families below
`minimum_relative_intensity` (default 0.001 of the strongest) are never offered for matching:
assigning a strong observed peak to a line that should be invisible is not an explanation.

The assignment is the Hungarian algorithm on $|2\theta_{\text{obs}} - 2\theta_{\text{calc}}|$, as in
{doc}`precise_lattice_parameter_determination` §5 — global and one-to-one, because a greedy
nearest-line pass can assign two peaks to one reflection and strand the true partner, and neither
failure is visible in the result it produces.

**A candidate that predicts no line in range is recorded, not raised.** It is scored zero with the
reason carried in `rejection`, and the ranking of the others proceeds. This matters specifically
because the candidates are user-supplied: one unreadable or implausible CIF among five must cost
the user that one, not the whole comparison.

## 5. Stage 5 — the four criteria

Defined and derived in {doc}`../theory/phase_identification_from_powder_patterns` §2. In summary:

| Criterion | Expression | Undefined when |
| --- | --- | --- |
| `explained_intensity_fraction` | $\sum_{\text{indexed}} A_p \big/ \sum_{\text{all}} A_p$ | — |
| `completeness` | observed strong lines / predicted strong lines, inside the measured span | — |
| `position_score` | $\max(0,\, 1 - \langle|\Delta 2\theta|\rangle / \varepsilon)$ | nothing indexed |
| `intensity_agreement` | $1 - \tfrac12 \sum_i |\hat o_i - \hat c_i|$ | fewer than two indexed lines |

An undefined criterion is renormalised out of the weighted mean rather than scored zero.

The weighting is a parameter. The workbench offers three presets, stated as specimen situations
rather than as four numbers, because the situation is what the operator knows:

| Preset | $(E, C, P, S)$ | For |
| --- | --- | --- |
| Balanced | $(0.40, 0.25, 0.20, 0.15)$ | a well-prepared random powder |
| Textured specimen | $(0.40, 0.30, 0.30, 0.00)$ | rolled sheet, coatings, anything with a rolling or fibre texture |
| Positions only | $(0.50, 0.00, 0.50, 0.00)$ | the strictest reading: line positions and unexplained intensity alone |

## 6. Reading a failure

The point of reporting four criteria rather than one score is that the pattern of failure names the
fault. Worked through on a synthetic nickel scan offered four candidates:

| Candidate | Score | $E$ | $C$ | $P$ | $s - 1$ | Reading |
| --- | --- | --- | --- | --- | --- | --- |
| Nickel (fcc) | 0.962 | 1.00 | 1.00 | 0.94 | $+0.26\%$ | correct; the dilation recovers the scan's deliberate one |
| Copper (fcc) | 0.617 | 0.78 | 0.40 | 0.34 | $-2.00\%$ | right centring, wrong cell size — pinned at the search-range edge |
| Ferrite (bcc) | 0.582 | 0.49 | 0.25 | 0.94 | $+0.64\%$ | good positions, quarter of its lines seen: wrong centring |
| Halite (NaCl) | 0.407 | 0.15 | 0.20 | 0.94 | $-0.94\%$ | explains a seventh of the intensity: not present |

Copper and ferrite score similarly and fail completely differently — copper on position, ferrite on
completeness — which the total alone would hide. The general readings:

| Symptom | Reading |
| --- | --- |
| High $P$, low $C$ | right cell metric, wrong centring or basis |
| High $C$ and $P$, low $E$ | the candidate is present and so is something else: a second phase |
| Low $P$, everything else moderate | wrong cell dimensions, or $\varepsilon$ narrower than the instrument's aberrations |
| High $E, C, P$, low $S$ | right framework, wrong basis — *or* a textured specimen; check that first |
| $s$ pinned at $\pm\delta$ | stretched as far as permitted and still not fitting |

## 7. Stage 6 — the two qualifications

`is_conclusive` and `is_decisive` are separate booleans because they fail for different reasons and
have different remedies:

| | Test | If false |
| --- | --- | --- |
| `is_conclusive` | best score $\ge$ `minimum_score` (0.55) | none of the candidates offered accounts for this pattern: widen the list, suspect a mixture, or check $\varepsilon$ against the aberrations |
| `is_decisive` | best $-$ runner-up $\ge$ `decisive_margin` (0.05) | this *scan* does not tell the top two apart: count longer at high angle, change wavelength, or use chemistry |

`describe()` states whichever holds in prose, with its remedy. `margin` is `nan` for a single
candidate rather than zero, because a comparison of one has no margin and reporting zero would read
as a tie; `is_decisive` is then true by convention, since nothing was chosen between.

## 8. Failure modes, and what they look like

| Symptom | Cause | Fix |
| --- | --- | --- |
| Every candidate scores near zero | detection found no peaks, or the radiation is wrong | check the returned peak table first, then `radiation` |
| The true phase indexes only its low-angle lines | `cell_scale_range` set to 0 with a real specimen | leave the default; the effect grows as $\tan\theta$ |
| Two candidates tie at a high score | genuinely indistinguishable on this angular range | §7: the remedy is measurement |
| A wrong candidate scores above the threshold | the range admits many lines and the candidate predicts many | read $C$ and $E$, not the total; consider `positions_only` |
| Ranking changes when candidates are added | it cannot: scores are computed per candidate and independently | if observed, that is a defect — the property is pinned by a test |
| A strong peak is unindexed by the winner | a second phase | multi-phase {doc}`rietveld_refinement`, not a further search |

## 9. Verification

- `tests/unit/test_xrd_phase_identification.py` — 44 tests, including that the generating phase of
  a synthetic pattern ranks first; that a centring is separated by completeness rather than by
  position; that the refined dilation recovers an imposed one; that the ranking is independent of
  the order candidates were offered in; and that a candidate predicting nothing is recorded rather
  than raised.
- `tests/unit/test_app_xrd.py` — the service operation, including candidates supplied as uploaded
  CIF text.
- {doc}`../examples/index` — four executable worked examples with independent provenance,
  including the algebraic identity that makes the cell dilation safe.

## 10. References

### Normative

- Hanawalt, J. D., Rinn, H. W. & Frevel, L. K. (1938). Chemical analysis by X-ray diffraction.
  *Industrial & Engineering Chemistry Analytical Edition* **10**, 457-512.
  <https://doi.org/10.1021/ac50125a001>
- Smith, G. S. & Snyder, R. L. (1979). $F_N$: a criterion for rating powder diffraction patterns.
  *Journal of Applied Crystallography* **12**, 60-65.
  <https://doi.org/10.1107/S002188987901178X>
- de Wolff, P. M. (1968). A simplified criterion for the reliability of a powder pattern indexing.
  *Journal of Applied Crystallography* **1**, 108-113.
  <https://doi.org/10.1107/S002188986800508X>
- Dollase, W. A. (1986). Correction of intensities for preferred orientation in powder
  diffractometry. *Journal of Applied Crystallography* **19**, 267-272.
  <https://doi.org/10.1107/S0021889886089458>

### Informative

- Cullity, B. D. & Stock, S. R. (2001). *Elements of X-Ray Diffraction*, 3rd ed., Ch. 14.
  Prentice Hall.
- Gates-Rector, S. & Blanton, T. (2019). The Powder Diffraction File. *Powder Diffraction* **34**,
  352-360. <https://doi.org/10.1017/S0885715619000812>
- Kuhn, H. W. (1955). The Hungarian method for the assignment problem. *Naval Research Logistics
  Quarterly* **2**, 83-97. <https://doi.org/10.1002/nav.3800020109>

## See also

- {doc}`../theory/phase_identification_from_powder_patterns` — the derivations
- {doc}`precise_lattice_parameter_determination` — detection, fitting and indexing in detail
- {doc}`rietveld_refinement` — the quantitative step once the phases are settled
- {doc}`../workflows/xrd_generation` — the workbench workflow
- {doc}`../workflows/phases_and_cif` — how a candidate arrives from a CIF
