# Rietveld Refinement: Fitting The Whole Powder Profile

**Surface:** `pytex.diffraction.rietveld.refine_rietveld`, `RietveldResult`,
`RefinedParameter`, `REFINABLE_PARAMETERS`, `DEFAULT_REFINEMENT_SET`, with
`pytex.diffraction.xrd_instrument.InstrumentBroadening`,
`xrd_measurement.MeasuredPowderPattern`, and the workbench operation
`xrd.rietveld`.

Rietveld's insight was to stop extracting peak intensities and fit the
**measured profile point by point** instead. Overlapped peaks then need no
deconvolution — they are modelled as what they are, a sum — and every point of
the pattern constrains the model. This page states the forward model, the
parameters, the refinement strategy, and the statistics, with particular
attention to the ways a refinement produces an excellent $R_{wp}$ and a wrong
answer.

```{figure} ../../figures/rietveld_refinement_algorithm.svg
:alt: Four-lane flow sheet. Lane 1 takes the raw profile with its background
  intact and enumerates reflections once over a padded window. Lane 2
  recomputes positions, intensities, profiles and background every evaluation.
  Lane 3 minimises the weighted residual by bounded trust-region least squares.
  Lane 4 reports the R factors, Durbin-Watson, and the residual curve.
:width: 100%

The refinement, with the bounds and the practice that govern it.
```

## 1. The forward model

At each measured $2\theta_i$ the calculated intensity is

$$
y_i^{\text{calc}} \;=\;
  b(2\theta_i) \;+\; s\sum_{hkl}
     m_{hkl}\,\lvert F_{hkl}\rvert^{2}\,L(\theta)\,P_{hkl}\,
     \Omega\!\left(2\theta_i - 2\theta_{hkl}\right),
$$

each factor a separate modelling decision:

| Factor | What it is | Model used |
| --- | --- | --- |
| $b$ | background | Chebyshev polynomial, degree 4-8 |
| $s$ | scale | one refined scalar |
| $m_{hkl}$ | multiplicity | from the phase point group |
| $\lvert F_{hkl}\rvert^{2}$ | structure factor | from the atomic basis and $B_{\text{iso}}$ |
| $L(\theta)$ | Lorentz-polarisation | geometry of the diffractometer |
| $P_{hkl}$ | preferred orientation | March-Dollase |
| $\Omega$ | peak shape | Thompson-Cox-Hastings pseudo-Voigt |
| $2\theta_{hkl}$ | peak position | from the dilated cell plus the zero shift |

### 1.1 Peak widths: the Caglioti form

The Gaussian width follows

$$
\mathrm{FWHM}^{2} = U\tan^{2}\theta + V\tan\theta + W,
$$

with a Lorentzian component $Y$ carrying the size and strain broadening. $U$,
$V$, $W$ and $Y$ start from the supplied `InstrumentBroadening` — the
resolution function of the diffractometer, ideally calibrated on a standard —
so a refinement separates specimen broadening from instrumental only when the
instrument function is known independently. Refining all four against one
pattern with no standard does not measure microstrain; it redistributes it.

### 1.2 Reflection enumeration, and the discontinuity it avoids

Reflection families are enumerated **once**, from the starting cell, over an
angular window **padded** beyond the fitted range. The padding matters: as the
cell dilates during refinement a reflection near the edge would otherwise move
into or out of the fitted set, changing the model discontinuously and giving the
least-squares solver a step in its objective where it expects a derivative.

## 2. The refinement

```text
input : measured profile, phase, refine-list, instrument, background degree

1  enumerate reflection families once, over a padded angular window
2  assemble the parameter vector from `refine`, each with physical bounds
3  trust-region least squares (scipy least_squares) on the weighted residual:
       r_i = sqrt(w_i) * (y_i_obs - y_i_calc(p))
   each evaluation recomputing:
       positions   from the dilated cell + zero shift
       |F|^2       with the current B_iso
       texture     March-Dollase
       profiles    Thompson-Cox-Hastings pseudo-Voigt
       background  Chebyshev, refined jointly
4  standard uncertainties from the Jacobian at the solution
5  report R_p, R_wp, R_exp, GoF, R_B, Durbin-Watson, and the residual curve
```

### 2.1 Every parameter is bounded, and that is deliberate

`_PARAMETER_BOUNDS` constrains each refined quantity, because **an unbounded
profile refinement will happily reach a lower $R_{wp}$ at a physically
impossible cell**. Two bounds worth quoting:

- `zero_shift_deg` $\in [-1, 1]$. A zero error beyond a degree is a broken
  diffractometer, not a refinable offset — and left free it swaps places with
  the cell, since both shift peak positions.
- `lattice_scale` $\in [0.9, 1.1]$. Ten per cent is far outside any real
  thermal or compositional cell change; a refinement that wants more has
  misidentified the phase.

These are not numerical guards. They encode the difference between a parameter
and an excuse.

### 2.2 The default refinement set, and why it is small

`DEFAULT_REFINEMENT_SET` is `scale`, `zero_shift_deg`, `lattice_scale`,
`caglioti_w` — plus the background. These are the parameters almost always worth
refining and least likely to run away *together*. The full
`REFINABLE_PARAMETERS` adds `caglioti_u`, `caglioti_v`, `lorentzian_y`,
`b_iso_overall`, `march_coefficient`.

Rietveld refinement is conventionally done **incrementally**: add parameters a
few at a time, watching the residual curve rather than only $R_{wp}$. Turning
everything on at once is the classic way to reach a low $R_{wp}$ at a meaningless
minimum, because strongly correlated parameters — zero shift against cell,
background against scale, $B_{\text{iso}}$ against scale — trade against each
other freely.

### 2.3 Do not subtract the background first

The background is refined **jointly**. Subtracting it beforehand removes the
correlation between background and scale that the parameter uncertainties depend
on, so the refinement then reports standard uncertainties that are too small.
The result is a number that looks better determined than it is, which is worse
than a slightly worse fit.

A high Chebyshev degree absorbs genuine broad features — including an amorphous
halo you may want to *see*. Degree is a modelling choice, not a knob to minimise.

## 3. Reading the result

| Statistic | Definition | Read it for |
| --- | --- | --- |
| $R_p$ | $\sum\lvert y_o - y_c\rvert / \sum y_o$ | unweighted profile agreement |
| $R_{wp}$ | $\sqrt{\sum w (y_o-y_c)^2 / \sum w y_o^2}$ | the quantity actually minimised |
| $R_{\text{exp}}$ | the value $R_{wp}$ would take if only counting statistics remained | the floor set by the data |
| **GoF** | $R_{wp}/R_{\text{exp}}$ | the honest headline; its square is reduced $\chi^2$ |
| $R_B$ | from Rietveld-partitioned integrated intensities | how well the *structure* fits, not the profile |
| **Durbin-Watson** | serial correlation of weighted residuals | near 2 uncorrelated; well below 1 means systematic misfit |

Returned as **fractions, not percentages**, so a reported $R_{wp}$ of $0.087$ is
$8.7\%$.

### 3.1 The two statistics that matter most are not $R_{wp}$

- **GoF** approaching 1 means the model explains the data down to counting
  noise. A GoF far below 1 means the weights are wrong, not that the fit is
  superb.
- **Durbin-Watson** catches what $R_{wp}$ hides. A refinement can reach a low
  $R_{wp}$ with residuals that are *systematically* wrong — a wrong peak shape
  leaves a sinusoidal residual through every peak — and DW well under 1 detects
  exactly that pattern of correlated misfit.

And above both: **the residual curve is the most informative single output of a
refinement**. Structure in it — a derivative-shaped residual at every peak (wrong
positions), a symmetric residual at every peak (wrong widths), a residual only at
high angle (wrong $B_{\text{iso}}$ or absorption) — diagnoses what to refine
next in a way no scalar does. `RietveldResult` carries
`residual_intensity` for this reason, and the workbench plots it beneath the
profile rather than as an optional extra.

## 4. Preferred orientation

March-Dollase models texture with one coefficient about one declared axis:
$r = 1$ is a random powder, $r \ne 1$ a preferred orientation of $(hkl)$ normals
along the specimen axis. `preferred_orientation_plane` is **required** before
`march_coefficient` may be refined, because a texture strength with no declared
axis is not a physical statement — and refining it blind lets it absorb intensity
errors from any source.

For real texture, a one-parameter model is a correction and not a measurement:
see {doc}`pole_figure_inversion`.

## 5. Constraints and failure modes

| Situation | What happens | What to do |
| --- | --- | --- |
| No atomic basis on the phase | only geometry and multiplicity contribute | supply a structure, or accept a profile-only fit |
| Background pre-subtracted | uncertainties too small | pass the raw profile |
| All parameters refined at once | plausible $R_{wp}$ at a meaningless minimum | refine incrementally |
| Widths refined with no instrument standard | specimen and instrument broadening confounded | calibrate `InstrumentBroadening` first |
| `converged=False` | ran out of evaluations | do not read the parameters |
| Wrong phase | bounds bite; `lattice_scale` pinned at 0.9 or 1.1 | the bound is telling you something |

## 6. Where it sits

Rietveld fits *everything at once*, which is its strength and its risk. When only
the cell is wanted, {doc}`precise_lattice_parameter_determination` extracts it
from peak positions with the systematic error extrapolated away — fewer
assumptions, and a better cell. Use Rietveld when phase fractions, structure or
profile parameters are the target.

## Verification

- The refinement recovering a known cell and scale from a synthetic pattern,
  and the R-factor definitions, in
  {doc}`../examples/generated/lattice-parameters`.

## See also

- {doc}`../theory/powder_xrd_and_saed` — the powder pattern forward model.
- {doc}`../theory/preferred_orientation_in_powder_intensities` — March-Dollase.
- {doc}`precise_lattice_parameter_determination` — the cell without the rest.
- {doc}`pole_figure_inversion` — texture properly measured.

## References

### Normative

- Rietveld, H. M. (1969). A profile refinement method for nuclear and magnetic
  structures. *Journal of Applied Crystallography* **2**, 65-71.
  <https://doi.org/10.1107/S0021889869006558>
- Thompson, P., Cox, D. E. & Hastings, J. B. (1987). Rietveld refinement of
  Debye-Scherrer synchrotron X-ray data from Al2O3. *Journal of Applied
  Crystallography* **20**, 79-83.
  <https://doi.org/10.1107/S0021889887087090>
- Dollase, W. A. (1986). Correction of intensities for preferred orientation in
  powder diffractometry: application of the March model. *Journal of Applied
  Crystallography* **19**, 267-272.
  <https://doi.org/10.1107/S0021889886089458>

### Informative

- Toby, B. H. (2006). R factors in Rietveld analysis: how good is good enough?
  *Powder Diffraction* **21**, 67-70.
  <https://doi.org/10.1154/1.2179804>
- Caglioti, G., Paoletti, A. & Ricci, F. P. (1958). Choice of collimators for a
  crystal spectrometer for neutron diffraction. *Nuclear Instruments* **3**,
  223-228. <https://doi.org/10.1016/0369-643X(58)90029-X>
- McCusker, L. B. et al. (1999). Rietveld refinement guidelines. *Journal of
  Applied Crystallography* **32**, 36-50.
  <https://doi.org/10.1107/S0021889898009856>
