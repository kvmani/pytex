# Residual Stress by the sin²ψ Method

**Surface:** `pytex.diffraction.xrd_residual_stress` —
`DiffractionElasticConstants`, `locate_stress_peak`, `fit_sin2psi_lines`,
`determine_residual_stress`, `residual_stress_pipeline`, `simulate_sin2psi_measurement`,
`parse_stress_scans`, `parse_stress_peak_positions`; in the workbench,
**XRD → Residual stress (sin²ψ)** (`xrd.residual_stress`), with its report drawn by
`pytex.app.services.xrd_stress_report`.

The sin²ψ method measures the macroscopic residual stress in the surface layer of a crystalline
specimen from the shift of one Bragg reflection as the specimen is tilted. At azimuth $\varphi$ and
tilt $\psi$ the planes normal to $\mathbf{m}(\varphi,\psi)$ diffract, and their strain is

$$
\varepsilon_{\varphi\psi} = \frac{d_{\varphi\psi}-d_{0}}{d_{0}}
= \tfrac{1}{2}S_{2}\,\mathbf{m}\cdot\boldsymbol{\sigma}\cdot\mathbf{m} + S_{1}\operatorname{tr}\boldsymbol{\sigma},
$$

so that under plane stress $d$ is linear in $\sin^{2}\psi$ with slope $d_{0}\tfrac{1}{2}S_{2}\sigma_{\varphi}$.
The derivations are in {doc}`../theory/residual_stress_sin2psi`; this page states how the result is
computed, what each stage of the report shows, and how to judge it.

```{figure} ../../figures/sin2psi_geometry.svg
:alt: The sin-squared-psi measurement geometry.
:width: 80%

The measurement geometry: specimen frame $S_1S_2S_3$, scattering vector $\mathbf{m}$ at azimuth
$\varphi$ and tilt $\psi$, and the diffracting planes normal to it.
```

## 1. Algorithmic architecture

| Stage | Function | Output | Purpose |
| --- | --- | --- | --- |
| 1. Reference | `DiffractionElasticConstants`, the phase cell | $d_{0}$, $u(d_{0})$, $S_{1}$, $\tfrac{1}{2}S_{2}$ | The stress-free spacing and the constants of the reflection |
| 2. Peak location | `locate_stress_peak` | $2\theta$, $u(2\theta)$ per $(\varphi,\psi)$ | LPA correction; Kα1/Kα2 pseudo-Voigt fit, or parabola/centroid after Rachinger stripping |
| 3. Spacing and strain | `determine_residual_stress` | $d$, $u(d)$, $\varepsilon$, $u(\varepsilon)$ | Bragg's law and $u(d) = d\cot\theta\,u(\theta)$ |
| 4. Per-azimuth lines | `fit_sin2psi_lines` | slope, $\sigma_{\varphi}$, $\tau_{\varphi}$, curvature $t$ | The classical evaluation and its linearity and ψ-splitting tests |
| 5. Tensor | `determine_residual_stress` | $\sigma_{ij}$, covariance, $\chi^{2}_{\nu}$ | One weighted linear least-squares over all measurements |
| 6. Budget | `StressTensorFit` | statistical, $d_{0}$, elastic-constant parts; Monte Carlo | GUM combination and its Monte Carlo check |
| 7. Derived | `StressTensorFit.in_plane_principal`, `von_mises` | $\sigma_{\mathrm{I}}, \sigma_{\mathrm{II}}, \varphi_{\mathrm{I}}$, $\sigma_{\mathrm{vM}}$ | Propagated by the Jacobian of each function |

## 2. Stage 1: the stress-free spacing and the elastic constants

$d_{0}$ is computed from the phase's cell with the stress-free lattice parameter(s) $a_{0}$ (and
$c_{0}$) substituted, through the reciprocal metric: $d_{0} = 1/\lvert\mathbf{B}^{*}\mathbf{h}\rvert$. Its
uncertainty is taken as $u(d_{0})/d_{0} = u(a_{0})/a_{0}$. With **Determine d₀ from the data** it is
refined instead (section 6).

The diffraction elastic constants come from single-crystal stiffness (tabulated for ferrite,
nickel, aluminium, copper, tungsten and α-titanium, or entered as $C_{11}$, $C_{12}$, $C_{44}$) under
the Kröner (default), Reuss, Voigt or Neerfeld–Hill model, or from $(E,\nu)$, or entered directly.
The relative uncertainty declared for them enters the budget. The **Diffraction elastic
constants** stage plots $\tfrac{1}{2}S_{2}$ of low-index reflections under all four models against
the cubic orientation parameter $\Gamma$ and marks the one used; the Reuss–Voigt spread at that
reflection is a fair guide to $u(S)/S$.

## 3. Stage 2: locating the peaks

In each scan the window of width **Fit window** about the stress-free Bragg angle is used. With the
LPA correction on, counts and their Poisson uncertainties are divided by
$(1+\cos^{2}2\theta)/\sin^{2}\theta\times(1-\tan\psi\cot\theta)$ (the last factor under
$\omega$-tilting only), normalized at the window centre.

- **Pseudo-Voigt profile fit** (default): Kα1 + Kα2 at the Bragg-law separation with the tabulated
  ratio, shared width and mixing, linear background; $u(2\theta)$ from the covariance scaled by the
  window's $\chi^{2}_{\nu}$.
- **Parabola**: weighted fit through the points above 80 % of the Kα2-stripped net maximum;
  vertex $-b/2c$ and its propagated uncertainty.
- **Centroid**: first moment of the Kα2-stripped net profile over a continuous window at half
  maximum, centred on itself; uncertainty with the fixed-point factor $1/(1-\partial F/\partial c)$.

The Kα2 stripping is essential for the last two at stress angles: an unstripped doublet biases a
parabola by more than 100 MPa on ferrite (211) with Cr Kα. The profile fit models the doublet
instead and is the only one of the three that is unbiased over many noise realizations; the others
are cross-checks.

## 4. Stages 3 and 4: strains and the per-azimuth lines

$d = \lambda_{1}/(2\sin\theta)$, $u(d) = d\cot\theta\,u(\theta)$, $\varepsilon = d/d_{0}-1$. At each
azimuth, weighted least squares of $d = c_{0}+c_{1}\sin^{2}\psi\,[+c_{2}\sin2\psi]$ (the last term when
both signs of $\psi$ were measured) gives $\sigma_{\varphi}-\sigma_{33} = c_{1}/(d_{0}\tfrac{1}{2}S_{2})$ and
$\tau_{\varphi} = c_{2}/(d_{0}\tfrac{1}{2}S_{2})$, with uncertainties scaled by $\sqrt{\chi^{2}_{\nu}}$ when that
exceeds one. A refit with an added $\sin^{4}\psi$ term tests linearity.

## 5. Stages 5 to 7: the tensor, its budget and derived quantities

All measurements enter $\boldsymbol{\varepsilon} = \mathbf{A}\mathbf{x}$ at once, with one row per
$(\varphi,\psi)$ and one column per free component of the chosen stress state:

| Stress state | Free components | Needs |
| --- | --- | --- |
| Biaxial | $\sigma_{11}, \sigma_{22}, \sigma_{12}$ | three azimuths not $180^{\circ}$ apart |
| Biaxial with shear | $+\,\sigma_{13}, \sigma_{23}$ | tilts of both signs (identifiable without, but poorly) |
| Triaxial | all six | an exact $d_{0}$; $d_{0}$ cannot be refined |

Identifiability is checked from the singular values of $\mathbf{A}$ before solving. The budget is
$\mathbf{V} = \mathbf{V}_{\mathrm{int}}\max(1,\chi^{2}_{\nu}) + \mathbf{s}_{d_{0}}\mathbf{s}_{d_{0}}^{\mathsf T}u^{2}(d_{0}) + \mathbf{V}_{S}$,
cross-checked by redrawing every input and re-solving in one batched solve.

## 6. Determining d₀ under plane stress

$d_{k} = d_{0} + \mathbf{a}_{k}\cdot(d_{0}\mathbf{x})$ is linear in $(d_{0}, d_{0}\mathbf{x})$; one solve gives
both. The information comes from making the intercepts agree with the slopes through $S_{1}$ — the
strain-free-direction method made exact. Refused for the triaxial state, where $\sigma_{33}$ and
$d_{0}$ are inseparable.

## 7. Reading the workbench report

The result opens with the stress and how far to trust it, then the evidence, the diagnostics and
the method — the order of `pytex.app.results.REPORT_SECTIONS`. **Report + figures** downloads
the whole of it: `report.md`, every figure as `figures/<key>.svg`, and `result.json` with every
number at full precision, including the library result under `data.library`.

### Warnings

Each warning is a condition the evaluation found in its own data:

- strain-fit $\chi^{2}_{\nu} > 3$ (unexplained scatter; the statistical uncertainty has already been
  widened) or $< 0.3$ (overstated peak uncertainties);
- two components correlated beyond $\pm0.95$ **by the measurement design** (judged on the
  statistical covariance, since $d_{0}$ correlates $\sigma_{11}$ and $\sigma_{22}$ by construction);
- $u(d_{0})$ dominating the budget;
- curvature of $d$ in $\sin^{2}\psi$ beyond $3u$ at any azimuth;
- significant ψ-splitting in a biaxial evaluation (choose the one with shear);
- tilts reaching only $\sin^{2}\psi < 0.3$; a reflection below $2\theta = 120^{\circ}$;
- peak locations that did not converge; a generated rather than measured data set.

### Result — the stress tensor

The table lists every free component with its combined standard uncertainty and the part from
each source, the Monte Carlo value beside them, and, for the demonstration, the stress it was
generated with. **Normal stress along each azimuth** plots $\sigma_{\varphi}$ from each slope against
the curve $\sigma_{11}\cos^{2}\varphi + \sigma_{12}\sin2\varphi + \sigma_{22}\sin^{2}\varphi$ of the tensor with its
uncertainty band; **Mohr's circle** shows the in-plane state and its principal stresses.

### Result — uncertainty budget

Bars per component and source. The tallest coloured bar is the one worth reducing; the hatched
Monte Carlo bar matching the black combined bar confirms the linear propagation.

### Evidence — the measurement

The scans as measured, the **measured directions** on the specimen hemisphere (open markers are
$\psi<0$, drawn where that direction points), and **the peak moves with tilt**: the normalized
profiles of each azimuth coloured by $\psi$, the raw signal of the method.

### Evidence — peak positions

Every located peak with $u(2\theta)$, width, height and peak-fit $\chi^{2}_{\nu}$; **every scan and its
located peak** (a grid, one column per azimuth and one row per tilt, the fitted profile drawn on the
counts on the scale it was fitted) and **peak width, height and precision against tilt**. The
peak-fit $\chi^{2}_{\nu}$ judges one profile; it is not the strain-fit $\chi^{2}_{\nu}$.

### Evidence — d against sin²ψ at each azimuth

The table of per-azimuth lines, the **d against sin²ψ** figure the method is named for (filled
$\psi\ge0$, open $\psi<0$, the weighted line in orange, the tensor's prediction dashed for each branch),
and **lattice strain against sin²ψ** for all azimuths together, whose absolute level carries $d_{0}$.

### Diagnostics

**Residuals of the stress fit** in units of $u$ against $\sin^{2}\psi$ with $\pm2$ and $\pm3$ bands, the
**correlation** matrix, and **linearity and ψ-splitting** with the branch mean $a_{1}$ and
half-difference $a_{2}$ plots.

### Method

The theory, the elastic constants (with the model comparison figure) and the algorithm as steps,
so that the downloaded report is self-contained.

## Verification

- `tests/unit/test_xrd_residual_stress.py` — the library, against closed forms and exact data.
- `tests/unit/test_app_xrd_stress.py` — the operation: recovery of the demonstration's generating
  stress, the report's derived quantities (normalized residuals reproduce $\chi^{2}_{\nu}$, the tensor
  curve is $\sigma_{\varphi}$, the budget combines in quadrature and agrees with Monte Carlo), every
  input route and refusal, and the contents of the downloadable bundle.
- `tests/unit/test_residual_stress_figures.py` — the geometry figure is the generator's output.
- Executable worked examples: {doc}`../examples/generated/residual-stress`.

## See also

- {doc}`../theory/residual_stress_sin2psi` — the derivations.
- {doc}`precise_lattice_parameter_determination` — precise spacings, and the stress-free $d_{0}$.
- {doc}`elastic_homogenization` — the Voigt, Reuss and Hill aggregate moduli.

## References

### Normative

- Macherauch, E. & Müller, P. (1961). *Zeitschrift für angewandte Physik* **13**, 305–312.
- Noyan, I. C. & Cohen, J. B. (1987). *Residual Stress: Measurement by Diffraction and
  Interpretation*. Springer. <https://doi.org/10.1007/978-1-4613-9570-6>
- Dölle, H. (1979). *Journal of Applied Crystallography* **12**, 489–501.
  <https://doi.org/10.1107/S0021889879013169>
- Welzel, U. et al. (2005). *Journal of Applied Crystallography* **38**, 1–29.
  <https://doi.org/10.1107/S0021889804029516>
- JCGM 100:2008 and JCGM 101:2008 — the GUM and its Monte Carlo supplement.

### Informative

- Kröner, E. (1958). *Zeitschrift für Physik* **151**, 504–518.
  <https://doi.org/10.1007/BF01337948>
- Fitzpatrick, M. E. et al. (2005). *Determination of Residual Stresses by X-ray Diffraction*,
  NPL Good Practice Guide No. 52, issue 2.
- SAE HS-784 (2003). *Residual Stress Measurement by X-Ray Diffraction*.
- Hauk, V. (ed.) (1997). *Structural and Residual Stress Analysis by Nondestructive Methods*.
  Elsevier. <https://doi.org/10.1016/B978-0-444-82476-9.X5000-2>
