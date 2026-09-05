# Determining A Lattice Parameter Precisely

**Surface:** `pytex.diffraction.xrd_peaks.detect_and_fit_peaks`,
`pytex.diffraction.xrd_indexing.index_peaks`,
`pytex.diffraction.xrd_lattice_parameter.determine_lattice_parameters`,
`determine_lattice_parameters_le_bail`,
`determine_lattice_parameters_from_pattern`, with
`pytex.diffraction.xrd_corrections` for the aberrations and the display
transforms.

The derivations are in
{doc}`../theory/precise_lattice_parameter_determination`. This page is the
implementation: the stages in order, the settings that matter with the values
they take, the numbers each stage produces on cases with known answers, and the
ways it fails.

## 1. The pipeline, and where each stage can go wrong

| Stage | Surface | Produces | Fails as |
| --- | --- | --- | --- |
| 1. Detect | `detect_peaks` | candidate angles | too many candidates (threshold low), or none (threshold high, or wrong expected width) |
| 2. Fit | `fit_peaks` | positions **with ESDs** | reduced $\chi^2 \gg 1$: wrong shape, wrong weights, or a neighbour inside the window |
| 3. Index | `index_peaks` | $(hkl)$ per peak, $M_N$, $F_N$ | reflections dropped when the starting cell is wrong; a strong unindexed peak means a second phase |
| 4. Determine | `determine_lattice_parameters` | cell + ESDs + drift term | reduced $\chi^2 \gg 1$: an aberration remains, or the wrong $f(\theta)$ |
| 4'. Determine | `determine_lattice_parameters_le_bail` | cell + refined aberration | reduced $\chi^2 \gg 1$: the doublet, the weights, or the background |

`determine_lattice_parameters_from_pattern` runs 1–4 and **repeats 3–4**, which
is not optional; see §5.

## 2. Stage 1 — detection

```
background  <- SNIP(pattern, half_window_deg)
y           <- 2 sqrt(max(I - background, 0) + 3/8)     if unit is counts
y           <- y - median(y)
for each scale s in geomspace(s_min, s_max, 8):
    k_s     <- ricker(s); k_s <- k_s - mean(k_s); k_s <- k_s / ||k_s||_2
    R_s     <- convolve(y, k_s)
R(i)        <- interpolate R_s(i) at s(i) = FWHM(2*theta_i) / (step * 2 sqrt(2 ln 2))
noise       <- 1.4826 * MAD(R)
candidates  <- local maxima of R with R > prominence_sigma * noise
            -> parabolic sub-step refinement on R
            -> thin: reject within one FWHM of a stronger peak,
                     and at the K-alpha2 partner angle of a stronger peak
```

Three properties are worth stating because each removes a free parameter:

- **Zero-mean kernel** — the response ignores a constant offset, and being even
  it nearly ignores a linear ramp.
- **Unit $L^2$ norm** — white noise of unit variance gives a response of unit
  variance, so `prominence_sigma` is in noise standard deviations, not counts.
- **Anscombe transform** — Poisson noise becomes homoscedastic, so one threshold
  works across a pattern whose peaks differ by decades.

| Setting | Default | Choose it by |
| --- | --- | --- |
| `prominence_sigma` | 5.0 | 5 is conservative; drop to 3–4 to chase weak lines and expect more rejected at stage 2 |
| `expected_fwhm_deg` / `instrument` | 0.12° fallback | only needs to be right within about a factor of two |
| `background_half_window_deg` | 2.0° | wider than the broadest peak, narrower than the background's curvature |
| `suppress_kalpha2` | `True` | leave on unless inspecting raw detection |

## 3. Stage 2 — fitting

Each candidate gets a window of `window_fwhm` expected widths either side, fitted
by bounded trust-region least squares over
$[c,\ h,\ w_L,\ (w_R,)\ \eta,\ b_0,\ b_1]$:

$$y(2\theta) = h\,P(2\theta; c, w, \eta)
+ h\,r\,P\!\left(2\theta; c_2(c), w, \eta\right)
+ b_0 + b_1 (2\theta - c)$$

with $c_2$ **not free**: $\sin\theta_2 = (\lambda_2/\lambda_1)\sin\theta_1$.
Weights are $1/\sigma$ from the pattern's own uncertainties, else Poisson, else
unity. The reported position ESD is
$\sqrt{[(J^{\mathsf T}WJ)^{-1}]_{11}}$ scaled by the window's reduced $\chi^2$.

Fits that converge to the same angle within `merge_tolerance_fwhm` (default 0.25)
of their mean FWHM are merged, keeping the lower reduced $\chi^2$. Two fits at
one angle are one reflection, and leaving both would let the indexer give the
copies different $(hkl)$ — a reflection at a spacing nothing diffracted from.

**Calibrated behaviour.** On a Poisson-noised synthetic Ni pattern
(30 000 peak counts, 0.12° FWHM, 0.01° step, Cu K$\alpha_1$/K$\alpha_2$):

| Quantity | Value |
| --- | --- |
| Positions recovered | all 7, within $10^{-3}$ degrees of the generating reflection list |
| Position ESD | $\approx 2 \times 10^{-4}$ degrees |
| Reduced $\chi^2$ | 0.93 – 1.32 |
| Cost of *not* modelling the doublet | every centre displaced $> 5$ millidegrees upward |

## 4. Stage 3 — indexing

Reflections are enumerated from the candidate phase over the measured range
padded by $\max(2\,\varepsilon, 1^{\circ})$, filtered at
`minimum_relative_intensity` (default 0.001) of the strongest, and matched by
`scipy.optimize.linear_sum_assignment` on $|2\theta_{\text{obs}} -
2\theta_{\text{calc}}|$ with out-of-tolerance pairs given a prohibitive cost.

Forbidden pairs are made *expensive*, not infinite: the rectangular assignment
needs a feasible complete solution, which an infinity makes inexpressible.

`tolerance_deg` (default 0.3°) must be **wider than any uncorrected zero or
displacement error and narrower than the spacing between neighbouring calculated
lines**. Those two constraints can conflict on a low-symmetry pattern; when they
do, use Le Bail.

## 5. Re-indexing is part of the algorithm

A starting cell wrong by a fraction $e$ misplaces a reflection by

$$\Delta(2\theta) = 2 e \tan\theta$$

which diverges towards back-reflection. Worked: $e = 3 \times 10^{-3}$ — an
ordinary discrepancy between a tabulated cell and a real alloy — puts a
$121^{\circ}$ reflection $0.6^{\circ}$ out, far outside any usable tolerance.
Those reflections are silently dropped, and they carry almost all the precision.

Measured on the workbench demonstration scan (cell dilated by 1.003):

| Passes | Reflections indexed | $a$ error |
| --- | --- | --- |
| 1 | 3 of 6 | $+9.4\ \mu\text{Å}$ (on 3 low-angle lines) |
| 2 | 6 of 6 | $-12.0\ \mu\text{Å}$ |
| 3 | 6 of 6 (converged) | $-12.0\ \mu\text{Å}$ |

The angular floor `minimum_two_theta_deg` is applied **after** the passes
converge, never during them: the first pass indexes the low-angle reflections,
which are exactly the ones a floor discards, so applying both at once can leave
nothing at all.

## 6. Stage 4 — Cohen least squares

```
X[i, :n]    <- (lambda^2 / 4) * [h^2 k^2 l^2 2hk 2hl 2kl]_i @ C        # C: 6 x n
X[i, n]     <- sin^2(theta_i) * f(theta_i)                             # if f != none
y[i]        <- sin^2(theta_i)
w[i]        <- 1 / (sin(2 theta_i) sigma(2 theta_i) / 2)
p           <- lstsq(X * w[:, None], y * w)
chi2_red    <- sum(((y - X p) / sigma)^2) / (m - n - 1)
cov         <- inv(X^T W X) * chi2_red
G*          <- C @ p[:n];  G <- inv(G*);  a = sqrt(G_11), cos(alpha) = G_23/(bc)
sigma(cell) <- sqrt(diag(J cov J^T)),  J numerically differentiated
```

The constraint matrix $C$ per system:

| System | $n$ | Free parameters |
| --- | --- | --- |
| cubic | 1 | $a^{*2}$ |
| tetragonal, trigonal, hexagonal | 2 | $a^{*2}, c^{*2}$ |
| orthorhombic | 3 | $a^{*2}, b^{*2}, c^{*2}$ |
| monoclinic | 4 | $+\ a^{*}c^{*}\cos\beta^{*}$ |
| triclinic | 6 | all of $\mathbf{G}^{*}$ |

Choosing $f(\theta)$:

| Dominant aberration | Use | Note |
| --- | --- | --- |
| Detector zero | `cot_theta` | fitted $D$ **is** the zero, in radians |
| Specimen displacement | `cos_squared_over_sin` | exact for Bragg–Brentano |
| Both / unknown | `nelson_riley` | the usual compromise |
| Camera absorption | `bradley_jay` | Cohen's classical column |
| None | `none` | run it to see what the correction was worth |

**Calibrated behaviour.** Synthetic Ni, 100 µm specimen displacement injected:

| Method | $f(\theta)$ | Relative error | Reduced $\chi^2$ |
| --- | --- | --- | --- |
| average | — | $4.0 \times 10^{-4}$ | $2.5 \times 10^{5}$ |
| cohen | none | $8.7 \times 10^{-5}$ | $3.7 \times 10^{4}$ |
| cohen | bradley_jay | $4.0 \times 10^{-5}$ | $4.3 \times 10^{3}$ |
| cohen | nelson_riley | $4.8 \times 10^{-6}$ | 12.9 |
| cohen | cos_squared_over_sin | $9.6 \times 10^{-8}$ | **0.96** |

The reduced $\chi^2$ column is the one to read on real data, where no known
answer exists: it tracks the accuracy across four orders of magnitude.

## 7. Stage 4' — Le Bail

```
observed <- raw - SNIP(raw)                       # not clipped at zero
P[k, :]  <- pseudo_voigt(centre_k) + r * pseudo_voigt(kalpha2 partner of k)
P[k, :]  <- P[k, :] / sum(P[k, :])                # unit sum: I_k is INTEGRATED
I[k]     <- sum(observed) / K                     # equal start
repeat `cycles` times:
    share    <- I[:, None] * P / sum_j (I[j] * P[j])          # Le Bail partition
    I        <- share @ max(observed - background_line, 0)
    minimise sum(((I @ P + background_line) - observed) * w)^2
             over [cell (n), systematic (1), U, V, W, eta, bg0, bg1]
             with w = 1 / sqrt(raw)  and  x_scale = "jac"
```

Four points are load-bearing, and each produced a **wrong cell** rather than a
poor fit when got wrong:

1. **Model the whole K$\alpha$ multiplet.** Alpha1 only, against doublet data,
   cost 39 µÅ on a clean pattern.
2. **Weight from the measured counts**, never the subtracted profile.
   Subtraction removes signal, not variance.
3. **Normalise each profile to unit sum.** The partition sums counts, so it
   returns an integrated intensity; using it as a unit-height amplitude
   overstated peaks about twentyfold at a 0.01° step.
4. **Refine a stiff residual background.** SNIP removes the background's shape
   and leaves a level offset; unmodelled, it was ~90% of the total misfit, all
   of it between the peaks.

`systematic` refines **exactly one** of zero or displacement. Refining both from
one specimen scan is ill-conditioned: they differ only as constant against
$\cos\theta$.

**Calibrated behaviour.** 100 µm displacement injected, `systematic="displacement"`:

| Phase | Refined displacement | $a$ error | $c$ error | Reduced $\chi^2$ |
| --- | --- | --- | --- | --- |
| Ni (cubic) | 0.10048 mm | $-3.2\ \mu\text{Å}$ | — | 1.15 |
| Ti (hexagonal) | 0.09996 mm | $-2.2\ \mu\text{Å}$ | $-1.4\ \mu\text{Å}$ | 1.23 |
| Ni, term omitted | — | $+228\ \mu\text{Å}$ | — | **18.9** |
| Ti, term omitted | — | $+142\ \mu\text{Å}$ | $+831\ \mu\text{Å}$ | **32.1** |

Returning the aberration in millimetres is a stronger check than a plausible
cell: a wrong model can reach a right-looking cell by compensating errors, but
will not also reproduce a physical quantity nobody fitted for.

## 8. When to use which

| Situation | Use |
| --- | --- |
| Cubic, well-resolved peaks | `cohen` with the matching $f$ |
| Hexagonal or lower, resolvable peaks | `cohen` — the joint solution is the only kind available |
| Peaks overlap badly | `le_bail` |
| Teaching the failure of averaging | `average`, cubic only |
| Both a zero **and** a displacement present | calibrate the zero against a standard, then refine displacement only |

On alpha-uranium (orthorhombic, heavily overlapped) the two routes agree on the
cell to $10^{-4}$ while their goodness of fit does not: 75 for the
peak-position route against 1.7 for the whole-pattern one. That gap is what Le
Bail is for.

## 9. Constraints and limits

- **`average` is cubic-only, by construction.** Outside the cubic system a
  lattice parameter *per reflection* does not exist; the surface raises rather
  than returning a number.
- **One drift term removes one aberration.** A pattern carrying both a zero and
  a displacement cannot be fully corrected by any single $f$.
- **Le Bail intensities are extracted, not measured.** For two completely
  overlapped reflections the partition is whatever ratio the iteration started
  with. They describe the profile; they are unfit for structural work. This is
  also exactly why the method cannot be biased by texture.
- **$R_{wp}$ here is computed on the background-subtracted profile** and is
  systematically higher than a Rietveld program's. The two must not be compared.
- **This is not a stress.** A symmetric scan measures one strain component,
  normal to the specimen surface. See §8 of the theory note.
- **Unknown-cell autoindexing is not implemented.** ITO, TREOR and DICVOL solve
  a different and much harder problem.

## Verification

`tests/unit/test_xrd_peaks.py`, `test_xrd_corrections.py`,
`test_xrd_indexing.py`, `test_xrd_lattice_parameter.py` and
`test_app_xrd.py`. Every accuracy assertion compares against the cell that
generated the pattern or the aberration that was injected, never against a
stored output. Worked examples with independent provenance are in
{doc}`../examples/generated/lattice-parameters`; the tutorial is
{doc}`../tutorials/notebooks/34_precise_lattice_parameters`.

## See also

- {doc}`../theory/precise_lattice_parameter_determination` — the derivations.
- {doc}`../theory/powder_xrd_and_saed` — the forward simulation this consumes.

## References

Cullity, B. D. & Stock, S. R., *Elements of X-Ray Diffraction*, 3rd ed.,
Prentice Hall (2001), Ch. 11.

Cohen, M. U., *Rev. Sci. Instrum.* **6** (1935) 68,
[doi:10.1063/1.1751937](https://doi.org/10.1063/1.1751937).

Nelson, J. B. & Riley, D. P., *Proc. Phys. Soc.* **57** (1945) 160,
[doi:10.1088/0959-5309/57/3/302](https://doi.org/10.1088/0959-5309/57/3/302).

Le Bail, A., Duroy, H. & Fourquet, J. L., *Mater. Res. Bull.* **23** (1988) 447,
[doi:10.1016/0025-5408(88)90019-0](https://doi.org/10.1016/0025-5408(88)90019-0).

de Wolff, P. M., *J. Appl. Crystallogr.* **1** (1968) 108,
[doi:10.1107/S002188986800508X](https://doi.org/10.1107/S002188986800508X).

Smith, G. S. & Snyder, R. L., *J. Appl. Crystallogr.* **12** (1979) 60,
[doi:10.1107/S002188987901178X](https://doi.org/10.1107/S002188987901178X).

Wilson, A. J. C., *Mathematical Theory of X-ray Powder Diffractometry*, Philips
Technical Library (1963).
