# Correcting The Ghosts In A Pole-Figure ODF

An orientation distribution reconstructed from pole figures is missing a half of itself. Friedel's
law makes a diffraction pole figure blind to the sign of a plane normal, the forward operator
annihilates every odd-degree harmonic exactly, and the least-squares solution therefore returns the
even part with the odd part silently set to zero. That is not a neutral default: it puts false
maxima where the specimen is empty and depresses the true maxima to pay for them, which is what the
literature calls a *ghost*.

This page states how `pytex.texture.correct_ghosts` recovers an odd part, what each setting is
calibrated against, and where the method fails. The derivation is in
{doc}`../theory/ghost_problem_and_odd_harmonics`; the numbers are computed live in
{doc}`../examples/index`.

## 1. The pipeline

| Stage | Input | Output | Where it goes wrong |
| --- | --- | --- | --- |
| 1. Even solution | measured pole figures | even-degree coefficients | ill-posed if the data are fewer than the coefficients |
| 2. Odd basis | crystal and specimen symmetry, bandlimit | orthonormal odd functions on the quadrature | the symmetry may admit none at all |
| 3. Minimization | even density, odd basis | odd coefficients | the feasible set may be empty (truncation, not ghosts) |
| 4. Report | both parts | corrected ODF and its cost | a correction quoted without its cost reads as a measurement |

## 2. Stage 1 — the even part is the data, and is held fixed

`HarmonicODF.invert_pole_figures` solves

$$
\min_{\mathbf{c}} \; \lVert \mathbf{A}\mathbf{c} - \mathbf{p} \rVert^{2}
+ \lambda \lVert \mathbf{c} \rVert^{2},
$$

with $\mathbf{A}$ the pole-density response of each symmetry-projected even basis function at each
measured direction, and $\mathbf{p}$ the measured intensities in multiples of a random
distribution. Two properties of $\mathbf{A}$ matter here and are enforced rather than assumed:

- **it folds opposite normals** whenever the pole figure declares itself antipodal, so odd-degree
  functions lie exactly in its null space. Without the folding an odd function produces a visible
  predicted density and the correction could not claim to leave the fit alone;
- **it is normalized by the random level of the same folded kernel.** The folded random level is
  asymptotically twice the unfolded one, so mismatching them is a factor-of-two error in every
  density the ODF reports.

Everything downstream treats $\tilde f$, the even density on the quadrature, as fixed. It is what
the measurement determined, and a correction is not entitled to change it.

## 3. Stage 2 — the odd basis, and when there is none

The odd terms up to the bandlimit are enumerated, projected onto the crystal and specimen
symmetries by averaging over the group, and orthonormalized against the quadrature weights by the
eigendecomposition of the Gram matrix, with an eigenvalue floor (`basis_tolerance`, default
$10^{-10}$) discarding the numerically null directions.

**A symmetry admits odd terms only where it has an odd-degree invariant.** By character theory the
dimension of the degree-$\ell$ invariant subspace is the group average of the SO(3) character
$\chi_\ell(\theta) = \sin((\ell + \tfrac{1}{2})\theta)/\sin(\theta/2)$, which gives:

| Rotation group | First odd degree with an invariant |
| --- | --- |
| 432 (cubic) | 9 |
| 622 (hexagonal) | 7 |
| 222 (orthorhombic) | 3 |
| 1 (triclinic) | 1 |

So a cubic ODF expanded to degree 6 or 8 has **no ghost part to correct**, and the correction
returns an empty odd basis and says so, rather than reporting a correction of size zero as though
one had been made. This is checked live in the worked example
`ghost-cubic-first-odd-invariant-is-degree-nine`.

## 4. Stage 3 — the minimization

Write the corrected density on the quadrature as $f = \tilde f + \mathbf{O}\hat{\mathbf{c}}$, with
$\mathbf{O}$ the orthonormal odd basis. Define the *inadmissible part*

$$
v(f) = \begin{cases}
f & \text{inside a declared zero range},\\
\min(f, 0) & \text{elsewhere},
\end{cases}
$$

and minimize

$$
\Phi(\hat{\mathbf{c}}) = \tfrac{1}{2}\sum_q w_q \, v\big(f_q\big)^{2}
+ \tfrac{\mu}{2}\,\hat{\mathbf{c}}^{\mathsf{T}}\hat{\mathbf{c}} .
$$

```text
start with c = 0
repeat (L-BFGS-B, analytic gradient):
    f    = even_density + O @ c
    v    = violation(f)                        # min(f, 0), or f in the zero range
    Phi  = 0.5 * sum(w * v**2) + 0.5 * mu * c.c
    grad = O.T @ (w * v) + mu * c
until |grad| < tolerance or max_iterations
```

$\Phi$ is convex and continuously differentiable in $\hat{\mathbf{c}}$, so the minimizer is unique
and a quasi-Newton method reaches it in tens of iterations. It is the same point the classical
Dahms–Bunge alternating projection converges to — projection onto the non-negative densities, then
back onto the densities with the measured even part — without that iteration's long tail: measured
on the repository's demonstration case, alternating projection needed 2667 iterations to the same
answer the minimizer reaches in 4.

The cost of a correction is dominated by stage 2 rather than by this minimization. The half-angle
powers of the Wigner $d$ functions are tabulated once per basis evaluation and shared across every
term, which is what makes a degree-9 odd basis affordable: it cut that basis from 16.4 s to 2.9 s,
and a degree-9 ghost-corrected inversion in the workbench from 57 s to 14 s.

**Why the second term is not optional.** Positivity *bounds* the odd part; it does not determine
it. Once the density is admissible, every remaining direction in the odd subspace is free, and an
unregularized minimizer stops at whichever admissible point it happens to reach first. On the
demonstration case that point had a ghost amplitude ratio of 0.198 and a distance-to-truth of
0.137, against 0.145 and 0.0945 for the minimum-norm solution — a larger odd part than the data
force, reported as though the data forced it.

## 5. The settings, and what each is calibrated against

| Setting | Default | Calibrated against |
| --- | --- | --- |
| `method` | `"positivity"` | The constraint physics guarantees. `"zero_range"` additionally asserts that a range the data show as empty *is* empty, which is an assumption about the specimen. |
| `zero_range_threshold` | 0.05 m.r.d. | The density below which the even solution is read as declaring an empty range. Only meaningful for `"zero_range"`; a threshold of zero makes it weaker than plain positivity, not stronger. |
| `odd_regularization` ($\mu$) | $10^{-6}$ | Six orders below a typical violation norm, so it selects the minimum-norm solution without biasing admissibility. |
| `degree_bandlimit` | the ODF's own | An odd part resolved more finely than the even part it corrects would put detail into the answer that no data constrain. |
| `max_iterations` | 500 | Reaching it is reported, not raised: a correction stopped early is a valid lower bound provided the reader is told. |
| `tolerance` | $10^{-12}$ | Gradient tolerance. Read `infeasibility_after` for what the residual violation actually is. |
| `basis_tolerance` | $10^{-10}$ | Eigenvalue floor of the odd Gram matrix; below it the direction is numerically absent from the symmetry-projected span. |

## 6. Calibrated behaviour

Measured on a single broad orthorhombic component whose answer is known by construction — the case
pinned by `tests/unit/test_ghost_correction.py`, with a degree-4 expansion broad enough that
truncation is not in play:

| Quantity | Even-only solution | Positivity-corrected | True distribution |
| --- | --- | --- | --- |
| minimum density (m.r.d.) | −0.465 | ≈ 0 | 0.049 |
| maximum density (m.r.d.) | 3.74 | 4.24 | 4.06 |
| negative fraction of SO(3) | 9.3 % | 0.06 % | 0 |
| weighted distance to truth | 0.190 | 0.092 | — |
| ghost amplitude ratio | — | 0.145 | — |
| change in predicted pole densities | — | $8\times10^{-4}$ m.r.d. | — |

The last row is the check that matters: the measured intensities run from about 1 to 4 m.r.d., so
the correction moved the fit by 0.02 % — the quadrature error, and nothing more.

## 7. Failure modes

1. **Truncation masquerading as ghosts.** A texture too sharp for the bandlimit rings, and the
   ringing is negative. Positivity cannot repair it with odd terms of the same bandlimit, and the
   minimization then converges to a point that is still infeasible. `infeasibility_after` is what
   exposes this: it is near zero when the correction succeeded and stays comparable to
   `infeasibility_before` when the problem was never the ghosts.
2. **No odd terms at all.** Stage 2 returns an empty basis; the report says so and the ODF is
   returned unchanged.
3. **An under-determined even part.** If the even solution is itself a picture of the regularizer,
   correcting it produces a self-consistent picture of the regularizer. The correction says nothing
   about this; the reconstruction report's `matrix_rank`, `condition_number` and observation count
   do.
4. **Reading the odd part as a measurement.** It is an inference from positivity, and no
   pole-figure experiment can confirm or refute it. `describe()` says so on every report, and the
   corrected distribution is reached through `final_odf` rather than by replacing `odf`, so that a
   reader can always see what the data alone gave.

## Verification

- `tests/unit/test_ghost_correction.py` — 17 tests: the invisibility of odd harmonics to a
  Friedel-symmetric operator, the known-answer case above, mean-density preservation, the empty-odd
  basis for cubic symmetry, and the refusal to correct an ODF that already carries odd degrees.
- Worked examples `ghost-cubic-first-odd-invariant-is-degree-nine`,
  `ghost-correction-restores-a-non-negative-density` and
  `ghost-correction-leaves-the-measured-fit-untouched`, whose expected values are analytic
  identities and a cited standard result rather than recorded program output.

## See also

- {doc}`../theory/ghost_problem_and_odd_harmonics` — why the odd part is unmeasurable at all.
- {doc}`../theory/harmonic_odf_reconstruction` — the even solution this corrects.
- {doc}`../workflows/pole_figure_presentation` — drawing the result honestly.
- {doc}`../workflows/workbench_application` — the same correction in the application.

## References

- H.-J. Bunge, *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths (1969).
- M. Dahms and H.-J. Bunge, *The iterative series-expansion method for quantitative texture
  analysis. I. General outline*, Journal of Applied Crystallography **22** (1989) 439–447.
  DOI: <https://doi.org/10.1107/S0021889889005261>.
- S. Matthies, *On the reproducibility of the orientation distribution function of texture samples
  from pole figures (ghost phenomena)*, Physica Status Solidi (b) **92** (1979) K135–K138.
  DOI: <https://doi.org/10.1002/pssb.2220920253>.
- R. Hielscher and H. Schaeben, *A novel pole figure inversion method: specification of the MTEX
  algorithm*, Journal of Applied Crystallography **41** (2008) 1024–1037.
  DOI: <https://doi.org/10.1107/S0021889808030112>.
