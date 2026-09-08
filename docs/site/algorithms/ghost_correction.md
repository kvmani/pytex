# Ghost Correction in Harmonic ODF Inversion

**Surface:** `pytex.texture.ghosts.correct_ghosts`, `HarmonicODFInversionReport`,
`pytex.texture.harmonics.HarmonicODF`, and workbench operation
`texture.measured_pole_figures`.

In kinematic diffraction, Friedel's law ($I_{\mathbf{h}} = I_{-\mathbf{h}}$) renders
diffraction pole figures intrinsically centrosymmetric. Consequently, experimental
pole figures represent projections of the centrosymmetric (even-degree) portion of
the orientation distribution function (ODF), leaving the odd-degree harmonic
components entirely unconstrained by diffraction data. Inverting experimental pole
figures with odd coefficients set to zero introduces systematic mathematical artifacts
termed **ghosts** (Matthies, 1979): artificial local maxima arise in regions where the
true orientation density is zero, while genuine texture peaks are proportionally
attenuated to conserve unit integral probability.

This page explains the mathematical basis of the ghost phenomenon, describes how
`pytex.texture.correct_ghosts` reconstructs the unmeasured odd harmonic subspace
by enforcing physical non-negativity ($f(g) \ge 0$), specifies the role of numerical
regularization parameters, and delineates the physical boundary conditions where
reconstruction succeeds or fails.

## 1. Algorithmic architecture

The ghost correction workflow operates as a decoupled post-processing stage
following even-degree harmonic inversion:

| Stage | Input Data | Output Structure | Governing Mathematical Condition |
| --- | --- | --- | --- |
| 1. Even-degree inversion | Measured pole figures $P_{hkl}$ | Even harmonic coefficients $\tilde{\mathbf{c}}$ | Least-squares fit to centrosymmetric projection operator |
| 2. Odd basis construction | Crystal and specimen symmetry, $L$ | Orthonormalized odd basis $\mathbf{O}$ | Group character projection and Gram matrix truncation |
| 3. Constrained optimization | Even density $\tilde{f}(g)$, odd basis $\mathbf{O}$ | Odd coefficients $\hat{\mathbf{c}}$ | Convex penalized minimization enforcing $f(g) \ge 0$ |
| 4. Solution reporting | Even and odd coefficients | Corrected ODF and diagnostic metrics | Diagnostic residual and admissibility verification |

## 2. Mathematical formulation

### 2.1 Decoupling the even and odd subspaces

An orientation distribution function expanded in generalized spherical harmonics
$T_\ell^{\mu\nu}(g)$ decomposes orthogonally into even- and odd-degree components:

$$
f(g) = \tilde{f}(g) + f_{\text{odd}}(g)
= \sum_{\substack{\ell=0 \\ \ell \text{ even}}}^{L} \sum_{\mu,\nu} C_\ell^{\mu\nu}\,\dot{T}_\ell^{\mu\nu}(g)
+ \sum_{\substack{\ell=1 \\ \ell \text{ odd}}}^{L} \sum_{\mu,\nu} C_\ell^{\mu\nu}\,\dot{T}_\ell^{\mu\nu}(g).
$$

Because pole figures generated under Friedel symmetry satisfy $P_{hkl}(\mathbf{y}) = P_{\bar{h}\bar{k}\bar{l}}(\mathbf{y})$,
the forward projection operator $\mathbf{A}$ annihilates $f_{\text{odd}}(g)$ identically:
$\mathbf{A} f_{\text{odd}} = \mathbf{0}$. Therefore, standard regularized least-squares
inversion of experimental pole figures minimizes $\lVert \mathbf{A}\tilde{\mathbf{c}} - \mathbf{p} \rVert^2$,
which determines the even coefficients $\tilde{\mathbf{c}}$ while providing zero gradient
information for odd degrees. The even component $\tilde{f}(g)$ is retained as fixed
experimental input during ghost correction.

### 2.2 Group character projection and symmetry invariants

Odd-degree harmonic functions must satisfy the point-group symmetries of both the
crystal ($\mathcal{G}_{\text{xtal}}$) and the specimen ($\mathcal{G}_{\text{spec}}$).
By character theory, the dimension $d_\ell$ of the invariant subspace at harmonic
degree $\ell$ equals the group average of the $\mathrm{SO}(3)$ rotation character:

$$
d_\ell = \frac{1}{|\mathcal{G}|} \sum_{g \in \mathcal{G}} \frac{\sin\left[(\ell + \frac{1}{2})\omega(g)\right]}{\sin\left(\frac{1}{2}\omega(g)\right)},
$$

where $\omega(g)$ denotes the rotation angle of symmetry operation $g$. For groups
possessing high rotational symmetry, the lowest permissible odd-degree invariant is
substantially greater than unity:

| Point Group Symmetry | Minimum Degree $\ell_{\text{min}}$ for Non-Trivial Odd Invariants |
| --- | --- |
| Cubic (432, $O$) | 9 |
| Hexagonal (622, $D_6$) | 7 |
| Orthorhombic (222, $D_2$) | 3 |
| Triclinic (1, $C_1$) | 1 |

Consequently, for cubic materials reconstructed at series bandlimits $L < 9$, the
odd harmonic invariant subspace is identically empty ($d_\ell = 0$ for all odd $\ell \le 8$).
In such cases, `correct_ghosts` reports an empty basis and returns the distribution
without alteration, preserving exact group-theoretical invariants.

### 2.3 Penalized optimization objective

Let $\mathbf{O} \in \mathbb{R}^{Q \times M}$ represent the orthonormalized basis of
admissible odd harmonics evaluated on a discrete quadrature grid of $Q$ orientations
with positive integration weights $w_q$. The total density is parameterized as
$\mathbf{f} = \tilde{\mathbf{f}} + \mathbf{O}\hat{\mathbf{c}}$, where $\hat{\mathbf{c}} \in \mathbb{R}^M$
are the odd-basis expansion coefficients.

To enforce physical admissibility while resolving the indeterminacy of the odd
subspace, PyTex defines an infeasibility penalty function $v(f)$:

$$
v(f) = \begin{cases}
f, & \text{inside a declared zero-range domain}, \\
\min(f, 0), & \text{elsewhere (positivity constraint)}.
\end{cases}
$$

The optimization objective balances constraint satisfaction against a Tikhonov
minimum-norm regularizer:

$$
\Phi(\hat{\mathbf{c}}) = \frac{1}{2}\sum_{q=1}^Q w_q \left[v(f_q)\right]^2 + \frac{\mu}{2}\,\lVert\hat{\mathbf{c}}\rVert_2^2.
$$

The objective function $\Phi(\hat{\mathbf{c}})$ is convex and continuously differentiable
($C^1$) with analytical gradient:

$$
\nabla \Phi(\hat{\mathbf{c}}) = \mathbf{O}^{\mathsf{T}}\left(\mathbf{w} \odot v(\mathbf{f})\right) + \mu\hat{\mathbf{c}},
$$

where $\odot$ denotes elementwise multiplication. The system is solved using the
quasi-Newton L-BFGS-B algorithm, converging within tens of iterations. This formulation
avoids the slow asymptotic convergence of classical alternating projection methods
(Dahms & Bunge, 1989).

The regularization parameter $\mu > 0$ ensures a strictly convex problem and selects
the minimum-norm odd distribution among all admissible non-negative solutions,
preventing the introduction of arbitrary high-amplitude odd fluctuations unconstrained
by the non-negativity boundary.

## 3. Algorithmic parameters and configuration

| Parameter | Default | Physical / Numerical Interpretation |
| --- | --- | --- |
| `method` | `"positivity"` | Selection of physical constraint. `"positivity"` enforces $f(g) \ge 0$. `"zero_range"` additionally forces density to zero in orientations where $\tilde{f}(g)$ falls below threshold. |
| `zero_range_threshold` | 0.05 m.r.d. | Density cutoff below which orientations are assigned to the zero-range domain. Applicable only when `method="zero_range"`. |
| `odd_regularization` ($\mu$) | $10^{-6}$ | Weight of the minimum-norm Tikhonov term. Balances strict non-negativity against minimal odd energy. |
| `degree_bandlimit` | ODF bandlimit | Maximum series degree for odd harmonics. Restricting odd degrees to match the even bandlimit prevents spurious high-frequency artifacts. |
| `max_iterations` | 500 | Upper bound on L-BFGS-B iterations. Early termination yields a partial correction and is reported in the diagnostic summary. |
| `tolerance` | $10^{-12}$ | Projected gradient termination tolerance for the optimizer. |
| `basis_tolerance` | $10^{-10}$ | Truncation threshold for singular values during Gram matrix orthogonalization. |

## 4. Benchmark performance and numerical verification

The quantitative performance of `correct_ghosts` is verified using synthetic textures
with known analytical distributions (Bunge, 1982; Matthies et al., 1987). For an
orthorhombic specimen with an isolated Gaussian orientation component evaluated at
bandlimit $L=4$:

| Metric | Even-Degree ODF ($f_{\text{even}}$) | Positivity-Corrected ODF ($f_{\text{corr}}$) | Ground Truth ($f_{\text{true}}$) |
| --- | --- | --- | --- |
| Minimum orientation density (m.r.d.) | −0.465 | $\approx 0.000$ | 0.049 |
| Peak orientation density (m.r.d.) | 3.740 | 4.240 | 4.060 |
| Negative volume fraction of $\mathrm{SO}(3)$ | 9.30 % | 0.06 % | 0.00 % |
| Weighted $L_2$ distance to true ODF | 0.190 | 0.092 | 0.000 |
| Ghost amplitude ratio $\lVert f_{\text{odd}} \rVert / \lVert \tilde{f} \rVert$ | 0.000 | 0.145 | 0.151 |
| Residual perturbation on pole figures | — | $< 10^{-3}$ m.r.d. | — |

The reconstructed odd harmonics eliminate negative density artifacts, restore peak
amplitudes toward ground-truth values, and alter forward-projected pole figure
intensities by less than $0.02\%$, consistent with numerical quadrature precision.

## 5. Physical assumptions and failure modes

1. **Series truncation ringing versus ghost phenomena.** When an orientation
   distribution exhibits sharp components that exceed the series truncation limit $L$,
   Gibbs ringing induces negative lobes in both even and odd subspaces. In such cases,
   no combination of odd harmonics of order $\le L$ can satisfy non-negativity, and
   the post-correction infeasibility (`infeasibility_after`) remains elevated.
2. **Symmetry-imposed absence of odd invariants.** As shown in Section 2.2, high crystal
   symmetries preclude odd harmonic invariants at low series degrees. When $L < \ell_{\text{min}}$,
   the ghost correction algorithm correctly identifies an empty odd basis and returns
   the input ODF without modification.
3. **Underdetermined even-degree inversion.** If input pole figures are sparse or
   poorly conditioned, the even-degree coefficients $\tilde{\mathbf{c}}$ reflect the
   inversion regularizer rather than true sample texture. Ghost correction regularizes
   the resulting distribution toward non-negativity, but cannot compensate for
   insufficient experimental projection data.
4. **Epistemic status of the odd component.** The reconstructed odd harmonics
   represent a mathematically regularized inference derived from non-negativity, not
   a direct diffraction measurement. PyTex exposes the corrected ODF as `final_odf`
   while preserving the uncorrected `odf` on the report object to maintain full
   analytical provenance.

## Verification

- `tests/unit/test_ghost_correction.py`: Verifies null-space projection invariance,
  orthonormalization stability, group character dimensions, and non-negativity restoration.
- Executable worked examples:
  - {doc}`../examples/generated/ghost-problem`
  - {doc}`../examples/generated/pole-figure-arithmetic`

## See also

- {doc}`../theory/ghost_problem_and_odd_harmonics` — Comprehensive mathematical derivation of the ghost phenomenon.
- {doc}`../theory/harmonic_odf_reconstruction` — Symmetrized spherical harmonic basis functions and series expansion.
- {doc}`pole_figure_inversion` — Inversion of experimental diffraction pole figures.

## References

### Normative

- Bunge, H. J. (1982). *Texture Analysis in Materials Science: Mathematical Methods*.
  Butterworths. <https://doi.org/10.1016/C2013-0-11769-2>
- Dahms, M. & Bunge, H. J. (1989). The iterative series-expansion method for quantitative
  texture analysis. I. General outline. *Journal of Applied Crystallography* **22**,
  439–447. <https://doi.org/10.1107/S0021889889005261>
- Matthies, S. (1979). On the reproducibility of the orientation distribution function
  of texture samples from pole figures (ghost phenomena). *Physica Status Solidi (b)*
  **92**, K135–K138. <https://doi.org/10.1002/pssb.2220920253>

### Informative

- Hielscher, R. & Schaeben, H. (2008). A novel pole figure inversion method:
  specification of the MTEX algorithm. *Journal of Applied Crystallography* **41**,
  1024–1037. <https://doi.org/10.1107/S0021889808030112>
- Matthies, S., Vinel, G. W. & Helming, K. (1987). *Standard Distributions in Texture
  Analysis*. Akademie-Verlag.
