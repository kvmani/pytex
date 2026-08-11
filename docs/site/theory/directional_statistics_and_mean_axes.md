# Directional Statistics: Mean Axes And The Orientation Tensor

Averaging directions looks like averaging vectors and is not. A crystal direction $[uvw]$ and its
negative $[\bar{u}\bar{v}\bar{w}]$ are the *same physical axis*, and for data of that kind the
arithmetic mean is not merely inaccurate — it is identically zero. This note gives the estimator
that works, the eigenvalue test that says whether the answer means anything, and the closed forms
that let both be checked.

It covers `SphericalVectorSet.orientation_tensor` and `SphericalVectorSet.mean_direction`.

## Why The Vector Mean Fails

If a set of axes is closed under negation — every $\hat{\mathbf{v}}$ appearing with
$-\hat{\mathbf{v}}$ — its resultant is exactly zero:

$$
\sum_{i} \hat{\mathbf{v}}_{i} = \mathbf{0} .
$$ (eq-ds-cancel)

Numerically confirmed at $5.2\times10^{-14}$ for 5000 axis pairs. This is not an edge case that
careful data avoids: measured axes carry an arbitrary sign, so any real axial dataset is *randomly*
signed, and the cancellation is partial and unpredictable rather than complete.

The failure is worth seeing concretely. Take 3000 axes tightly clustered about $\hat{\mathbf{z}}$ —
a strong, unambiguous texture — and assign each an independent random sign:

| Estimator | Result | True axis |
| --- | --- | --- |
| normalized resultant $\sum_i \hat{\mathbf{v}}_i / \lVert\cdot\rVert$ | $(-0.244,\ -0.157,\ -0.957)$ | $(0,0,1)$ |
| principal eigenvector of $\boldsymbol{\Theta}$ | $(0.002,\ -0.006,\ 1.000)$ | $(0,0,1)$ |

The naive answer is wrong in sign and tilted several degrees off axis, and it would change if the
signs were redrawn. Nothing about the data is ambiguous; the estimator is simply inapplicable.

## The Orientation Tensor

The fix is to use a statistic that cannot see the sign. The outer product does not:
$(-\hat{\mathbf{v}})(-\hat{\mathbf{v}})^{\mathsf{T}} = \hat{\mathbf{v}}\hat{\mathbf{v}}^{\mathsf{T}}$.
So the second moment

$$
\boldsymbol{\Theta} = \frac{1}{n}\sum_{i} \hat{\mathbf{v}}_{i}\hat{\mathbf{v}}_{i}^{\mathsf{T}}
$$ (eq-ds-tensor)

is the lowest-order moment that survives antipodal identification, and it is the natural summary of
axial data. It is symmetric positive semi-definite with

$$
\operatorname{tr}\boldsymbol{\Theta} = \frac{1}{n}\sum_{i}\lVert\hat{\mathbf{v}}_{i}\rVert^{2} = 1,
\qquad
\lambda_{1} + \lambda_{2} + \lambda_{3} = 1,
\qquad
\lambda_{1} \le \lambda_{2} \le \lambda_{3} .
$$ (eq-ds-trace)

The unit trace is automatic for unit vectors and is a free check on any implementation.

## The Eigenvalues Classify The Distribution

Because the eigenvalues are non-negative and sum to one, they live on a triangle whose corners are
the three limiting distributions — and each corner is an exact closed form:

| Distribution | $(\lambda_{1}, \lambda_{2}, \lambda_{3})$ | Computed |
| --- | --- | --- |
| uniform over the sphere | $(\tfrac{1}{3}, \tfrac{1}{3}, \tfrac{1}{3})$ | $(0.33196, 0.33258, 0.33546)$, $n=2\times10^{5}$ |
| perfect girdle (all in one plane) | $(0, \tfrac{1}{2}, \tfrac{1}{2})$ | $(0, 0.498037, 0.501963)$ |
| perfect cluster (all parallel) | $(0, 0, 1)$ | $(0, 0, 1)$ exactly |

The uniform case is exact in expectation: for directions uniform on the sphere,
$\mathbb{E}[\hat{v}_{i}\hat{v}_{j}] = \delta_{ij}/3$, so
$\mathbb{E}[\boldsymbol{\Theta}] = \mathbf{I}/3$. The girdle case follows from
$\langle\cos^{2}\rangle = \langle\sin^{2}\rangle = 1/2$ around a circle, and the cluster case is
immediate.

This is the practical value of {eq}`eq-ds-tensor`: **three numbers say what kind of texture the
data has**, without contouring anything. $\lambda_{3}$ near 1 is a fibre; $\lambda_{1}$ near 0 with
$\lambda_{2}\approx\lambda_{3}$ is a girdle — poles spread around a great circle, which is what a
rolling texture's plane normals do; all three near $1/3$ is no texture at all.

## The Mean Direction, And When It Does Not Exist

For antipodal data the mean axis is the eigenvector of $\lambda_{3}$, sign-canonicalized into the
upper hemisphere so that a rerun cannot return the opposite vector — a determinism requirement, not
a physical statement, since the sign was never meaningful. For non-antipodal data, where sign *is*
data, the estimator is the normalized resultant, and PyTex raises rather than returning a direction
when that resultant is numerically zero.

**The eigenvalue gap decides whether the answer is meaningful, and it must be checked separately.**
The principal eigenvector is identified only when $\lambda_{3} > \lambda_{2}$. As
$\lambda_{3} \to \lambda_{2}$ the principal direction becomes arbitrary within the corresponding
eigenplane, and for uniform data — where all three eigenvalues coincide — it is arbitrary
altogether. The estimator still returns a perfectly good unit vector: in the $\pm$-symmetric test
above it returned $(0.813, -0.569, 0.126)$ with norm exactly 1, a confidently reported direction
carrying no information whatever.

Nothing in the return type reveals this. A mean direction should therefore be quoted with its
eigenvalues, and a girdle distribution — $\lambda_{2}\approx\lambda_{3}$ — has no mean axis at all,
only a mean *plane*, whose normal is the eigenvector of $\lambda_{1}$ and is well determined
precisely when the axis is not.

## Assumptions And Limits

- $\boldsymbol{\Theta}$ is a second moment and sees only the ellipsoidal part of a distribution. It
  cannot distinguish a single cluster from two antipodal clusters, nor a girdle from a
  four-fold-symmetric arrangement in the same plane.
- The tensor is unweighted here: every direction counts once. Weighting by solid angle is required
  when the directions come from a raster rather than from a population — see
  {doc}`pole_figure_arithmetic_and_mrd`, where the same neglect costs exactly 50%.
- Eigenvalue uncertainty is not estimated. With small $n$ the gap $\lambda_{3}-\lambda_{2}$ can be
  an artefact of sampling, and no test for that is provided.

## References

### Informative

- N. I. Fisher, T. Lewis and B. J. J. Embleton, *Statistical Analysis of Spherical Data*, Cambridge
  University Press (1987). The orientation tensor, axial data, and why the resultant fails.
- N. H. Woodcock, *Specification of fabric shapes using an eigenvalue method*, Geological Society of
  America Bulletin **88** (1977) 1231–1236.
  DOI: <https://doi.org/10.1130/0016-7606(1977)88<1231:SOFSUA>2.0.CO;2>. The eigenvalue
  classification of cluster and girdle fabrics.
- K. V. Mardia and P. E. Jupp, *Directional Statistics*, Wiley (2000).

## See Also

- {doc}`pole_figure_arithmetic_and_mrd` — solid-angle weighting, needed whenever these directions
  come from a measured raster.
- {doc}`ipf_color_keys` — the other place where antipodal identification changes the answer.
