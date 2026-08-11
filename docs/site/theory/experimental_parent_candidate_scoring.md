# Experimental Parent Candidate Scoring

This note records the current experimental PyTex method for ranking candidate parent orientations
against an observed child-orientation set under a fixed orientation relationship.

## Scope

The current experimental implementation supports:

- a fixed `PhaseTransformationRecord`,
- optional explicit variant assignments for each observed child orientation,
- vectorized scoring of multiple candidate parent orientations, and
- mean, median, or max reduction across child residual angles.

## Predicted Child Orientations

For a candidate parent orientation $g_p^{(m)}$ and a child-specific variant rotation
$\Delta g_i$, the predicted child orientation is

$$
\hat{g}_{c,i}^{(m)} = \Delta g_i g_p^{(m)}
$$

If explicit variant assignments are absent, the record's base orientation relationship is reused
for every child.

## Residual Metric

The candidate residual is computed from the child-child misorientation between the observed child
orientation $g_{c,i}$ and the predicted child orientation $\hat{g}_{c,i}^{(m)}$:

$$
\omega_i^{(m)} = \omega\!\left(g_{c,i}, \hat{g}_{c,i}^{(m)}\right)
$$

When symmetry-aware mode is enabled, $\omega$ is reduced to the disorientation angle under the
child symmetry on both sides.

## Candidate Score

PyTex currently reduces the per-child residuals by one of

$$
\operatorname{mean}_i \omega_i^{(m)}, \qquad
\operatorname{median}_i \omega_i^{(m)}, \qquad
\max_i \omega_i^{(m)}
$$

The resulting score is a ranking aid for research workflows, not a stable parent-reconstruction
claim.

## Current Limits

- This method is intentionally experimental and does not claim full parent reconstruction.
- It does not yet include spatial regularization, variant-compatibility priors, or grainwise clustering.
- Literature-backed dataset breadth remains ahead of the current in-repo benchmark surface.
