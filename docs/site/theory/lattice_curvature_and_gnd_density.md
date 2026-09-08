# Lattice Curvature And GND Density

This note fixes how PyTex converts a measured orientation map into lattice
curvature and geometrically necessary dislocation density, and states precisely
what a two-dimensional surface measurement can and cannot determine. It is the
theory behind `pytex.ebsd.gnd`:
`lattice_curvature_tensor`,
`nye_dislocation_density_tensor`, and
`geometrically_necessary_dislocation_density`.

## Why Curvature Carries Dislocation Content

Dislocations whose Burgers vectors cancel over a region — statistically stored
dislocations — leave the lattice orientation unchanged on average. Dislocations
that do *not* cancel leave a net Burgers vector, and the only way a
continuous lattice can accommodate that is by bending. The bending is therefore a
direct measurement of the non-cancelling, or geometrically necessary,
dislocation content.

This is what makes the quantity useful: it is measurable from orientation data
alone, in physical units, and it is the part of the dislocation population that
governs long-range internal stress and the stored energy driving
recrystallization.

## Curvature And The Nye Tensor

For a lattice rotation field with rotation vector
$\boldsymbol{\omega}(\mathbf{x})$, the lattice curvature tensor is

$$
\kappa_{ij} = \frac{\partial \omega_i}{\partial x_j}
$$

and Nye's dislocation density tensor follows as

$$
\alpha_{ij} = \kappa_{ji} - \delta_{ij}\,\kappa_{kk}
$$

For the small elastic strains of a deformed but unloaded specimen the curvature
is dominated by lattice rotation, which is exactly what an orientation map
records.

PyTex forms $\boldsymbol{\omega}$ between neighbouring points in the
*specimen* frame, as the rotation vector of
$\mathbf{g}_{2}\mathbf{g}_{1}^{\mathsf{T}}$, and differentiates by central
differences along the two map axes, falling back to one-sided differences at the
map edges.

**Symmetry is deliberately not reduced.**

Between neighbouring points of a continuously bent lattice the true relative
rotation is small. Applying a disorientation reduction would replace it with a
symmetry-equivalent representative whose *axis* is unrelated to the
physical curvature, destroying the directional information the Nye tensor is
built from. The consequence is that points genuinely straddling a grain boundary
produce large, meaningless curvature, and must be masked by a grain segmentation
before the map is interpreted. This is a real limitation, and it is stated
rather than hidden.

## What A Surface Map Cannot Determine

A two-dimensional map supplies gradients along two directions only, so the third
column of $\kappa$ is unmeasured. Six of the nine components of $\kappa$ are
therefore available, and hence five of the nine components of $\alpha$:

- $\alpha_{01}, \alpha_{02}, \alpha_{10}, \alpha_{12}$ are fully determined, since off the diagonal $\alpha_{ij} = \kappa_{ji}$ and those require only the first two columns of $\kappa$;
- $\alpha_{20}$ and $\alpha_{21}$ require the unmeasured third column and are not determined;
- every diagonal component contains the trace, which carries the unmeasurable $\kappa_{22}$; only the *difference* $\alpha_{00} - \alpha_{11} = \kappa_{00} - \kappa_{11}$ is determined.

PyTex reports the unmeasured entries as `NaN` rather than as zero, so
“measured to be zero” and “not measured” remain distinguishable. The
implementation subtracts the trace on the diagonal directly rather than as
$\mathrm{tr}(\kappa)\,\mathbf{I}$, because in IEEE arithmetic
$\mathrm{NaN} \times 0$ is $\mathrm{NaN}$ and the identity form would propagate
the unmeasured trace into every off-diagonal component — destroying exactly the
information the tensor exists to carry. This was a real defect, caught by a
planted-gradient test.

## Scalar Density

The default estimate sums the absolute values of the four fully determined
components and divides by the Burgers vector magnitude,

$$
\rho_{\mathrm{GND}} = \frac{1}{b}
    \bigl( |\alpha_{01}| + |\alpha_{02}| + |\alpha_{10}| + |\alpha_{12}| \bigr)
$$

PyTex also offers the widely quoted scalar estimate

$$
\rho_{\mathrm{GND}} = \frac{2\theta}{b\,u}
$$

where $\theta$ denotes the kernel average misorientation in radians, $u$ is the spatial step
size, and $b$ is the Burgers vector magnitude. While this scalar formulation omits gradient
directionality, it is standard in the metallurgical literature for comparative characterization.
For an idealized single-axis tilt boundary, the scalar model and the curvature tensor formulation
yield identical results: on a four-connected grid, the interior KAM equals half the per-step
misorientation, so $2\theta/(bu)$ reduces to $(\partial\theta/\partial x)/b$. PyTex verifies this
exact equivalence in unit tests.

## Methodological Assumptions and Physical Constraints

**Lower-bound property.** Dislocation dipoles and multipoles that produce no net lattice
curvature do not induce spatial orientation gradients, and statistically stored dislocations (SSDs)
cancel over short correlation lengths. Consequently, $\rho_{\mathrm{GND}}$ represents a rigorous
lower bound on the total dislocation density ($\rho_{\text{total}} = \rho_{\mathrm{GND}} + \rho_{\mathrm{SSD}}$).

**Resolution dependence.** Dislocation density derived from orientation gradient mapping is not an
intrinsic material constant independent of measurement scale. Finer acquisition step sizes resolve
higher local curvature gradients, systematically increasing computed $\rho_{\mathrm{GND}}$.
Quantitative comparisons across specimens are valid only when acquired at identical spatial step
sizes and evaluated with identical kernel topologies. PyTex requires explicit declaration of the
step scale argument to ensure transparency.

## Units

Curvature is returned in radians per metre and density in reciprocal square
metres. Map coordinates carry no declared unit, so the conversion is an explicit
argument, `step_scale_m`, defaulting to $10^{-6}$ — micrometres, the
EBSD convention. Burgers vectors are given in nanometres. A wrong unit
declaration rescales every result, which is why it is a named parameter rather
than an assumption.

## References

Nye, J. F., *Acta Metall.* **1**, 153–162 (1953),
`DOI: 10.1016/0001-6160(53)90054-6`.

Pantleon, W., *Scripta Mater.* **58**, 994–997 (2008),
`DOI: 10.1016/j.scriptamat.2008.01.050`.

Kysar, J. W., Saito, Y., Oztop, M. S., Lee, D. and Huh, W. T.,
*Int.\ J.\ Plasticity* **26**, 1097–1123 (2010),
`DOI: 10.1016/j.ijplas.2010.03.009`.

Calcagnotto, M., Ponge, D., Demir, E. and Raabe, D.,
*Mater.\ Sci.\ Eng.\ A* **527**, 2738–2746 (2010),
`DOI: 10.1016/j.msea.2010.01.004`.
