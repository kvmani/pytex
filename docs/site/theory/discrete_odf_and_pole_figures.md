# Discrete ODF and Pole-Figure Foundations in PyTex

## Scope

The current PyTex texture layer implements the foundational algorithms needed before full inversion and harmonic reconstruction:

- pole-figure synthesis from orientation sets,
- inverse-pole-figure synthesis from orientation sets,
- kernel-weighted ODF evaluation,
- volume-fraction queries around a reference orientation,
- discrete pole-figure inversion over an explicit orientation dictionary,
- contour pole-figure plotting from projected density grids,
- kernel-smoothed Bunge-section plotting for discrete ODF inspection.
- non-negative fitting of named texture-component mixtures to normalized ODF density.

## Pole-Figure Construction

Given a pole normal in the crystal frame and an orientation set $\{g_i\}$, the corresponding specimen directions are obtained by applying each orientation to the pole family. The current implementation uses explicit crystal-symmetry expansion of the pole family and stores the resulting specimen directions together with their weights.

## ODF Evaluation

PyTex currently treats the ODF as a weighted orientation set together with an angular kernel. This is a foundational representation rather than a final inversion framework. Two kernels are presently supported:

- a halfwidth-calibrated cosine-power kernel exposed as `de_la_vallee_poussin`,
- a von Mises–Fisher style kernel.

The current ODF representation is therefore a discrete support together with non-negative weights.
This supports explicit estimation, explicit inversion, and explicit plotting without pretending that
the present implementation is already a harmonic expansion on $SO(3)$.

## Current Limits

## Discrete Pole-Figure Inversion

The current inversion path is deliberately explicit. Given measurement directions $\mathbf{s}_m$, an orientation dictionary $\{g_j\}$, a pole family $\mathcal{H}$, and an angular kernel $K$, PyTex builds

$$
A_{mj} = \frac{1}{|\mathcal{H}|}\sum_{h \in \mathcal{H}} K\!\left(\angle(\mathbf{s}_m, g_j h)\right)
$$

then solves a regularized nonnegative least-squares problem

$$
\min_{\mathbf{w} \ge 0}\ \frac{1}{2}\lVert A\mathbf{w} - \mathbf{b}\rVert_2^2 + \frac{\lambda}{2}\lVert \mathbf{w}\rVert_2^2
$$

followed by normalization of the recovered weights so they define a discrete ODF.

This fits the present PyTex architecture because the orientation support remains explicit and scientifically inspectable.

## Named-Component Mixture Fitting

Given normalized ODF density $f(g_i)$ on explicit evaluation orientations and normalized kernel
responses $K(g_i,c_j)$ around named ideal components $c_j$, PyTex solves

$$
\min_{a_j,a_r\geq 0}\left\|\sum_j a_j K(g_i,c_j)+a_r-f(g_i)\right\|_2^2,
\qquad \sum_j a_j+a_r=1.
$$

The constant $a_r$ column is the optional random-texture contribution because a random ODF has
density one multiple of random everywhere. `fit_odf_components(...)` uses SciPy's constrained
SLSQP solver, checks that the design has full column rank, and returns `ODFComponentFit` with the
fractions, observed and predicted densities, RMS/maximum residuals, $R^2$, `describe()`, and a JSON
contract. The default component peaks use the ODF kernel; a different declared `KernelSpec` can be
supplied.

This is a basis fit, not an automatic phase-discovery method. Fractions depend on the named
components offered, the kernel halfwidth, and the orientation sampling used for the unweighted
residual. The default evaluates on the ODF support; a sparse or strongly non-uniform support should
be replaced with an explicit approximately uniform `evaluation_orientations` grid. A deficient
support raises instead of returning arbitrary fractions. The random term absorbs only uniform
density, not an omitted broad component, ghost correction, or uncertainty.

## Plotting And Inspection

PyTex uses the same discrete texture model for computation and plotting.

- Pole-figure contours are rendered from a smoothed density grid in the chosen projection plane.
- Euler-space ODF contours are rendered from the weighted support directly.
- Classical Bunge sections are rendered as kernel-smoothed inspection views through the discrete support.

These plots are intended to be faithful inspection surfaces for the implemented discrete model. They
should not yet be read as a claim of full harmonic ODF reconstruction.

## Current Limits

Bootstrap uncertainty and automatic component-centre/halfwidth refinement remain future work. The
current named-component surface fits declared ideal centres and keeps its residual explicit.

## Normative References

- H.-J. Bunge, *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths, 1969. DOI: <https://doi.org/10.1016/C2013-0-11769-2>.
- R. A. Fisher, “Dispersion on a Sphere”, *Proceedings of the Royal Society A*, 217(1130), 1953, 295–305. DOI: <https://doi.org/10.1098/rspa.1953.0064>.
