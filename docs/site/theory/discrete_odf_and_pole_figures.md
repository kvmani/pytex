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

## Plotting And Inspection

PyTex uses the same discrete texture model for computation and plotting.

- Pole-figure contours are rendered from a smoothed density grid in the chosen projection plane.
- Euler-space ODF contours are rendered from the weighted support directly.
- Classical Bunge sections are rendered as kernel-smoothed inspection views through the discrete support.

These plots are intended to be faithful inspection surfaces for the implemented discrete model. They
should not yet be read as a claim of full harmonic ODF reconstruction.

## Current Limits

Full harmonic ODF expansion, broad experimental PF inversion doctrine, and rigorous kernel normalization on $SO(3)$ remain future work. The current representation is deliberately explicit and testable so later algorithms can be layered on top of a stable semantic core.

## Normative References

- H.-J. Bunge, *Texture Analysis in Materials Science: Mathematical Methods*, Butterworths, 1969. DOI: <https://doi.org/10.1016/C2013-0-11769-2>.
- R. A. Fisher, “Dispersion on a Sphere”, *Proceedings of the Royal Society A*, 217(1130), 1953, 295–305. DOI: <https://doi.org/10.1098/rspa.1953.0064>.
