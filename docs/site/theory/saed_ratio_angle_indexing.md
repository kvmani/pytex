# Ratio/Angle Indexing Of A Measured SAED Pattern

The kinematic-spot note fixes the forward problem: given a phase and a zone axis, where do the
spots fall. This note fixes the inverse problem as implemented in
`pytex.diffraction.solving.solve_saed_pattern`: given picked spot positions and enough
calibration to scale them, determine the phase, the zone axis, the crystal orientation in the
pattern frame, and the indices of every spot.

## Scope, And Why Two Indexing Surfaces Exist

PyTex carries two indexing paths that consume different evidence and must not be conflated.
`pytex.diffraction.models.index_saed_pattern` starts from a calibrated
`DiffractionGeometry` — detector distance, pattern centre, tilt, pixel pitch — and works
in detector pixels; it suits an automated pipeline attached to an instrument. The surface treated
here starts from a bare list of spot coordinates and a camera constant, which is what a reader of
a printed micrograph has. Neither subsumes the other.

## Calibration

A picked position becomes a reciprocal-space length through the camera constant $L\lambda$,

$$
\lVert \mathbf{g} \rVert = \frac{r}{L\lambda},
\qquad
d = \frac{1}{\lVert \mathbf{g} \rVert}
$$ (eq-saed-camera-constant)

with $r$ the distance from the transmitted beam. Pixel coordinates are scaled by the pixel pitch
before {eq}`eq-saed-camera-constant`; the camera constant is either supplied directly or formed
as the product of a camera length and the relativistic electron wavelength.

The transmitted beam is the calibration centre and not a listed spot; a spot coinciding with it
carries no direction and is rejected. Coordinates in pixels or millimetres without a camera
constant are rejected at construction rather than at first use, an uncalibrated length being
unrecoverable.

## Admissibility: The Ratio And Angle Tests

Two non-collinear reflections determine a zone, so the two shortest non-collinear measured
vectors seed the solution; shortest because they are the best determined relative to the picking
error, which is approximately constant in absolute terms.

For a candidate phase, all reflections permitted by its lattice centring are enumerated to a
bound. A calculated pair $\left(\mathbf{g}^{c}_{1}, \mathbf{g}^{c}_{2}\right)$ is admissible for
an observed pair $\left(\mathbf{g}^{o}_{1}, \mathbf{g}^{o}_{2}\right)$ when

$$
\frac{\bigl| \lVert \mathbf{g}^{c}_{i} \rVert
- \lVert \mathbf{g}^{o}_{i} \rVert \bigr|}{\lVert \mathbf{g}^{o}_{i} \rVert}
\ \le\ \varepsilon_{\ell}
\quad (i = 1,2),
\qquad
\bigl| \theta^{c} - \theta^{o} \bigr| \ \le\ \varepsilon_{\theta}
$$ (eq-saed-admissible)

where $\theta$ denotes the interplanar angle within the pair. Condition {eq}`eq-saed-admissible`
is the classical ratio-and-angle test of single-crystal pattern indexing; it is evaluated as one
pairwise cosine matrix over the two admissible reflection pools rather than as a nested loop.

The default tolerances are calibrated rather than arbitrary. A relative length tolerance of
$0.03$ covers the centring error of a hand-picked spot at typical camera constants while still
separating the $\{111\}$ and $\{200\}$ reflections of a face-centred cubic metal, whose lengths
differ by approximately fifteen percent. The default angular tolerance is $2^{\circ}$.

**Intensities are excluded by design.** A kinematic intensity model is not reliable enough
to index against, and a printed pattern seldom carries calibrated intensities. Geometry alone
decides an assignment; intensity is carried through for plotting and record-keeping only.

## Zone Axis And Orientation

The zone axis follows from

$$
\mathbf{z} \ \propto\ \mathbf{g}^{c}_{1} \times \mathbf{g}^{c}_{2}
$$

converted to direct-lattice components through the metric tensor and rationalized to the nearest
primitive integer triple.

The crystal-to-pattern rotation is obtained by matching two right-handed orthonormal triads. Let

$$
\hat{\mathbf{e}}_{1} = \frac{\mathbf{g}^{c}_{1}}{\lVert \mathbf{g}^{c}_{1} \rVert},
\qquad
\hat{\mathbf{e}}_{2} = \frac{\mathbf{g}^{c}_{2}
- \left( \mathbf{g}^{c}_{2} \cdot \hat{\mathbf{e}}_{1} \right) \hat{\mathbf{e}}_{1}}
{\bigl\lVert \mathbf{g}^{c}_{2}
- \left( \mathbf{g}^{c}_{2} \cdot \hat{\mathbf{e}}_{1} \right) \hat{\mathbf{e}}_{1} \bigr\rVert},
\qquad
\mathbf{E} = \bigl[\, \hat{\mathbf{e}}_{1} \ \ \hat{\mathbf{e}}_{2} \ \
\hat{\mathbf{e}}_{1} \times \hat{\mathbf{e}}_{2} \,\bigr]
$$

and let $\mathbf{F}$ be the same construction applied to the observed pair embedded in the
detector plane. Then

$$
\mathbf{R} = \mathbf{F}\,\mathbf{E}^{\mathsf{T}}
$$ (eq-saed-rotation)

carries crystal Cartesian vectors into the pattern frame. Equation {eq}`eq-saed-rotation` is a
proper rotation by construction, as a product of two orthonormal right-handed bases, so no
re-orthogonalization or determinant repair is required.

## Multi-Spot Assignment And Residual Error Minimization

Once the seed rotation $\mathbf{R}$ is formed from two non-collinear spots, all permitted reflections
$\mathbf{g}_j^{c}$ for the candidate phase are projected into the detector frame:

$$
\mathbf{g}_j^{\text{det}} = \mathbf{R}\,\mathbf{g}_j^{c}
$$

Reflections lying in the zero-order Laue zone satisfy $\lvert (\mathbf{g}_j^{\text{det}})_z \rvert \le \varepsilon_{\ell} \max_i \lVert \mathbf{g}_i^o \rVert$.
Their in-plane components $\mathbf{g}_j^{\text{proj}} = \left( (\mathbf{g}_j^{\text{det}})_x, (\mathbf{g}_j^{\text{det}})_y \right)$
form the pool of predicted spot positions.

Each measured spot $\mathbf{g}_i^o$ ($i = 1, \dots, N$) searches for the nearest unclaimed prediction within
its adaptive match radius:

$$
r_{\text{match}, i} = \varepsilon_{\ell}\,\lVert \mathbf{g}_i^o \rVert
$$

A predicted reflection $\mathbf{g}_j^{\text{proj}}$ is assigned to spot $i$ if:

$$
\delta_i \equiv \lVert \mathbf{g}_i^o - \mathbf{g}_j^{\text{proj}} \rVert \le r_{\text{match}, i}
$$ (eq-saed-residual)

and $\mathbf{g}_j^{\text{proj}}$ has not been claimed by any closer measured spot. Reflections are not reused,
enforcing a strictly bijective (one-to-one) correspondence between physical reflections and observed spots.

The residual error across all $N_{\text{indexed}}$ assigned spots is quantified by the mean reciprocal residual:

$$
\bar{\delta} = \frac{1}{N_{\text{indexed}}} \sum_{i \in \text{indexed}} \delta_i
$$ (eq-saed-mean-residual)

and the matched fraction:

$$
f_{\text{match}} = \frac{N_{\text{indexed}}}{N}
$$ (eq-saed-matched-fraction)

When users pick four or more spots on the interactive plate, the planar lattice is over-determined.
PyTex provides linear least-squares refinement of the beam centre $\mathbf{c}$ and lattice basis vectors
$\mathbf{a}, \mathbf{b}$ prior to indexing (`pytex.diffraction.lattice_fit`), eliminating manual
beam-picking bias and isolating outlier picks. See [Fitting The Pattern Lattice And Scoring The Solutions](lattice_fit_and_solution_scoring.md)
for the full least-squares derivation.

## Ranking And Deduplication

Solutions are ordered by matched fraction first and mean residual second. The ordering is
deliberate: a solution accounting for every spot with moderate residuals is preferable to one
accounting for half of them exactly, the latter usually indicating a coincidence on a sub-lattice.

Many seed assignments are related by a crystal symmetry operation and represent one physical
answer expressed through different bookkeeping. Two solutions are identified when

$$
\mathbf{R}_{1} \simeq \mathbf{R}_{2}\,\mathbf{S},
\qquad \mathbf{S} \in G
$$

and the retained representative is rewritten into the conventional description — fewest negative
indices, then lowest — so that a cubic pattern down a cube axis is reported as $[001]$ rather
than the equally valid $[0\bar{1}0]$ that the seed search happened to encounter first.

## The Intrinsic Zone-Sense Ambiguity

A single selected-area pattern cannot distinguish a zone axis from its reverse when the
reflection set is centrosymmetric. Inversion of the crystal maps $\mathbf{g} \mapsto -\mathbf{g}$
and, by Friedel's law in the absence of anomalous scattering, leaves both the positions and the
kinematic intensities unchanged. The two senses correspond to genuinely distinct proper rotations
that index the pattern equally well.

This is a property of the experiment and not a deficiency of the algorithm, so it is reported as
such rather than resolved by an arbitrary choice. Resolution requires additional evidence: a
second zone axis obtained by tilting, or a dynamical observation such as convergent-beam
diffraction.

## Assumptions And Limits

- The spots are assumed to lie in one zero-order Laue zone about a low-index axis. A crystal tilted off zone — for instance a transformation variant viewed along a *parent* zone axis, whose own child zone axis is generally irrational — is only partly indexed. The partial match is the correct outcome; a complete match would indicate that reflections had been invented.
- Systematic absences are taken from each candidate phase's space group, so a phase supplied without one is treated as primitive and may be offered reflections its true symmetry forbids.
- Higher-order Laue zone rings and double-diffraction spots are outside the model.
- Spot detection from image data is outside the scope; the algorithm consumes picked or listed coordinates.

## Normative references

International Tables for Crystallography, Vol. A (reflection conditions and metric relations).

## Informative references

Edington, J. W., *Practical Electron Microscopy in Materials Science*, Monograph 2:
Electron Diffraction in the Electron Microscope.
Williams, D. B., Carter, C. B., *Transmission Electron Microscopy*, 2nd ed., Springer, 2009.
De Graef, M., *Introduction to Conventional Transmission Electron Microscopy*, Cambridge
University Press, 2003.
