# Reciprocal Space And Kinematic Spots

PyTex now includes an explicit reciprocal-space and zone-axis layer for minimal diffraction spot prediction.

## Reciprocal Vectors

Given Miller indices $(h,k,l)$, PyTex constructs a reciprocal-lattice vector

$$
\mathbf{g} = h \mathbf{a}^\ast + k \mathbf{b}^\ast + l \mathbf{c}^\ast
$$

## Zone Axis Condition

For a zone axis $\mathbf{z}$ in direct space, a reciprocal vector belongs to the zone when

$$
\mathbf{z} \cdot \mathbf{g} = 0
$$

The current implementation uses this relation as an explicit filter when a zone axis is provided.

## Kinematic Spot Construction

With incident wavevector $\mathbf{k}_{\mathrm{in}}$ and reciprocal vector $\mathbf{g}$ in the laboratory frame, the candidate outgoing wavevector is

$$
\mathbf{k}_{\mathrm{out}} = \mathbf{k}_{\mathrm{in}} + \mathbf{g}
$$

PyTex then evaluates the excitation error as the scalar mismatch

$$
s = \lVert \mathbf{k}_{\mathrm{out}} \rVert - \lVert \mathbf{k}_{\mathrm{in}} \rVert
$$

and accepts candidates within a configured tolerance.

## Reflection Families

PyTex groups reflections into symmetry-aware families in crystal reciprocal space. Let $\mathbf{g}_c$ be the reciprocal vector in the crystal basis and let $\mathcal{S}$ denote the proper crystal symmetry group. The family direction key is constructed from a canonicalized representative

$$
\widehat{\mathbf{g}}_{\mathrm{fam}} =
  \operatorname{canon}_{\mathcal{S}}\!\left(
    \frac{\mathbf{g}_c}{\lVert \mathbf{g}_c \rVert}
  \right)
$$

and the family key retains the reciprocal magnitude $\lVert \mathbf{g}_c \rVert$ as a separate invariant so that $(100)$ and $(200)$ do not collapse into the same family.

## Detector Projection

Accepted outgoing directions are normalized and intersected with the detector plane defined by the camera length, detector basis, and pattern center.

## Detector Acceptance Masks

After detector projection, PyTex may apply an explicit acceptance mask. The current implementation supports a rectangular inset and an optional radial bound about the pattern center. This separates detector-plane intersection, detector containment, and workflow-specific acceptance into distinct states.

## Proxy Intensity Weighting

PyTex currently provides a minimal ranking intensity rather than a physical structure-factor model. In the default proxy mode, the intensity is

$$
I_{\mathrm{proxy}} =
  \frac{1}{1 + \left(s / \sigma_s\right)^2}
  \cdot
  \frac{1}{1 + \lVert \mathbf{g} \rVert^2}
$$

where $s$ is the excitation error, $\sigma_s$ is the configured excitation scale, and $\mathbf{g}$ is the reciprocal vector in the laboratory frame. The first factor penalizes off-Ewald candidates and the second suppresses higher-order reflections in a simple teaching-grade way.

## Double Diffraction And Forbidden Reflections

A kinematic calculation reports a reflection as absent whenever its structure factor vanishes. Recorded zone-axis patterns routinely show such reflections anyway, because a beam diffracted by $\mathbf{g}_1$ is itself an incident beam inside the crystal and can diffract again by $\mathbf{g}_2$, leaving the specimen along

$$
\mathbf{g} = \mathbf{g}_1 + \mathbf{g}_2
$$

The set of reflections reachable this way is therefore the set of pairwise algebraic sums of the excited reflections. This selection rule is purely geometric: it follows from the additivity of scattering vectors and is independent of the multi-beam dynamical theory that determines how much intensity each path actually carries. PyTex implements the rule exactly and the intensity only as an estimate.

Two structural consequences follow, and both are checked by tests.

First, *centring absences can never be revived*. The reflections allowed by a centring condition form a sublattice $\Lambda$ of the reciprocal lattice, and a sublattice is closed under addition, so $\mathbf{g}_1, \mathbf{g}_2 \in \Lambda$ implies $\mathbf{g}_1 + \mathbf{g}_2 \in \Lambda$. Only absences imposed by the *basis* — a glide plane, a screw axis, or the arrangement of the motif — can be produced by double diffraction. This matches the experimental record: the forbidden $002$ of diamond-cubic silicon appears along $[110]$, produced by $(111) + (\bar{1}\bar{1}1)$, while an $F$-centred $100$ never does.

Second, the added reflections are not kinematic observations and must not be reported as if they were. PyTex assigns them

$$
I_{\mathrm{dd}}(\mathbf{g}) = c \sum_{\mathbf{g}_1 + \mathbf{g}_2 = \mathbf{g}} I(\mathbf{g}_1) \, I(\mathbf{g}_2)
$$

where the product follows from the two-step amplitude scaling as $F(\mathbf{g}_1) F(\mathbf{g}_2)$, the sum runs over unordered pairs of excited reflections, and the coupling constant $c$ absorbs everything a kinematic treatment cannot supply — beam coupling strength and specimen thickness. Path phases are unknown at this level of theory, which is why contributions are summed in intensity rather than in amplitude. The result is an observability estimate. Every such reflection is flagged in the spot table, named with the path that produced it, and rendered distinctly.

## Detector-Space Indexing Association

Observed detector coordinates can be clustered in detector space and then associated with simulated spots. Let $\mathbf{p}_i$ denote a detector-space observation cluster center and let $\mathbf{s}_j$ denote a simulated detector coordinate. PyTex currently evaluates the detector residual

$$
r_{ij} = \lVert \mathbf{p}_i - \mathbf{s}_j \rVert
$$

accepts matches under a configured residual threshold, and summarizes the result as an indexing candidate with match fraction, mean residual, and a simple aggregate score.

## Local Candidate Refinement

PyTex now also provides a deterministic local refinement loop around a seed orientation. The current implementation parameterizes the neighborhood in Bunge Euler space, evaluates a local Cartesian grid, and retains the best-scoring candidate after each shrink step. This is intentionally a transparent local search rather than a hidden optimizer.

## Family-Level Indexing Reports

Because spot-by-spot matches can be hard to interpret, PyTex also aggregates matches to reflection families. The current family report records the representative Miller indices, multiplicity, represented simulated-spot count, matched count, matched fraction, total family intensity, and mean detector residual.

## Current Limits

- This is a minimal teaching-grade and geometry-grade kinematic spot workflow.
- Excitation error is currently a simple magnitude mismatch, not a full refinement metric.
- Proxy intensity weighting is not a substitute for structure factors, polarization terms, or dynamical scattering.
- Local refinement now exists, but continuous gradient-style or probabilistic refinement remains future work.
- Double diffraction supplies the correct *set* of extra reflections, but their intensities are a coupling estimate, not a solved multi-beam calculation.
