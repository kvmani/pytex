# Misorientation, Disorientation, And Boundary Statistics

**Surface:** `Orientation.misorientation_to`, `Misorientation.disorientation`,
`OrientationSet.misorientation_angles_to`,
`pytex.core.misorientation_distribution.MisorientationDistribution`,
`random_disorientation_angles_deg`, with
`pytex.ebsd.csl.classify_misorientations` classifying the result and the
workbench operations `ebsd.distribution` and `ebsd.or_from_grains` consuming it.

Grain boundary properties depend critically on relative crystallographic
orientation. Because crystal symmetry generates multiple equivalent descriptions
of any boundary, quantitative boundary analysis requires the unique canonical
representation known as **disorientation**. This page details the misorientation
symmetry orbit, the algorithm for reducing misorientations into the fundamental
zone, and the theoretical random baseline (Mackenzie distribution) used for
statistical comparison.

## 1. Symmetry equivalence and the misorientation orbit

For two orientations $g_1, g_2$ (crystal-to-specimen, Bunge) the misorientation
is the rotation carrying one crystal frame onto the other:

$$
\Delta g \;=\; g_1^{-1} g_2 .
$$

Each crystal is defined only up to its own symmetry, so $g_i$ and $S g_i$
describe the same crystal for every operator $S$ in the point group. The
misorientation is therefore not one rotation but an **orbit**:

$$
\bigl\{\, S_1 \,\Delta g\, S_2^{\mathsf{T}} \;:\; S_1 \in G_1,\; S_2 \in G_2 \,\bigr\}.
$$

For cubic-cubic bicrystals, that orbit contains $24 \times 24 = 576$ members.
They are all physically equivalent descriptions of the same boundary, with
rotation angles ranging from small values up to $180^\circ$. Consequently,
the scalar misorientation angle is ill-defined unless a specific symmetry
representative is specified. Unreduced or arbitrarily selected orbit
representatives produce inconsistent angular distributions.

## 2. The disorientation: one representative, chosen canonically

The **disorientation** is the orbit member with the smallest rotation angle —
the unique representative lying in the misorientation fundamental zone.

```text
input : misorientation m, symmetry groups G_left, G_right

1  enumerate every candidate  S_l . m . S_r^T          (576 for cubic-cubic)
2  convert each to a quaternion
3  score each with the canonical fundamental-region key
4  take the minimum key
```

Step 3 provides deterministic tie-breaking. Selecting a representative based
solely on rotation angle is insufficient when multiple orbit members attain the
minimum angle simultaneously (common in symmetric tilt and twist boundaries).
PyTex scores candidate quaternions using a canonical fundamental-region key,
guaranteeing that identical boundary misorientations always yield the same
rotation axis and angle across independent evaluations.

**Fundamental bounds.** For cubic-cubic symmetry, the disorientation angle
never exceeds $62.8^\circ$ (the Mackenzie cutoff). A reported misorientation
above this bound represents an unreduced orbit member rather than a disorientation.

### Cost

Naively $O(|G_1| \cdot |G_2|)$ rotations per pair — 576 matrix products for a
cubic pair, and an EBSD map has millions of pairs. The batch surfaces
(`OrientationSet.misorientation_angles_to`, and the pair routines behind
segmentation and KAM) evaluate the orbit vectorised over all pairs at once and
reduce with a scalar projection rather than forming every candidate matrix,
which is what makes a full-map KAM tractable.

## 3. Theoretical random baseline: the Mackenzie distribution

A measured distribution requires comparison against a theoretical reference.
For randomly textured polycrystals, the disorientation angle distribution follows
the Mackenzie distribution. PyTex samples this reference distribution directly:
`random_disorientation_angles_deg` draws Haar-uniform quaternions (uniform on
$SO(3)$, rather than uniform in Euler angles) and reduces each to its disorientation.

Key angular characteristics for cubic–cubic boundaries:

| Quantity | Cubic value | Definition |
| --- | --- | --- |
| **Mode** | $\approx 45^\circ$ | peak of the probability density function |
| **Mean** | $\approx 40.7^\circ$ | expected value across the distribution |
| **Maximum** | $62.8^\circ$ | upper boundary of the cubic misorientation fundamental zone |

The mean angle lies below the mode because the distribution is left-skewed:
a low-angle tail extends to $0^\circ$ while the sharp cutoff at $62.8^\circ$ limits
the upper tail. PyTex documents this distinction explicitly at `MisorientationDistribution.mean_angle_deg`.

### Correlated versus uncorrelated

`MisorientationDistribution` carries a `correlated` flag, and the distinction is
physical:

- **Correlated** — misorientations between *neighbouring* grains. This is the
  boundary population, and it is what twinning and variant selection modify.
- **Uncorrelated** — misorientations between randomly chosen grain pairs. This
  reflects the *texture* alone, with no information about which grains touch.

A correlated distribution departing from its uncorrelated counterpart is
evidence of preferential boundary formation. Comparing a correlated measurement
against the Mackenzie curve instead conflates texture with boundary selection:
a strongly textured material has a non-Mackenzie uncorrelated distribution
before any boundary preference exists at all.

## 4. Reading a distribution

`histogram` bins the angles; the shape is then read against the two baselines
above. The characteristic signals:

| Feature | Usually means |
| --- | --- |
| Sharp spike at $60^\circ$ | $\Sigma 3$ annealing twins in an fcc material |
| Excess below $\approx 15^\circ$ | subgrain structure, or a segmentation threshold set too high |
| Depletion at high angle | strong texture — neighbouring grains are similar |
| Peak matching a transformation OR | variant selection; hand it to {doc}`orientation_relationship_determination` |

## 5. From a distribution to boundary character

The spike at $60^\circ$ is where this page hands over to
{doc}`csl_boundaries`, which classifies individual boundaries against the
coincidence-site-lattice registry rather than reading a histogram.

## Verification

- The random-disorientation baseline against the published Mackenzie mean, in
  {doc}`../examples/generated/random-disorientation`.
- Disorientation reduction and its cubic bound, in
  {doc}`../examples/generated/orientation`.

## See also

- {doc}`../theory/orientation_space_and_disorientation` — the fundamental zone
  and why the reduction is well defined.
- {doc}`../theory/random_disorientation_baseline` — the Mackenzie distribution.
- {doc}`../theory/fundamental_region_reduction` — the canonical key.
- {doc}`csl_boundaries` — classifying a boundary once it is reduced.
- {doc}`ebsd_grains_and_local_misorientation` — where these pairs come from.

## References

### Normative

- Mackenzie, J. K. (1958). Second paper on statistics associated with the
  random disorientation of cubes. *Biometrika* **45**, 229-240.
  <https://doi.org/10.1093/biomet/45.1-2.229>
- Morawiec, A. (2004). *Orientations and Rotations: Computations in
  Crystallographic Textures*. Springer.
  <https://doi.org/10.1007/978-3-662-09156-2>

### Informative

- Randle, V. (2004). Twinning-related grain boundary engineering. *Acta
  Materialia* **52**, 4067-4081.
  <https://doi.org/10.1016/j.actamat.2004.05.031>
