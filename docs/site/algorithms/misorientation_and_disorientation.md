# Misorientation, Disorientation, And Boundary Statistics

**Surface:** `Orientation.misorientation_to`, `Misorientation.disorientation`,
`OrientationSet.misorientation_angles_to`,
`pytex.core.misorientation_distribution.MisorientationDistribution`,
`random_disorientation_angles_deg`, with
`pytex.ebsd.csl.classify_misorientations` classifying the result and the
workbench operations `ebsd.distribution` and `ebsd.or_from_grains` consuming it.

Everything that distinguishes one grain boundary from another is a statement
about **misorientation**, and almost every such statement is wrong unless it is
made about the **disorientation** instead. This page states the difference, the
reduction that turns one into the other, the baseline a measured distribution
must be compared against, and the three numbers that are routinely quoted for
each other.

## 1. Misorientation is not unique, and that is the whole problem

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

For cubic-cubic that orbit has $24 \times 24 = 576$ members. They are all equally
valid descriptions of the same physical boundary, and their rotation angles range
from a few degrees to nearly $180^\circ$. **Any statement about "the
misorientation angle" that does not say which representative was taken is
undefined**, and a histogram built from arbitrary representatives is noise.

## 2. The disorientation: one representative, chosen canonically

The **disorientation** is the orbit member with the smallest rotation angle — the
one lying in the misorientation fundamental zone.

```text
input : misorientation m, symmetry groups G_left, G_right

1  enumerate every candidate  S_l . m . S_r^T          (576 for cubic-cubic)
2  convert each to a quaternion
3  score each with the canonical fundamental-region key
4  take the minimum key
```

Step 3 is worth dwelling on. Selecting purely on angle is ambiguous, because
distinct orbit members can share the minimum angle exactly — a symmetric
boundary has several equally small representatives. Choosing among them by
whichever floating-point comparison happens to win makes the answer depend on
rounding, and a disorientation *axis* that flips between runs breaks every
downstream statistic that groups by axis. The canonical key breaks ties
deterministically, so the same boundary always returns the same representative.

**Bounds worth memorising.** For cubic-cubic symmetry the disorientation angle
never exceeds $62.8^\circ$. A reported "misorientation" above that is not a
disorientation, and is either an unreduced representative or a different
convention.

### Cost

Naively $O(|G_1| \cdot |G_2|)$ rotations per pair — 576 matrix products for a
cubic pair, and an EBSD map has millions of pairs. The batch surfaces
(`OrientationSet.misorientation_angles_to`, and the pair routines behind
segmentation and KAM) evaluate the orbit vectorised over all pairs at once and
reduce with a scalar projection rather than forming every candidate matrix,
which is what makes a full-map KAM tractable.

## 3. The baseline: Mackenzie, and the three numbers people confuse

A measured distribution means nothing on its own. The reference is the
distribution of a **randomly textured** aggregate — the Mackenzie distribution —
and PyTex generates it by sampling rather than by transcribing a curve:
`random_disorientation_angles_deg` draws Haar-uniform quaternions (uniform on
$SO(3)$, not uniform in Euler angles, which is a different and wrong
distribution) and reduces each to its disorientation.

Three numbers are routinely quoted for one another:

| Quantity | Cubic value | What it is |
| --- | --- | --- |
| **Mode** | $\approx 45^\circ$ | where the distribution peaks |
| **Mean** | $\approx 40.7^\circ$ | the average angle |
| **Maximum** | $62.8^\circ$ | the hard cutoff of the fundamental zone |

The mean is *below* the mode because the distribution is left-skewed: a long
low-angle tail pulls the mean down while the hard cutoff at $62.8^\circ$ stops
the upper side compensating. **Quoting $45^\circ$ as the mean conflates the mode
with the mean**, and `MisorientationDistribution.mean_angle_deg` documents the
distinction at the point of use rather than leaving it to be rediscovered.

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
