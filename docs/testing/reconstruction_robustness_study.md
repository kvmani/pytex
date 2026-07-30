# Parent-Grain Reconstruction: Robustness Study

Operating envelope for `pytex.experimental.reconstruct_parent_grains`, measured against planted
ground truth. This is the evidence base for the reconstruction row of the
[Phase Transformation Validation Matrix](phase_transformation_validation_matrix.md), and it is the
characterization the surface must have before it can leave `experimental`.

Regenerate with:

```bash
python scripts/study_reconstruction_robustness.py
```

## What was measured

Two sweeps are run. The first isolates the noise and tolerance behaviour on sparse adjacency; the
second, reported further below, moves to map scale with a dense grain graph and reaches a
materially different conclusion about the safe tolerance.

The first sweep uses six randomly oriented austenite parents, each transformed through distinct
Kurdjumov-Sachs variants, with adjacency chaining the children of each parent and one contact edge
between consecutive parents — so every microstructure contains real cross-parent boundaries. Three axes
were swept over **25 seeds per cell, 48 cells**:

- orientation noise on the child grains: 0, 0.25, 0.5, 1.0, 2.0 deg (Gaussian, per grain),
- edge tolerance: 1, 2, 3, 5 deg (the one free parameter of the edge test),
- children per parent: 2, 5, 12.

Trials whose planted ground truth is **genuinely ambiguous** at the tolerance under test — a
cross-parent boundary that really does land on the same-parent fingerprint — are counted
separately and excluded from the accuracy statistics. No algorithm can resolve those from
orientations alone, and scoring them as failures would misattribute a physical limit to the
implementation.

## Headline result: no false merges

**Across all 48 cells and roughly 700 judged trials, the false-link rate was exactly zero.** The
edge test never linked a cross-parent boundary whose ground truth was separable. This is the
property that matters most for reconstruction, because a false link merges two parent grains
irreversibly and silently.

For contrast, the angle-only edge test this replaced accepted 52.8% of uniformly random unrelated
boundaries at the same 3 deg tolerance, and recovered only 4 to 7 of 12 planted parents on an
equivalent fixture. See the CHANGELOG entry for that measurement.

## The remaining failure mode is splitting, not merging

With false links eliminated, the only way the partition breaks is a **missed** link: a genuine
same-parent boundary rejected because noise pushed it off the fingerprint. That is governed by the
tolerance relative to the noise:

| noise (deg) | tolerance (deg) | ratio | exact partitions | missed links |
| --- | --- | --- | --- | --- |
| 0.25 | 1.0 | 4x | 83–100% | 0.0–0.3% |
| 0.50 | 1.0 | 2x | 0–48% | 14–17% |
| 0.50 | 2.0 | 4x | 82–100% | 0.0–0.3% |
| 1.00 | 2.0 | 2x | 0–48% | 13–17% |
| 1.00 | 3.0 | 3x | 33–95% | 0.8–1.9% |
| 1.00 | 5.0 | 5x | 80–100% | 0.0–0.3% |
| 2.00 | 5.0 | 2.5x | 0–80% | 0.8–10% |

**Practical rule: set `tolerance_deg` to at least four times the per-grain orientation scatter.**
At 2x the partition collapses; at 4x it is essentially always exact.

That is a *lower* bound only, and on this sparse adjacency it is the only bound — which is
precisely what makes the sparse sweep misleading on its own. The map-scale section below adds an
**upper** bound that is far more restrictive, and it is the binding one for real grain graphs. Do
not read the rule above as endorsing a large tolerance.

Note the counter-intuitive interaction with grain count: **more children per parent makes exact
partition recovery harder at marginal tolerance**, because each additional intra-parent edge is
another chance to miss one and split the cluster (at 1.0 deg noise and 3.0 deg tolerance: 95%
exact with 2 children, 79% with 5, 33% with 12).

## Parent-orientation accuracy improves as sqrt(n)

The quaternion-eigen-mean refinement averages per-member noise rather than inheriting one
member's. Mean parent error at 0.5 deg noise, 3.0 deg tolerance:

| children per parent | measured error (deg) | sigma / sqrt(n) (deg) |
| --- | --- | --- |
| 2 | 0.316 | 0.354 |
| 5 | 0.200 | 0.224 |
| 12 | 0.135 | 0.144 |

The measured errors track $\sigma/\sqrt{n}$ closely, and sit slightly below it. Parent estimates
are therefore *more* accurate than the individual measurements they are built from, which is the
expected behavior of a correct symmetry-aware average.

## The cost of raising tolerance

Tolerance cannot be raised without limit. As it grows, genuinely ambiguous cross-parent boundaries
become common — at 5 deg, roughly 20 of every 25 random microstructures contained at least one
cross-parent boundary sitting on the same-parent fingerprint (versus 0 to 1 of 25 at 1 deg). Those
are real physical coincidences, not algorithmic failures, but they bound how far the tolerance can
usefully go and are why the study reports them separately.

This is the genuine trade-off of the method: too tight and parents split, too loose and distinct
parents become indistinguishable. The window is wide for clean data and narrows as noise rises.

## Map scale changes the answer: the tolerance has an upper bound too

The sweep above chains each parent's grains and gives each parent pair a single contact edge. A
measured grain graph is far denser, and the second sweep models that: parents tile a square grid,
each holding a 3x3 patch of child grains, with four-connected adjacency over all of them — so
every shared parent boundary contributes *several* edges.

That difference is decisive, because the failure is asymmetric. **One chance link anywhere along a
shared boundary merges two parents irreversibly**, so many edges per boundary means many chances.
With 100 parents (900 grains, ~1740 edges):

| relationship | tolerance (deg) | parents recovered of 100 | ambiguous cross-edges |
| --- | --- | --- | --- |
| Kurdjumov-Sachs | 0.5 | 99.7 | 0.06% |
| Kurdjumov-Sachs | 1.0 | 98.0 | 0.37% |
| Kurdjumov-Sachs | 2.0 | 88.7 | 2.28% |
| Kurdjumov-Sachs | **3.0 (the default)** | **69.7** | 6.11% |
| Burgers (bcc→hcp) | 0.5 | 100.0 | 0.00% |
| Burgers (bcc→hcp) | 1.0 | 100.0 | 0.00% |
| Burgers (bcc→hcp) | 2.0 | 98.3 | 0.49% |
| Burgers (bcc→hcp) | 3.0 | 97.0 | 0.74% |

(All at zero added noise, 3 seeds per cell. Adding 0.25 deg scatter does not change the picture:
Kurdjumov-Sachs recovers 98.7 at 1.0 deg and 91.3 at 2.0 deg.)

**The false-link rate among separable boundaries remained zero throughout.** Every one of these
merges came from a *genuinely* ambiguous boundary — two unrelated parents that really do share a
fingerprint-consistent misorientation. This is a property of the relationship and the tolerance,
not a defect in the edge test, and no algorithm working from orientations and a binary edge test
can avoid it.

Two consequences follow.

**The default `tolerance_deg=3.0` is not safe for dense grain graphs under Kurdjumov-Sachs.** It
loses about 30% of parents at map scale. The default is retained because it is appropriate for the
sparse-adjacency and noisier cases in the first sweep, but map-scale work should set it explicitly.
Combined with the noise rule above, the usable window for cubic-cubic Kurdjumov-Sachs on a dense
graph is roughly $4\sigma \le$ tolerance $\lesssim 1^{\circ}$, which means the method needs
orientation scatter below about $0.25^{\circ}$ to be reliable there.

**The window is relationship-dependent, and fingerprint size is why.** Burgers is dramatically more
forgiving than Kurdjumov-Sachs at the same tolerance because it has 12 variants rather than 24, so
its admissible set is far smaller (about 2 800 distinct elements against about 10 700) and occupies
correspondingly less of orientation space. Relationships with more variants are intrinsically
harder to reconstruct from boundary evidence alone.

Because this is a property of the relationship and tolerance rather than of the data,
`ParentGrainReconstructionResult` now reports it directly as `chance_link_probability`, and
`describe()` warns when the expected number of chance links across the tested edges reaches one:

> WARNING: at this tolerance 7.3% of unrelated grain pairs fall within the same-parent
> fingerprint, so about 2 of the 29 tested edge(s) are expected to link by chance alone.

That diagnostic needs no ground truth, so it is available on real data.

## Limitations

This study is synthetic and deliberately narrow. It does **not** establish:

- behavior on measured EBSD data, where noise is neither Gaussian nor independent per grain, and
  grain size varies;
- behavior for non-cubic parents or children (only cubic-cubic Kurdjumov-Sachs was swept);
- behavior on measured grain topology — the map-scale sweep uses square parent blocks with
  four-connected adjacency, which captures the density that matters but not the irregular grain
  shapes, size distribution, or boundary lengths of a real map;
- MTEX parity. The `or_transformation_v1` parity campaign is defined and its PyTex side generated,
  but no MTEX results exist yet because no MATLAB installation was available. No parity claim may
  be made until that campaign has been run and compared.

Closing the first and last of those is the remaining work before reconstruction can be promoted
out of `experimental`.

## References

### Normative

- [Phase Transformation Validation Matrix](phase_transformation_validation_matrix.md)
- [Orientation Relationship Analysis Foundation](../architecture/orientation_relationship_analysis_foundation.md)
- [Benchmark And Tolerance Governance](../standards/benchmark_and_tolerance_governance.md)

### Informative

- Morito, Tanaka, Konishi, Furuhara and Maki, *Acta Materialia* 51 (2003) 1789 — the
  Kurdjumov-Sachs intervariant table underlying the same-parent fingerprint.
