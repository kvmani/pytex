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

That is a *lower* bound. There is an upper bound too, from a different mechanism — genuinely
ambiguous cross-parent boundaries — which the sparse topology here barely exercises. The map-scale
section below shows it dominating, and shows the clustering change that recovers most of the loss.
Do not read the rule above on its own as endorsing a large tolerance.

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

## Map scale exposed a ceiling, and it needed an algorithm change

The sweep above chains each parent's grains and gives every parent pair a single contact edge. A
measured grain graph is far denser, and the second sweep models that: parents tile a square grid,
each holding a 3x3 patch of child grains, with four-connected adjacency over all of them — so every
shared parent boundary contributes *several* edges.

That difference is decisive, because the failure is asymmetric. **One chance link anywhere along a
shared boundary merges two parents irreversibly**, so many edges per boundary means many chances.
With connectivity-only clustering (union-find over linked edges), 100 parents / 900 grains /
~1740 edges gave:

| relationship | tolerance (deg) | parents recovered of 100 | ambiguous cross-edges |
| --- | --- | --- | --- |
| Kurdjumov-Sachs | 1.0 | 98.0 | 0.37% |
| Kurdjumov-Sachs | 2.0 | 88.7 | 2.28% |
| Kurdjumov-Sachs | **3.0 (the default)** | **69.7** | 6.11% |
| Burgers (bcc→hcp) | 3.0 | 97.0 | 0.74% |

**The false-link rate among separable boundaries was exactly zero throughout.** Every one of those
merges came from a *genuinely* ambiguous boundary — two unrelated parents that really do share a
fingerprint-consistent misorientation. The edge test was not wrong; **no** edge test can reject
those, because they are indistinguishable from same-parent boundaries by construction.

The window is also strongly relationship-dependent, and fingerprint size is the mechanism. Burgers
is far more forgiving than Kurdjumov-Sachs because both its variant count and its child point group
are smaller: the deduplicated double coset $G_c \left( R\,G_p\,R^{\mathsf{T}} \right) G_c$ holds
**684** distinct elements for Burgers against **10 584** for Kurdjumov-Sachs — a factor of 15 — so
it occupies correspondingly less of orientation space. Relationships with larger admissible sets
are intrinsically harder to reconstruct from boundary evidence alone.

*(Corrected 2026-08-04. This paragraph previously quoted "about 2 800" against "about 10 700", and
both numbers were wrong. The Burgers estimate had treated the hexagonal child as contributing 24
proper operators rather than 12. The Kurdjumov-Sachs figure came from a deduplication that
canonicalized a quaternion's sign on its largest-magnitude component; two components tie in
magnitude for the 90 and 180 degree elements of a crystal point group, so numerically identical
rotations were canonicalized to* $q$ *and* $-q$ *and counted twice — 81 of them. That routine now
deduplicates on the matrices, which carry no sign ambiguity. The correction strengthens this
paragraph's argument rather than weakening it. Both counts, and their independence from lattice
parameters, are pinned in* `tests/unit/test_transformation.py`*.)*

### The fix: connectivity proposes, consistency disposes

Connectivity is not sufficient evidence of a shared parent, so it is no longer treated as such.
Each connected cluster is now split by agreement: every member proposes the parent it implies,
each proposal is scored by how many members it explains, and the best-supported proposal claims its
supporters; unexplained members repeat the vote. A cluster spanning two parents separates because
no single orientation explains all of it, while a genuine single-parent cluster is returned whole.

Re-running the identical sweep:

| relationship | tolerance (deg) | parents recovered of 100 | grains in pure clusters | exact partitions |
| --- | --- | --- | --- | --- |
| Kurdjumov-Sachs | 0.5 | 100.0 | 100.00% | 100% |
| Kurdjumov-Sachs | 1.0 | 100.0 | 100.00% | 100% |
| Kurdjumov-Sachs | 2.0 | 99.7 | 97.81% | 0% |
| Kurdjumov-Sachs | **3.0 (the default)** | **99.7** | **95.19%** | 0% |
| Burgers (bcc→hcp) | 1.0 | 100.0 | 100.00% | 100% |
| Burgers (bcc→hcp) | 3.0 | 99.7 | 98.93% | 33% |

At the default tolerance under Kurdjumov-Sachs this turns **69.7 recovered parents into 99.7**, with
95% of grains landing in clusters drawn from a single true parent. At 1 deg and below the partition
is recovered *exactly* — every grain, every parent — for both relationships, including with 0.25 deg
added scatter. The sparse sweep above is unaffected: its cells are already single-parent clusters,
which the vote returns untouched.

Note the honest gap between the last two columns. "Exact partition" demands that all 900 grains be
grouped perfectly, so a single misplaced grain fails the whole trial; at 2–3 deg a handful of grains
still land in mixed clusters even though the parent *count* is essentially right. Tightening the
tolerance remains the way to close that.

### Small parents are not the weak case, contrary to expectation

A majority vote invites an obvious worry: a parent contributing only one grain has one vote, so a
large neighbour should be able to swallow it. Measured against that hypothesis, it does not happen.
Conditioning on the junction being *genuinely* ambiguous — the only situation where the question
arises, since otherwise the edge test separates the parents outright — a small parent is absorbed
at these rates:

| grains in the small parent | absorbed at 2 deg | absorbed at 3 deg |
| --- | --- | --- |
| 1 | 5.2% | 6.0% |
| 2 | 6.9% | 6.1% |
| 3 | 7.8% | 8.0% |
| 5 | 6.0% | 10.7% |

(Against a 9-grain neighbour, 3000 seeds per row, of which 58–228 produced an ambiguous junction.)

The rate is essentially flat in parent size, and if anything mildly *worse* for larger small
parents. The reason is that vote counts only decide which proposal is considered first; they never
let a majority claim grains it does not explain. Claiming is a per-grain consistency test —
$\mathbf{C}_j^{\mathsf{T}}\mathbf{P}$ near the variant-description set — so a grain from a genuinely
different parent is simply not claimed, however many votes the neighbour has. Absorption requires
the *additional* coincidence that the absorbing parent's own hypothesis also explains the absorbed
grains, which is independent of how many grains either parent contributes.

Because absorption needs both coincidences, and ambiguous junctions are themselves uncommon (that
is what `chance_link_probability` measures), the compound rate is small — consistent with the
95–100% grain purity measured at map scale.

### What still limits it

The residual failures at loose tolerance are the irreducible ones: when two parents share a boundary
inside the fingerprint *and* the vote cannot separate them because the ambiguity extends to the
parent hypotheses themselves. Reporting remains the defence — because the coincidence rate depends
only on the relationship and the tolerance, never on the data,
`ParentGrainReconstructionResult` exposes `chance_link_probability`, and `describe()` warns when the
expected number of chance links across the tested edges reaches one:

> WARNING: at this tolerance 7.3% of unrelated grain pairs fall within the same-parent
> fingerprint, so about 2 of the 29 tested edge(s) are expected to link by chance alone.

That diagnostic needs no ground truth, so it is available on real data.

## Limitations

This study is synthetic and deliberately narrow. It does **not** establish:

- behavior on measured EBSD data, where noise is neither Gaussian nor independent per grain, and
  grain size varies;
- behavior beyond the two relationships swept — the noise/tolerance sweep is cubic-cubic
  Kurdjumov-Sachs only, and the map-scale sweep adds Burgers (bcc→hcp); other symmetries and
  variant counts are untested, and fingerprint size is known to matter;
- behavior on measured grain topology — the map-scale sweep uses square parent blocks with
  four-connected adjacency, which captures the density that matters but not the irregular grain
  shapes, size distribution, or boundary lengths of a real map;
- (Small parents were expected to be the weak case and measured otherwise — see below. The
  singleton caveat in `describe()` still applies for a different reason: a parent with no
  orientation-consistent neighbour at all is symmetry-ambiguous regardless of the clustering.)
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
