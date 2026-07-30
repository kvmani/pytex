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

Six randomly oriented austenite parents, each transformed through distinct Kurdjumov-Sachs
variants, with adjacency chaining the children of each parent and one contact edge between
consecutive parents — so every microstructure contains real cross-parent boundaries. Three axes
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
At 2x the partition collapses; at 4x it is essentially always exact. The default of 3.0 deg is
therefore appropriate for data with scatter up to about 0.75 deg, which covers well-indexed EBSD,
and should be raised for noisier maps.

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

## Limitations

This study is synthetic and deliberately narrow. It does **not** establish:

- behavior on measured EBSD data, where noise is neither Gaussian nor independent per grain, and
  grain size varies;
- behavior for non-cubic parents or children (only cubic-cubic Kurdjumov-Sachs was swept);
- behavior on realistic map topology — adjacency here is a chain plus single contacts, whereas a
  real grain graph is far more densely connected, which changes both how easily a missed link
  splits a cluster and how many chances there are to make a false one;
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
