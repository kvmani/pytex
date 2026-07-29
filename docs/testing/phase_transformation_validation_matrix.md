# Phase Transformation Validation Matrix

This document is the authoritative validation ledger for PyTex phase-transformation workflows.

## Status Keys

- `implemented`: automated coverage and validation notes exist for the current category
- `foundational`: the implementation exists and is scientifically structured, but the external or
  literature-backed validation surface is not yet complete
- `planned`: the category is accepted but not yet validated adequately
- `n/a`: not applicable to current PyTex scope, with explanation

## Matrix

| Area | Baseline | Status | Notes |
| --- | --- | --- | --- |
| Orientation-relationship semantic contracts | Core invariant tests and canonical data-model doctrine | implemented | Parent and child phases, frames, and provenance are checked through automated tests. |
| Orientation-relationship construction from explicit plane-direction correspondence | Deterministic right-handed basis construction plus unit tests | implemented | `OrientationRelationship.from_parallel_plane_direction(...)` is covered for matrix recovery and phase mismatch rejection. |
| Named Bain correspondence helper | Deterministic correspondence construction plus cubic-phase guards | implemented | `OrientationRelationship.from_bain_correspondence(...)` is covered for the stated `(001)_p || (001)_c`, `[110]_p || [100]_c` mapping and rejects non-cubic parents. |
| Named Nishiyama-Wassermann helper | Deterministic correspondence construction plus cubic-phase guards | implemented | `OrientationRelationship.from_nishiyama_wassermann_correspondence(...)` is covered for the stated `(111)_p || (011)_c`, `[1-10]_p || [100]_c` mapping and rejects non-cubic children. |
| Variant generation and uniqueness | Internal deterministic tests plus symmetry-backed invariants | implemented | `TransformationVariant` generation is covered by unit tests and benchmark manifests. |
| Variant-indexed predicted child orientations | Internal deterministic tests and manifest-backed workflow identity | implemented | `PhaseTransformationRecord.predicted_child_orientations()` now respects explicit variant assignments. |
| Transformation manifest schema | Stable JSON schema plus round-trip tests | implemented | `TransformationManifest` now records dedicated transformation workflow context. |
| Experimental parent-candidate scoring | Internal deterministic tests and benchmark identity | implemented | Candidate-parent ranking is staged under `pytex.experimental` with explicit non-stable status. |
| Literature-tracked starter transformation families | Textbook Bain and Nishiyama-Wassermann correspondences | foundational | PyTex now tracks named literature-facing helpers and tests, but broader family breadth and curated datasets remain ahead. |
| Index correspondence (parent-child Miller mapping) | Analytic identities: defining parallelisms (KS, Bain, Burgers), inverse-transpose and zone-law invariants, round trips, KS close-packed-group counts across all 24 variants | implemented | `correspondence_direct`/`correspondence_reciprocal`, `map_plane_to_child`/`map_direction_to_child` and parent-inverses, variant tables; rationalized indices with atan2 angular residuals; hexagonal Miller-Bravais covered; worked examples in the transformation gallery. |
| OR misorientation representation | Literature axis/angle representatives: KS 42.85 deg <0.968 0.178 0.178>, NW 45.99 deg, Bain 45 deg <100> | implemented | `OrientationRelationship.misorientation()` returns the deterministic disorientation representative; pinned in tests and a worked example. |
| OR deviation metric | Analytic zero on exact synthetic children with planted-variant recovery; documented KS-GT (2.40 deg) and NW-GT (2.86 deg) separations reproduced | implemented | `or_deviation(...)` returns per-pair minimal deviations and best-variant indices with aggregate statistics; the entry point for OR fitting. |
| Parallelism finders over symmetry families | KS close-packed pairing: each of the 24 variants pairs exactly one {111} member with a {011} child plane (and one <110> member with a <111> child direction) at zero deviation | implemented | `find_parallel_planes` / `find_parallel_directions` enumerate integer symmetry orbits and report per-variant matches within tolerance as a typed `ParallelismReport`. |
| Explainable OR reports (`describe()`) | Substring-validated prose: conventions, defining parallelisms, misorientation representative, deviation statistics, variant frequencies, ambiguity flags | implemented | `describe()` on `OrientationRelationship`, correspondence results, `ORDeviationReport`, `ParallelismReport`, `VariantSelectionReport`, `ParentReconstructionReport` per the development-guide explainable-results doctrine. |
| OR fitting from measured pairs | Analytic recovery: exact GT pairs refit GT identically from both a GT and a KS starting nominal (reported nominal distance reproduces the documented 2.40 deg KS-GT separation); seeded 0.5-deg-noise fits land within 0.15 deg of truth with residuals at the noise level | implemented | `fit_orientation_relationship(...)` performs symmetry-aligned quaternion eigen-mean averaging with iterative realignment; MTEX `calcParent2Child` parity comparison and a worked example remain queued. |
| Variant packet classification (close-packed groups) | Morito lath-martensite hierarchy: KS {111} yields four packets of six variants; Burgers {110} yields six groups of two | implemented | `variant_close_packed_groups(...)` labels variants by the parent family member each carries into exact parallelism; validated with the full 24-variant single-parent fixture below. |
| Lath-martensite structure fixture (single parent, all 24 KS variants) | Morito et al. block/packet structure: reconstruction gathers all 24 children into one parent recovered exactly; variant selection recovers every planted index; packets 4x6 | implemented | End-to-end literature-structure validation of reconstruction + selection + packet classification on one fixture. External measured-data fixtures remain queued. |
| Variant pole figures (predicted overlays) | Packet-plane coincidence: every KS variant's predicted {011} pole set contains the specimen-frame normal of its packet's parent {111} member (all 24 variants, exact) | implemented | `variant_pole_figure(...)` computes specimen-frame poles per variant under the canonical composition; `plot_variant_pole_figure(...)` renders the color-per-variant stereographic overlay. |
| OR-fitting worked example | Gallery example `or-fit-recovers-gt-from-ks-nominal`: zero residual + documented 2.40 deg KS-GT separation, computed live | implemented | Closes the queued documentation gap for `fit_orientation_relationship`. |
| Transformation deformation gradients (Bain-strain family) | Textbook Bain principal stretches (sqrt(2) a_c/a_p twice, a_c/a_p once) and volume ratio 2(a_c/a_p)^3 exact; KS/NW share Bain stretches with polar rotations 11.06/9.74 deg (the literature rigid rotations from Bain) | implemented | `deformation_gradient()` uses the nearest-integer lattice correspondence and polar decomposition; describe() reports principal strains, volume change, and residual rotation. |
| Named Shoji-Nishiyama (fcc-hcp) correspondence | Literature variant count (4, one per {111} parent plane), defining parallelism (111)->(0001) exact, one-variant-per-packet structure | implemented | `from_shoji_nishiyama_correspondence(...)` with cubic/hexagonal guards; `standard_fcc_hcp_relationships` catalog. |
| OR identification from boundaries (experimental, F7 first stage) | Synthetic microstructures: the generating relationship wins with mean fingerprint distance 0 (KS and GT cases) with >1 deg margins over all other fcc-bcc candidates; 0.3-deg noise preserves the ranking | foundational | `pytex.experimental.identify_orientation_relationship(...)` scores boundary misorientations against each candidate's double-coset intervariant fingerprint; no parent orientations required. The second stage (rotation refinement) is a separate row below. |
| OR rotation refinement from boundaries (experimental, F7 second stage) | Exact GT boundaries + KS nominal: converges to the true GT rotation (symmetry-reduced distance ~0 at the 1e-6 deg round-trip noise floor) with the 2.404 deg KS-GT update reported; 0.3-deg-noise KS boundaries + NW nominal recover KS to <0.5 deg | foundational | `pytex.experimental.refine_orientation_relationship_from_boundaries(...)`: alternates nearest-double-coset-element assignment with least-squares rotation updates on smooth $2\sin(\theta/2)$ chordal residuals (SciPy LM); rotation identifiable up to coset symmetry, distances reported symmetry-reduced. |
| Named Pitsch-Schrader (hcp-bcc) correspondence | Defining parallelism (0001)->(110) exact; literature 5.26 deg separation from the inverse Burgers relationship; 3 variants (internally derived orbit count) | implemented | `from_pitsch_schrader_correspondence(...)` with hexagonal-parent/cubic-child guards; `standard_hcp_bcc_relationships` catalog bundles PS + inverse Burgers. |
| Named Bagaryatsky and Isaichev (ferrite-cementite) correspondences | All three Bagaryatsky axis parallelisms exact in the Pnma setting ([1-1-1]->[100], [211]->[010], (0-11)->(001), Bhadeshia MST 34 (2018)); Isaichev (101)->(031) exact; symmetry-reduced Bagaryatsky-Isaichev separation 3.586 deg about exactly the cementite a-axis (literature ~3.8 deg, axial-ratio-dependent); 12 and 24 variants (internally derived orbit counts) | implemented | `from_bagaryatsky_correspondence(...)` / `from_isaichev_correspondence(...)` with cubic-parent/orthorhombic-child guards (first orthorhombic-child constructors); `standard_ferrite_cementite_relationships` catalog. |
| Named Potter (hcp-bcc) correspondence | Defining parallelisms exact: (01-11)->{110} and <2-1-10>-><111> (Potter 1973, V-N precipitates); basal-plane image residual equals the symmetry-reduced separation from inverse Burgers (1.370 deg at c/a=4.68/2.95; literature "~2 deg" is c/a-dependent); 12 variants (internally derived orbit count) | implemented | `from_potter_correspondence(...)` with hexagonal-parent/cubic-child guards; appended to the `standard_hcp_bcc_relationships` catalog. |
| Map-scale parent-grain reconstruction (experimental v1) | Synthetic planted-parent recovery: three KS parents with five children each are re-grouped exactly (cross-parent edges rejected by the intervariant fingerprint), parents recovered to 0 deg; 0.3-deg-noise case keeps the partition with residuals at the noise level | foundational | `pytex.experimental.reconstruct_parent_grains(...)`: intervariant-fingerprint edge test, union-find clustering, per-cluster candidate-parent scoring; `describe()` states ambiguity caveats. Literature fixtures (martensite→austenite, alpha→beta Ti) and EBSD grain-graph wiring remain ahead before stabilization. |

## Current Posture

PyTex now has explicit transformation primitives, a dedicated manifest contract, named literature
starter correspondences, and a bounded experimental reconstruction-scoring surface. The semantics
are no longer ad hoc, but the validation posture is still foundational because broad
literature-backed families and curated reconstruction datasets are not yet in place.

## Evidence Hardening Queue

Before broader transformation algorithms move out of experimental scope, the next validation pass
should add:

- a curated orientation-relationship family ledger that records source convention, plane-direction
  correspondences, and PyTex canonical mapping
- literature-backed variant-count and variant-orbit checks beyond the starter Bain and
  Nishiyama-Wassermann helpers
- at least one small parent-reconstruction benchmark with explicit ambiguity and failure-mode notes

## References

### Normative

- `strategy.md`
- `../architecture/phase_transformation_foundation.md`

### Informative

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*
- `../../benchmarks/transformation/variant_prediction_benchmark_manifest.json`
