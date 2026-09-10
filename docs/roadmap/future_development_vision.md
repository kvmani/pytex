# PyTex Future Development Vision

**Repository-wide assessment and proposed research program — 10 September 2026.**
Baseline: commit `e497c8a`, package version 0.9.0. Status: **proposal for prioritization**;
the features and acceptance targets below are not claims of delivered capability.

PyTex should become a reference environment for **evidence-backed phase-transformation
crystallography**: identify a relationship, reconstruct the parent, quantify what remains
ambiguous, connect the result to interfaces and diffraction, and export an analysis another
laboratory can reproduce. Its next advantage should come from trustworthy answers across
measurements, rather than the number of individually available calculations.

The recommended first program is **measured OR validation, robust uncertainty-aware OR fitting,
reconstruction reliability, and quantitative interface geometry**. Texture uncertainty and
multiphase diffraction follow as a second connected program. Instrument-aware experiment design
is the strongest longer-term research bet. Continue useful smaller improvements, but do not
start all of these programs at once.

## 1. Authority, coverage, and how to use this document

This vision is subordinate to the [mission](https://github.com/kvmani/pytex/blob/e497c8a/mission.md),
[specifications](https://github.com/kvmani/pytex/blob/e497c8a/specifications.md), and
[governing development guide](critical_review_and_development_guide.md). It extends the
[OR foundation](../architecture/orientation_relationship_analysis_foundation.md) and reconciles
the proposals in the [interface and texture vision](vision_interface_crystallography_and_texture_quantification.md)
against current implementation. Existing foundation feature IDs are unchanged; **V01–V24** are
portfolio IDs for this document. A proposed API name is not a reservation of stable public API.

The review inventoried every tracked package module, read public declarations and module
contracts across every scientific subsystem, and inspected representative implementations,
tests, fixtures, schemas, application services, CI, executable examples, and roadmap records.
It is a repository-wide capability and development assessment, **not a line-by-line correctness
proof or an independent replication of every scientific result**. Repository counts below were
computed from tracked files and Python syntax trees; historical test outcomes remain attributed
to their recorded run. New verification is recorded in the
[progress ledger](https://github.com/kvmani/pytex/blob/main/docs/development/active_task_progress.md).

Primary external literature and maintained tool documentation were checked on 10 September
2026. This is a targeted prior-art check, not a systematic literature review. Every proposed
methods paper needs a renewed search before its claim and experiment plan are frozen.

For planning, read the capability assessment, choose the next stage in Section 7, and turn each
selected feature card into an implementation plan. Each card states the **additional scope,
use cases, impact, difficulty, publication opportunity, tests, and validation**. Section 6 gives
shared dataset protocols and measurable acceptance targets. Adopt targets before seeing held-out
results; never move a threshold merely to pass a benchmark.

## 2. What the repository actually provides

The tracked baseline contains **172 Python package modules, 137,029 Python source lines,
621 entries in the root `__all__`, 166 Python test files, 36 notebooks, and 11 JSON schema
files**. Counts include private modules and package initializers; source-line counts include
comments and docstrings. These are inventory measures, not scientific quality scores.

Source links below are pinned to the reviewed commit so later implementations do not silently
change the meaning of the assessment. “Implemented” means an identifiable implementation with
repository test evidence; it does not mean externally validated in every regime.

| Area | Implemented baseline and evidence | Additional opportunity |
| --- | --- | --- |
| Canonical crystallography | Frames and frame graphs, direct/reciprocal bases, point/space-group distinction, Miller-Bravais forms, orientation representations, semantic batches, symbols and notation. [Core source](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/core). | Shared uncertainty semantics, tighter API organization, and controlled growth of result contracts. |
| OR analysis | Correspondences, variants and groups, weighted measured-pair fitting, seedless characterization, deformation gradients, dossier, stereograms and composite views. [Transformation](https://github.com/kvmani/pytex/blob/e497c8a/src/pytex/core/transformation.py), [weighted tests](https://github.com/kvmani/pytex/blob/e497c8a/tests/unit/test_or_weighted_fitting.py). | Robust fitting, calibrated confidence and alternatives; interface calculations. Scalar evidence weights and residual scatter are already present but are not confidence intervals. |
| Parent reconstruction | Experimental boundary-based OR identification/refinement, full-rotation fingerprints, parent clustering, consistency splitting, chance-link and singleton diagnostics. [Experimental modules](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/experimental), [reconstruction tests](https://github.com/kvmani/pytex/blob/e497c8a/tests/unit/test_parent_grain_reconstruction.py). | Measured benchmark, external comparison, twin-aware joint inference, explicit abstention and reliable stabilization. This is not a greenfield reconstruction project. |
| Texture | Discrete/harmonic ODFs, kernels, fibres, named-component fractions and constrained fixed-component fitting, PF inversion, random-standard defocus correction, ghost correction, Kearns estimators. [Texture source](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/texture). | Axial specimen symmetry, component-shape refinement, uncertainty and cross-measurement consistency. A named fibre object does not implement continuous specimen symmetry. |
| EBSD | Canonical multiphase maps, square/hex topology, KAM, grains, boundary graphs, CSL/twins, local metrics, curvature/GND and texture workflows. [Map model](https://github.com/kvmani/pytex/blob/e497c8a/src/pytex/ebsd/models.py), [GND](https://github.com/kvmani/pytex/blob/e497c8a/src/pytex/ebsd/gnd.py). | Out-of-core execution, hex/irregular curvature, spatial uncertainty, 3-D and time-series identity. Existing surface-map GND is a lower-bound analysis. |
| Import and interoperability | ANG, CTF, EDAX OIM HDF5, XRDML and LaboTex, CIF/structure, orix, KikuchiPy and abTEM/ASE bridges. [Adapters](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/adapters). | Explicit additional HDF5 dialects, large-array transport, 4D-STEM result import and frame-calibrated cross-modal registration. Readable data are not automatically experimentally calibrated data. |
| Measured powder XRD | XY/XRDML measurement I/O, background, peaks, indexing, candidate phase ranking, corrections, size/strain tools, lattice determination including Le Bail, and bounded single-phase profile refinement. [Diffraction source](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/diffraction), [Rietveld tests](https://github.com/kvmani/pytex/blob/e497c8a/tests/unit/test_xrd_rietveld.py). | Joint multiphase fractions, texture/scale identifiability and calibrated residual stress. Do not propose “add measured XRD” or “add Rietveld” as absent features. |
| Electron diffraction | Geometric SAED solving, lattice fitting, composite OR patterns, Kikuchi maps, finite-thickness kinematics, HOLZ, CBED, Bloch-wave foundations and diffraction-group reasoning. [Solving](https://github.com/kvmani/pytex/blob/e497c8a/src/pytex/diffraction/solving.py), [dynamical model](https://github.com/kvmani/pytex/blob/e497c8a/src/pytex/diffraction/dynamical.py). | Automated image metrology, rings, model-discrepancy validation and shared calibration uncertainty. Existing dynamical support is not validation of all specimen regimes. |
| TEM operation | Stage models, calibration, orientation recovery, reachability, routing, zone-axis atlas and symmetry-ambiguity analysis. [TEM source](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/tem). | Select the next measurement by expected discrimination, stage feasibility, time and uncertainty; current routing is a substantial starting point. |
| HREM | Pure-Python phase-object simulation, abTEM multislice dispatch, microscope aberrations, azimuthal CTF, atomic snapshots and a workbench. [HREM](https://github.com/kvmani/pytex/blob/e497c8a/src/pytex/diffraction/hrem.py), [abTEM adapter](https://github.com/kvmani/pytex/blob/e497c8a/src/pytex/adapters/abtem.py). | Validated ensembles, sampling convergence and measured-image comparison. Optical-coefficient transfer parity does not establish full image equivalence. |
| Material properties | Elastic tensors, directional moduli, homogenization, slip/Schmid and Taylor factors. [Properties](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/properties). | HCP slip/twinning breadth, boundary transmission and simulation adapters; avoid recreating a full crystal-plasticity solver. |
| Visualization and application | Semantic plot builders, themes, composite crystal scenes, dossiers, shared offline desktop/browser workbench, typed operation registry, exports, progress and logbook. [Plotting](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/plotting), [application](https://github.com/kvmani/pytex/tree/e497c8a/src/pytex/app). | Uncertainty-aware views, assumption comparison, replayable studies and ADP geometry. A new panel alone seldom creates a scientific contribution. |
| Documentation and quality | Sphinx/MyST, executable examples, notebooks, canonical SVGs, parity/validation ledgers, strict typing, warning controls and browser tests. [CI](https://github.com/kvmani/pytex/blob/e497c8a/.github/workflows/ci.yml), [test tree](https://github.com/kvmani/pytex/tree/e497c8a/tests). | Evidence freshness, contract completeness and performance regression coverage. Extend the existing infrastructure instead of proposing it from scratch. |

### 2.1 The most consequential gaps

1. **External evidence lags the flagship.** The tracked campaign-result directory for MTEX
   contains only `.gitkeep`; older core parity fixtures do exist elsewhere. The OR foundation
   explicitly withholds measured reconstruction and OR-fitting parity claims. No campaign
   generator, synthetic fixture, or accurate-looking map closes that gap by itself.
2. **Uncertainty is fragmented.** Calibration metadata, measurement uncertainties, fit scatter
   and some least-squares parameter errors exist. There is no reviewed common path that
   propagates correlated measurement uncertainty through OR, reconstruction and texture results.
3. **The interface is still absent.** `ORDossier.interface` is explicitly `None`. Composite
   drawing, index correspondence and deformation-gradient calculation must not be described as
   computed interface misfit or PTMC habit-plane prediction.
4. **Scale is not solved by a few bounded kernels.** Weighted OR alignment is chunked, while
   the scan/map workflow remains predominantly in-memory. Global graph reconciliation and
   halo handling are separate requirements from reading an HDF5 file in blocks.
5. **Roadmap freshness is a scientific governance issue.** Old documents still list shipped
   XRD, HREM, GUI and CI capabilities as missing. Treat dated historical sections as history.
   The August guide's higher-level priorities still govern, with this assessment supplying the
   current scope of the proposed extensions.
6. **Maintenance pressure is visible.** `transformation.py` has 4,916 lines and the root exports
   number 621. This motivates staged internal decomposition and curated domain entry points,
   not a disruptive blanket rewrite or removal of existing imports.

The preceding release ledger reports 91.5180669597191% coverage, 67 browser tests passing
without retries, zero Sphinx warnings and eight green CI jobs for release commit `233f638`.
Those are **historical release results**. The current CI coverage threshold is 87%, below the
observed release coverage; propose ratcheting it only against a reproducible baseline and
documenting the aggregation method. Neither high coverage nor tool parity proves physical truth.

## 3. Prior-art boundaries and the publication thesis

The novelty ratings below concern a **testable contribution beyond known methods**. They are
not predictions of acceptance, and implementation difficulty is not a proxy for novelty.

| Proposed theme | Existing work to acknowledge | Defensible PyTex contribution to investigate |
| --- | --- | --- |
| Parent reconstruction and OR discovery | [Niessen et al. (2022)](https://journals.iucr.org/j/issues/2022/01/00/nb5309/) already describe generic reconstruction, OR analysis and several transformation systems; [Hielscher et al.'s variant graph](https://arxiv.org/abs/2201.02103) is established prior art. | Calibrated ambiguity, robustness to sampling/segmentation and independently measured recovery, with complete convention-preserving interchange. Generic reconstruction itself is not novel. |
| Measured habit planes and interface matching | [Habit-plane determination from reconstructed maps](https://arxiv.org/abs/2303.07750) and [pymatgen interface matching](https://pymatgen.org/pymatgen.analysis.interfaces.html) cover substantial neighboring ground. | Connect measured OR uncertainty, variant-resolved interface candidates and forward diffraction tests; quantify what a 2-D trace cannot identify. |
| Joint texture and powder analysis | [GSAS-II](https://doi.org/10.1107/S0021889813003531) is an established refinement platform; [quartz/MAUD coordinate validation](https://journals.iucr.org/j/issues/2023/06/00/xx5024/) demonstrates texture refinement and cross-frame concerns. | Diagnose cross-modality inconsistency and confounding between texture, calibration and phase fraction, with independently measured reference mixtures. Texture correction alone is not novel. |
| Electron simulation and mapping | [abTEM multislice](https://abtem.github.io/doc/user_guide/walkthrough/multislice.html) and [pyxem orientation mapping](https://arxiv.org/abs/2111.07347) establish simulation and 4D-STEM capabilities. | Link uncertainty, OR hypotheses and instrument-constrained experimental discrimination through canonical PyTex contracts. Do not claim a new multislice method for an adapter. |
| Mechanical predictions | [DAMASK](https://damask-multiphysics.org/) supplies multiphysics crystal plasticity; [MTEX slip transmission](https://mtex-toolbox.github.io/SlipTransmission.html) is adjacent prior art. | Trace reconstructed-parent/variant uncertainty into useful property predictions using adapters and independently measured outcomes. |
| Scientific interchange | [NeXus application definitions](https://manual.nexusformat.org/classes/applications/index.html) already standardize experiment-specific information. | A tested semantic mapping and portable study replay for this scientific workflow; JSON or HDF5 serialization alone is not a new method. |

Three distinct publication routes are realistic: a **software paper** about a useful,
validated research instrument; a **methods paper** proving a new estimator or diagnostic;
and an **application paper** answering a materials question with independent experiments.
All three require evidence, but they require different evidence. Keep these routes separate.

## 4. Feature portfolio and assessment scale

Difficulty estimates cover implementation, scientific tests, examples and documentation for
the bounded first version, assuming an experienced scientific Python developer with access to
a domain reviewer. They exclude experiment scheduling and journal review. **S: 1–2 person-weeks;
M: 3–6; L: 7–12; XL: 13–24+**. These are planning ranges, not delivery promises; an unsuccessful
research investigation may stop without a stable feature. Data access is shown separately.

Impact: **5** changes a flagship research workflow; **4** removes a major practical barrier;
**3** materially improves an adjacent workflow. Novelty: **E** established method/engineering,
**I** integration or application opportunity, **R** research hypothesis requiring a distinct
method and comparison. R is not a certification of originality. “Now” means the first program;
“Next” means after its prerequisites; “Explore” requires a bounded feasibility study.

| ID | Feature | Priority | Impact | Difficulty | Novelty |
| --- | --- | --- | --- | --- | --- |
| V01 | Measured OR and reconstruction benchmark | Now | 5 | M + external data | I |
| V02 | Robust OR inference with calibrated uncertainty | Now | 5 | L | R |
| V03 | Ambiguity-aware joint parent reconstruction | Now, after V01/V02 | 5 | XL | R |
| V04 | Variant-resolved interface misfit and coincidence | Now | 5 | L | I/R |
| V05 | Habit-plane inference from traces and sections | Next | 4 | L | R |
| V06 | PTMC and invariant-plane compatibility | Explore | 5 | XL | I/R |
| V07 | Stress-assisted variant selection | Next | 4 | L | I/R |
| V08 | Axial specimen symmetry | Next; bounded early win | 4 | M | E |
| V09 | Texture uncertainty and component refinement | Next | 5 | L | I/R |
| V10 | Joint EBSD–diffraction texture consistency | Next | 5 | XL | R |
| V11 | Texture-aware multiphase powder quantification | Next | 5 | L | I/R |
| V12 | Calibrated diffraction residual-stress analysis | Next | 4 | L | I |
| V13 | Out-of-core EBSD and grain-graph execution | Now, infrastructure track | 5 | L | E/I |
| V14 | Hex/irregular curvature and GND uncertainty | Next | 4 | M–L | I |
| V15 | Vendor and 4D-STEM semantic bridges | Next, demand-led | 4 | M per bridge | E/I |
| V16 | Image-to-SAED and ring-pattern metrology | Next | 4 | L | I |
| V17 | Information-guided TEM measurement planning | Explore | 5 | XL | R |
| V18 | Electron-model validation and physical ensembles | Next | 4 | L | I/R |
| V19 | 3-D and time-resolved transformation tracking | Explore | 5 | XL | R |
| V20 | HCP mechanics and simulation-ready microstructures | Next | 4 | L | I |
| V21 | Neutron/TOF texture workflow | Explore, partner-led | 4 | XL | I/R |
| V22 | Replayable study bundles and assumption comparison | Now, infrastructure track | 5 | M | I |
| V23 | API sustainability and performance evidence | Now, continuous | 4 | M initial | E |
| V24 | Anisotropic displacement and quantitative crystal scenes | Later | 3 | M | E/I |

### V01 — A measured OR and reconstruction benchmark

**Scope and use cases.** Extend existing manifests and campaigns into a published benchmark
covering partially transformed steel, fully transformed steel, and a Burgers Ti/Zr system.
Separate measured parent truth, expert annotations and simulated truth. Researchers use it to
choose tolerances and compare algorithms; maintainers use it as the stabilization gate.

**Impact, difficulty, novelty.** Impact 5, M plus data/MTEX access. Primarily a validation and
dataset contribution (I), potentially supporting a software or benchmark paper. It cannot claim
a new reconstruction algorithm. Dependencies: an experimental partner and comparator operator.

**Tests and validation.** Protocol A below; test frame/variant label remapping, hashes,
deterministic metrics and label-permutation invariance. Execute the prepared MTEX campaign
with version and settings recorded. Compare full boundary rotations and parent partitions,
not only angle histograms. Release the exact comparison scripts and licensed fixtures or
documented acquisition instructions. Promotion remains blocked if measured truth is unavailable.

### V02 — Robust OR inference with calibrated uncertainty

**Scope and use cases.** Build on the existing weighted Markley fit: explicit outlier models,
multiple initializations, symmetry-aligned bootstrap or likelihood-based uncertainty regions,
candidate alternatives and an abstention rule. Preserve excluded observations. Distinguish
user quality weights from calibrated measurement variances and grain-level dependence.
Use for noisy paired EBSD data and discrimination between nearby nominal ORs.

**Impact, difficulty, novelty.** Impact 5, L, research hypothesis R. A publishable result would
show calibrated uncertainty and fewer confidently wrong OR labels than the existing fit and
a matched external baseline. A generic robust loss or bootstrap implementation alone is E.
Depends on V01 and a shared covariance/frame contract introduced within this feature.

**Tests and validation.** Exact rotations, quaternion sign/symmetry invariance, zero-weight
exclusion and common weight rescaling; contamination, minority variants and correlated samples.
Use Protocol B with clustered resampling and held-out specimen-level evaluation. Report interval
coverage, angular error and abstention versus error; retain multimodal answers when necessary.

### V03 — Joint reconstruction with visible alternatives

**Scope and use cases.** Extend the experimental graph engine to alternate OR refinement,
parent assignment and consistency checks, with retained-parent anchors, twins and per-grain
alternative candidates. Provide ambiguity maps and an unassigned state. Do not force every
grain into a single parent. Use for consumed-austenite steel and transformed Ti/Zr maps.

**Impact, difficulty, novelty.** Impact 5, XL, R. The research claim is reliable joint inference
under ambiguity, tested against current PyTex and MTEX variant-graph reconstruction. Existing
single-parent splitting and chance-link estimates must remain credited. Dependencies: V01,
V02; V13 before a large-data claim.

**Tests and validation.** Protocol A; exact planted partitions, parent twins, missing variants,
singleton parents, mixed phases and segmentation perturbations. Require improvement in
false-merge/false-split and selective-risk metrics, not merely more assigned pixels. Verify
that input order and graph relabeling do not change scientifically equivalent results.
Stabilize only the measured and documented regimes; unsupported regimes remain experimental.

### V04 — Interface misfit and coincidence for each variant

**Scope and use cases.** Populate the currently empty dossier interface block from two phases,
an OR variant and declared interface planes. Compute in-plane bases, strain/misfit, bounded
supercell matches, residuals and candidate rankings; render a measured-scale interface patch.
Expose area/index limits and alternate correspondences. Use for precipitates, martensite
interfaces, epitaxy and teaching. Defer atomistic relaxation and interface-energy prediction.

**Impact, difficulty, novelty.** Impact 5, L, I/R. Matching itself is established; the opportunity
is uncertainty-aware, variant-resolved connection from measurement to interface to diffraction.
Compare against pymatgen's Zur–McGill implementation rather than claiming the first matcher.
Depends on existing lattice/OR/dossier surfaces; V02 adds uncertainty later.

**Tests and validation.** Protocol C; identical lattices give zero misfit, known rational
supercells match exactly, in-plane unimodular basis changes preserve invariant measures, and
phase exchange follows the declared inverse transformation. Verify bounded exhaustive cases
independently. Geometry alone must never be reported as a measured habit plane or interface energy.

### V05 — Habit-plane evidence from measured traces

**Scope and use cases.** Combine map orientations, specimen surface normals and segmented
interface traces to rank compatible plane families and report uncertainty. Support multiple
nonparallel sections or 3-D normals when available. A single surface trace constrains a family
of 3-D normals and must not yield artificial uniqueness. Use for lath/plate morphology studies.

**Impact, difficulty, novelty.** Impact 4, L, R only if ambiguity/selection-bias treatment improves
on published trace-distribution methods. Depends on V01/V03 for reconstructed orientations,
V04 for comparison, and independently annotated microscopy.

**Tests and validation.** Analytic plane/surface intersections, sign-equivalent normals, tilted
sections and near-parallel ill-conditioning; validate against independent 3-D microscopy or
serial sections in Protocol C. Compare the published MTEX habit-plane approach. With only
one 2-D map, report compatibility or a population model conditional on its sampling assumptions.

### V06 — PTMC and invariant-plane compatibility

**Scope and use cases.** Starting from an explicitly selected lattice correspondence and
deformation, solve a bounded family of lattice-invariant shear/twin systems for invariant-plane
strain, shape strain and candidate OR/habit planes. Return all admissible solutions and
no-solution diagnostics. Use for martensitic compatibility and alloy-parameter sensitivity.

**Impact, difficulty, novelty.** Impact 5, XL, I/R; classical PTMC implementation is established.
The possible methods contribution is uncertainty/identifiability across correspondence choices.
Requires a transformation-crystallography reviewer, V04, and a signed scientific specification.
Do not treat nearest-integer correspondence as a unique physical mechanism.

**Tests and validation.** Exact compatible and incompatible stretches, rank-one residuals,
rigid-frame invariance and a cited literature case with all lattice parameters and shear
systems pinned. Compare independently derived solutions and measured habit-plane evidence.
Time-box initial feasibility to four person-weeks; stop if correspondence selection and
independent numerical reference cases cannot be fixed. Do not block the first paper on PTMC.

### V07 — Stress-assisted variant selection

**Scope and use cases.** Combine variant deformation tensors with a declared stress state to
report mechanical-work rankings and population predictions under an explicit statistical model.
Keep mechanical driving work separate from actual transformation probability. Use for loaded
transformations, additive-manufactured texture and heat-treatment studies.

**Impact, difficulty, novelty.** Impact 4, L, I/R. Stress-work ranking is established; a defensible
application or methods paper needs prediction on a held-out load path. Depends on V02/V04,
validated deformation correspondences and measured stress/temperature history.

**Tests and validation.** Zero stress removes mechanical preference; hydrostatic loading gives
equal work for variants with equal volume change; simultaneous frame rotation preserves work.
Pin stress and finite/small-strain conventions before implementation. Compare measured variant
fractions against uniform selection and mechanics-only baselines; report fitting versus prediction
separately. Failure to predict held-out fractions limits the output to a driving-work descriptor.

### V08 — Continuous axial specimen symmetry

**Scope and use cases.** Add an explicit continuous specimen-symmetry contract and exact harmonic
projection about a declared specimen axis; provide a discrete approximation only with stated
resolution and convergence evidence. Connect ODF, PF, IPF and Kearns reporting. Use for wires,
rods and nominally axisymmetric tubes; keep the assumption opt-in.

**Impact, difficulty, novelty.** Impact 4, M, E; valuable completeness, low standalone publication
novelty. Builds on specimen symmetry and harmonics, not a new private symmetry model.

**Tests and validation.** Projection idempotence, invariance to rotations around the declared
axis, unit ODF normalization, equal transverse Kearns factors and preservation of their sum.
Validate a uniform distribution, ideal axial fibre, and deliberately asymmetric data that must
change when the assumption is imposed. Pin a matched MTEX symmetry projection and quantify
discrete convergence. Rotating the requested axis and data together must preserve predictions.

### V09 — Texture uncertainty and refined component models

**Scope and use cases.** Add grain/block bootstrap intervals for Kearns factors, component
fractions and texture index; separately refine named-component centers and widths with
nonnegative fractions and model-selection diagnostics. Preserve current fixed-center fitting.
Expose sensitivity to smoothing, pole coverage, ghost correction and specimen symmetry.
Use for comparing processing routes and deciding whether apparent texture changes are real.

**Impact, difficulty, novelty.** Impact 5, L, I/R. Publication needs calibrated intervals or a
demonstrated bias/identifiability improvement. Bootstrap and mixture fitting themselves are
established. Depends on V08 where axial symmetry is used and V01-style specimen-level provenance.

**Tests and validation.** Protocol B/D; uniform texture, a known mixture, overlapping components,
grain-size bias and sharp fibres. Test conservation of fractions and density normalization.
Refit the whole pipeline within a bootstrap replicate where its uncertainty is claimed.
Report sensitivity separately from confidence: regularization does not recover unobserved
information, and ghost-corrected solutions are not uniquely measured ODFs.

### V10 — Joint EBSD–diffraction texture consistency

**Scope and use cases.** Fit or compare one declared texture model against EBSD orientations
and diffraction pole figures through their separate sampling, correction and calibration
models. Start with a consistency report before a joint optimizer. Explicitly represent surface
versus bulk sampling and reject an unjustified “same population” assumption. Use for rolled
sheet and zirconium tubes measured by more than one technique.

**Impact, difficulty, novelty.** Impact 5, XL, R. Candidate contribution: identify when a
discrepancy is sampling, frame calibration or model inadequacy, and quantify the evidence.
Depends on V09, calibrated measurements and paired specimens; V11 supplies powder integration.

**Tests and validation.** Protocol D; inject a frame error, pole-dependent defocus, incomplete
coverage and a surface texture different from the bulk. Compare separate fits, naive pooling
and the proposed consistency-aware method on held-out reflections/regions. Success is more
accurate prediction and correctly flagged incompatibility, not merely a lower combined residual.

### V11 — Texture-aware multiphase powder quantification

**Scope and use cases.** Extend the existing single-phase `refine_rietveld` into a bounded
multi-phase profile model: calibrated scales, physically defined mass/volume fractions,
shared instrument parameters, phase-specific texture and constraints. Retain fixed structural
coordinates initially. Add correlation/identifiability diagnostics and an unknown-phase
remainder warning. Use for retained austenite, transformation products and phase mixtures.

**Impact, difficulty, novelty.** Impact 5, L, I/R. Existing GSAS-II/MAUD comparisons are required;
novelty could be independently constrained texture with better fraction uncertainty. Merely
summing profiles is E. Depends on V09 for uncertain texture and Protocol E data access.

**Tests and validation.** Pure-phase limits, nonnegative fractions summing to one under the
declared phase-completeness model, scale/unit invariance, overlapping peaks and a missing phase.
Compare certified or gravimetrically prepared mixtures using matched structural models.
Convert mass and volume fractions with declared density information; an arbitrary intensity
scale fraction is not a mass fraction. Keep structural-coordinate/occupancy refinement out of v1.

### V12 — Calibrated residual-stress analysis

**Scope and use cases.** Build a multi-tilt, multi-reflection lattice-strain workflow using
existing peak/lattice tools, instrument corrections and elastic tensors. Begin with a clearly
declared plane-stress model and stress-free spacing; report uncertainty, excluded reflections
and sensitivity to the elastic model. Use for coatings, formed sheet and residual stress after
processing. A single peak shift must not determine a full stress tensor.

**Impact, difficulty, novelty.** Impact 4, L, I. Core methods are established; publication is
more likely an application or validated texture-dependent extension. Depends on V11's
instrument discipline, measured tilt geometry and independent mechanical loading/reference.

**Tests and validation.** Zero stress, known imposed elastic strain, rotated loading,
insufficient tilts and correlated zero-shift errors. Recover held-out applied loads within
combined measurement uncertainty; compare a matched external stress workflow. Texture-aware
diffraction elastic constants require a validated grain-interaction model, not automatic
substitution of a bulk elastic modulus.

### V13 — Out-of-core EBSD with correct global topology

**Scope and use cases.** Introduce a semantic array-storage boundary and one supported chunked
backend, then streamed scan conversion, halo-aware local metrics and global grain/parent graph
reconciliation. Record chunking, units and provenance. Use for scans that exceed RAM.
First deliver read/subset/KAM; global segmentation and reconstruction are separate increments.

**Impact, difficulty, novelty.** Impact 5, L, E/I. Practical adoption value is high; a methods
paper needs a demonstrably new graph/parallel algorithm. Existing bounded numerical kernels
remain useful. Depends on V22 for replay metadata and V23 for memory benchmarks.

**Tests and validation.** Protocol F: compare in-memory and chunked results across boundaries,
holes, phase interfaces and long grains spanning many chunks. Test chunk-order and chunk-size
invariance, interrupted writes and storage-version errors. Measure process peak memory including
native arrays, not just Python allocator traces. Document exact versus tolerance-equivalent
outputs; do not promise constant memory for a global graph without a demonstrated bound.

### V14 — Hexagonal/irregular-grid curvature with uncertainty

**Scope and use cases.** Extend existing square-grid curvature/GND using local geometry-aware
derivatives on hex and irregular neighborhoods, with rank/conditioning checks, phase-boundary
masks and noise propagation. Use for deformed HCP alloys and scans with missing pixels.

**Impact, difficulty, novelty.** Impact 4, M–L, I. A calibrated spatial estimator could support a
methods contribution; changing a finite-difference stencil alone has limited novelty.
Depends on current coordinate graphs and the uncertainty contract from V02/V09.

**Tests and validation.** Constant and analytic linear rotation fields sampled on square, hex
and perturbed grids, rigid specimen rotation, holes and step-size sweeps. Compare shared
identifiable quantities with a matched MTEX dislocation workflow and measured repeat scans.
Unmeasured derivative components remain unmeasured; the method must not reinterpret a 2-D GND
lower bound as total dislocation density or statistically stored dislocation content.

### V15 — Additional vendor and 4D-STEM bridges

**Scope and use cases.** Add Oxford H5OINA/Bruker readers when licensed examples and documented
semantics are available. Separately import pyxem/4D-STEM orientation, phase and calibration
results into canonical PyTex types. Choose one bridge by real user demand; a generic `.h5`
extension is not a file-format specification. Use for instrument interoperability.

**Impact, difficulty, novelty.** Impact 4, M per bridge, E/I. Adapter novelty is low; value comes
from joining measurements to OR analysis without convention loss. Depends on V13 for large
datasets and V22 for transport. Do not rebuild pyxem's raw-pattern engine.

**Tests and validation.** Vendor golden examples, independent orientation/frame fiducials,
non-square detector pixels, invalid units, phase IDs and scan ordering; round-trip semantic
equality where the target format supports it. Validate optional/missing metadata explicitly.
Use open, pinned upstream examples, with version and redistribution rights recorded; importer
agreement on a synthetic array alone does not establish vendor-frame correctness.

### V16 — Image-to-SAED and ring-pattern metrology

**Scope and use cases.** Convert detector images to spots with subpixel positions and covariance,
then call existing solvers. Add beam-stop masks, saturation/background handling and ellipse-aware
ring extraction for polycrystals. Preserve rejected features and manual corrections. Use for
precipitate identification and nanocrystalline films.

**Impact, difficulty, novelty.** Impact 4, L, I. Detection/ring indexing are established; a paper
needs demonstrated metrological reliability or a distinct OR-discrimination workflow.
Depends on calibrated detector models and an annotated measured-image corpus.

**Tests and validation.** Protocol G: planted spots with independent noise/background generators,
overlap, saturation, beam stops and elliptical distortion; blinded expert-picked images.
Score localization, missed/false detections, calibration bias and final indexing accuracy.
Do not report covariance from a saturated or unresolved peak as a trustworthy Gaussian error.

### V17 — Information-guided TEM measurement planning

**Scope and use cases.** Extend existing ambiguity analysis, discriminating-zone sweeps and
stage routing to rank the next feasible measurement by expected separation of remaining
phase/OR/orientation hypotheses. Include calibration uncertainty, acquisition time and an
explicit dose proxy. Use for selecting tilts or diffraction conditions with limited beam time.
The first version advises an operator; it does not control microscope hardware.

**Impact, difficulty, novelty.** Impact 5, XL, R. Strong research potential if fewer measurements
achieve the same reliable identification than geometric nearest-zone, random and expert
baselines. Depends on V02/V16/V18 and a microscope partner; an information score without
prospective experiment evidence is insufficient.

**Tests and validation.** Symmetry-indistinguishable hypotheses yield no claimed discrimination;
unreachable targets are never recommended. Protocol G tests noisy calibration and missing spots;
then prospective operator runs compare time, dose proxy, success and incorrect certainty.
Report model mismatch and inability to distinguish candidates, including equal-prior sensitivity.

### V18 — Electron-model validation and physical ensembles

**Scope and use cases.** Establish matched-input comparisons for PyTex kinematic/Bloch and
abTEM multislice paths; add explicit thickness, orientation, frozen-phonon and detector-response
ensembles through the adapter where physically supported. Separate epistemic model discrepancy
from stochastic variation. Use for measured HREM/CBED interpretation and quantitative simulations.

**Impact, difficulty, novelty.** Impact 4, L, I/R. abTEM already supplies multislice; adapter
plumbing is E. A calibrated applicability map or cross-model error estimator is a plausible
methods result. Depends on reproducible potentials, sampling and experimental calibration.

**Tests and validation.** Protocol H: vacuum/norm checks, thin/weak-scattering limits, sampling,
slice-thickness and beam-set convergence; compare physical arrays with aligned units and
normalization, not visually similar screenshots. Validate measured thickness/defocus series.
Agreement in a thin limit must not be extrapolated to strongly dynamical specimens.

### V19 — 3-D and time-resolved transformation tracking

**Scope and use cases.** Register volumes or sequential maps, track parent/child identities,
variant birth and grain split/merge events, and attach confidence to temporal links. Start with
an offline two-time-point workflow before streaming or full 3-D reconstruction. Use for in-situ
transformation studies and serial-section EBSD.

**Impact, difficulty, novelty.** Impact 5, XL, R. The research opportunity is a transformation-aware
tracking model, not simply stacking maps. Depends on V03/V13/V15 and independently annotated
time/3-D data; budget registration as a scientific problem with spatial distortion.

**Tests and validation.** Known rigid/affine motion, grain growth, transformation, split/merge,
missing slices and occlusion. Compare orientation-only, overlap-only and combined tracking on
held-out sequences. Require correct event identity and calibrated ambiguity; persistent IDs
cannot be validated by comparing arbitrary integer labels alone.

### V20 — HCP mechanics and simulation-ready microstructures

**Scope and use cases.** Add explicitly normalized HCP basal/prismatic/pyramidal slip and twin
families, boundary slip-transmission descriptors, and a narrow DAMASK export/return adapter.
Carry phase, orientation, units and grain identity through the boundary. Use for Ti/Zr
deformation and reconstructed-microstructure property studies. Defer a native constitutive solver.

**Impact, difficulty, novelty.** Impact 4, L, I. Existing mechanics algorithms are not new;
publication needs a demonstrated structure–property result or uncertainty propagation from OR
reconstruction into independently measured response. Depends on V03/V09 and calibrated material data.

**Tests and validation.** Analytic Schmid bounds under the declared stress normalization,
orthogonal slip direction/plane normals, exact Burgers-vector magnitudes, symmetry duplicates
and frame invariance. Compare MTEX descriptors and a single-crystal DAMASK case; then validate
aggregate response on a separate loading path. Do not treat geometric transmission as a full
prediction of slip activation or fracture.

### V21 — Neutron and time-of-flight texture workflow

**Scope and use cases.** Extend existing radiation/scattering and frame models into one
instrument-specific TOF/NeXus adapter, bank calibration and texture refinement workflow.
Use for bulk texture and phase fractions in thick components. Neutron intent in the mission
does not mean a complete instrument reduction workflow exists today.

**Impact, difficulty, novelty.** Impact 4, XL, I/R. Requires a beamline partner and V09–V11.
Established TOF reduction should be reused through adapters. A publication might establish
cross-instrument consistency or an identifiable joint refinement, not the first neutron texture code.

**Tests and validation.** Flight-path/time-zero conversion, bank frame rotations, wrong-bank
metadata and isotope/scattering conventions. Use a bank-calibration standard and the published
quartz coordinate-validation scenario as references, followed by a paired laboratory-XRD/neutron
specimen. A neutron instrument operator must approve the model and uncertainty budget.

### V22 — Replayable study bundles and assumption comparison

**Scope and use cases.** Extend existing manifests, dossier, result exports and logbook into
one study bundle: input hashes, calibrated frames, parameters, masks, versions, result contracts,
references and regenerable figure recipes. Add side-by-side results under different assumptions.
Use for paper supplements, review, laboratory handoff and teaching exercises.

**Impact, difficulty, novelty.** Impact 5, M, I; support for a software paper rather than a new
numerical method. New schemas must compose the existing contracts. Start with OR and then texture;
do not invent a second analysis engine or place scientific calculations in JavaScript.

**Tests and validation.** Replay in a clean environment without internet, reject changed hashes,
round-trip scientific meaning and masks, exercise schema migrations and compare core outputs.
Browser tests must cover imported studies, changed assumptions and incomplete evidence.
Large/raw data may be external with checksum and acquisition instructions; generated study
outputs stay outside repository history unless they are named canonical baselines.

### V23 — Sustainable APIs and performance evidence

**Scope and use cases.** Split large implementations behind compatible facades, curate domain
entry points, complete report `describe()`/JSON coverage and add versioned timing/peak-memory
benchmarks in a controlled CI lane. Use for contributor onboarding, safe feature growth and
realistic throughput planning. Preserve deprecation policy and existing imports.

**Impact, difficulty, novelty.** Impact 4, M initially, E. This supports every paper and release
but is not a standalone scientific innovation. Sequence small refactors by the feature touching
them; do not undertake a repo-wide rewrite before the research program.

**Tests and validation.** Import compatibility and isolation, schema reconstruction and
description/number agreement; benchmark fixed workloads after warm-up on recorded hardware.
Compare medians and process memory over repetitions, with a predeclared noise allowance.
Keep CPU workbench tests separate from heavy numerical runs when diagnosing timing effects.
Improve the coverage ratchet without lowering the measured baseline or hiding skips.

### V24 — Anisotropic displacement and quantitative crystal scenes

**Scope and use cases.** Extend `AtomicSite`/CIF boundaries with anisotropic displacement tensors,
then thermal ellipsoids, coordination-angle and polyhedron-distortion reports. Preserve the
isotropic representation and distinguish occupancy from displacement. Use for structure-quality
inspection and teaching. Electronic-density volumes remain outside the current mandate.

**Impact, difficulty, novelty.** Impact 3, M, E/I. Established VESTA-class functionality, low
standalone novelty; prioritize only if it serves an active interface or scattering study.
Depends on core tensor conventions and adapter validation before renderer changes.

**Tests and validation.** Positive-semidefinite tensors, isotropic sphere limit, probability-level
scaling, covariance under frame changes and CIF round-trip with documented units. Compare a
licensed anisotropic-CIF example and VESTA semiaxes/angles quantitatively. Structural assertions
and measurements take precedence over runtime SVG byte baselines.

## 5. Publication programs with go/no-go gates

### Paper A — A validated instrument for OR and interface analysis

**Likely form:** software/methods paper; the Computer Programs scope of *Journal of Applied
Crystallography* is a natural precedent, not an acceptance prediction. Minimum scope: existing
OR engine/dossier plus V01, V04 and V22. V02 strengthens it; V06 is not required.

**Research question:** can a convention-explicit workflow reproduce measured OR analyses and
interface candidates reliably across independent datasets and tools?

**Evidence:** one steel and one Burgers-system measured study, an external comparison, an
independent interface case, clear uncertainty limits, and a replayable versioned release.
Compare PyTex with an explicitly scripted multi-tool workflow. Measure scientific discrepancies,
successful reproduction and, if claiming usability, controlled task completion rather than
anecdotal convenience. Include failure/ambiguous examples in the paper.

**Go/no-go:** no broad “first integrated tool” claim; no parent-reconstruction performance claim
without Protocol A. If measured truth is missing, narrow to validated correspondence/interface
software or delay submission. A collection of rendered pictures is not the main result.

### Paper B — Reliable OR and parent inference under ambiguity

**Likely form:** methods paper built around V02/V03, with V01's benchmark. Compare weighted
PyTex, robust-only fitting, uncertainty-aware fitting, current graph reconstruction and the
matched MTEX method. Ablate priors, parent anchors, spatial consistency and resampling choice.

**Evidence:** calibration curves, false merges/splits, retained coverage versus assignment error,
angular error, scale/memory costs, and failure regimes on multiple independent specimens.
Split train/tuning/test by specimen; connected grains from one map cannot serve as independent
training and test populations. Report out-of-distribution material or instrument behavior.

**Go/no-go:** advance only if the method provides calibrated or demonstrably more reliable
decisions beyond the baselines. If the new method only repackages known inference, publish the
benchmark or software contribution honestly instead of inflating the algorithmic claim.

### Paper C — Texture-informed phase quantification across measurements

**Likely form:** methods/application paper joining V09–V11. Question: does independently
measured texture improve phase-fraction estimates, and can the model detect incompatible
surface and bulk populations?

**Evidence:** mixtures with independently known fractions, paired EBSD/PF/powder measurements,
instrument standards, texture-free and March–Dollase baselines, matched GSAS-II/MAUD analyses,
and held-out reflections/specimens. Ablate each modality and calibration component.

**Go/no-go:** phase fractions must improve against external truth, not simply against synthetic
patterns generated by the same code. Uncertainty must expose scale–texture confounding.
If all modalities do not sample the same population, report inconsistency rather than forcing
a common ODF. Do not bundle this paper into Paper A merely to increase feature count.

### Paper D — Choosing the next TEM observation

**Likely form:** higher-risk methods paper around V17, after V16/V18. Compare uncertainty-aware
planning against nearest reachable zone, random feasible selection and operator judgment.
Perform retrospective replay first, then a prospective study with fixed stopping rules.

**Evidence:** measurements/time/dose proxy to reach a specified error rate, calibration robustness,
unreachable/indistinguishable cases and blinded final identification. A microscope partnership
is essential. Stop at a planning prototype if prospective validation cannot be obtained.

## 6. Concrete testing and validation protocols

These are **proposed acceptance criteria**, not measurements of current performance or
universal physical tolerances. Apply the
[benchmark and tolerance standard](../standards/benchmark_and_tolerance_governance.md).
For numerical identities, prefer condition-aware tolerances tied to units, resolution and
floating-point propagation. For experimental comparisons, combine reference uncertainty,
measurement uncertainty, registration error and sampling variability before choosing a threshold.

| Protocol | Cases and independent reference | Proposed acceptance and reporting |
| --- | --- | --- |
| **A — OR reconstruction** | Exact synthetic KS/NW/GT and Burgers populations; 1, 10 and 100 parents; missing variants, twins, holes and 0/0.5/1/2/5 degree noise; measured partially transformed steel with independently observed parent orientations; Ti/Zr study with retained parent or independent reference; version-pinned MTEX output. | Exact noiseless identifiable cases recover the partition up to label permutation. For a pilot of well-resolved measured grains, target at least 95% of accepted parent orientations within 1 degree of reference and adjusted Rand index at least 0.95 on the reference-covered region. Freeze or justify those targets from instrument uncertainty before testing. Also report assigned fraction, false merges, false splits, boundary-location error and excluded regions; never raise apparent accuracy by hiding abstentions. |
| **B — Statistical calibration** | At least 500 independent synthetic datasets per declared regime with known OR/texture truth, clustered observations, heterogeneous precision and outliers; separate measured repeat scans/specimens. Use an independent data generator for at least one campaign. | A nominal 95% interval/region should show 92–98% empirical coverage in the declared benchmark regime, with a binomial interval reported for that coverage estimate. The band is a planning target, not a guarantee. Report interval width, bias, error after abstention and random seeds. Calibrate uncertainty at the actual independent unit, typically grain/block/specimen, not pixel. |
| **C — Interfaces/habit planes** | Identical cells, analytically commensurate 2:3 in-plane lengths, an oblique low-symmetry plane, a published lattice-matching case, and measured nonparallel traces or serial-section normals. Compare pymatgen matching only with the same search limits and metric. | Exact cases give zero mismatch within a justified numerical tolerance; symmetry/basis relabeling preserves invariant rankings. Record match area versus strain instead of hiding the tradeoff in one score. Predicted measured-plane uncertainty regions must be compatible with independent normals; single-trace data must retain the non-identifiable family. |
| **D — Texture/multimodal** | Uniform orientation distribution, ideal axial fibre, an independently generated known mixture; missing polar caps, defocus errors, frame rotations and conflicting surface/bulk textures; measured rolled sheet or tube with repeated EBSD and PF scans. | Density normalization and Kearns sum rule hold; noiseless fixed-component mixtures recover known fractions. For measurement, use Protocol B intervals, held-out pole-figure prediction and explicit residual structure. Quantify ghost/symmetry/width sensitivity separately. Same-code forward/inverse recovery is a consistency test, not external validation. |
| **E — Powder quantification/stress** | Pure end members and known 20/80 and 50/50 mixtures by mass; similar lattices, heavy overlap, unknown impurity, strong texture; line-position/shape calibration from [NIST SRM 660c](https://www.nist.gov/publications/certification-standard-reference-material-660c-powder-diffraction); independent applied-load stress case. | Initial mixture target: absolute fraction error at most 2 percentage points for a preregistered well-conditioned mixture set, with Protocol B interval coverage; report failures outside that regime. A line-profile standard validates calibration, not phase fractions. Stress estimates must agree within a predeclared combined uncertainty budget. Low profile residual alone is insufficient. |
| **F — Scale/storage** | Same maps at 100 thousand, 1 million and 10 million points; many chunk sizes, long grains crossing chunks, sparse/dense graphs, phase boundaries and interruptions; in-memory baseline and independently counted partitions. | Exact discrete labels after equivalence mapping where algorithm semantics permit; numerical metrics within the in-memory method's tolerance. Target the declared 10-million-point streamed local-metric workflow within 4 GiB process peak memory on a recorded host; separately measure global graph memory. Record median time, variability, peak resident memory, thread count and versions. No portable runtime promise until measured. |
| **G — Image metrology/planning** | Independently generated calibrated spots and rings, tilt/calibration perturbations, saturation, beam stops, missing/spurious spots; independently annotated detector images; held-out microscope sessions. | Pilot target: median localization error below 0.2 pixel for isolated unsaturated peaks in a fixed, declared noise/width regime. Report detection precision/recall and downstream wrong-phase rate, not localization alone. Planner target: at least 20% fewer observations at matched identification error than nearest-zone planning on held-out cases; confirm prospectively before a benefit claim. |
| **H — Electron simulation** | Vacuum, weak-phase and two-beam checks; Si plus an ionic or higher-scattering structure; thickness/defocus/orientation sweeps; matching abTEM potential/CTF and measured series with independent thickness estimates. | Converge detector observables to a predeclared budget before comparing models. Pilot target: under 1% relative integrated-intensity change on further sampling/beam/slice refinement for nonzero selected observables; use absolute tolerances near zero. Measure model error against experiment separately. A converged model can still be physically inappropriate. |

### 6.1 Dataset governance and acquisition backlog

| Needed resource | Scientific owner to appoint | Truth/metadata required | Acquisition gate |
| --- | --- | --- | --- |
| Steel and Ti/Zr parent-reference EBSD | OR lead + experimental collaborator | Parent truth source, phase structures, coordinate conventions, step size, indexing quality and segmentation policy | Permission to redistribute a canonical subset, or a stable public acquisition route with license and checksum. Expert reconstruction alone is an annotation, not independent parent truth. |
| Paired EBSD/PF/powder specimens | Texture lead + diffraction operator | Same-specimen identity, surface/bulk sampling, specimen axes, calibration and repeat measurements | Confirm comparable populations before a joint-fit study. |
| Powder standards and known mixtures | Diffraction lead | Standard certificate/version, mixture preparation uncertainty, wavelength, geometry and raw counts | Obtain experiment access; a certificate without the local measured scan does not validate local calibration. |
| Detector images and TEM series | TEM lead + microscope operator | Raw images, pixel response/geometry, voltage, tilt, thickness/defocus reference and annotations | Agree the prospective protocol and retain raw-data provenance. |
| 3-D/time and TOF data | Future partner-led programs | Spatial/time registration, event truth or bank/flight-path calibration | Defer implementation commitments until data and an independent evaluator exist. |

Appointments are proposed roles, not commitments by named people. Nothing in this document
authorizes contacting collaborators or obtaining restricted data. Include original file hashes,
license, source DOI/URL, instrument/software versions, transformations, exclusions and truth type
in each benchmark manifest. Synthetic and measured fixtures must remain clearly distinguishable.
Commit only sources and small referenced canonical fixtures; large data, plots, logs and generated
reports follow the repository's source-only rule.

### 6.2 Common feature exit gate

A stable feature needs the existing core semantic types or a reviewed extension, constructor
invariants, vectorized/bounded batch behavior, `describe()` and reconstruction-grade JSON where
applicable. It also needs deterministic/property tests, failure/ambiguity cases, independent
known answers, an executable worked example, theory and workflow documentation, a validation
ledger row and an appropriate benchmark. New symbols go into the central registry; geometry
requires a canonical SVG following the visualization guide. UI forms use the shared parameter
width/row/symbol declarations and show all mandatory inputs together.

Before a feature lands, run Ruff, strict mypy and the base suite with nondecreasing measured
coverage, no new warnings, repository integrity and the zero-warning Sphinx build. Run browser
journeys for any changed application behavior. External data and comparator runs have their own
versioned evidence gates; a green unit suite cannot waive them. Record each completed increment
and exact next action in the ledger, explicitly stage its paths, then commit and push to main.

## 7. Sequencing and scope control

Use milestones rather than calendar promises. With one scientific developer, reserve roughly
one-third of available effort for validation, documentation and maintainability. With additional
contributors, independent infrastructure work can proceed alongside science, but the gates and
dependencies remain. Do not sum all 24 estimates into a promised release schedule.

| Stage | Deliverables | Exit decision |
| --- | --- | --- |
| **0 — Evidence and contracts** | V01 dataset/comparator protocol and access plan; V22 minimal replay for an existing OR case; V23 reproducible benchmark/coverage baseline. | An independent reference can actually be obtained; metrics and hold-out splits are frozen. Start logistics now; no scientific novelty claim yet. |
| **1 — Flagship reliability** | V02 bounded robust/uncertainty fit, then V03 reconstruction pilot. Begin V13 streamed import/local metrics if map size is a demonstrated constraint. | Held-out evidence supports calibrated decisions. If not, retain experimental status and publish the failure boundary. |
| **2 — Interfaces and first paper** | V04 measured-scale interface report, dossier/replay integration, Paper A evaluation. V05 only when trace/3-D truth is available. | A reviewer can reproduce two measured material studies and an independent interface case. Do not wait for PTMC or a native mechanics solver. |
| **3 — Texture and diffraction** | V08, V09, bounded V11; then V10 consistency and joint inference. V12 only with a calibrated stress case. | Independent mixture/texture truth supports Paper C; otherwise retain separate analyses and diagnostics. |
| **4 — Partner-led research** | Four-week feasibility studies for V06, V17, V19 or V21, one at a time. V16/V18 support the TEM path. | A distinct question, suitable data, comparator and feasible pilot justify a funded implementation. Stop/defer otherwise. |

**First five implementable work packages:** (1) V01's benchmark specification and external-run
handoff; (2) V22's OR replay bundle; (3) V02's robust-fit prototype and coverage experiment;
(4) V04's identical/rational-interface kernel and report; (5) V13's streamed-read/KAM prototype.
The OR scientific track has priority over unrelated breadth when resources compete. V08 is a
good bounded alternative while measured-data access is pending, but it does not substitute for
the flagship's evidence gate.

Each approved work package should name a scientific owner and reviewer, its baseline commit,
feature IDs, out-of-scope items, input/output contract, evidence protocol, effort ceiling and
stop condition. Keep implementation status and publication readiness separate. Review this
portfolio after two development cycles, recording what shipped and retiring stale claims.

### Explicit deferrals

- A full native Rietveld structure-refinement engine, crystal-plasticity solver, electronic-density
  viewer or general raw-4D-STEM framework would compete with the flagship and mature external
  tools. Prefer narrow adapters until a demonstrated scientific need justifies otherwise.
- Generic AI labeling, opaque “confidence” scores and language-model-generated scientific
  interpretations are not substitutes for calibrated inference and deterministic `describe()`.
- GPU support is a response to measured bottlenecks, not a prerequisite for every new algorithm.
- Automatic microscope control, proprietary database redistribution and cloud collaboration
  infrastructure are separate projects, outside the scope of this vision's first programs.

## 8. Measures of progress

Track independent measured cases and material families, reproducible external comparisons,
calibrated uncertainty coverage, false confident conclusions, complete replay success,
time/peak memory at fixed scientific accuracy, and independently repeated analyses. Track paper
readiness by closed evidence gates, not manuscript length, feature count or plotted examples.

The first major success is a collaborator reproducing a measured OR/interface study from its
declared inputs and obtaining both the same result and the same justified limits. The next is
a demonstrated research improvement that survives held-out data and an honest prior-art comparison.

## 9. Governing references and implementation starting points

The external sources linked beside each claim are informative prior art and candidate validation
references. Numerical implementations must additionally pin the applicable normative equation,
units and convention in their theory note before coding. This vision introduces no new formula
or crystallographic convention.

- [Canonical data model](../architecture/canonical_data_model.md),
  [OR analysis foundation](../architecture/orientation_relationship_analysis_foundation.md),
  [multimodal foundation](../architecture/multimodal_characterization_foundation.md).
- {doc}`Testing strategy </validation/testing_strategy>`,
  {doc}`MTEX parity matrix </validation/mtex_parity_matrix>`,
  {doc}`phase-transformation validation </validation/phase_transformation_validation_matrix>`,
  {doc}`diffraction validation </validation/diffraction_validation_matrix>`,
  {doc}`VESTA parity </validation/vesta_parity_matrix>`.
- [Documentation architecture](../standards/documentation_architecture.md),
  [executable examples](../standards/executable_examples.md),
  [notation](../standards/notation_and_conventions.md),
  [symbol registry](../standards/terminology_and_symbol_registry.md),
  [visualization style](../standards/visualization_style_guide.md),
  [data contracts](../standards/data_contracts_and_manifests.md),
  [reference canon](../standards/reference_canon.md).
- [Previous September review](https://github.com/kvmani/pytex/blob/e497c8a/docs/development/repository_review_2026_09.md),
  [world-class roadmap](world_class_feature_roadmap.md),
  [MTEX/EBSD roadmap](https://github.com/kvmani/pytex/blob/e497c8a/docs/roadmap/mtex_parity_and_ebsd_feature_roadmap.md).
