# Implementation Roadmap

This roadmap is split into two layers:

1. the immediate development program for the next phase
2. the longer-horizon capability ladder that keeps the project aligned with the mission

The immediate phase is intentionally corrective. The repository already has substantial scientific
implementation depth, but the next phase must restore reproducibility guarantees, clean up the
rendered documentation surface, and raise the quality of validation evidence before new features
are added aggressively.

## Immediate Development Roadmap

### Phase 0: Trustworthiness And Documentation Hardening

Priority: highest

Goals:

- restore repository integrity and fixture pinning
- remove broken Sphinx/MyST link hygiene from the canonical docs
- turn validation pages into evidence-backed scientific references rather than policy summaries
- make the public docs easier to use for scientific understanding and implementation onboarding
- prove optional interoperability claims with live integration tests where dependencies are available

This phase is the default governing program for the next development cycle. New feature work should
only land when it is required to complete one of the workstreams below or when it clearly removes a
blocker inside this phase.

#### Phase 0 Operating Rules

- Do not add new stable public surface area unless it is required for a hardening task in this
  phase.
- Treat `pytest -q`, repository integrity checks, and Sphinx HTML builds as release-blocking gates
  during this phase.
- If a task changes scientific meaning, update code, tests, Sphinx docs, validation notes, and
  theory assets in the same change.
- Prefer changes that reduce ambiguity in claims, validation posture, or user interpretation over
  changes that merely broaden surface area.

#### Phase 0 Milestones

| Milestone | Purpose | Required outcome |
| --- | --- | --- |
| M0.1 | restore trust in pinned scientific assets | fixture corpus and integrity checks are green |
| M0.2 | make the docs build a reliable scientific surface | Sphinx build is warning-free for roadmap-targeted pages |
| M0.3 | upgrade validation from policy to evidence | high-value validation pages contain worked evidence blocks |
| M0.4 | align public claims with executable reality | optional adapters and workflow claims are backed by tests |
| M0.5 | make the learning path coherent | quickstart, concepts, notebooks, and validation form one teaching route |

#### Phase 0 Release Gates

The immediate roadmap is not complete until all of the following are true:

- `pytest -q` passes in the default repository environment
- repository integrity checks pass from a clean checkout
- the Sphinx HTML build completes without unresolved cross-reference warnings for the targeted docs
- the roadmap-targeted validation pages show source-backed numerical evidence rather than only
  policy text
- the public docs surface explains current limitations where validation or physics is still
  foundational
- adapter-facing claims in docs are no stronger than the available integration coverage

#### Workstream 0.1: Reproducibility Backbone

Tasks:

- refresh the phase fixture catalog hashes and keep the pinned CIF corpus synchronized with the
  actual repository contents
- keep the repo-integrity script as a hard gate for fixture, metadata, and manifest drift
- ensure every bundled phase fixture is referenced by tests, manifests, and at least one user-facing
  doc or workflow path

Primary repository surfaces:

- `fixtures/phases/catalog.json`
- `fixtures/phases/*/phase.cif`
- `fixtures/phases/*/metadata.json`
- `benchmarks/structure_import/`
- `scripts/check_repo_integrity.py`
- `scripts/regenerate_phase_fixture_catalog_hashes.py`

Verification gates:

- fixture hash tests pass
- repo-integrity script passes
- structure-import manifests refer only to pinned, present, and documented assets

Implementation notes:

- refresh digests only after confirming the underlying fixture content is intentional
- if fixture contents changed, update the corresponding metadata, benchmark references, and docs in
  the same change
- do not weaken integrity assertions to make the suite pass; fix the asset graph instead

Exit criteria:

- `pytest -q` passes with no integrity failures
- the phase fixture catalog and metadata digests are reproducible from a clean checkout
- repository integrity checks fail loudly when fixture contents drift

#### Workstream 0.2: Documentation Link Hygiene

Tasks:

- fix unresolved MyST/Sphinx cross-references in canonical architecture and standards docs
- prefer navigable repository-relative Markdown links or proper Sphinx `:doc:` links over bare file
  paths in prose
- expose the roadmap itself from the public Sphinx entry points so contributors can find the active
  development plan without hunting through source files

Primary repository surfaces:

- `docs/architecture/*.md`
- `docs/standards/*.md`
- `docs/site/index.md`
- `docs/site/reference/canonical_docs.md`
- `docs/site/validation/index.md`

Verification gates:

- `python -m sphinx -b html docs/site docs/_build/html` completes without unresolved cross-reference
  warnings for targeted pages
- roadmap, validation, and architecture pages are mutually reachable from the rendered HTML

Implementation notes:

- fix the canonical source documents rather than silencing warnings in build configuration
- convert source-doc links that currently point at `../site/...` or non-rendered relative targets
  into links that remain valid from the rendered site
- use this workstream to enforce the documentation architecture rule that required-reading docs are
  linked, not merely named

Exit criteria:

- the Sphinx HTML build completes without unresolved cross-reference warnings
- the public site exposes the roadmap as a discoverable navigation target
- readers can move from overview, concept, validation, and roadmap pages without dead links

#### Workstream 0.3: Validation Evidence And User Guidance

Tasks:

- expand documented test cases so major validation pages show source, formula, worked example,
  expected result, current output, and tolerance
- add clearer “what this object means” and “why this convention exists” introductions to the core
  concept pages
- revise quickstart and API guide pages so they teach the canonical scientific path instead of only
  listing symbols
- add explicit limitations and approximation callouts to validation and workflow pages where the
  implementation is geometric, kinematic, or pedagogical rather than physically complete

Primary repository surfaces:

- `docs/site/validation/automated_test_cases.md`
- `docs/site/validation/index.md`
- `docs/testing/automated_test_cases.md`
- `docs/site/tutorials/quickstart.md`
- `docs/site/api/index.md`
- the highest-value concept and workflow pages referenced by those docs

Verification gates:

- validation pages include at least one worked numerical evidence block for each targeted major area
- quickstart and API guide can be read in isolation by a technically competent new user
- targeted pages state current approximation boundaries explicitly

Implementation notes:

- start with the highest-risk scientific surfaces: structure import, diffraction, texture
  reconstruction, and frame/convention handling
- each evidence block should name the code surface, governing source, expected value, current code
  output, tolerance, and interpretation
- when a page explains a scientific object, lead with meaning and usage before exhaustive API lists

Exit criteria:

- a new contributor can follow the HTML docs from concept to workflow to validation without needing
  hidden repository context
- validation pages answer “what was checked?” and “what does the result mean?” directly
- onboarding pages explain frames, conventions, and failure modes before advanced API catalogs

#### Workstream 0.4: Optional Interoperability Reality Check

Tasks:

- add live integration coverage for optional adapters where the dependency is available in CI or a
  dedicated environment
- keep adapter docs honest about what is normalized, what is validated, and what is still
  dependency-limited
- verify EBSD and structure-import bridge paths against representative external payloads rather than
  only against object construction

Primary repository surfaces:

- `tests/integration/`
- `tests/unit/test_plotting_and_adapters.py`
- `tests/unit/test_orix_miller_adapter.py`
- adapter docs and workflow pages that currently describe interoperability

Verification gates:

- each claimed stable adapter path has at least one executable integration test or an explicit
  documented limitation
- adapter docs distinguish optional runtime dependency from canonical PyTex semantics

Implementation notes:

- prefer narrow, high-signal integration tests over broad but brittle environment-dependent suites
- if a dependency cannot be made reliable in default CI, keep the doc claim conservative and note
  the test environment required
- use representative external payloads where they validate boundary normalization rather than only
  constructor round-trips

Exit criteria:

- optional adapter claims are backed by executable integration tests
- adapter documentation distinguishes canonical PyTex semantics from source-library semantics
- dependency-limited paths are called out explicitly instead of implied by stable API exports

#### Workstream 0.5: Teaching-Surface Consolidation

Tasks:

- refresh the highest-value tutorials, especially the quickstart, core-model path, validation path,
  and notebook index
- make the notebook atlas mirror the documented learning sequence instead of acting as a parallel
  browse-only tree
- standardize visual and explanatory callouts so pedagogical examples are clearly labeled as such

Primary repository surfaces:

- `docs/site/tutorials/quickstart.md`
- `docs/site/tutorials/notebooks.md`
- `docs/site/concepts/*.md`
- priority notebook files under `docs/site/tutorials/notebooks/`

Verification gates:

- the recommended reading path matches the actual notebook and workflow sequence
- priority notebooks still smoke-execute after doc and example updates
- terminology and notation match the registry and concept pages

Implementation notes:

- keep the teaching path anchored to the implemented, validated surfaces rather than aspirational
  future features
- add “pedagogical approximation” callouts anywhere plots, intensities, or reconstructions are
  intentionally simplified
- prefer one coherent path from phase/structure to orientation/texture to diffraction/validation
  rather than many lightly explained entry points

Exit criteria:

- the HTML site presents one coherent learning path from core model to workflow to validation
- notebooks, concepts, workflows, and theory notes agree on terminology and notation
- the main tutorials explain the scientific model clearly enough for teaching and self-service use

### Phase 0 Suggested Execution Sequence

The workstreams above are coupled. The recommended implementation order is:

1. Workstream 0.1 first, because fixture integrity and repository trust must be restored before
   validation or docs can make stronger claims.
2. Workstream 0.2 second, because warning-heavy docs builds make every later documentation task
   harder to verify.
3. Workstream 0.3 third, because once the docs graph is healthy, validation pages can be upgraded
   into evidence-backed scientific guidance.
4. Workstream 0.4 fourth, because interoperability claims should be tightened only after the
   validation and docs posture is clearer.
5. Workstream 0.5 fifth, because teaching-surface consolidation should consume the corrected docs,
   validation, and interoperability posture rather than racing ahead of them.

### Phase 0 Deliverables By End Of Cycle

The repository should have all of the following by the end of this immediate phase:

- a passing fixture-integrity and repo-integrity baseline
- a warning-free or materially warning-reduced Sphinx build for the roadmap-targeted docs
- upgraded validation pages with explicit evidence blocks and limitations
- a public roadmap visible from the documentation entry points
- conservative, test-backed adapter claims
- a coherent teaching path across quickstart, concepts, workflows, and notebook atlas

## Capability Ladders

### Semantic Foundation

| Subsystem | Status | Notes |
| --- | --- | --- |
| Documentation and standards | implemented | Mission, standards, Sphinx, LaTeX, SVG, and validation doctrine are in place. |
| Canonical core model | implemented | Frames, transforms, symmetry, lattice, provenance, orientations, and reciprocal primitives exist. |
| Multimodal characterization doctrine | implemented | Shared cross-modality semantics are documented and the core acquisition/calibration/quality primitive family now exists. |
| Phase-transformation doctrine | implemented | Foundation docs exist and the first stable transformation primitive family now exists in the core model. |

### Validated Foundational Implementation

| Subsystem | Status | Notes |
| --- | --- | --- |
| Orientation and texture | implemented | Rotations, misorientation, PF/IPF, symmetry reduction, discrete ODF foundations, band-limited harmonic ODF reconstruction, and runtime contour or section plotting surfaces exist with tests and MTEX-backed ledgers. |
| EBSD map workflows | implemented | KAM, segmentation, GROD, boundaries, cleanup, graph aggregation, multiphase phase selection, and manifest-backed normalization now cover both regular-grid and graph-backed coordinate modes. |
| Diffraction foundations | implemented | Geometry, powder XRD generation, SAED spot workflows, and diffraction plotting now have pinned external-baseline cases for one powder XRD path and one SAED path; broader coverage still remains ahead. |
| Scientific visualization | implemented | Shared YAML style themes, runtime diffraction plotting, and VESTA-like crystal visualization now include structural plotting validation for the highest-value publication-facing surfaces without tracked runtime SVG baselines. |
| CIF and structure import | implemented | Core-model construction, space-group semantics, hash-pinned fixture integrity, and a manifest-backed fixture audit workflow now form the default reproducible structure-import baseline. |
| Multimodal acquisition core | foundational | Acquisition geometry, calibration, quality, scattering, experiment manifests, and workflow entry points now exist, but broader modality coverage remains ahead. |
| Phase-transformation foundation | foundational | Core transformation primitives, dedicated transformation manifests, and experimental parent-candidate scoring now exist, but literature-backed validation and richer algorithms remain ahead. |

### Research-Grade Algorithmic Expansion

| Subsystem | Status | Notes |
| --- | --- | --- |
| Exact orientation-space boundary catalogs | planned | Required for broader class-by-class parity claims. |
| Harmonic ODF inversion and richer reconstruction | foundational | A stable band-limited harmonic reconstruction surface now exists; broader experimental correction doctrine, higher-order validation breadth, and stronger external parity are still ahead. |
| Rich diffraction refinement and intensity models | planned | Current implementation is geometric and kinematic, not full physical modeling. |
| Phase transformation and parent reconstruction | planned | Stable transformation semantics now exist; algorithmic breadth and validation remain ahead. |

### Teaching-Grade Explanatory Surface

| Subsystem | Status | Notes |
| --- | --- | --- |
| Sphinx concepts and workflows | implemented | Public entry point is live and buildable, and the site now renders concepts, workflows, architecture, standards, validation, and benchmark context from one browsable surface. |
| Canonical LaTeX theory notes | implemented | Major foundation notes exist and are cross-linked from the site. |
| SVG geometry figures | implemented | Core orientation, diffraction, and EBSD figures exist. |
| Multimodal and transformation teaching notes | foundational | Architectural prose is now defined; broader workflow coverage remains ahead. |

## Immediate Default Development Order

Unless explicitly redirected, immediate development should follow this order:

1. engineering hygiene: restore fixture integrity, keep tests green, and make the docs build warning-free
2. documentation link hygiene: fix cross-references and expose the roadmap in the public site
3. validation evidence expansion: upgrade high-value validation pages and documented test cases
4. reproducibility and manifest hardening: tighten pinned assets, schemas, and workflow contracts
5. interoperability verification: add live adapter coverage where dependencies allow it
6. teaching-surface consolidation: improve quickstart, concepts, notebooks, and validation navigation
7. only then expand algorithmic breadth or new stable feature surface area

## Near-Term Non-Goals

The current roadmap intentionally deprioritizes the following until the hardening program above is further advanced:

- broad new modality expansion
- major stable adapter-surface expansion
- full physical diffraction intensity modeling before geometric and benchmark validation are stronger
- large new algorithmic breadth without theory, benchmark, and validation scaffolding

## References

### Normative

- [Engineering Governance](../standards/engineering_governance.md)
- [Reference Canon](../standards/reference_canon.md)

### Informative

- [Repository Review 2026 Foundation Audit](../architecture/repo_review_2026_foundation_audit.md)
