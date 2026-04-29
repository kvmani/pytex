# Implementation Roadmap

This roadmap is split into two layers:

1. the immediate development program for the next phase
2. the longer-horizon capability ladder that keeps the project aligned with the mission

The immediate phase is still intentionally corrective, but it is no longer a foundation-rescue
phase. Repository integrity, the default test suite, the rendered Sphinx surface, and the full
scientific lane are now green. The next phase must broaden evidence-backed validation before new
features are added aggressively.

## Immediate Development Roadmap

### Phase 1: Validation And Interoperability Hardening

Priority: highest

Goals:

- keep repository integrity and docs-build gates green while broadening evidence coverage
- formalize one base contributor lane and one full scientific lane
- turn validation pages into evidence-backed scientific references rather than only policy summaries
- expand structure-import, diffraction, and transformation evidence before broadening claims
- prove optional interoperability claims with live tests in the full scientific lane

This phase is the default governing program for the next development cycle. New feature work should
only land when it is required to complete one of the workstreams below or when it clearly removes a
blocker inside this phase.

#### Phase 1 Operating Rules

- Do not add new stable public surface area unless it is required for a hardening task in this
  phase.
- Treat the base lane and full scientific lane as first-class release gates during this phase.
- If a task changes scientific meaning, update code, tests, Sphinx docs, validation notes, and
  theory assets in the same change.
- Prefer changes that reduce ambiguity in claims, validation posture, or user interpretation over
  changes that merely broaden surface area.

#### Phase 1 Milestones

| Milestone | Purpose | Required outcome |
| --- | --- | --- |
| M1.1 | formalize the executable quality lanes | base and full scientific environments are documented and exercised in CI |
| M1.2 | keep the public docs phase-accurate | roadmap, README, and audit memo reflect current green gates and current risks |
| M1.3 | upgrade validation from policy to evidence | high-value validation pages contain worked evidence blocks for structure, diffraction, and transformation |
| M1.4 | align optional claims with executable reality | adapter and CIF-backed workflow claims are backed by full-lane tests |
| M1.5 | make the learning path follow the validated surface | quickstart, concepts, notebooks, and validation navigation reflect the current evidence path |

#### Phase 1 Release Gates

The current baseline has now closed the Phase 1 release gates below:

- `python scripts/check_repo_integrity.py` passes in both lanes
- `python -m pytest -q` passes in the base lane
- `python -m pytest -q -rs` passes in the full scientific lane without the previous `pymatgen`-gated skips
- `python -m ruff check .` and `python -m mypy src` are runnable in the documented contributor environment
- `python -m sphinx -b html docs/site docs/_build/html` completes without unresolved cross-reference warnings for the targeted docs
- the roadmap-targeted validation pages show source-backed numerical evidence rather than only policy text
- the public docs surface explains current limitations where validation or physics is still foundational
- adapter-facing claims in docs are no stronger than the available integration coverage

#### Workstream 1.1: Executable Quality Lanes

Tasks:

- standardize the base contributor lane around `.[dev,docs]`
- standardize the full scientific lane around `.[dev,docs,adapters]`
- keep CI, local-development docs, and README synchronized with those two lanes
- ensure the full scientific lane covers the optional structure and diffraction paths that depend on `pymatgen`

Primary repository surfaces:

- `pyproject.toml`
- `.github/workflows/ci.yml`
- `README.md`
- `docs/development/local_development.md`
- `docs/site/tutorials/installation_and_build.md`

Verification gates:

- the documented base lane can run integrity, lint, mypy, tests, and docs builds
- the documented full scientific lane can run the optional scientific coverage without skips
- CI reflects the same lane split as the docs

#### Workstream 1.2: Validation Evidence And Status Re-Baselining

Tasks:

- re-baseline the roadmap, README, and repo audit memo so they reflect the current green infrastructure gates
- expand the structure-import, diffraction, and phase-transformation validation pages with stronger evidence framing
- keep validation manifests, benchmark manifests, and documented evidence blocks synchronized
- continue tightening public limitation statements where coverage remains foundational

Primary repository surfaces:

- `docs/roadmap/implementation_roadmap.md`
- `docs/architecture/repo_review_2026_foundation_audit.md`
- `docs/testing/*.md`
- `docs/site/validation/index.md`

Verification gates:

- roadmap and audit docs no longer imply that integrity or docs-build recovery is still unfinished
- validation pages clearly separate hard claims, foundational claims, and current limitations
- the rendered validation section points users from theory to evidence to workflow limits

#### Workstream 1.3: External-Baseline Expansion

Tasks:

- prioritize external-baseline expansion for structure import, diffraction, and phase transformation
- keep the pinned fixture corpus, benchmark manifests, and validation manifests aligned
- use literature-tracked or open-tool-backed starter cases before claiming broader equivalence
- keep diffraction angle or shell geometry claims separate from stronger intensity-model claims

Primary repository surfaces:

- `docs/testing/automated_test_cases.md`
- `docs/testing/structure_validation_matrix.md`
- `docs/testing/diffraction_validation_matrix.md`
- `docs/testing/phase_transformation_validation_matrix.md`
- `benchmarks/validation/*.json`
- `tests/unit/`

Verification gates:

- each targeted subsystem has at least one stronger evidence block or literature-tracked acceptance statement
- each targeted subsystem keeps its manifests and fixtures synchronized with the documented claim
- transformation validation grows before transformation algorithm breadth grows

#### Workstream 1.4: Optional Interoperability Reality Check

Tasks:

- add live integration coverage for optional adapters where the dependency is available in CI or a dedicated environment
- keep adapter docs honest about what is normalized, what is validated, and what is still dependency-limited
- verify EBSD and structure-import bridge paths against representative external payloads rather than only object construction

Primary repository surfaces:

- `tests/integration/`
- `tests/unit/test_plotting_and_adapters.py`
- `tests/unit/test_orix_miller_adapter.py`
- adapter docs and workflow pages that currently describe interoperability

Verification gates:

- each claimed stable adapter path has at least one executable integration test or an explicit documented limitation
- adapter docs distinguish optional runtime dependency from canonical PyTex semantics

#### Workstream 1.5: Teaching-Surface Consolidation

Tasks:

- refresh the highest-value tutorials, especially the quickstart, core-model path, validation path, and notebook index
- make the notebook atlas mirror the documented learning sequence instead of acting as a parallel browse-only tree
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

### Phase 1 Suggested Execution Sequence

The workstreams above are coupled. The recommended implementation order is:

1. Workstream 1.1 first, because contributor and CI expectations must be explicit before stronger scientific claims can rely on them.
2. Workstream 1.2 second, because status docs and validation framing should reflect the already-green infrastructure baseline.
3. Workstream 1.3 third, because external-baseline growth is now the highest-value scientific hardening task.
4. Workstream 1.4 fourth, because interoperability claims should be tightened only after the validation posture is clearer.
5. Workstream 1.5 fifth, because teaching-surface consolidation should consume the corrected docs, validation, and interoperability posture rather than racing ahead of them.

### Phase 1 Deliverables By End Of Cycle

The repository should have all of the following by the end of this immediate phase:

- a documented and CI-exercised base lane plus full scientific lane
- a passing fixture-integrity and repo-integrity baseline
- a warning-free or materially warning-reduced Sphinx build for the roadmap-targeted docs
- upgraded validation pages with explicit evidence blocks and limitations
- conservative, test-backed adapter and CIF-backed workflow claims
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

1. executable quality lanes: keep the base lane green and formalize the full scientific lane
2. status and validation re-baselining: keep roadmap, README, and audit prose aligned with actual repo state
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
