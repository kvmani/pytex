# AGENTS.md

This file defines how automation agents and human contributors must work in the PyTex repository.

PyTex is not a throwaway prototype. It is being built as a long-horizon scientific software library for crystallographic texture and diffraction, with equal support for research and teaching.

## Primary References

Read these first when working on the repository:

- `mission.md`
- `specifications.md`
- `docs/roadmap/critical_review_and_development_guide.md` (the governing development guide:
  current priorities, quality bars, and the explainable-results doctrine)
- `docs/architecture/orientation_relationship_analysis_foundation.md` (the OR-analysis flagship
  program)
- `docs/README.md`
- `docs/architecture/overview.md`
- `docs/architecture/canonical_data_model.md`
- `docs/testing/strategy.md`
- `docs/testing/mtex_parity_matrix.md`
- `docs/testing/vesta_parity_matrix.md`
- `docs/standards/engineering_governance.md`
- `docs/standards/notation_and_conventions.md`
- `docs/standards/documentation_architecture.md`
- `docs/standards/executable_examples.md`
- `docs/standards/latex_and_figures.md`
- `docs/standards/visualization_style_guide.md`
- `docs/standards/terminology_and_symbol_registry.md`
- `docs/standards/development_principles.md`
- `docs/standards/data_contracts_and_manifests.md`

If implementation choices conflict with these documents, stop and reconcile the conflict before continuing.

## Project Mandate

- Build a pure-Python-first, GPL-compatible texture and diffraction library.
- Treat the canonical crystallographic data model as a first-class product surface.
- Treat orientation-relationship analysis as the flagship capability; when priorities compete
  within a horizon, OR-analysis work per the OR foundation document takes precedence.
- Keep research-grade rigor and teaching-grade clarity in the same repo.
- Use external libraries through adapters where practical, but do not leak their raw domain semantics into PyTex public APIs.

## Non-Negotiable Rules

- Sphinx is the primary browsable and searchable documentation surface.
- LaTeX is the canonical source for major scientific notes.
- Stable public numerical surfaces must be documented with executable worked examples whose outputs
  are computed live from the code and checked against cited reference values. Documentation numbers
  must not be hand-transcribed. See `docs/standards/executable_examples.md`.
- Every documented method must explain its purpose, when and where to use it, its expected inputs,
  and its expected outputs, with at least one computed example for verifiable numerical behavior.
- Nomenclature and symbols are governed by one registry
  (`docs/standards/terminology_and_symbol_registry.md`); new symbols must be registered before use.
- Publication-quality SVG figures are mandatory where reference frames, geometry, or conventions matter.
- Canonical architecture, process-flow, validation, workflow, and teaching diagrams must follow the
  central visualization style guide.
- MTEX parity is the validation floor for relevant functionality, not the ceiling.
- Major scientific docs must cite normative sources explicitly.
- The preferred source hierarchy is: IUCr and International Tables, other formal standards, canonical textbooks, peer-reviewed papers, maintained tool documentation, then vendor notes.
- No stable public API may rely on naked arrays when frame, symmetry, or basis meaning would be ambiguous.
- No subsystem may define its own private frame or symmetry model.
- Stable features are incomplete until docs, figures, tests, and validation notes all exist.
- Stable terminology and symbol meaning must be fixed centrally and reused across docs, theory notes, notebooks, code explanations, and figures.
- Correctness, provenance, and interpretability take priority over premature optimization.
- Stable report and result objects must be explainable: they carry a `describe()` surface
  producing convention-explicit, citation-backed scientific prose, tested like any other output.
  See the explainable-results doctrine in
  `docs/roadmap/critical_review_and_development_guide.md`.
- New test or runtime warnings are defects; matplotlib figures opened by tests must be closed;
  measured coverage must not decrease.
- Where a literature convention exists (variant numbering, axis/angle representatives, section
  conventions), conform and pin it in tests, or document the deviation explicitly.
- Construction-time invariant checks are preferred over downstream error recovery.
- Any stable workflow that crosses a tool boundary must eventually have a machine-readable manifest and schema.

## Engineering Priorities

1. correctness
2. traceability
3. maintainability
4. interoperability
5. speed

Speed matters, but only after semantics are explicit and scientifically defensible.

## Expected Repository Shape

- `src/pytex/core/`: canonical data model, conventions, and low-level transformations
- `src/pytex/texture/`: PF, IPF, ODF, fibers, and texture-domain behavior
- `src/pytex/ebsd/`: crystal maps, grain workflows, and EBSD-specific semantics
- `src/pytex/diffraction/`: diffraction geometry, stereonets, and simulation-facing models
- `src/pytex/adapters/`: optional bridges to ORIX, KikuchiPy, PyEBSDIndex, pymatgen, diffsims
- `src/pytex/experimental/`: unstable research methods
- `docs/site/`: Sphinx-facing concepts, tutorials, worked examples, workflows, and curated API docs
- `docs/site/examples/`: generated executable-worked-example gallery (do not hand-edit)
- `docs/tex/`: canonical scientific documentation
- `docs/figures/`: canonical SVG figures
- `worked_examples/`: source of truth for executable worked examples (framework, registry, examples)

## Working Process

### Durable progress and resumability

- For every substantial multi-step task, keep a repository-local running progress note that records
  the objective, decisions, completed work, verification results, current worktree state, and exact
  next actions.
- Update the note throughout the task, especially before long-running verification, commits,
  pushes, or any likely interruption point, so another agent or a later session can resume without
  reconstructing context from chat history.
- Treat resumability as a primary working principle: leave the repository in an intelligible state,
  distinguish completed work from planned work, and record blockers or unverified assumptions
  explicitly.
- Use `docs/development/active_task_progress.md` for the current task unless a more specific tracked
  handoff document already exists. Archive or reset it only after the task is fully verified and its
  durable outcomes have been incorporated into canonical documentation.

### Before coding

- Read the closest authoritative docs for the affected subsystem.
- Confirm whether the task changes:
  - scientific semantics
  - public APIs
  - documentation standards
  - validation or test obligations
  - figure or notation requirements
- Prefer expanding the shared core model over encoding one-off conversions locally.

### During coding

- Keep public types explicit and strongly named.
- Prefer immutable metadata objects and contiguous NumPy-backed arrays for vectorized data.
- Write numerical routines as highly vectorized array operations (NumPy `einsum`/broadcasting, SciPy batch primitives); avoid Python per-element loops on the hot path. SciPy is a core dependency for linear algebra, optimization, spatial queries, and special functions. See `docs/standards/development_principles.md` principle 9a.
- Add tests with implementation, not afterward.
- For any stable public numerical surface, add or update an executable worked example in
  `worked_examples/`, regenerate the gallery with `python scripts/generate_worked_examples.py`, and
  keep `tests/unit/test_worked_examples.py` green. The example's expected value must have independent
  provenance (analytic identity or cited standard), never a copied prior program output.
- Update docs when behavior, conventions, or surface area changes — including `docs/README.md`
  (the index must stay complete) and any foundation document whose claims the change affects;
  stale foundational claims are treated as defects.
- Give new stable report objects a `describe()` method per the explainable-results doctrine, and
  keep JSON contracts and `describe()` in lockstep.
- Add or update cross-links when a page relies on terms, conventions, or workflows defined elsewhere in the docs.
- Add local module indexes or README files when a subsystem grows enough to need them.
- Treat repository artifact hygiene as part of the implementation task: generated outputs, local inspection assets, caches, build products, screenshots, notebooks checkpoints, benchmark scratch files, and similar non-canonical artifacts must be excluded in `.gitignore` before or alongside the change that creates them.
- Only commit generated files when they are intentional canonical repository assets referenced by docs, tests, manifests, validation workflows, or pinned regression baselines.
- Runtime plotting validation must prefer structural and semantic assertions over repo-tracked SVG byte baselines. Canonical SVG tracking is reserved for documentation figures in `docs/figures/`, not for routine runtime-regression fixtures.

### When touching algorithms

- Document the theory path in `docs/tex/`.
- Record assumptions, normalization rules, and failure modes.
- Record normative and informative citations in the corresponding Markdown and LaTeX docs.
- Update `docs/testing/mtex_parity_matrix.md` if the algorithm overlaps MTEX functionality.
- Add benchmark fixtures or explicit placeholders if the benchmark cannot land yet.
- Add a computed worked example that demonstrates the algorithm on a case with a known answer.

## Documentation And Executable-Examples Standard

PyTex documentation is held to a teaching-and-research standard in which every documented surface
answers four questions and, for numerical surfaces, proves the fourth:

1. What does it do (theory, mathematics, and algorithm)?
2. When and where is it used (scenarios, workflows, and the modules/methods involved)?
3. What are the expected inputs and outputs?
4. Does it actually produce the expected output? Demonstrate with an executable worked example whose
   value is computed live from the code and compared with a cited reference value.

The mechanics are fixed centrally:

- Executable worked examples live in `worked_examples/`, are rendered into `docs/site/examples/` by
  `scripts/generate_worked_examples.py`, and are validated by `tests/unit/test_worked_examples.py`.
- Theory and derivations remain in `docs/tex/`; interactive narratives remain in the notebook
  tutorials; canonical geometry remains in `docs/figures/` SVGs. Worked examples are the verifiable
  numerical bridge between them, not a replacement.
- Nomenclature is governed by `docs/standards/terminology_and_symbol_registry.md`; symbols must be
  registered before use across prose, math, figures, notebooks, worked examples, and docstrings.
- The docstring contract in `docs/standards/documentation_architecture.md` (purpose, when-to-use,
  parameters, returns, cross-references) applies to stable public surfaces; the `angle_*` helpers in
  `pytex.core.miller` are the reference exemplar.

## Anti-Goals

- Do not chase feature breadth before the core data model is coherent.
- Do not bury frame conventions in comments or hidden helper code.
- Do not create silent import-time coupling to heavy optional scientific stacks.
- Do not claim scientific equivalence to MTEX, ORIX, KikuchiPy, or other tools without tests or validation notes.
