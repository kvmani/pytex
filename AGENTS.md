# AGENTS.md

This file defines how automation agents and human contributors must work in the PyTex repository.

PyTex is not a throwaway prototype. It is being built as a long-horizon scientific software library for crystallographic texture and diffraction, with equal support for research and teaching.

## Primary References

Read these first when working on the repository:

- `mission.md`
- `specifications.md`
- `docs/README.md`
- `docs/architecture/overview.md`
- `docs/architecture/canonical_data_model.md`
- `docs/testing/strategy.md`
- `docs/testing/mtex_parity_matrix.md`
- `docs/standards/engineering_governance.md`
- `docs/standards/notation_and_conventions.md`
- `docs/standards/documentation_architecture.md`
- `docs/standards/latex_and_figures.md`
- `docs/standards/development_principles.md`
- `docs/standards/data_contracts_and_manifests.md`

If implementation choices conflict with these documents, stop and reconcile the conflict before continuing.

## Project Mandate

- Build a pure-Python-first, GPL-compatible texture and diffraction library.
- Treat the canonical crystallographic data model as a first-class product surface.
- Keep research-grade rigor and teaching-grade clarity in the same repo.
- Use external libraries through adapters where practical, but do not leak their raw domain semantics into PyTex public APIs.

## Non-Negotiable Rules

- Sphinx is the primary browsable and searchable documentation surface.
- LaTeX is the canonical source for major scientific notes.
- Publication-quality SVG figures are mandatory where reference frames, geometry, or conventions matter.
- MTEX parity is the validation floor for relevant functionality, not the ceiling.
- Major scientific docs must cite normative sources explicitly.
- The preferred source hierarchy is: IUCr and International Tables, other formal standards, canonical textbooks, peer-reviewed papers, maintained tool documentation, then vendor notes.
- No stable public API may rely on naked arrays when frame, symmetry, or basis meaning would be ambiguous.
- No subsystem may define its own private frame or symmetry model.
- Stable features are incomplete until docs, figures, tests, and validation notes all exist.
- Correctness, provenance, and interpretability take priority over premature optimization.
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
- `docs/site/`: Sphinx-facing concepts, tutorials, workflows, and curated API docs
- `docs/tex/`: canonical scientific documentation
- `docs/figures/`: canonical SVG figures

## Working Process

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
- Add tests with implementation, not afterward.
- Update docs when behavior, conventions, or surface area changes.
- Add local module indexes or README files when a subsystem grows enough to need them.

### When touching algorithms

- Document the theory path in `docs/tex/`.
- Record assumptions, normalization rules, and failure modes.
- Record normative and informative citations in the corresponding Markdown and LaTeX docs.
- Update `docs/testing/mtex_parity_matrix.md` if the algorithm overlaps MTEX functionality.
- Add benchmark fixtures or explicit placeholders if the benchmark cannot land yet.

## Anti-Goals

- Do not chase feature breadth before the core data model is coherent.
- Do not bury frame conventions in comments or hidden helper code.
- Do not create silent import-time coupling to heavy optional scientific stacks.
- Do not claim scientific equivalence to MTEX, ORIX, KikuchiPy, or other tools without tests or validation notes.
