# PyTex

PyTex is a GPL-compatible, pure-Python-first library for crystallographic texture and diffraction with a deliberate focus on materials-science research and teaching.

The repository is being built around three non-negotiable foundations:

- a canonical crystallographic data model for frames, symmetry, orientations, maps, and diffraction geometry
- hybrid scientific documentation: Sphinx for the primary browsable/searchable docs surface, LaTeX for authoritative scientific notes, and SVG for canonical figures
- MTEX-plus validation, where MTEX parity is the baseline and PyTex-specific interoperability and provenance checks extend beyond it
- explicit normative-source hierarchy, so conventions are fixed from authoritative references and not re-litigated locally

## Start Here

- `mission.md`
- `specifications.md`
- `AGENTS.md`
- `docs/README.md`
- `docs/architecture/overview.md`
- `docs/architecture/canonical_data_model.md`
- `docs/testing/strategy.md`
- `docs/testing/mtex_parity_matrix.md`
- `docs/roadmap/implementation_roadmap.md`
- `docs/standards/notation_and_conventions.md`
- `docs/standards/latex_and_figures.md`
- `docs/standards/documentation_architecture.md`
- `docs/standards/scientific_citation_policy.md`
- `docs/standards/benchmark_and_tolerance_governance.md`
- `docs/standards/hexagonal_and_trigonal_conventions.md`
- `docs/standards/development_principles.md`
- `docs/standards/data_contracts_and_manifests.md`

## Current Status

The repository currently implements Phases 0 through 2 of the foundation roadmap:

- modern Python packaging and CI skeleton
- canonical core data structures under `src/pytex/`
- documentation governance and hybrid scientific doc scaffold
- validation strategy and MTEX parity matrix scaffold
- baseline tests for the foundational data model
- symmetry-aware orientation and disorientation foundations
- PF/IPF containers, class-specific IPF sector reduction, and discrete kernel-based ODF foundations
- minimal regular-grid EBSD neighborhood, KAM, grain segmentation, GROD, grain-boundary, and cleanup workflow foundations

Exact orientation-space fundamental-region polyhedra, harmonic ODF inversion, richer EBSD workflows, diffraction engines, and external-tool adapters are intentionally staged after this foundation so they do not invent conflicting conventions.

## Quick Start

Install the package in editable mode with development tools:

```bash
python -m pip install -e '.[dev]'
python scripts/check_repo_integrity.py
pytest
```

Inspect the documentation inventory from the CLI:

```bash
pytex info
pytex docs inventory
pytex core demo
```

## Repository Layout

```text
pytex/
├─ src/pytex/
│  ├─ core/
│  ├─ texture/
│  ├─ ebsd/
│  ├─ diffraction/
│  ├─ adapters/
│  ├─ plotting/
│  └─ experimental/
├─ tests/
├─ docs/
│  ├─ architecture/
│  ├─ testing/
│  ├─ roadmap/
│  ├─ standards/
│  ├─ development/
│  ├─ site/
│  ├─ tex/
│  └─ figures/
├─ fixtures/
├─ benchmarks/
├─ schemas/
├─ examples/
└─ scripts/
```

## Design Direction

- Own the domain model instead of leaking raw arrays through public APIs where frame or symmetry meaning would be ambiguous.
- Reuse proven projects such as ORIX, KikuchiPy, PyEBSDIndex, pymatgen, and diffsims through adapters instead of coupling the whole library to any single external representation.
- Treat documentation, figures, and validation artifacts as product deliverables rather than release polish.
- Keep research-grade depth and teaching-grade clarity in the same repository.

## License

PyTex is released under the GPL-3.0-or-later license. See `LICENSE` for the repository license notice. The licensing posture is intentional so GPL-compatible scientific dependencies can be integrated cleanly where that makes technical sense.
