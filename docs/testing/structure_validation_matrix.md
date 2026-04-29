# Structure Import Validation Matrix

This document is the authoritative validation ledger for PyTex structure-import and CIF-backed
phase-construction workflows.

## Status Keys

- `implemented`: automated coverage and validation notes exist for the current category
- `foundational`: the implementation exists and is scientifically structured, but the
  external-baseline surface is not yet complete
- `planned`: the category is accepted but not yet validated adequately
- `n/a`: not applicable to current PyTex scope, with explanation

## Matrix

| Area | Baseline | Status | Notes |
| --- | --- | --- | --- |
| CIF-backed phase creation | PyTex import tests and IUCr-backed CIF semantics | implemented | `Phase.from_cif()` and `Phase.from_cif_string()` are covered by automated tests. |
| Point-group and space-group preservation | IUCr-style crystallographic conventions and analyzer-backed invariants | implemented | `SpaceGroupSpec`, `Phase.space_group`, and point-group preservation are directly tested. |
| Unit-cell and atomic-site normalization | Internal invariant tests and structure-backed examples | implemented | Lattice, unit-cell, and site normalization are checked through CIF and `pymatgen`-backed cases. |
| Hash-pinned fixture-corpus audit | Catalog digests plus bundled fixture metadata and audit summaries | implemented | The bundled fixture corpus now records conventional and primitive expectations across all pinned starter phases. |
| Full-scientific-lane fixture loading | `pymatgen`-backed full-lane tests against the pinned fixture corpus | implemented | The controlling environment for CIF-backed claims is the full scientific lane, where the pinned fixture-audit coverage now executes without skips. |
| Cross-reference against literature or independent datasets | Textbook or literature-backed structure examples | foundational | Current evidence is still dominated by invariant-driven tests plus pinned fixture metadata rather than a broad literature ledger. |
| External-tool or external-dataset parity | Curated CIF corpora and pinned comparison baselines | foundational | More pinned benchmark inputs should be added before stronger equivalence claims are made. |

## Current Posture

The structure-import layer is semantically explicit, reproducible, and stronger than simple
parser-success testing. The main remaining gap is not basic correctness scaffolding. It is broader
external-baseline breadth beyond the current pinned starter corpus.

The strongest current claim is:

- the full scientific lane can reconstruct the bundled phase-fixture corpus into canonical PyTex
  objects with pinned space-group, point-group, lattice, and site-count expectations

The repo does not yet claim:

- broad parity against large external CIF corpora
- comprehensive literature-backed structure families beyond the starter set

## References

### Normative

- `strategy.md`
- `../standards/reference_canon.md`

### Informative

- `../tex/theory/crystal_structures_and_cif_import.tex`
- `../../benchmarks/structure_import/foundation_benchmark_manifest.json`
