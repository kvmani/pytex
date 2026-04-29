# Built-In Phase Fixtures

This directory contains the pinned structure-fixture corpus used for structure-import tests, benchmark manifests, validation manifests, and small reproducible demos.

## Starter Corpus

- `fe_bcc`
  Alpha iron in the body-centered cubic structure.
- `zr_hcp`
  Zirconium in the hexagonal close-packed structure.
- `ni_fcc`
  Nickel in the face-centered cubic structure.
- `nicl`
  Nickel chloride, stored as `NiCl2` while keeping the short fixture id requested for the starter corpus.
- `diamond`
  Diamond cubic carbon.

## Layout

Each fixture directory contains:

- `phase.cif`
  The pinned CIF used by tests and demos.
- `metadata.json`
  Machine-readable provenance, citation, redistribution, symmetry, lattice, and usage metadata.

The top-level [`catalog.json`](catalog.json) file is the authoritative manifest for this corpus.
Catalog entries are hash-pinned so the repo can verify fixture integrity as part of the reproducibility contract.

PyTex also exposes this corpus through the public loader helpers:

- `phase_fixture_catalog_path()`
- `list_phase_fixtures()`
- `get_phase_fixture(...)`

## Authoring Policy

- These fixtures are intentionally compact and suitable for in-repo use.
- The CIFs are pinned PyTex copies synthesized from canonical public crystallographic records and standard prototype symmetry so the test corpus remains stable.
- Metadata records point back to the upstream COD information cards and original cited publications.
- New bundled phase fixtures must be added to the catalog, include complete metadata, be hash-pinned, and be referenced by at least one benchmark or validation manifest.
