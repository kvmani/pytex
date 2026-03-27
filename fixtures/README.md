# Fixtures

This directory contains canonical datasets and placeholders used by tests, examples, and validation workflows.

The initial scaffold focuses on semantic fixtures and compact, pinned benchmark assets rather than large experimental datasets.

## Current Fixture Families

- `fixtures/mtex_parity/`
  MTEX-aligned parity fixtures for rotations, symmetry reduction, and EBSD behavior.
- `fixtures/phases/`
  Built-in crystallographic phase fixtures backed by pinned CIF files and explicit provenance metadata.

## Fixture Policy

- Bundled scientific fixtures must be small, redistributable, and provenance-documented.
- Phase fixtures must include:
  - a pinned CIF file
  - machine-readable metadata with source, citation, redistribution note, expected symmetry, and intended uses
  - a stable fixture id referenced from tests, docs, and benchmark or validation manifests
- The preferred source family for bundled phase fixtures is the Crystallography Open Database or similarly canonical open crystallographic sources.
- Large benchmark corpora should remain manifest-addressable and be fetched separately in later phases rather than committed directly to the repo.
