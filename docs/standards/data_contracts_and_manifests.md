# Data Contracts And Manifests

PyTex will eventually exchange data across importers, adapters, benchmarks, examples, reports, and CLI workflows. Those boundaries need a formal contract before they become stable.

## Policy

- Stable workflow interchange must use versioned, machine-readable manifests.
- Schemas belong under `schemas/`.
- Manifests must preserve:
  - source provenance
  - canonicalized conventions
  - original source conventions when relevant
  - frame and symmetry metadata
  - software versions
  - benchmark or fixture identity when applicable
- Stable schemas must not be invented independently inside adapters or notebooks.

## Immediate Scope

The current repository does not yet expose stable manifest schemas, but the next adapter-heavy phases should add at least:

- import manifests for vendor and external-tool normalization
- benchmark fixture manifests
- workflow result manifests for reproducibility and reporting

## Minimum Manifest Fields

Every stable manifest family should include:

- schema identifier and version
- creation timestamp
- PyTex version
- source-system provenance
- canonical convention set identifier
- frame and symmetry declarations
- referenced files or fixture identifiers
- validation or benchmark context when applicable

## Design Rule

If a workflow output would be difficult to interpret without hidden context, it is not ready to be a stable interchange artifact until that context is captured by a manifest and schema.
