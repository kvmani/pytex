# Data Contracts And Manifests

## Policy

- Stable workflow interchange must use versioned, machine-readable manifests.
- Major stable scientific outputs should also use versioned, machine-readable JSON contracts when they must cross workflow, storage, or tool boundaries.
- Schemas belong under `schemas/`.
- Manifests must preserve:
  - source-system provenance
  - canonical convention identifiers
  - frame and symmetry declarations
  - references to external files or fixtures
  - validation and benchmark context where relevant
- Stable schemas must not be invented independently inside adapters or notebooks.

## Broader Contract Rule

PyTex distinguishes between:

- `workflow manifests`
  Scientific context, provenance, validation, and workflow-bound evidence records.
- `object or result contracts`
  Canonical JSON serializations for major stable outputs that must be reconstructible into PyTex objects or stable result surfaces.

Both are part of the machine-readable interoperability posture of the project.

The expectation is:

- manifests describe workflows, evidence, and reproducibility context
- JSON object contracts describe major scientific outputs themselves

When a stable output is likely to be consumed by another tool, validation pipeline, or research agent, PyTex should prefer an explicit JSON contract over an informal ad hoc dictionary.

## Immediate Scope

The current repository exposes several stable manifest schemas today and a broader roadmap that must not be implemented ad hoc:

- import manifests for vendor and external-tool normalization
- benchmark fixture manifests
- workflow result manifests for reproducibility and reporting
- experiment manifests for multimodal acquisition and calibration state
- validation manifests for diffraction and structure-import baselines

## Minimum Manifest Fields

- schema identifier and version
- creation timestamp
- PyTex version
- source-system provenance
- canonical convention set identifier
- frame and symmetry declarations
- referenced files or fixture identifiers
- validation or benchmark context when applicable

## Design Rule

Stable manifests must describe scientific semantics explicitly enough that downstream workflows do not need hidden side channels to recover frame, symmetry, or provenance meaning.

The same rule applies to major stable JSON output contracts: they must preserve enough scientific semantics for reconstruction and downstream interpretation without hidden side channels.

## Manifest Family Roadmap

The repo now treats the following manifest families as first-class roadmap items:

- `import manifests`
  For vendor or external-tool normalization into PyTex primitives.
- `experiment manifests`
  For acquisition geometry, calibration state, detector metadata, and modality-specific context.
- `benchmark manifests`
  For pinned fixtures, comparison baselines, tolerances, and regeneration provenance.
- `workflow result manifests`
  For reproducible outputs, figure generation, and reporting.
- `validation manifests`
  For parity or literature-backed validation campaigns that span more than one script or test file.

The currently stable manifest families are:

- `EBSDImportManifest`
- `ExperimentManifest`
- `BenchmarkManifest`
- `ValidationManifest`
- `WorkflowResultManifest`

Stable manifests should be reused through the shared manifest surface rather than redefined inside subsystems.

## References

### Normative

- `reference_canon.md`
- Hall, S. R. and McMahon, B. (eds.), *International Tables for Crystallography, Volume G*

### Informative

- `../architecture/multimodal_characterization_foundation.md`
