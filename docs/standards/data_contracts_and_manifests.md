# Data Contracts And Manifests

## Policy

- Stable workflow interchange must use versioned, machine-readable manifests.
- Schemas belong under `schemas/`.
- Manifests must preserve:
  - source-system provenance
  - canonical convention identifiers
  - frame and symmetry declarations
  - references to external files or fixtures
  - validation and benchmark context where relevant
- Stable schemas must not be invented independently inside adapters or notebooks.

## Immediate Scope

The current repository exposes one stable manifest schema today and a broader roadmap that must not be implemented ad hoc:

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

Only the EBSD import manifest is currently stable. The others are required architectural targets and should not be implemented ad hoc inside subsystems.

## References

### Normative

- `reference_canon.md`
- Hall, S. R. and McMahon, B. (eds.), *International Tables for Crystallography, Volume G*

### Informative

- `../architecture/multimodal_characterization_foundation.md`
