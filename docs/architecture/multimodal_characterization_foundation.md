# Multimodal Characterization Foundation

This document defines the shared semantic backbone that PyTex must use across EBSD, XRD, neutron diffraction, and TEM-facing workflows.

## Why This Exists

PyTex already has meaningful texture, EBSD, and diffraction implementations. The next risk is not lack of capability. It is semantic fragmentation: each modality inventing its own frame chain, calibration model, and measurement metadata.

The role of this document is to prevent that.

## Shared Backbone

All modality-specific workflows must compose from one shared backbone:

- crystallographic structure and phase semantics
- explicit reference frames and frame transforms
- orientation and misorientation semantics
- reciprocal-space semantics
- acquisition geometry and calibration state
- provenance and measurement-quality metadata

Texture remains the semantic nucleus, but the surrounding modalities must bind to the same vocabulary rather than creating private interpretations.

## Canonical Frame Chain

The canonical frame chain for multimodal workflows is:

`crystal -> specimen -> map -> detector -> laboratory -> reciprocal`

Not every modality uses every link, but no stable workflow may bypass this vocabulary when the distinction matters.

- EBSD:
  Usually exercises `crystal -> specimen -> map`, and may later connect to detector geometry for pattern-level workflows.
- XRD and neutron diffraction:
  Require explicit `specimen -> laboratory -> reciprocal` reasoning and, where detector geometry matters, `detector` semantics as well.
- TEM diffraction:
  Requires explicit detector and reciprocal-space semantics while preserving the same crystal and specimen foundations.

## Required Future Primitive Families

The following concepts are now required architectural targets:

- `AcquisitionGeometry`
  Shared acquisition metadata that is broader than any single diffraction or EBSD container.
- `CalibrationRecord`
  Explicit calibration state, not hidden inside loose metadata.
- `MeasurementQuality`
  A stable place for confidence, masking, uncertainty, and quality metrics.
- `ScatteringSetup`
  Shared diffraction experiment context that can support XRD, neutron, and TEM variants without coupling them prematurely.

These types now exist as foundational core primitives, but broader modality-specific workflow adoption is still ahead.

## Modality Boundaries

PyTex should keep these boundaries explicit:

- `texture/orientation`
  Owns the mathematical representation of orientation, symmetry reduction, and texture-domain objects.
- `ebsd`
  Owns scan-topology, neighborhood, segmentation, and map-domain semantics.
- `diffraction`
  Owns detector, beam, scattering, and indexing semantics.
- `adapters`
  Own vendor or third-party normalization boundaries, not the core scientific meaning.

## Calibration And Quality Doctrine

Acquisition geometry is not complete without calibration and quality state.

Future multimodal APIs must be able to express:

- nominal geometry
- calibrated geometry
- confidence or uncertainty metadata
- masking and acceptance state
- source-system provenance for those states

These concepts must not be trapped inside vendor-specific adapters.

## Current Limits

- The shared acquisition, calibration, quality, and scattering primitives now exist, and stable experiment manifests now exist with them, but adoption across every modality surface is still incomplete.
- The current diffraction and EBSD layers are foundational, not yet a full multimodal platform.

## References

### Normative

- `../standards/reference_canon.md`
- `../standards/notation_and_conventions.md`
- `canonical_data_model.md`

### Informative

- `ebsd_foundation.md`
- `diffraction_foundation.md`
- Nolze et al. (2023), DOI: <https://doi.org/10.1107/S1600576723009275>
