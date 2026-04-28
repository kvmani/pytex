# Diffraction Foundation

This document records the initial Phase 4 implementation posture for diffraction-facing workflows.

## Implemented

- `DiffractionGeometry` as the canonical detector/specimen/laboratory geometry container
- electron wavelength from accelerating voltage
- detector-plane coordinates in millimeters relative to an explicit pattern center
- outgoing ray directions in the laboratory frame
- scattering-vector computation in reciprocal-length units
- detector-space `2θ` and azimuth evaluation
- Bragg-angle and ring-radius prediction from `d` spacing or `CrystalPlane`
- powder XRD reflection enumeration and `2θ` spectrum generation with configurable radiation and broadening
- explicit `RadiationSpec`, `PowderReflection`, and `PowderPattern` objects for XRD workflows
- explicit `ReciprocalLatticeVector` and `ZoneAxis` core-model objects
- specimen-to-laboratory rotation as an explicit diffraction-geometry contract
- minimal Ewald-style kinematic spot simulation with excitation-error filtering and detector projection
- explicit SAED spot-pattern generation from a `ZoneAxis` with detector coordinates in a named detector frame
- validation of detector projection edge cases, integer Miller inputs, and off-detector spot semantics
- symmetry-aware reflection-family grouping with explicit multiplicity records
- explicit detector acceptance masks for workflow-level detector gating
- minimal proxy intensity weighting for spot ranking and family representation
- detector-space clustering and simulated/observed indexing-candidate association
- local orientation-candidate ranking and deterministic local refinement
- family-level indexing reports built from matched reflection families

## Deliberate Current Limits

- no full detector-to-specimen transform calibration workflow yet
- no physically complete structure-factor or dynamical intensity simulation yet
- no continuous or probabilistic orientation-refinement workflow yet
- no adapter-backed bridges to diffsims or related diffraction stacks yet

## Why This Is The Right First Step

Phase 4 should begin by making reciprocal-space and detector-space semantics explicit before attempting simulation breadth. A usable diffraction foundation needs detector coordinates, wave-vector transfer, Bragg geometry, and clear frame ownership before more elaborate indexing or pattern-generation work can stay interpretable.

## References

### Normative

- [Canonical Data Model](canonical_data_model.md)
- [Multimodal Characterization Foundation](multimodal_characterization_foundation.md)
- [Reference Canon](../standards/reference_canon.md)

### Informative

- <a href="../tex/algorithms/diffraction_geometry_and_bragg_rings.tex">Diffraction Geometry And Bragg Rings</a>
- <a href="../tex/algorithms/powder_xrd_and_saed.tex">Powder XRD And SAED</a>
- <a href="../tex/algorithms/reciprocal_space_and_kinematic_spots.tex">Reciprocal Space And Kinematic Spots</a>
