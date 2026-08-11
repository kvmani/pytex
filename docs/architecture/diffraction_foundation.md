# Diffraction Foundation

This document records the initial Phase 4 implementation posture for diffraction-facing workflows.

## Implemented

- `DiffractionGeometry` as the canonical detector/specimen/laboratory geometry container
- electron wavelength from accelerating voltage
- detector-plane coordinates in millimeters relative to an explicit pattern center
- outgoing ray directions in the laboratory frame
- scattering-vector computation in reciprocal-length units
- detector-space $2\theta$ and azimuth evaluation
- Bragg-angle and ring-radius prediction from $d$ spacing or `CrystalPlane`
- powder XRD reflection enumeration and $2\theta$ spectrum generation with configurable radiation and broadening
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
- electron structure factors on an **absolute** scale (angstrom), from the Mott-Bethe
  conversion of the tabulated X-ray form factors with the relativistic correction applied,
  and the two-beam extinction distances that follow from them
- convergent-beam diffraction: disc geometry from the convergence semi-angle, the excitation
  error across each disc, the two-beam dynamical rocking curve, the Kossel-Moellenstedt
  versus Kossel regime test, and HOLZ ring radii
- foil-thickness determination from CBED fringe minima, returning the extinction distance in
  the same fit so the thickness does not depend on a tabulated constant

## Deliberate Current Limits

- no full detector-to-specimen transform calibration workflow yet
- dynamical intensity is **two-beam only**, and each CBED disc is computed independently, so
  the discs of one simulated pattern are not mutually consistent and their relative
  intensities carry no information; there is no many-beam or Bloch-wave solver
- no absorption model, so simulated thickness fringes do not decay as real ones do
- no HOLZ *line* simulation inside the bright-field disc, and no diffraction-group symmetry
  determination, so CBED cannot yet supply a point group or decide centrosymmetry
- the absolute structure-factor scale rests on a fitted parametrization: good to about
  1.5 percent for light elements (aluminium, validated against published extinction
  distances) and only indicative for heavy ones
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

- <a href="../site/theory/diffraction_geometry_and_bragg_rings.md">Diffraction Geometry And Bragg Rings</a>
- <a href="../site/theory/powder_xrd_and_saed.md">Powder XRD And SAED</a>
- <a href="../site/theory/reciprocal_space_and_kinematic_spots.md">Reciprocal Space And Kinematic Spots</a>
- <a href="../site/theory/convergent_beam_electron_diffraction.md">Convergent-Beam Electron Diffraction</a>
- <a href="../site/theory/dynamical_cbed_and_symmetry_determination.md">Dynamical CBED And Symmetry Determination</a>
