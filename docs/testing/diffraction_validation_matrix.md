# Diffraction Validation Matrix

This document is the authoritative validation ledger for PyTex diffraction-facing workflows.

## Status Keys

- `implemented`: automated coverage and validation notes exist for the current category
- `foundational`: the implementation exists and is scientifically structured, but the external-baseline surface is not yet complete
- `planned`: the category is accepted but not yet validated adequately
- `n/a`: not applicable to current PyTex scope, with explanation

## Matrix

| Area | Baseline | Status | Notes |
| --- | --- | --- | --- |
| Detector and beam geometry invariants | PyTex geometry tests and literature-backed conventions | implemented | Wavelength, detector projection, `2θ`, azimuth, and Bragg-ring geometry are directly tested. |
| Reciprocal-space primitives | IUCr-style crystallographic relations and internal invariant tests | implemented | `ReciprocalLatticeVector`, `CrystalPlane`, and `ZoneAxis` consistency is unit-tested. |
| Kinematic spot generation | Internal geometric invariants and current workflow examples | foundational | Spot simulation, acceptance masks, and family grouping exist; broader external comparison remains ahead. |
| Reflection family aggregation | Internal invariant tests | implemented | Multiplicity and grouping behavior are tested against symmetry-aware family keys. |
| Orientation candidate ranking and local refinement | Internal workflow tests | foundational | Deterministic local ranking exists, but continuous or statistically calibrated refinement is not yet in scope. |
| Intensity modeling | literature-backed physical models | planned | Current intensity is a proxy ranking model, not a full physical simulation. |
| External package or literature parity | diffsims and literature-backed spot geometry comparisons | planned | Needed before stronger public equivalence claims are made. |

## Current Posture

The diffraction layer is scientifically meaningful and internally tested, but it is still a foundational implementation rather than a finished external-baseline program.

## References

### Normative

- `strategy.md`
- `../standards/reference_canon.md`

### Informative

- `../tex/algorithms/diffraction_geometry_and_bragg_rings.tex`
- `../tex/algorithms/reciprocal_space_and_kinematic_spots.tex`
