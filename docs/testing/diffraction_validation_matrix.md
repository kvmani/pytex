# Diffraction Validation Matrix

This document is the authoritative validation ledger for PyTex diffraction-facing workflows.

## Status Keys

- `implemented`: automated coverage and validation notes exist for the current category
- `foundational`: the implementation exists and is scientifically structured, but the
  external-baseline surface is not yet complete
- `planned`: the category is accepted but not yet validated adequately
- `n/a`: not applicable to current PyTex scope, with explanation

## Matrix

| Area | Baseline | Status | Notes |
| --- | --- | --- | --- |
| Detector and beam geometry invariants | PyTex geometry tests and literature-backed conventions | implemented | Wavelength, detector projection, `2θ`, azimuth, and Bragg-ring geometry are directly tested. |
| Powder XRD reflection enumeration | Internal invariant tests and Bragg-law checks | implemented | `d` spacing, `2θ` filtering, and wavelength configuration are exercised through automated tests. |
| Powder XRD spectrum construction | Internal deterministic tests plus pinned `pymatgen` peak-position baselines | implemented | Broadening and plotting are stable, and pinned `ni_fcc` and `fe_bcc` Cu Kα reference cases now check peak positions and multiplicities. |
| Reciprocal-space primitives | IUCr-style crystallographic relations and internal invariant tests | implemented | `ReciprocalLatticeVector`, `CrystalPlane`, and `ZoneAxis` consistency is unit-tested. |
| SAED zone-axis spot generation | Internal geometric invariants and detector-coordinate tests | implemented | Zone-axis filtering, reciprocal construction, and detector mapping are directly tested. |
| Kinematic spot generation | Internal geometric invariants plus pinned `diffsims` shell-geometry baselines | implemented | Spot simulation, acceptance masks, and family grouping exist, and pinned `ni_fcc` and `fe_bcc` `[001]` shell-geometry cases now check external shell coverage. |
| Reflection family aggregation | Internal invariant tests | implemented | Multiplicity and grouping behavior are tested against symmetry-aware family keys. |
| Orientation candidate ranking and local refinement | Internal workflow tests | foundational | Deterministic local ranking exists, but continuous or statistically calibrated refinement is not yet in scope. |
| Intensity modeling | literature-backed physical models | planned | Current intensity is a proxy ranking model, not a full physical simulation. |
| XRD and SAED plotting | Runtime plotting tests and style-config tests | implemented | Plotters return Matplotlib figures and reuse the shared YAML style system. |
| External package or literature parity | pinned `pymatgen` XRD and `diffsims` SAED reference artifacts | foundational | Compact external-baseline cases are now pinned in-repo for two starter materials, but broader material and orientation coverage remains ahead. |

## Current Posture

The diffraction layer is scientifically meaningful, internally tested, and no longer limited to
purely internal checks. It now has pinned external-baseline cases for multiple starter materials
across one powder XRD workflow family and one SAED workflow family, but broader external coverage
remains ahead.

## References

### Normative

- `strategy.md`
- `../standards/reference_canon.md`

### Informative

- `../tex/algorithms/diffraction_geometry_and_bragg_rings.tex`
- `../tex/algorithms/powder_xrd_and_saed.tex`
- `../tex/algorithms/reciprocal_space_and_kinematic_spots.tex`
