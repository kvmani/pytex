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
| Detector and beam geometry invariants | PyTex geometry tests and literature-backed conventions | implemented | Wavelength, detector projection, $2\theta$, azimuth, and Bragg-ring geometry are directly tested. |
| Powder XRD reflection enumeration | Internal invariant tests and Bragg-law checks | implemented | $d$ spacing, $2\theta$ filtering, and wavelength configuration are exercised through automated tests. |
| Powder XRD spectrum construction | Internal deterministic tests plus pinned `pymatgen` peak-position baselines | implemented | Broadening and plotting are stable; pinned `ni_fcc`, `fe_bcc`, and non-cubic `zr_hcp` Cu K-alpha cases check peak positions, multiplicities, and one emitted representative per symmetry family. |
| Reciprocal-space primitives | IUCr-style crystallographic relations and internal invariant tests | implemented | `ReciprocalLatticeVector`, `CrystalPlane`, and `ZoneAxis` consistency is unit-tested. |
| SAED zone-axis spot generation | Internal geometric invariants and detector-coordinate tests | implemented | Zone-axis filtering, reciprocal construction, and detector mapping are directly tested. |
| Kinematic spot generation | Internal geometric invariants plus pinned `diffsims` shell-geometry baselines | implemented | Spot simulation, acceptance masks, and family grouping exist, and pinned `ni_fcc` and `fe_bcc` `[001]` shell-geometry cases now check external shell coverage. |
| Reflection family aggregation | Internal invariant tests | implemented | Multiplicity and grouping behavior are tested against symmetry-aware family keys. |
| Full-scientific-lane diffraction baselines | `pymatgen` and `diffsims`-backed full-lane execution against pinned fixture payloads | implemented | The controlling environment for these claims is the full scientific lane, where the pinned diffraction baseline tests now execute without skips. |
| Orientation candidate ranking and local refinement | Internal workflow tests | foundational | Deterministic local ranking exists, but continuous or statistically calibrated refinement is not yet in scope. |
| Intensity modeling | Literature-backed physical models | planned | Current intensity is a proxy ranking model, not a full physical simulation. |
| XRD and SAED plotting | Runtime plotting tests and style-config tests | implemented | Plotters return Matplotlib figures and reuse the shared YAML style system. |
| Composite SAED reflection-table export (TX3) | Column-level identities: d = 1/\|g\| to machine precision on every row; detector radius equals the camera constant times the *in-plane* part of g exactly (and never exceeds camera constant x \|g\|); the table lists exactly the pattern's own spots; Friedel symmetry I(g) = I(-g) with mirrored detector positions for a centrosymmetric phase | implemented | `pytex.diffraction.export.composite_reflection_table(...)`. Every value is read from the engine's `SpotTable`, so the table cannot drift from the rendered figure. CSV column order is a declared public contract (`REFLECTION_TABLE_COLUMNS`). Worked example in the composite-diffraction gallery. |
| Child-zone-anchored composite geometry (TX4) | Consistency identity: anchoring on variant k's own view of a parent zone reproduces the parent-anchored pattern to 1e-13 mm on every sub-pattern, with identical hkl ordering, for k = 1..4. Burgers anchored on the alpha (0001) basal zone recovers a beta <110> parent direction at <1e-9 deg, and the anchor variant's own zone axis is exactly rational | implemented | `simulate_composite_saed_from_child_zone(...)` maps the child zone back through R_k^T and delegates to the parent-anchored engine, so there is one detector-geometry definition. `align_child_g` is expressed in the child's own indices. Worked example in the composite-diffraction gallery. |
| Kinematic spot-order stability | The same zone axis reached exactly and via a 1e-15 perturbation gives identical row order; repeated simulation is bitwise stable | implemented | Sort keys are quantized (1 pm of detector radius, 1e-12 of full-scale intensity) before `lexsort`, so symmetry-equivalent reflections fall through to the exact hkl tie-break instead of being ordered by floating-point noise. This was a real defect: the two anchoring routes previously produced correctly-positioned but permuted spot tables. |
| Composite SAED manifest and file export | Manifest validates against `schemas/composite_saed_manifest.schema.json`; its file inventory equals the set of files actually written; parent + per-variant reflection counts sum to the pattern's total; every variant's exact and nearest-rational child zone axis round-trips from the pattern | implemented | `export_composite_saed(...)` writes table, figure(s), coincidence table and manifest; rendering leaves no open matplotlib figures (asserted). |
| Lattice-centering audit | A phase declaring `Im-3m` reports centering `I` as *declared*; a phase with no space group reports `P` as *assumed* and triggers an explicit warning in `describe()` and the reflection table; disabling absences is stated explicitly. Body centring forbids h + k + l odd, and no such reflection survives for a declared bcc phase | implemented | `CompositeSAEDPattern.centering_audit()` / `phase_centering_is_declared(...)`. Closes a silent failure mode: an undeclared body-centred phase was simulated as primitive and listed forbidden reflections. |
| External package or literature parity | Pinned `pymatgen` XRD and `diffsims` SAED reference artifacts | foundational | Compact cubic and HCP external-baseline cases are pinned in-repo, but broader low-symmetry material and orientation coverage remains ahead. |

## Current Posture

The diffraction layer is scientifically meaningful, internally tested, and backed by real
external-baseline artifacts. The strongest current claim is geometric and shell-structure agreement
for the pinned starter cases, not full physical parity across intensity models or material space.

In practical terms:

- peak-position and multiplicity agreement are the hard powder-XRD claims
- shell geometry and family coverage are the hard SAED claims
- intensity differences remain informative rather than normative

## Evidence Hardening Queue

Before stronger diffraction claims are made, the next validation pass should add:

- additional open material fixtures beyond the new HCP starter case that exercise lower symmetry
  and overlapping multi-family reflection behavior
- separate ledgers for geometric agreement, shell/family agreement, and intensity-model limitations
- pinned external-baseline regeneration metadata for each new `pymatgen` or `diffsims` comparison
  artifact

## References

### Normative

- `strategy.md`
- `../standards/reference_canon.md`

### Informative

- `../tex/algorithms/diffraction_geometry_and_bragg_rings.tex`
- `../tex/algorithms/powder_xrd_and_saed.tex`
- `../tex/algorithms/reciprocal_space_and_kinematic_spots.tex`
