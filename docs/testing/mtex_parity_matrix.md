# MTEX Parity Matrix

This document is the authoritative ledger for PyTex parity against public MTEX tests and documented example categories.

Reference baseline:

- MTEX `6.0.0`
- public repository `tests/` tree
- public documentation examples where they expose behavior not captured by tests

The campaign-based parity runner targets MTEX `6.0.0`, released in November 2024, so validation
does not depend on the newest MTEX release line. Regeneration scripts live under
`scripts/mtex_generators/`, and shared campaign inputs live under `fixtures/mtex_parity/campaigns/`.

## Status Keys

- `implemented`: PyTex has equivalent or stronger automated coverage
- `foundational`: PyTex has a correct foundational implementation, but not yet full parity with the MTEX behavior category
- `planned`: parity target is accepted but not implemented yet
- `n/a`: not applicable to current PyTex scope, with explanation

## Matrix

| Area | MTEX reference examples/tests | PyTex status | Notes |
| --- | --- | --- | --- |
| Euler/quaternion conversions | `check_eulerquat.m`; `fixtures/mtex_parity/campaigns/orientation_core_cases.json` | implemented | Baseline conversion and normalization tests are present, and the public rotation surface now also covers vectorized axis-angle, matrix, and Rodrigues / Rodrigues--Frank conversions. The campaign runner adds shared JSON inputs and MTEX/PyTex result JSON for Euler, matrix, quaternion, axis-angle, Miller construction, composition, inverse, and misorientation cases. |
| Orientation construction and grids | `orientation.byEuler`, `orientation.byAxisAngle`, `orientation.byMiller`, `regularSO3Grid`, `equispacedSO3Grid` | foundational | Scalar and batch orientation constructors now cover Euler, axis-angle, matrix, quaternion, Rodrigues / Rodrigues-Frank, Miller plane-direction correspondences, regular Bunge SO3 grids, and deterministic quaternion SO3 grids with shared frame and symmetry checks. Campaign coverage exists for scalar construction and operations; external grid parity breadth remains ahead. |
| Miller planes, directions, and family expansion | MTEX `Miller` object workflows and indexing examples; `fixtures/mtex_parity/campaigns/miller_geometry_cases.json` | foundational | First-class `MillerPlane` / `MillerDirection` scalar and batch objects now cover vectorized reciprocal/direct vectors, d-spacings, hexagonal index conversion, angle relations, projection, symmetry-family expansion, JSON contracts, orix interop, and Miller-backed orientation construction from three- or four-index inputs. The campaign runner adds compact cubic and hexagonal geometry/family checks. Broader MTEX parity for plotting idioms and every crystal-class helper remains ahead. |
| Fundamental regions | `check_fundamentalRegion.m` | foundational | Exact orbit-based minimum-angle reduction is implemented and now regression-tested across the supported proper point groups, but broader class-by-class closed-form boundary catalogs are still ahead. |
| Symmetry operators and SO(3) basics | `SO3FunTests`, `check_WignerD.m`, symmetry-related checks | foundational | Common proper point-group generation, exact orbit reduction, symmetry actions, and a first Wigner-basis harmonic reconstruction surface are implemented; broader external parity is still ahead. |
| Spherical projections and stereonets | public spherical-projection examples and plotting workflows | foundational | Wulff-net plotting, stereographic direction/plane projection, great-circle traces, and rotational symmetry-axis symbols are implemented with deterministic regression coverage, but full MTEX visual-parity claims are still ahead. The core model now also provides `vector3d`-style semantic sphere primitives: `SphericalVectorSet` (polar constructors, antipodal-aware angles, spherical/orientation-tensor means, folding) and `S2Grid` (equispaced and regular sphere grids with normalized band-area quadrature weights), with the canonical polar-angle convention centralized in `core/sphere.py` and reused by the stereonet surface. |
| EBSD container basics | `check_ebsd.m` | implemented | `CrystalMap` plus fixture-backed regular-grid segmentation coverage are implemented. |
| KAM-related behavior | `testKAM2.m` | implemented | Fixture-backed regular-grid KAM support covers order, threshold, and max-style aggregation. |
| GROD and grain-local orientation metrics | public EBSD workflow examples | implemented | Fixture-backed GROD relative to a representative grain orientation is implemented. |
| Grain-boundary and cleanup workflows | public EBSD segmentation examples | implemented | Fixture-backed boundary extraction and adjacency-based small-grain merging are implemented for regular grids. |
| IPF color coding | `checkIpfColorCoding.m`; `fixtures/mtex_parity/campaigns/ipf_color_cases.json` | foundational | `IPFColorKey` exists and is symmetry-aware. The campaign runner now emits comparable RGB and reduced-direction JSON for cubic and hexagonal cases, but full MTEX color-key visual parity is not yet claimed. |
| ODF and PF reconstruction | `check_FourierODF.m`, PF reconstruction examples; `fixtures/mtex_parity/campaigns/odf_discrete_cases.json` | foundational | Discrete/kernel ODF evaluation, PF/IPF synthesis, explicit dictionary-based PF inversion, and a band-limited harmonic ODF inversion surface are implemented. The campaign runner currently covers sampled discrete ODF density from weighted Euler orientations. XRDML PF and reconstruction campaigns are defined as pending until cubic and hexagonal fixtures are supplied. |
| Interfaces/import-export | `checkInterfaces.m`; `fixtures/mtex_parity/campaigns/xrdml_pole_figure_cases.json` | foundational | Stable EBSD import manifests, manifest IO, object-backed vendor bridge adapters, CIF-backed phase creation, XRDML pole-figure import, LaboTex PPF/EPF pole-figure import, and lightweight orix / KikuchiPy bridge surfaces now exist. Shared XRDML parity campaign inputs are present but pending real cubic and hexagonal XRDML fixtures. |

## PyTex-Only Extensions

The following categories must exceed MTEX coverage:

- vendor reference-frame normalization
- provenance retention and manifest integrity
- adapter interoperability with ORIX, KikuchiPy, PyEBSDIndex, pymatgen, and diffsims
- LaTeX/SVG documentation asset integrity
- workflow-level reproducibility from CLI entry points

## Current Hardening Priorities

The highest-value remaining MTEX-facing hardening work is concentrated in the foundational rows
that underpin multiple downstream subsystems:

- orientation construction and deterministic SO(3) grid coverage
- Miller plane and direction workflows beyond the current starter families
- broader class-by-class fundamental-region coverage
- richer symmetry and harmonic parity breadth
- broader ODF/PF reconstruction breadth
- stronger interface and import or export parity at the adapter boundary

## Campaign Regeneration Workflow

The detailed step-by-step workflow is documented in:

- `docs/site/validation/mtex_regeneration.md`
- `scripts/mtex_generators/README.md`

Minimal example on a MATLAB system with MTEX started:

```matlab
addpath("scripts/mtex_generators")
run_mtex_parity_campaign("fixtures/mtex_parity/campaigns/orientation_core_cases.json", ...
                         "fixtures/mtex_parity/results/mtex")
```

Repeat for each campaign file, bring `fixtures/mtex_parity/results/mtex/` back to the PyTex
machine, then run:

```powershell
python scripts/generate_pytex_parity_campaign.py fixtures/mtex_parity/campaigns fixtures/mtex_parity/results/pytex
python scripts/compare_parity_results.py fixtures/mtex_parity/results/mtex fixtures/mtex_parity/results/pytex
```

## References

### Normative

- `strategy.md`
- `../standards/reference_canon.md`

### Informative

- MTEX documentation: <https://mtex-toolbox.github.io/>
