# MTEX Parity Matrix

This document is the authoritative ledger for PyTex parity against public MTEX tests and documented example categories.

Reference baseline:

- MTEX `6.1.1`
- public repository `tests/` tree
- public documentation examples where they expose behavior not captured by tests

## Status Keys

- `implemented`: PyTex has equivalent or stronger automated coverage
- `foundational`: PyTex has a correct foundational implementation, but not yet full parity with the MTEX behavior category
- `planned`: parity target is accepted but not implemented yet
- `n/a`: not applicable to current PyTex scope, with explanation

## Matrix

| Area | MTEX reference examples/tests | PyTex status | Notes |
| --- | --- | --- | --- |
| Euler/quaternion conversions | `check_eulerquat.m` | implemented | Baseline conversion and normalization tests are present, and the public rotation surface now also covers vectorized axis-angle, matrix, and Rodrigues / Rodrigues--Frank conversions. |
| Orientation construction and grids | `orientation.byEuler`, `orientation.byAxisAngle`, `orientation.byMiller`, `regularSO3Grid`, `equispacedSO3Grid` | foundational | Scalar orientation constructors, SO2 grids, regular Bunge SO3 grids, and deterministic quaternion SO3 grids now exist with PyTex frame and symmetry metadata. MTEX numerical parity is intentionally not claimed while external validation is paused. |
| Miller planes, directions, and family expansion | MTEX `Miller` object workflows and indexing examples | foundational | First-class `MillerPlane` / `MillerDirection` scalar and batch objects now cover vectorized reciprocal/direct vectors, d-spacings, hexagonal index conversion, angle relations, projection, symmetry-family expansion, JSON contracts, and orix interop. Broader MTEX parity for plotting idioms and every crystal-class helper remains ahead. |
| Fundamental regions | `check_fundamentalRegion.m` | foundational | Exact orbit-based minimum-angle reduction is implemented for the supported proper point groups, but broader class-by-class boundary catalogs are still ahead. |
| Symmetry operators and SO(3) basics | `SO3FunTests`, `check_WignerD.m`, symmetry-related checks | foundational | Common proper point-group generation, exact orbit reduction, symmetry actions, and a first Wigner-basis harmonic reconstruction surface are implemented; broader external parity is still ahead. |
| Spherical projections and stereonets | public spherical-projection examples and plotting workflows | foundational | Wulff-net plotting, stereographic direction/plane projection, great-circle traces, and rotational symmetry-axis symbols are implemented with deterministic regression coverage, but full MTEX visual-parity claims are still ahead. |
| EBSD container basics | `check_ebsd.m` | implemented | `CrystalMap` plus fixture-backed regular-grid segmentation coverage are implemented. |
| KAM-related behavior | `testKAM2.m` | implemented | Fixture-backed regular-grid KAM support covers order, threshold, and max-style aggregation. |
| GROD and grain-local orientation metrics | public EBSD workflow examples | implemented | Fixture-backed GROD relative to a representative grain orientation is implemented. |
| Grain-boundary and cleanup workflows | public EBSD segmentation examples | implemented | Fixture-backed boundary extraction and adjacency-based small-grain merging are implemented for regular grids. |
| IPF color coding | `checkIpfColorCoding.m` | foundational | `IPFColorKey` exists and is symmetry-aware, but full MTEX color-key parity is not yet claimed. |
| ODF and PF reconstruction | `check_FourierODF.m`, PF reconstruction examples | foundational | Discrete/kernel ODF evaluation, PF/IPF synthesis, explicit dictionary-based PF inversion, and a band-limited harmonic ODF inversion surface are implemented; broader MTEX-class correction workflows and parity breadth remain ahead. |
| Interfaces/import-export | `checkInterfaces.m` | foundational | Stable EBSD import manifests, manifest IO, object-backed vendor bridge adapters, CIF-backed phase creation, XRDML pole-figure import, LaboTex PPF/EPF pole-figure import, and lightweight orix / KikuchiPy bridge surfaces now exist, but dependency-pinned live-package parity is still ahead. |

## PyTex-Only Extensions

The following categories must exceed MTEX coverage:

- vendor reference-frame normalization
- provenance retention and manifest integrity
- adapter interoperability with ORIX, KikuchiPy, PyEBSDIndex, pymatgen, and diffsims
- LaTeX/SVG documentation asset integrity
- workflow-level reproducibility from CLI entry points

## References

### Normative

- `strategy.md`
- `../standards/reference_canon.md`

### Informative

- MTEX documentation: <https://mtex-toolbox.github.io/>
