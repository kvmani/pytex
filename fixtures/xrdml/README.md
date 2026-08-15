# XRDML Fixtures

This directory stores pinned XRDML pole-figure fixtures used to validate PyTex's texture-import boundary.

Current contents:

- `polefig_Ge113_xrayutilities.xrdml.bz2`
  Source: `xrayutilities` example data
  Upstream URL: `https://github.com/dkriegner/xrayutilities`
  Original example path: `examples/data/polefig_Ge113.xrdml.bz2`
  License context: GPL-2.0-or-later as distributed by the upstream project
  Role in PyTex: real-file regression coverage for PANalytical-style pole-figure import
- `synthetic_random_standard.xrdml`
  Source: constructed PyTex validation fixture; explicitly not an experimental measurement
  Analytic definition: after subtracting 10 counts, the 0°, 30°, and 60° rings are
  100, 80, and 50 counts, so their exact defocusing factors are 1, 0.8, and 0.5
  Role in PyTex: end-to-end XRDML import and random-standard defocusing calibration

The fixture is used for import-shape and metadata validation only. PyTex does not infer the crystallographic pole from this file; callers must supply the corresponding `CrystalPlane` explicitly because the reflection indices are not encoded in a stable machine-readable field in this example file.

The current validation surface also checks that vendor axis names are normalized onto the PyTex specimen-direction convention explicitly rather than implicitly, and that any intensity normalization or random-standard correction applied later is chosen by the caller instead of being silently inferred during import.
