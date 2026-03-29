# Visual Regression Fixtures

This directory contains pinned SVG outputs used by the publication-grade visualization regression
tests.

Current coverage:

- `xrd_ni_fcc_journal.svg`
  Powder XRD output for the pinned `ni_fcc` fixture using the `journal` theme.
- `saed_ni_fcc_dark.svg`
  SAED output for the pinned `ni_fcc` `[001]` case using the `dark` theme.
- `crystal_zr_hcp_journal.svg`
  Hexagonal crystal-scene output for the pinned `zr_hcp` fixture, including unit-cell, plane,
  and direction overlays.
- `ipf_ni_fcc_journal.svg`
  Inverse pole figure output for a pinned FCC orientation set.

The SVGs are generated with a fixed Matplotlib SVG hash salt and metadata date so the automated
tests can compare them deterministically while still keeping human-reviewable artifacts in-repo.
