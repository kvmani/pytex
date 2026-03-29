# Diffraction External Baselines

This directory contains compact pinned diffraction reference artifacts derived from the built-in
phase fixture corpus and independent open-source reference implementations.

Current scope:

- `ni_fcc_pymatgen_xrd_cuka.json`
  Peak-position and multiplicity baseline for powder XRD using `pymatgen` and the pinned `ni_fcc`
  fixture.
- `ni_fcc_diffsims_saed_001_200kev.json`
  Indexed shell-geometry baseline for the `ni_fcc` `[001]` SAED case using `diffsims`.

These artifacts are intentionally small and deterministic so the default test suite can validate
PyTex against pinned external baselines without needing to install the generator packages during
ordinary test execution.
