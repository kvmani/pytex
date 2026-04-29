# Diffraction External Baselines

This directory contains compact pinned diffraction reference artifacts derived from the built-in
phase fixture corpus and independent open-source reference implementations.

Current scope:

- `ni_fcc_pymatgen_xrd_cuka.json`
  Peak-position and multiplicity baseline for powder XRD using `pymatgen` and the pinned `ni_fcc`
  fixture.
- `fe_bcc_pymatgen_xrd_cuka.json`
  Peak-position and multiplicity baseline for powder XRD using `pymatgen` and the pinned `fe_bcc`
  fixture.
- `ni_fcc_diffsims_saed_001_200kev.json`
  Indexed shell-geometry baseline for the `ni_fcc` `[001]` SAED case using `diffsims`.
- `fe_bcc_diffsims_saed_001_200kev.json`
  Indexed shell-geometry baseline for the `fe_bcc` `[001]` SAED case using `diffsims`.

These artifacts are intentionally small and deterministic so the full scientific lane can validate
PyTex against pinned external baselines without needing to regenerate the reference payloads during
ordinary contributor work.
