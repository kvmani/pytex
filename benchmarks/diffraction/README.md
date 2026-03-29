# Diffraction Benchmarks

This directory holds benchmark manifests for implemented diffraction workflows.

Current scope:

- detector and beam geometry invariants
- kinematic spot simulation
- reflection-family grouping
- orientation candidate ranking and local refinement
- compact phase-backed diffraction demos and validation inputs derived from `fixtures/phases/`
- pinned external-baseline reference artifacts for powder XRD and SAED derived from the `ni_fcc` fixture

The external-baseline program currently uses compact pinned artifacts generated from independent
open-source reference implementations and versioned in-repo so the default test suite can validate
geometry agreement without fetching data at runtime.
