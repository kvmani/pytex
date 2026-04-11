# Validation

PyTex validates stable work through multiple layers, not just unit tests. A feature is not
considered scientifically settled merely because it executes without error.

```{toctree}
:maxdepth: 1

testing_strategy
automated_test_cases
mtex_parity_matrix
diffraction_validation_matrix
structure_validation_matrix
plotting_validation_matrix
phase_transformation_validation_matrix
```

## What Validation Means Here

- deterministic unit tests for invariants and public API behavior
- human-auditable documented test cases for important formulas and conversions
- parity ledgers where MTEX is the correct scientific floor
- benchmark manifests that pin workflow identity and tolerances
- validation manifests that record what is and is not scientifically claimed
- theory notes that document assumptions, approximations, and current limits

## How To Use This Section

Read the validation surface in this order:

1. start here to identify which categories PyTex currently validates
2. open the relevant matrix or documented test-case page
3. check whether the result is `implemented`, `foundational`, or still limited
4. read the linked workflow or theory note before treating a numerical match as a broad scientific claim

The key question is not only “does a test pass?” The key question is “what exact scientific claim
does this passing test justify?”

## Explicitly Checked Today

- Euler, quaternion, matrix, axis-angle, and Miller-surface semantics
- harmonic ODF reconstruction invariance, retained-basis diagnostics, and synthetic PF refit behavior
- regular-grid and graph-backed EBSD neighborhood workflows
- multiphase EBSD normalization, phase selection, and phase-aware texture extraction
- optional ORIX adapter transfer for Miller, rotation, orientation, symmetry, and phase semantics
- stable import, experiment, transformation, benchmark, validation, and workflow-result manifests
- detector geometry, reciprocal-space, and kinematic diffraction invariants
- CIF-backed phase creation, point-group preservation, and space-group preservation
- runtime plotting semantics for the highest-value publication-facing surfaces
- transformation-variant prediction and bounded experimental parent-candidate scoring

## What Is Still Limited

- Some validated surfaces are foundational rather than broad parity claims.
- Diffraction coverage is still stronger on geometry and kinematics than on full physical intensity
  modeling.
- Adapter support is not automatically equivalent to broad live-package interoperability coverage.
- A passing validation page does not imply every neighboring workflow variant is equally mature.

For optional interoperability specifically, the controlling question is: which semantic boundary has
been tested, under which dependency assumptions, and with which explicit limits?

Where those limits matter, the detailed validation pages should be treated as the controlling
scientific statement.

## Related Material

- {doc}`../benchmarks/index`
- {doc}`../theory/index`
- {doc}`../architecture/phase_transformation_foundation`
- {doc}`automated_test_cases`
