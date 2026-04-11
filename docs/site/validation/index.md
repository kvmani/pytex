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

## Explicitly Checked Today

- Euler, quaternion, matrix, axis-angle, and Miller-surface semantics
- harmonic ODF reconstruction invariance, retained-basis diagnostics, and synthetic PF refit behavior
- regular-grid and graph-backed EBSD neighborhood workflows
- multiphase EBSD normalization, phase selection, and phase-aware texture extraction
- stable import, experiment, transformation, benchmark, validation, and workflow-result manifests
- detector geometry, reciprocal-space, and kinematic diffraction invariants
- CIF-backed phase creation, point-group preservation, and space-group preservation
- runtime plotting semantics for the highest-value publication-facing surfaces
- transformation-variant prediction and bounded experimental parent-candidate scoring

## Related Material

- {doc}`../benchmarks/index`
- {doc}`../theory/index`
- {doc}`../architecture/phase_transformation_foundation`
