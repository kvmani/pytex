# Validation

PyTex validates stable work through multiple layers, not just unit tests. A feature is not considered scientifically settled merely because it executes without error.

## Validation Surfaces

- `docs/testing/strategy.md`
- `docs/testing/mtex_parity_matrix.md`
- `docs/testing/diffraction_validation_matrix.md`
- `docs/testing/structure_validation_matrix.md`
- `docs/testing/plotting_validation_matrix.md`
- [../../tex/validation/validation_program.tex](../../tex/validation/validation_program.tex)

## What Fixture-Backed Validation Means Here

For parity-covered areas, PyTex checks behavior against pinned MTEX-derived fixture bundles stored in `fixtures/mtex_parity/`. Those fixtures are versioned alongside the code and exercised through the parity suite in `tests/parity/`.

This is important for two reasons:

- the expected behavior is explicit and reviewable
- parity regressions become ordinary test failures rather than vague scientific suspicions

## MTEX Baseline

- pinned target: MTEX `6.1.1`
- pinned parity fixtures: `fixtures/mtex_parity/`
- parity test tree: `tests/parity/`

## Explicitly Checked

- Euler and quaternion convention handling for Bunge and Matthies or ABG-facing workflows
- quaternion, matrix, and axis-angle round trips against pinned fixture values
- vector and orientation projection into currently supported fundamental-region surfaces
- regular-grid EBSD KAM, GROD, grain segmentation, boundary extraction, and cleanup workflows
- stable EBSD import-manifest schema and normalization validation
- stable experiment, benchmark, validation, and workflow-result manifest schemas
- detector geometry, reciprocal-space, and kinematic diffraction invariants
- CIF-backed phase creation, point-group preservation, and space-group preservation
- semantic plotting inputs, contour plot builders, and ODF section rendering behavior

## Intentional Current Differences

- exact polyhedral orientation-region boundaries are not yet implemented for every crystal class
- parity fixtures are pinned and repo-local; future MTEX regeneration is prepared but not automated yet
- adapter-backed import normalization remains a separate follow-on phase
- diffraction now has an explicit validation ledger, but its external-baseline program is still foundational rather than mature
- structure import now has an explicit validation ledger, but its broader external-baseline program is still foundational rather than mature
- plotting now has an explicit validation ledger, but parity-oriented visual baselines are still ahead

## Reading The Validation Posture Correctly

“Parity-backed” does not mean “identical in every public API detail.” PyTex often uses clearer or narrower public surface names than MTEX while checking the underlying scientific behavior against MTEX-derived fixtures.

Likewise, “foundational” in the parity ledger does not mean “hand-wavy.” It means the implementation exists and is scientifically structured, but the automated parity surface is not yet complete enough to declare full coverage for that category.

## Current Posture

- core model invariants are covered by automated unit tests
- orientation and texture foundations have dedicated unit and parity tests
- EBSD parity coverage now includes KAM, GROD, segmentation, boundaries, and cleanup on regular grids
- diffraction foundations have dedicated unit tests and a separate validation ledger
- powder XRD and SAED foundations now sit inside that diffraction validation surface with dedicated unit tests
- structure-import foundations have dedicated unit tests and a separate validation ledger
- plotting foundations have dedicated unit tests and a separate validation ledger
- exact orientation-space polyhedral parity remains ahead of the current build

## Related Material

- `docs/testing/mtex_parity_matrix.md`
- `docs/testing/diffraction_validation_matrix.md`
- `docs/testing/structure_validation_matrix.md`
- `docs/testing/plotting_validation_matrix.md`
- [../../tex/validation/validation_program.tex](../../tex/validation/validation_program.tex)

## References

### Normative

- `../../testing/strategy.md`
- `../../testing/mtex_parity_matrix.md`
- `../../testing/diffraction_validation_matrix.md`
- `../../testing/structure_validation_matrix.md`
- `../../testing/plotting_validation_matrix.md`

### Informative

- MTEX documentation: [Kernel Average Misorientation](https://mtex-toolbox.github.io/EBSDKAM.html)
- MTEX documentation: [Grain Reference Orientation Deviation](https://mtex-toolbox.github.io/EBSDGROD.html)
- MTEX documentation: [Grain Orientation Parameters](https://mtex-toolbox.github.io/GrainOrientationParameters.html)
