# Repository review, September 2026

This review starts at `60ffbff` (version 0.8.1). It covers package structure,
scientific contracts, executable examples, application workflows, validation and release
metadata. It is a repository-wide engineering and scientific-method review, not an external
validation of every algorithm. Verification is recorded in the
[active progress ledger](active_task_progress.md).

## Findings and disposition

| Area | Evidence from the current tree | Action or remaining boundary |
| --- | --- | --- |
| Core orientation relationships | Both fitting entry points used equal weights only; all pair/symmetry matrices were retained simultaneously. Iteration settings were unchecked, characterization accepted incompatible nominal/catalog phases, and a nonconverged fit could be conclusive. | Add declared pair weights with retained exclusions, bounded alignment, accurate final residuals, phase validation and convergence-aware naming. Analytic circular-mean tests supplement existing multi-variant tests. |
| Core architecture | `core.parent_reconstruction` imported `experimental.phase_transformation.score_parent_orientations`, despite the old guide identifying it as a priority. | Shared implementation moved to private core; the experimental facade preserves callers. Import-isolation and facade identity tests enforce the boundary. Map reconstruction retains experimental status. |
| Texture | Kernels, harmonic reconstruction, corrections, Kearns analysis and component/fibre workflows exist; old missing-feature lists are not reliable implementation inventories. | Retain the distinction between synthetic tests and measured/external parity in the MTEX ledger. No new ODF uncertainty or inversion-completeness claim. |
| EBSD | Canonical maps, vendor readers, segmentation, local metrics and measured-pair OR services exist. | Improve measured-pair evidence handling. Larger map backing and measured parent-reconstruction validation remain separate work. |
| Diffraction | XRD fitting, phase identification, SAED, CBED and recent HREM code exist. HREM tests assumed optional ASE/abTEM were installed; the adapter omitted coma/trefoil coefficients. | Test capability detection independently of the machine; forward all represented optical coefficients and compare phase transfer against installed abTEM. Preserve the distinction between phase-object and multislice models. |
| TEM and material properties | Stage/tilt geometry, calibration, slip, Taylor and elastic tensor surfaces have domain tests and examples. Abstract stage methods are intentional, not unfinished concrete operations. | Retain their existing contracts and numerical checks; no speculative feature additions. |
| Plotting | Shared themes, geometry builders, SVG governance and introspected class diagrams exist. The class atlas was stale after the preceding HREM addition. | Regenerate named canonical assets and reconcile the atlas counts; no runtime screenshot baselines added. |
| GUI | OR paste parser silently discarded rows with a non-six-column shape. Euler labels were handwritten despite the central symbol registry. | Reject a malformed paste atomically and preserve the pasted text; accept six angles plus an optional weight, display included/excluded evidence, and use registered symbols. Verify actual browser behavior. |
| Scientific documentation | OR theory called the Markley eigen-mean a Karcher mean, and some prose described independently reducing all pairs, contrary to the implementation and its Bain regression. | State the chordal objective, single-positive-weight seed, assumptions and limits. Add computed weighted examples and practical measurement guidance. |
| Notebook documentation | Notebook 33 has 38 cells without IDs; Sphinx's warning counter missed the resulting `nbformat` warning. | Add stable source cell IDs without outputs or execution counts; verify the next build's complete log. |
| CLI and release | The examples command lost its integer return during the preceding HREM addition, failing strict typing. Release policy still tells contributors to edit a literal version in `pyproject.toml`, although `_version.py` is authoritative. | Restore the return and align release instructions before the version bump. |

## What the review does not establish

Weighting records a user's evidence choices; it is not robust outlier inference, a confidence
interval, or proof that rows are independent. A converged local fit is not a global optimum.
MTEX measured-pair fitting and measured-data reconstruction validation remain outstanding.
Existing synthetic tests and literature identities must not be relabeled as those comparisons.

The first baseline browser run passed 60 tests directly and two on retry. Those retries are
recorded as timing sensitivities, not a clean first-attempt result. Final browser verification
must be run without concurrent CPU-intensive scientific checks before deciding whether fixes
are needed.

## References

- [Governing development guide](../roadmap/critical_review_and_development_guide.md)
- [OR foundation](../architecture/orientation_relationship_analysis_foundation.md)
- [Testing strategy](../testing/strategy.md)
- [MTEX parity ledger](../testing/mtex_parity_matrix.md)
- Markley et al., *Averaging Quaternions* (2007),
  [doi:10.2514/1.28949](https://doi.org/10.2514/1.28949).
