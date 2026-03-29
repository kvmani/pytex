# Combined Structure, Diffraction, And Visualization Workflow

PyTex is designed so a single `Phase` can anchor structure visualization, powder diffraction,
SAED reasoning, and later multimodal workflows without redefining the crystallographic semantics
in each subsystem.

The immediate roadmap teaching path now uses the pinned `ni_fcc` fixture for this workflow, so the
same example is shared by docs, notebooks, structural plotting-validation cases, diffraction
baselines, and manifest-backed validation notes.

## Workflow Pattern

1. load or import a pinned `Phase`
2. visualize its unit cell or supercell
3. generate powder XRD reflections or a broadened spectrum
4. generate a SAED spot pattern for an explicit zone axis
5. render each output through the shared plotting and style system
6. trace the same case through workflow-result and validation manifests

## Why This Matters

This is the core PyTex architecture in action:

- one `Phase`
- one lattice and reciprocal basis
- one point-group reduction surface
- one shared plotting style layer

The structure does not need to be translated into a second private format just to move from a crystal viewer to an XRD or SAED workflow.

## Reproducibility Trail

The immediate roadmap path is also pinned to machine-readable artifacts:

- `benchmarks/structure_import/foundation_workflow_result_manifest.json`
- `benchmarks/diffraction/external_baseline_workflow_result_manifest.json`
- `benchmarks/validation/diffraction_validation_manifest.json`

Those files connect the fixture corpus, the first external diffraction baselines, and the
validation posture for the same teaching example.

## Related Material

- {doc}`phases_and_cif`
- {doc}`xrd_generation`
- {doc}`saed_generation`
- {doc}`crystal_visualization`
- {doc}`style_customization`
- {doc}`../tutorials/notebooks/15_structure_diffraction_visualization_pipeline`
- `../../testing/diffraction_validation_matrix.md`

## References

### Normative

- `../../architecture/canonical_data_model.md`
- `../../standards/data_contracts_and_manifests.md`

### Informative

- [../../tex/theory/crystal_visualization_geometry.tex](../../tex/theory/crystal_visualization_geometry.tex)
- [../../tex/algorithms/powder_xrd_and_saed.tex](../../tex/algorithms/powder_xrd_and_saed.tex)
