# Tutorial Notebooks

PyTex now treats executable notebooks as a first-class tutorial layer inside the documentation system.

These notebooks are meant to do two things at once:

- teach the scientific meaning behind the APIs
- provide runnable end-to-end examples that users can adapt for their own work

The notebooks are intentionally aligned with the architecture and theory notes rather than acting as informal side material.

They are also expected to track the current runtime API closely. When plotting, inversion, or
batch behavior changes in code, the notebook generator and the built notebooks must be updated in
the same change so the executable examples stay trustworthy.

The priority roadmap notebooks are also smoke-executed in the default test suite. PyTex therefore
treats them as stable teaching artifacts rather than as unverified examples.

## Notebook Atlas

```{toctree}
:maxdepth: 1

notebooks/01_reference_frames_and_transforms
notebooks/02_rotations_orientations_and_batch_primitives
notebooks/03_symmetry_and_fundamental_regions
notebooks/04_phases_lattices_space_groups_and_cif
notebooks/05_multimodal_acquisition_and_manifests
notebooks/06_texture_odf_and_pole_figure_inversion
notebooks/07_ebsd_regular_grid_workflows
notebooks/08_diffraction_geometry_and_kinematic_spots
notebooks/09_phase_transformation_foundations
notebooks/10_plotting_semantic_primitives
notebooks/11_powder_xrd_workflows
notebooks/12_saed_workflows
notebooks/13_crystal_visualization_workflows
notebooks/14_yaml_style_customization
notebooks/15_structure_diffraction_visualization_pipeline
notebooks/16_ebsd_to_texture_outputs
```

## How To Use These

- Read them in the built Sphinx site when you want concept-plus-code explanation.
- Open the raw `.ipynb` files locally when you want to execute and modify the examples.
- Use the texture and plotting notebooks to see the currently implemented contour pole-figure and
  Bunge-section ODF surfaces, not just older scatter-style examples.
- Use the XRD, SAED, crystal-visualization, and style notebooks to see the current diffraction and
  structure-view surfaces rather than relying on stale standalone scripts.
- Use the linked concept pages and LaTeX notes when you need the deeper formal derivation behind the tutorial flow.

## Immediate Roadmap Path

The main near-term teaching path is now organized around the pinned fixture corpus and the same
validation artifacts used elsewhere in the repository:

1. `04_phases_lattices_space_groups_and_cif`
2. `13_crystal_visualization_workflows`
3. `11_powder_xrd_workflows`
4. `12_saed_workflows`
5. `15_structure_diffraction_visualization_pipeline`

That sequence takes one pinned phase from CIF-backed construction through structure visualization,
diffraction generation, and manifest-backed reproducibility notes.

## References

### Normative

- `../../standards/documentation_architecture.md`
- `../../standards/development_principles.md`

### Informative

- `../../architecture/overview.md`
