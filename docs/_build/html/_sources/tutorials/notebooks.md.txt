# Tutorial Notebooks

PyTex now treats executable notebooks as a first-class tutorial layer inside the documentation system.

These notebooks are meant to do two things at once:

- teach the scientific meaning behind the APIs
- provide runnable end-to-end examples that users can adapt for their own work

The notebooks are intentionally aligned with the architecture and theory notes rather than acting as informal side material.

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
```

## How To Use These

- Read them in the built Sphinx site when you want concept-plus-code explanation.
- Open the raw `.ipynb` files locally when you want to execute and modify the examples.
- Use the linked concept pages and LaTeX notes when you need the deeper formal derivation behind the tutorial flow.

## References

### Normative

- `../../standards/documentation_architecture.md`
- `../../standards/development_principles.md`

### Informative

- `../../architecture/overview.md`
