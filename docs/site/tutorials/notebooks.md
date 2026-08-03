# Tutorial Notebooks

PyTex now treats executable notebooks as a first-class tutorial layer inside the documentation system.

These notebooks are meant to do two things at once:

- teach the scientific meaning behind the APIs
- provide runnable end-to-end examples that users can adapt for their own work

The notebooks are intentionally aligned with the architecture and theory notes rather than acting as informal side material.

They are also expected to track the current runtime API closely. The notebooks are hand-authored
and committed *executed*; when plotting, inversion, or batch behavior changes in code, the
affected notebooks must be re-run (`python scripts/execute_notebooks.py --only <prefix>`) in the
same change so the executable examples stay trustworthy. Every notebook computes its results live
and checks them against analytic or literature values, so a divergence between code and tutorial
surfaces as a test failure rather than a silent documentation lie.

The priority roadmap notebooks are smoke-executed in the default suite when they stay within the
lightweight surface, and the heavier structure or diffraction notebooks are controlled by the full
scientific lane. PyTex therefore treats them as stable teaching artifacts rather than as unverified
examples.

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
notebooks/17_miller_vectorized_workflows
notebooks/18_orientation_relationships_fundamentals
notebooks/19_lattice_correspondence_and_transformation_strain
notebooks/20_or_catalogs_identification_and_reconstruction
notebooks/21_composite_or_diffraction_patterns
notebooks/22_burgers_beta_to_alpha_zirconium
notebooks/23_transformation_crystallography_end_to_end
```

## Orientation-Relationship Teaching Track

The orientation-relationship notebooks are committed *executed*: their pole
figures, spectra, and reports render directly in this site. They pair with the
scientific diagrams under `_static/or/` and the OR concept and foundation
pages.

1. `18_orientation_relationships_fundamentals` — conventions, the KS
   relationship, variants, intervariant spectrum, packets, variant pole
   figures, OR deviation and fitting.
2. `19_lattice_correspondence_and_transformation_strain` — index
   correspondence with rationalized indices and residuals, the Bain strain
   computed and rendered in 3D, KS/NW/GT polar rotations, and the Burgers
   basal-on-(110) overlay computed from the actual rotation.
3. `20_or_catalogs_identification_and_reconstruction` — the standard catalogs
   with computed variant counts and literature separations, then the
   experimental pipeline end to end: OR identification from boundaries,
   rotation refinement, and parent-grain reconstruction.
4. `21_composite_or_diffraction_patterns` — the diffraction observable of an
   OR: the Ewald/excitation-error theory behind the kinematic engine, the
   shared parent-anchored detector geometry, and composite SAED patterns for
   both canonical cases — Kurdjumov-Sachs (cubic, with the KS-NW 5.26 deg
   separation appearing as a child-zone deviation) and Burgers beta->alpha
   (hexagonal, with the six-fold basal view, four-index Miller-Bravais
   labels, and the {110}_beta / (0002)_alpha superposition at 0.1545 mm).
5. `22_burgers_beta_to_alpha_zirconium` — the single-system deep dive, on the
   alloy Burgers himself worked on. Builds the whole relationship from the
   distorted hexagon on {110}_beta, then derives every result twice: the
   12 variants from group theory, the five-valued intervariant spectrum
   against the literature table, all three principal strains in closed form,
   the half-integer correspondence that exposes the shuffle, OR-placed 3D
   unit cells, and composite SAED down [110], [111], [001] and [112]_beta.
   Shows why zirconium is not titanium with different numbers: the
   {110}_beta / (0002)_alpha coincidence splits about eight times wider.
6. `23_transformation_crystallography_end_to_end` — the five questions a
   transformation study actually asks, answered on Burgers beta->alpha in one
   pass: the orientation relationship determined from measured Euler angles with
   no nominal supplied, the variant correspondence table for an arbitrary parent
   plane, the composite SAED down [110]_beta with its reflection table and
   manifest, the same composite re-anchored on [0001]_alpha of one variant, and
   a measured pattern solved back to phase, zone axis, spot indices and variant.

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

## Stable Teaching Sequence

For contributors and new users, the default reading and execution path during the current hardening
phase is:

1. `01_reference_frames_and_transforms`
2. `02_rotations_orientations_and_batch_primitives`
3. `04_phases_lattices_space_groups_and_cif`
4. `05_multimodal_acquisition_and_manifests`
5. `06_texture_odf_and_pole_figure_inversion`
6. `07_ebsd_regular_grid_workflows`
7. `08_diffraction_geometry_and_kinematic_spots`
8. `10_plotting_semantic_primitives`

That path mirrors the intended public learning route: canonical frames and orientations first,
structure and provenance second, then texture, EBSD, diffraction, and plotting.

For the validated fixture-to-diffraction route in the full scientific lane, continue with:

9. `11_powder_xrd_workflows`
10. `12_saed_workflows`
11. `13_crystal_visualization_workflows`
12. `15_structure_diffraction_visualization_pipeline`

## How This Relates To Validation

- The lightweight notebooks are smoke-executed in the base lane, and the fixture-backed
  structure/diffraction notebooks are controlled by the full scientific lane.
- Notebook examples are expected to agree with the public concept, workflow, and validation pages.
- If a notebook shows a pedagogical simplification, the corresponding workflow page should say so
  explicitly rather than leaving the distinction implicit in code cells.

## References

### Normative

- {doc}`../standards/documentation_architecture`
- {doc}`../standards/development_principles`

### Informative

- {doc}`../architecture/overview`
