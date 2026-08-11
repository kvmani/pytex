# Theory And Algorithm Notes

These are the canonical scientific notes for PyTex: the derivations, the convention fixing, and
the assumptions and failure modes behind the implemented surfaces. Where the workflow pages
explain *how to get a result* and the {doc}`../algorithms/index` pages explain *how a result is
computed*, these notes explain *why the scientific contracts are defined the way they are*.

They are authored as MyST Markdown and render here in full, mathematics included. There is no
separate LaTeX source: a single source keeps the derivation on the page identical to the
derivation of record. A typeset PDF of this documentation is produced by Sphinx itself with
`sphinx -b latexpdf`, so the print and web forms cannot drift apart.

## How To Use This Section

- read the concept and workflow pages first if you need orientation to the library surface;
- read these notes when you need convention fixing, a derivation, or an algorithm's assumptions;
- treat these notes as canonical when a workflow page deliberately summarizes rather than repeats
  a full derivation.

Symbols and nomenclature follow {doc}`../standards/terminology_and_symbol_registry`, which is
normative for every note here.

## Foundations

```{toctree}
:maxdepth: 1

project_philosophy
```

## Theory

```{toctree}
:maxdepth: 1

reference_frames
canonical_data_model
euler_convention_handling
orientation_representations
orientation_space_and_disorientation
fundamental_region_reduction
hexagonal_conventions
crystal_structures_and_cif_import
crystal_visualization_geometry
```

## Algorithms

```{toctree}
:maxdepth: 1

vectorized_miller_planes_and_directions
orientation_representations_and_plane_direction_construction
discrete_odf_and_pole_figures
harmonic_odf_reconstruction
stereographic_projections_and_xrdml_texture_import
preferred_orientation_in_powder_intensities
ebsd_kam_parameterization
ebsd_local_misorientation
ebsd_grain_segmentation_and_grod
ebsd_boundaries_and_cleanup
lattice_curvature_and_gnd_density
multiphase_ebsd_graph_workflows
orientation_relationship_determination
orientation_relationship_index_correspondence
phase_transformation_relationship_construction
experimental_parent_candidate_scoring
diffraction_geometry_and_bragg_rings
reciprocal_space_and_kinematic_spots
powder_xrd_and_saed
saed_ratio_angle_indexing
kikuchi_bands_and_gnomonic_projection
stereographic_kikuchi_maps
tem_specimen_tilt_navigation
convergent_beam_electron_diffraction
dynamical_cbed_and_symmetry_determination
foundation_feature_priorities
```

## Validation

```{toctree}
:maxdepth: 1

validation_program
```

## Related Sections

Several algorithms are documented at two levels of detail, deliberately. The note here carries the
derivation; the {doc}`../algorithms/index` page carries the implementation as reimplementable
steps, the calibrated tolerances, the worked cubic and hexagonal numbers, and the failure modes.

| Note | Implementation page |
| --- | --- |
| {doc}`orientation_relationship_determination` | {doc}`../algorithms/orientation_relationship_determination` |
| {doc}`orientation_relationship_index_correspondence` | {doc}`../algorithms/variant_correspondence` |
| {doc}`saed_ratio_angle_indexing` | {doc}`../algorithms/saed_pattern_indexing` |
| {doc}`tem_specimen_tilt_navigation` | {doc}`../algorithms/tem_tilt_navigation` |
| {doc}`reciprocal_space_and_kinematic_spots` | {doc}`../algorithms/composite_saed_assembly` |
