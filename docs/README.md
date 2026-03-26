# Documentation Index

This folder is the navigation hub for implementation-facing and science-facing documentation.

## Core Documents

- `../mission.md`
- `../specifications.md`
- `../AGENTS.md`

## Architecture

- `architecture/overview.md`
- `architecture/canonical_data_model.md`
- `architecture/orientation_and_texture_foundation.md`
- `architecture/ebsd_foundation.md`
- `architecture/diffraction_foundation.md`

## Testing And Validation

- `testing/strategy.md`
- `testing/mtex_parity_matrix.md`

## Standards

- `standards/engineering_governance.md`
- `standards/notation_and_conventions.md`
- `standards/documentation_architecture.md`
- `standards/latex_and_figures.md`
- `standards/scientific_citation_policy.md`
- `standards/benchmark_and_tolerance_governance.md`
- `standards/hexagonal_and_trigonal_conventions.md`
- `standards/development_principles.md`
- `standards/data_contracts_and_manifests.md`

## Development

- `development/local_development.md`
- `site/README.md`

## Roadmap

- `roadmap/implementation_roadmap.md`

## LaTeX Source Tree

- `tex/README.md`
- `tex/foundations/project_philosophy.tex`
- `tex/theory/canonical_data_model.tex`
- `tex/theory/euler_convention_handling.tex`
- `tex/theory/fundamental_region_reduction.tex`
- `tex/theory/hexagonal_conventions.tex`
- `tex/theory/orientation_space_and_disorientation.tex`
- `tex/theory/reference_frames.tex`
- `tex/algorithms/discrete_odf_and_pole_figures.tex`
- `tex/algorithms/ebsd_kam_parameterization.tex`
- `tex/algorithms/ebsd_local_misorientation.tex`
- `tex/algorithms/ebsd_grain_segmentation_and_grod.tex`
- `tex/algorithms/ebsd_boundaries_and_cleanup.tex`
- `tex/algorithms/diffraction_geometry_and_bragg_rings.tex`
- `tex/algorithms/reciprocal_space_and_kinematic_spots.tex`
- `tex/validation/validation_program.tex`

## Figures

- `figures/reference_frames.svg`
- `figures/reference_frames_vectors.svg`
- `figures/orientation_conventions.svg`
- `figures/bunge_euler_geometry.svg`
- `figures/orientation_mapping_semantics.svg`
- `figures/active_passive_rotation.svg`
- `figures/crystal_symmetry_actions.svg`
- `figures/ipf_sector_reduction.svg`
- `figures/disorientation_fundamental_region.svg`
- `figures/diffraction_geometry.svg`
- `figures/zone_axis_ewald_geometry.svg`
- `figures/kinematic_spot_projection.svg`
- `figures/hcp_reference_frame.svg`
- `figures/pole_figure_construction.svg`

## Documentation Rules

- Root Markdown documents provide discoverable guidance and links.
- Sphinx is the primary user-facing documentation layer.
- Major scientific notes are authored canonically in LaTeX.
- Scientific geometry diagrams are maintained canonically as SVG.
- Stable features are not considered complete until docs, figures, examples, and validation notes all exist.
- Foundational conventions, frame mappings, symmetry reductions, and major algorithms must be explained through prose, explicit mathematics, and annotated figures together.
