# Documentation Index

## Core Documents

- [Mission](../mission.md)
- [Specifications](../specifications.md)
- [Agent Instructions](../AGENTS.md)

## Architecture

- [Architecture Overview](architecture/overview.md)
- [Canonical Data Model](architecture/canonical_data_model.md)
- [Orientation And Texture Foundation](architecture/orientation_and_texture_foundation.md)
- [EBSD Foundation](architecture/ebsd_foundation.md)
- [Diffraction Foundation](architecture/diffraction_foundation.md)
- [Multimodal Characterization Foundation](architecture/multimodal_characterization_foundation.md)
- [Phase Transformation Foundation](architecture/phase_transformation_foundation.md)
- [Repository Review 2026 Foundation Audit](architecture/repo_review_2026_foundation_audit.md)

## Testing And Validation

- [Testing Strategy](testing/strategy.md)
- [MTEX Parity Matrix](testing/mtex_parity_matrix.md)
- [Diffraction Validation Matrix](testing/diffraction_validation_matrix.md)
- [Structure Validation Matrix](testing/structure_validation_matrix.md)
- [Plotting Validation Matrix](testing/plotting_validation_matrix.md)

## Standards

- [Engineering Governance](standards/engineering_governance.md)
- [Notation And Conventions](standards/notation_and_conventions.md)
- [Documentation Architecture](standards/documentation_architecture.md)
- [LaTeX And Figures](standards/latex_and_figures.md)
- [Terminology And Symbol Registry](standards/terminology_and_symbol_registry.md)
- [Scientific Citation Policy](standards/scientific_citation_policy.md)
- [Benchmark And Tolerance Governance](standards/benchmark_and_tolerance_governance.md)
- [Hexagonal And Trigonal Conventions](standards/hexagonal_and_trigonal_conventions.md)
- [Development Principles](standards/development_principles.md)
- [Data Contracts And Manifests](standards/data_contracts_and_manifests.md)
- [Reference Canon](standards/reference_canon.md)

## Development

- [Local Development](development/local_development.md)
- [Sphinx Site README](site/README.md)

## Roadmap

- [Implementation Roadmap](roadmap/implementation_roadmap.md)

## LaTeX Source Tree

- [LaTeX README](tex/README.md)
- [Project Philosophy](tex/foundations/project_philosophy.tex)
- [Canonical Data Model Theory Note](tex/theory/canonical_data_model.tex)
- [Euler Convention Handling](tex/theory/euler_convention_handling.tex)
- [Fundamental Region Reduction](tex/theory/fundamental_region_reduction.tex)
- [Hexagonal Conventions](tex/theory/hexagonal_conventions.tex)
- [Orientation Space And Disorientation](tex/theory/orientation_space_and_disorientation.tex)
- [Reference Frames](tex/theory/reference_frames.tex)
- [Crystal Visualization Geometry](tex/theory/crystal_visualization_geometry.tex)
- [Discrete ODF And Pole Figures](tex/algorithms/discrete_odf_and_pole_figures.tex)
- [EBSD KAM Parameterization](tex/algorithms/ebsd_kam_parameterization.tex)
- [EBSD Local Misorientation](tex/algorithms/ebsd_local_misorientation.tex)
- [EBSD Grain Segmentation And GROD](tex/algorithms/ebsd_grain_segmentation_and_grod.tex)
- [EBSD Boundaries And Cleanup](tex/algorithms/ebsd_boundaries_and_cleanup.tex)
- [Diffraction Geometry And Bragg Rings](tex/algorithms/diffraction_geometry_and_bragg_rings.tex)
- [Powder XRD And SAED](tex/algorithms/powder_xrd_and_saed.tex)
- [Reciprocal Space And Kinematic Spots](tex/algorithms/reciprocal_space_and_kinematic_spots.tex)
- [Validation Program](tex/validation/validation_program.tex)

## Figures

- [Reference Frames](figures/reference_frames.svg)
- [Reference Frames Vectors](figures/reference_frames_vectors.svg)
- [Orientation Conventions](figures/orientation_conventions.svg)
- [Bunge Euler Geometry](figures/bunge_euler_geometry.svg)
- [Orientation Mapping Semantics](figures/orientation_mapping_semantics.svg)
- [Active Passive Rotation](figures/active_passive_rotation.svg)
- [Crystal Symmetry Actions](figures/crystal_symmetry_actions.svg)
- [IPF Sector Reduction](figures/ipf_sector_reduction.svg)
- [Disorientation Fundamental Region](figures/disorientation_fundamental_region.svg)
- [Diffraction Geometry](figures/diffraction_geometry.svg)
- [Zone Axis Ewald Geometry](figures/zone_axis_ewald_geometry.svg)
- [Kinematic Spot Projection](figures/kinematic_spot_projection.svg)
- [HCP Reference Frame](figures/hcp_reference_frame.svg)
- [Pole Figure Construction](figures/pole_figure_construction.svg)

## Documentation Rules

- Root Markdown documents provide discoverable guidance and links.
- When a document points to another repository document, use a clickable Markdown link rather than a plain backticked path.
- Sphinx is the primary user-facing documentation layer.
- Major scientific notes are authored canonically in LaTeX.
- Scientific geometry diagrams are maintained canonically as SVG.
- Stable features are not considered complete until docs, figures, examples, and validation notes all exist.
- Foundational conventions, frame mappings, symmetry reductions, major algorithms, and modality boundaries must be explained through prose, explicit mathematics, and annotated figures together.
