# Documentation Index

## Core Documents

- [Mission](../mission.md)
- [Changelog](../CHANGELOG.md)
- [Specifications](../specifications.md)
- [Agent Instructions](../AGENTS.md)

## Architecture

- [Architecture Overview](architecture/overview.md)
- [Class And Object Model Atlas](site/architecture/class_model_atlas.md) — generated class hierarchy and per-domain object-model diagrams
- [Canonical Data Model](architecture/canonical_data_model.md)
- [Reference Frame Foundation](architecture/reference_frame_foundation.md)
- [Orientation And Texture Foundation](architecture/orientation_and_texture_foundation.md)
- [EBSD Foundation](architecture/ebsd_foundation.md)
- [Diffraction Foundation](architecture/diffraction_foundation.md)
- [Multimodal Characterization Foundation](architecture/multimodal_characterization_foundation.md)
- [Phase Transformation Foundation](architecture/phase_transformation_foundation.md)
- [Orientation Relationship Analysis Foundation](architecture/orientation_relationship_analysis_foundation.md)
- [Transformation Crystallography And Composite Diffraction Program](architecture/transformation_crystallography_and_diffraction_program.md)
- [TEM Tilt Navigation Foundation](architecture/tem_tilt_navigation_foundation.md)
- [Repository Review 2026 Foundation Audit](architecture/repo_review_2026_foundation_audit.md)

## Testing And Validation

- [Testing Strategy](testing/strategy.md)
- [Automated Test Cases](testing/automated_test_cases.md)
- [MTEX Parity Matrix](testing/mtex_parity_matrix.md)
- [Diffraction Validation Matrix](testing/diffraction_validation_matrix.md)
- [Structure Validation Matrix](testing/structure_validation_matrix.md)
- [Plotting Validation Matrix](testing/plotting_validation_matrix.md)
- [VESTA Parity Matrix](testing/vesta_parity_matrix.md)
- [Phase Transformation Validation Matrix](testing/phase_transformation_validation_matrix.md)
- [Reconstruction Robustness Study](testing/reconstruction_robustness_study.md)

## Standards

- [Engineering Governance](standards/engineering_governance.md)
- [Notation And Conventions](standards/notation_and_conventions.md)
- [Documentation Architecture](standards/documentation_architecture.md)
- [Scientific Notes And Figures](standards/scientific_notes_and_figures.md)
- [Visualization Style Guide](standards/visualization_style_guide.md)
- [Terminology And Symbol Registry](standards/terminology_and_symbol_registry.md)
- [Scientific Citation Policy](standards/scientific_citation_policy.md)
- [Benchmark And Tolerance Governance](standards/benchmark_and_tolerance_governance.md)
- [Hexagonal And Trigonal Conventions](standards/hexagonal_and_trigonal_conventions.md)
- [Development Principles](standards/development_principles.md)
- [Data Contracts And Manifests](standards/data_contracts_and_manifests.md)
- [Reference Canon](standards/reference_canon.md)
- [Executable Examples](standards/executable_examples.md)
- [API Stability And Deprecation](standards/api_stability_and_deprecation.md)

## Development

- [Local Development](development/local_development.md)
- [Active Task Progress](development/active_task_progress.md) — the durable handoff record for the current task
- [Sphinx Site README](site/README.md)

## Roadmap

- [Critical Review And Development Guide](roadmap/critical_review_and_development_guide.md) — the governing development guide
- [Implementation Roadmap](roadmap/implementation_roadmap.md)
- [World-Class Feature And Foundation Roadmap](roadmap/world_class_feature_roadmap.md)
- [MTEX Parity And EBSD Feature Roadmap](roadmap/mtex_parity_and_ebsd_feature_roadmap.md)
- [Working Notes: Transformation Crystallography And Diffraction Program](roadmap/working_notes_transformation_diffraction_program.md) — the TX phase ledger
- [Working Notes: Algorithm Documentation Program](roadmap/working_notes_algorithm_documentation_program.md) — the active TD phase ledger

## Scientific Notes

Canonical theory, algorithm, and validation notes, authored as MyST Markdown and rendered in
full on the Sphinx site. See [Theory And Algorithm Notes](site/theory/index.md) for the grouped
index with cross-links to the matching implementation pages.

### Foundations

- [PyTex Project Philosophy](site/theory/project_philosophy.md)

### Theory

- [Reference Frames and Conventions in PyTex](site/theory/reference_frames.md)
- [Canonical Data Model for PyTex](site/theory/canonical_data_model.md)
- [Euler Convention Handling](site/theory/euler_convention_handling.md)
- [Orientation Representations And The Equal-Volume Maps](site/theory/orientation_representations.md)
- [Orientation Space, Symmetry Reduction, and Disorientation in PyTex](site/theory/orientation_space_and_disorientation.md)
- [Random Disorientation And The Mackenzie Baseline](site/theory/random_disorientation_baseline.md)
- [Fundamental Region Reduction](site/theory/fundamental_region_reduction.md)
- [Hexagonal and Trigonal Conventions in PyTex](site/theory/hexagonal_conventions.md)
- [Crystal Structures And CIF Import](site/theory/crystal_structures_and_cif_import.md)
- [Crystal Visualization Geometry](site/theory/crystal_visualization_geometry.md)

### Algorithms

- [Vectorized Miller Planes and Directions in PyTex](site/theory/vectorized_miller_planes_and_directions.md)
- [Orientation Representations and Plane–Direction Construction in PyTex](site/theory/orientation_representations_and_plane_direction_construction.md)
- [Discrete ODF and Pole-Figure Foundations in PyTex](site/theory/discrete_odf_and_pole_figures.md)
- [Harmonic ODF Reconstruction in PyTex](site/theory/harmonic_odf_reconstruction.md)
- [Stereographic Projections and XRDML Texture Import in PyTex](site/theory/stereographic_projections_and_xrdml_texture_import.md)
- [Preferred Orientation In Powder Intensities](site/theory/preferred_orientation_in_powder_intensities.md)
- [Inverse-Pole-Figure Colour Keys](site/theory/ipf_color_keys.md)
- [EBSD KAM Parameterization](site/theory/ebsd_kam_parameterization.md)
- [EBSD Local Misorientation Foundations](site/theory/ebsd_local_misorientation.md)
- [EBSD Grain Segmentation And GROD Foundations](site/theory/ebsd_grain_segmentation_and_grod.md)
- [EBSD Grain Boundaries And Cleanup Foundations](site/theory/ebsd_boundaries_and_cleanup.md)
- [Lattice Curvature And GND Density](site/theory/lattice_curvature_and_gnd_density.md)
- [Multiphase EBSD Graph Workflows](site/theory/multiphase_ebsd_graph_workflows.md)
- [Determining An Orientation Relationship From Measured Orientations](site/theory/orientation_relationship_determination.md)
- [Orientation-Relationship Index Correspondence](site/theory/orientation_relationship_index_correspondence.md)
- [Phase-Transformation Relationship Construction From Plane-Direction Correspondence](site/theory/phase_transformation_relationship_construction.md)
- [Experimental Parent Candidate Scoring](site/theory/experimental_parent_candidate_scoring.md)
- [Diffraction Geometry And Bragg Rings](site/theory/diffraction_geometry_and_bragg_rings.md)
- [Reciprocal Space And Kinematic Spots](site/theory/reciprocal_space_and_kinematic_spots.md)
- [Powder XRD And SAED Foundations](site/theory/powder_xrd_and_saed.md)
- [Ratio/Angle Indexing Of A Measured SAED Pattern](site/theory/saed_ratio_angle_indexing.md)
- [Kikuchi Bands And The Gnomonic Projection](site/theory/kikuchi_bands_and_gnomonic_projection.md)
- [Stereographic Kikuchi Maps And Zone-Axis Routing](site/theory/stereographic_kikuchi_maps.md)
- [TEM Specimen Tilt Navigation](site/theory/tem_specimen_tilt_navigation.md)
- [Convergent-Beam Electron Diffraction](site/theory/convergent_beam_electron_diffraction.md)
- [Dynamical CBED: Many-Beam Coupling, Absorption, HOLZ Lines, and Point-Group Determination](site/theory/dynamical_cbed_and_symmetry_determination.md)
- [Foundation Feature Priorities](site/theory/foundation_feature_priorities.md)

### Validation

- [PyTex Validation Program](site/theory/validation_program.md)

## Figures

- [Reference Frames](figures/reference_frames.svg)
- [Reference Frame Catalog](figures/reference_frame_catalog.svg) (generated by `scripts/generate_reference_frame_figures.py`)
- [Sample Frame RD TD ND](figures/sample_frame_rd_td_nd.svg) (generated by `scripts/generate_reference_frame_figures.py`)
- [Reference Frames Vectors](figures/reference_frames_vectors.svg)
- [Orientation Conventions](figures/orientation_conventions.svg)
- [Bunge Euler Geometry](figures/bunge_euler_geometry.svg)
- [Orientation Mapping Semantics](figures/orientation_mapping_semantics.svg)
- [Active Passive Rotation](figures/active_passive_rotation.svg)
- [Crystal Symmetry Actions](figures/crystal_symmetry_actions.svg)
- [Harmonic ODF Symmetry Projection](figures/harmonic_odf_symmetry_projection.svg)
- [IPF Sector Reduction](figures/ipf_sector_reduction.svg)
- [Disorientation Fundamental Region](figures/disorientation_fundamental_region.svg)
- [Diffraction Geometry](figures/diffraction_geometry.svg)
- [Zone Axis Ewald Geometry](figures/zone_axis_ewald_geometry.svg)
- [Kinematic Spot Projection](figures/kinematic_spot_projection.svg)
- [OR Determination Algorithm](figures/or_determination_algorithm.svg) (generated by `scripts/generate_algorithm_figures.py`)
- [Variant Correspondence Algorithm](figures/variant_correspondence_algorithm.svg) (generated by `scripts/generate_algorithm_figures.py`)
- [Composite SAED Algorithm](figures/composite_saed_algorithm.svg) (generated by `scripts/generate_algorithm_figures.py`)
- [SAED Indexing Algorithm](figures/saed_indexing_algorithm.svg) (generated by `scripts/generate_algorithm_figures.py`)
- [HCP Reference Frame](figures/hcp_reference_frame.svg)
- [Pole Figure Construction](figures/pole_figure_construction.svg)
- [PyTex System Structure](figures/pytex_system_structure.svg)
- [PyTex Scientific Data Flow](figures/pytex_scientific_data_flow.svg)
- [PyTex Governance And Completion Model](figures/pytex_governance_completion_model.svg)
- [PyTex Current State And Planned Expansion](figures/pytex_current_state_expansion.svg)
- [Core Foundation Map](figures/core_foundation_map.svg)
- [Texture Foundation Flow](figures/texture_foundation_flow.svg)
- [EBSD Foundation Flow](figures/ebsd_foundation_flow.svg)
- [Diffraction Foundation Flow](figures/diffraction_foundation_flow.svg)

## Documentation Rules

- Root Markdown documents provide discoverable guidance and links.
- When a document points to another repository document, use a clickable Markdown link rather than a plain backticked path.
- Sphinx is the primary user-facing documentation layer.
- Major scientific notes are authored canonically as MyST Markdown under `site/theory/`.
- Scientific geometry diagrams are maintained canonically as SVG.
- Stable features are not considered complete until docs, figures, examples, and validation notes all exist.
- Foundational conventions, frame mappings, symmetry reductions, major algorithms, and modality boundaries must be explained through prose, explicit mathematics, and annotated figures together.
