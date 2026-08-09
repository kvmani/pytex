# API Guide

PyTex keeps the stable API centered on named scientific primitives. This page is intentionally curated rather than exhaustive: the goal is to show the public objects users should build around, then point back to the concept and workflow pages that explain their meaning.

```{toctree}
:maxdepth: 2

full_reference
```

## How To Read This Page

This is not a replacement for the concept and workflow docs.

Use it in this order:

1. identify the scientific object family you need
2. read the linked concept or workflow page that defines its meaning
3. come back here to find the stable constructor or helper surface

If a type below would be ambiguous to you without knowing its frame, symmetry, provenance, or
reduction rules, do not start from the API list alone. Follow the linked concept page first.

If you want to see how these types *fit together* rather than what each one is, read
{doc}`../architecture/class_model_atlas`. It carries generated UML-style diagrams of the class
hierarchy and of the composition relations between the objects listed below, one per domain.

## Recommended Entry Points

- Start with `ReferenceFrame`, `SymmetrySpec`, `Phase`, and `Orientation` if you are building core
  crystallographic objects from scratch.
- Start with `Phase.from_cif(...)` and the phase or diffraction workflows if your data begins from
  structure definitions.
- Start with `normalize_ebsd(...)` and `CrystalMap` if your data begins from EBSD tooling.
- Start with `DiffractionGeometry`, `RadiationSpec`, `PowderPattern`, or `SAEDPattern` if your
  work begins in detector or reciprocal-space reasoning.
- Start with the manifest family when workflow context, provenance, validation, or interchange must
  survive beyond one Python call boundary.

## Current Interpretation Rule

The stable PyTex surface prefers semantically explicit objects over raw arrays. In practice that
means:

- use domain types when frame or symmetry meaning matters
- use semantic batch types when vectorized data shares one scientific interpretation
- use manifests or JSON contracts when results must remain reconstructible outside in-memory use

## Core

### Reference Frames

See {doc}`../architecture/reference_frame_foundation` and
{doc}`../concepts/reference_frames_and_conventions` for the model these build on.

- `ReferenceFrame` - named, domain-typed frame carrying axis labels and axis geometry
- `FrameTransform` - typed, invertible, composable map between two frames
- `FrameGraph` - registry that resolves the transform between any two connected frames
- Catalog constants: `CARTESIAN_FRAME`, `SPECIMEN_FRAME`, `SAMPLE_RD_TD_ND_FRAME`,
  `CRYSTAL_FRAME`, `MAP_FRAME`, `DETECTOR_FRAME`, `LABORATORY_FRAME`
- Catalog builders: `cartesian_frame`, `specimen_frame`, `sample_frame`, `crystal_frame`,
  `map_frame`, `detector_frame`, `laboratory_frame`, `reciprocal_frame_for`,
  `rolling_frame_graph`, `get_standard_frame`, `list_standard_frames`

### Other Core Primitives

- `VectorSet`
- `AcquisitionGeometry`
- `CalibrationRecord`
- `MeasurementQuality`
- `SymmetrySpec`
- `IPFSectorBoundary`
- `OrientationFundamentalRegion`
- `SpaceGroupSpec`
- `Lattice`
- `AtomicSite`
- `UnitCell`
- `Phase`
- `MillerPlane`
- `MillerDirection`
- `MillerPlaneSet`
- `MillerDirectionSet`
- `reduce_indices(...)`
- `canonicalize_sign(...)`
- `antipodal_keys(...)`
- `angle_plane_plane_rad(...)`
- `angle_dir_dir_rad(...)`
- `project_directions_onto_planes(...)`
- `ReciprocalLatticeVector`
- `ZoneAxis`
- `EulerSet`
- `EulerConventionTransform`
- `QuaternionSet`
- `Rotation`
- `RotationSet`
- `Orientation`
- `Orientation.from_euler(...)`
- `Orientation.from_axis_angle(...)`
- `Orientation.from_matrix(...)`
- `Orientation.from_quaternion(...)`
- `Orientation.from_miller(...)`
- `OrientationSet`
- `OrientationSet.from_axes_angles(...)`
- `OrientationSet.from_equispaced_so3_grid(...)`
- `OrientationSet.from_matrices(...)`
- `OrientationSet.from_plane_direction(...)`
- `OrientationSet.from_quaternions(...)`
- `OrientationSet.from_regular_so3_grid(...)`
- `OrientationSet.from_so2_grid(...)`
- `specimen_direction_vector(...)`
- `ScatteringSetup`
- `format_miller_indices`
- `format_plane_indices`
- `format_direction_indices`
- `OrientationRelationship`
- `OrientationRelationshipCatalog`
- `OrientationRelationship.from_parallel_plane_direction(...)`
- `OrientationRelationship.from_bain_correspondence(...)`
- `OrientationRelationship.from_nishiyama_wassermann_correspondence(...)`
- `TransformationVariant`
- `PhaseTransformationRecord`
- `ParentReconstructionConfig`
- `ParentReconstructionReport`
- `VariantSelectionReport`
- `reconstruct_parent_orientation(...)`

See {doc}`../concepts/core_model`, {doc}`../concepts/how_pytex_differs`,
{doc}`../concepts/miller_planes_directions`, and {doc}`../concepts/orientation_texture`.

## Texture

- `PoleFigure`
- `InversePoleFigure`
- `ODF`
- `HarmonicODF`
- `ODFInversionReport`
- `HarmonicODFReconstructionReport`
- `ODFReconstructionConfig`
- `PoleFigureCorrectionSpec`
- `PoleFigureResidualReport`
- `residual_reports_for_pole_figures(...)`
- `KernelSpec`
- `XRDMLPoleFigureMeasurement`
- `LaboTexPoleFigureMeasurement`
- `read_xrdml_pole_figure`
- `load_xrdml_pole_figure`
- `invert_xrdml_pole_figures`
- `read_labotex_pole_figures`
- `load_labotex_pole_figures`
- `invert_labotex_pole_figures`

See {doc}`../concepts/orientation_texture`.

## EBSD

- `CrystalMap`
- `CrystalMapPhase`
- `CoordinateNeighborGraph`
- `Grain`
- `GrainSegmentation`
- `GrainBoundarySegment`
- `GrainBoundaryNetwork`
- `GrainGraph`
- `EBSDImportManifest`
- `ExperimentManifest`
- `BenchmarkManifest`
- `TransformationManifest`
- `ValidationManifest`
- `WorkflowResultManifest`
- `NormalizedEBSDDataset`
- `TextureReport`
- `OrientationQualityWeights`
- `EBSDTextureWorkflow`
- `EBSDTextureWorkflowResult`
- `normalize_ebsd(...)`
- `index_hough(...)`
- `refine_orientations(...)`

Adapter-boundary utilities:

- `to_orix_phase(...)`
- `to_orix_miller_plane(...)`
- `to_orix_miller_direction(...)`
- `from_orix_miller(...)`
- `to_orix_rotation(...)`
- `from_orix_rotation(...)`
- `to_orix_orientation(...)`
- `from_orix_orientation(...)`

See {doc}`../workflows/ebsd_kam`, {doc}`../workflows/ebsd_grains`, and {doc}`../workflows/ebsd_to_texture_outputs`.

Those adapter helpers are intentionally grouped as boundary utilities rather than as the center of
the EBSD API. They preserve PyTex semantics at the edge of optional ORIX or KikuchiPy
interoperability, but they should not be read as a blanket parity claim for the external packages.

## Experimental

The following surfaces are intentionally outside the stable API contract but are documented because
they are useful for research workflows:

- `pytex.experimental.score_parent_orientations(...)`
- `pytex.experimental.ParentReconstructionResult`

See {doc}`../workflows/phase_transformation_manifests_and_scoring`.

## Diffraction

- `DetectorAcceptanceMask`
- `DiffractionGeometry`
- `DiffractionPattern`
- `RadiationSpec`
- `ScatteringFactorTable`
- `StructureFactor`
- `ReflectionCondition`
- `DiffractionIntensityModel`
- `PowderReflection`
- `PowderPattern`
- `SAEDSpot`
- `SAEDPattern`
- `KinematicSpot`
- `KinematicSimulation`
- `ReflectionFamily`
- `IndexingCandidate`
- `FamilyIndexingReport`
- `OrientationRefinementResult`

See {doc}`../workflows/diffraction_geometry`, {doc}`../workflows/diffraction_spots`,
{doc}`../workflows/xrd_generation`, and {doc}`../workflows/saed_generation`.

## Plotting

### Reference-Frame Visualization

- `frame_triad(...)` / `frame_triad_primitives(...)` - triad primitives for 3D scenes
- `add_frame_indicator(axes, frame, ...)` - embeddable corner gizmo for any 2D figure
  (diffractograms, pole figures, IPF maps, crystal-viewer panels)
- `plot_reference_frame(...)` / `plot_frame_relationship(...)` - standalone 3D figures
- `reference_frame_svg(...)` / `frame_catalog_svg(...)` - documentation SVG, no matplotlib
- `project_orthographic(...)`, `FrameTriad`, `TRIAD_AXIS_COLORS`

### Other Plotting Surfaces

- `IPFColorKey`
- `ipf_color(...)`
- `ipf_colors(...)`
- `CrystalCellOverlay`
- `CrystalPlaneOverlay`
- `CrystalDirectionOverlay`
- `plot_vector_set`
- `plot_wulff_net`
- `plot_crystal_directions`
- `plot_crystal_planes`
- `plot_symmetry_orbit`
- `plot_symmetry_elements`
- `plot_euler_set`
- `plot_quaternion_set`
- `plot_rotations`
- `plot_orientations`
- `plot_pole_figure`
- `plot_inverse_pole_figure`
- `plot_ipf_map`
- `plot_kam_map`
- `plot_odf`
- `plot_xrd_pattern`
- `plot_saed_pattern`
- `CrystalScene`
- `build_crystal_scene`
- `plot_crystal_structure_3d`
- `list_style_themes`
- `load_style_theme`
- `read_style_yaml`
- `resolve_style`
- `save_documentation_figure_svg`

See {doc}`../workflows/ipf_colors`, {doc}`../workflows/plotting_primitives`,
{doc}`../workflows/stereographic_projections`, {doc}`../workflows/crystal_visualization`,
{doc}`../workflows/ebsd_to_texture_outputs`, and {doc}`../workflows/style_customization`.

`plot_pole_figure(...)` supports scatter, histogram, and contour rendering. `plot_odf(...)`
supports scatter, contour, and classical Bunge-section rendering for both the discrete
`ODF` surface and the harmonic `HarmonicODF` surface.

`plot_xrd_pattern(...)` and `plot_saed_pattern(...)` return ordinary Matplotlib figures using the
shared YAML style system. `plot_crystal_structure_3d(...)` provides a publication-oriented 3D
structure view while preserving PyTex lattice and plane semantics, including optional unit-cell and
hexagonal-prism overlays where scientifically appropriate.

For architectural context, see {doc}`../architecture/overview`.

## Limits Of This Guide

- This page is curated, not exhaustive API documentation.
- The exhaustive module, class, method, and function reference is available at {doc}`full_reference`.
- The structural view — which class holds which, and what inherits from what — is
  {doc}`../architecture/class_model_atlas`.
- A symbol appearing here does not mean every downstream workflow built on it is equally validated.
- Use {doc}`../validation/index` to check current parity, evidence, and limitations before relying
  on a surface for stronger scientific claims.
