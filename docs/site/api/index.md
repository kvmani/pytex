# API Guide

PyTex keeps the stable API centered on named scientific primitives. This page is intentionally curated rather than exhaustive: the goal is to show the public objects users should build around, then point back to the concept and workflow pages that explain their meaning.

## Core

- `ReferenceFrame`
- `FrameTransform`
- `VectorSet`
- `AcquisitionGeometry`
- `CalibrationRecord`
- `MeasurementQuality`
- `SymmetrySpec`
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
- `QuaternionSet`
- `Rotation`
- `RotationSet`
- `Orientation`
- `OrientationSet`
- `OrientationSet.from_axes_angles(...)`
- `OrientationSet.from_matrices(...)`
- `OrientationSet.from_plane_direction(...)`
- `OrientationSet.from_quaternions(...)`
- `ScatteringSetup`
- `format_miller_indices`
- `format_plane_indices`
- `format_direction_indices`
- `OrientationRelationship`
- `TransformationVariant`
- `PhaseTransformationRecord`

See {doc}`../concepts/core_model`, {doc}`../concepts/how_pytex_differs`,
{doc}`../concepts/miller_planes_directions`, and {doc}`../concepts/orientation_texture`.

## Texture

- `PoleFigure`
- `InversePoleFigure`
- `ODF`
- `ODFInversionReport`
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
- `normalize_ebsd(...)`
- `index_hough(...)`
- `refine_orientations(...)`
- `to_orix_phase(...)`
- `to_orix_miller_plane(...)`
- `to_orix_miller_direction(...)`
- `from_orix_miller(...)`
- `to_orix_rotation(...)`
- `from_orix_rotation(...)`
- `to_orix_orientation(...)`
- `from_orix_orientation(...)`

See {doc}`../workflows/ebsd_kam`, {doc}`../workflows/ebsd_grains`, and {doc}`../workflows/ebsd_to_texture_outputs`.

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

- `IPFColorKey`
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
supports scatter, contour, and classical Bunge-section rendering.

`plot_xrd_pattern(...)` and `plot_saed_pattern(...)` return ordinary Matplotlib figures using the
shared YAML style system. `plot_crystal_structure_3d(...)` provides a publication-oriented 3D
structure view while preserving PyTex lattice and plane semantics, including optional unit-cell and
hexagonal-prism overlays where scientifically appropriate.

For architectural context, see `docs/architecture/overview.md`.
