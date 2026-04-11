# API Guide

PyTex keeps the stable API centered on named scientific primitives. This page is intentionally curated rather than exhaustive: the goal is to show the public objects users should build around, then point back to the concept and workflow pages that explain their meaning.

## How To Read This Page

This is not a replacement for the concept and workflow docs.

Use it in this order:

1. identify the scientific object family you need
2. read the linked concept or workflow page that defines its meaning
3. come back here to find the stable constructor or helper surface

If a type below would be ambiguous to you without knowing its frame, symmetry, provenance, or
reduction rules, do not start from the API list alone. Follow the linked concept page first.

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
- `OrientationRelationship.from_parallel_plane_direction(...)`
- `OrientationRelationship.from_bain_correspondence(...)`
- `OrientationRelationship.from_nishiyama_wassermann_correspondence(...)`
- `TransformationVariant`
- `PhaseTransformationRecord`

See {doc}`../concepts/core_model`, {doc}`../concepts/how_pytex_differs`,
{doc}`../concepts/miller_planes_directions`, and {doc}`../concepts/orientation_texture`.

## Texture

- `PoleFigure`
- `InversePoleFigure`
- `ODF`
- `HarmonicODF`
- `ODFInversionReport`
- `HarmonicODFReconstructionReport`
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
supports scatter, contour, and classical Bunge-section rendering for both the discrete
`ODF` surface and the harmonic `HarmonicODF` surface.

`plot_xrd_pattern(...)` and `plot_saed_pattern(...)` return ordinary Matplotlib figures using the
shared YAML style system. `plot_crystal_structure_3d(...)` provides a publication-oriented 3D
structure view while preserving PyTex lattice and plane semantics, including optional unit-cell and
hexagonal-prism overlays where scientifically appropriate.

For architectural context, see {doc}`../architecture/overview`.

## Limits Of This Guide

- This page is curated, not exhaustive API documentation.
- A symbol appearing here does not mean every downstream workflow built on it is equally validated.
- Use {doc}`../validation/index` to check current parity, evidence, and limitations before relying
  on a surface for stronger scientific claims.
