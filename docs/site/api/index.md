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
- `UnitCell`
- `Phase`
- `ReciprocalLatticeVector`
- `ZoneAxis`
- `EulerSet`
- `QuaternionSet`
- `Rotation`
- `RotationSet`
- `Orientation`
- `OrientationSet`
- `ScatteringSetup`
- `OrientationRelationship`
- `TransformationVariant`
- `PhaseTransformationRecord`

See {doc}`../concepts/core_model`, {doc}`../concepts/how_pytex_differs`, and {doc}`../concepts/orientation_texture`.

## Texture

- `PoleFigure`
- `InversePoleFigure`
- `ODF`
- `ODFInversionReport`
- `KernelSpec`

See {doc}`../concepts/orientation_texture`.

## EBSD

- `CrystalMap`
- `Grain`
- `GrainSegmentation`
- `GrainBoundarySegment`
- `GrainBoundaryNetwork`
- `GrainGraph`
- `EBSDImportManifest`
- `ExperimentManifest`
- `BenchmarkManifest`
- `ValidationManifest`
- `WorkflowResultManifest`
- `NormalizedEBSDDataset`

See {doc}`../workflows/ebsd_kam` and {doc}`../workflows/ebsd_grains`.

## Diffraction

- `DetectorAcceptanceMask`
- `DiffractionGeometry`
- `DiffractionPattern`
- `KinematicSpot`
- `KinematicSimulation`
- `ReflectionFamily`
- `IndexingCandidate`
- `FamilyIndexingReport`
- `OrientationRefinementResult`

See {doc}`../workflows/diffraction_geometry` and {doc}`../workflows/diffraction_spots`.

## Plotting

- `IPFColorKey`
- `plot_vector_set`
- `plot_symmetry_orbit`
- `plot_symmetry_elements`
- `plot_euler_set`
- `plot_quaternion_set`
- `plot_rotations`
- `plot_orientations`
- `plot_pole_figure`
- `plot_inverse_pole_figure`
- `plot_odf`
- `save_documentation_figure_svg`

See {doc}`../workflows/ipf_colors` and {doc}`../workflows/plotting_primitives`.

`plot_pole_figure(...)` supports scatter, histogram, and contour rendering. `plot_odf(...)`
supports scatter, contour, and classical Bunge-section rendering.

For architectural context, see `docs/architecture/overview.md`.
