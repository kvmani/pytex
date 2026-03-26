# API Guide

PyTex keeps the stable API centered on named scientific primitives. This page is intentionally curated rather than exhaustive: the goal is to show the public objects users should build around, then point back to the concept and workflow pages that explain their meaning.

## Core

- `ReferenceFrame`
- `FrameTransform`
- `SymmetrySpec`
- `Lattice`
- `UnitCell`
- `Phase`
- `ReciprocalLatticeVector`
- `ZoneAxis`
- `Rotation`
- `Orientation`
- `OrientationSet`

See {doc}`../concepts/core_model`, {doc}`../concepts/how_pytex_differs`, and {doc}`../concepts/orientation_texture`.

## Texture

- `PoleFigure`
- `InversePoleFigure`
- `ODF`
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

See {doc}`../workflows/ipf_colors`.

For architectural context, see `docs/architecture/overview.md`.
