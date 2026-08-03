# Phase Transformation Foundation

This document defines the architectural contracts that anchor PyTex phase-transformation workflows.

## Why This Exists

Phase transformation work cannot be bolted onto an orientation library by passing around pairs of phases and a few matrices. Stable transformation workflows need explicit semantics for:

- parent and child phase identities
- orientation relationships
- variant generation and indexing
- habit-plane and direction correspondences
- provenance for how the relationship was derived or selected
- explicit crystallographic correspondence when the relationship is derived from parallel planes and
  directions rather than only from a precomputed matrix

Without that foundation, transformation features would create private local semantics and fragment the core model.

## Stable Primitive Families

PyTex now treats the following as stable foundational primitives:

- `OrientationRelationship`
  Defines the named mapping between parent and child phase crystallographic objects.
- `TransformationVariant`
  Represents one generated or selected variant within a parent-child transformation family.
- `PhaseTransformationRecord`
  Stores provenance, assumptions, and workflow context for a transformation analysis.
- `TransformationManifest`
  Records dedicated machine-readable transformation workflow context rather than forcing that
  context into generic result metadata.

## Parent-Child Semantics

Stable transformation workflows must express:

- parent phase
- child phase
- the reference frames in which the relationship is defined
- whether the relationship is exact, fitted, literature-adopted, or inferred
- the source used to justify the relationship

The stable surface must not rely on unnamed arrays or undocumented variant numbering.

PyTex now also supports constructing an `OrientationRelationship` directly from one explicit
parent-plane / child-plane plus parent-direction / child-direction correspondence pair. That keeps
the crystallographic evidence attached to the relationship instead of forcing users to preserve only
the derived rotation matrix.

## Variant Doctrine

Variant generation states:

- which symmetry groups are acting (parent operators generate candidates; child operators define
  the equivalence orbit)
- which equivalence relation defines a distinct variant (child-symmetry orbit reduction by
  default, with the raw operator-product enumeration available explicitly)
- how variants are indexed and reproduced (enumeration order of
  `OrientationRelationship.generate_variants()`)

This doctrine is implemented and validated in `pytex.core.transformation` (literature-correct
variant counts for Bain, KS, NW, GT, Pitsch, Burgers). Alignment of variant numbering with
published conventions (e.g. Morito V1-V24) is tracked in the
[Orientation Relationship Analysis Foundation](orientation_relationship_analysis_foundation.md),
which is the normative program for all further OR-analysis work.

## Habit Planes And Direction Correspondence

If PyTex exposes habit-plane or direction-correspondence features, the stable API must keep direct-versus-reciprocal meaning explicit and link them to the parent and child phase semantics already in the core model.

The current stable constructor is intentionally narrow: one plane-normal correspondence plus one
in-plane direction correspondence define a right-handed parent-to-child crystal mapping. Broader
literature family catalogs and ambiguity handling remain future work.

PyTex now also includes a small named-helper layer on top of the explicit correspondence
constructor, starting with `OrientationRelationship.from_bain_correspondence(...)` and
`OrientationRelationship.from_nishiyama_wassermann_correspondence(...)`. The goal is explicitness:
each builder encodes one stated correspondence rather than hiding an unnamed matrix.

## Determining The Relationship From Measurements

A relationship is not always known in advance — usually it is what the measurement is *for*.
`characterize_orientation_relationship(parents, children)` (and its Euler-angle entry point)
takes paired parent/child orientations and returns the fitted rotation, the matching named
relationship, the parallel planes and directions that define it, and an explicit verdict on
whether the identification is trustworthy. `describe_orientation_relationship(...)` performs the
last step alone, turning any rotation back into the crystallographic statement the literature
uses.

This closes the loop with the constructors above: a relationship built from an explicit
correspondence records the parallelisms that define it, and those same parallelisms are what the
determination surface recovers and reports. The full specification is the
[Transformation Crystallography And Composite Diffraction Program](transformation_crystallography_and_diffraction_program.md),
feature TX1.

## Current Limits (updated 2026-07-17)

The full current-state statement lives in the
[Orientation Relationship Analysis Foundation](orientation_relationship_analysis_foundation.md)
(§1 implemented, §5 limits) — index correspondence, misorientation representation, deviation,
fitting, parallelism finders, packet classification, and variant pole figures are stable and
literature-pinned; map-scale parent-grain reconstruction is experimental pending
measured-data fixtures and MTEX parity; deformation gradients, habit-plane/PTMC analysis,
OR-from-boundaries determination, and the broader named-OR catalog remain future work.

## References

### Normative

- [Reference Canon](../standards/reference_canon.md)
- [Canonical Data Model](canonical_data_model.md)
- [Orientation Relationship Analysis Foundation](orientation_relationship_analysis_foundation.md)

### Informative

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*
