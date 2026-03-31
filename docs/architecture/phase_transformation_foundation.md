# Phase Transformation Foundation

This document defines the architectural contracts that anchor PyTex phase-transformation workflows.

## Why This Exists

Phase transformation work cannot be bolted onto an orientation library by passing around pairs of phases and a few matrices. Stable transformation workflows need explicit semantics for:

- parent and child phase identities
- orientation relationships
- variant generation and indexing
- habit-plane and direction correspondences
- provenance for how the relationship was derived or selected

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

## Variant Doctrine

Variant generation must eventually state:

- which symmetry groups are acting
- which equivalence relation defines a distinct variant
- how variants are indexed and reproduced
- what happens when multiple literature conventions exist

Until that doctrine is implemented, transformation code should remain experimental.

## Habit Planes And Direction Correspondence

If PyTex exposes habit-plane or direction-correspondence features, the stable API must keep direct-versus-reciprocal meaning explicit and link them to the parent and child phase semantics already in the core model.

## Current Limits

- The primitive family now exists and has a dedicated manifest schema, but literature-backed
  validation is still foundational rather than broad.
- Variant generation and variant-indexed prediction are now benchmarked and validated in-repo, but
  broader transformation-family coverage remains ahead.
- Parent reconstruction remains outside the stable algorithmic surface. PyTex now stages bounded
  candidate-parent scoring under `pytex.experimental` so research workflows can proceed without
  overstating stability.

## References

### Normative

- [Reference Canon](../site/standards/reference_canon.md)
- [Canonical Data Model](../site/architecture/canonical_data_model.md)

### Informative

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*
