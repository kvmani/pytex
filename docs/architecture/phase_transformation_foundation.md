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

- The primitive family now exists, but variant-generation validation is still weak.
- No dedicated transformation manifest schema exists yet; current workflow-result and validation manifests must carry that context until a richer transformation contract is justified.
- Parent reconstruction should remain outside the stable algorithmic surface until the primitive family has stronger literature-backed validation.

## References

### Normative

- `../standards/reference_canon.md`
- `canonical_data_model.md`

### Informative

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*
