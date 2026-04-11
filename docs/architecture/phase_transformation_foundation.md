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

Variant generation must eventually state:

- which symmetry groups are acting
- which equivalence relation defines a distinct variant
- how variants are indexed and reproduced
- what happens when multiple literature conventions exist

Until that doctrine is implemented, transformation code should remain experimental.

## Habit Planes And Direction Correspondence

If PyTex exposes habit-plane or direction-correspondence features, the stable API must keep direct-versus-reciprocal meaning explicit and link them to the parent and child phase semantics already in the core model.

The current stable constructor is intentionally narrow: one plane-normal correspondence plus one
in-plane direction correspondence define a right-handed parent-to-child crystal mapping. Broader
literature family catalogs and ambiguity handling remain future work.

PyTex now also includes a small named-helper layer on top of the explicit correspondence
constructor, starting with `OrientationRelationship.from_bain_correspondence(...)` and
`OrientationRelationship.from_nishiyama_wassermann_correspondence(...)`. The goal is explicitness:
each builder encodes one stated correspondence rather than hiding an unnamed matrix.

## Current Limits

- The primitive family now exists and has a dedicated manifest schema, but literature-backed
  validation is still foundational rather than broad.
- Relationship construction from one explicit plane-direction correspondence is now implemented, but
  broader named literature families and competing convention catalogs remain ahead.
- A small named cubic transformation-helper layer now exists, but this is still the start of the
  catalog, not the catalog itself.
- Variant generation and variant-indexed prediction are now benchmarked and validated in-repo, but
  broader transformation-family coverage remains ahead.
- Parent reconstruction remains outside the stable algorithmic surface. PyTex now stages bounded
  candidate-parent scoring under `pytex.experimental` so research workflows can proceed without
  overstating stability.

## References

### Normative

- [Reference Canon](../standards/reference_canon.md)
- [Canonical Data Model](canonical_data_model.md)

### Informative

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*
