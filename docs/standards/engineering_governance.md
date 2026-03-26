# Engineering Governance

## Core Rules

- Stable APIs should expose scientific meaning explicitly.
- Major scientific claims require validation notes.
- Documentation changes are part of behavior changes.
- Major scientific conventions and algorithms must be documented with explicit mathematics and annotated graphics, not prose alone.
- Optional adapters must not weaken stable-core semantics.
- Code should prefer correctness and clarity before optimization.
- Normative conventions must be documented once and reused everywhere.
- Reference-frame, symmetry, and notation conversions must be centralized, not reimplemented ad hoc.
- Every major scientific document must declare its normative sources.
- Stable domain objects must reject inconsistent scientific metadata early.
- Schema-backed manifests are required before workflow interchange is considered stable.

## Repository Expectations

- Add module-local indexes when a subsystem grows materially.
- Keep docs synchronized with actual behavior.
- Use deterministic tests where possible.
- Do not add hidden or silent convention conversions.
- Do not depend on external scientific libraries for basic importability of the core package.
- When a subsystem depends on a scientific convention, update the corresponding standards document in the same change.
- When a change affects a foundational frame, symmetry, reduction, or algorithm contract, update the Sphinx concept or workflow page, the LaTeX note, and the canonical SVG figures in the same change.
