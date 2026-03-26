# Architecture Overview

PyTex is organized as a modular scientific library built around a canonical internal data model.

## Architectural Priorities

1. one shared model for frames, symmetry, and provenance
2. stable APIs that communicate scientific meaning explicitly
3. optional adapter boundaries for heavy external tools
4. documentation and validation artifacts treated as first-class outputs
5. foundational semantics explained through mathematics and graphics, not only code and prose

## Module Map

```text
src/pytex/
  core/
    conventions.py
    provenance.py
    frames.py
    symmetry.py
    lattice.py
    orientation.py
  texture/
    models.py
  ebsd/
    models.py
  diffraction/
    models.py
  adapters/
  plotting/
  experimental/
```

## Layering Rules

- `core/` owns canonical primitives and low-level math semantics.
- `texture/`, `ebsd/`, and `diffraction/` build on `core/`.
- `adapters/` may depend on external libraries, but stable core types must remain usable without them.
- `experimental/` can depend on stable types but must not weaken or bypass stable semantics.

## Data Flow Principle

Imports from vendors or external scientific libraries should be normalized into PyTex canonical primitives at the boundary. Internal algorithms should then work on PyTex types, not on source-specific representations.
