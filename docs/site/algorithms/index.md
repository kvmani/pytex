# Algorithms

These pages state **how** each scientific surface computes what it computes: the
mathematics, the algorithm as steps a reader could reimplement, the constraints
and tolerances with what each is calibrated against, the complexity, and the
failure modes.

They sit between the other documentation layers rather than duplicating them:

| Layer | Answers | Example |
| --- | --- | --- |
| concepts | what the objects mean | {doc}`../concepts/orientation_relationships` |
| workflows | how to get a result | {doc}`../workflows/composite_or_diffraction` |
| **algorithms** (here) | how the result is computed, and what constrains it | this section |
| worked examples | the numbers, computed live and checked | {doc}`../examples/index` |
| theory notes | the canonical LaTeX derivations | {doc}`../theory/index` |

## Conventions used throughout

Every page works two systems side by side so nothing hexagonal is left implicit:

- **cubic** — Kurdjumov-Sachs, fcc → bcc (austenite → martensite);
- **hexagonal** — Burgers, bcc → hcp (β → α titanium and zirconium).

Angles are in degrees on every public surface, lengths on a detector in
millimetres, reciprocal lengths in Å⁻¹. Orientation matrices are crystal-to-
specimen in the Bunge convention, so a parent/child pair shows the rotation
$\mathbf{V} = \mathbf{C}^{\mathsf{T}}\mathbf{P}$. Symbols follow
{doc}`../standards/terminology_and_symbol_registry`.

Each figure on these pages is **generated** by
`scripts/generate_algorithm_figures.py` rather than drawn by hand, so a diagram
cannot drift from the algorithm it illustrates, and each is held to the
repository's figure layout guards.

```{toctree}
:maxdepth: 1

orientation_relationship_determination
variant_correspondence
composite_saed_assembly
saed_pattern_indexing
```
