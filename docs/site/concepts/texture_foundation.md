# Texture Foundation

This page drills into the texture subsystem as a scientific layer on top of the canonical core.

## What The Texture Layer Owns

- rotation and misorientation semantics
- pole figures and inverse pole figures
- symmetry-aware reduction into inverse-pole-figure sectors
- orientation-space reduction and fundamental-region keys
- kernel-based ODF evaluation
- discrete pole-figure inversion over an explicit orientation dictionary
- band-limited harmonic ODF reconstruction over the same canonical texture objects
- IPF color-key generation

## Texture Flow

:::{figure} ../../figures/texture_foundation_flow.svg
:alt: Texture foundation flow from canonical core semantics through rotations, symmetry reduction, PF/IPF, ODF, harmonic reconstruction, and teaching outputs.
:class: architecture-poster-figure
:::

## Why This Layer Matters

The texture layer is where PyTex demonstrates that its canonical core is not abstract governance. It has to produce scientifically useful outputs:

- orientations must reduce correctly under symmetry
- pole figures must reflect the chosen crystal directions and specimen directions
- ODF evaluation must remain deterministic and interpretable
- pole-figure inversion must stay explicit about its dictionary, kernel, harmonic bandlimit, and convergence report
- IPF color mappings must stay tied to explicit symmetry semantics

## Current State

- rotation import/export is implemented
- symmetry-aware misorientation and disorientation are implemented
- PF/IPF and ODF foundations are implemented
- discrete dictionary-based PF inversion is implemented
- band-limited harmonic PF-to-ODF inversion is implemented
- exhaustive orientation-region boundaries and broader experimental correction doctrine remain ahead

## Related Material

- {doc}`orientation_texture`
- {doc}`../architecture/orientation_and_texture_foundation`
- {doc}`../standards/reference_canon`

## References

### Normative

- {doc}`../architecture/orientation_and_texture_foundation`
- {doc}`../standards/reference_canon`

### Informative

- {doc}`../workflows/ipf_colors`
- {doc}`symmetry_and_fundamental_regions`
