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
- statistical sample symmetry is implemented for triclinic, monoclinic, orthorhombic and **axial**
  (fibre) specimens; axial symmetry is imposed on measured pole figures exactly, as a
  trapezoid-weighted ring average ({func}`~pytex.texture.sample_symmetry.impose_sample_symmetry`)
- ODF sections at constant $\varphi_2$, $\varphi_1$ or $\sigma$ are implemented for both ODF
  representations over the Euler box the crystal and sample symmetries require
  ({func}`~pytex.texture.sections.odf_sections`, {func}`~pytex.texture.sections.euler_section_ranges`)
- recalculated and difference pole figures, with the RP factor, are the acceptance test of every
  inversion in the workbench's measured-texture analysis
- component volume fractions integrate a continuous ODF over misorientation balls and are reported
  against their Haar random reference; a hexagonal ideal-orientation catalogue is provided
- the Kearns triad is computed from three measured theta-2theta section scans, with every peak,
  intensity and quadrature node reported — see {doc}`../workflows/texture_analysis_workbench`
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
