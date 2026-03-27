# Texture Foundation

This page drills into the texture subsystem as a scientific layer on top of the canonical core.

## What The Texture Layer Owns

- rotation and misorientation semantics
- pole figures and inverse pole figures
- symmetry-aware reduction into inverse-pole-figure sectors
- orientation-space reduction and fundamental-region keys
- kernel-based ODF evaluation
- discrete pole-figure inversion over an explicit orientation dictionary
- IPF color-key generation

## Texture Flow

```{mermaid}
flowchart LR
    core["Canonical core<br/>Frames, symmetry, lattice, orientation"]
    rotations["Rotation and Orientation<br/>Quaternion-backed semantics"]
    symmetry["Symmetry reduction<br/>Fundamental sector and orbit logic"]
    pfipf["PF / IPF<br/>PoleFigure and InversePoleFigure"]
    odf["ODF<br/>Kernel evaluation and density queries"]
    colors["IPFColorKey<br/>Symmetry-aware presentation semantics"]

    core --> rotations
    rotations --> symmetry
    symmetry --> pfipf
    rotations --> pfipf
    symmetry --> odf
    symmetry --> colors
    rotations --> colors

    classDef core fill:#f0f4f8,stroke:#334e68,color:#102a43,stroke-width:1.5px;
    classDef domain fill:#d9f2e6,stroke:#2d6a4f,color:#1b4332,stroke-width:1.5px;
    class core core;
    class rotations,symmetry,pfipf,odf,colors domain;
```

## Why This Layer Matters

The texture layer is where PyTex demonstrates that its canonical core is not abstract governance. It has to produce scientifically useful outputs:

- orientations must reduce correctly under symmetry
- pole figures must reflect the chosen crystal directions and specimen directions
- ODF evaluation must remain deterministic and interpretable
- pole-figure inversion must stay explicit about its dictionary, kernel, and convergence report
- IPF color mappings must stay tied to explicit symmetry semantics

## Current State

- rotation import/export is implemented
- symmetry-aware misorientation and disorientation are implemented
- PF/IPF and ODF foundations are implemented
- discrete dictionary-based PF inversion is implemented
- exact harmonic inversion and exhaustive orientation-region boundaries remain ahead

## Related Material

- {doc}`orientation_texture`
- `docs/architecture/orientation_and_texture_foundation.md`
- `docs/standards/reference_canon.md`

## References

### Normative

- `docs/architecture/orientation_and_texture_foundation.md`
- `docs/standards/reference_canon.md`

### Informative

- `docs/site/workflows/ipf_colors.md`
- `docs/site/concepts/symmetry_and_fundamental_regions.md`
