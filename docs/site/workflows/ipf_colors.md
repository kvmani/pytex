# Texture: IPF Color Keys

PyTex now includes a first stable plotting-side surface for inverse pole figure color assignment
through `IPFColorKey`, plus convenience helpers `ipf_color(...)` and `ipf_colors(...)`.

## Scope

- explicit IPF color-key objects rather than hidden plotting defaults
- symmetry-aware reduction of crystal directions into the supported IPF sector
- deterministic RGB assignment from sector-vertex barycentric weights
- color generation directly from `Orientation`, `OrientationSet`, or crystal directions
- named specimen directions through `RD`, `TD`, `ND`, `x`, `y`, and `z`

## Why This Exists

IPF coloring is not just decoration. It is a compact encoding of orientation information, and its meaning depends on:

- which specimen direction is being colored
- which crystal symmetry is active
- whether antipodal identification is in effect
- how the chosen IPF sector is mapped into color space

PyTex keeps those choices explicit by making the color key a named public object.

## Example

```python
import numpy as np

from pytex import IPFColorKey, ipf_colors

color_key = IPFColorKey(
    crystal_symmetry=orientations.symmetry,
    specimen_direction="ND",
)

rgb = color_key.colors_from_orientations(orientations)
rgb_direct = ipf_colors(orientations, direction="ND")
```

## Current Interpretation

The current implementation is a teaching-grade and workflow-grade color key:

- directions are reduced into the supported IPF sector
- sector vertices act as canonical color anchors
- RGB values are built from barycentric weights inside the sector cone

This keeps the semantics explicit and deterministic, even though it does not yet claim full MTEX color-key parity.

## Current Limits

- no full MTEX-parity IPF color-key reproduction yet
- IPF rendering now exists, but richer map-style plotting and publication presets are still ahead
- color-key semantics are ahead of full MTEX visual parity

## Related Material

- {doc}`../concepts/orientation_texture`
- {doc}`plotting_primitives`
- `docs/testing/mtex_parity_matrix.md`
