# SAED Generation

PyTex now includes a kinematic selected-area electron diffraction workflow built on explicit
reciprocal-space and detector semantics.

![SAED Example](../../figures/saed_demo.svg)

## Scope

- reciprocal-lattice construction from `Phase`
- explicit `ZoneAxis` handling in crystal coordinates
- reflection filtering by zone condition
- detector-space projection through a camera-constant abstraction
- spot labeling and styling through the shared runtime plotting system

## Scientific Model

The current SAED workflow is a geometric and kinematic foundation:

1. enumerate candidate Miller indices
2. convert them into reciprocal-lattice vectors
3. apply the zone condition with the explicit direct-space zone axis
4. project in-zone reciprocal vectors into a detector basis orthogonal to the zone axis
5. assign a proxy intensity for ranking and plotting

The detector map is controlled by `camera_constant_mm_angstrom`, which acts as a simple
camera-length style scale factor between reciprocal-length units and detector millimeters.

## Example

```python
import numpy as np

from pytex import (
    FrameDomain,
    Handedness,
    ReferenceFrame,
    ZoneAxis,
    generate_saed_pattern,
    get_phase_fixture,
    plot_saed_pattern,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
phase = get_phase_fixture("ni_fcc").load_phase(crystal_frame=crystal)

pattern = generate_saed_pattern(
    phase,
    ZoneAxis(indices=np.array([0, 0, 1]), phase=phase),
    camera_constant_mm_angstrom=180.0,
    max_index=5,
    max_g_inv_angstrom=3.0,
)
figure = plot_saed_pattern(pattern, theme="dark")
figure.savefig("ni_fcc_saed.png", dpi=200)
```

## Coordinate Semantics

The current SAED workflow keeps three coordinate meanings separate:

- crystal direct-space coordinates for the `ZoneAxis`
- reciprocal-space coordinates for reflection construction
- detector-plane coordinates in millimeters for plotting

This is important because zone-axis reasoning is defined in direct space, while diffraction spots
live in reciprocal space and are finally rendered in detector coordinates.

`SAEDPattern` is the stable pattern-level container carrying the generated `SAEDSpot` collection,
named detector and reciprocal frames, the camera constant, and the crystal-basis information used
for the detector projection.

The first pinned external-baseline case for this workflow now uses the built-in `ni_fcc` fixture
for a `[001]` zone-axis pattern and records shell geometry against a `diffsims` reference result
under `fixtures/diffraction/`.

## Current Limits

- the current intensity is a proxy, not a dynamical diffraction model
- no Ewald-sphere curvature treatment for high-angle electron diffraction yet
- external-baseline coverage currently validates shell geometry for a pinned case rather than a
  broad orientation library

## Related Material

- {doc}`../concepts/technical_glossary_and_symbols`
- {doc}`phases_and_cif`
- {doc}`xrd_generation`
- {doc}`../tutorials/notebooks/12_saed_workflows`
- {doc}`style_customization`
- [../../tex/algorithms/powder_xrd_and_saed.tex](../../tex/algorithms/powder_xrd_and_saed.tex)

## References

### Normative

- `../../standards/reference_canon.md`
- `../../standards/notation_and_conventions.md`

### Informative

- `../../testing/diffraction_validation_matrix.md`
