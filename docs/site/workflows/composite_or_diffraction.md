# Composite OR Diffraction

PyTex simulates and renders **composite selected-area electron diffraction (SAED)
patterns** — a parent phase together with any subset of the transformation
variants of an orientation relationship (OR) — for an arbitrary parent zone
axis. This is the quantitative tool behind OR verification in the TEM: the
composite pattern shows exactly which parent and child reflections superimpose
and how the variants distribute in the diffraction plane.

The engine is fully vectorized and kinematic. It is a deliberate,
publication-oriented successor to the ad-hoc composite-pattern plotting many
groups maintain in scripts: configurable spot shapes, sizes, colors and
symbols; variant subsetting; in-plane rotation; and non-overlapping,
coincidence-merging spot annotation are all first-class, typed configuration.

## Scope

- vectorized kinematic zone-axis simulation ({func}`~pytex.diffraction.simulate_zone_axis_spots`)
- composite parent + OR-variant assembly for any parent zone axis
  ({func}`~pytex.diffraction.simulate_composite_saed`)
- a shared, parent-anchored detector frame so every phase and variant overlays
  one physically consistent pattern
- exact treatment of the (generally irrational) child zone axes, with
  nearest-rational `[uvw]` labels and their angular deviation
- quantitative spot-coincidence analysis
  ({func}`~pytex.diffraction.find_spot_coincidences`)
- layered, highly configurable publication rendering with crowding-aware
  annotations ({func}`~pytex.plotting.render_composite_saed`)

## Scientific Model

The workflow is strictly kinematic (no dynamical / multi-beam effects):

1. **Beam and zone axis.** The beam travels antiparallel to the zone axis; the
   zone-axis unit vector points toward the gun. The shared detector basis
   `(u, v, z)` is anchored in the **parent** crystal frame.
2. **Variant mapping.** Each variant rotation `V_i` re-expresses
   parent-crystal-frame Cartesian vectors in the child crystal frame. The same
   physical beam direction becomes `z_c = V_i z_p`, and each child pattern is
   simulated on the rotated basis `V_i (u, v, z)` — algebraically identical to
   pulling child reciprocal vectors back to the parent frame before projection.
   All sub-patterns therefore share one detector.
3. **Reflection selection.** A reflection `g` is kept when its excitation error
   `s_g = g_z - g^2 λ / 2` satisfies `|s_g| ≤ s_max`, where `g_z` is the
   zone-axis component of `g` and `λ` is the relativistic electron wavelength.
   This small-angle kinematic criterion treats rational parent zones and
   irrational child zones uniformly and stays honest about Ewald-sphere
   curvature (`s_g = -g^2 λ / 2 ≤ 0` for exact zero-order-Laue-zone spots).
4. **Intensity.** `I ∝ |F_hkl|²` from the atomic-number electron
   structure-factor proxy (with isotropic Debye-Waller damping), optionally
   relrod-damped, with lattice-centering systematic absences applied. Each
   sub-pattern is max-normalized; kinematic cross-phase intensity ratios are
   undefined at this level of theory and are a rendering choice.
5. **Detector map.** `r_mm = (camera constant) · g_⊥`, the standard SAED
   small-angle scale between reciprocal length and detector millimeters.

The relativistic wavelength `λ(V) = h / sqrt(2 m₀ e V (1 + eV / 2 m₀ c²))`
reproduces the standard tabulated values (0.02508 Å at 200 kV); see the
[composite diffraction worked examples](../examples/generated/composite-diffraction.md).

## Example

```python
import numpy as np

from pytex import (
    FrameDomain,
    Handedness,
    ReferenceFrame,
    ZoneAxis,
    get_phase_fixture,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction import find_spot_coincidences, simulate_composite_saed
from pytex.plotting import CompositeSAEDPlotConfig, render_composite_saed

parent_frame = ReferenceFrame("parent", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
child_frame = ReferenceFrame("child", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
austenite = get_phase_fixture("ni_fcc").load_phase(crystal_frame=parent_frame, phase_name="austenite")
ferrite = get_phase_fixture("fe_bcc").load_phase(crystal_frame=child_frame, phase_name="ferrite")

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)

composite = simulate_composite_saed(
    ks,
    ZoneAxis(np.array([0, 1, -1]), phase=austenite),
    variant_indices=(1, 2, 3, 4),
)

# Which child reflections superimpose on parent reflections?
report = find_spot_coincidences(composite, tolerance_mm=2.5)
print(report.describe())

figure = render_composite_saed(composite, config=CompositeSAEDPlotConfig())
figure.savefig("ks_composite.png", dpi=200, bbox_inches="tight")
```

`composite.describe()` and `report.describe()` produce convention-explicit,
citation-ready prose summaries of the geometry and the superposition counts.

## Configuration Highlights

- **Variant selection.** Pass `variant_indices` to
  {func}`~pytex.diffraction.simulate_composite_saed`, or restrict at render time
  with `CompositeSAEDPlotConfig(variant_indices=...)`; call
  `composite.select_variants(...)` to derive a subset pattern.
- **Spot styling.** {class}`~pytex.plotting.SpotStyle` controls marker shape,
  color, fill, edge, opacity, z-order and the intensity→size mapping
  (`intensity_area`, `intensity_radius` or `constant`). Variants cycle a
  colorblind-aware palette by default, or take explicit per-variant overrides.
- **In-plane rotation and alignment.** `in_plane_rotation_deg` rotates the whole
  composite; `align_parent_g` places a chosen parent reflection along `+u`.
- **Axes units.** Render in calibrated detector `mm` or reciprocal `Å⁻¹`.
- **Annotations.** {class}`~pytex.plotting.SpotAnnotationConfig` merges
  coincident reflections into phase-tagged multi-line labels and places labels
  greedily without overlapping each other or covering spots, with an intensity
  floor, a label budget, and optional leader lines. Pass
  `return_annotations=True` to receive an `AnnotationResult` report.

## Coordinate Semantics

The workflow keeps three coordinate meanings separate, exactly as the
single-phase {doc}`saed_generation` workflow does:

- crystal direct-space coordinates for the parent `ZoneAxis`
- reciprocal-space coordinates for reflection construction (in each phase)
- shared parent-anchored detector-plane coordinates (mm) for overlay and plotting

The child zone axes are stored as exact `CrystalDirection` objects;
{func}`~pytex.diffraction.rationalize_zone_axis` supplies nearest-integer
`[uvw]` labels with an honest angular deviation, so a Kurdjumov-Sachs child
zone that lands exactly on `⟨111⟩` reports `0°` while an off-zone variant
reports its true tilt.

## Current Limits

- kinematic intensities only — no dynamical (Bloch-wave / multi-beam) model,
  no double-diffraction excitation of kinematically forbidden spots
- no higher-order-Laue-zone (HOLZ) ring construction yet
- cross-phase intensity scaling is a rendering choice, not physics

## Related Material

- {doc}`saed_generation`
- {doc}`../concepts/orientation_relationships`
- {doc}`../concepts/diffraction_foundation`
- {doc}`phases_and_cif`
- {doc}`style_customization`
- [Composite diffraction worked examples](../examples/generated/composite-diffraction.md)

## References

### Normative

- `../../standards/reference_canon.md`
- `../../standards/notation_and_conventions.md`

### Informative

- De Graef, *Introduction to Conventional Transmission Electron Microscopy*,
  Cambridge University Press, 2003 (electron wavelength, kinematic SAED geometry).
- Kurdjumov and Sachs, *Z. Physik* **64** (1930) 325; Morito et al.,
  *Acta Materialia* **51** (2003) 1789 (variant conventions).
- `../../testing/diffraction_validation_matrix.md`
