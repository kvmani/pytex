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
   $(\mathbf{u}, \mathbf{v}, \mathbf{z})$ is anchored in the **parent** crystal frame.
2. **Variant mapping.** Each variant rotation $\mathbf{V}_i$ re-expresses
   parent-crystal-frame Cartesian vectors in the child crystal frame. The same
   physical beam direction becomes $\mathbf{z}_c = \mathbf{V}_i \mathbf{z}_p$, and
   each child pattern is simulated on the rotated basis
   $\mathbf{V}_i (\mathbf{u}, \mathbf{v}, \mathbf{z})$ — algebraically identical to
   pulling child reciprocal vectors back to the parent frame before projection.
   All sub-patterns therefore share one detector.
3. **Reflection selection.** A reflection $\mathbf{g}$ is kept when its excitation
   error $s_g = g_z - g^{2}\lambda / 2$ satisfies $|s_g| \le s_{\max}$, where $g_z$
   is the zone-axis component of $\mathbf{g}$ and $\lambda$ is the relativistic
   electron wavelength. This small-angle kinematic criterion treats rational parent
   zones and irrational child zones uniformly and stays honest about Ewald-sphere
   curvature ($s_g = -g^{2}\lambda / 2 \le 0$ for exact zero-order-Laue-zone spots).
4. **Intensity.** $I \propto |F_{hkl}|^{2}$ from the atomic-number electron
   structure-factor proxy (with isotropic Debye-Waller damping), optionally
   relrod-damped, with lattice-centering systematic absences applied. Each
   sub-pattern is max-normalized; kinematic cross-phase intensity ratios are
   undefined at this level of theory and are a rendering choice.
5. **Detector map.** $r_{\mathrm{mm}} = \Lambda\, g_{\perp}$ with camera constant
   $\Lambda$, the standard SAED
   small-angle scale between reciprocal length and detector millimeters.

The relativistic wavelength

$$
\lambda(V) = \frac{h}{\sqrt{2 m_0 e V \left(1 + \dfrac{eV}{2 m_0 c^{2}}\right)}}
$$

reproduces the standard tabulated values ($0.02508\ \text{\AA}$ at $200\ \mathrm{kV}$); see the
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
- **Axes units.** Render in calibrated detector $\mathrm{mm}$ or reciprocal
  $\text{\AA}^{-1}$.
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
zone that lands exactly on $\langle 111 \rangle$ reports $0^{\circ}$ while an off-zone variant
reports its true tilt.

## Anchoring On A Product Zone Instead

`simulate_composite_saed` takes the **parent** zone axis, which matches how the
crystallography is derived but not always how the microscope is used. In
practice the operator tilts to a low-index zone of the *product* — say
$[0001]$ of one alpha variant — and then wants to know what the matrix and the
other variants contribute to that same pattern.

```python
pattern = simulate_composite_saed_from_child_zone(
    burgers,
    ZoneAxis(np.array([0, 0, 1]), phase=alpha),
    anchor_variant_index=3,
    variant_indices=(1, 2, 3),
)
```

The anchor variant's rotation $\mathbf{R}_k$ carries parent Cartesian vectors
into that child's frame, so the requested child zone $\mathbf{z}_c$ corresponds
to the parent direction $\mathbf{R}_k^{\mathsf{T}} \mathbf{z}_c$. That direction
is **generally irrational**, and the result reports it exactly alongside its
nearest rational label — the same honesty child zone axes already get in the
parent-anchored case, so neither crystal is privileged in the output.

The detector basis is then built from that parent direction through the same
`zone_basis_from_axis` call the parent-anchored path uses, so there is exactly
one geometry definition. The consequence is a testable identity: **anchoring on
variant $k$'s image of a parent zone reproduces the parent-anchored pattern for
that zone exactly.** `align_child_g` works in the child's own indices and is
mapped to the parent frame internally.

`variant_indices` chooses which variants are drawn and is independent of
`anchor_variant_index`; the anchor variant is not added automatically.

Every export states which crystal defined the geometry: the pattern carries
`anchor_variant_index`, `describe()` and `anchor_description()` name it, and the
manifest records both it and the nearest rational parent zone.

### A note on spot order

Symmetry-equivalent reflections have mathematically equal intensity and detector
radius that differ in the last few ULPs depending on how the basis was built.
The engine's sort quantizes those two continuous keys (1 pm of radius, $10^{-12}$
of full-scale intensity) before ordering, so ties fall through to the exact
lexicographic `hkl` comparison rather than to floating-point noise. Without
that, the two anchoring routes produced correctly-positioned but *permuted*
spot tables, and any exported table or pinned figure inherited the permutation.

## Exporting The Pattern

A simulated pattern usually has to leave the process — as a table for a paper, a
figure for a slide, and a record of how both were produced.
`pytex.diffraction.export` is that boundary.

`composite_reflection_table(pattern)` returns one row per rendered spot carrying
its source (parent, or variant $k$), phase, Miller indices and formatted label,
$d$-spacing, $|g|$, detector position and radius in millimetres, excitation
error, structure-factor amplitude and relative intensity. Every value is read
from the `SpotTable` objects the engine produced, so the table and the figure
cannot disagree. It exports through `to_csv`, `to_markdown`, `to_records` and
`to_json_dict`.

Two things about those numbers are worth stating plainly, and `describe()`
states both:

- **Intensities are max-normalized within each sub-pattern separately.** Compare
  intensities within one source, never across sources: kinematic theory defines
  no intensity ratio between two different phases.
- **The detector radius uses the in-plane part of $\mathbf{g}$**, so it is
  slightly smaller than $(L\lambda)\,|\mathbf{g}|$. The difference is the
  out-of-plane component that the excitation error records; $d = 1/|g|$ uses the
  full vector.

`export_composite_saed(pattern, directory)` writes the table, the rendered
figure in the requested formats, the parent/child coincidence table, and a JSON
manifest validated by `schemas/composite_saed_manifest.schema.json` — the
relationship, both phases with their applied lattice centering, every variant's
exact and nearest-rational child zone axis, and the full simulation
configuration. Figures are closed after writing, so calling it in a loop leaks
nothing.

### The centering trap

`ReflectionCondition.from_phase` reads the lattice centering from the first
letter of a phase's space-group symbol, and falls back to primitive when the
phase carries none. A body-centred phase supplied **without** that metadata is
therefore simulated as primitive, and its pattern shows reflections the real
structure forbids — with nothing in the spot list to say so.

`pattern.centering_audit()` reports, per phase, the centering applied and
whether it was *declared* or *assumed*; `describe()` and the manifest carry the
same statement, and the reflection table's `describe()` raises a warning when
anything was assumed. If a simulated bcc pattern shows a $\{100\}$ reflection,
this is why.

## Current Limits

- kinematic intensities only — no dynamical (Bloch-wave / multi-beam) model,
  no double-diffraction excitation of kinematically forbidden spots
- no higher-order-Laue-zone (HOLZ) ring construction yet
- cross-phase intensity scaling is a rendering choice, not physics

## Stating The Parent Frame In The Pattern

Every sub-pattern of a composite figure shares one parent-anchored detector basis, so the detector
axes are trivially the page axes and say nothing useful. What a reader actually needs is where the
**parent crystal** axes land on that detector, which the pattern's zone basis supplies:

```python
render_composite_saed(
    pattern, config=CompositeSAEDPlotConfig(show_frame_indicator=True)
)
```

See {doc}`../architecture/reference_frame_foundation`.

## Related Material

- {doc}`../tutorials/notebooks/21_composite_or_diffraction_patterns` — the
  executed teaching notebook: Ewald/excitation-error theory, the shared
  detector construction, and both canonical cases end to end
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
