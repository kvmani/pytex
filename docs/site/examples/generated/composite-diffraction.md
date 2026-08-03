<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Composite OR diffraction

Numerical cornerstones of composite orientation-relationship SAED simulation: the relativistic electron wavelength against the standard 200 kV value, the exactness of the Kurdjumov-Sachs child-zone mapping, and the two defining Burgers beta->alpha signatures (exact basal zone and the {110}/(0002) near-coincidence), plus the identities the exported reflection table must satisfy.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Relativistic electron wavelength at 200 kV

Every kinematic TEM computation starts from the electron wavelength, which fixes the Ewald-sphere radius k = 1/lambda and hence every excitation error. The relativistic formula lambda = h / sqrt(2 m0 e V (1 + e V / (2 m0 c^2))) must reproduce the standard tabulated value at a 200 kV accelerating voltage.

**Symbols**

- $\lambda$ &mdash; Radiation wavelength.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.diffraction.kinematic import electron_wavelength_angstrom
```

:::

**Compute**

```python
result = electron_wavelength_angstrom(200.0)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `composite-electron-wavelength-200kv` | 0.02508 | 0.02508 | angstrom | 6.60e-07 | 5e-06 | ✅ pass |

**Why this value**: The standard relativistic electron wavelength at 200 kV is 2.508 pm = 0.02508 angstrom.

**Citation**: De Graef, Introduction to Conventional Transmission Electron Microscopy, Cambridge University Press, 2003, Table 2.2.

**See also**: {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## KS maps the parent [0 1 -1] zone exactly onto a <1 1 1> child zone

The Kurdjumov-Sachs relationship is defined by the parallelism <-1 0 1>_fcc || <-1 -1 1>_bcc. When the composite SAED simulator maps a parent [0 1 -1] zone axis (a member of the <-1 0 1> family) through all 24 variants, at least one variant's child zone axis must land exactly on a <1 1 1>-type direction: the minimal angular deviation between mapped and rational child zones over the variants is zero.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import simulate_composite_saed

parent_frame = ReferenceFrame(
    "parent_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
child_frame = ReferenceFrame(
    "child_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
austenite = Phase(
    "austenite",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
    crystal_frame=parent_frame,
)
martensite = Phase(
    "martensite",
    lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
    crystal_frame=child_frame,
)
ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=martensite
)
```

:::

**Compute**

```python
zone = ZoneAxis(np.array([0, 1, -1]), phase=austenite)
composite = simulate_composite_saed(ks, zone, include_parent=False)
result = min(
    pattern.nearest_zone_axis.deviation_deg
    for pattern in composite.variant_patterns
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `composite-ks-exact-child-zone` | 0.0000 | 0.0000 | deg | 9.00e-15 | 1e-09 | ✅ pass |

**Why this value**: The defining KS direction parallelism makes the mapped child zone rational, so the deviation of the best variant is exactly 0 degrees.

**Citation**: Kurdjumov and Sachs, Z. Physik 64 (1930) 325; Morito et al., Acta Materialia 51 (2003) 1789 (variant conventions).

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## Burgers maps the parent <110> zone exactly onto the hcp [0001] basal zone

The Burgers relationship governing the beta->alpha transformation of titanium, zirconium and hafnium is defined by the plane parallelism {110}_bcc || (0001)_hcp. Viewing a beta crystal down a <110> zone axis must therefore look straight down the hcp c-axis for the variants whose basal plane is that particular {110}: the minimal angular deviation between the mapped child zone and a rational [0001] zone must be exactly zero.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import simulate_composite_saed

beta_frame = ReferenceFrame(
    "beta_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
alpha_frame = ReferenceFrame(
    "alpha_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
# Beta-titanium (bcc) and alpha-titanium (hcp), room-temperature parameters.
beta_ti = Phase(
    "beta-titanium",
    lattice=Lattice(3.3065, 3.3065, 3.3065, 90.0, 90.0, 90.0, crystal_frame=beta_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=beta_frame),
    crystal_frame=beta_frame,
    space_group_symbol="Im-3m",
)
alpha_ti = Phase(
    "alpha-titanium",
    lattice=Lattice(2.9508, 2.9508, 4.6855, 90.0, 90.0, 120.0, crystal_frame=alpha_frame),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=alpha_frame),
    crystal_frame=alpha_frame,
    space_group_symbol="P6_3/mmc",
)
burgers = OrientationRelationship.from_burgers_correspondence(
    parent_phase=beta_ti, child_phase=alpha_ti
)
```

:::

**Compute**

```python
zone = ZoneAxis(np.array([1, 1, 0]), phase=beta_ti)
composite = simulate_composite_saed(burgers, zone, include_parent=False)
result = min(
    pattern.nearest_zone_axis.deviation_deg
    for pattern in composite.variant_patterns
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `composite-burgers-exact-basal-zone` | 0.0000 | 0.0000 | deg | 4.80e-15 | 1e-09 | ✅ pass |

**Why this value**: The defining Burgers plane parallelism {110}_bcc || (0001)_hcp makes the mapped child zone exactly rational, so the deviation of the best variant is 0 degrees.

**Citation**: Burgers, Physica 1 (1934) 561.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Composite OR diffraction workflow <../../workflows/composite_or_diffraction>`

## Burgers {110}_bcc and (0002)_hcp reflections nearly superimpose

The practical TEM signature of the Burgers relationship is that the beta {110} reflection lands almost exactly on the alpha (0002) reflection, because the plane parallelism pairs two nearly equal interplanar spacings: d(110)_bcc = a/sqrt(2) = 2.3381 angstrom against d(0002)_hcp = c/2 = 2.3428 angstrom. At a 180 mm*angstrom camera constant the residual detector separation is well under a spot diameter, so the composite pattern reads as a single decorated pattern. This computes that separation from the simulated composite.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import simulate_composite_saed

beta_frame = ReferenceFrame(
    "beta_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
alpha_frame = ReferenceFrame(
    "alpha_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
# Beta-titanium (bcc) and alpha-titanium (hcp), room-temperature parameters.
beta_ti = Phase(
    "beta-titanium",
    lattice=Lattice(3.3065, 3.3065, 3.3065, 90.0, 90.0, 90.0, crystal_frame=beta_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=beta_frame),
    crystal_frame=beta_frame,
    space_group_symbol="Im-3m",
)
alpha_ti = Phase(
    "alpha-titanium",
    lattice=Lattice(2.9508, 2.9508, 4.6855, 90.0, 90.0, 120.0, crystal_frame=alpha_frame),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=alpha_frame),
    crystal_frame=alpha_frame,
    space_group_symbol="P6_3/mmc",
)
burgers = OrientationRelationship.from_burgers_correspondence(
    parent_phase=beta_ti, child_phase=alpha_ti
)
```

:::

**Compute**

```python
from pytex.diffraction.composite import find_spot_coincidences

zone = ZoneAxis(np.array([1, 1, 0]), phase=beta_ti)
composite = simulate_composite_saed(burgers, zone)
report = find_spot_coincidences(composite, tolerance_mm=1.0)
result = report.coincidences[0].separation_mm
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `composite-burgers-110-0002-coincidence` | 0.15450 | 0.15450 | mm | 2.02e-06 | 1e-04 | ✅ pass |

**Why this value**: Analytically the separation is (sqrt(2)/a_bcc - 2/c_hcp) * camera_constant = (1.414214/3.3065 - 2/4.6855) * 180 = 0.15450 mm.

**Citation**: Burgers, Physica 1 (1934) 561; lattice parameters from standard Ti data.

**See also**: {doc}`Composite OR diffraction workflow <../../workflows/composite_or_diffraction>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## The exported reflection table obeys d = 1/|g|, body centring, and the {110} spacing

Before an exported reflection table can serve as a measurement reference, it must satisfy the identities its own columns imply. This tabulates a Burgers composite viewed along beta [110] and checks four things: that the table lists exactly the pattern's own spots and no others, that every row's d-spacing is the reciprocal of its reported |g| rather than a separately computed quantity that could drift, that no body-centring-forbidden beta reflection survived, and that the strongest beta reflection is the {110} whose spacing is a / sqrt(2).

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import simulate_composite_saed

beta_frame = ReferenceFrame(
    "beta_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
alpha_frame = ReferenceFrame(
    "alpha_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
# Beta-titanium (bcc) and alpha-titanium (hcp), room-temperature parameters.
beta_ti = Phase(
    "beta-titanium",
    lattice=Lattice(3.3065, 3.3065, 3.3065, 90.0, 90.0, 90.0, crystal_frame=beta_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=beta_frame),
    crystal_frame=beta_frame,
    space_group_symbol="Im-3m",
)
alpha_ti = Phase(
    "alpha-titanium",
    lattice=Lattice(2.9508, 2.9508, 4.6855, 90.0, 90.0, 120.0, crystal_frame=alpha_frame),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=alpha_frame),
    crystal_frame=alpha_frame,
    space_group_symbol="P6_3/mmc",
)
burgers = OrientationRelationship.from_burgers_correspondence(
    parent_phase=beta_ti, child_phase=alpha_ti
)
```

:::

**Compute**

```python
from pytex.diffraction.export import composite_reflection_table

zone = ZoneAxis(np.array([1, 1, 0]), phase=beta_ti)
composite = simulate_composite_saed(burgers, zone, variant_indices=(1,))
table = composite_reflection_table(composite)
parent_rows = table.rows_for_source("parent")
# The strongest beta reflection along [110] is a {110}: d = a / sqrt(2).
strongest = parent_rows[0]
result = [
    # The table is a view of the pattern, so it lists every spot and no others.
    len(table) - composite.spot_count(),
    # d and |g| are one quantity reported two ways.
    max(abs(row.d_angstrom - 1.0 / row.g_inv_angstrom) for row in table.rows),
    # Body centring forbids h + k + l odd, so no such beta row may survive.
    sum(1 for row in parent_rows if sum(row.hkl) % 2 != 0),
    strongest.d_angstrom,
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `composite-burgers-reflection-table-identities` | [0.0000, 0.0000, 0.0000, 2.3380] | [0.0000, 0.0000, 0.0000, 2.3380] | counts, angstrom, counts, angstrom | 1.43e-06 | 1e-05 | ✅ pass |

**Why this value**: The first two are identities that hold to machine precision: the table is a view of the simulated pattern, and d = 1/|g| is a definition, so any nonzero value is an export defect. The third is the body-centred reflection condition h + k + l = 2n (International Tables Vol. A), which the beta phase's Im-3m space group imposes. The fourth is the analytic bcc {110} interplanar spacing a / sqrt(2) = 3.3065 / 1.414214 for the standard beta-Ti lattice parameter.

**Citation**: Burgers, Physica 1 (1934) 561; lattice parameters from standard Ti data. Reflection conditions and interplanar spacings: International Tables for Crystallography, Vol. A.

**See also**: {doc}`Composite OR diffraction workflow <../../workflows/composite_or_diffraction>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`

## Anchoring on a variant's own zone reproduces the parent-anchored pattern

A composite pattern can be set up two ways: choose the parent zone axis, or choose a zone axis of one product variant and let the parent direction follow. The two must agree, because the anchor variant's rotation R_k satisfies R_k^T (R_k z_p) = z_p, so both routes build the detector basis about the same parent direction. This simulates a Burgers composite along beta [110], re-anchors it on variant 2's own view of that zone, and measures the largest detector displacement of any spot — for the variants and for the parent.

:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    ZoneAxis,
)
from pytex.core.transformation import OrientationRelationship
from pytex.diffraction.composite import simulate_composite_saed

beta_frame = ReferenceFrame(
    "beta_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
alpha_frame = ReferenceFrame(
    "alpha_crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT
)
# Beta-titanium (bcc) and alpha-titanium (hcp), room-temperature parameters.
beta_ti = Phase(
    "beta-titanium",
    lattice=Lattice(3.3065, 3.3065, 3.3065, 90.0, 90.0, 90.0, crystal_frame=beta_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=beta_frame),
    crystal_frame=beta_frame,
    space_group_symbol="Im-3m",
)
alpha_ti = Phase(
    "alpha-titanium",
    lattice=Lattice(2.9508, 2.9508, 4.6855, 90.0, 90.0, 120.0, crystal_frame=alpha_frame),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=alpha_frame),
    crystal_frame=alpha_frame,
    space_group_symbol="P6_3/mmc",
)
burgers = OrientationRelationship.from_burgers_correspondence(
    parent_phase=beta_ti, child_phase=alpha_ti
)
```

:::

**Compute**

```python
from pytex.diffraction.composite import simulate_composite_saed_from_child_zone

selection = (1, 2, 3, 4)
parent_zone = ZoneAxis(np.array([1, 1, 0]), phase=beta_ti)
reference = simulate_composite_saed(burgers, parent_zone, variant_indices=selection)
# Anchor on variant 2's own view of that same parent zone.
recovered = simulate_composite_saed_from_child_zone(
    burgers,
    reference.variant_pattern(2).zone_axis_child,
    anchor_variant_index=2,
    variant_indices=selection,
)
largest_shift = max(
    float(
        np.max(
            np.abs(
                reference.variant_pattern(index).spots.detector_mm
                - recovered.variant_pattern(index).spots.detector_mm
            )
        )
    )
    for index in selection
)
parent_shift = float(
    np.max(
        np.abs(
            reference.parent_spots.detector_mm - recovered.parent_spots.detector_mm
        )
    )
)
result = [largest_shift, parent_shift]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `composite-child-anchored-geometry-consistency` | [0.0000, 0.0000] | [0.0000, 0.0000] | mm | 5.68e-14 | 1e-09 | ✅ pass |

**Why this value**: An exact identity of the construction, not a measured agreement: the child-anchored entry point maps the requested child zone back through R_k^T and then delegates to the parent-anchored engine, so the shared detector basis is the same object built the same way. Any nonzero displacement would mean the two paths had diverged. The 1e-9 mm tolerance is the floating-point round trip through the rotation, not a physical margin.

**Citation**: Burgers, Physica 1 (1934) 561.

**See also**: {doc}`Composite OR diffraction workflow <../../workflows/composite_or_diffraction>`, {doc}`Diffraction foundation <../../concepts/diffraction_foundation>`
