<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Orientation-relationship correspondence

Index-correspondence identities for named orientation relationships: mapping parent planes and directions to their product-phase counterparts, with rationalized indices and angular residuals.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Kurdjumov-Sachs maps (111) austenite onto (011) ferrite exactly

Given the Kurdjumov-Sachs relationship, find which ferrite plane corresponds to the austenite close-packed plane (111). Because {111}_fcc || {011}_bcc is the defining parallelism of the relationship, the mapped plane must rationalize to (011) with zero angular residual — the residual is the verifiable quantity.

**Symbols**

- $\mathbf{M}^{*}$ &mdash; Plane-index correspondence matrix mapping parent (hkl) to child (hkl).
- $(hkl)$ &mdash; Miller plane indices.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationRelationship,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

parent_frame = ReferenceFrame(
    name="austenite_crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
child_frame = ReferenceFrame(
    name="ferrite_crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
austenite = Phase(
    "austenite",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
    crystal_frame=parent_frame,
)
ferrite = Phase(
    "ferrite",
    lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
    crystal_frame=child_frame,
)
```

:::

**Compute**

```python
ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
mapped = ks.map_plane_to_child(
    CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=austenite), phase=austenite)
)
result = np.concatenate(
    [mapped.rational_indices.astype(float), [mapped.angular_residual_deg]]
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-ks-plane-correspondence-identity` | [0.0000, 1.0000, 1.0000, 0.0000] | [0.0000, 1.0000, 1.0000, 0.0000] | indices, deg | 3.79e-15 | 1e-09 | ✅ pass |

**Why this value**: The Kurdjumov-Sachs relationship is constructed from the parallelism {111}_fcc || {011}_bcc, so mapping the defining parent plane must recover the defining child plane identically (analytic identity).

**Citation**: Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## Bain maps [110] austenite onto [100] ferrite exactly

Given the Bain correspondence, find which ferrite direction corresponds to the austenite [110] direction. The Bain construction fixes [110]_fcc || [100]_bcc, so the mapped direction must rationalize to [100] with zero angular residual.

**Symbols**

- $\mathbf{M}$ &mdash; Direction-index correspondence matrix mapping parent [uvw] to child [uvw].
- $[uvw]$ &mdash; Miller direction indices.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    CrystalDirection,
    CrystalPlane,
    FrameDomain,
    Handedness,
    Lattice,
    MillerIndex,
    OrientationRelationship,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

parent_frame = ReferenceFrame(
    name="austenite_crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
child_frame = ReferenceFrame(
    name="ferrite_crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
austenite = Phase(
    "austenite",
    lattice=Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=parent_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=parent_frame),
    crystal_frame=parent_frame,
)
ferrite = Phase(
    "ferrite",
    lattice=Lattice(2.87, 2.87, 2.87, 90.0, 90.0, 90.0, crystal_frame=child_frame),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=child_frame),
    crystal_frame=child_frame,
)
```

:::

**Compute**

```python
bain = OrientationRelationship.from_bain_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
mapped = bain.map_direction_to_child(
    CrystalDirection([1.0, 1.0, 0.0], phase=austenite)
)
result = np.concatenate(
    [mapped.rational_indices.astype(float), [mapped.angular_residual_deg]]
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-bain-direction-correspondence-identity` | [1.0000, 0.0000, 0.0000, 0.0000] | [1.0000, 0.0000, 0.0000, 0.0000] | indices, deg | 1.03e-14 | 1e-09 | ✅ pass |

**Why this value**: The Bain correspondence is constructed from (001)_fcc || (001)_bcc with [110]_fcc || [100]_bcc, so mapping the defining parent direction must recover the defining child direction identically (analytic identity).

**Citation**: Bain, Trans. AIME 70 (1924) 25.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`
