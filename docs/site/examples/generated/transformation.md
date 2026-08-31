<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Orientation-relationship correspondence

Index-correspondence identities for named orientation relationships: mapping parent planes and directions to their product-phase counterparts, with rationalized indices and angular residuals, the misorientation representation used for EBSD comparison, and the recovery of a relationship and its parallelism statement from measured parent/child orientation pairs.

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
| `or-ks-plane-correspondence-identity` | [0.0000, 1.0000, 1.0000, 0.0000] | [0.0000, 1.0000, 1.0000, 0.0000] | indices, deg | < 1e-11 | 1e-09 | ✅ pass |

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
| `or-bain-direction-correspondence-identity` | [1.0000, 0.0000, 0.0000, 0.0000] | [1.0000, 0.0000, 0.0000, 0.0000] | indices, deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: The Bain correspondence is constructed from (001)_fcc || (001)_bcc with [110]_fcc || [100]_bcc, so mapping the defining parent direction must recover the defining child direction identically (analytic identity).

**Citation**: Bain, Trans. AIME 70 (1924) 25.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## Kurdjumov-Sachs as a misorientation: 42.85 deg about <0.968 0.178 0.178>

Express the Kurdjumov-Sachs relationship the way it is measured from EBSD boundary data: as the minimal-angle symmetry-reduced misorientation. The published representative is a rotation of 42.85 deg about an axis with components <0.968 0.178 0.178>; the computed angle and sorted absolute axis components are compared against that tabulated value.

**Symbols**

- $(\mathbf{n}, \omega)$ &mdash; Axis-angle pair of the symmetry-reduced misorientation representative.


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
misorientation = ks.misorientation()
axis = np.sort(np.abs(misorientation.rotation.axis))[::-1]
result = np.concatenate([[misorientation.angle_deg], axis])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-ks-misorientation-representation` | [42.8478, 0.9679, 0.1776, 0.1776] | [42.8500, 0.9679, 0.1776, 0.1776] | deg, axis components | 2.24e-03 | 5e-03 | ✅ pass |

**Why this value**: The Kurdjumov-Sachs disorientation representative is tabulated as a 42.85 deg rotation about <0.968 0.178 0.178> in standard thermo-mechanical processing references.

**Citation**: Verlinden, Driver, Samajdar, Doherty, Thermo-Mechanical Processing of Metallic Materials (2007); Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## OR fitting recovers Greninger-Troiano from a Kurdjumov-Sachs start

Twenty parent/child pairs are generated with the Greninger-Troiano relationship, then fitted starting from a Kurdjumov-Sachs nominal. The fit must land on the operative relationship exactly (zero mean residual) while reporting the documented 2.40 deg separation between the fitted relationship and the assumed KS nominal.

**Symbols**

- $(\mathbf{n}, \omega)$ &mdash; Axis-angle pair of the symmetry-reduced misorientation representative.


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

from pytex import (
    Orientation,
    OrientationSet,
    fit_orientation_relationship,
)

specimen = ReferenceFrame(
    name="specimen",
    domain=FrameDomain.SPECIMEN,
    axes=("x", "y", "z"),
    handedness=Handedness.RIGHT,
)
gt = OrientationRelationship.from_greninger_troiano_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
gt_variants = gt.generate_variants()
rng = np.random.default_rng(11)
eulers = rng.uniform(0.0, 60.0, size=(20, 3))
parents = OrientationSet.from_orientations(
    [
        Orientation.from_euler(
            *euler, specimen_frame=specimen, symmetry=austenite.symmetry, phase=austenite
        )
        for euler in eulers
    ]
)
picks = rng.integers(0, len(gt_variants), size=20)
children = OrientationSet(
    quaternions=np.stack(
        [
            parents[index]
            .rotation.compose(gt_variants[int(picks[index])].parent_to_child_rotation.inverse())
            .quaternion
            for index in range(20)
        ],
        axis=0,
    ),
    crystal_frame=ferrite.crystal_frame,
    specimen_frame=specimen,
    symmetry=ferrite.symmetry,
    phase=ferrite,
)
```

:::

**Compute**

```python
ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
report = fit_orientation_relationship(parents, children, ks)
result = np.array([report.deviation_from_nominal_deg, report.mean_residual_deg])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-fit-recovers-gt-from-ks-nominal` | [2.4037, 0.0000] | [2.4037, 0.0000] | deg | < 5e-05 | 5e-03 | ✅ pass |

**Why this value**: Exact GT-generated pairs must refit GT identically (zero residual is an analytic identity), and the reported distance from the KS nominal is the documented KS-GT representative separation of 2.40 deg.

**Citation**: Greninger and Troiano, Trans. AIME 185 (1949) 590; Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## The Sigma3 twin is an admissible Kurdjumov-Sachs same-parent boundary

Deciding whether two neighbouring martensite grains descend from one austenite grain means asking whether their boundary misorientation is one the relationship can actually produce. That admissible set is ``G_c (R G_p R^T) G_c``, because two children of one parent satisfy ``C_i^T C_j = V_i V_j^T``. Two identities are checked: the published Kurdjumov-Sachs intervariant table contains a 60 deg rotation about <111> — the Sigma3 twin relation, Morito's V1-V20 pair — so the exact Sigma3 rotation must sit at zero distance from the fingerprint; and every one of the 276 distinct variant-pair boundaries of a common parent must sit at zero distance too, since they generate the set by construction.

**Symbols**

- $G_c \left(R G_p R^{\mathsf{T}}\right) G_c$ &mdash; Same-parent boundary fingerprint: the admissible child-child misorientations of one parent grain.


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
from pytex import (
    Rotation,
    boundary_fingerprint_distances_deg,
    intervariant_boundary_fingerprint,
)

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
fingerprint = intervariant_boundary_fingerprint(ks)

sigma3 = Rotation.from_axis_angle(
    np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0), np.deg2rad(60.0)
).as_matrix()
sigma3_distance = float(
    boundary_fingerprint_distances_deg(sigma3[None, :, :], fingerprint)[0]
)

variants = ks.generate_variants()
children = np.stack(
    [variant.parent_to_child_rotation.inverse().as_matrix() for variant in variants]
)
left, right = np.triu_indices(len(variants), k=1)
boundaries = np.einsum(
    'nji,njk->nik', children[left], children[right], optimize=True
)
worst_variant_pair = float(
    boundary_fingerprint_distances_deg(boundaries, fingerprint).max()
)
result = [sigma3_distance, worst_variant_pair]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-ks-same-parent-boundary-fingerprint` | [0.0000, 0.0000] | [0.0000, 0.0000] | deg | 1.21e-06 | 1e-05 | ✅ pass |

**Why this value**: Both values are identities, not fitted numbers. The Kurdjumov-Sachs intervariant table published by Morito et al. lists a 60 deg / <111> variant pair (V1-V20), which is exactly the Sigma3 coincidence-site relation, so the Sigma3 rotation belongs to the admissible set. The variant-pair boundaries generate the set by construction, so their distance to it is identically zero. The 1e-5 deg tolerance is the arccos and quaternion/matrix round-trip noise floor, not a physical margin.

**Citation**: Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) 1789 (KS intervariant table); Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## Kurdjumov-Sachs recovered from measured parent/child orientation pairs

The everyday EBSD question: a parent grain and several child grains were indexed, and the operative orientation relationship is wanted. Children are synthesized here through six known Kurdjumov-Sachs variants of one parent, and characterization runs with no nominal relationship supplied, so the answer comes from the data alone. Two quantities are checked: the deviation of the fitted rotation from catalog Kurdjumov-Sachs, and its deviation from Nishiyama-Wassermann. The first must be zero because the data were built from that relationship; the second must be the published separation between the two relationships, which is what makes them distinguishable at all.

**Symbols**

- $\mathbf{R}$ &mdash; Parent-to-child rotation of an orientation relationship.
- $(\mathbf{n}, \omega)$ &mdash; Axis-angle pair of the symmetry-reduced misorientation representative.


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
from pytex import (
    OrientationSet,
    Rotation,
    characterize_orientation_relationship,
    specimen_frame,
)

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
variants = ks.generate_variants()
parent_matrix = Rotation.from_axis_angle([1.0, 2.0, 3.0], 0.7).as_matrix()
# Canonical crystal->specimen convention: C = P V^T.
child_matrices = np.stack(
    [
        parent_matrix @ variants[k].parent_to_child_rotation.as_matrix().T
        for k in (0, 4, 8, 13, 17, 22)
    ]
)
frame = specimen_frame()
parents = OrientationSet.from_matrices(
    np.stack([parent_matrix] * 6), specimen_frame=frame, phase=austenite
)
children = OrientationSet.from_matrices(
    child_matrices, specimen_frame=frame, phase=ferrite
)
report = characterize_orientation_relationship(parents, children)
deviations = dict(zip(report.catalog_names, report.catalog_deviations_deg, strict=True))
result = [deviations["kurdjumov_sachs"], deviations["nishiyama_wassermann"]]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-ks-identified-from-measured-orientations` | [0.0000, 5.2644] | [0.0000, 5.2600] | deg | 4.39e-03 | 1e-02 | ✅ pass |

**Why this value**: The first value is an analytic identity: the children were generated from exact Kurdjumov-Sachs variants, so the fitted rotation must coincide with the relationship it was built from. The second is the tabulated 5.26 deg angular separation between the Kurdjumov-Sachs and Nishiyama-Wassermann relationships. Neither number is copied from a previous program output.

**Citation**: Kurdjumov and Sachs, Z. Phys. 64 (1930) 325; Nishiyama, Sci. Rep. Tohoku Univ. 23 (1934) 637; Wassermann, Arch. Eisenhuettenwes. 16 (1933) 647.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## Reading the Kurdjumov-Sachs parallelisms back out of its rotation

An orientation relationship is stored as a rotation, but the literature reports it as parallel planes and directions. This example recovers that statement from the rotation alone and checks it against the defining Kurdjumov-Sachs parallelisms: the parent plane must belong to {111} and its child image to {011}, the parent direction to <110> and its child image to <111>, all at zero angular deviation. Sorted absolute indices are compared because any member of a family is an equally correct statement of the same relationship.

**Symbols**

- $\mathbf{R}$ &mdash; Parent-to-child rotation of an orientation relationship.
- $(hkl)$ &mdash; Miller plane indices.
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
from pytex import describe_orientation_relationship

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
planes, directions = describe_orientation_relationship(ks)
plane, direction = planes[0], directions[0]
result = np.concatenate(
    [
        np.sort(np.abs(plane.parent_indices)),
        np.sort(np.abs(plane.child_indices)),
        np.sort(np.abs(direction.parent_indices)),
        np.sort(np.abs(direction.child_indices)),
        [plane.deviation_deg, direction.deviation_deg],
    ]
)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-ks-parallelism-statement-from-rotation` | [1.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0000, 0.0000] | [1.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 0.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0000, 0.0000] | indices, deg | < 1e-06 | 1e-04 | ✅ pass |

**Why this value**: The Kurdjumov-Sachs relationship is defined by {111}_fcc || {011}_bcc and <110>_fcc || <111>_bcc, so recovering the statement from the rotation must reproduce exactly those families at zero deviation (analytic identity). The 1e-4 deg tolerance is the matrix-quaternion round-trip noise floor, not a physical margin.

**Citation**: Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## The (111) variant correspondence table is the four Kurdjumov-Sachs packets

Ask what one austenite plane becomes in every martensite variant. Mapping (111) through all 24 Kurdjumov-Sachs variants and grouping the images by index family must reproduce the packet structure of lath martensite: four crystallographically distinct answers, six variants each, of which exactly one group — six variants — carries (111) onto a {011} ferrite plane at zero residual. The computed values are the row count, the number of distinct images, the number of exactly parallel variants, and the smallest and largest group sizes.

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
from pytex import variant_correspondence_table

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
table = variant_correspondence_table(
    ks, CrystalPlane(MillerIndex(np.array([1, 1, 1]), phase=austenite), phase=austenite)
)
exact = table.exact_rows()
group_sizes = {}
for row in table.rows:
    group_sizes[row.equivalence_group] = group_sizes.get(row.equivalence_group, 0) + 1
result = [
    len(table.rows),
    table.distinct_image_count((1, 1, 1)),
    len(exact),
    min(group_sizes.values()),
    max(group_sizes.values()),
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-ks-variant-correspondence-packets` | [24, 4, 6, 6, 6] | [24, 4, 6, 6, 6] | counts | exact | exact | ✅ pass |

**Why this value**: Crystallographic identity, not a measured coincidence. Kurdjumov-Sachs has 24 variants; each carries exactly one member of the four-member {111} family onto its {011} close-packed child plane, so any nominated member is the close-packed plane of exactly 24/4 = 6 of them. Those six are one packet in the sense of Morito et al., and the parent symmetry acts transitively on the remaining images, so every group holds six.

**Citation**: Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) 1789 (packet structure); Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## Every Kurdjumov-Sachs variant carries its own {111}-{011} pair, not variant 1's

A variant is generated as ``V = S_c R S_p^T``, so the parallelism it realizes is the defining pair carried by those operators, ``(S_p n_parent) || (S_c n_child)`` -- not the nominal pair the relationship was written with. ``TransformationVariant.parallel_planes`` returns that per-variant pair. Two numbers check it. First, the worst angle over all 24 variants between the variant rotation applied to its own parent normal and its own child normal must be zero, which is an identity. Second, the number of distinct parent {111} members named across the variants must be 4 -- the four close-packed planes of the fcc parent, one per Morito packet. Substituting the nominal pair instead would name a single member and open a non-zero angle, which is exactly the figure that looks right and is wrong.

**Symbols**

- $\mathbf{S}_{p}$ &mdash; Parent point-group operator generating a transformation variant.
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
worst_deviation_deg = 0.0
parent_members = set()
for variant in ks.generate_variants():
    rotation = variant.parent_to_child_rotation.as_matrix()
    for parent_plane, child_plane in variant.parallel_planes:
        rotated = rotation @ parent_plane.normal
        # atan2 of the chord against the sum, not arccos of the dot product:
        # these normals are meant to be exactly parallel, and arccos cannot
        # report an angle below about 1e-6 deg however exact the pair is.
        deviation = np.rad2deg(
            2.0
            * np.arctan2(
                float(np.linalg.norm(rotated - child_plane.normal)),
                float(np.linalg.norm(rotated + child_plane.normal)),
            )
        )
        worst_deviation_deg = max(worst_deviation_deg, float(deviation))
    indices = np.asarray(variant.parallel_planes[0][0].miller.indices)
    canonical = indices if indices[0] >= 0 else -indices
    parent_members.add(tuple(int(value) for value in canonical))
result = [worst_deviation_deg, float(len(parent_members))]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-ks-variant-parallelisms-are-per-variant` | [0.0000, 4.0000] | [0.0000, 4.0000] | deg, count | < 1e-08 | 1e-06 | ✅ pass |

**Why this value**: Both values are identities. ``V = S_c R S_p^T`` maps ``S_p n_parent`` onto ``S_c n_child`` by construction, so the deviation is exactly zero up to floating-point noise. The fcc {111} family has four members and the 24 Kurdjumov-Sachs variants distribute over them six apiece, which is the packet structure of lath martensite reported by Morito et al.

**Citation**: Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) 1789 (packet structure); Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## The OR dossier reports the numbers its own functions report

The dossier aggregates an orientation-relationship declaration into one serializable document. Its whole value rests on a rule -- that it calls the existing functions and never reimplements them -- so what is worth checking is not any single number but the agreement. Six values are computed here: the difference between the cell volume the dossier reports and the one the lattice reports; the difference between that volume and the cube of the cubic edge; the variant count; the packet count; the number of distinct intervariant angles; and the largest of them. The first two are identities and must be exactly zero; the last four are the published Kurdjumov-Sachs figures -- 24 variants, 4 packets, 10 distinct angles, the largest being the 60 degree Sigma3 twin relation.

**Symbols**

- $V$ &mdash; Unit-cell volume of a phase.
- $(\mathbf{n}, \omega)$ &mdash; Axis-angle pair of the symmetry-reduced misorientation representative.


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
from pytex import or_dossier

ks = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
dossier = or_dossier(ks)
payload = dossier.to_json()

# Each number beside the function a reader would check it against.
volume = payload['parent']['volume_angstrom3']
spectrum = payload['misorientation']
result = [
    volume - austenite.lattice.volume_angstrom3(),
    volume - 3.6 ** 3,
    float(spectrum['variant_count']),
    float(spectrum['packet_count']),
    float(len(spectrum['intervariant_angles_deg'])),
    max(spectrum['intervariant_angles_deg']),
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-dossier-agrees-with-its-sources` | [0.0000, 0.0000, 24.0000, 4.0000, 10.0000, 60.0000] | [0.0000, 0.0000, 24.0000, 4.0000, 10.0000, 60.0000] | angstrom^3, angstrom^3, counts, deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: The first two entries are identities: the dossier reads the volume from the lattice rather than recomputing it, and the volume of a cubic cell is the cube of its edge. The remaining four are the published Kurdjumov-Sachs figures -- 24 crystallographically distinct variants, four packets on the four members of the parent {111} family, and the ten distinct intervariant disorientations of Morito et al., whose largest is the 60 degree rotation about <111> -- the Sigma3 twin relation.

**Citation**: Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) 1789 (intervariant table and packet structure); Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`

## Writing Greninger-Troiano in low indices costs the 2.40 deg to Kurdjumov-Sachs

Six exact Greninger-Troiano children of one austenite grain are characterized, then the fitted relationship is restated in integers with the index bound held at two. Greninger-Troiano has no low-index direction pair, so the tidiest statement available at that bound is the Kurdjumov-Sachs one -- {111} parallel to {110} with <110> parallel to <111> -- and writing it is not free. Four numbers are computed: the fit residual, which is zero because the data are exact; whether the rationalized plane family is {111} and the direction family <110>; and the price of the idealization. That price must be the published Kurdjumov-Sachs to Greninger-Troiano separation. An idealization returned without it would read as a measurement of Kurdjumov-Sachs.

**Symbols**

- $\Delta\omega$ &mdash; Symmetry-reduced angle between a measured relationship and the integer statement it is idealized to.
- $(hkl)$ &mdash; Miller plane indices.
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
from pytex import (
    OrientationSet,
    Rotation,
    characterize_orientation_relationship,
    specimen_frame,
)

gt = OrientationRelationship.from_greninger_troiano_correspondence(
    parent_phase=austenite, child_phase=ferrite
)
variants = gt.generate_variants()
parent_matrix = Rotation.from_axis_angle([1.0, 2.0, 3.0], 0.7).as_matrix()
picks = (0, 4, 8, 13, 17, 22)
# Canonical crystal->specimen convention: C = P V^T.
children = np.stack(
    [parent_matrix @ variants[k].parent_to_child_rotation.as_matrix().T for k in picks]
)
frame = specimen_frame()
report = characterize_orientation_relationship(
    OrientationSet.from_matrices(
        np.stack([parent_matrix] * len(picks)), specimen_frame=frame, phase=austenite
    ),
    OrientationSet.from_matrices(children, specimen_frame=frame, phase=ferrite),
)

rationalized = report.as_rational_relationship(max_index=2)
plane = sorted(abs(int(v)) for v in rationalized.plane_statement.parent_indices)
direction = sorted(abs(int(v)) for v in rationalized.direction_statement.parent_indices)
result = [
    report.mean_residual_deg,
    float(plane == [1, 1, 1]),
    float(direction == [0, 1, 1]),
    rationalized.residual_rotation_deg,
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-rationalization-costs-the-ks-gt-separation` | [0.0000, 1.0000, 1.0000, 2.4037] | [0.0000, 1.0000, 1.0000, 2.4037] | deg, booleans, deg | < 5e-05 | 5e-03 | ✅ pass |

**Why this value**: The zero residual is an analytic identity: the pairs were generated from the relationship being fitted. The 2.40 deg is the documented separation between the Greninger-Troiano and Kurdjumov-Sachs representatives, which is exactly what it costs to write the former with the latter's indices.

**Citation**: Greninger and Troiano, Trans. AIME 185 (1949) 590; Kurdjumov and Sachs, Z. Phys. 64 (1930) 325.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`Transformation API <../../api/index>`
