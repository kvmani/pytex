<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Crystal geometry: angles, spacings, and multiplicities

Interplanar and interdirection angles, interplanar spacings, and symmetry multiplicities for cubic and hexagonal phases, and the naming of a direction or a plane that arrived as a Cartesian vector. Each result is checked against an analytic identity for the relevant crystal system.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Angle between (100) and (110) in a cubic crystal

You have indexed two poles as {100} and {110} on a cubic phase and want to confirm the geometry of a pole figure or a Kikuchi band intersection. In a cubic system the answer is exactly 45 degrees, independent of the lattice parameter, so this is the first check that your frame and symmetry wiring is correct.

**Symbols**

- $\angle(\mathbf{n}_1, \mathbf{n}_2)$ &mdash; Angle between two plane normals.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
result = float(np.degrees(angle_plane_plane_rad(
    MillerPlane.from_hkl([1, 0, 0], phase=cubic),
    MillerPlane.from_hkl([1, 1, 0], phase=cubic),
)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `cubic-angle-100-110` | 45.0000 | 45.0000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: For cubic metrics the interplanar angle is arccos(h1 h2 + k1 k2 + l1 l2 over norms); for (100) and (110) this is arccos(1/sqrt(2)) = 45 degrees, independent of a.

**Citation**: Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed., Appendix 3.

**See also**: {doc}`Miller planes and directions <../../concepts/miller_planes_directions>`, {doc}`angle_plane_plane_rad / angle_dir_dir_rad <../../api/index>`

## Angle between (100) and (111) in a cubic crystal

The 100-to-111 angle sets the classic stereographic-triangle geometry and appears whenever you relate a rolling-plane normal to an octahedral slip plane. The exact value is arccos(1/sqrt(3)).

**Symbols**

- $\angle(\mathbf{n}_1, \mathbf{n}_2)$ &mdash; Angle between two plane normals.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
result = float(np.degrees(angle_plane_plane_rad(
    MillerPlane.from_hkl([1, 0, 0], phase=cubic),
    MillerPlane.from_hkl([1, 1, 1], phase=cubic),
)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `cubic-angle-100-111` | 54.7356 | 54.7356 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: arccos(1/sqrt(3)) = 54.7356 degrees for cubic (100) vs (111).

**Citation**: Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed., Appendix 3.

**See also**: {doc}`Miller planes and directions <../../concepts/miller_planes_directions>`, {doc}`angle_plane_plane_rad / angle_dir_dir_rad <../../api/index>`

## Angle between [110] and [111] directions in a cubic crystal

Slip-system and Schmid-factor calculations repeatedly need the angle between a slip direction such as [110] and a loading or plane-normal direction such as [111]. In cubic metrics the direction angle equals the same-index plane angle.

**Symbols**

- $\angle(\mathbf{d}_1, \mathbf{d}_2)$ &mdash; Angle between two lattice directions.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
result = float(np.degrees(angle_dir_dir_rad(
    MillerDirection.from_uvw([1, 1, 0], phase=cubic),
    MillerDirection.from_uvw([1, 1, 1], phase=cubic),
)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `cubic-angle-dir-110-111` | 35.2644 | 35.2644 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: arccos(sqrt(2/3)) = 35.2644 degrees for cubic [110] vs [111].

**Citation**: Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed., Appendix 3.

**See also**: {doc}`Miller planes and directions <../../concepts/miller_planes_directions>`, {doc}`angle_plane_plane_rad / angle_dir_dir_rad <../../api/index>`

## Interplanar spacing of (111) in a cubic crystal (a = 4 angstrom)

Interplanar spacing is the bridge from crystallography to diffraction: it fixes the Bragg angle for a reflection. For a cubic lattice d_hkl = a / sqrt(h^2 + k^2 + l^2), so d_111 = a / sqrt(3).

**Symbols**

- $d_{hkl}$ &mdash; Interplanar spacing of the (hkl) family.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
result = MillerPlane.from_hkl([1, 1, 1], phase=cubic).d_spacing_angstrom
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `cubic-dspacing-111` | 2.30940 | 2.30940 | angstrom | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: d_111 = a / sqrt(3) = 4 / sqrt(3) = 2.30940 angstrom.

**Citation**: Cullity and Stock, Elements of X-Ray Diffraction, 3rd ed., Eq. 3-10.

**See also**: {doc}`Miller planes and directions <../../concepts/miller_planes_directions>`, {doc}`Diffraction geometry worked examples <diffraction>`

## Symmetry multiplicity of {100}, {110}, {111}, {321} under m-3m

Powder-diffraction intensities and pole-figure normalization both depend on how many symmetry-equivalent planes a family contains. PyTex treats plane families with antipodal equivalence (a plane and its opposite normal are the same plane), so the reduced multiplicities are half the full point-group orbit: 3, 6, 4, and 24.

**Symbols**

- $m_{\{hkl\}}$ &mdash; Symmetry multiplicity of a plane family.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
result = [
    MillerPlane.from_hkl(hkl, phase=cubic).symmetry_equivalent_indices()[0].shape[1]
    for hkl in ([1, 0, 0], [1, 1, 0], [1, 1, 1], [3, 2, 1])
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `cubic-plane-multiplicity` | [3, 6, 4, 24] | [3, 6, 4, 24] | &mdash; | exact | exact | ✅ pass |

**Why this value**: Full m-3m orbit sizes are 6, 12, 8, 48; antipodal folding halves them to 3, 6, 4, 24.

**Citation**: Hahn (ed.), International Tables for Crystallography Vol. A, point group m-3m.

**See also**: {doc}`Symmetry and fundamental regions <../../concepts/symmetry_and_fundamental_regions>`, {doc}`angle_plane_plane_rad / angle_dir_dir_rad <../../api/index>`

## Angle between (0001) and (10-10) in a hexagonal crystal

In HCP metals the basal plane (0001) and the prismatic planes {10-10} are the dominant slip and texture planes. Their normals are exactly perpendicular for any c/a ratio, which makes this a robust convention check for the hexagonal metric and the four-index handling.

**Symbols**

- $\angle(\mathbf{n}_1, \mathbf{n}_2)$ &mdash; Angle between two plane normals.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
hexagonal = Phase(
    "hcp-demo",
    lattice=Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
result = float(np.degrees(angle_plane_plane_rad(
    MillerPlane.from_hkl([0, 0, 1], phase=hexagonal),
    MillerPlane.from_hkl([1, 0, 0], phase=hexagonal),
)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hex-angle-basal-prism` | 90.0000 | 90.0000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: The basal-plane normal is c*; prismatic normals lie in the basal plane, so the angle is 90 degrees.

**Citation**: Partridge, The crystallography and deformation modes of HCP metals, Metall. Rev. 12 (1967).

**See also**: {doc}`Hexagonal and trigonal conventions <../../standards/hexagonal_and_trigonal_conventions>`, {doc}`angle_plane_plane_rad / angle_dir_dir_rad <../../api/index>`

## Angle between (10-10) and (01-10) in a hexagonal crystal

Adjacent first-order prismatic planes in HCP are separated by 60 degrees. This confirms that the 120-degree gamma angle of the hexagonal cell is handled correctly when reasoning about prismatic slip variants.

**Symbols**

- $\angle(\mathbf{n}_1, \mathbf{n}_2)$ &mdash; Angle between two plane normals.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
hexagonal = Phase(
    "hcp-demo",
    lattice=Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
result = float(np.degrees(angle_plane_plane_rad(
    MillerPlane.from_hkl([1, 0, 0], phase=hexagonal),
    MillerPlane.from_hkl([0, 1, 0], phase=hexagonal),
)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hex-angle-prism-prism` | 60.0000 | 60.0000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: First-order prismatic normals are separated by 60 degrees in the hexagonal basal plane.

**Citation**: Partridge, The crystallography and deformation modes of HCP metals, Metall. Rev. 12 (1967).

**See also**: {doc}`Hexagonal and trigonal conventions <../../standards/hexagonal_and_trigonal_conventions>`, {doc}`angle_plane_plane_rad / angle_dir_dir_rad <../../api/index>`

## Naming a direction that arrived as geometry recovers its own indices

An inverse pole figure, a stereogram readout and an interactive viewer all face the same problem: the direction is known as a vector and the reader wants a label. The naming is a search over low-index triples, so the first thing to establish is that it is exact where it should be — a direction built from [3 2 1] must be named [3 2 1], with an angular residual of zero rather than of the search's step size.

**Symbols**

- $\angle(\mathbf{d}_1, \mathbf{d}_2)$ &mdash; Angle between two lattice directions.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
cubic = Phase(
    "cubic-demo",
    lattice=Lattice(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("m-3m", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
vector = MillerDirection.from_uvw([3, 2, 1], phase=cubic).unit_vector_cartesian
indices, residual_deg = nearest_low_index_direction(vector, phase=cubic)
result = float(residual_deg)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `cubic-nearest-low-index-round-trip` | 0.0000 | 0.0000 | deg | < 1e-11 | 1e-09 | ✅ pass |

**Why this value**: The vector is the Cartesian image of [3 2 1] under the direct basis, and gcd(3, 2, 1) = 1 so the triple is already primitive. The nearest primitive triple to a direction that is one is therefore itself, at zero angle.

**Citation**: International Tables for Crystallography, Volume A, section 1.3 (direct and reciprocal bases and their index conventions).

**See also**: {doc}`Miller planes and directions <../../concepts/miller_planes_directions>`, {doc}`angle_plane_plane_rad / angle_dir_dir_rad <../../api/index>`

## In a hexagonal phase a plane normal and the like-indexed direction are not parallel

Miller indices are reciprocal-basis components and direction indices are direct-basis components, so naming a vector as a plane and naming it as a direction are different questions with different answers in every non-cubic phase. The pyramidal normal is the standard demonstration: name the (10-11) normal as a plane and it comes back exactly, name it as a direction and the angle to [10-11] is the metric difference itself.

**Symbols**

- $\angle(\mathbf{n}_1, \mathbf{n}_2)$ &mdash; Angle between two plane normals.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameDomain,
    Handedness,
    Lattice,
    MillerDirection,
    MillerPlane,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
    angle_dir_dir_rad,
    angle_plane_plane_rad,
    nearest_low_index_direction,
    nearest_low_index_plane,
)

crystal = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)
hexagonal = Phase(
    "hcp-demo",
    lattice=Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=crystal),
    symmetry=SymmetrySpec.from_point_group("6/mmm", reference_frame=crystal),
    crystal_frame=crystal,
)
```

:::

**Compute**

```python
normal = MillerPlane.from_hkl([1, 0, 1], phase=hexagonal).normal_cartesian
indices, residual_deg = nearest_low_index_plane(normal, phase=hexagonal)
result = [int(value) for value in indices]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `hex-nearest-low-index-plane-vs-direction` | [1, 0, 1] | [1, 0, 1] | &mdash; | exact | exact | ✅ pass |

**Why this value**: The vector is the unit normal of (10-11) by construction, and the plane search runs against the reciprocal basis, so it must return the indices the normal was built from. The direction search, which runs against the direct basis, returns a different triple for the same vector because c/a is not 1.

**Citation**: International Tables for Crystallography, Volume A, section 1.3; Partridge, The crystallography and deformation modes of HCP metals, Metall. Rev. 12 (1967).

**See also**: {doc}`Hexagonal and trigonal conventions <../../standards/hexagonal_and_trigonal_conventions>`, {doc}`angle_plane_plane_rad / angle_dir_dir_rad <../../api/index>`
