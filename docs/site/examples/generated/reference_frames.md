<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Reference Frames And Frame Transforms

Creating standard frames, declaring frame relationships in words, and letting the frame graph compose multi-step chains — with the rotation angles, components, and invariants checked against exact analytic values. The last two examples pin the IUCr notation convention: the reciprocal star marks the basis, never the indices.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Rotation angle implied by a declared axis correspondence

EBSD vendors and analysis tools disagree about which specimen axis is called what. Rather than hand-writing a permutation matrix, you declare the correspondence in words: specimen x is the sample TD axis, specimen y is the reversed RD axis, specimen z is ND. That declaration is a 90-degree rotation about the shared third axis, and this example checks that PyTex builds exactly that rotation.

**Symbols**

- $\mathbf{R}$ &mdash; Rotation matrix mapping source-frame to target-frame components.
- $\omega$ &mdash; Rotation angle of a frame-to-frame transform.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameGraph,
    FrameTransform,
    crystal_frame,
    rolling_frame_graph,
    sample_frame,
    specimen_frame,
)

specimen = specimen_frame()
sample = sample_frame()
crystal = crystal_frame()
```

:::

**Compute**

```python
transform = FrameTransform.from_axis_correspondence(
    specimen, sample, {"x": "TD", "y": "-RD", "z": "ND"}
)
result = transform.rotation_angle_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `frame-axis-correspondence-angle` | 90.0000 | 90.0000 | deg | < 1e-12 | 1e-09 | ✅ pass |

**Why this value**: The declaration fixes R e_x = e_TD, R e_y = -e_RD, R e_z = e_ND, i.e. the signed permutation [[0,-1,0],[1,0,0],[0,0,1]]. Its trace is 1, and the rotation angle follows from cos(omega) = (trace - 1) / 2 = 0, so omega = 90 degrees exactly.

**Citation**: Hahn, Th. (ed.), International Tables for Crystallography, Volume A, DOI: 10.1107/97809553602060000100 (right-handed axis conventions); trace identity for a proper rotation, Bunge, Texture Analysis in Materials Science, DOI: 10.1016/C2013-0-11769-2.

**See also**: {doc}`Reference frames and conventions <../../concepts/reference_frames_and_conventions>`, {doc}`ReferenceFrame / FrameTransform / FrameGraph <../../api/index>`

## Components of the specimen x axis in a relabelled sample frame

The point of a typed frame transform is that it converts components, not just angles. Having declared that specimen x is the sample TD axis, a direction lying along specimen x must come back with sample components (0, 1, 0): purely TD, no RD, no ND. This is the check that catches a reversed or transposed convention immediately.

**Symbols**

- $\mathbf{R}$ &mdash; Rotation matrix mapping source-frame to target-frame components.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameGraph,
    FrameTransform,
    crystal_frame,
    rolling_frame_graph,
    sample_frame,
    specimen_frame,
)

specimen = specimen_frame()
sample = sample_frame()
crystal = crystal_frame()
```

:::

**Compute**

```python
transform = FrameTransform.from_axis_correspondence(
    specimen, sample, {"x": "TD", "y": "-RD", "z": "ND"}
)
result = np.asarray(transform.apply_to_directions(np.array([1.0, 0.0, 0.0])))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `frame-axis-correspondence-components` | [0.0000, 1.0000, 0.0000] | [0.0000, 1.0000, 0.0000] | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: By definition of a basis, the source frame's x axis has source components e_x. The declaration 'x is TD' therefore forces its target components to be the TD basis vector (0, 1, 0).

**Citation**: Hahn, Th. (ed.), International Tables for Crystallography, Volume A, DOI: 10.1107/97809553602060000100.

**See also**: {doc}`Reference frames and conventions <../../concepts/reference_frames_and_conventions>`, {doc}`ReferenceFrame / FrameTransform / FrameGraph <../../api/index>`

## Composing a two-hop frame chain with the frame graph

A rolled sheet mounted 30 degrees off the stage axis gives two declared relationships: the canonical Cartesian reference to the specimen frame, and the specimen frame to the RD/TD/ND sample frame. You never declared the direct Cartesian-to-sample relationship, but you need it. The frame graph composes the shortest declared chain for you, and the result must be the 30-degree mounting rotation.

**Symbols**

- $\omega$ &mdash; Rotation angle of a frame-to-frame transform.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameGraph,
    FrameTransform,
    crystal_frame,
    rolling_frame_graph,
    sample_frame,
    specimen_frame,
)

specimen = specimen_frame()
sample = sample_frame()
crystal = crystal_frame()
```

:::

**Compute**

```python
graph = rolling_frame_graph(rd_offset_deg=30.0)
transform = graph.transform_between("cartesian", "sample_rd_td_nd")
result = transform.rotation_angle_deg
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `frame-graph-multi-hop-angle` | 30.0000 | 30.0000 | deg | < 1e-12 | 1e-09 | ✅ pass |

**Why this value**: The Cartesian-to-specimen edge is the identity and the specimen-to-sample edge is the declared 30-degree mounting rotation about ND, so the composition R = R2 R1 = R2 has rotation angle 30 degrees exactly.

**Citation**: Bunge, H.-J., Texture Analysis in Materials Science: Mathematical Methods, DOI: 10.1016/C2013-0-11769-2 (specimen-frame conventions for rolled sheet).

**See also**: {doc}`Reference frames and conventions <../../concepts/reference_frames_and_conventions>`, {doc}`ReferenceFrame / FrameTransform / FrameGraph <../../api/index>`

## Round-tripping components through a frame chain and back

Any chain of frame transforms must be exactly invertible: converting a direction from the specimen frame into the sample frame and back has to return the original components. This is the invariant that guarantees no convention is silently lost when data crosses several module boundaries, so the residual is checked against exact zero.

**Symbols**

- $\mathbf{R}$ &mdash; Rotation matrix mapping source-frame to target-frame components.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameGraph,
    FrameTransform,
    crystal_frame,
    rolling_frame_graph,
    sample_frame,
    specimen_frame,
)

specimen = specimen_frame()
sample = sample_frame()
crystal = crystal_frame()
```

:::

**Compute**

```python
graph = rolling_frame_graph(rd_offset_deg=37.5)
direction = np.array([0.3, -0.7, 0.5])
forward = graph.convert(
    direction, source="specimen", target="sample_rd_td_nd", directions=True
)
back = graph.convert(
    forward, source="sample_rd_td_nd", target="specimen", directions=True
)
result = float(np.max(np.abs(np.asarray(back) - direction)))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `frame-round-trip-residual` | 0.0000 | 0.0000 | &mdash; | < 1e-12 | 1e-14 | ✅ pass |

**Why this value**: A frame transform is a proper rotation, so R^-1 = R^T and R^T R = I exactly. The round-trip residual is therefore zero up to floating-point rounding.

**Citation**: Orthogonality of proper rotations; Bunge, Texture Analysis in Materials Science, DOI: 10.1016/C2013-0-11769-2.

**See also**: {doc}`Reference frames and conventions <../../concepts/reference_frames_and_conventions>`, {doc}`ReferenceFrame / FrameTransform / FrameGraph <../../api/index>`

## A right-handed frame has axis-vector determinant +1

PyTex refuses to build a frame whose declared handedness contradicts its axis geometry, because a silently mirrored frame turns every downstream chirality result — variant selection, twin sense, pole-figure handedness — inside out. This example shows the invariant being reported for the standard RD/TD/ND sample frame.

**Symbols**

- $\mathbf{R}$ &mdash; Rotation matrix mapping source-frame to target-frame components.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    FrameGraph,
    FrameTransform,
    crystal_frame,
    rolling_frame_graph,
    sample_frame,
    specimen_frame,
)

specimen = specimen_frame()
sample = sample_frame()
crystal = crystal_frame()
```

:::

**Compute**

```python
result = sample.determinant
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `frame-right-handed-determinant` | 1.0000 | 1.0000 | &mdash; | < 1e-12 | 1e-12 | ✅ pass |

**Why this value**: For a right-handed orthonormal triad the determinant of the matrix whose columns are the axis vectors is exactly +1; a left-handed triad gives -1.

**Citation**: Hahn, Th. (ed.), International Tables for Crystallography, Volume A, DOI: 10.1107/97809553602060000100 (right-handed axial-frame convention).

**See also**: {doc}`Reference frames and conventions <../../concepts/reference_frames_and_conventions>`, {doc}`ReferenceFrame / FrameTransform / FrameGraph <../../api/index>`

## Every reciprocal-frame axis carries the IUCr star

In a workflow holding both direct and reciprocal quantities, the single most valuable safeguard is that a reciprocal-space vector cannot be mistaken for a direct-space one. PyTex enforces that by starring every axis of a reciprocal-domain frame. This example counts the starred axes on both frames: exactly three on the reciprocal frame, because the star belongs to the basis, and none on the direct crystal frame.

**Symbols**

- $\mathbf{a}^{*}$ &mdash; Reciprocal basis vector, dual to the direct basis vector a.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    crystal_frame,
    format_plane_family_indices,
    format_plane_indices,
    reciprocal_frame_for,
)

crystal = crystal_frame()
reciprocal = reciprocal_frame_for(crystal)
```

:::

**Compute**

```python
starred = sum(1 for axis in reciprocal.axes if axis.endswith("*"))
direct_starred = sum(1 for axis in crystal.axes if axis.endswith("*"))
result = np.array([starred, direct_starred])
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `reciprocal-frame-star-count` | [3, 0] | [3, 0] | &mdash; | exact | exact | ✅ pass |

**Why this value**: The reciprocal basis of a three-dimensional lattice has exactly three vectors a*, b*, c*, each conventionally starred; the direct basis a, b, c is never starred.

**Citation**: Hahn, Th. (ed.), International Tables for Crystallography, Volume A, DOI: 10.1107/97809553602060000100; reciprocal-space definitions, International Tables Volume C.

**See also**: {doc}`Reference frames and conventions <../../concepts/reference_frames_and_conventions>`, {doc}`ReferenceFrame / FrameTransform / FrameGraph <../../api/index>`

## Miller indices are never starred, in any bracket form

The star marks the basis, not the indices. Miller indices are already reciprocal-basis components by definition, so starring them would name a different quantity - a mistake that is easy to make when 'reciprocal quantities are starred' is applied too broadly. This example counts the stars produced by every bracket form, which must be zero.

**Symbols**

- $\mathbf{g}_{hkl}$ &mdash; Reciprocal-lattice vector of the (hkl) reflection.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import (
    crystal_frame,
    format_plane_family_indices,
    format_plane_indices,
    reciprocal_frame_for,
)

crystal = crystal_frame()
reciprocal = reciprocal_frame_for(crystal)
```

:::

**Compute**

```python
forms = [
    format_plane_indices((1, 1, 1), style="plain"),
    format_plane_family_indices((1, 1, 1), style="plain"),
    format_plane_indices((1, -1, 0), style="plain"),
]
result = float(sum(text.count("*") for text in forms))
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `miller-indices-carry-no-star` | 0.0000 | 0.0000 | &mdash; | exact | exact | ✅ pass |

**Why this value**: By definition g_hkl = h a* + k b* + l c*: the indices (h, k, l) are the scalar coefficients of the starred basis vectors, so the indices carry no star themselves.

**Citation**: Hahn, Th. (ed.), International Tables for Crystallography, Volume A, DOI: 10.1107/97809553602060000100.

**See also**: {doc}`Reference frames and conventions <../../concepts/reference_frames_and_conventions>`, {doc}`ReferenceFrame / FrameTransform / FrameGraph <../../api/index>`
