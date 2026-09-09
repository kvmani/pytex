<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Weighted orientation-relationship fitting

Retain every measured pair while controlling its influence with declared evidence weights. These checks use a common-axis circular identity and an exact exclusion.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## A weighted quaternion fit agrees with the circular mean

Two pairs differ by two degrees about one axis. Give the first three times the evidence weight and compare the fitted displacement against the analytic circular mean.

**Symbols**

- $w_i$ &mdash; Normalized scalar evidence weight of measured pair i.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import OrientationRelationship, OrientationSet, Rotation, specimen_frame
from pytex.core import fit_orientation_relationship
from pytex.app.phases import builtin_phase

parent = builtin_phase('austenite_fcc').to_phase()
child = builtin_phase('fe_bcc').to_phase()
nominal = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=parent, child_phase=child,
)
angles = np.radians([0.0, 2.0])
measured = np.stack([
    Rotation.from_axis_angle([1, 0, 0], angle).as_matrix()
    @ nominal.parent_to_child_rotation.as_matrix() for angle in angles
])
parents = OrientationSet.from_matrices(
    np.repeat(np.eye(3)[None], 2, axis=0), specimen_frame=specimen_frame(), phase=parent,
)
children = OrientationSet.from_matrices(
    measured.transpose(0, 2, 1), specimen_frame=specimen_frame(), phase=child,
)
```

:::

**Compute**

```python
report = fit_orientation_relationship(parents, children, nominal, pair_weights=[3, 1])
analytic_angle = np.degrees(np.arctan2(np.sin(angles[1]), 3 + np.cos(angles[1])))
result = abs(report.residuals_deg[0] - analytic_angle)
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-weighted-fit-circular-identity` | 0.0000 | 0.0000 | deg | < 1e-12 | 1e-10 | ✅ pass |

**Why this value**: For common-axis rotations the weighted chordal optimum has angle arg(sum(w*exp(i*angle))); the two independent formulas must agree.

**Citation**: Markley et al. (2007), Averaging Quaternions, doi:10.2514/1.28949.

**See also**: {doc}`OR determination <../../theory/orientation_relationship_determination>`

## Excluding a pair leaves its residual visible

Exclude the perturbed second pair while retaining it for review. The first exact pair fixes the rotation; the second must still report its imposed two-degree residual.

**Symbols**

- $w_i$ &mdash; Normalized scalar evidence weight of measured pair i.


:::{dropdown} Setup (imports and object construction)

```python
import numpy as np
from pytex import OrientationRelationship, OrientationSet, Rotation, specimen_frame
from pytex.core import fit_orientation_relationship
from pytex.app.phases import builtin_phase

parent = builtin_phase('austenite_fcc').to_phase()
child = builtin_phase('fe_bcc').to_phase()
nominal = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=parent, child_phase=child,
)
angles = np.radians([0.0, 2.0])
measured = np.stack([
    Rotation.from_axis_angle([1, 0, 0], angle).as_matrix()
    @ nominal.parent_to_child_rotation.as_matrix() for angle in angles
])
parents = OrientationSet.from_matrices(
    np.repeat(np.eye(3)[None], 2, axis=0), specimen_frame=specimen_frame(), phase=parent,
)
children = OrientationSet.from_matrices(
    measured.transpose(0, 2, 1), specimen_frame=specimen_frame(), phase=child,
)
```

:::

**Compute**

```python
report = fit_orientation_relationship(parents, children, nominal, pair_weights=[1, 0])
assert report.weighted_mean_residual_deg < 1e-10
assert report.effective_pair_count == 1.0
result = report.residuals_deg[1]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `or-weighted-fit-excluded-residual` | 2.0000 | 2.0000 | deg | < 1e-12 | 1e-10 | ✅ pass |

**Why this value**: The excluded pair was constructed by a two-degree rotation relative to the only positive-weight pair. It contributes no term to the scatter matrix.

**Citation**: Markley et al. (2007), Averaging Quaternions, doi:10.2514/1.28949.

**See also**: {doc}`OR determination <../../theory/orientation_relationship_determination>`
