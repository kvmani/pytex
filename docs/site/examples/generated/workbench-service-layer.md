<!-- GENERATED FILE. Do not edit by hand.
     Source of truth: worked_examples/ (rendered by scripts/generate_worked_examples.py).
     Run `python scripts/generate_worked_examples.py` to regenerate. -->

# Workbench service layer

The three quantitative claims the workbench user guide makes, each checked against a value fixed independently of this code: the Kurdjumov-Sachs packet structure and intervariant spectrum from Morito et al., the closure of the m.r.d. scale as an exact identity, and the assertion a Miller component label makes about where its poles land.

```{note}
Every number on this page is computed live from the public PyTex API when the documentation is regenerated, then checked against an independently known reference value by `tests/unit/test_worked_examples.py`. The code shown is exactly the code that produced the computed value, so you can copy any snippet and reproduce the tabulated output.
```

## Kurdjumov-Sachs gives four packets of six variants

One austenite grain transforming under Kurdjumov-Sachs produces 24 child orientations, and they are not unstructured: each variant carries exactly one member of the parent {111} family into exact parallelism with a child {110}, and the variants sharing that member form a packet. The {111} family has four members, so there are four packets, and 24 variants divided among them evenly gives six each. That grouping is what a lath martensite micrograph shows as a block, and it is why one parent grain gives 24 orientations but only four apparent plate directions.

:::{dropdown} Setup (imports and object construction)

```python
from pytex.app import REGISTRY

AUSTENITE_TO_FERRITE = {
    "phase": {"builtin": "austenite_fcc"},
    "child_phase": {"builtin": "fe_bcc"},
    "relationship": "kurdjumov_sachs",
    "packet_plane": [1, 1, 1],
}
```

:::

**Compute**

```python
response = REGISTRY.call(
    'variants.pole_figure',
    dict(AUSTENITE_TO_FERRITE, pole=[1, 0, 0], projection='stereographic',
         include_parent=False),
)
sizes = response['data']['packet_sizes']
result = sorted(sizes.values())
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `workbench-ks-packet-size` | [6, 6, 6, 6] | [6, 6, 6, 6] | variants per packet | exact | exact | ✅ pass |

**Why this value**: Morito et al. report the 24 Kurdjumov-Sachs variants of lath martensite as four packets of six, one packet per member of the parent {111} family.

**Citation**: Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) 1789, Table 2.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`The PyTex Workbench <../../workflows/workbench_application>`

## The 276 variant pairs fall on ten disorientations

Two child grains that grew from the same parent cannot meet at an arbitrary misorientation: the admissible set is fixed by the relationship and the two point groups, and it is discrete. That discreteness is what makes a measured misorientation histogram a test — peaks away from these ten angles are boundaries between different parent grains, which is the reasoning parent-grain reconstruction rests on. The 24 variants give 276 unordered pairs, and they collapse onto ten values.

:::{dropdown} Setup (imports and object construction)

```python
from pytex.app import REGISTRY

AUSTENITE_TO_FERRITE = {
    "phase": {"builtin": "austenite_fcc"},
    "child_phase": {"builtin": "fe_bcc"},
    "relationship": "kurdjumov_sachs",
    "packet_plane": [1, 1, 1],
}
```

:::

**Compute**

```python
response = REGISTRY.call(
    'variants.intervariant_misorientations',
    dict(AUSTENITE_TO_FERRITE, merge_equal_angles=True),
)
result = [row['angle_deg'] for row in response['table']['rows']]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `workbench-ks-intervariant-spectrum` | [10.53, 14.88, 20.61, 21.06, 47.11, 49.47, 50.51, 51.73, 57.21, 60.00] | [10.53, 14.88, 20.61, 21.06, 47.11, 49.47, 50.51, 51.73, 57.21, 60.00] | deg | 5.00e-03 | 1e-02 | ✅ pass |

**Why this value**: The ten distinct intervariant disorientation angles of the Kurdjumov-Sachs relationship, as tabulated for lath martensite by Morito et al.

**Citation**: Morito, Tanaka, Konishi, Furuhara and Maki, Acta Materialia 51 (2003) 1789, Table 2.

**See also**: {doc}`Orientation relationships <../../concepts/orientation_relationships>`, {doc}`The PyTex Workbench <../../workflows/workbench_application>`

## The area-weighted mean of any pole figure is 1 m.r.d.

Multiples of a random distribution is not a display convenience but a normalisation with an exact consequence: since 1 m.r.d. is by definition what a texture-free material gives everywhere, the area-weighted mean over the hemisphere must be 1 for *every* texture, however sharp. The mean therefore carries no information, which is precisely why every feature of a pole figure is a departure from it — and why a figure whose mean is not 1 has not been normalised, and its numbers mean nothing outside itself.

The identity is checked here on a sharp single component and on a random texture at once, because holding for both is the whole claim. Note that it is an *area*-weighted mean: the unweighted average of the same grid is biased, because an equispaced grid on a hemisphere is not an equal-area one.

**Symbols**

- $m.r.d.$ &mdash; Multiples of a random distribution.


:::{dropdown} Setup (imports and object construction)

```python
from pytex.app import REGISTRY

def pole_figure(model):
    return REGISTRY.call(
        "texture.pole_figure",
        {
            "phase": {"builtin": "ni_fcc"},
            "model": model,
            "spread_deg": 10.0,
            "grain_count": 400,
            "halfwidth_deg": 10.0,
            "seed": 7,
            "pole": [1, 1, 1],
            "projection": "equal_area",
            "resolution_deg": 5.0,
        },
    )
```

:::

**Compute**

```python
result = [
    pole_figure('random')['data']['mean_mrd'],
    pole_figure('goss')['data']['mean_mrd'],
    pole_figure('fcc_rolling')['data']['mean_mrd'],
]
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `workbench-mrd-mean-is-one` | [1.000, 1.000, 1.000] | [1.000, 1.000, 1.000] | m.r.d. | 2.97e-04 | 1e-02 | ✅ pass |

**Why this value**: Definitional: a pole figure normalised to multiples of a random distribution integrates to the sphere's area, so its area-weighted mean is exactly 1 for any texture.

**Citation**: Randle and Engler, Introduction to Texture Analysis, 2nd ed., chapter 5 (pole figure normalisation).

**See also**: {doc}`Texture foundation <../../concepts/texture_foundation>`, {doc}`Pole figure arithmetic and the m.r.d. scale <../../theory/pole_figure_arithmetic_and_mrd>`, {doc}`The PyTex Workbench <../../workflows/workbench_application>`

## A Miller label is a testable claim: Goss puts (011) on ND

The Goss component is written {011}<100>, and the first half of that notation asserts that the {011} plane lies in the sheet plane. So the (011) pole of a Goss texture must point along the sheet normal — polar angle zero, the centre of the pole figure. Checking it needs no reference figure at all, only the notation, which makes it the sharpest available end-to-end test of the whole chain: the Euler convention, the crystal-to-specimen mapping, the symmetry family of the pole, and the projection.

:::{dropdown} Setup (imports and object construction)

```python
from pytex.app import REGISTRY

def pole_figure(model):
    return REGISTRY.call(
        "texture.pole_figure",
        {
            "phase": {"builtin": "ni_fcc"},
            "model": model,
            "spread_deg": 6.0,
            "grain_count": 400,
            "halfwidth_deg": 10.0,
            "seed": 7,
            "pole": [1, 1, 1],
            "projection": "equal_area",
            "resolution_deg": 5.0,
        },
    )
```

:::

**Compute**

```python
rows = REGISTRY.call(
    'texture.pole_figure',
    {
        'phase': {'builtin': 'ni_fcc'}, 'model': 'goss', 'spread_deg': 6.0,
        'grain_count': 400, 'halfwidth_deg': 10.0, 'seed': 7,
        'pole': [0, 1, 1], 'projection': 'equal_area', 'resolution_deg': 5.0,
    },
)['table']['rows']
result = max(rows, key=lambda row: row['mrd'])['polar_deg']
```

**Result**

| Quantity | Computed (live) | Expected (reference) | Unit | Deviation | Tolerance | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `workbench-goss-pole-at-nd` | 0.0 | 0.0 | deg | 0.00e+00 | 6e+00 | ✅ pass |

**Why this value**: The Goss component {011}<100> places {011} in the sheet plane, so its (011) pole lies along ND at a polar angle of zero. The tolerance is the grid spacing, not a fitted margin.

**Citation**: Randle and Engler, Introduction to Texture Analysis, 2nd ed., chapter 5 (ideal orientations and their Miller descriptions).

**See also**: {doc}`Texture foundation <../../concepts/texture_foundation>`, {doc}`The PyTex Workbench <../../workflows/workbench_application>`
