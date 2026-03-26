# EBSD: Kernel-Average Misorientation

Kernel-average misorientation (KAM) is a local orientation-gradient metric on an EBSD map. PyTex implements a regular-grid KAM workflow that is intentionally explicit about neighborhood order, thresholding, and symmetry handling.

## Definition

For an orientation field \(o_{i,j}\) on a regular grid and a neighborhood set \(N(i,j)\), PyTex uses the standard local-average form

```{math}
\mathrm{KAM}_{i,j} =
\frac{1}{|N(i,j)|}
\sum_{(k,l) \in N(i,j)}
\omega\!\left(o_{i,j}, o_{k,l}\right),
```

where \(\omega\) is the misorientation or disorientation angle, depending on whether symmetry-aware reduction is enabled.

## What PyTex Exposes

- regular 2D `CrystalMap` validation through `grid_shape`
- deterministic neighbor-pair generation for 4- and 8-connectivity
- neighbor order control for first-, second-, and higher-order neighborhoods
- optional misorientation thresholding
- symmetry-aware or symmetry-disabled evaluation
- mean-style and max-style aggregation

## Why The Parameters Matter

### Neighbor Order

Larger neighborhoods suppress some noise but also smooth away strongly local structure. PyTex therefore makes the order an explicit parameter rather than silently fixing it.

### Threshold

A threshold can suppress very large cross-boundary misorientations when the user wants an intra-grain local-gradient view instead of a boundary-dominated field.

### Symmetry Awareness

For crystallographic interpretation, the relevant quantity is often the disorientation angle after crystal symmetry reduction. For debugging imported data or checking raw rotational differences, symmetry-disabled behavior can also be useful.

![EBSD Grain Workflow](../../figures/ebsd_grain_workflow.svg)

## Example

```python
import numpy as np

from pytex import (
    CrystalMap,
    FrameDomain,
    Handedness,
    Orientation,
    OrientationSet,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)

orientations = OrientationSet.from_orientations(
    [
        Orientation(Rotation.identity(), crystal, specimen, symmetry=symmetry),
        Orientation(Rotation.from_bunge_euler(3.0, 0.0, 0.0), crystal, specimen, symmetry=symmetry),
        Orientation(Rotation.from_bunge_euler(6.0, 0.0, 0.0), crystal, specimen, symmetry=symmetry),
        Orientation(Rotation.from_bunge_euler(9.0, 0.0, 0.0), crystal, specimen, symmetry=symmetry),
    ]
)

crystal_map = CrystalMap(
    coordinates=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
    orientations=orientations,
    map_frame=specimen,
    grid_shape=(2, 2),
    step_sizes=(1.0, 1.0),
)

kam = crystal_map.kernel_average_misorientation_deg(
    order=1,
    threshold_deg=5.0,
    symmetry_aware=True,
    statistic="mean",
)
print(kam)
```

## Interpretation Notes

- high KAM often marks sharp local orientation gradients, but it is not itself a grain-reconstruction algorithm
- KAM responds strongly to noise if the orientation field is noisy
- KAM maps become much more interpretable when shown together with grain boundaries or after well-documented denoising

## Related Material

- `docs/architecture/ebsd_foundation.md`
- [../../tex/algorithms/ebsd_local_misorientation.tex](../../tex/algorithms/ebsd_local_misorientation.tex)
- [../../tex/algorithms/ebsd_kam_parameterization.tex](../../tex/algorithms/ebsd_kam_parameterization.tex)
- [../../figures/ebsd_grain_workflow.svg](../../figures/ebsd_grain_workflow.svg)

## References

### Informative

- MTEX documentation: [Kernel Average Misorientation (KAM)](https://mtex-toolbox.github.io/EBSDKAM.html)
- MTEX documentation: [Grain Orientation Parameters](https://mtex-toolbox.github.io/GrainOrientationParameters.html)
