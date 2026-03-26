# Quickstart

This example constructs the minimal objects needed to work with orientations in PyTex.

```python
from pytex import FrameDomain, Handedness, Orientation, ReferenceFrame, Rotation, SymmetrySpec

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)

orientation = Orientation(
    rotation=Rotation.from_bunge_euler(45.0, 35.0, 15.0),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
)

print(orientation.as_matrix())
```

## Next Steps

- Use `OrientationSet.from_euler_angles()` for vectorized orientation input.
- Use `Phase.from_cif(...)` or `Phase.from_cif_string(...)` when your phase definition starts from crystallographic structure data.
- Use `InversePoleFigure.from_orientations()` to reduce crystal directions into symmetry sectors.
- Use `CrystalMap.kernel_average_misorientation_deg()` for a minimal EBSD neighbor metric on regular grids.
- Use `CrystalMap.segment_grains()` for thresholded grain segmentation and `GrainSegmentation.grod_map_deg()` for GROD.
