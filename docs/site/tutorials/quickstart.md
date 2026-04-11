# Quickstart

This example constructs the minimal objects needed to work with orientations in PyTex.

The important point is not the matrix at the end. The important point is the semantic contract:

- the crystal frame is explicit
- the specimen frame is explicit
- the active symmetry model is explicit
- the orientation object retains that meaning instead of collapsing into an anonymous array

If those four pieces are missing, most texture and EBSD workflows become ambiguous at the first
tool boundary.

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

## What This Example Means

- `ReferenceFrame("crystal", ...)` defines the frame in which crystal directions, planes, and
  structure semantics live.
- `ReferenceFrame("specimen", ...)` defines the frame used by texture and EBSD-facing orientation
  results.
- `SymmetrySpec.from_point_group("m-3m", ...)` fixes the reduction and equivalence rules that will
  apply later.
- `Orientation(...)` is not just a rotation. It is a crystal-to-specimen mapping with attached
  frame and symmetry meaning.

## What To Check First In Your Own Data

Before building larger workflows, verify these points explicitly:

- Are your Euler angles really Bunge angles, and are they in degrees?
- Does your source define the mapping as crystal-to-specimen or specimen-to-crystal?
- Is the active symmetry a point-group reduction surface or a space-group structure surface?
- Do you need a single orientation or a vectorized `OrientationSet`?

If you cannot answer those questions directly, start with the concept pages before loading a full
dataset. PyTex is designed to force that ambiguity into the open early.

## Current Limits

- This quickstart uses a minimal orientation-only path. It does not yet introduce provenance,
  manifests, or batch workflows.
- The printed matrix is useful for inspection, but the semantically stable object is
  `Orientation`, not the bare matrix output.
- If your workflow starts from CIFs, EBSD imports, or diffraction geometry, use the dedicated
  workflow pages instead of extending this example ad hoc.

## Next Steps

1. Read {doc}`../concepts/core_model` and {doc}`../concepts/reference_frames_and_conventions` if
   you need the frame and mapping semantics first.
2. Use `OrientationSet.from_euler_angles()` when your input is naturally vectorized.
3. Use `Phase.from_cif(...)` or `Phase.from_cif_string(...)` when the workflow starts from
   crystallographic structure data rather than from orientations alone.
4. Move to {doc}`../workflows/ebsd_import_normalization` if your data comes from EBSD tooling.
5. Move to {doc}`../validation/index` if you want to inspect what PyTex has actually validated so
   far before building on a workflow surface.
