# FIB Lamella Planning for a Target Zone Axis

**Workbench:** EBSD → **FIB lamella**. **Library:** `pytex.fib`. **Theory:**
{doc}`../theory/fib_lamella_zone_axis_geometry`. **Protocol:** {doc}`fib_azimuth_calibration`.

This page takes you from an EBSD scan of a bulk surface to the one page an operator carries to the
FIB: which grain to cut, where in the SEM field of view to put the lamella, at what azimuth, and what
tilt the TEM holder will then have to supply to reach the zone axis you want.

## When to use it

- You need a TEM lamella that shows a particular zone axis — `<011>` of an fcc phase for a
  dislocation study, `<0001>` of a hexagonal phase, a low-index axis for CBED — and you have an EBSD
  map of the surface you will cut from.
- You want to choose the grain *before* milling rather than tilting a finished lamella in the hope
  that the axis is within reach.
- You want the choice written down: the conventions, the calibration it relies on, and how sure it
  is.

It is not for tilted or wedge milling, bicrystal lamellae that need an axis in two crystals, or
serial-section data; those are outside version 1 and the module says so.

## The idea in one paragraph

A lamella milled straight down has its faces perpendicular to the surface, so the TEM beam — which
enters along the lamella normal — can only lie along a direction *in* the surface. For each grain
PyTex finds the symmetry equivalent of your zone axis that comes closest to the surface plane. The
angle by which it misses, **ε\***, is the tilt the holder must supply. Grains with a small ε\*, room
for the rectangle, and a well-defined orientation are the good sites.

## Step by step in the workbench

1. **Open the scan** in the EBSD workspace (any tab: the scan is shared) — an EDAX/TSL `.ang`, an
   Oxford `.ctf` or an EDAX `.oh5`/`.h5`. The file must already carry the 70° tilt correction; PyTex
   does not apply a second one. Without a file, the practice maps are used, drawn at a coarser step
   so that real lamellae fit.
2. Go to the **FIB lamella** sub-tab. It opens on its first example so there is always something to
   read.
3. **Target.** Enter the zone axis `[u v w]`. Any symmetry equivalent, in either sense, is
   accepted. To plan one grain you have already chosen, enter its id in *Plan one grain*;
   otherwise every grain is ranked.
4. **Frames and registration.** Say which way the scan's specimen `z` axis points (out of or into
   the surface) and which way its rows run. These two settings are vendor-dependent and no header
   states them; check them once per instrument. Choose how the SEM image is registered: identity
   (the SEM image *is* the scan), a declared rotation/scale/flip, or an affine fit to three or more
   control points — whose misfit then widens the uncertainty.
5. **Chamber.** Pick the ion-column angle (54° by default, or 52°) and enter the calibrated rotation
   sense and offset. Tick *Sense and offset come from the fiducial calibration* only after running
   {doc}`fib_azimuth_calibration`; until then every result carries an UNCALIBRATED warning.
6. **Lamella and holder.** Length, width (trench to trench), depth and final thickness; the holder's
   α and β limits (±30° by default) and the safety margin.
7. **Uncertainty and ranking.** EBSD accuracy, mount repeatability and the coverage factor; the
   ranking weights are in a closed section, stated back in the result.
8. Press **Plan the lamella**.

## Reading the result

The headline gives the grain, the residual tilt with its expanded uncertainty and feasibility class,
the lamella-normal azimuth in the sample frame and in the SEM image, and the FIB rotation with its
calibration status. Below it, in report order:

- **The recommended lamella** — the plan view (the grain outlined on the image-quality map, the
  rectangle, the normal, a scale bar and the sample axes) and the **stereogram**, which is the reason
  for the answer: its outer circle is the surface plane, every orbit member is a dot, and the ringed
  member's distance from the circle is ε\*.
- **The unknown mounting rotation** — a polar plot of the lamella's in-plane rotation φ on the grid.
  The shaded region is what the holder can remove at each φ; the arcs show where the target is
  actually reached. *Guaranteed* means reachable for every φ at the upper uncertainty bound with the
  margin to spare; *probabilistic* means for the stated fraction only.
- **The other candidates** and the **preparability map**: ε\* at every point, low is good, the
  planned rectangles numbered by rank.
- **Section and predicted pattern**, the **uncertainty budget**, **alternative azimuths** in the same
  grain, the **frames and calibration** in words, the **work order**, and the whole answer as prose.

The table lists every planned candidate with its residual, azimuths in all three frames,
feasibility, fraction of mounting rotations that work, clearance, spread and score, and exports to
CSV and XLSX like any result.

**Print work order** opens the one-page work order — grain, site coordinates in scan and image
units, stage tilt and rotation, the FIB rotation with its calibration caveat, the rectangle, the
expected TEM residual with its uncertainty, the feasibility class and every warning — ready for A4.
**Download work order** saves the same page.

## Doing the same from Python

```python
from pytex.adapters.scan_files import read_scan
from pytex.fib import ChamberGeometry, SurfaceGeometry, rank_grains, work_order_html

scan = read_scan("site.ang").dataset.crystal_map
report = rank_grains(
    scan.segment_grains(max_misorientation_deg=5.0),
    (0, 1, 1),
    surface=SurfaceGeometry(normal_sign=1, scan_y_sign=1),
    chamber=ChamberGeometry(column_angle_deg=54.0, rotation_sense=-1,
                            rotation_offset_deg=90.0, calibrated=True),
)
best = report.best()
print(best.describe())
open("work_order.html", "w", encoding="utf-8").write(work_order_html(best))
```

`report.to_json_dict()` validates against `schemas/fib_lamella_plan.schema.json`; every number in
`describe()` is also in the JSON.

## What the result does not know

- **The mounting rotation φ and the flip.** Read them off the first TEM pattern and solve the tilt
  with the TEM tilt-navigation workspace ({doc}`../theory/tem_specimen_tilt_navigation`).
- **The subsurface.** The grain is assumed columnar over the milling depth; check the lamella's first
  image for a buried boundary.
- **A real-instrument validation.** The chain is validated analytically and against PyTex's TEM
  solver, not yet against a lamella cut on a real instrument with its achieved tilt recorded.

## Examples

The panel ships four runnable examples, each executed by the test suite: a polycrystal ranked for
`<011>`; the same map for the harder `<111>`, where the mounting rotation starts to matter; one chosen
grain on a calibrated 52° instrument; and an affine SEM registration from five control points.
Executable worked examples with analytic answers are in {doc}`../examples/index`.
