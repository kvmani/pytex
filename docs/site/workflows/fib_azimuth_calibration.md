# Protocol: Calibrating the FIB Pattern Rotation (s_R and R0)

**Purpose.** Measure, once per instrument, the two numbers that turn a lamella azimuth found in
the sample frame into the rotation typed into the FIB pattern box,

$$
\theta_{\mathrm{ion}} = s_{R}\,(\theta_{S} + R_{\mathrm{st}}) + R_{0},
$$

where $s_{R} = \pm 1$ is the sense of the pattern rotation relative to the sample azimuth and $R_{0}$
is the offset. PyTex does not know either for any vendor and does not assume them: an uncalibrated
chamber gives every plan an UNCALIBRATED warning, because a wrong sign mills the lamella at the
mirror azimuth and wastes the site. See {doc}`../theory/fib_lamella_zone_axis_geometry`, section 7.

**When.** Before the first planned lamella on an instrument; again after any stage service, stage
exchange, software update, or change to how the pattern box's rotation is defined.

**Time.** About thirty minutes, most of it imaging.

## What you need

- A flat, polished specimen with an EBSD map (any material) or, failing that, any flat specimen whose
  SEM image you can register to the sample frame you use for planning.
- The same `SurfaceGeometry` and SEM registration you will use for planning. The calibration is only
  valid for that chain.

## Procedure

1. **Mount and record.** Mount the specimen as for a real lift-out. Note the stage rotation
   $R_{\mathrm{st}}$ you will mill at (keep it the same for all fiducials; a compucentric stage is
   assumed — otherwise re-find the area after any rotation).
2. **Tilt** the stage to the ion-column angle $T$ (52° or 54°), so the ion beam is normal to the
   surface.
3. **Mill three fiducial trenches** in an area with no other features, each a thin line box about
   10 µm × 0.5 µm and 0.5 µm deep, with pattern rotations $\rho_{1} = 0^{\circ}$, $\rho_{2} = 30^{\circ}$ and
   $\rho_{3} = 60^{\circ}$ as typed into the pattern box. A fourth at 100° adds redundancy. Do **not**
   use only 0° and 90°: those cannot tell the two senses apart, and the fit refuses them.
4. **Image** the trenches at zero tilt with the SEM (or map them with EBSD) in the frame you use for
   planning, through the same registration.
5. **Measure** for each trench the azimuth $\theta_{i}$ of its **normal** (perpendicular to the long
   side) in the sample frame, from $X_s$ towards $Y_s$. Only its value modulo 180° matters.
6. **Fit** with PyTex:

   ```python
   from pytex.fib import FiducialObservation, ChamberGeometry, calibrate_chamber_from_fiducials

   chamber = calibrate_chamber_from_fiducials(
       [
           FiducialObservation(pattern_rotation_deg=0.0, sample_azimuth_deg=92.1, stage_rotation_deg=0.0),
           FiducialObservation(pattern_rotation_deg=30.0, sample_azimuth_deg=61.8, stage_rotation_deg=0.0),
           FiducialObservation(pattern_rotation_deg=60.0, sample_azimuth_deg=32.3, stage_rotation_deg=0.0),
       ],
       base=ChamberGeometry(column_angle_deg=54.0),
   )
   print(chamber.rotation_sense, chamber.rotation_offset_deg, chamber.calibration_residual_deg)
   ```

   The fit tries both senses, takes the circular mean (period 180°) of the implied offsets, and keeps
   the sense with the smaller residual.
7. **Judge the residual.** An RMS residual below about 0.5° is a good calibration. Above 1–2°,
   suspect a misread trench, an uncorrected registration rotation, or a stage that moved between
   trenches, and repeat.
8. **Record** $s_{R}$, $R_{0}$, the residual, the date, the instrument and the software version in
   the instrument log. In the workbench, enter $s_{R}$ and $R_{0}$ under **Chamber** and tick
   *Sense and offset come from the fiducial calibration*.

## Checking a calibration before trusting it

Mill one more short trench using a *planned* azimuth — for instance the $\theta_{\mathrm{ion}}$ PyTex
gives for a trench normal at $\theta_{S} = 45^{\circ}$ — and measure its normal. It must come out at
45° within the residual. This is the check that would have caught a sign error before it reached a
real lamella.

## Failure modes

| Symptom | Likely cause |
| --- | --- |
| Measured normals change in the opposite sense to $\rho$ | $s_{R} = -1$: expected, the fit handles it |
| Residual ~90° on one trench | its long side and normal were swapped when measuring |
| Residual grows with each trench | the stage drifted or rotated between trenches |
| Calibration reproduces on one specimen but not another | the specimen was mounted in a different orientation relative to the EBSD reference; the chain, not the chamber, changed |
