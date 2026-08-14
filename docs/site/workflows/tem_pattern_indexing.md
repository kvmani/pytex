# Indexing A Pattern And Choosing Where To Tilt Next

This page is the microscope session, in order: what am I looking at, what is it,
where should I go next, and how do I get there. It is the workflow the TEM panel
of the workbench implements, and every step of it is available from Python
without the application.

The companion pages are
[Solving a measured SAED pattern](saed_pattern_solving.md) for the indexing
algorithm itself, [SAED generation](saed_generation.md) for the kinematic
simulation, and
[TEM specimen tilt navigation](../theory/tem_specimen_tilt_navigation.md) for the
holder geometry.

## 1. Open a pattern

Two kinds of pattern enter the workflow, and from the first click they are
treated identically.

**A micrograph of your own.** In the workbench the image never leaves the
browser; only the coordinates you click are sent. From Python, you supply the
coordinates directly.

**A practice plate from the gallery.** Three simulated patterns ship with the
application, chosen to be the three cases a microscopist meets first:

| Entry | Phase and axis | What it is for |
| --- | --- | --- |
| `fcc_al_001` | Aluminium down [001] | The square four-fold reference case for ratio-and-angle indexing; every reflection has unmixed indices. |
| `bcc_fe_110` | Ferrite down [110] | A centred rectangle whose two shortest vectors are perpendicular but in the ratio √2 — the pattern most often misread as a cubic ⟨001⟩ square. |
| `hcp_zr_2-1-10` | Zirconium down [2̄110] | A rectangle whose aspect ratio is √3·a/c, so it measures the axial ratio with no calibration at all. |

These are calculations, not pictures. Every spot sits where the lattice, the
zone axis and the camera constant put it. What they are *not* is a substitute for
a real plate: intensities are kinematic and therefore indicative, and double
diffraction is not modelled, so a reflection a real plate shows through double
diffraction is absent here. Positions are exact; brightnesses are a guide.

The realism that *is* modelled is the realism that breaks workflows. The beam is
not always at the centre of the frame, the pattern is rolled about the beam so it
does not line up with the detector axes, and each spot carries a sub-pixel
centroiding scatter. A workflow that only works on an idealised pattern has not
been tested.

```python
from pytex.core.lattice import ZoneAxis
from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

image = synthesize_saed_image(
    phase,
    ZoneAxis(indices=(0, 0, 1), phase=phase),
    camera_constant_mm_angstrom=10.0317,   # 200 kV, 400 mm camera length
    raster=DetectorRaster(width_px=1024, height_px=1024, pixel_size_mm=0.024),
    in_plane_rotation_deg=17.0,
    position_jitter_px=0.9,
    rng_seed=1001,
)
```

### The pixel convention

The raster is the recorded image: column index increases along the detector's
`+X` axis and row index increases along its `+Y` axis, which is how a camera
reads out and how every viewer displays the result. **No handedness flip is
applied.** A picked `(column, row)` pair is therefore the detector `(X, Y)` pair
divided by the pixel pitch and offset to the beam centre, and a pattern
synthesized here is interpreted by the solver under exactly the convention it was
built with. The round-trip tests hold both sides to that one convention.

## 2. Calibrate

The camera constant is the single number that turns a picture into a
measurement:

$$ r = \frac{L\lambda}{d} $$

with $r$ the distance of a reflection from the transmitted beam. In the gallery
this is *computed*, not typed: the panel takes a camera length and an
accelerating voltage and forms $L\lambda$ from the relativistic electron
wavelength, so the number in the calibration field and the number the geometry
used cannot drift apart. Shorten the camera and more of reciprocal space fits on
the plate; lengthen it and the outer reflections leave the frame.

Note what the calibration does *not* affect. The ratio of two lengths on the same
plate, and the angle between them, are calibration-free — which is why the
ratio-and-angle method identifies a zone axis without it, and why the hcp
prism-zone aspect ratio measures $c/a$ on an uncalibrated instrument. The camera
constant enters only when an absolute spacing is wanted. That is exactly why a
wrong camera constant is dangerous: it produces a self-consistent pattern of the
wrong material rather than an obviously broken one.

## 3. Fit the lattice, and settle the beam centre

Before indexing, impose the one constraint the pattern already satisfies: its
spots lie on a plane lattice. `pytex.diffraction.lattice_fit.fit_planar_lattice`
does that, and two things follow.

**The beam centre is solved for rather than clicked.** Four or more spots
over-determine it, so least squares gives the centre that best explains all of
them at once. This is worth more than it sounds: the camera equation measures
every radius *from the beam*, so an error there biases every spacing in the same
direction and produces a self-consistent answer for the wrong material.

**A mis-picked spot becomes visible.** The workbench draws the fitted lattice
over the pattern as two families of ruled lines, with the two basis vectors as
labelled arrows from the beam to the picks that generate them. Move the spot an
arrow points at and every line in the grid turns with it — so the two picks worth
being careful about are the two the arrows are on, and a spot clicked one node
out, or on a dust particle, stops matching the grid.

The panel offers both a directional pad for nudging the beam by a chosen step and
a *Refine from the spots* button, because the two are needed for different
reasons. The refinement is better than any single click; but a centre wrong by an
exact lattice vector fits perfectly, and only a person looking at which spot is
brightest can settle that one.

The full method, its three failure modes and its two irreducible limits are in
[Fitting the pattern lattice, and scoring the solutions](../theory/lattice_fit_and_solution_scoring.md).
**This step is geometry, not indexing**: a lattice that fits says the picks are
mutually consistent, which is necessary for a correct indexing and far from
sufficient.

## 4. Pick and index

Click the transmitted beam first. It is not a reflection; it is the origin every
spot is measured from, so an error there biases every d-spacing in the pattern at
once.

Then pick reflections. Two non-collinear ones are the minimum that fixes a zone
axis, and *non-collinear* is the operative word: the brightest spots of a pattern
almost always include a Friedel pair, $g$ and $-g$, which is collinear through
the beam and cannot seed anything. The workbench's auto-pick walks the brightness
order but skips a spot whose direction from the beam duplicates one already
taken; `SyntheticSAEDImage.independent_seed_spots` is the same rule in the
library.

Indexing itself is geometric and never uses intensity: see
[Solving a measured SAED pattern](saed_pattern_solving.md).

### Checking the answer

When the axis is already known — a practice plate, a reference specimen, a
pattern indexed before — pass it as `expected_zone_axis` and the result states
whether the indexing agrees and by how many degrees. **The comparison is made up
to symmetry.** A bcc [110] pattern is indistinguishable from a [101] one, because
the crystal symmetry maps one onto the other, so comparing index triples would
call a correct answer wrong. What is compared is the smallest angle between the
two directions over the symmetry orbit.

### Reading the candidates

Every candidate carries three things beyond its indices.

**Deviations**, which are measurements: the measured d-spacing against the
calculated one for each spot, and the measured angle against the calculated one
for each pair. The *same* relative deviation on every spot is the signature of a
wrong camera constant; a scatter of them is the signature of a wrong indexing.

**A fused accuracy score**, which is a policy — length agreement, angle agreement
and coverage, weighted and combined into one number in `[0, 1]`. The weights are
configurable, documented, and travel with every score, because a number whose
policy is invisible is an assertion rather than a measurement. Angles are
weighted above lengths by default: a wrong camera constant scales every length
and leaves every angle alone, so an angular disagreement is evidence about the
crystallography while a length disagreement may only be evidence about the
instrument. Candidates are ranked by this score rather than by the solver's own
sort key, and when the two orders disagree the result says so — that disagreement
is itself a sign the pattern does not settle the answer.

**Its calculated pattern**, in the picking coordinates. Superimposing what a
solution predicts on what was measured turns accepting it into a judgement made
by looking: a calculated pattern uniformly too large is a camera constant, one
turned is a roll, one with rows the plate does not show is the wrong phase. The
prediction is bounded by the index limit, so a plate spot with no calculated
ring beside it means *check the index limit* before doubting the solution.

Selecting a candidate draws it. **Accepting** one is a separate, deliberate act,
and it is what carries the phase and axis into the steps below — a tilt planned
from a solution nobody chose is a tilt planned from a guess.

## 5. Choose where to go next

Tilt planning answers *can I reach the axis I named*. At the column the question
usually comes the other way round, and `pytex.tem.atlas.zone_axis_atlas` answers
that one: which axes are worth going to at all.

```python
from pytex.tem.atlas import zone_axis_atlas

atlas = zone_axis_atlas(
    phase,
    current_zone_axis=ZoneAxis(indices=(0, 0, 1), phase=phase),
    max_index=2,
    max_angle_deg=60.0,
)
```

Each entry is one symmetry-distinct zone-axis family with the four things that
decide the choice:

- **how far** it is, measured to the nearest member of the family rather than to
  the one you typed, because every member gives the same pattern;
- **how many members** it has — a family of twelve offers twelve chances that one
  of them lies inside the holder's range;
- **how many reflections** its pattern shows inside a fixed cut-off, which is
  what the trip buys;
- **the rotational symmetry** of the pattern, which is what you recognise on the
  screen when you arrive, and the first confirmation that you arrived where you
  intended.

The symmetry is *measured on the simulated spot set*, not deduced from the point
group, so it reports what is actually there — including the Friedel centre of
symmetry a kinematic pattern always has, whether or not the crystal does.

The answer is usually not the nearest axis. From ferrite ⟨110⟩ the nearest family
within an index limit of 3 is ⟨320⟩ at 11.31° with eight reflections; ⟨111⟩ is
three times farther at 35.26° and shows thirty-six, with six-fold symmetry that
is unmistakable on arrival. The cost of a tilt is a few minutes and the risk of
losing the grain; the value is the information the new pattern carries.

The default index limit of 2 is deliberate: it gives the axes a standard
stereogram labels — ⟨100⟩, ⟨110⟩, ⟨111⟩, ⟨210⟩, ⟨211⟩, ⟨221⟩. Raising it admits
nearer but sparser families, which is a trade-off worth making consciously rather
than by default.

**Reflection counts are kinematic and exclude double diffraction**, so a real
plate of a diamond-structure or hexagonal phase shows a few more spots than the
count states. The count understates richness; it never invents it.

## 6. Tilt

The chosen destination goes to `plan_tilt_to_zone_axis`, which is documented in
[TEM specimen tilt navigation](../theory/tem_specimen_tilt_navigation.md). Two
things are worth repeating here because they surprise people:

**Reachability depends on the rotation about the beam, which one indexed pattern
cannot give you.** Every roll about the beam produces identical spot positions,
so a single pattern leaves it undetermined. Which member of a target family is
nearest, and how far the crystal must turn, do not depend on it; how that turn
divides between alpha and beta does, and therefore so does whether the holder can
make the move at all. Fix it from a second pattern at a different tilt.

**A zone axis corresponds to essentially one (alpha, beta) pair**, not to a menu
of routes that can be traded off. There is no combining a large alpha with a
large beta to reach further, which is why a ±30° holder cannot make a 54.7° move
however the tilts are divided.

## Exporting, at every stage

Every result in every panel exports in four formats, and the buttons are
generated from the manifest rather than written into the browser, so a format
added in Python appears everywhere at once.

| Format | For |
| --- | --- |
| **CSV** | One row per entity, at full precision, for another program. |
| **Excel** | The same table plus a sheet recording the inputs it came from. |
| **JSON** | The complete result, round-trippable back into the application. |
| **Report** | A readable Markdown page: the answer in prose, the caveats, the data, the exact inputs, and the citations. |

The first three are machine-readable to varying degrees; the last is the one to
paste into a notebook entry, because it says what was computed, from what, and on
whose authority. A result with no table still exports as a report — the prose and
the provenance are the point.

## What the numbers rest on

The three quantities this workflow depends on are pinned as executable worked
examples that recompute on every test run:
[simulated SAED plates and the zone-axis atlas](../examples/generated/saed_practice_patterns.md).
The camera-constant identity is checked against $L\lambda/d$ from the lattice
parameter; the hcp aspect ratio against $\sqrt{3}\,a/c$ from the hexagonal
reciprocal metric; and the basal-to-prism angle against the exactly 90° that the
hexagonal cell setting guarantees for any axial ratio.
