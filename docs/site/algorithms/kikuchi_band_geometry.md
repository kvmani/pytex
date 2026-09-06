# Kikuchi Band Geometry And Pattern Simulation

**Surface:** `pytex.diffraction.kikuchi.simulate_kikuchi_pattern`,
`KikuchiPattern`, `KikuchiBand`, `KikuchiZoneAxis`, `GnomonicProjection`, with
`pytex.diffraction.kikuchi_map` for the stereographic band network and the
workbench operations `ebsd.simulate_kikuchi_pattern` and `ecci.solve_workflow`.

An EBSD pattern is the projection of the crystal's own diffraction geometry onto
a detector, and every orientation measurement in the technique is the inversion
of that projection. This page states the forward model — which planes give
visible bands, where they land, how wide they are — because the forward model is
what indexing inverts, what a calibration fits, and what a simulated pattern is
compared against.

```{figure} ../../figures/kikuchi_geometry_algorithm.svg
:alt: Three-lane flow sheet. Lane 1 enumerates candidate planes, reduces
  antipodal pairs, and filters by centring and a kinematic intensity threshold.
  Lane 2 carries each normal into the laboratory frame and computes the Bragg
  angle from the relativistic wavelength. Lane 3 projects gnomonically and
  reports band centre lines, edges at twice the Bragg angle, and zone axes.
:width: 100%

The forward model that orientation determination inverts.
```

## 1. Why bands, and why they are straight

A Kikuchi band comes from a lattice plane, not from a reflection. Inelastically
scattered electrons travel in all directions inside the crystal; those meeting a
plane at exactly the Bragg angle diffract. The locus of directions at a fixed
angle $\theta_B$ from a plane is a **cone** — the Kossel cone — and there are two
of them, one either side of the plane. A band is the pair of cone traces on the
detector, and its centre line is the trace of the plane itself.

Because $\theta_B$ is small for electrons (fractions of a degree to a few
degrees), the cones are very flat and their traces are nearly straight lines. In
**gnomonic** coordinates they are exactly straight, which is why that projection
is the natural frame for the whole subject.

### 1.1 The gnomonic projection

`GnomonicProjection` maps a direction to the point where the ray from the source
along it meets the detector, in units of the detector distance. Its defining
property:

> **Great circles map to straight lines.**

So a lattice-plane trace is a straight line in gnomonic coordinates and a curve
in pixel coordinates on a tilted detector. Two further properties matter in
practice:

- The origin is the **pattern centre**, and one gnomonic unit is one detector
  distance — so a coordinate of $1.0$ lies at $45^\circ$ from the detector
  normal.
- Coordinates are **detector-size independent**. The same crystal gives the same
  gnomonic pattern on any detector, which is what makes gnomonic space the right
  frame for comparing patterns, fitting a calibration, or transferring a
  solution between microscopes.

Only the forward hemisphere projects. Directions travelling away from the
detector, or parallel to its plane, have **no intersection**, and are reported
invalid rather than given a spurious coordinate — a silent wrap would place a
band on the opposite side of the pattern.

## 2. Band width measures the plane spacing, inversely

The band's angular width is exactly $2\theta_B$, and

$$
\sin\theta_B = \frac{\lambda}{2d},
$$

so **band width is a direct measurement of interplanar spacing** — an *inverse*
one. Wide bands are low-$d$ planes; narrow bands are high-$d$.

This inverts an intuition that catches people out and is worth stating
explicitly: `min_d_spacing_angstrom` excludes the **widest** bands, not the
narrowest, because it drops the small-$d$ planes. The weak high-order bands that
clutter a pattern are the wide ones.

$\lambda$ is the **relativistically corrected** electron wavelength. At 20 kV
the correction is about 2 %, and at 200 kV about 30 % — omitting it puts every
band width wrong by that factor, which a calibration then absorbs into a wrong
detector distance.

## 3. The simulation algorithm

```text
input : geometry (detector, energy, frames), phase, orientation,
        max_index, min_d_spacing, min_relative_intensity, max_bands

1  enumerate candidate planes (hkl) up to max_index
2  reduce to one representative per antipodal pair
       -- (hkl) and (-h-k-l) are the same plane and would draw the same band twice
3  apply the lattice-centring reflection condition
       -- a systematically absent reflection produces no band
4  compute kinematic |F|^2; drop bands below min_relative_intensity
5  carry each plane normal into the laboratory frame:
       n_lab = T_specimen->lab . R_crystal->specimen . n_crystal
6  Bragg angle from sin(theta_B) = lambda / 2d, relativistic lambda
       -- drop planes whose spacing cannot satisfy it at all
7  keep the strongest max_bands, if a cap is given
8  enumerate zone axes up to zone_axis_max_index; a band belongs to a zone
   when the zone law h u + k v + l w = 0 holds
```

Step 5 is where the frames must be right, and the implementation **checks rather
than assumes**: the geometry's specimen frame must match the orientation's, and
the orientation's crystal frame and phase are checked against the phase
argument. A frame mismatch here produces a plausible, wrong pattern — the most
expensive kind of error in this domain.

### 3.1 Zone axes organise the pattern

Bands intersect at **zone axes**: the direction $[uvw]$ common to every plane
whose zone law $hu + kv + lw = 0$ holds. On a pattern these are the bright
intersections where several bands cross, and they are what a human uses to index
by eye. `KikuchiZoneAxis` carries them with the bands that meet there, so a
solved pattern can name the zone rather than only the bands.

## 4. What the intensities are, and are not

The intensities are a **kinematic $\lvert F\rvert^2$ proxy**. That is enough to
decide which bands are present and roughly how strong, which is what a geometric
model needs. It is *not* a photometric prediction:

| Present in a real pattern | Modelled here |
| --- | --- |
| Band positions and widths | **yes, exactly** |
| Which planes give bands | yes, via $\lvert F\rvert^2$ and centring |
| Excess/deficiency asymmetry across a band | **no** |
| Dynamical contrast, band profile shape | **no** |
| Inelastic background | **no** |
| Detector response | **no** |

So this is a map of *where the bands are*, not a photograph of one. For a
dynamical treatment see {doc}`cbed_thickness_and_symmetry` and the Bloch-wave
solver; for the band network on the crystal sphere rather than a detector, see
the stereographic Kikuchi map.

## 5. Cost, and the parameter that controls it

Plane enumeration is $O(\text{max\_index}^3)$, so raising `max_index` admits
narrower higher-order bands at cubic cost. Everything downstream is vectorised
over the surviving planes. `zone_axis_max_index` is enumerated separately and is
usually kept lower, since high-index zones are not visible intersections in
practice.

| Parameter | Effect |
| --- | --- |
| `max_index` | more, narrower bands; cubic cost |
| `min_d_spacing_angstrom` | drops the **widest** bands (section 2) |
| `min_relative_intensity` | drops weak bands |
| `max_bands` | legibility cap, strongest kept |
| `zone_axis_max_index` | how many intersections are labelled |

## 6. What this is the forward model for

| Consumer | Inverts or uses it for |
| --- | --- |
| Hough/Radon indexing | detects band centre lines, matches against a simulated look-up |
| Dictionary indexing | compares a measured pattern against simulated ones |
| Pattern-centre calibration | fits detector distance and centre so simulated bands land on measured ones |
| ECCI workflow | predicts the pattern at a stage position, to find a two-beam condition |
| Teaching | the same geometry, drawn rather than solved |

Because the same simulation serves all of these, a change in the geometry model
propagates to indexing, calibration and the workbench together, and cannot drift
between them.

## Verification

- Band positions and widths against the analytic Bragg relation, and the
  gnomonic straight-line property, in {doc}`../examples/generated/ebsd`.

## See also

- {doc}`../theory/kikuchi_bands_and_gnomonic_projection` — the derivation, the
  Kossel cone, and the detector frames in full.
- {doc}`../theory/stereographic_kikuchi_maps` — the band network on the crystal
  sphere, for planning a tilt.
- {doc}`cbed_thickness_and_symmetry` — what a dynamical treatment adds.
- {doc}`ebsd_grains_and_local_misorientation` — what the indexed orientations
  become.

## References

### Normative

- Kikuchi, S. (1928). Diffraction of cathode rays by mica. *Japanese Journal of
  Physics* **5**, 83-96.
- Schwartz, A. J., Kumar, M., Adams, B. L. & Field, D. P., eds. (2009).
  *Electron Backscatter Diffraction in Materials Science*, 2nd ed. Springer.
  <https://doi.org/10.1007/978-0-387-88136-2>

### Informative

- Britton, T. B. et al. (2016). Tutorial: Crystal orientations and EBSD - or
  which way is up? *Materials Characterization* **117**, 113-126.
  <https://doi.org/10.1016/j.matchar.2016.04.008>
- Winkelmann, A., Trager-Cowan, C., Sweeney, F., Day, A. P. & Parbrook, P.
  (2007). Many-beam dynamical simulation of electron backscatter diffraction
  patterns. *Ultramicroscopy* **107**, 414-421.
  <https://doi.org/10.1016/j.ultramic.2006.10.006>
