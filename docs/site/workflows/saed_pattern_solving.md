# Solving A Measured SAED Pattern

This workflow answers the question asked in front of a microscope: *here are the
spots I picked — which phase is this, down which zone axis, and what is every
spot?*

The input is deliberately minimal: spot positions relative to the transmitted
beam, plus enough calibration to turn them into reciprocal-space lengths. That
is all a printed pattern gives you, and it is all this needs.

## Which surface to use

PyTex has two indexing paths, and they solve different problems.

| | `pytex.diffraction.solving.solve_saed_pattern` | `pytex.diffraction.models.index_saed_pattern` |
| --- | --- | --- |
| Starts from | a list of picked spot coordinates plus a camera constant | a calibrated `DiffractionGeometry` — detector distance, pattern centre, tilt, pixel size |
| Works in | reciprocal angstroms | detector pixels |
| Suits | a person reading a micrograph or a printed pattern | an automated pipeline attached to an instrument |

Use `index_saed_pattern` when the detector model is known. Use
`solve_saed_pattern` when only the pattern is.

## The measured-pattern file

The YAML file is the reproducibility boundary. A pattern solved from a committed
file gives the same answer on any machine, whether the spots were clicked or
typed. It is validated by `schemas/measured_saed_pattern.schema.json`.

```yaml
schema: pytex.measured_saed_pattern/1
name: alpha_zone_01
calibration:
  units: px                 # px | mm | reciprocal_angstrom
  centre: [512.0, 512.0]    # transmitted beam, in the coordinates' own units
  pixel_size_mm: 0.014
  camera_constant_mm_angstrom: 180.0
  # or, equivalently: {camera_length_mm: 800, beam_energy_kev: 200}
spots:
  - {x: 612.0, y: 498.0, label: A}
  - {x: 545.0, y: 631.0, label: B}
```

**The transmitted beam is not a spot.** It is the calibration's `centre`, and
every position is taken relative to it; a spot sitting on the centre is
rejected, because it has no direction.

Calibration is validated at construction, not at first use: coordinates in `px`
or `mm` without a camera constant fail immediately, because an uncalibrated
length is not a recoverable state.

## Picking spots

`pytex.plotting.saed_picker.SAEDSpotPicker` displays an image and collects
clicks — left to add, right to remove the nearest, middle to set the beam
centre, `u` to undo, `c` to clear.

```python
picker = SAEDSpotPicker(image, calibration=calibration).show()
# ... click ...
picker.save_yaml("alpha_zone_01.yaml", name="alpha_zone_01")
```

The picking *logic* lives in `SpotPickerState`, a plain object with no
Matplotlib dependency, and the picker is a thin event adapter over it. That is
deliberate: an interactive tool that cannot be tested is a liability in a
scientific library, so the state machine is exercised headlessly and the GUI is
not on the critical path.

## Solving

```python
report = solve_saed_pattern(pattern, [alpha, beta], max_index=6)
print(report.describe())
best = report.best()
```

The algorithm is classical ratio/angle indexing. Two non-collinear reflections
fix the zone, so the two shortest measured vectors seed the solution; a
calculated pair is admissible when both lengths match within
`length_tolerance_relative` and their interplanar angle within
`angle_tolerance_deg`. The zone axis follows from the pair, the
crystal-to-pattern rotation from aligning calculated onto observed, and every
remaining spot is indexed by projection. Solutions are ranked by matched
fraction first, then by residual — a solution that explains every spot with
moderate residuals beats one that explains half of them perfectly, which is
usually a coincidence on a sub-lattice.

**Intensities are never used.** A kinematic intensity model is not reliable
enough to index against, and a printed pattern rarely carries calibrated
intensities at all. Geometry alone decides.

Systematic absences come from each phase's space group, so a phase supplied
without one is treated as primitive and may be offered reflections its real
structure forbids — the same trap described in
{doc}`composite_or_diffraction`.

## Reading the answer honestly

Three things the report tells you that a bare indexing would not:

**The zone-sense ambiguity is intrinsic.** A single SAED pattern cannot
distinguish a zone axis from its reverse when the reflection set is
centrosymmetric, because inverting the crystal leaves the pattern unchanged. The
report names this and does not count the two senses as competing answers.

**Symmetry-equivalent descriptions are one answer, not several.** Solutions are
deduplicated by crystal symmetry, and the survivor is rewritten into the
description a crystallographer would write — fewest negative indices, then
lowest — so a cubic cube-axis pattern reports $[001]$ rather than the equally
valid $[0\bar{1}0]$ that the seed search happened to find first.

**Not solving is a legitimate answer.** A report may hold no solutions, meaning
no candidate phase explains the pattern at these tolerances. `best()` raises
rather than returning a guess, and `is_conclusive` is `False` whenever the best
solution leaves spots unindexed or a genuinely different candidate explains the
pattern equally well.

An unindexed spot is often just an index bound: a reflection beyond the solver's
`max_index` is never offered a match, so raise it before widening tolerances.
`describe()` says so.

## Naming the transformation variant

When the pattern comes from a product phase and the parent's orientation in the
same pattern frame is known — from solving the parent's own spots, or from EBSD
— `assign_transformation_variant` names which variant it is: the child
orientation a variant predicts is $\mathbf{P}\mathbf{V}_k^{\mathsf{T}}$, and the
closest prediction wins. The reported deviation is symmetry-reduced under the
child point group, so a large value means the pattern does not belong to that
relationship at all — worth checking rather than assuming.

## Current limits

- The solver assumes the spots form a **zone-axis pattern about a low-index
  zone**. A crystal tilted off zone — for instance a transformation variant seen
  from a *parent* zone axis, whose own child zone axis is generally irrational —
  produces spots that do not all lie in one zero-order Laue zone, and will be
  only partly indexed. That partial match is the honest outcome; a full match
  would mean the solver was inventing reflections.
- Zero-order Laue zone only: no HOLZ rings, and no double-diffraction spots.
- Spot *detection* from image data is out of scope. The solver consumes picked
  or listed coordinates; `DiffractionPattern.cluster_observations` is the
  nearest available detection helper.

## Related Material

- {doc}`composite_or_diffraction` — simulating the patterns this solves
- {doc}`saed_generation`
- {doc}`../concepts/diffraction_foundation`
- {doc}`../concepts/orientation_relationships`

## References

### Normative

- {doc}`../architecture/diffraction_foundation`
- {doc}`../standards/data_contracts_and_manifests`

### Informative

- Edington, *Practical Electron Microscopy in Materials Science*, Monograph 2 —
  ratio/angle indexing of single-crystal patterns.
- Williams and Carter, *Transmission Electron Microscopy*, 2nd ed. — camera
  constant calibration and SAED indexing practice.
