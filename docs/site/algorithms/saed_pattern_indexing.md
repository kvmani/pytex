# Solving A Measured SAED Pattern

**Surface:** `pytex.diffraction.solving.solve_saed_pattern`,
`solve_saed_pattern_file`, `assign_transformation_variant`, with
`pytex.plotting.saed_picker` for the picking front end.

Given spot positions relative to the transmitted beam and enough calibration to
scale them, this determines the phase, the zone axis, the crystal orientation in
the pattern frame, and the Miller indices of every spot — with residuals,
ranked alternatives, and an explicit verdict on whether the answer is
unambiguous.

```{figure} ../../figures/saed_indexing_algorithm.svg
:alt: Three-lane flow sheet. Lane 1 converts picked spot positions into
  reciprocal-space vectors through the camera constant. Lane 2 seeds from the two
  shortest non-collinear vectors, admits a calculated pair only when both lengths
  and the interplanar angle match, and builds the zone axis and rotation. Lane 3
  indexes every spot by projection, deduplicates symmetry-equivalent solutions,
  and ranks with a conclusiveness verdict.
:width: 100%

The algorithm, with the constraint governing each stage.
```

## 1. Which surface, and why there are two

| | `solve_saed_pattern` (here) | `models.index_saed_pattern` |
| --- | --- | --- |
| Starts from | picked spot coordinates plus a camera constant | a calibrated `DiffractionGeometry` |
| Works in | reciprocal ångström | detector pixels |
| Knows | nothing about the detector | distance, pattern centre, tilt, pixel size |
| Suits | a person reading a micrograph | a pipeline attached to an instrument |

Use `index_saed_pattern` when the detector model is known; use this one when only
the pattern is. They are not redundant — they consume different evidence.

## 2. Calibration

A picked position becomes a reciprocal-space vector through the camera constant
$L\lambda$:

$$\lVert \mathbf{g} \rVert = \frac{r}{L\lambda}, \qquad d = \frac{1}{\lVert \mathbf{g} \rVert},$$

with $r$ measured from the transmitted beam. Pixels scale by the pixel pitch
first. The camera constant may be given directly or derived from a camera length
and accelerating voltage through the relativistic electron wavelength.

The **transmitted beam is not a spot** — it is the calibration's centre, and a
spot coinciding with it is rejected, because it has no direction. Coordinates in
pixels or millimetres without a camera constant are rejected **at construction**,
not at first use: an uncalibrated length is not a recoverable state, so failing
early is the only honest option.

## 3. Seeding: ratio and angle

Two non-collinear reflections fix a zone, so the two shortest non-collinear
measured vectors seed the solution — shortest because they are the best
determined relative to picking error.

For a candidate phase, all reflections its lattice centring allows are enumerated
to `max_index`. A calculated pair $(\mathbf{g}_1^{c}, \mathbf{g}_2^{c})$ is
**admissible** for an observed pair $(\mathbf{g}_1^{o}, \mathbf{g}_2^{o})$ when
both tests pass:

$$\frac{\bigl| \lVert\mathbf{g}_i^{c}\rVert - \lVert\mathbf{g}_i^{o}\rVert \bigr|}
{\lVert\mathbf{g}_i^{o}\rVert} \le \varepsilon_{\text{len}}, \qquad
\bigl| \theta^{c} - \theta^{o} \bigr| \le \varepsilon_{\text{ang}} ,$$

where $\theta$ is the interplanar angle. This is the classical ratio/angle
indexing test (Edington; Williams and Carter), evaluated as a single pairwise
cosine matrix over the two admissible reflection pools rather than a loop.

The defaults are calibrated, not arbitrary:

- $\varepsilon_{\text{len}} = 0.03$ — three percent covers the centring error of
  a hand-picked spot at typical camera constants while still separating the
  $\{111\}$ and $\{200\}$ rings of an fcc metal, which differ by 15%.
- $\varepsilon_{\text{ang}} = 2°$.

:::{admonition} Constraint: intensities are never used
:class: important

A kinematic intensity model is not reliable enough to index against, and a
printed pattern rarely carries calibrated intensities at all. **Geometry alone
decides.** Intensity is carried through for plotting and record-keeping only.
:::

## 4. Zone axis and orientation

The zone axis follows from the cross product of the two calculated reflections,
$\mathbf{z} \propto \mathbf{g}_1^{c} \times \mathbf{g}_2^{c}$, converted to a
direct-lattice direction and rationalized.

The crystal-to-pattern rotation is built by matching two right-handed orthonormal
triads. Gram–Schmidt on the calculated pair gives
$\mathbf{E} = [\hat{\mathbf{e}}_1\;\hat{\mathbf{e}}_2\;\hat{\mathbf{e}}_1\times\hat{\mathbf{e}}_2]$
and on the observed pair (embedded in the detector plane, $z=0$) gives
$\mathbf{F}$. Then

$$\mathbf{R} = \mathbf{F}\,\mathbf{E}^{\mathsf{T}} ,$$

which is a proper rotation by construction — no re-orthogonalization, no
determinant repair.

Every allowed zone reflection is then projected through $\mathbf{R}$ and each
measured spot claims its nearest free prediction within its match radius
$\varepsilon_{\text{len}} \lVert \mathbf{g} \rVert$. Reflections are not reused,
so two spots can never claim one reflection.

## 5. Ranking, deduplication, and the verdict

Solutions are ranked by **matched fraction first**, then by mean residual. That
ordering is deliberate: a solution explaining every spot with moderate residuals
is a better answer than one explaining half of them perfectly, which is usually a
coincidence on a sub-lattice.

:::{admonition} Symmetry-equivalent descriptions are one answer, not several
:class: tip

Many seed assignments are related by a crystal symmetry operation and give the
same physical answer through different bookkeeping. Two solutions are the same
when their rotations differ by an element of the point group:
$\mathbf{R}_1 \simeq \mathbf{R}_2 \mathbf{S}$.

They are deduplicated on that test, and the survivor is rewritten into the
description a crystallographer would write — fewest negative indices, then lowest
— so a cubic cube-axis pattern reports $[001]$ rather than the equally valid
$[0\bar{1}0]$ the seed search happened to find first. Without this, an
unambiguous solve reported five competing "100% matched" solutions and
`is_conclusive` was `False`.
:::

:::{admonition} The zone-sense ambiguity is intrinsic, not a failure
:class: warning

A single SAED pattern **cannot** distinguish a zone axis from its reverse when
the reflection set is centrosymmetric: inverting the crystal leaves the pattern
unchanged (Friedel). The two senses are genuinely different proper rotations that
index equally well, and the report names the ambiguity rather than presenting one
sense as the answer. Resolving it needs a second zone axis or a
convergent-beam/dynamical observation.
:::

`is_conclusive` requires the best solution to index **every** spot and to face no
genuinely different competitor that does the same. The zone-sense pair is not
counted as a competitor.

## 6. Worked behaviour

### Cubic: simulate, then solve

An fcc austenite pattern simulated at 200 kV with a 180 mm·Å camera constant, its
spot positions handed back as if picked, and solved with both fcc and bcc offered
as candidates:

| zone | spots | recovered | indexed | mean residual (Å⁻¹) | conclusive |
| --- | --- | --- | --- | --- | --- |
| $[001]$ | 20 | $[001]$ | 100% | $3\times10^{-16}$ | yes |
| $[011]$ | 34 | $[011]$ | 100% | $6\times10^{-16}$ | yes |
| $[111]$ | 12 | $[111]$ | 100% | $6\times10^{-16}$ | yes |
| $[112]$ | 16 | $[112]$ | 100% | $7\times10^{-18}$ | yes |
| $[013]$ | 12 | $[013]$ | 100% | $10^{-15}$ | yes |

The residuals are at the floating-point floor because the round trip is exact by
construction. The bcc decoy is rejected on systematic absences.

Down $[001]$ the solver's own output reproduces the cubic identities
$\lVert\mathbf{g}_{220}\rVert / \lVert\mathbf{g}_{200}\rVert = \sqrt{2}$ exactly
and a $(200)$–$(220)$ angle of exactly $45°$ — independent geometry, not stored
values.

### Degradation under picking noise

Gaussian noise added to each picked position, fcc $[011]$:

| noise | outcome |
| --- | --- |
| 0.0 mm | 100% indexed, zero residual |
| 0.5 mm | 100% indexed, residual $0.0077$ Å⁻¹ — the noise, visible |
| 2.0 mm | **no solution** |
| 25.0 mm | **no solution** |

The failure mode is refusal, not a confident wrong answer. `best()` raises rather
than guessing, and `describe()` says no candidate phase explains the pattern at
these tolerances.

### Hexagonal: four-index labels

$\alpha$-Ti down $[0001]$, 54 spots, all indexed, with labels emitted in
four-index Miller-Bravais form — $(\bar{2}110)$, $(\bar{1}\bar{1}20)$,
$(\bar{1}2\bar{1}0)$ — because a three-index hexagonal label hides the symmetry
of the family.

## 7. Naming the transformation variant

When the pattern comes from a product phase and the parent's orientation in the
same pattern frame is known, `assign_transformation_variant` names the variant.
Each variant predicts the child orientation
$\mathbf{P}\mathbf{V}_k^{\mathsf{T}}$, and the closest prediction wins, with the
distance symmetry-reduced under the child point group:

$$k^{\star} = \arg\min_k \;\min_{\mathbf{S}_c \in G_c}
\angle\!\left( \mathbf{S}_c \mathbf{P}\mathbf{V}_k^{\mathsf{T}},\;
\mathbf{R}_{\text{solved}} \right).$$

A large deviation means the pattern does not belong to that relationship at all —
worth checking rather than assuming.

## 8. Constraints and limits

:::{admonition} Constraint: a zone-axis pattern is assumed
:class: warning

The spots must lie in one zero-order Laue zone. A crystal tilted off zone — for
instance a transformation variant seen from a **parent** zone axis, whose own
child zone axis is generally irrational — produces spots that do not all lie in
one ZOLZ and is only **partly** indexed. That partial match is the honest
outcome; a full match there would mean the solver was inventing reflections, and
a test pins the partial result rather than leaving it to be discovered.
:::

- Systematic absences come from each phase's space group, so a phase supplied
  without one is treated as primitive — the same trap as in simulation, with the
  same consequence.
- An unindexed spot is often just an index bound: a reflection beyond `max_index`
  is never offered a match. `describe()` names `max_index` as the first thing to
  raise, before the tolerances.
- No HOLZ rings and no double-diffraction spots.
- Spot **detection** from image data is out of scope. The solver consumes picked
  or listed coordinates.

## 9. Reproducibility: the file is the boundary

Spots may be clicked, but the measured-pattern YAML — validated by
`schemas/measured_saed_pattern.schema.json` — is what makes a solve reproducible:
a pattern solved from a committed file gives the same answer on any machine.

The picking *logic* lives in `SpotPickerState`, a plain object with no Matplotlib
dependency, and `SAEDSpotPicker` is a thin event adapter over it. An interactive
tool that cannot be tested is a liability in a scientific library, so the state
machine is exercised headlessly and the GUI is not on the critical path.

## Verification

| Claim | Where it is checked |
| --- | --- |
| Simulate-then-solve closure on five fcc zones and hexagonal $\alpha$-Ti | `tests/unit/test_saed_solving.py` |
| Cubic $\sqrt{2}$ ratio and $45°$ angle read back from the solver's output | same |
| fcc not solved as bcc; true phase outranks the decoy | same |
| An unexplainable pattern returns no solution; `best()` raises | same |
| Noise robustness at 0.5 mm, degradation beyond | same |
| Off-zone variant only partly indexed | same |
| Variant assignment recovers the planted variant | same |
| Calibration in all three unit systems; YAML round trip and schema | same |
| Live numerical demonstration | [worked examples](../examples/generated/composite-diffraction.md) |

## See also

- {doc}`composite_saed_assembly` — the forward problem
- {doc}`../workflows/saed_pattern_solving` — how to drive it
- {doc}`../tutorials/notebooks/23_transformation_crystallography_end_to_end`

## References

### Normative

- {doc}`../architecture/diffraction_foundation`
- {doc}`../standards/data_contracts_and_manifests`

### Informative

- Edington, *Practical Electron Microscopy in Materials Science*, Monograph 2 —
  ratio/angle indexing of single-crystal patterns.
- Williams and Carter, *Transmission Electron Microscopy*, 2nd ed. — camera
  constant calibration and indexing practice.
