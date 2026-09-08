# Navigating To A Zone Axis On A Real Holder

**Surface:** `pytex.tem.plan_tilt_to_zone_axis`, with `CurrentState` for the
orientation, `DoubleTiltStage` and the `TiltEnvelope` family for the holder,
`analyze_ambiguity` for identifiability, and
`pytex.plotting.tilt_stereogram.plot_tilt_stereogram` for the figure.

Navigating between crystallographic zone axes in transmission electron microscopy
requires determining the goniometer tilt angles $(\alpha, \beta)$ that bring a
target zone axis into coincidence with the electron beam. This page outlines
stage kinematics, reachability within coupled tilt envelopes, closed-form
inversion, and orientation determination from single patterns versus two-zone
orientations.

The theory note is {doc}`/theory/tem_specimen_tilt_navigation`; the architecture
foundation with full identifiability analysis is
{doc}`../architecture/tem_tilt_navigation_foundation`.

## 1. Physical and mathematical considerations

While ideal geometric alignment can be inverted analytically, practical tilt
planning is governed by three physical constraints:

| Nominal problem | Physical reality |
| --- | --- |
| Two-angle inversion | Multiple mathematical solution branches, crystal symmetry equivalents (up to 48 in cubic crystals), and physical ambiguity families. |
| Zone-sense ambiguity | Friedel's law introduces antipodal symmetry; determining true orientation additionally requires accounting for microscope lens rotation (diffraction rotation). |
| Tilt stage limits | Real pole-piece clearances impose coupled $(\alpha, \beta)$ envelopes rather than independent Cartesian angular limits. |

## 2. Coordinate frame definitions

In PyTex TEM workflows, the **holder coordinate system serves directly as the
specimen-domain reference frame**, with the crystal-to-holder rotation $\mathbf{U}$
represented as a standard `Orientation`.

Separating the specimen frame from the holder frame by an intermediate mounting
rotation introduces an unobservable degree of freedom: diffraction measurements
observe only the combined crystal-to-holder rotation $\mathbf{U}$. Defining
$\mathbf{U}$ directly eliminates parameter redundancy and ensures traceability
throughout tilt planning.

## 3. Stage kinematics

The $\alpha$ axis is the holder rod, fixed in the microscope column. The $\beta$
axis is a cradle carried *inside* the rod, so its orientation rotates with $\alpha$.
Composing the $\beta$ rotation about its instantaneous laboratory axis with $\alpha$:

$$\mathbf{R}_{\text{stage}} = \bigl[\mathbf{R}_x(\alpha)\mathbf{R}_y(\beta)\mathbf{R}_x(\alpha)^{\mathsf T}\bigr]\mathbf{R}_x(\alpha) = \mathbf{R}_x(\alpha)\,\mathbf{R}_y(\beta),$$

so the moving-axis factors cancel algebraically. This simplification is specific
to an orthogonal double-tilt gimbal where the inner $\beta$ axis rotates about
the outer $\alpha$ axis. For non-orthogonal or non-ideal stage configurations,
the forward kinematics must be composed explicitly via rotation operators.

Transposing for the beam axis gives the form everything else uses:

$$\hat{\mathbf{b}}_H(\alpha,\beta) = \bigl(-\cos\alpha\sin\beta,\ \sin\alpha,\ \cos\alpha\cos\beta\bigr).$$

This is spherical coordinates whose **pole is the $\beta$ axis**, and three
useful things follow at once.

- The Jacobian is $\cos\alpha$, so a tilt rectangle reaches
  $\Omega = \Delta\beta\,(\sin\alpha_{\max} - \sin\alpha_{\min})$ steradians.
  A $\pm30°$ holder gets $\pi/3 = 1.047\ \mathrm{sr}$ — **8.3 % of all beam
  directions.** That number is why symmetry equivalents matter so much.
- Constant-$\beta$ curves are great circles through the $\beta$ pole and
  constant-$\alpha$ curves are small circles about it, so the reachable region
  draws as exact arcs rather than a sampled outline.
- $\lVert\partial\hat{\mathbf{b}}/\partial\beta\rVert = \cos\alpha$: a degree of
  $\beta$ buys only $\cos\alpha$ degrees of crystal rotation. This is the
  double-tilt holder's gimbal lock, and the origin of the conditioning factor
  each solution reports.

## 4. The closed form

With $\hat{\mathbf{w}} = \mathbf{U}\hat{\mathbf{t}}_C$ the target in holder
coordinates,

$$\beta^{*} = \operatorname{atan2}(-w_1, w_3), \qquad \alpha^{*} = \operatorname{atan2}\!\left(w_2, \sqrt{w_1^2+w_3^2}\right).$$

Four branches satisfy alignment; two need $|\beta| > 90°$ and one needs about
$180°$ of $\alpha$, so on a real holder exactly one is relevant. The engine
enumerates all four anyway, so a reachability report can say *why* nothing is
available rather than returning an empty list.

**Degeneracy.** If $w_1 = w_3 = 0$ the target lies along the $\beta$ axis,
$\beta$ is indeterminate, and $\alpha = \pm90°$. Detected and reported as a
one-parameter family, never as an arbitrary `atan2(0, 0)`.

A calibrated non-ideal stage has no closed form, so the analytic result seeds a
Gauss–Newton refinement. Because the seed is analytic and the residual is a
physical angle, a convergence failure is diagnostic rather than a numerical
accident.

## 5. Reconstruction: prefer two zones to one pattern

| | Single pattern | **Two zone axes** |
| --- | --- | --- |
| Needs the diffraction rotation | **yes** | no |
| Needs the image parity | **yes** | no |
| Needs a detector model | no | no |
| Self-checking | `det U = +1` only | interzonal-angle test |
| Constructor | `CurrentState.from_pattern_solution` | `CurrentState.from_two_zone_axes` |

The two-zone approach relies strictly on crystallographic zone indices and recorded
stage tilt readouts: $\mathbf{U}\hat{\mathbf{n}}_i = \hat{\mathbf{b}}_H(\alpha_i,\beta_i)$,
solved as a classic two-vector attitude determination problem. This method bypasses
the need for the microscope diffraction lens rotation calibration constant.

Furthermore, it provides an intrinsic physical consistency check: the angle between
the two crystallographic zone axes $\hat{\mathbf{n}}_1 \cdot \hat{\mathbf{n}}_2$ is
fixed by lattice geometry and must equal the angle between the two calibrated beam
directions $\hat{\mathbf{b}}_{H,1} \cdot \hat{\mathbf{b}}_{H,2}$. A discrepancy
between these angles immediately diagnoses stage miscalibration, zone indexing errors,
or localized specimen bending.

A residual two-fold ambiguity persists when inverting both zone senses simultaneously,
corresponding to a $180°$ rotation about $\hat{\mathbf{n}}_1\times\hat{\mathbf{n}}_2$.
This ambiguity is physically inconsequential whenever the two-fold operation coincides
with an element of the crystal point group, which the engine evaluates and reports.
For example, between cubic $[001]$ and $[110]$, the cross-product axis is
$\langle 1\bar{1}0\rangle$, which is an exact diad in $m\bar{3}m$.

## 6. Identifiability: three separate ambiguities

### 6.1 Friedel/Laue — usually harmless, and the engine says so

Kinematic intensities obey $|F(\mathbf{g})| = |F(-\mathbf{g})|$, so the pattern
is centrosymmetric even when the crystal is not, and the reconstruction is
determined only up to the rotations of the **Laue class** that map the zone
plane to itself.

Enumerated over all 32 point groups, that group exceeds the crystal's own proper
group for exactly **ten** of them — those with improper operations other than
inversion — and by a factor of exactly two:

| Class | Point groups | Index |
| --- | --- | --- |
| Centrosymmetric (11) | $\bar{1}$, $2/m$, $mmm$, $4/m$, $4/mmm$, $\bar{3}$, $\bar{3}m$, $6/m$, $6/mmm$, $m\bar{3}$, $m\bar{3}m$ | 1 |
| Enantiomorphic (11) | $1$, $2$, $222$, $4$, $422$, $3$, $32$, $6$, $622$, $23$, $432$ | 1 |
| Improper, acentric (10) | $m$, $mm2$, $\bar{4}$, $4mm$, $\bar{4}2m$, $3m$, $\bar{6}$, $6mm$, $\bar{6}m2$, $\bar{4}3m$ | **2** |

So for every cubic and hexagonal metal the much-feared $180°$ ambiguity is
absorbed entirely by crystal symmetry, and the engine reports no ambiguity
rather than warning about one. Note also that **"non-centrosymmetric" is the
wrong test**: quartz ($32$) is acentric yet completely ambiguity-free.

### 6.2 The diffraction rotation — the one that costs a session

An error $\delta\varphi$ is a rotation about the beam axis, is **not** absorbed
by symmetry, and lands the specimen

$$\Delta = 2\arcsin\!\left(\sin\tfrac{\delta\varphi}{2}\sin\theta_c\right) \approx \delta\varphi\,\sin\theta_c$$

from the target, where $\theta_c$ is the angle between current and target zones.

| $\delta\varphi$ | 5° hop | 30° hop | 60° hop | 90° hop |
| --- | --- | --- | --- | --- |
| 2° | 0.17° | 1.00° | 1.73° | 2.00° |
| 5° | 0.44° | 2.50° | 4.33° | 5.00° |
| 180° | 10° | 60° | 120° | 180° |

Two consequences. **At $180°$ both tilt angles are negated**: the specimen goes
exactly the wrong way while the calculation reports a clean zero residual. And
the error **scales with the hop**, so long excursions are routed through
intermediate low-index zones with re-indexing at each — several short hops are
self-correcting where one long hop is open-loop.

### 6.3 Parity

Image mirroring and readout signs *reflect* the trajectory rather than rotating
it, so unlike §6.2 they cannot be absorbed into any angle. Caught by
$\det\mathbf{U} = -1$ and by the interzonal-angle test.

## 7. Calibration: two exposures

Apply a known positive $\alpha$, record the azimuth $\psi_\alpha$ at which a
tracked Kikuchi feature moved; repeat for $\beta$. A rigid rotation carries a
pole at the pattern centre along $-\hat{\mathbf{y}}_L$ for $+\alpha$ and
$+\hat{\mathbf{x}}_L$ for $+\beta$, so

$$\varphi_D = -\psi_\beta, \qquad \psi_\alpha - \psi_\beta = \begin{cases} -90° & \text{unmirrored} \\ +90° & \text{mirrored} \end{cases}$$

Three results from two exposures — the rotation numerically, the parity from a
sign, and a redundancy check, since the azimuths must be $90°$ apart. If they
are not, either the axes are non-orthogonal or the tracked feature was
misidentified, so the procedure diagnoses itself.

$\varphi_D$ is **hysteretic** in the lens settings, so a calibration carries the
camera length and voltage it was measured at and `check_applicable` refuses a
mismatch rather than interpolating.

## 8. Reachability: four distinct verdicts

| Verdict | Meaning |
| --- | --- |
| `EXACT` | Forward-validated, comfortably inside the envelope |
| `WITHIN_TOLERANCE` | Reachable, with a residual worth stating |
| `NEAREST_APPROACH` | No orbit member reachable; this is the closest — **and it is guaranteed to be a position the holder can actually reach** |
| `UNREACHABLE` | Not close enough for a partial tilt to help |

These are qualitatively different answers, not gradations of one, and the API
keeps them apart so a nearest approach cannot be mistaken for a hit. A
single-tilt holder reaches a set of *measure zero* — a great circle — so an
exact zone axis there is coincidence only, and the correct output is always a
nearest approach.

## 8a. The reverse direction: orientation *from* an indexed pattern

Everything above inverts the master equation for the tilts. The same equation
inverts for the **orientation**, which is what texture and microstructure work
wants out of indexing: given the crystal-to-pattern rotation that indexing
produced and the holder angles at which the pattern was recorded,

$$\mathbf{U} = \mathbf{R}_{\text{stage}}(\alpha,\beta)^{\mathsf T}\ \mathbf{F}\ \mathbf{R}_z(\varphi_D)\ \mathbf{R}_{P \leftarrow C}.$$

Because the holder frame is the specimen-domain frame, $\mathbf{U}$ *is* an
ordinary `Orientation` — reportable as Bunge $(\varphi_1, \Phi, \varphi_2)$ and
directly comparable with an EBSD measurement, with no convention conversion in
between.

**Surface:** `pytex.tem.orientation_from_indexed_pattern` for one pattern,
`orientation_from_indexed_patterns` for several, and
`orientations_from_pattern_report` to map over a solver's ranked solutions.

### Why one pattern is not enough

That composition needs $\varphi_D$, and an error in it is **not** absorbed by
crystal symmetry: it rotates the reported orientation bodily about the beam
axis, degree for degree. A $180°$ error yields a clean, self-consistent
orientation that is wrong by $180°$. The single-pattern surface therefore
**raises** rather than guessing when no calibration exists.

### Two patterns determine the orientation *and* the calibration

The way out needs no extra equipment. With two indexed patterns at different
stage positions:

1. The zone axes alone fix the orientation —
   $\mathbf{U}\hat{\mathbf{n}}_i = \hat{\mathbf{b}}_H(\alpha_i,\beta_i)$, two
   non-parallel correspondences, Wahba's problem. No calibration enters.
2. Each pattern's in-plane indexing then over-determines $\varphi_D$: the
   residual
   $\mathbf{M}_i = \mathbf{R}_{\text{stage},i}\,\mathbf{U}\,\mathbf{R}_i^{\mathsf T}$
   must be a rotation **about the beam axis**, and its angle *is* the diffraction
   rotation. It is read off, not searched for, so recovery is exact.

Two independent diagnostics fall out:

| Quantity | What it catches |
| --- | --- |
| Scatter of $\varphi_D$ across patterns | Patterns disagreeing about one instrument constant |
| **Beam-axis deviation of $\mathbf{M}_i$** | The part of the residual that **no** $\varphi_D$ can absorb — a mirrored stored pattern, a mis-indexed reflection, a reversed stage sign, or a bent specimen |

The second is the one to read. `MultiPatternOrientation.as_calibration()` refuses
to return a calibration when it is large, because propagating a constant the data
already contradict is worse than having none.

### How a mirrored pattern actually presents

Not as an improper matrix. Indexing builds right-handed triads from both the
observed and the calculated vectors, so it *always* returns a proper rotation.
Mirroring instead flips the stored $y$ axis **and** reverses the derived pattern
normal — together a $180°$ rotation about the pattern $x$ axis, which is proper.
It therefore slips past any determinant check and is caught only by the
beam-axis deviation above. An improper matrix arriving at this surface means a
caller composed a mirror in by hand, and is rejected as an input error.

## 9. Path planning: the geodesic is the Kikuchi band

A straight line in $(\alpha,\beta)$ is not a straight line on the crystal
sphere. The geodesic between two zone axes lies in the plane they span, whose
normal rationalizes to a low-index $(hkl)$ — and the Kikuchi band of that
reflection is exactly the band joining the two poles.

**The mathematically optimal path is the one experienced operators already
follow by eye**, so the planner names the band:

> *Follow the $(1\bar{1}0)$ Kikuchi band from $[001]$ toward $[111]$; total
> travel 54.7°.*

Every sampled point is checked against the envelope with a working clearance,
because a path that grazes a mechanical limit will hit it once backlash is
included.

## 10. Forward validation

No solution reaches the caller without being re-derived through the calibrated
forward model:

1. take the estimated $\mathbf{U}$;
2. build $\mathbf{R}_{\text{stage}}$ from the **calibrated** model;
3. compute $\mathbf{v}_L = \mathbf{R}_{\text{stage}}\mathbf{U}\hat{\mathbf{t}}_C$;
4. measure its angle from the beam axis;
5. **reject** above tolerance.

Rejection happens *before* ranking, never as a down-weight, so a solver bug
cannot present itself as a low-scoring answer. Note what the residual does and
does not mean: it measures how well the solver did, **not** whether the inputs
were true. A mis-calibrated plan reports a residual of zero. That is what the
propagated $\sigma$ values are for.

## 11. Complexity and failure modes

Solving is $O(F \times N \times 4)$ closed-form evaluations for $F$ ambiguity
families and $N$ orbit members — 192 for a cubic general direction, each a
handful of flops. The nearest-approach search is a vectorized grid over the
envelope. Neither is a bottleneck at interactive scale.

Failure modes, in order of practical likelihood:

1. **$\varphi_D$ wrong** — a clean, confident, exactly-backwards answer. Detected
   only by calibration or by the two-zone path.
2. **Parity or sign convention wrong** — reflected trajectory.
3. **Specimen bent or the area polycrystalline** — the rigid assumption fails;
   bounded by the interzonal-angle residual across pairs.
4. **Target orbit outside the envelope** — not a failure but an answer, and the
   most common useful output.
5. **Target near the $\beta$-axis pole** — ill-conditioned by $1/\cos\alpha$,
   bounded at 1.31 within a real envelope.
6. **One of the ten affected point groups with polarity-sensitive purpose** —
   irreducible from kinematic SAED. Reported, not guessed.

## See also

- {doc}`../architecture/tem_tilt_navigation_foundation` — the design record.
- {doc}`../tutorials/notebooks/24_tem_tilt_navigation` — the workflow end to end.
- {doc}`../examples/generated/tem_tilt_navigation` — computed, checked numbers.
- {doc}`saed_pattern_indexing` — producing the indexed pattern this consumes.
