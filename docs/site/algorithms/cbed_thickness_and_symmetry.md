# Reading A CBED Pattern: Thickness And Point Group

**Surface:** `pytex.diffraction.cbed.simulate_cbed_pattern`,
`thickness_from_fringe_minima`, `extinction_distance_angstrom`,
`fringe_minimum_excitation_errors`, `two_beam_rocking_curve`,
`holz_ring_radii_inv_angstrom`, `CBEDPattern.symmetry_observations`,
`CBEDPattern.determine_point_group`, with
`pytex.diffraction.diffraction_groups` carrying the group theory and
`pytex.diffraction.dynamical` the Bloch-wave solver.

A convergent-beam pattern answers two questions no other single measurement
answers as directly: **how thick is the foil here**, and **what is the point
group of the crystal, centre of symmetry included**. This page states how each
is computed, as steps a reader could reimplement, and what constrains the
answer. The derivations are in
{doc}`../theory/convergent_beam_electron_diffraction` and
{doc}`../theory/dynamical_cbed_and_symmetry_determination`; the numbers are
computed live in {doc}`../examples/index`.

Silicon runs through the page as the worked case, because it is the specimen the
technique is taught and calibrated on, and because both of its answers are known
in advance: the diamond glide extinguishes $(200)$ exactly, and $m\bar{3}m$ fixes
the diffraction group at every zone by group theory alone. An answer that can be
checked against something other than another run of this code is the only kind
worth demonstrating.

## 1. Why a convergent beam carries the information

A parallel beam sets one incident direction, so each reflection is one spot at
one excitation error. Converging the beam onto the specimen illuminates a *cone*
of incident directions at once. Every reflection becomes a **disc**, and
position within a disc is incident direction — so a disc is a continuous scan of
the rocking curve, recorded in one exposure.

Two consequences drive everything below:

- **Intensity varies across a disc**, and that variation is the thickness
  fringe pattern. Section 2 inverts it.
- **The pattern has a symmetry** that is a property of the crystal and the beam
  direction together, and — unlike kinematic diffraction — it *is* sensitive to
  the centre of symmetry. Section 3 reads it.

The disc half-angle is the convergence semi-angle $\alpha$; discs of radius
$\alpha/\lambda$ sit at the reciprocal-lattice positions. Whether neighbouring
discs overlap is the Kossel-Möllenstedt (separate) versus Kossel (overlapping)
distinction, and `ConvergentBeamConfig` reports which regime it is in rather
than leaving it to be inferred from the picture.

## 2. Thickness from the two-beam fringes

### 2.1 The physics, in one line

In the two-beam approximation the diffracted intensity oscillates with depth at
a rate set by the **effective** extinction distance, which stiffens as the
reflection is tilted away from the Bragg condition:

$$
I_g(s, t) = \frac{\sin^{2}(\pi t s_{\mathrm{eff}})}{(\xi_g s_{\mathrm{eff}})^{2}},
\qquad
s_{\mathrm{eff}}^{2} = s^{2} + \xi_g^{-2}.
$$

Dark fringes are the zeros of the numerator, at $t\,s_{\mathrm{eff},n} = n$ for
integer $n$. Both the foil thickness $t$ and the extinction distance $\xi_g$ are
unknown, and one fringe cannot separate them.

### 2.2 The Kelly linearisation

Squaring the fringe condition and dividing by $n^{2}$ removes the coupling:

$$
\left(\frac{s_n}{n}\right)^{2}
  = \frac{1}{t^{2}} - \frac{1}{\xi_g^{2}}\cdot\frac{1}{n^{2}}.
$$

Plot $(s_n/n)^2$ against $1/n^{2}$ and the points fall on a **straight line**
whose intercept is $t^{-2}$ and whose slope is $-\xi_g^{-2}$. One weighted least
-squares fit returns both unknowns, and the extinction distance comes out as a
by-product rather than having to be assumed — which matters, because $\xi_g$
depends on the accelerating voltage, the reflection, and the Debye-Waller factor,
and an assumed value propagates straight into the thickness.

This is Kelly *et al.* (1975), and `thickness_from_fringe_minima` implements it.

### 2.3 The algorithm

```text
input : s_1..s_N        (excitation errors of the dark fringes, any sign)
        first_order     (optional; the n to assign to the innermost fringe)

1  take |s_i|, sort ascending          -- innermost fringe gets the lowest order
2  for each candidate n0 in 1..max_first_order:      (skipped if first_order given)
3      assign orders n_i = n0, n0+1, ...
4      least-squares fit  y_i = (s_i/n_i)^2  against  x_i = 1/n_i^2
5      reject unless intercept > 0 and slope < 0      -- a physical fit
6      keep the assignment with the best R^2
7  t     = 1/sqrt(intercept)
8  xi_g  = 1/sqrt(-slope)
9  report t, xi_g, the orders used, R^2, and the residuals
```

Step 5 is what makes the search safe rather than a fishing expedition: a wrong
order assignment does not merely shift the line, it **curves** the plot, so the
non-physical fits are rejected on sign before $R^2$ is consulted.

### 2.4 The failure mode, stated plainly

The order assignment is the method's known weakness. If the innermost *visible*
minimum is not truly $n = 1$ — common when the disc is small or the foil thin,
because the $n=1$ fringe can lie outside the illuminated cone — then assuming
$n = 1$ inflates the thickness. The search in steps 2-6 exists for exactly this
case and usually recovers the right assignment from the curvature, but it is a
recovery, not a guarantee.

Two failures are refused outright rather than answered:

- **Fewer than two minima.** Two unknowns need two equations.
- **No physically valid assignment.** When nothing in $1..n_{\max}$ gives a
  positive intercept and a negative slope, the data are not two-beam fringes,
  and a `ValueError` is a more useful outcome than a plausible number.

### 2.5 What is not modelled

The two-beam expression is exact only when one reflection is strongly excited
and the rest are negligible. Near a zone axis that is false, and the fringe
positions shift. The Bloch-wave solver (`method="bloch"`) computes the
many-beam intensities, but `thickness_from_fringe_minima` still inverts the
two-beam relation — so a thickness read off a many-beam disc carries a
systematic error the fit cannot see. Take the fringes at a **two-beam
condition**, tilted away from the zone axis, which is what the technique
prescribes anyway.

Absorption damps the outer fringes and shifts nothing to first order, so it
degrades precision rather than accuracy.

## 3. Point group from the pattern's symmetry

### 3.1 What makes this possible at all

Kinematic diffraction obeys Friedel's law: $I_{\mathbf{g}} = I_{-\mathbf{g}}$
whether or not the crystal has a centre of symmetry. So a kinematic pattern
cannot distinguish the 32 point groups — only the 11 Laue classes. **Dynamical
scattering breaks Friedel's law**, and CBED is the standard method that
recovers the missing information. This is not a refinement of the kinematic
answer; it is a different observable.

The theory note derives the propagator theorem behind it. The consequence used
here is exact: an operation acting as $-1$ on both the beam direction and the
transverse plane is the inversion and nothing else, so Buxton's $2_R$ element is
present at **every** beam direction of a centrosymmetric crystal and at none of
an acentric one.

### 3.2 Diffraction groups: the forward direction

Fix the beam direction $\mathbf{b}$. Each operator $S$ of the crystal point group
either fixes $\mathbf{b}$, reverses it, or relates different patterns. The first
two contribute their restriction $T = S|_\perp$ to the plane normal to the beam,
the second tagged with a reciprocity flag because such an operation is a symmetry
only in combination with the reciprocity theorem. The resulting subgroups of
$G_2 \times \mathbb{Z}_2$ are **Buxton's 31 diffraction groups**, and
`pytex.diffraction.diffraction_groups` constructs them rather than transcribing
a table, then checks the count and membership against the published set.

Two of the three observables are rigid symmetries of one recorded pattern:

| Observable | Definition | Read from |
| --- | --- | --- |
| **whole pattern (WP)** | untagged elements only | the full pattern, disc centres and all |
| **bright field (BF)** | $\varphi(D)$, tagged elements acting as $-T$ | the transmitted disc alone |
| $\pm\mathbf{g}$ two-fold | tagged, compares $+g$ and $-g$ dark fields | *two exposures*; see 3.5 |

### 3.3 Silicon, and why the comparison is with gallium arsenide

Silicon is $m\bar{3}m$ (centrosymmetric); gallium arsenide is the same lattice
with a two-atom motif of different species, $\bar{4}3m$ (acentric). They are the
canonical pair because the *lattice* is identical, so anything that separates
them separates the centre of symmetry and nothing else. Computed from the group
theory alone:

| Zone | Si ($m\bar{3}m$) | GaAs ($\bar{4}3m$) |
| --- | --- | --- |
| $[001]$ | $4mm1_R$ — BF $4mm$, **WP $4mm$** | $4_Rmm_R$ — BF $4mm$, **WP $2mm$** |
| $[111]$ | $6_Rmm_R$ — BF $3m$, WP $3m$ | $3m$ — BF $3m$, WP $3m$ |
| $[110]$ | $2mm1_R$ — BF $2mm$, WP $2mm$ | $m1_R$ — BF $2mm$, WP $m$ |

Down $[001]$ the two crystals share a bright-field $4mm$ and differ in the
**whole-pattern** symmetry, $4mm$ against $2mm$. That single observation settles
centrosymmetry. Down $[111]$ both give BF $3m$ and WP $3m$ — the diffraction
groups differ ($6_Rmm_R$ against $3m$) but only in tagged elements, so the two
recorded observables do not separate them. **The zone axis is part of the
measurement**, and choosing one at which the candidate groups differ in BF or WP
is a step of the method, not a detail.

Note also the $[111]$ *projection* symmetry of silicon: $6mm$, higher than the
crystal's own $3m$. Projection symmetry is generally higher than the truth, which
is why section 3.4 refuses to report without higher-order Laue zones.

### 3.4 The algorithm

```text
input : phase, zone axis, ConvergentBeamConfig(method="bloch", laue_zones=(0,1,-1))

1  simulate the pattern              -- Bloch-wave, HOLZ beams included
2  build candidate plane operations from the pattern itself:
       rotations of order 2, 3, 4, 6
       mirror lines at the azimuths of the disc centres, and their perpendiculars
3  BF test : resample the transmitted disc under each candidate; keep those
             whose mismatch is under `tolerance` x the map's mean |deviation|
4  WP test : keep candidates that both permute the disc centres AND map each
             disc's intensity onto its partner's
5  close each surviving set under multiplication, and name it as a plane group
6  look up the diffraction groups consistent with (BF, WP)
7  report the point groups those imply, and whether all of them are centric
```

Step 2 matters more than it looks. Testing a dense sweep of angles would find
spurious mirrors in any smooth intensity map; taking mirror candidates only from
the disc-centre azimuths restricts them to orientations at which a mirror could
actually permute the discs. Step 5 is needed because $\{1, R_2, R_3, R_6\}$ is
four matrices but a six-fold group, and naming it before closure would understate
the symmetry.

### 3.5 Three refusals, and why each is right

The implementation declines to answer in three situations rather than returning
a number that would sometimes be wrong.

- **`method="two-beam"` is refused.** Each disc there is an independent rocking
  curve, symmetric in $s$ by construction, so every $\pm\mathbf{g}$ pair matches
  and the pattern reports a centre of symmetry whatever the crystal is. The
  symmetry would be the method's, not the specimen's.
- **A zeroth-Laue-zone-only beam set is refused** unless `require_holz=False`.
  Confined to the ZOLZ, the calculation samples the potential *projected* along
  the beam, and that projection is frequently centrosymmetric when the crystal is
  not — zincblende down $[111]$ is the standard example. Without HOLZ beams, GaAs
  down $[001]$ reports silicon's $4mm$ whole-pattern symmetry.
- **The $\pm\mathbf{g}$ two-fold is not measured**, and is passed through from
  the caller. Buxton's $2_R$ observation compares the $+g$ and $-g$ dark-field
  discs, each recorded with *its own* reflection at the Bragg condition — two
  exposures at different specimen tilts, related by reciprocity. It is not a
  two-fold rotation of a single zone-axis pattern. Treating it as one gives a
  test that provably fails: the excitation errors satisfy
  $s_{-g}(-\theta) - s_g(\theta) = -2g_z$, which vanishes only in the zeroth
  Laue zone, so once HOLZ beams are admitted the two-fold is broken for a
  centrosymmetric crystal too. The residual *grows* with the beam set for centric
  and acentric structures alike, so it is physics rather than truncation. The
  determination does not need it: BF and WP settle centrosymmetry at any zone
  whose candidate groups differ in them.

### 3.6 Tolerance

`tolerance` is the largest accepted mismatch as a fraction of each map's mean
absolute deviation, default $0.05$. Operations that are not grid-aligned —
three- and six-fold rotations and their mirrors — carry resampling error that
falls with `config.disc_samples`, so a trigonal or hexagonal zone needs at least
101 samples. The silicon-versus-GaAs separation down $[001]$ is not marginal:
residuals of $0.00$ against $0.32$.

## 4. HOLZ rings and the repeat along the beam

The first-order Laue zone appears as a ring of radius

$$
R_1 \approx \sqrt{\frac{2H}{\lambda}},
$$

with $H$ the reciprocal-lattice spacing along the beam, so measuring the ring
gives the **real-space repeat along the beam direction** — information a
zone-axis spot pattern does not carry at all, since it samples one plane of
reciprocal space. `holz_ring_radii_inv_angstrom` returns the radii for the
zones requested.

HOLZ *lines* inside the discs are sharp, and their intersections are sensitive
to lattice parameter at the $10^{-4}$ level. One trap is worth stating because it
is exact rather than approximate: **strain and accelerating voltage are
degenerate** in their effect on HOLZ line positions, so a lattice parameter
refined from HOLZ geometry is only as good as the voltage calibration.

## 5. Constraints and limits

| | |
| --- | --- |
| Two-beam thickness | needs a genuine two-beam condition; inverts the two-beam relation even on many-beam data |
| Order assignment | inferred from curvature; wrong assignment inflates thickness |
| Symmetry | needs `method="bloch"` and HOLZ beams; ZOLZ-only reports projection symmetry |
| $\pm\mathbf{g}$ | not measured; supply it from two-exposure experiment if wanted |
| Intensities | dynamical via Bloch waves, with an imaginary optical potential for absorption |
| Not modelled | inelastic background, detector response, beam divergence within a probe position |

## Verification

- Extinction distances of aluminium at 100 kV against the published table, and
  the Kelly plot recovering both $t$ and $\xi_g$ from fringe positions, in
  {doc}`../examples/generated/convergent-beam-diffraction`.
- The 31 diffraction groups constructed and checked against the published set;
  the $2_R \Leftrightarrow$ centrosymmetry correspondence verified over all 32
  point groups, in
  {doc}`../examples/generated/dynamical-cbed-and-symmetry`.
- Silicon's $(200)$ extinction and its diffraction groups against $\bar{4}3m$,
  in the silicon example of the same group.

## See also

- {doc}`../theory/convergent_beam_electron_diffraction` — disc geometry, the
  rocking curve, the absolute scale of $\xi_g$.
- {doc}`../theory/dynamical_cbed_and_symmetry_determination` — the coupled beam
  equations, absorption, and the construction of the 31 groups.
- {doc}`../tutorials/notebooks/28_convergent_beam_diffraction` — thickness
  measurement worked interactively, including the failure mode triggered on
  purpose.
- {doc}`../tutorials/notebooks/29_dynamical_cbed_and_point_groups` — the
  determination end to end, on silicon.

## References

### Normative

- Buxton, B. F., Eades, J. A., Steeds, J. W. & Rackham, G. M. (1976). The
  symmetry of electron diffraction zone axis patterns. *Philosophical
  Transactions of the Royal Society A* **281**, 171-194.
  <https://doi.org/10.1098/rsta.1976.0024>
- Kelly, P. M., Jostsons, A., Blake, R. G. & Napier, J. G. (1975). The
  determination of foil thickness by scanning transmission electron microscopy.
  *Physica Status Solidi (a)* **31**, 771-780.
  <https://doi.org/10.1002/pssa.2210310251>

### Informative

- Spence, J. C. H. & Zuo, J. M. (1992). *Electron Microdiffraction*. Plenum.
  <https://doi.org/10.1007/978-1-4899-2353-0>
- Williams, D. B. & Carter, C. B. (2009). *Transmission Electron Microscopy*,
  2nd ed. Springer. <https://doi.org/10.1007/978-0-387-76501-3>
