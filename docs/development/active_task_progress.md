# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history. Governed by the cardinal rule in `AGENTS.md`: ledger plus commit-and-push to `main`
after every substantial increment.

## Stereographic Kikuchi Maps For Zone-Axis Navigation — CAPABILITY COMPLETE (2026-08-11)

**Objective.** Give PyTex the TEM operator's road atlas: the Kikuchi bands and zone axes of a
phase, drawn on a stereographic projection of the crystal sphere, with routing along bands
from the current zone axis to a target one. Then tutorial 30, with inline cubic and hexagonal
graphics.

**Why this is not already covered.** Three neighbouring surfaces exist and none of them is
this:

- `pytex.diffraction.kikuchi` simulates an EBSD-style pattern in **gnomonic** projection on a
  flat detector of finite extent. A gnomonic projection cannot show a hemisphere — a band at
  90 degrees from the pattern centre is at infinity — so it is structurally unable to be an
  atlas.
- `pytex.tem.path.connecting_band` already knows the key fact that the geodesic between two
  zone axes *is* a Kikuchi band, but it answers only for one pair of poles given in advance.
- `pytex.plotting.tilt_stereogram` draws poles and the stage envelope stereographically, but
  has no bands, so it shows the destinations without the roads between them.

The map is the missing middle: the global band network, from which the routes are read.

**Where the code goes.** `src/pytex/diffraction/kikuchi_map.py`, beside `kikuchi.py` and
`stereonets.py`. Same reasoning as the dynamical CBED entry below: `pytex/tem/` is scoped to
instrument operation, and band geometry on the crystal sphere is diffraction geometry. The
*routing* result is consumed by `pytex.tem` navigation, not defined by it.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 0 | Survey the three neighbouring surfaces; open this entry and the round-2 notebook rubric | done | (step 0) |
| 1 | Notebooks 02 and 06 to the round-2 rubric, proving the standard is reachable | done | (step 1) |
| 2 | `pytex.diffraction.kikuchi_map`: bands, zone axes, routing, `describe()`, JSON | done | (step 2) |
| 3 | `pytex.plotting.kikuchi_map`: the rendered atlas, cubic and hexagonal | done | (step 3) |
| 4 | Theory note, worked example, symbol registry, docs index, parity matrix | done | (this commit) |
| 5 | Tutorial 30 | done | (this commit) |
| 6 | Remaining round-2 notebooks, per `docs/development/notebook_improvement_progress.md` | pending | |

### Steps 2-5 outcome (2026-08-11)

`src/pytex/diffraction/kikuchi_map.py` (~1200 lines) and
`src/pytex/plotting/kikuchi_map.py`, with 42 tests across
`tests/unit/test_kikuchi_map.py` and `tests/unit/test_kikuchi_map_plotting.py`, the theory note
`docs/tex/algorithms/stereographic_kikuchi_maps.tex`, the worked example
`diffraction-kikuchi-map-zone-axis-tilt-angles`, four new registry symbols and four new terms, a
parity-matrix row claiming no MTEX parity with the reason, and tutorial 30.

**Zone axes are integer cross products, and that is not a micro-optimization.** The first
implementation crossed the *Cartesian* band normals, then searched a 7^3 index grid for the nearest
low-index direction per crossing, then deduplicated by an O(n^2) sweep over float vectors with an
angular tolerance. On a phase carrying an atomic basis that is tolerable, because the intensity
filter keeps the band count near 70. On a phase carrying **only a lattice** nothing can be filtered
on intensity, the band count goes to several hundred, and the sweep hangs — the worked example, which
builds its phase inline and therefore has no unit cell, never returned. Since the direct and
reciprocal bases are dual, the zone axis of two planes is the integer cross product of their index
triples: exact, deduplicated by a set lookup on tuples, no grid search. The same map now builds in
0.31 s instead of not at all.

**Three sign traps, all the same trap.** A zone axis is a line, not a direction, and the two senses
are the same axis. (a) Route legs reported each endpoint in the map's canonical sense, so consecutive
legs did not join and the drawn arc ended on the opposite side of the disc from the target marker.
(b) An equatorial axis has z of order 1e-17 of either sign, and a one-hemisphere projection folds on
the sign of z, so which side it landed on was decided by round-off. (c) Negating an exact zero gives
-0.0, whose sign bit is set, so the fold fired on it. Fixed by orienting the route chain explicitly,
snapping numerically-zero z to +0.0 at construction, and adding zero after negation. Tutorial 30
demonstrates the underlying fact as a failure mode, because any code comparing zone axes by index
equality will eventually hit it.

**A bare lattice degrades rather than raising.** `electron_structure_factor_angstrom` refuses a phase
without atom positions, correctly. The map now falls back to uniform intensities and records
`has_intensity_model=False`, which `describe()` states, because every geometric quantity — traces,
widths, zone axes, routes — is fully determined without a structure factor and that is most of what a
map is for.

### Step 0 outcome (2026-08-10)

Two things found and one fixed. Tutorial 29 had been committed **with its outputs and per-cell
execution timings**, leaving `tests/unit/test_notebooks.py` red on `main` and the file at
1.1 MB; cleared to 55 KB and pushed as the first commit of this task. And the round-1 notebook
overhaul, recorded as complete, is complete only against its own much lower bar — the explicit
round-2 rubric now lives in `docs/development/notebook_improvement_progress.md` and measures
the round-1 notebooks as short by roughly a factor of three in every dimension.

## Dynamical CBED: Many-Beam Coupling, Absorption, HOLZ Lines, Diffraction Groups — COMPLETE (2026-08-10)

**Objective.** Close the four gaps that the CBED step of the previous task explicitly listed as
not implemented, and that its `describe()` still advertises as missing:

1. **Many-beam coupling.** Every disc is currently its own two-beam calculation, so the discs of
   one pattern are not mutually consistent. Replace that with a Bloch-wave solution of the full
   coupled system.
2. **Absorption.** Without an imaginary part of the crystal potential the fringes never decay and
   the Borrmann (anomalous absorption) asymmetry — the thing that makes a real bright-field disc
   look the way it does — cannot appear.
3. **HOLZ lines.** Only the ring *radii* are given; the sharp deficiency lines inside the
   bright-field disc, which are the lattice-parameter metrology instrument, are absent.
4. **Diffraction-group symmetry determination.** CBED's most celebrated capability: the point
   group *including the presence or absence of a centre of symmetry*, which Friedel's law hides
   from kinematic SAED.

These are one capability, not four. The chain is: many-beam coupling makes the pattern a single
mutually consistent object; absorption makes it a physically realizable one; HOLZ beams in the
same beam set break the projection (ZOLZ-only) symmetry that would otherwise make every pattern
look centrosymmetric; and the symmetry that survives that breaking *is* the diffraction group.

**Where the code goes, and why not in `pytex.tem`.** The request said "in the TEM module". The
work lands in `pytex/diffraction/` alongside `cbed.py`, because `pytex/tem/` is scoped by
`docs/architecture/tem_tilt_navigation_foundation.md` and its own package docstring to
*instrument operation* — stage models, tilt solving, reachability — and explicitly not to new
crystallography. Dynamical electron scattering is diffraction physics and belongs with the
diffraction engine it extends. Recorded here rather than decided silently.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 0 | Survey `cbed.py`, the symmetry surface, and the docs hooks; open this entry | done | (step 1) |
| 1 | `pytex.diffraction.dynamical`: Bloch waves + absorption, tests, theory note | done | e671b8f |
| 2 | `pytex.diffraction.holz`: HOLZ line loci, chords, metrology sensitivity | done | 660ff94 |
| 3 | `pytex.diffraction.diffraction_groups`: the 31 groups, forward and inverse | done | 660ff94 |
| 4 | Wire all three into `CBEDPattern`; retire the "not implemented" limits text | done | f12addc |
| 5 | Notebook 29, docs index, symbol registry, worked examples, parity matrix | done | (this commit) |

### Step 1 outcome (2026-08-09)

`src/pytex/diffraction/dynamical.py` (≈900 lines), 20 tests in
`tests/unit/test_dynamical.py`, and the theory note
`docs/tex/algorithms/dynamical_cbed_and_symmetry_determination.tex`.

**The scale is inherited, not re-asserted.** The off-diagonal coupling is
`nu_g = lambda F_g / (pi V_c cos theta_g)`, chosen so that `|nu_g| = 1/xi_g` for the
extinction distance already validated against Williams & Carter Table 23.1. The two-beam limit
of the many-beam solver then reproduces `two_beam_rocking_curve` to 2e-15, which pins the
diagonal convention (`2 s_g`), the off-diagonal scale and the `i pi` in the propagator
simultaneously. A many-beam module that introduced its own absolute scale would have needed a
second validation; this one does not.

**Three exact properties are asserted rather than assumed:** the two-beam limit above;
unitarity (`sum I_g = 1` to 1e-12 at every tilt and thickness with absorption off, which is
what catches the classic error of obtaining the Bloch-wave excitation amplitudes by
projection instead of by solving — the eigenvectors are not orthogonal); and the fact that
normal absorption factors exactly out of the matrix exponential as `exp(-2 pi t / xi'_0)`, so
the phenomenological `mean_ratio` provably cannot contaminate any statement about shape,
position or symmetry. `normal_absorption_factor` exposes that scalar so a caller can divide
it out.

**Absorption is structural, its magnitude is not, and the docstring says which is which.**
The imaginary optical potential enters the structure matrix, so anomalous absorption emerges
from the eigenvector structure rather than being applied to the output. The test is the
Hashimoto-Howie-Whelan theorem: with absorption the bright-field rocking curve becomes
asymmetric (>10 percent) while the dark-field one stays symmetric to 1e-10. The ratios
themselves are the customary 1/10 of Hirsch et al.; Bird-King absorptive form factors are not
implemented and `AbsorptionModel.describe()` says so. `reflection_ratio > mean_ratio` is
rejected because it would give a Bloch wave that gains intensity with depth.

**The centrosymmetry mechanism was derived, and it changed the design.** `A` is Hermitian for
any real potential but *symmetric* only when every included `nu_g` is real. Relabelling
`g -> -g`, `theta -> -theta` turns `A` into `A^T` (exactly, for a ZOLZ set), so the propagator
becomes `M^T` and `I_g(theta) = I_-g(-theta)` holds **iff** `M` is symmetric **iff** the
sampled structure is centrosymmetric. Friedel's law is therefore a theorem about the
propagator, not a kinematic accident.

What that derivation exposed is the thing the module now warns about loudly: a ZOLZ-only beam
set samples the *projected* potential, and for zincblende down `[111]` every ZOLZ coefficient
is real. A projection calculation reports Friedel's law to 1e-14 for GaAs and cannot see the
polarity at all. Admitting the first-order Laue zone breaks it by 15 percent absolute
(26 percent relative). The control that makes this a measurement rather than a coincidence is
a rocksalt structure on the *same* lattice with the *same* two species, differing only by
where the second sublattice sits: its violation stays below 1e-3. Three orders of magnitude
separate them. `BeamSet.holz_mask` and `BeamSet.describe()` exist so that a caller cannot
draw a symmetry conclusion from a projection calculation without being told.

**Consequently the beam selection had to change.** A HOLZ reflection is far from Bragg on axis
and exactly at Bragg somewhere inside the bright-field disc, so selecting on the zero-tilt
excitation error discards every HOLZ beam and with it the whole mechanism. Selection is on
`min |s_g|` over the illumination cone, which is available in closed form because `s_g` is
affine in the tilt.

**Not implemented, and stated in the module docstring:** Bethe perturbation of weak beams
(so a full HOLZ ring costs `O(m n^3)` in earnest — the working economy is a tighter
excitation window, not coarser tilt sampling), wedge/bent/strained specimens, and probe
aberrations.

**Side effect.** Nine new public names moved the class-model atlas counts from 250/233 to
253/236, so `docs/figures/class_model_*.svg` were regenerated and the atlas page's prose
counts updated.

### Step 2 outcome (2026-08-10)

`src/pytex/diffraction/holz.py`, 16 tests in `tests/unit/test_holz.py`.

**The geometry is exact and is checked against the other module.** Because `s_g` is affine in the
incident tilt, `s_g = 0` is a straight line, so line positions are closed-form rather than sampled.
The test takes points on a line and asks `pytex.diffraction.dynamical` — which derived the
excitation error for a different purpose — what `s_g` is there: zero to 1e-15. Two modules agreeing
to machine precision is worth more than either agreeing with a stored number.

**The metrology trap is a theorem, not a caveat.** Scaling the lattice by `1+eps` and scaling the
wavelength by the same factor return every line to its original position **exactly**, at every
reflection simultaneously (asserted to 1e-16 for strains from 1e-4 to 2e-2). So HOLZ line positions
cannot separate a strained lattice from a mis-set accelerating voltage, which is why quantitative
HOLZ work calibrates the voltage on a standard first. `offset_at(lattice_strain=..., 
wavelength_angstrom=...)` takes both arguments precisely so the degeneracy cannot be papered over.

**The numbers explain the practice.** For Ni [001] at 200 kV in a 1000 A foil, the best *single*
line resolves a strain of 3.6e-3 — far short of the 1e-4 the technique is known for. The best
*intersection* of two near-parallel lines resolves 6.3e-5, because a crossing moves as
`1/sin(phi)` times faster than its lines. That is why HOLZ measurements read intersections, and the
module reports the amplification rather than leaving the reader to wonder where the sensitivity
comes from. Line half-width `1/(t |g_perp|)` falls as `1/t`, so the resolvable strain does too:
HOLZ metrology wants a *thick* foil, the opposite of the usual instinct.

**Deliberately geometry only.** No intensities: line positions are exact, line contrast needs the
coupled dynamical solution and inherits all of its approximations. `holz_line_pattern` accordingly
does not require a unit cell, and a test asserts that — atoms decide whether a line is visible, not
where it is.

### Step 3 outcome (2026-08-10)

`src/pytex/diffraction/diffraction_groups.py`, 33 tests in
`tests/unit/test_diffraction_groups.py`.

**The 31 diffraction groups are derived, not transcribed.** Each crystal operator `S` is classified
by its action on the beam direction — fixes it, reverses it, or neither — and the first two classes
contribute their transverse restriction `S|_perp`, tagged with the reciprocity flag when `S`
reverses the beam. That map is a homomorphism onto a subgroup of (plane point group) x Z2, and
scanning all 32 point groups over their characteristic directions yields exactly Buxton's 31, whose
membership the test compares symbol by symbol. The point-group-to-diffraction-group table
(Buxton Table 2) is likewise computed by inversion rather than copied.

**Two of the observables are derived too.** Whole-pattern symmetry is the untagged subgroup, because
a tagged element needs reciprocity, which maps a point in one disc to an incident direction outside
the illumination cone. Bright-field symmetry is `phi(D)` with `phi(T, tagged) = -T`, because inside
the direct disc the reciprocity displacement is proportional to `g_perp` and therefore vanishes,
leaving only reciprocity's own inversion of the incident direction. Every canonical entry checks
out: `m-3m [001] -> 4mm1_R`, `-43m [001] -> 4_Rmm_R (BF 4mm, WP 2mm)`, `432 [001] -> 4m_Rm_R`,
`-6m2 [001] -> 3m1_R` with a **six**-fold bright-field disc over a `3m` whole pattern, `m-3m [111]
-> 6_Rmm_R` with `3m` whole pattern (not the 6mm a kinematic pattern appears to show).

**Centrosymmetry is an exact correspondence.** `2_R` requires an operator acting as `-1` on the beam
direction and as `-1` on the transverse plane, which is the inversion and nothing else. So `2_R` is
in the diffraction group at *every* beam direction of a centrosymmetric crystal and at none of an
acentric one — asserted over all 32 point groups at every characteristic direction. Supplying only
the `+-g` observation therefore splits the 32 point groups into exactly 21 and 11, which is the
arithmetic of the whole technique and is asserted as such.

**And `2_R` is invisible in BF and WP** (`phi(2, tagged) = -2 = 1`), so no amount of disc symmetry
decides the centre. That is why `SymmetryObservations` carries `friedel_pair_two_fold` as a separate
field and why leaving it unknown leaves the verdict `None` by construction rather than by accident.

**Not implemented, and `describe()` says so:** Buxton's dark-field and `+-g` observations for
reflections lying on symmetry lines, recorded at their own Bragg condition. They would narrow cases
such as `4_Rmm_R -> {-42m, -43m}` that the three implemented observations leave open; the report
recommends a second zone axis instead, and names the tool that finds it.

**Side effect.** Sixteen more public names moved the class-model atlas counts from 253/236 to
260/243; figures regenerated and the atlas prose updated.

### Step 4 outcome (2026-08-10)

`ConvergentBeamConfig` gains `method`, `absorption`, `laue_zones`, `holz_max_index` and
`holz_g_max_inv_angstrom`; `CBEDPattern` gains `beam_set`, `holz_lines`,
`predicted_diffraction_group()`, `symmetry_observations()` and `determine_point_group()`. Nine
new tests in `tests/unit/test_cbed.py`. The `describe()` text that advertised all four features
as unimplemented is gone.

**The import cycle was removed structurally, not with a lazy import.**
`electron_structure_factor_angstrom` moved from `cbed.py` to `scattering.py`, which is its
proper home next to `electron_scattering_factors` and which nothing in the diffraction package
depends on. `cbed` now imports `dynamical`, `holz` and `diffraction_groups` at module level and
none of them imports `cbed`. The public import path is preserved by re-export, and
`tests/unit/test_cbed.py` still imports it from `cbed` unchanged.

**`method="two-beam"` stays the default.** It is the model `thickness_from_fringe_minima`
inverts, and switching the default would silently change every existing result. `"bloch"` is
the coupled path; the config refuses combinations that would be silently ignored — absorption
without `bloch` (the closed form has no absorptive term), HOLZ zones without `bloch` (each disc
is independent there, so a HOLZ beam could not change it), and `laue_zones` without `0` (there
would be no discs to draw).

**The headline result, end to end.** Zincblende GaAs and a rocksalt structure on the *same*
lattice with the *same* two species, differing only by where the second sublattice sits:

| | predicted | measured BF | measured WP | determination |
| --- | --- | --- | --- | --- |
| GaAs `[001]`, ZOLZ+FOLZ | `4_Rmm_R` | `4mm` | `2mm` | `{-42m, -43m}`, **not centrosymmetric** |
| control `[001]`, ZOLZ+FOLZ | `4mm1_R` | `4mm` | `4mm` | includes `m-3m`, centre not excluded |
| GaAs `[001]`, **ZOLZ only** | `4_Rmm_R` | `4mm` | `4mm` | looks exactly like the control |

The residuals are `0.00` against `0.32` — no tolerance judgement is involved. The third row is
the point of the whole exercise: same crystal, same code, one flag, and the missing centre
becomes invisible, because the projected potential of zincblende down `[001]` *is*
centrosymmetric. `symmetry_observations()` therefore refuses a ZOLZ-only pattern unless asked
twice, and refuses a two-beam pattern outright.

**Three measurement bugs found and fixed, each of which produced a plausible wrong answer:**

1. **The surviving operations must be closed into a group before naming.** They are tested as
   generators — one rotation per order — so `{1, R2, R3, R6}` is four matrices and was being
   counted as a four-fold axis where the crystal has a six-fold.
2. **A worst-case criterion is unusable.** HOLZ lines are narrower than the tilt sampling can
   resolve, so resampling at a rotation that is not grid-aligned produces large errors along a
   few thin loci while the map as a whole is symmetric. The residual is now the mean absolute
   deviation, which weights those loci by their area.
3. **Per-disc normalization needs a floor.** A systematically absent reflection has an
   identically zero disc whose floating-point noise has a mean absolute deviation of order
   1e-30; dividing by it turned rounding error into a catastrophic symmetry violation. Silicon
   down `[001]` has four such discs, the absent `{200}`, and they destroyed the four-fold.
   But normalizing by the *brightest* disc instead is equally wrong — GaAs breaks its four-fold
   in the near-forbidden `{200}` discs, whose contrast is half a percent of the strongest, and
   that normalization hid it. The floor is 1e-6 of the strongest disc's contrast.

**The `+-g` observation is deliberately not measured, and this was the hardest finding.**
Buxton's `2_R` compares the `+g` and `-g` *dark-field* discs, each recorded with its own
reflection at the Bragg condition — two exposures at different specimen tilts, related by
reciprocity. It is not a two-fold rotation of a single zone-axis pattern. Treating it as one
gives a test that fails, and the derivation says why:
`s_-g(-theta) - s_g(theta) = -2 g_z`, which vanishes only in the zeroth Laue zone. The
numerical check settles it: the residual *grows* with the beam set for centric and acentric
structures alike (Si 0.06 -> 0.37, GaAs 0.13 -> 0.55 as the excitation window widens), so it is
physics and not truncation. `symmetry_observations` therefore leaves `friedel_pair_two_fold` to
the caller rather than reporting a number that would sometimes be wrong. The determination does
not need it: at `[001]` the bright-field and whole-pattern symmetries settle the centre outright.

### Step 5 outcome (2026-08-10)

`docs/site/tutorials/notebooks/29_dynamical_cbed_and_point_groups.ipynb` (54 cells, 24 code,
runs in 23 s), six worked examples in `worked_examples/examples/dynamical_cbed.py`, four new
rows in `docs/testing/diffraction_validation_matrix.md`, and a CHANGELOG entry.

**The control was changed while writing the notebook, and that was a real finding.** The
rocksalt structure used in the unit tests — same lattice, same two species, offset by 1/2
instead of 1/4 — is an excellent control for the whole-pattern test and *useless* for the
`+-g` one, because gallium and arsenic are neighbours in the periodic table and its
higher-order reflections go as `f_Ga - f_As`, which is nearly zero. It would have shown a
residual near zero and made the naive `+-g` test look reliable. Diamond silicon is the right
control: identical structure type to zincblende, identical site geometry, differing only in
whether the two sublattices carry the same element. Notebook section 5.3 now makes that
point explicitly — a control that is perfect for one measurement can be uninformative for
another.

**Every number in the notebook is computed, including the ones in the prose.** A sentence
quoting the symmetry residuals as "0.00 against 0.32" was replaced by a cell that reads the
symmetry back at four tolerances spanning two decades and shows the answer is flat, which is
the claim that actually matters and is checked rather than transcribed.

**Worked-example provenance.** All six expected values are analytic identities or published
counts, never a copied program output: the two-beam closed form (exact, tolerance 1e-12),
unitarity (exact), the strain/wavelength cancellation (exact, 1e-15), Buxton's 31 diffraction
groups, the 21/11 split of the point groups from the International Tables, and the published
diffraction-group assignments for `-43m` and `m-3m` down a four-fold axis.

## Dynamical CBED — COMPLETE (2026-08-10)

All five steps are landed and pushed. The four gaps the previous task listed as unimplemented
are closed, and the limits that remain are named in `describe()` on every object that has them:
Bethe perturbation of weak beams, computed absorptive form factors, Buxton's special dark-field
and `+-g` observations, and specimen realism.

## Repository Content Rule And PDF History Purge — COMPLETE (2026-08-09)

**Objective.** Make "the repository holds sources and canonical assets only" a cardinal rule that
is enforced rather than remembered, and resolve the tracked-PDF exception it exposed.

**The rule.** `AGENTS.md` gains a second cardinal rule and a `Repository content` section stating
one test: if a command in this repository regenerates a file, it is committed only when
documentation, a test, a manifest, or a pinned baseline names it. Both lists — what is committed
(sources, canonical `docs/figures/` SVGs including the generated ones, generated galleries and
fixtures tests load) and what never is (build output, caches, notebook outputs, inspection renders,
logs, scratch, reference PDFs) — are explicit, because the previous phrasing left "generated" and
"canonical" to judgement. Cross-referenced from `docs/standards/engineering_governance.md` and
`.gitignore`.

**The enforcement.** `scripts/check_repo_integrity.py::_check_repository_content` fails on any
tracked path in the excluded categories; it runs in the base lane via
`tests/unit/test_repo_integrity.py`, which tests one path per category and the canonical generated
assets it must *not* reject. Prose alone would not have held.

**What it found immediately.** Eight reference PDFs, 107 MB of third-party textbooks, tracked since
before the no-PDFs rule existed. On the maintainer's decision they were untracked and purged from
history with `git filter-repo`, and `main` was force-pushed. Result: `.git` 131 MB -> 35 MB
(pack 122.7 MiB -> 34.2 MiB). Local copies were untouched; `references/reference_index.md` now
states that no PDF is tracked, which files are expected locally, and that a filename is not a
citation. There is no grandfather list left — the check rejects every tracked PDF.

**Backups taken before the rewrite** (outside the repo, safe to delete once the rewrite is trusted):
`../pytex-backup-2026-08-09/pytex-full-history.bundle` (124 MB, the complete pre-rewrite history)
and `../pytex-backup-2026-08-09/references/` (all nine PDFs).

**Consequence for other clones.** Every commit hash changed. The maintainer confirmed on
2026-08-09 that no other clone or fork existed, so nothing had to be re-cloned and no follow-up is
outstanding. Any clone taken before that date would need re-cloning.

## Completed Task: Orientation Representations, TEM Round-Trip Indexing, And CBED — COMPLETE (2026-08-09)

**Objective.** Three things that belong together because they share one conversion spine:

1. **A core orientation-representation surface.** PyTex can already build a `Rotation` or an
   `Orientation` from a matrix, a quaternion, axis-angle, Rodrigues, three Euler conventions, and
   `(hkl)[uvw]`, and can emit most of those again. What is missing is (a) the two representations
   with no constructor at all — **homochoric** and **cubochoric** — and (b) the operation users
   actually ask for: *give me one orientation and print every form of it at once*, in readable,
   standard notation. That becomes `pytex.core.representations`, with a vectorized batch path.
2. **A TEM zone-axis round trip.** Generate a pattern with `generate_saed_pattern`, feed it back to
   `solve_saed_pattern`, and show the original zone axis and orientation come back. This is the
   honest self-consistency proof for the indexing chain, and it currently exists nowhere as a
   documented artifact.
3. **CBED.** Convergent-beam diffraction is absent from the library. It is the natural companion to
   SAED — same geometry, one extra parameter (the convergence semi-angle) — and it is what gives
   thickness measurement and point-group determination.

Each gets a tutorial notebook held to the standard of notebook 25 (pole-figure arithmetic): live
computation, rendered mathematics, algorithm boxes, and the failure modes shown rather than
described. The two TEM notebooks each cover **Ni (FCC)** and **Zr (HCP)**, the fixtures already in
`pytex.core.fixtures`.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 0 | Survey the existing conversion surface; open this entry | done | (this commit) |
| 1 | `pytex.core.representations` + tests + theory note | done | (this commit) |
| 2 | Notebook 26: rotation and orientation representations | done | (this commit) |
| 3 | CBED module + tests + algorithm note | done | (this commit) |
| 4 | Notebook 27: TEM zone-axis indexing round trip (Ni, Zr) | done | (this commit) |
| 5 | Notebook 28: CBED analysis (Ni, Zr) | done | (this commit) |
| 6 | Docs index, symbol registry, worked examples, parity matrix | done | (this commit) |

### Step 1 outcome (2026-08-09)

`src/pytex/core/representations.py` (≈1150 lines) plus 39 tests in
`tests/unit/test_orientation_representations.py`, and the theory note
`docs/tex/theory/orientation_representations.tex`.

**The cube-to-ball map was derived, not transcribed.** The Roşca–Morawiec–De Graef map is
usually quoted as a block of constants, and a mis-remembered constant produces a map that still
looks plausible. It is instead derived here from the two conditions that determine it —
sub-cubes map onto spheres of equal volume, giving `r = z (6/π)^(1/3)`; and each cube face maps
to its own spherical *sector*, not a cap, because six caps of solid angle 2π/3 overlap and
cannot tile. Those force an area-preserving square-to-sector map, which factors into a planar
wedge and a Lambert lift, and fix the wedge prefactor analytically at `2^(1/4) √(6/π)` —
confirmed by the fact that the same value then satisfies the sector-boundary condition at every
azimuth, not just the one it was solved at.

**The tests assert the defining properties rather than stored numbers.** Jacobian determinant
1 everywhere (finite differences, off the branch kinks), the cube corner landing exactly on the
ball surface, every face of every sub-cube landing on one sphere, continuity across the pyramid
boundary and the face diagonal, and — the sharpest one — a uniform cube sample reproducing the
analytic mean rotation angle of SO(3), `π/2 + 2/π = 126.4756 deg`. A wrong constant fails all of
them; a table of copied outputs would have failed none.

**Also landed:** vectorized Euler conversions that match the per-object path exactly (the
existing `RotationSet.as_euler_set` loops in Python); `convert_orientations`, routing all ten
representations through quaternions so there are ten conversions rather than ninety;
`ideal_orientation_indices`, the inverse of `Orientation.from_miller`, reporting the nearest
`(hkl)[uvw]` **with its residual angle** and the four-index form for hexagonal phases; and the
`OrientationRepresentations` / `OrientationRepresentationSet` reports with `describe()` and
`to_json_dict()`.

**Side effect.** The 21 new public names moved the class-model atlas counts from 242/226 to
246/229, so `docs/figures/class_model_*.svg` were regenerated and the atlas page's prose count
updated — the test that forbids hand-transcribed counts caught it, as designed.

### Step 2 outcome (2026-08-09)

`docs/site/tutorials/notebooks/26_orientation_representations.ipynb`: 54 cells, 25 of them
code, every number computed live. Its spine is the *measure problem* — section 6 puts three
samplers against the analytic law `(1 - cos w)/pi` and shows uniform Euler angles producing a
visibly non-uniform set of orientations, which is the fact that makes the equal-volume charts
worth having rather than exotic. Sections 8.1-8.3 verify the cube-to-ball map in the notebook
itself (nested spheres to nine decimals, unit Jacobian over 400 probes, and the grid picture),
so a reader is shown the evidence rather than asked to trust it.

**Two real inaccuracies found while writing it**, both fixed in this commit:

1. **The Rodrigues-Frank docstrings claimed the form "stays finite" at `omega = pi`.** It does
   not: the magnitude is `tan(omega/2)`, and PyTex stores `inf` there. The claim was wrong but
   the *behaviour* is right and better than the claim — the magnitude is a projective
   coordinate, so 180 degrees is exactly representable and exactly invertible (axis preserved),
   whereas the plain 3-vector overflows to ~1e16 and loses its axis in the product. Five
   docstrings across `orientation.py` and `batches.py` now say that instead.
2. **`OrientationRepresentations` documented a non-negative quaternion scalar part it did not
   enforce.** `Rotation` keeps whichever sign the arithmetic produced, so a batch row and a
   single report of the same rotation could differ by a global sign. Reporting surfaces now
   canonicalize through the new public `canonical_quaternions`, which also resolves the
   180-degree tie where both signs have `w = 0`.

Also corrected in the notebook draft: `Rotation.distance_to` returns **radians**, not degrees
(three cells were mislabelling it), and the S component was written `(123)[6 -3 4]`, which does
not satisfy the zone law — `(123)[6 3 -4]` does, and now round-trips exactly.

### Step 3 outcome (2026-08-09)

`src/pytex/diffraction/cbed.py`, 31 tests in `tests/unit/test_cbed.py`, and
`docs/tex/algorithms/convergent_beam_electron_diffraction.tex`.

**The hard part was the absolute scale, not the geometry.** Disc geometry and rocking curves
are easy to make *look* right; an extinction distance that is wrong by a constant produces
perfectly plausible fringes at the wrong spacing. Two things fix it and both are pinned by
tests:

1. **Mott-Bethe.** Electrons scatter from the potential, so `f_e = (Z - f_x)/(8 pi^2 a0 s^2)`.
   PyTex's X-ray table is stored as `Z - 41.78214 s^2 sum a_i exp(-b_i s^2)`, and 41.78214
   *is* `8 pi^2 a0` — so the inversion returns exactly `sum a_i exp(-b_i s^2)` and introduces
   no new constant. That identity is asserted directly, which is what pins the scale.
2. **The relativistic factor** `gamma = 1 + E/m0c^2` (1.39 at 200 kV). Omitting it lengthens
   every extinction distance by 39 percent, and nothing else in the pattern would look wrong.

Validated against Williams & Carter Table 23.1 for aluminium at 100 kV: `{111}` 555 vs 556,
`{200}` 664 vs 673, `{220}` 1063 vs 1057 — within 1.4 percent. Aluminium was chosen because
the fitted parametrization is most accurate for light elements; for Ni it runs ~11 percent
high, which is stated in the docstring rather than hidden, and is exactly why the thickness
fit *measures* the extinction distance instead of taking it from a table.

**The capability claim is a round trip**, tested through the public path: simulate a pattern
at a known thickness, read the fringe minima off a disc, and recover both the thickness and
the extinction distance from the Kelly `(s_n/n)^2` versus `1/n^2` line.

**Two things the tests taught, both now documented:** forcing the classic wrong assumption
`n = 1` on the innermost minimum does not merely bias the thickness — it tilts the fitted
line the wrong way and yields a *negative* `1/xi^2`, so the routine can and does refuse.
And a disc is centred at `s = -lambda g^2/2`, so exact Bragg lies inside it only when the
convergence angle exceeds the Bragg angle; below that the two fringe branches are unequal
and the richer one is the one to measure.

**Not implemented, and said so in `describe()`:** many-beam coupling (each disc is its own
two-beam calculation), absorption, HOLZ lines within the bright-field disc, and the
diffraction-group symmetry determination that would give the point group including its
centrosymmetry.

### Step 4 outcome (2026-08-09)

`docs/site/tutorials/notebooks/27_tem_pattern_indexing_round_trip.ipynb`: 31 cells, 12 code.
Simulate a zone-axis pattern, hand the bare spot positions to `solve_saed_pattern` as if
picked off a micrograph, and check what comes back — for four Ni zones and four Zr zones.

**The notebook states what a round trip does not prove.** It tests the forward and inverse
models for *mutual* consistency; a convention error shared by both would round-trip perfectly.
So it is an internal-consistency test that complements, and does not replace, the external
pymatgen baselines. Saying which is which was worth a paragraph.

**Three findings the run produced, each now a section rather than a footnote:**

1. **The recovered zone axis is a symmetry orbit, not a triple.** Zr `[10-10]` comes back as a
   different member of its 3-element orbit. Checking a round trip by string comparison would
   call that a failure; checking up to symmetry calls it correct, which it is.
2. **Ni `[112]` has a genuine residual ambiguity, and the others do not.** Every pattern is
   invariant under a half turn about the beam, because Friedel's law makes the spot set
   centrosymmetric. For `[001]`, `[011]` and `[111]` that half turn *is* a cubic symmetry
   operation, so the disorientation is zero. `<112>` is not a two-fold axis of m-3m, so the
   two returned solutions are physically different orientations with identical patterns — the
   true one is in the report but is not ranked first. This is what `pytex.tem.ambiguity`
   enumerates, and the notebook now connects the two.
3. **Comparing orientations without symmetry reduction looks like a bug.** The raw angle
   between the true and recovered rotations is routinely 90-180 degrees while the
   disorientation is zero. Shown as a column rather than described.

**The inconclusive example had to be made honest.** Two spots at the default 3 percent length
tolerance still identify nickel; the first draft claimed otherwise. At 5 percent — realistic
for hand-picked spots — a zirconium zone explains them equally well and `is_conclusive` flips.
The notebook shows both tolerances and draws the real lesson: set the tolerance to what the
picking achieves, because a tighter one buys confidence the data does not support.

### Step 5 outcome (2026-08-09)

`docs/site/tutorials/notebooks/28_convergent_beam_diffraction.ipynb`: 29 cells, 13 code, every
number computed live for both Ni and Zr.

**The geometric chain it makes explicit, and which is not obvious from the literature:**
discs touch when `alpha = theta_B` of the innermost reflection (the Kossel-Moellenstedt
threshold), and exact Bragg lies inside a disc only when `alpha > theta_B`. These are the
*same* condition, so **at a zone axis in the KM regime a disc necessarily shows one wing of
the rocking curve**, never the symmetric two-wing curve textbooks draw. That is not a
simulation artefact — real thickness measurement is done at a two-beam condition reached by
tilting off the zone axis — and the Kelly fit does not care, because it needs the minimum
positions and their orders, not both wings. The notebook derives this rather than working
around it.

**Verified live in the notebook:** aluminium extinction distances against Williams & Carter
(within 1.4 percent); thickness round trips for Ni (200) at 1500 and 2000 A and Zr (10-10) at
2000 and 3000 A, all in the KM regime and all recovering `t` to better than 0.1 percent and
`xi_g` to better than 0.1 percent; and HOLZ ring radii inverted to recover Ni `a` down [001],
`a*sqrt(3)` down [111], Zr `c` down [0001] and Zr `a` down [11-20] — the Zr `c` case being the
point, since the [0001] spot pattern cannot see `c` at all.

**Both failure modes shown running, not described:** the Zr case at 10 mrad falls into the
Kossel regime, where R^2 drops to 0.993 and both fitted numbers are wrong (the R^2 is the
warning); and forcing `first_order=1` on a clean case raises, because the wrong assignment
implies a negative `1/xi^2`.

### Step 6 outcome (2026-08-09) — TASK COMPLETE

Four executable worked examples in two new groups, `orientation-representations` and
`convergent-beam-diffraction`, each with independent provenance rather than a copied program
output: the equal-volume identity (ball and cube both enclosing `pi^2`, and the cube corner
landing on the ball surface), the `from_miller` / `ideal_orientation_indices` inversion, the
aluminium extinction distances against Williams & Carter Table 23.1, and the Kelly inversion
of analytically generated fringe minima. Gallery regenerated; `test_worked_examples.py` green.

Registry, matrices, and foundations updated: nine new symbols registered before use
(`rho`, `rho_F`, `h`, `c`, `R_1`, `a_p`, `s_g`, `xi_g`, `f_e`, `alpha`, `t`, `H`); seven CBED
rows added to the diffraction validation matrix including an explicit `planned` row for the
dynamical breadth that is *not* implemented; two MTEX rows added that state plainly that MTEX
has **no** equal-volume chart, so no parity is claimed and the comparable surface is
ORIX/EMsoft; and both the diffraction and orientation foundation documents corrected, since
the diffraction one still claimed "no dynamical intensity simulation".

**The Sphinx build was broken before this task started, and is now fixed.** The notebook policy
rests on the site executing every notebook (`nb_execution_mode = "cache"`,
`nb_execution_raise_on_error`), so "a notebook that no longer runs fails the docs build" is the
guarantee that replaces stored outputs. But `docs/site/tutorials/notebooks/15_*.ipynb` read
three manifests by repository-relative path, and myst-nb executes a notebook with *its own
directory* as the working directory — so the build had been failing since the policy was
switched on, while `tests/unit/test_notebooks.py` passed because it `chdir`s to the repo root.
Notebook 15 now resolves the repository root by walking up to `pyproject.toml`, which works
under all three execution contexts. With that fixed the full site builds, which is the first
confirmation that notebooks 26, 27 and 28 execute under the documentation build rather than
only under a manual runner.

Two smaller defects the build surfaced and this commit fixes: the CBED worked example linked to
a non-existent concept page (`concepts/diffraction_geometry`; the page is
`concepts/diffraction_foundation`), and notebook 27 jumped from an H1 to an H3.

### Task status

**COMPLETE.** All six steps landed and pushed. Base lane green (`ruff`, `mypy`, `pytest`), the
repository-integrity check passes, and `python -m sphinx -b html docs/site docs/_build/html`
completes with only the pre-existing duplicate-object-description warnings from
`pytex.adapters.*`, which predate this task.

### Next task

None claimed. Candidates in the capability review's recommended order: roadmap reconciliation
plus `windows-latest` in CI, then the defocus model and ghost correction. A natural follow-on
to this task specifically would be many-beam (Bloch-wave) CBED intensities, which is the
prerequisite for the diffraction-group symmetry determination that the current module
explicitly does not attempt.

## Completed Task: Class & Object Model Atlas — COMPLETE (2026-08-09)

**Objective.** Give the Sphinx site a Class & Object Model Atlas: an overview of the library
architecture, an honest class-hierarchy view, and UML-style object-model diagrams per domain
(core, texture, EBSD, diffraction, TEM), generated from the source so they cannot drift.

**Key finding from the first inspection pass (drives the whole design).** PyTex has 241 public
classes, 225 of them dataclasses, and only **6 internal inheritance edges** (`TiltEnvelope` ->
4 envelopes; `ElasticTensor` -> stiffness/compliance) plus 10 `StrEnum` vocabularies and 2
`Protocol`s. A Doxygen-style inheritance atlas would therefore be nearly empty and would
misrepresent the design. The rich, real structure is **composition**: 390 typed dataclass-field
references between public classes. The atlas leads with the object model and states the
composition-first fact explicitly, with the small hierarchy shown in full rather than padded.

**Tooling constraint.** No Graphviz `dot` binary and no `pylint`/`pyreverse` in this environment,
so `sphinx.ext.inheritance_diagram` and pyreverse are both unavailable — and adding a
system-binary dependency to the docs build would be a regression. The atlas is therefore rendered
by PyTex's own SVG stack (`pytex.plotting.svg_primitives`, canonical style tokens, Arial advance
metrics), the same path the reference-frame and algorithm figures already use.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 0 | Inspect codebase + docs architecture, choose approach | done | (this entry) |
| 1 | Introspection module: build the class model from the source | done | (this commit) |
| 2 | `pytex.plotting.class_diagrams` layered UML renderer + tests | done | (this commit) |
| 3 | Generator script + canonical SVG assets in `docs/figures/` | done | (this commit) |
| 4 | Sphinx page, cross-links from architecture/API/library-structure | done | (this commit) |
| 5 | Build the site, inspect the rendered diagrams, iterate | done | (this commit) |

### Next actions

Step 1: `scripts/class_model.py` — walk `pytex`, resolve dataclass field type hints, emit nodes
(name, module, stereotype, key fields) and typed relations (composition / association /
inheritance), with per-domain view selection.

## Current Task: Pole-Figure Arithmetic (PFA) — COMPLETE (2026-08-08)

**Objective.** Make pole-figure arithmetic possible. It is today structurally blocked, not merely
unwritten: `PoleFigure` holds *scattered* specimen directions, so two figures generally share no
support and cannot be combined at all; there is no resampling onto a common grid, no spherical
interpolation, and the XRDML/LaboTex adapters normalize by `max` or `sum` rather than to multiples
of a random distribution, so magnitudes are not physically comparable between figures. Arithmetic
dunders cannot be added until resampling and m.r.d. normalization exist.

**Secondary outcome.** ODF inversion gains the residual-figure QC product it currently lacks:
measured-minus-recalculated pole figures, plotted.

### Dependency order (this is why the work is sequenced, not parallel)

```
sampling semantics  ->  on_grid (resampling)  ->  normalize_to_mrd  ->  arithmetic  ->  residual QC product
```

Nothing downstream is meaningful without the step before it.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 0 | Cardinal rule (ledger + commit/push cadence) into `AGENTS.md` | done | `ee9591c` |
| 1 | `PoleFigure.sampling` semantics field + contract round-trip | done | `6c925b2` |
| 2 | `PoleFigure.on_grid` spherical resampling | done | `82dc073` |
| 3 | `PoleFigure.normalize_to_mrd` + `raster_solid_angle_weights` + adapter `mrd` mode | done | `4b1bdec` |
| 4 | Arithmetic, `PoleFigureDifference`, `symmetrize`/`rotate`/`restrict_polar_range` | done | `ce3146f` |
| 5 | Residual QC product + difference plotting | done | `920e6c4` |
| 6 | Exports, worked example, docs, CHANGELOG, review-scorecard update | done | this commit |

### Design decisions worth not re-deriving

1. **`PoleFigure` carries two different scientific meanings today, and arithmetic forces the
   distinction into the open.** `from_orientations` produces a *cloud of poles* whose `intensities`
   are per-pole weights; the XRDML/LaboTex adapters and `HarmonicODF.reconstruct_pole_figure`
   produce a *sampled density field* whose `intensities` are pole densities at given directions.
   Resampling a cloud means kernel **density estimation** (a weighted sum); resampling a field
   means kernel **interpolation** (a weighted mean, Nadaraya-Watson). Applying the wrong estimator
   to a regular latitude-longitude grid biases the result towards the poles, because such a grid
   oversamples there. So the reading is recorded on the object as `sampling`, in the same spirit as
   the existing `antipodal` and `includes_symmetry_family` flags, rather than assumed at the call
   site.
2. **Subtraction cannot return a `PoleFigure`.** The class enforces non-negative intensities, which
   is correct — a pole density is non-negative. A difference is signed. So `__sub__` returns a
   distinct `PoleFigureDifference` type carrying the signed field plus its residual statistics;
   `__add__`, `__mul__` and `__truediv__` return `PoleFigure`. The asymmetry is deliberate and is
   what keeps the non-negativity invariant honest.
3. **m.r.d. needs solid-angle weights, which scattered directions do not carry.** Where the
   integration weights are known (an `S2Grid`, or an adapter's own polar/azimuth raster) they are
   used directly. Otherwise the spherical mean is estimated by resampling onto an equal-area grid
   and integrating with that grid's weights — approximate, so it is documented as such rather than
   presented as exact.

4. **A latent equality defect blocked the whole thing and is now fixed.** `MillerIndex`,
   `CrystalDirection`, `ZoneAxis` and `ReciprocalLatticeVector` inherited the dataclass `__eq__`,
   which compares their index arrays with `==` and raises "truth value of an array ... is
   ambiguous" for every distinct-but-equal pair. Comparing two poles — the first check any
   arithmetic operator must make — was therefore impossible. `SymmetrySpec` already carried the
   custom `__eq__`/`__hash__` for exactly this reason; those four now do too. `CrystalPlane` is
   still unhashable because `Phase` holds a dict; that is pre-existing and out of this sprint's
   scope, but it is why `PoleFigureDifference` compares poles rather than hashing them.

### Verification

Full suite green: 4340 tests, strict mypy over 106 source files, ruff clean. The identities that
calibrate the new numerics are asserted directly rather than by tolerance, and are the first
things to re-run if any of this is touched:

- Nadaraya-Watson reproduces a constant field to machine precision, for any halfwidth and any
  query direction (partition of unity — distinguishes a weighted mean from a weighted sum).
- Density estimation of a random texture converges to 1 m.r.d. at **second order** under cloud
  refinement. The test asserts the convergence *rate*, not a tolerance, because a wrong
  normalizing constant would leave a floor that refinement could not remove.
- `raster_solid_angle_weights` reproduces the analytic solid-angle mean of `cos(polar)` over a
  cap to 5e-4 at a 2.5 degree step, where the unweighted mean is wrong by more than 1e-2.
- On the real Ge(113) XRDML fixture the m.r.d. grid has weighted mean exactly 1 and a peak of
  1.2e4; the `max` mode reports a peak of 1.0 for a single crystal, which says nothing.

### What was deliberately not done

- No S2 kernel library to match the SO(3) one; the resampling kernel is one fixed von
  Mises-Fisher shape.
- No contoured rendering of a difference figure — scatter is honest for a scattered support but
  less readable on a dense grid.
- `CrystalPlane`/`Phase` remain unhashable (`Phase` holds a dict). Pre-existing, out of scope,
  and the reason `PoleFigureDifference` compares poles rather than hashing them.
- Ghost correction and a defocus model remain absent; they are the next texture gap, tracked in
  `docs/roadmap/feature_capability_review_2026_08.md` section 5.

## Follow-on: Notebook Policy And Tutorial 25 — COMPLETE (2026-08-08)

**Objective.** A hand-authored tutorial on pole-figure arithmetic, and a guarantee that no
notebook output is ever committed.

**The conflict that had to be resolved first.** `scripts/execute_notebooks.py` was not a
generator — the notebooks have always been hand-authored. It existed to bake outputs *into* the
committed file, because `docs/site/conf.py` set `nb_execution_mode = "off"` and myst-nb then
publishes only stored outputs. Stripping outputs without changing that would have turned all 24
tutorial pages into bare code listings. So the site now executes notebooks itself
(`nb_execution_mode = "cache"`, `nb_execution_raise_on_error = True`), which is a *stronger*
guarantee: a notebook that no longer runs fails the build, whereas a stored output proves only
that it ran once against whatever the library looked like then. Committed notebooks went from
13.3 MB to 0.45 MB.

**Two real defects found while writing the tutorial**, each fixed with its own commit and tests:

1. `HarmonicODF` crashed for any `degree_bandlimit >= 7` (`88eb2eb`). The Wigner small-`d`
   coefficient multiplies four factorials; at degree 7 the product exceeds int64 and became a
   Python big integer that NumPy could hold only as an object, so `np.sqrt` raised. The default
   is 6, so it broke exactly when a user raised the bandwidth for a sharp texture. Now evaluated
   in log-gamma.
2. `PoleFigureResidualReport.from_odf` compared two different scales (`05bbd0e`). A discrete
   ODF's `evaluate_pole_density` returns a kernel *response*, not m.r.d.; a random texture
   returns the kernel's spherical mean (~0.016 at 12 deg). A perfect fit reported a relative
   residual of 0.99. Now 0.010 on the same case, via the new public
   `pytex.texture.models.random_pole_density`, which also replaced the private duplicate in
   `diffraction/preferred_orientation.py` — the same constant had been discovered twice.

**Notebook 25** (`docs/site/tutorials/notebooks/25_pole_figure_arithmetic.ipynb`): 54 cells,
25 of them code, every result computed live on simulated rolling textures. Its spine is the
dependency chain — support, then scale, then arithmetic — and its most useful figures are the
two failure modes: the wrong resampling estimator returning 0.66 to 6.3 on a field that is
constant 1.0, and an over-smoothed ODF whose residual is *organized* (−4.17 m.r.d. at the peaks,
+0.74 in the background) rather than noisy.

### Next task

None claimed. The capability review's recommended order puts roadmap reconciliation plus
`windows-latest` in CI next, then the defocus model and ghost correction.
