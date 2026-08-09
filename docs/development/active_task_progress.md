# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history. Governed by the cardinal rule in `AGENTS.md`: ledger plus commit-and-push to `main`
after every substantial increment.

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

## Current Task: Orientation Representations, TEM Round-Trip Indexing, And CBED — IN PROGRESS (2026-08-09)

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
| 2 | Notebook 26: rotation and orientation representations | pending | |
| 3 | CBED module + tests | pending | |
| 4 | Notebook 27: TEM zone-axis indexing round trip (Ni, Zr) | pending | |
| 5 | Notebook 28: CBED analysis (Ni, Zr) | pending | |
| 6 | Docs index, symbol registry, worked examples, parity matrix | pending | |

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

### Next actions

Step 2: notebook 26, the rotation and orientation representation tutorial.

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
