# FIB Lamella Planning Foundation

**Status:** implemented in `pytex.fib` (PyTex 0.12.0); the end-to-end chain is validated analytically and
against `pytex.tem.navigation`, and remains **unvalidated against a real instrument** until a lamella cut
from a known EBSD scan has its achieved TEM tilt recorded (section 12, *Real case*).
**Target module:** `src/pytex/fib/`
**Sibling foundation documents:** `docs/architecture/tem_tilt_navigation_foundation.md`,
`docs/architecture/reference_frame_foundation.md`, `docs/architecture/ebsd_foundation.md`

---

## 1. Purpose

Answer one operational question, reproducibly and with its conventions stated:

> Given an EBSD scan of a bulk surface and a target TEM zone axis `<uvw>`, **which grain should be
> cut, where exactly in the SEM field of view should the lamella be placed, at what in-plane
> azimuth, and what residual tilt will the TEM holder then have to supply?**

The module is a **forward-geometry and selection layer**. It sits upstream of
`pytex.tem.navigation`, which already solves stage angles once a beam direction is known. It must
not duplicate that solver, `pytex.core.frames`, `pytex.core.symmetry`, or
`pytex.ebsd.models.GrainSegmentation`.

## 2. The governing physical constraint

The lamella is milled **strictly vertically** — the ion beam travels along the inward surface
normal. The two large faces of the slab are therefore **perpendicular to the surface**, so the
lamella plane normal is forced to lie **in the surface plane**.

In the TEM the electron beam enters along the lamella normal. Hence:

> The set of achievable lamella normals is the **great circle of in-plane directions** of the
> sample frame. A grain is usable for zone axis `<uvw>` **iff** at least one symmetry-equivalent of
> `<uvw>` lies within `eps_max` of that great circle. The angular distance of that equivalent from
> the great circle **is** the residual tilt the TEM holder must supply.

Everything in this specification follows from that statement. The crystallography is small; the
**convention bookkeeping across five frames is the product**.

## 3. Frozen decisions

These were settled by the maintainer and are not to be re-litigated during implementation.

| # | Decision |
|---|---|
| D1 | Readers: `.ang`, `.ctf`, `.oh5`/`.h5` via existing PyTex EBSD readers. **Assume the orientation file already carries the 70 deg tilt correction.** Do not apply a second correction. |
| D2 | Sample frame, SEM frame, SEM stage frame and TEM goniometer frame must all be **user-definable and changeable in the GUI**, not hard-coded. |
| D3 | SEM/EBSD registration **default**: SEM `x` along EBSD scan `x`, SEM `y` along EBSD scan `y`, unit scale, no flip. This default is configurable. An **optional, user-activated affine registration step** is provided. |
| D4 | FIB ion-column angle is a **user input**, accepting 52 deg and 54 deg, **defaulting to 54 deg**. |
| D5 | Milling is **strictly vertical**. Tilted or wedge milling is out of scope for v1. |
| D6 | There is **no deterministic lift-out / attach convention**. The in-plane mounting rotation `phi` about the lamella normal, and the front/back flip, are **uncontrolled and unknown** at planning time. This is a first-class modelling requirement, not a caveat — see section 7.4. |
| D7 | TEM holder envelope is a **user input**, defaulting to a rectangular envelope `alpha = +/-30 deg`, `beta = +/-30 deg`. |
| D8 | **Any symmetry-equivalent of the target `<uvw>` is acceptable.** Evaluate the full orbit, both senses, and report the best member. |
| D9 | Results must be presented **comprehensively, legibly and visually**, with geometry figures. A numeric answer without its figure is an incomplete deliverable. |

## 4. Explicit non-goals for v1

State these in the module docstring so scope creep is visible:

- tilted / wedge / pre-tilted milling (D5)
- bicrystal or interface lamellae where a zone axis is required in **both** crystals
- "boundary plane edge-on" as a co-equal objective (an *optional secondary score* is allowed, see
  section 7.6, but it must not constrain the primary solve)
- 3D / serial-section EBSD; the subsurface is assumed columnar and that assumption is **reported**,
  never silently made
- writing instrument-native pattern files for any FIB vendor; v1 emits human-readable work orders
  plus a JSON manifest

## 5. Reference frames and the transform chain

No private frame model. Register every frame below in `pytex.core.frames` / the frame catalogue and
route all conversions through `FrameGraph`. Register every new symbol in
`docs/standards/terminology_and_symbol_registry.md` and `pytex.core.symbols` **before** first use.

| Symbol | Frame | Definition |
|---|---|---|
| `C` | crystal | phase crystal frame, as PyTex already defines it |
| `S` | sample / specimen | `X_s`, `Y_s` in the surface plane; `Z_s` = outward surface normal. Tilt-corrected EBSD frame (D1). |
| `I` | SEM image | 2D pixel frame of the SEM or IQ image; related to `S` by in-plane rotation `omega`, optional handedness flip, scale, and optional affine (D3) |
| `P` | FIB stage / ion-beam view | the plan view seen by the ion column when the stage is tilted to the column angle `T` (D4) |
| `L` | lamella | `n_L` = lamella plane normal (in the surface plane), `t_L` = lamella long axis, `Z_s` = depth axis (mill direction is `-Z_s`) |
| `H` | TEM holder / goniometer | as already defined by `pytex.tem.stage` |

Chain:

```
  EBSD (.ang/.ctf/.oh5)        SEM / IQ image            FIB ion view            TEM holder
  orientation g(x,y) in S  ->  registration S<->I   ->   stage T, R_stage   ->   alpha, beta
        |                            |                         |                      |
        |  section 6 geometry        |  D3 + optional affine   |  D4 + calibrated    |  D6 phi unknown
        v                            v                         v  offset R0           v  D7 envelope
   eps, n_L, t_L            rectangle (x, y, L, W)      pattern azimuth        feasibility + risk
```

### 5.1 The FIB azimuth is calibrated, never asserted

At stage tilt equal to the column angle `T`, the **ion view is the undistorted plan view** and the
SEM view is foreshortened by `cos T` perpendicular to the tilt axis. The azimuth to enter into the
FIB pattern box is:

```
theta_ion = s_R * ( theta_S + R_stage ) + R0
```

where `s_R` in `{+1, -1}` is the chamber rotation handedness and `R0` is a **per-instrument offset
calibrated once with a fiducial**. Both are user inputs with documented defaults. The
implementation **must not hard-code a sign it cannot verify**; the docstring and `describe()` must
state that `s_R` and `R0` require a one-time fiducial calibration, and the module must ship a
documented calibration procedure. Getting this wrong destroys a real lamella, so the honest
position is a calibrated parameter, not a guessed convention.

## 6. Normative geometry

Let `g` be the measured orientation and `u_C` the unit target zone axis in the crystal frame.

**Verify against `pytex.core.orientation` whether PyTex's `g` maps sample -> crystal or
crystal -> sample, and write the correct one. Do not assume; a sign error here is undetectable in a
symmetric test case.**

```
  1. orbit          U = { s_i . u_C : s_i in point group }  union  { -U }      (D8)
  2. to sample      d_S = g^-1 . u   for each u in U
  3. out-of-plane   eps(u) = arcsin( | d_S . Z_s | )
  4. best member    u* = argmin eps ;  eps* = eps(u*) ;  d* = d_S(u*)
  5. lamella normal n_L = normalize( d* - (d* . Z_s) Z_s )
  6. long axis      t_L = Z_s x n_L
  7. surface azimuth  theta_S = atan2( n_L . Y_s , n_L . X_s )
```

The lamella footprint drawn on the surface is a rectangle of length `L` along `t_L` and width `W`
along `n_L`, extending to depth `D` along `-Z_s`.

Defaults (all user-editable): `L = 15 um`, `W = 2 um` (trench-to-trench including the protective
strap), `D = 8 um`, final thickness `t = 80 nm`.

Vectorize steps 1-4 over (grains x orbit members) with `einsum`/broadcasting. No per-grain Python
loop on the hot path.

## 7. Required behaviour

### 7.1 Single-grain plan

Input: a grain (or a point), a phase, a target `<uvw>`, geometry settings. Output: a
`LamellaPlan` carrying `eps*`, the chosen orbit member, `n_L`, `t_L`, `theta_S`, `theta_I`,
`theta_ion`, the placement rectangle, the predicted TEM residual, the feasibility verdict, the risk
flags and the figures.

### 7.2 Map-wide ranking (the primary workflow)

Given a segmented map and a target `<uvw>`, return grains ranked by a transparent, documented score
combining: `eps*` (primary), whether the footprint fits inside the grain, grain orientation spread
(GOS) as an orientation-reliability proxy, distance from the map edge, and neighbour context. The
score's weights are user-visible inputs and the `describe()` output must state them. Also produce a
per-point **preparability raster** of `eps*` for display over the SEM/IQ image.

### 7.3 Footprint fit inside the grain

The azimuth `theta_S` is **not free** once the orientation is fixed, so the only placement freedom
is translation. Find the placement maximising clearance: rotate the grain mask by `-theta_S`, then
locate the largest inscribed axis-aligned rectangle of the required `L x W` using a summed-area
table or `scipy.ndimage` distance transform. Pure NumPy/SciPy only — **no new dependencies**
(air-gapped deployment constraint). Report the achieved margin in micrometres, and if the footprint
does not fit, say so and report the largest `L` that does.

### 7.4 The unknown mounting rotation — a first-class output

Because `phi` is uncontrolled (D6), `eps*` is a known **angle** whose decomposition into
`(alpha, beta)` is unknown until the lamella is physically on the holder. The module must therefore
report **three** distinct things and never conflate them:

1. **Guaranteed reachable** — `eps* <= min(alpha_max, beta_max) - margin` (default margin 5 deg).
   Reachable for every `phi`. This is the recommended selection criterion.
2. **Probabilistically reachable** — sweep `phi` over `[0, 2pi)`, solve the true stage angles at
   each `phi` with the existing `pytex.tem.navigation.plan_tilt_to_zone_axis` and
   `pytex.tem.stage.TiltEnvelope`, and report the **fraction of `phi` for which the target is
   inside the envelope**. Report it as a number and as a polar feasibility figure.
3. **Unreachable** — no `phi` works.

The small-angle decomposition `alpha ~ eps cos phi`, `beta ~ eps sin phi` may appear in explanatory
prose as a stated approximation. It must **not** be the computational path; use the exact solver.

The front/back **flip** `s = +/-1` reverses the beam sense, giving `-[uvw]`. Report both branches
explicitly. Note in `describe()` that for centrosymmetric diffraction conditions the distinction is
usually immaterial, but for CBED and polar/non-centrosymmetric phases it is not.

### 7.5 Uncertainty budget

Propagate, and report as a combined standard uncertainty on `eps*`:

- EBSD angular accuracy (user input, default 0.5 deg)
- intragranular orientation spread for the chosen grain (GOS/KAM, from the map — already available)
- registration residual from the affine step, when used
- stage/mount repeatability (user input, default 2 deg)

Report `eps* +/- u(eps*)` and let the feasibility verdict use the upper bound, not the point
estimate. A plan that is feasible only at the point estimate must be flagged.

### 7.6 Optional secondary score (non-constraining)

Allow the user to name a secondary plane (for example a grain-boundary plane or a habit plane) and
report, for information only, the angle between that plane and the beam direction — i.e. how close
to edge-on it would be in the resulting lamella. It **ranks** candidates; it must never alter the
primary solve. Full multi-objective optimisation is deferred.

## 8. Public API sketch

```
src/pytex/fib/
  __init__.py
  frames.py      # frame registrations, SEM<->EBSD registration, chamber geometry
  geometry.py    # section 6, vectorized, frame-typed, no naked arrays on the public surface
  placement.py   # section 7.3 footprint fit
  planning.py    # LamellaPlan, single-grain solve, phi sweep, uncertainty budget
  selection.py   # map-wide ranking, preparability raster
  report.py      # LamellaPlanReport, describe(), to_json_dict()
src/pytex/plotting/fib_figures.py
src/pytex/app/services/fib_lamella.py
```

Immutable dataclasses with construction-time invariant checks:

- `SurfaceGeometry` — sample frame axes, surface normal sense
- `ImageRegistration` — rotation, scale, flip, optional affine; must expose `residual_deg`
- `ChamberGeometry` — column angle `T` (D4), `s_R`, `R0`, compucentric flag
- `LamellaSpec` — `L`, `W`, `D`, final thickness
- `MountModel` — `phi` unknown (D6), flip branch handling
- `LamellaPlan` — one candidate, with `eps`, orbit member, azimuths, rectangle, TEM residual,
  feasibility class, uncertainty, risk flags
- `LamellaPlanReport` — ranked plans plus what is undetermined and what to do about it

Every stable report object carries `describe()` producing convention-explicit, citation-backed
prose, and `to_json_dict()` kept in lockstep with it, per the explainable-results doctrine. Mirror
the existing style of `pytex.tem.navigation.TiltPlanReport` exactly.

## 9. Mandatory figures

Publication-quality SVG, via `pytex.plotting`, following the visualization style guide. All of
these are part of the deliverable, not decoration (D9):

1. **Frame-convention figure** (canonical, committed to `docs/figures/`): sample, SEM image, FIB
   stage and lamella frames with the mill direction and the in-plane constraint drawn.
2. **Plan-view overlay** — SEM or IQ image with grain outline, placement rectangle, azimuth
   annotation, scale bar, and the lift-out direction.
3. **Section schematic** — side view of the trench, slab, depth `D`, thickness `t`, showing the
   normal lying in the surface plane.
4. **Stereographic projection** — orbit members of `<uvw>` plotted against the in-plane great
   circle, with `eps*` annotated for the selected member. This is the figure that makes the
   criterion in section 2 self-evident and it should be treated as the module's signature figure.
5. **Polar `phi`-feasibility plot** — reachable / unreachable arcs of mounting rotation (section
   7.4), with the guaranteed-reachable disc drawn.
6. **Preparability map** — per-point `eps*` raster over the ROI with a sequential colormap and an
   explicit colour scale.
7. **Predicted SAED** at the achieved axis, reusing the existing diffraction code.

## 10. Workbench panel

One panel, generated from a self-describing operation manifest like every other panel. It must obey
the cardinal form rule: everything needed on one screen, declared field widths, registry symbols on
labels, related quantities grouped by `row`, every mandatory input reachable without a disclosure.

Input groups: file + phase; target `<uvw>`; frames and registration (D2, D3); chamber (D4); lamella
dimensions; holder envelope (D7); uncertainty inputs; ranking weights.

Output: ranked candidate table, the figures of section 9, the `describe()` prose, and a
**printable work order** stating grain id, stage coordinates, `theta_ion` with its calibration
caveat, rectangle dimensions, expected TEM residual with its uncertainty, and the feasibility class.

Ship runnable examples the test suite executes, as with every other panel.

## 11. Contracts

Add `schemas/fib_lamella_plan.schema.json` and validate `to_json_dict()` against it in tests. Follow
`docs/standards/data_contracts_and_manifests.md`. Preserve existing schemas; add, do not restructure.

## 12. Validation strategy

**There is no MTEX equivalent for this workflow, so MTEX parity is not available as the floor.**
State that explicitly in `docs/testing/strategy.md` rather than leaving a silent gap. The
substitute, in descending order of strength:

| Lane | Test | Why it is trustworthy |
|---|---|---|
| Analytic | Orientations constructed so `eps` is exactly known (target axis placed exactly in-plane -> `eps = 0`; exactly along `Z_s` -> `eps = 90`; at a constructed 12.5 deg) | closed-form ground truth, not a prior program output |
| Invariance | Applying any symmetry operator to `g` leaves `eps*`, `n_L` and `t_L` invariant up to the symmetry | catches orbit and sense bugs |
| Round trip | forward geometry composed with the inverse recovers the input to 1e-9 | catches frame-composition errors |
| Registration | identity registration is a no-op; a known rotation+flip round-trips; affine residual is reported | catches the S<->I bugs that ruin real lamellae |
| Cross-check | For a plan with `eps*`, feeding `n_L` to `plan_tilt_to_zone_axis` must return the same residual angle | ties the new module to the validated existing solver |
| Envelope | `phi` sweep feasibility fraction is 1.0 for `eps < min(alpha_max, beta_max)` on a rectangular envelope, and matches a closed-form arc computation on an asymmetric one | catches section 7.4 errors |
| Golden prose | snapshot tests on `describe()` | the doctrine requires the prose be tested like any other output |
| Panel | Playwright, as for every other panel | |
| Real case | **Requires maintainer data**: one lamella cut from a known EBSD scan with the achieved TEM tilt recorded. Until it exists, mark the end-to-end chain **unvalidated against a real instrument** in the foundation document and in `describe()`. |

Pytest and Playwright only; introduce no other framework. Warnings are defects; close matplotlib
figures; coverage must not decrease.

## 13. Documentation obligations

- this document committed as `docs/architecture/fib_lamella_planning_foundation.md`, linked from
  `docs/README.md` and `README.md`
- a theory note under `docs/site/theory/` deriving section 6 with the mathematics rendered, citing
  normative sources per the citation policy and source hierarchy
- a workflow page under `docs/site/workflows/` covering the panel end to end
- executable worked examples in `worked_examples/`, values computed at build time with independent
  analytic provenance
- every new symbol registered centrally before use
- the FIB `s_R`/`R0` calibration procedure written as a standalone, followable instrument protocol

## 14. Implementation phases and success criteria

Land each phase on `main` with the progress ledger updated in the same commit; each phase must leave
`ruff`, `mypy` and `pytest` green.

**Phase 1 — geometry core.** `fib/frames.py` + `fib/geometry.py`, single-grain solve, no I/O, no
plotting.
*Success:* analytic, invariance, round-trip and cross-check lanes green; `eps` for a constructed
12.5 deg case matches to 1e-9; the orientation-convention question of section 6 is resolved **with
a test that would fail under the opposite convention**.

**Phase 2 — plan and report.** `LamellaPlan`, `LamellaPlanReport`, `describe()`, `to_json_dict()`,
JSON schema, `phi` sweep, uncertainty budget, flip branches.
*Success:* schema validation green; golden prose tests green; feasibility fraction matches the
closed-form arc on an asymmetric envelope; every number in `describe()` also present in
`to_json_dict()`.

**Phase 3 — map integration and placement.** EBSD readers, registration (default + optional
affine), footprint fit, map-wide ranking, preparability raster.
*Success:* identity registration is a proven no-op; footprint fit reproduces a hand-computed
rectangle on a synthetic mask; ranking is stable and its weights appear in `describe()`; runs on a
full-size map without a per-grain Python loop.

**Phase 4 — figures.** All seven figures of section 9; the frame-convention SVG committed as a
canonical asset.
*Success:* byte-comparison test on the canonical SVG; every figure renders in both a light and a
dark theme; no open-figure or warning leaks.

**Phase 5 — panel, docs, worked examples.**
*Success:* manifest test green including the field-width and `row` rules; Playwright green; worked
examples regenerate and pass; `docs/README.md` index complete; work order prints legibly on one
page.

## 15. Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| FIB azimuth sign / offset wrong | **destroys a real lamella** | calibrated `s_R`, `R0` (section 5.1); caveat carried in `describe()` and on the printed work order; fiducial protocol documented |
| EBSD -> sample frame convention wrong for a vendor | high, silent | per-reader convention test with an asymmetric fixture that fails under a flip |
| `phi` unknown treated as zero | high | section 7.4 is a hard requirement; a plan may not report a single `(alpha, beta)` pair as if `phi` were known |
| Subsurface orientation assumed | medium, unavoidable | always reported as an assumption with a risk flag; never silent |
| Orientation convention `g` vs `g^-1` | high | Phase 1 success criterion requires a convention-discriminating test |
| Over-trusting the point estimate of `eps` | medium | feasibility verdict uses the upper bound of the uncertainty budget |

## 16. Novelty classification

Required by project practice; keep this section in the committed document.

- **Established methods** — orientation-to-sample-direction transform, symmetry orbits,
  angle-to-plane, zone-axis tilt solving, inscribed-rectangle placement. Nothing here is new.
- **Useful engineering improvements** — the convention-explicit EBSD -> FIB -> TEM chain in one
  reproducible tool; map-wide preparability ranking instead of per-grain eyeballing; the printable
  work order.
- **Potentially original research directions** — (a) the **uncertainty budget** of section 7.5
  propagating EBSD accuracy and intragranular spread into a predicted residual tilt with a
  confidence interval; (b) the **`phi`-marginalised feasibility** treatment of section 7.4, which
  turns an uncontrolled mounting rotation into a quantified probability rather than an ignored
  variable. Both are modest but, as far as is known here, not standard practice.
- **Requires literature verification, before any novelty is claimed in writing** — EBSD-guided
  site-specific lamella preparation is described in the microscopy literature and parts of this
  chain exist in vendor software. Whether an open tool already covers the full chain, and whether
  (a) or (b) above have prior art, has **not** been verified. Run that search and record the result
  in `references/reference_index.md` before the module's documentation asserts novelty.
