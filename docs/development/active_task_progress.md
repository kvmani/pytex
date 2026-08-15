# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history. Governed by the cardinal rule in `AGENTS.md`: ledger plus commit-and-push to `main`
after every substantial increment.

## Repository Governance And Five-Feature Delivery Program — IN PROGRESS (2026-08-15)

**Objective.** Repair the governing roadmap/ledger and executable quality gates, add a minimal
critical Playwright lane for the shared web/desktop workbench, then deliver five deliberately
bounded scientific features in order: measured powder-XRD I/O and comparison, random-standard
pole-figure defocus calibration, hex-grid EBSD support, a finite-thickness SAED shape factor, and
named-component fitting for measured ODFs. Each feature lands with tests, an independently known
numerical result, `describe()`, a portable JSON contract where scientifically appropriate, a
worked example, theory/workflow documentation, a parity-ledger update, and a benchmark case.
Where an open redistributable measurement is unavailable, the fixture is synthetic and labelled
as such in both its metadata and documentation.

### Governing decisions

1. Hardening precedes feature breadth: stale roadmap claims, the contradictory application-ledger
   state, Windows CI, the Sphinx-warning ratchet, and browser automation land before scientific
   feature 1.
2. The warning ratchet measures the ordinary CI Sphinx build first. It forbids warning growth
   immediately while the existing warning families are reduced deliberately; it is not a claim
   that the baseline warnings are acceptable.
3. Playwright covers critical user-visible behavior in a real Chromium browser without changing
   the workbench's zero-bundler, zero-third-party-runtime architecture.
4. Each scientific increment is independently revertible and leaves the full base lane green.

### Step ledger

| # | Increment | Status | Commit |
| --- | --- | --- | --- |
| G1 | Reconcile governing roadmaps and completed application ledger state | done | this commit |
| G2 | Add Windows base CI and a tested Sphinx-warning ratchet | done | this commit |
| G3 | Add minimal critical Playwright workbench tests and CI lane | pending | — |
| F1 | Measured powder-XRD I/O and comparison | pending | — |
| F2 | Random-standard defocus calibration | pending | — |
| F3 | Hex-grid EBSD support | pending | — |
| F4 | Finite-thickness SAED shape factor | pending | — |
| F5 | Named-component ODF fitting | pending | — |
| C | Full completion audit across code, contracts, docs, examples, ledgers, and benchmarks | pending | — |

### Current evidence and resume point

- Clean `main` at `33a4f45`; `origin/main` matches.
- Repository integrity, Ruff lint, and strict mypy over 138 source files are green at entry.
- The live tree contains roughly 82.5k Python source lines and 558 top-level exports; the July
  governing review still reports 31.7k lines, 804 tests, and no CI coverage.
- The ordinary Sphinx build passes at the enforced 602-warning ceiling. A nitpicky audit builds
  successfully but reports 2,300 warnings, including genuine unresolved current API references;
  those are tracked as documentation debt beyond the no-growth baseline.
- Windows/Python 3.11 now participates in the base CI lane; Linux and macOS retain the 3.11–3.13
  matrix. The ratchet parser has focused unit coverage and contributor commands use the same gate.
- Full-suite verification exposed one pre-existing generated-asset drift in
  `docs/figures/class_model_ebsd.svg`; the canonical generator refreshed it and its byte-exact test
  now passes.
- Governance verification: repository integrity, Ruff, strict mypy over 138 source files, the
  focused governance/class-atlas tests, the full pytest suite (two expected skips), and the
  602-warning Sphinx build pass.
- Next action: commit and push G1/G2 after the final suite passes, then start G3 by adding the
  test-only Playwright dependency, critical Chromium journeys, and the dedicated CI lane.

## Crystal Sphere Lighting And Depth — COMPLETE (2026-08-15)

**Objective.** Replace the Crystal Viewer's flat atom discs and flat bond strokes with a
professional sphere/cylinder presentation: configurable studio lighting, diffuse and specular
response, depth cueing and soft shadows in the byte-identical web/desktop viewer, with the same
lighting parameters carried into publication export.

**Initial evidence.** On clean `main` at `0fbf4b1`, browser inspection shows atoms as uniform SVG
`circle` fills with a thin outline; their radius and painter-order occlusion convey size/depth, but
there is no surface-normal cue, highlight, limb darkening or depth fade. Bonds are uniform SVG
lines. The Python publication renderer already has lit sphere/cylinder meshes and depth cueing, but
the shared appearance envelope exposes none of those theme keys, so the interactive and exported
looks cannot be tuned together.

**Design.** Keep orthographic geometry and crystallographic semantics unchanged. The interactive
renderer will use reusable SVG radial gradients for sphere lighting and perpendicular linear
gradients for cylindrical bonds, with bounded depth buckets to avoid one definition per atom.
Lighting controls will cover ambient/diffuse/specular strength, highlight size, azimuth/elevation,
depth cue and soft shadow. The export envelope will validate and map those same values to the
existing Blinn–Phong/depth-cue renderer; light direction will be transformed from screen space back
to the crystal frame using the exported camera matrix.

**Implemented.** The shared viewer now paints atoms with reusable, species-aware radial gradients,
bonds with layered cylindrical glyphs, and distant atoms with an opaque background veil that
preserves painter-order occlusion. The property rail exposes glossy/matte/flat finish, light
azimuth/elevation, ambient/diffuse/specular strengths, highlight sharpness, and depth cue. Export
normalizes and validates the appearance envelope, transforms the screen-space light through the
camera matrix, and maps it to the existing publication renderer. Publication specular response is
now camera-aware in both the single-crystal and composed-world render paths.

**Verification so far.** `node --check` passes. Targeted Ruff and 114 focused app/rendering tests
pass. Browser validation at `127.0.0.1:8765` confirms 27 sphere glyphs, two reused species
gradients, 54 cylindrical bonds, depth values that vary with projected distance, a genuinely flat
fallback with no gradients, and a highlight that moves from 35% to 65% across the glyph when light
azimuth changes. Light and system-dark screenshots keep the labels legible and give the crystal
the intended three-dimensional reading. A 150 dpi NaCl publication export was decoded and visually
inspected: its sphere highlights, cylindrical bonds, depth ordering, plane, labels, legend, and
frame indicator are all intact.

**Final verification.** The full base lane is green: `ruff check .`, mypy over all 138 source files,
and the complete pytest suite (two pre-existing skips) pass in 704 seconds. Native pywebview
inspection at 1536×816 confirms the same glossy sphere/cylindrical-bond presentation as the web
shell, a large uninterrupted crystal canvas, independently scrollable input/property rail, visible
completion status, and usable lighting controls. The native Surface finish selector was exercised
from glossy to flat and back to glossy; the scene redrew immediately and no calculation or shell
restart was triggered. No generated inspection image or runtime output was added to the repository.

**Outcome.** The flat-disc defect is resolved across both application shells and publication
export. All controls are presentation-only, documented, bounded by service validation, and covered
by interaction/source contracts plus camera-aware numerical shading tests. This task is complete;
no follow-up work remains in this ledger entry.

## Scientific Workbench Visual Interaction And Layout — COMPLETE (2026-08-14)

**Objective.** Bring the byte-identical desktop/web workbench closer to a professional texture and
crystallography application: mouse-wheel zoom and direct manipulation in every graphic area,
presentation-only property controls for scientific scene objects, real contour-line and filled-
contour texture plots with user-controlled levels, and a layout pass that reserves space for the
plot, controls, results, progress, and operational messages without hiding any of them.

**Governing boundary.** Desktop and web remain one application (`docs/architecture/application_platform.md`):
all changes belong in the shared static frontend or shared Python services, never in one launcher.
Presentation controls redraw existing scientific rows/scene geometry and do not masquerade as
calculation parameters. Scientific inputs and export provenance remain manifest-owned.

### Initial audit (clean `main` at `388ac5b`)

- `core/plotframe.js` supplies the mandated cursor and hover-detail surfaces, but has no viewport
  model. Only the Crystal Viewer implements wheel zoom, locally, so every 2-D SVG plot is static.
- The Crystal Viewer already receives explicit atom, bond, cell, plane, and direction glyphs and
  correctly separates camera arithmetic from crystallography, but presentation properties are
  hard-coded in `crystal.js` (atom radius, plane opacity, bond/cell/direction widths, labels).
- Texture calls pole-figure discs and ODF raster cells “filled contours”. They are mosaics over
  discrete samples: there is no interpolated isoline geometry, line-only/filled/combined mode, or
  adjustable contour-level policy. The service already returns the numerical grid/support needed
  to build those display layers without changing the ODF calculation.
- The shared shell already has a wide-screen control rail, plot stage, result cards, loading button
  state, error toasts, and responsive breakpoints. Interactive validation must check whether the
  new appearance sections and viewport toolbar keep those regions usable at desktop and tablet
  widths rather than assuming the existing grid is sufficient.

### Planned increments

| # | Increment | Status | Commit |
| --- | --- | --- | --- |
| 1 | Shared 2-D SVG viewport: wheel zoom, drag pan, reset/fit, cursor-correct coordinates | done | `579dd6e` |
| 2 | Crystal object appearance editor with live redraw and publication-export parity | done | `210d5fd` |
| 3 | Professional texture isolines/filled contours, adjustable level count/values and colour scale | done | `243acc5` |
| 4 | Homogeneous interaction/layout pass: plot allocation, property/result/log/progress space | done | `1598fd5` |
| 5 | Browser + native desktop interaction validation, automated quality lanes, docs closeout | done | (closeout commit) |

### Increment 1 result

`core/plotframe.js` now owns the 2-D viewport as part of the instrument contract: wheel zoom is
anchored at the pointer, Shift-left-drag and middle-drag pan without stealing ordinary plot clicks,
the header reports magnification and supplies Zoom in/out plus Fit, and viewBox coordinates remain
the source of the cursor conversion after every camera change. The optional preserved-viewport
path lets a presentation-only redraw retain the user's inspection location. Crystal explicitly
disables the 2-D viewport because its existing wheel/drag handlers are a 3-D camera.

**Verification so far.** `node --check` is green for the shared plot frame and crystal panel;
`pytest tests/unit/test_app_server.py tests/unit/test_app_desktop.py -q` passes (60 tests). Ruff
passes and strict mypy passes over 138 source files. The full unit lane reached 100% with one
unrelated stale canonical asset: the EBSD class-model generator produced a 1984 px canvas while the
committed SVG still declared 1975 px. Per the failing test's required repair path,
`scripts/generate_class_model_figures.py` regenerated the referenced atlas; only
`class_model_ebsd.svg` changed, and all eight generator/committed-asset comparisons now pass. Every
other unit test in the controlling run passed (two existing skips).

**Exact next action.** Commit and push increment 1 to `main`, then build the Crystal Viewer
appearance editor against the scene glyphs and publication renderer.

### Increment 1 landing

Committed and pushed to `main` as `579dd6e` (`Give every 2D workbench plot a shared viewport`).

### Increment 2 design

The editor will expose object-class visibility plus atom size/opacity and per-species colours;
bond, cell, plane and direction colour/opacity/width; annotation scale; and axis-gizmo visibility.
It redraws the JSON scene already returned by Python. `crystal.render` receives a validated nested
appearance object and translates the same settings into the shared YAML-style renderer; object
styles remain outside the scientific inputs/rows and never trigger scene recomputation in the UI.

**Exact next action.** Add the appearance model/control and renderer bindings, then pin both the
frontend-only redraw boundary and the validated publication-export mapping in focused tests.

### Increment 2 result

The Crystal Viewer now has one open **Object properties** section in the control rail. It toggles
atoms, bonds, cells, planes, directions, labels and the axis gizmo; scales/fades atoms; exposes a
colour picker for every species in the current scene; and controls colour, width and opacity for
bonds/cells/planes/directions plus annotation scale. Every input mutates presentation state and
calls `draw()` only—never the scientific service. The Figure path serializes the same state into a
strictly validated `appearance` object. Python maps it to the shared crystal style, filters hidden
overlay classes, and now supports theme-level per-species colour overrides.

**Verification.** JavaScript syntax passes. The focused Crystal/app-server/crystal-renderer suite
passes (96 tests). The first focused ruff call found one 102-character validation line; it was
wrapped immediately. The full base lane is green: ruff passes, strict mypy passes over 138 source
files, and the entire unit suite passes with the two existing skips.

**Exact next action.** Run the controlling quality lane, commit and push increment 2, then replace
the Texture panel's disc/raster mosaics with real configurable filled/line contour geometry.

### Increment 2 landing

Committed and pushed to `main` as `210d5fd` (`Add live crystal object property controls`).

### Increment 3 design

Pole-figure samples will be resampled for display only onto a clipped regular projection grid;
ODF sections already have a rectangular grid. One shared marching-squares implementation will draw
actual isolines over either grid. Filled contours quantize the same interpolated grid into the same
declared bands, so line, filled and combined views cannot disagree about levels. Controls will
offer automatic level count or explicit comma-separated levels, upper colour limit, palette, line
colour/weight, fill opacity and display-grid resolution. The source m.r.d. rows and their export
remain unchanged, and transparent hit regions retain hover readout of computed—not interpolated—
sample values. The current SVG becomes directly exportable through the shared desktop/web saver.

### Increment 3 verification

The Texture workbench now renders pole figures and ODF sections through one SVG contour pipeline:
filled bands and marching-squares isolines share the same interpolated field and declared levels.
Automatic and exact levels, colour upper limit and palette, line colour/weight, fill opacity and
display-grid density update live without recomputing scientific results. SVG export preserves the
configured presentation while CSV/workbook/JSON exports continue to carry computed source values.

Focused Texture/server tests pass. In the in-app browser, the Goss reference example produced all
eight requested automatic isolines; exact levels `0.5, 1, 2, 4` produced four matching SVG paths;
the plot zoomed to 228% under mouse-wheel input and Fit restored 100%; SVG export completed with no
browser warnings or errors.

**Exact next action.** Run the controlling quality lane, commit and push increment 3, then add a
shared operation-progress/activity surface and tune responsive plot/property allocation before the
final browser and native-desktop validation pass.

### Increment 3 landing

Committed and pushed to `main` as `243acc5` (`Add configurable texture contour plots`). The full
base lane passed: ruff, strict mypy over 138 source files, and the complete unit suite with the two
existing skips.

### Increment 4 design

The one shared API call path now emits uniquely identified start/finish events, including duration
and failure text. A persistent bottom activity bar reports current work without consuming plot
space; its bounded, expandable history overlays the stage only when requested. This makes progress
and recent operational messages available in every panel and both shells without duplicating panel
logic. The large-screen rail now scales between 21 and 25 rem. At mid-size desktop widths the
duplicated tagline and then action labels yield before tabs wrap; tablet behaviour still places a
height-bounded property sheet below the main graphic.

**Exact next action.** Reload the live app and verify current/finished/error activity, layout and
navigation at wide, laptop, tablet and phone widths; then source-test, commit and push increment 4.

### Increment 4 verification

At the 1164 × 655 in-app-browser viewport all seven workspaces remain on one 64 px masthead row;
the plot receives the flexible width, the rail receives 29% within its 21–25 rem bounds, and the
36 px activity bar remains visible below both independent scroll regions. Opening history showed
the catalogue and crystal operations with durations and did not reflow either region. Browser logs
contained no warnings or errors.

After installing the repository's declared `desktop` extra, Computer Use exercised the real
pywebview window rather than the browser fallback. The native Texture panel exposed all contour
controls; a real wheel gesture over the pole figure reached 192%, Fit returned it to 100%, and the
history panel listed pole-figure, crystal and catalogue calls with timings. Windows snap testing at
947 px confirmed all seven tabs and the below-plot property sheet remain reachable. The resulting
breakpoint pass gives the plot more vertical space (rail capped at 36 vh, 34 vh on phones), moves
all tabs to a deliberate full row at tablet width, and opens native windows maximized so the
activity bar cannot begin under the OS taskbar.

Focused app-server/desktop tests, JavaScript syntax checks, ruff and strict mypy are green. The
complete unit lane also passes: 6,002 passed and the two existing skips in 15:26. An initial
buffered invocation exceeded a 15-minute outer tool timeout without reporting a failure; an
unbuffered controlling rerun exposed continuous test-by-test progress and completed successfully.

**Exact next action.** Run the complete unit lane, commit and push increment 4, then relaunch from
that commit for final browser/native smoke coverage and close this ledger and goal.

### Increment 4 landing and final closeout

Committed and pushed to `main` as `1598fd5` (`Add shared calculation activity and responsive
layout`). A committed-build browser reload at 1164 × 655 rendered two random-baseline contour paths,
kept all tabs visible, allocated 826 px to the stage and 337 px to the rail, recorded the pole-figure
call, and emitted no warnings/errors. The committed native build relaunched maximized at 1536 × 816
with the activity bar inside the OS work area, all tabs visible, and the same Texture contour
controls. Earlier real native interaction on this exact viewport implementation reached 192% by
wheel and returned to 100% with Fit; the browser reached 228% and likewise returned to 100%.

Final quality evidence: JavaScript syntax checks pass; ruff passes; strict mypy passes over 138
source files; the focused app-server/desktop/integrity suites pass; and the complete unit suite
passes with 6,002 tests and the two existing skips. The optional `desktop` extra was installed only
into the local virtual environment to exercise pywebview; it added no repository artifact. The
implementation is complete and no follow-up code work is required for this objective.

## TEM Indexing: Lattice Overlay, Centre Refinement, Scored Solutions — COMPLETE (2026-08-14)

**Objective.** Close the loop between picking and trusting an answer. Five things, in the order a
microscopist meets them:

1. **A live 2D lattice overlay** fitted to the picked spots, drawn over the pattern, so a wrong
   beam centre or a mis-picked spot is visible rather than inferred.
2. **Centre refinement** — least-squares from the fit, plus small manual nudges — with the overlay
   updating live so the user can judge when it is right.
3. **Explicit deviations** for every candidate solution: measured against calculated d-spacings and
   interspot angles, not just a single residual.
4. **A fused, configurable, documented accuracy score** that ranks the candidates.
5. **The calculated pattern superimposed on the experimental one** for a chosen solution, an
   explicit accept, and the accepted solution carried into tilt planning.

Plus a cross-cutting requirement: **detailed export in both human-readable and machine-readable
form at every stage, across every panel.**

### Initial audit (worktree clean on `main` at `fb0df52`)

- Picking is blind. Nothing is drawn from the picks except the picks, so the beam centre is placed
  by eye and never checked. `PatternCalibration.centre` is taken as given.
- `PatternSolution` reports `residual_inv_angstrom` per spot and a `(matched_fraction, mean
  residual)` sort key. Ranking is documented as *a sort key, not a scalar quality* — deliberately,
  and the app says so. There is no per-spot d-spacing deviation, no interspot-angle deviation, and
  no single number a user can weigh.
- Nothing renders a calculated pattern over a measured one anywhere in the application.
- Export exists on every table-bearing result — CSV, XLSX, JSON — from one shared path, plus a
  "Copy summary" button. The gap is a genuinely *human-readable* document: JSON is machine-only,
  CSV is a grid. `EXPORT_FORMATS` is not published in the manifest, so the frontend hard-codes
  three buttons.

### Planned increments

| # | Increment | Status | Commit |
| --- | --- | --- | --- |
| 1 | `pytex.diffraction.lattice_fit`: 2D lattice fit with least-squares centre refinement | done | (this commit) |
| 2 | `pytex.diffraction.solution_scoring`: deviations and a configurable fused score | done | (this commit) |
| 3 | App: `tem.fit_lattice`, scoring and overlay data in `tem.solve_pattern` | done | (this commit) |
| 4 | Frontend: overlay, nudges, scored solution list, calculated-pattern fit, accept | done | (this commit) |
| 5 | Export: a human-readable report format on every panel, manifest-published | done | (this commit) |
| 6 | Docs, human-style pass, closeout | done | (this commit) |

### Increment 1 result

`pytex.diffraction.lattice_fit` imposes the one constraint a zone-axis pattern satisfies before any
crystallography — that its spots lie on a plane lattice — and gets two things out of it: a beam
centre that is over-determined by the picks rather than guessed from one of them, and a name for
any pick the lattice cannot explain. 39 tests, every case built from a lattice whose answer is known
before the fit runs.

**The first four drafts were each wrong in an instructive way, and driving them is what found it.**

1. **Seeding from offsets to the picked centre subdivided the cell.** A centre half a spacing out
   makes some offsets spuriously short; the shortest pair then generates a sub-lattice; and the fit
   explains every spot perfectly by halving the cell *around the wrong centre* — machine-precision
   residuals, an uncorrected centre, and nothing anywhere to say something was wrong. It failed in
   exactly the case the module exists to repair. Seeding from *differences between spots* cannot do
   that: a difference of two spots is a lattice vector however badly the centre was picked.
2. **One shortest difference is not enough either.** A spot clicked forty pixels off creates short
   differences that are not lattice vectors, and a fit seeded from one of those is confident about a
   lattice that does not exist. The seed is now a small search — six candidates, each fitted, the
   best kept.
3. **Ranking candidates by inlier count always prefers the finer lattice**, because halving a cell
   explains every spot it explained before plus the mis-picked one. Candidates are now ranked by
   evidence: a cell of area `A` puts about `π t² / A` of the plane within tolerance of some node, so
   an inlier is worth `log(A / π t²)` nats and a lattice half as coarse pays `log 4` per inlier for
   its extra nodes. The tolerance must be shared across candidates or the comparison is empty —
   which is why it is a fraction of the shortest *observed spot separation*, a stand-in for picking
   precision, and not a fraction of each candidate's own basis.
4. **The reported basis angle was not a property of the lattice.** A square lattice came back as two
   vectors 135° apart, which is correct and useless. The basis is now Gauss-reduced, so lengths and
   included angle are invariants a user can read and 90° means rectangular.

**Two honest limits, both now pinned by tests.** A centre wrong by an *exact lattice vector* is
undetectable from geometry alone — every spot is still an exact node — and what identifies the
transmitted beam is that it is the brightest thing on the plate, which is a judgement about
intensity. And a centre more than half a spacing out cannot be refined without relabelling which
node the origin is, so the fit is leashed there and says so rather than quietly choosing a different
origin.

Against the three practice plates: the true centre is recovered to within 3 px (the plates' own
sub-pixel scatter) from offsets of 0, 12, 28 and 30 px, and a 47 px mis-pick is flagged in all nine
placements tried.

### Increment 2 result

`pytex.diffraction.solution_scoring` keeps the evidence and the opinion apart, which is the whole
design.

**Deviations are measurements.** Per indexed spot, measured d against calculated d, absolute and
relative. Per pair of indexed spots, measured angle against calculated angle. No weighting, no
judgement — the numbers a user would read off the plate themselves.

**The score is a policy**, and it lives in `ScoringWeights` where it can be read and changed.
Length, angle and coverage each map to [0, 1] through `1 / (1 + (x/t)^s)`, which scores one half at
the stated tolerance; the fused score is their weighted mean, normalised so 1 means perfect
agreement on everything picked and 0.5 means disagreement at tolerance. Every default is documented
with its reasoning, and the score carries the weights that produced it — a number whose policy is
invisible is an assertion, not a measurement.

**Why angles outweigh lengths (1.5 against 1.0), and the test that proves it earns the weight.**
A wrong camera constant scales every length and leaves every angle untouched, so an angular
disagreement is evidence about the *crystallography* while a length disagreement may only be
evidence about the *instrument*. Driving it: a 5 percent camera-constant error moves the length
deviation from 0.20 percent to 5.01 percent and the angle deviation not at all — 0.214 degrees
before and after — dropping the fused score from 0.983 to 0.794. Coverage is weighted highest
(2.0) because an unindexed spot is unexplained evidence that precision elsewhere cannot answer.

**Two judgements worth recording.** A solution with one indexed spot has *no* pair to measure an
angle between; that is missing evidence rather than disagreement, so the term is held neutral at
0.5 instead of scoring zero and punishing a solution for a spot the user did not pick. And the
agreement curve is polynomial rather than Gaussian, so two badly wrong solutions stay comparable
instead of both underflowing to zero.

31 tests. The invariances are the load-bearing ones: the angle term does not move when the
calibration does, and neither term moves when the pattern is rolled about the beam — which one
pattern cannot fix, so scoring must not pretend otherwise.

### Increment 3 result

- **`tem.fit_lattice`** exposes the fit: one row per pick with the lattice node it was assigned to
  and its residual, the refined centre, the overlay nodes clipped to the frame, and the notes. The
  centre of every practice plate is recovered to within 3 px from offsets of 0, 14 and 25 px, and a
  44 px mis-pick is named in the table rather than averaged away.
- **`tem.solve_pattern` now scores every candidate and ranks by the score**, not by the solver's own
  sort key. That key orders by matched fraction then residual and is explicitly not a quality; the
  score is one, and a list sorted by something other than the number printed beside it would be a
  trap. When the two orders disagree, a note says so and explains that the disagreement is itself
  evidence the pattern does not settle the answer.
- **Every candidate carries its calculated pattern** in picking coordinates, so accepting a solution
  can be a judgement made by looking. Pinned hard: for all three plates with scatter off, every
  simulated spot has a predicted node within 1e-6 px of it. Getting the scale or the rotation wrong
  would look like a disagreement the crystallography never had, so the test is exact rather than
  approximate.
- The per-spot table gained a Δd column in percent. A wrong camera constant shows as *the same*
  deviation on every spot — driving a 5 percent error gives exactly +5.00 on all six — which is the
  signature that separates a calibration error from an indexing error, and a test pins it.
- Scoring weights and both tolerances are exposed as advanced parameters in a Scoring group, and the
  policy travels in `data.score.weights` and in `inputs`.

**One sign error found by driving it.** The library test scales the measured `g` up, which makes
measured `d` smaller; the app test scales the *camera constant* up, which makes measured `d`
larger. The first draft of the app test asserted the library's sign and failed. The code was right
both times; the expectation was not.

28 new app tests.

### Increment 4 result

The picking canvas now carries four layers, each answering a different question, each toggleable.

- **The fitted lattice**, refitted on every pick, nudge and auto-pick, debounced at 200 ms because
  picking is a burst of clicks and a request per click asks a question only the last one wants
  answered. Drawn as two families of ruled lines rather than dots — a grid of points is a second
  set of spots to confuse with the first.
- **The two basis vectors, as labelled arrows** from the beam to the picks that generate them
  (added at the user's request mid-increment). The grid shows *that* the picks are consistent; the
  arrows show *which two* are carrying the whole lattice, which is what a user needs while
  adjusting. Each arrow ends on the picked spot rather than on the ideal node, so the gap between
  head and node is the error in that pick, visible without reading anything. `tem.fit_lattice`
  returns `basis_vectors` with the pick each one lands on, and reports honestly when no pick sits
  on a unit node.
- **The beam-centre tool**: a readout of where the beam is and where the spots say it should be, a
  directional pad with a step size, "Refine from the spots", and "Undo refinement". A pad rather
  than two number fields because the judgement is visual — nudge, look, nudge again — and a number
  box breaks that loop by demanding a value before showing its effect.
- **The calculated pattern of the selected candidate**, as open rings wide enough for the measured
  spot to show through the middle.

Candidates are listed ranked by score with three bars — d, angle, spots — so "why is this one
lower" is answered without a click. Selecting one draws its calculated pattern; **accepting** is a
separate deliberate act that carries the phase and axis into the atlas and tilt panels. Looking
costs nothing and commits to nothing.

**Two defects found by rendering the live SVG and looking at it.**

12. **The calculated rings were drawn smaller than the spots they explain**, so they sat *inside*
    the bright core and read as part of the spot rather than as a prediction about it — the one
    thing a superimposed pattern must never do. They are now wider than the core.
13. **The overlays used `var(--teal)` and `var(--violet)`.** The plate is always dark whatever the
    interface theme, and those tokens resolve to deep saturated values in light mode — a dark line
    on a near-black ground for half the users. Overlay colours are now fixed bright values drawn
    over a dark halo, which is legible on a dark plate and on a light-ground micrograph both.

Also: the arrowheads at `marker * 1.5` were invisible against the spot and are now `2.6`.

Verified by rebuilding the page's own layers from the services and rasterising them: with the beam
displaced 22, -16 px and one pick moved 46, 30 px, the grid still passes through the true spots,
the arrows still point at the two generating picks, and the mis-picked spot sits visibly off the
lattice between two rows.

### Increment 5 result

Export already covered every table-bearing result in CSV, XLSX and JSON. The gap was a genuinely
*human-readable* document: CSV is a grid with no provenance, JSON is complete and unreadable, and a
workbook is those two in separate sheets. None is the thing to paste into a notebook entry.

- **`result_to_markdown`** writes a report in the order a reader needs it: the answer in prose, the
  caveats, the data, the exact inputs that produced it, and the citations. It works on results with
  no table at all, because the prose and the provenance are the point.
- **`EXPORT_FORMATS` is now published in the manifest**, and the browser builds its export buttons
  from it. The three formats were hard-coded in `result.js`; a fourth would have needed an edit in
  Python and an edit in JavaScript, which is exactly the drift the manifest exists to prevent. A
  format added in Python now appears on every result in every panel at once.
- Verified in the running application: all seven panels offer CSV, Excel, JSON and Report, and
  driving the Report button produces `aluminium-fcc-down-001.md` with the title, the summary and
  the notes intact.

### Increment 6 result

- **Theory note** `theory/lattice_fit_and_solution_scoring` carries the derivations: why the centre
  is over-determined and why solving for it matters, the three failure modes of the fit and what
  answers each, the evidence criterion with its `log(A / πt²)` per inlier, Gauss reduction, and the
  scoring policy with the calibration argument that sets the weights. Indexed in
  `theory/index` and `docs/README`.
- **Workflow page** gained the fit-and-settle-the-centre step, the reading-the-candidates section,
  and an export section; the user guide gained the same in the TEM walkthrough and a rewritten
  export paragraph naming the fourth format.
- **Two more executable worked examples**, both analytic: a beam centre displaced 30 px in each
  direction is recovered exactly from eight lattice nodes, and measured vectors stretched by 5
  percent report a relative length deviation of exactly `1/1.05 - 1` on every spot.

**One more defect, found by driving the hcp plate end to end.** The accepted-solution message read
"down [010]" while the card above it read "down [1̄21̄0]" — the same axis in two notations on one
screen, because the candidate list used `PatternSolution.zone_axis_label` (the solver's own
three-index rendering) instead of the phase-aware `direction_label` everything else goes through.

**Final verification.** `ruff check .`, `mypy src` (138 files), the full `pytest tests/unit` lane
and `scripts/check_repo_integrity.py` are green. `python -m sphinx -b html docs/site` exits 0 with
the new theory page rendering; the seven non-autodoc warnings are the pre-existing docstring
formatting in `cbed`, `holz`, `models` and `plotting.runtime`. Driving the running application: the
hcp plate opened, auto-picked, the beam displaced 22 px and refined back to 0.0, indexed with the
correct-axis verdict, a candidate accepted, and the zone-axis list produced from the accepted
orientation. All seven panels export CSV, Excel, JSON and Report, and the Report button produces a
readable Markdown file.

### What was deliberately not done

- **No automatic spot detection.** Auto-pick works from the simulated truth, not from image
  analysis; finding spots in a real micrograph is local maxima, background estimation and a beam
  stop, and pretending otherwise would be worse than not offering it.
- **No intensity term in the score.** Relative intensities in a real pattern are dynamical and vary
  with thickness and tilt, so a solution scored on them would be scored on the specimen rather than
  on the crystallography.
- **The overlay is bounded by the index limit**, so a plate spot with no calculated ring beside it
  means the limit is too low, not that the solution is wrong. Stated in the help text rather than
  worked around.
- **The report export is Markdown, not PDF.** PDF would need a rendering dependency the zero-build
  deployment rule excludes; Markdown is readable as plain text and converts anywhere.

### Next task

None claimed. This goal is complete.

## TEM Module: Practice SAED Gallery And Zone-Axis Navigation — COMPLETE (2026-08-14)

**Objective.** Make the TEM panel of both shells usable for a researcher's day-to-day workflow
without needing a micrograph in hand. Ship a gallery of synthetic but scientifically exact SAED
patterns — fcc [001], hcp [2̄110]-family, bcc [110] — that can be picked and indexed exactly like a
real plate; add the infrastructure a real pattern needs (calibration from a known reflection,
auto-detected picks, an answer-check against the simulated truth); and add a zone-axis navigation
atlas that answers "what else can I reach from here, and how" rather than only "can I reach the one
axis I typed".

**Architecture.** Science stays in the library (`pytex.tem`, `pytex.diffraction`); the app layer
holds only the curated catalogue and the JSON contract; the browser holds presentation only. No new
JavaScript dependency, no build step. The synthetic pattern is transmitted as coordinates plus
intensities, not as a raster, so the same picking code serves both a real uploaded image and a
gallery entry.

### Initial audit (worktree clean on `main` at `ff3cfcc`)

- The TEM panel has two operations: `tem.solve_pattern` (pick → index) and `tem.plan_tilt`
  (orientation → alpha/beta). Both are sound and well documented.
- **The panel is unusable without a file.** Every one of its four examples is a tilt plan; the
  picking canvas shows a placeholder until the user uploads an image, so the indexing half of the
  panel — the half the microscope actually starts with — cannot be tried at all.
- `pytex.diffraction.saed.generate_saed_pattern` already produces exact zone-axis spot geometry in
  detector millimetres. Nothing converts that into the pixel-coordinate, finite-spot-size,
  noise-bearing thing a user clicks on, and nothing renders it.
- Tilt planning requires the user to already know which target axis they want. There is no surface
  that enumerates the low-index axes of a phase, their angles from the current one, and which of
  them the holder can reach — which is the question a microscopist actually has at the column.

### Planned increments

| # | Increment | Status | Commit |
| --- | --- | --- | --- |
| 1 | `pytex.tem.synthetic`: detector-pixel synthetic SAED patterns, with tests | done | (this commit) |
| 2 | `pytex.tem.atlas`: symmetry-reduced zone-axis atlas with angles from a current axis | done | (this commit) |
| 3 | App: gallery catalogue, `tem.gallery_pattern`, `tem.zone_axis_atlas`, calibration helper | done | (this commit) |
| 4 | Frontend: gallery, synthetic rendering, auto-pick, answer check, atlas navigation | done | (this commit) |
| 5 | Docs: workflow page, theory cross-links, user guide, worked examples | done | (this commit) |
| 6 | Human-style browser + desktop pass, defect fixes, closeout | done | (this commit) |

### Increments 1 and 2 result

- `pytex.tem.synthetic` projects a simulated zone-axis pattern onto a stated detector raster and
  returns pixel coordinates, relative intensities and display radii — the last step between
  `generate_saed_pattern` and something a user can click. `DetectorRaster` allows an off-centre
  beam, because a real plate has one and a workflow that assumes the middle of the frame has not
  been tested. Deterministic seeded jitter emulates centroiding error so residuals are realistic
  rather than machine-epsilon, and `independent_seed_spots` skips Friedel pairs, which are
  collinear through the beam and cannot seed an index.
- `pytex.tem.atlas` enumerates the symmetry-distinct zone-axis families of a phase with the angle
  from the axis currently on the beam, the family size, the reflection count inside a stated
  cut-off, and the pattern's apparent n-fold symmetry — the last measured on the simulated spot
  set rather than deduced from the point group, so it reports what the operator will actually
  recognise on the screen.
- **Convention fixed and documented:** the raster is the recorded image, with no handedness flip,
  so a picked `(column, row)` is the detector `(X, Y)` scaled by the pixel pitch. The round-trip
  tests hold the construction and `solve_saed_pattern` to that one convention.
- 44 tests. Radii are checked against `r = (camera constant)/d`; angles against the cubic closed
  forms (45°, 54.7356°) and the hcp basal-to-prism 90°; family sizes against the point-group orbit
  sizes (3, 6, 4); and every gallery pattern is round-tripped through the real indexer at three
  different rolls about the beam.
- **Two assumptions were wrong and the tests now say so.** A roll about the beam does *not*
  preserve the set of visible reflections, because a square frame clips differently as the pattern
  turns — which is what a real plate does. And an indexed pattern fixes the zone axis only up to
  symmetry, so a bcc [110] pattern legitimately indexes as [101]; the round-trip check compares
  families, not index triples.
- `ruff check .`, `mypy src` (135 files) and the new test module are green.

### Increment 3 result

- `pytex.app.tem_gallery` holds the three curated plates: aluminium fcc [001], ferrite bcc [110],
  and zirconium hcp [2̄110]. Each carries the lesson it teaches, three suggested next axes with
  reasons, and its own instrument setting — 200 kV, 400 mm camera length — from which the camera
  constant is *computed* as `L·λ`, never typed, so the field a user reads and the geometry the
  pattern used cannot drift apart.
- `tem.gallery_pattern` returns the pattern, the answer, the calibration and a set of suggested
  picks in one payload, so the browser never has to transcribe a calibration between panels.
- `tem.zone_axis_atlas` answers the question the panel could not: not "can I reach the axis I
  typed" but "which axis should I name". Reachability is computed by the same planner
  `tem.plan_tilt` uses, against the same envelope and roll, and a test pins the two to name the
  same destination and the same Δα.
- `_family_label` in the calculator service was promoted to a public `family_label`, since two
  services now need ⟨uvw⟩ and {hkl} rendering and reaching into a private name is how two
  renderings of one convention start.
- 23 new app tests. The important one is the end-to-end round trip: every gallery entry is opened,
  its suggested picks are handed to the real `tem.solve_pattern` with the calibration it reported,
  and the indexed axis must be the one the entry was built from.

**Three defects found and fixed while driving the new surfaces.**

1. **⟨110⟩ vanished from a 45° search.** The angle is computed as an arccos through a basis
   product, which lands a few ulps either side of exactly 45°, so a bare `>` comparison dropped
   the single most-wanted target about half the time. The filter now carries a tolerance.
2. **The axis already on the beam was offered as somewhere to tilt to.** The same arccos is
   square-root-behaved near 1, turning 1e-16 of cosine error into ~1e-6 of a degree — enough to
   pass a `> 1e-6` test. For ferrite [110] the summary read "the nearest reachable one is [110] at
   0.00°". The threshold is now a thousandth of a degree, which no distinct family can be inside.
3. **The default index limit buried the useful axes.** At `max_index = 3` a dozen high-index
   families with six-spot patterns crowded ⟨111⟩ — the six-fold axis anyone would actually want —
   off the end of a twelve-row table. The default is now 2, which is exactly the set a standard
   stereogram labels, and the help text explains what raising it admits and what it costs.

### Increments 4 and 5 result

- The TEM panel is now laid out as the session is: open a pattern, calibrate and index, choose
  where to go next, tilt. The gallery cards are read from the manifest — from the `pattern`
  parameter's options — so a fourth practice plate appears in the browser the moment it is added
  in Python, with no edit to the panel.
- A simulated plate is drawn as SVG from coordinates and brightnesses, a few kilobytes rather than
  a raster, crisp at any zoom, on a dark ground in both themes because that is what a diffraction
  pattern is. Auto-pick places the beam and six mutually non-collinear strong reflections; Show
  answer labels every spot; and the indexed result is checked against the construction, with the
  symmetry comparison done in Python where the symmetry lives.
- Atlas rows are actionable rather than transcribable: choosing one sets the tilt target below.
  Reachability is carried by a border stripe *and* by words, never by colour alone.
- Docs: a new workflow page `workflows/tem_pattern_indexing`, a step-by-step TEM section in the
  user guide, cross-links from `saed_pattern_solving`, and three executable worked examples whose
  expected values are analytic — `L·λ/d` for the calibration identity, `√3·a/c` for the hcp
  prism-zone aspect ratio, and the exactly-90° basal-to-prism angle.

**Four more defects found by driving the running application.**

4. **The transmitted beam was indistinguishable from a strong reflection.** Rendering the live SVG
   showed the direct beam as a 9 px core against 7.4 px for a 200 spot — while the panel's first
   instruction is "click the transmitted beam". It is now sized against the *nearest* reflection,
   so it is unmistakable without swallowing the inner spots of a dense hexagonal zone.
5. **The plate had no scale.** A diffraction pattern without one is a picture. A reciprocal-space
   bar in Å⁻¹, chosen from a 1-2-5 sequence, now makes the calibration visible: change the camera
   length and the bar changes with the pattern.
6. **A catalogue phase carried into the panel was renamed "(edited)".** The gallery sent an
   expanded `PhaseSpec`, which the phase picker reads as a user-edited phase; the indexed result
   was titled "Aluminium (fcc) (edited)" on a phase nobody had touched.
7. **A tilt plan counted "members of [012]".** A specific direction has no members. The sentence
   now writes the family form ⟨012⟩ per the notation standard, with a test pinning it.

**Also removed:** `GalleryEntry.beam_energy_kev` and `camera_length_mm` were dead — the operation
always used the request values — while the module docstring claimed each entry carried "its own
instrument setting". The voltage and camera length belong to the microscope, not the specimen, and
the docstring now says so.

**Verification.** `ruff check .`, `mypy src` (136 files) and the full `pytest tests/unit` lane are
green. The class-model atlas figures and the count stated on its page were regenerated, since two
new modules changed the model (272 public classes, 255 dataclasses, still only 6 inheritance
relations). Browser driving covered: gallery load, auto-pick, index, the correct-answer verdict,
the atlas, choosing a destination, and the resulting tilt plan; zero console messages; no
page-level horizontal overflow at 375 px or 1280 px; and the reachable/unreachable stripes checked
in dark mode.

**Known limitation of this session's verification.** The Browser pane would not composite frames,
so `screenshot` was unavailable throughout. Appearance was checked by extracting the live SVG's
own circle, line and text attributes from the DOM and rasterising *those* — a view of the real
output, not of a parallel implementation — which is how defects 4 and 5 were found.

### Increment 6 result

**Driven end to end in the browser**, repeatedly and from a cold reload: choose a plate, auto-pick,
index, read the verdict, list the zone axes, choose a destination, plan the tilt. Also: manual
picking at real screen coordinates (the beam plus three spots, indexed correctly), uploading an
image of one's own, the Show answer overlay, undo and clear, the example menu, the command palette,
the help drawer, 375 px and 1280 px layout, and dark mode. Zero console messages throughout; no
page-level horizontal overflow at either width (the only elements wider than the viewport are
inside `.table-wrap`, which scrolls).

**Driven in the real desktop shell.** `python -m pytex.app desktop` opens a native `pywebview`
window over a loopback server, reports `native_window: true`, serves the byte-identical frontend,
and answers `tem.gallery_pattern` through the same route. The window's *contents* could not be
scripted, so the panel itself was exercised through the browser shell, which loads the same files.

**A programmatic sweep stands in for the clicks nobody has time to make.** Every gallery plate at
three camera lengths, indexed through the real operation and checked against its construction; then
every atlas row at three index depths pushed back through the planner, which is exactly the pairing
a user makes by clicking a destination. 0 failures — after the sweep found the one below.

**Four further defects, all found by driving rather than by reading.**

8. **Pressing "Index the pattern" with nothing picked did nothing at all.** The error is real and
   well worded, but its field is the spot picker, which the panel hides because the value comes
   from the canvas — so the message landed on an invisible row. Errors on hidden fields now go to
   the toast and the plot status.
9. **A stale "Correct — that is the axis" card stood beside a failed index.** A verdict answers one
   attempt; a failed attempt now removes it, while the rest of the previous result stays.
10. **Loading your own micrograph kept the practice plate's calibration silently.** A camera
    constant from another exposure is the one error this panel cannot detect — it indexes to a
    plausible, self-consistent, wrong material rather than failing — so the plot now says so.
11. **`tem.plan_tilt` returned a 500 on a legitimate high-index hexagonal target.**
    `TiltSolution.orbit_member_indices` is `None` when the member has no low-index integer form
    within the navigation module's bound; both the tilt panel and the new atlas read it as an
    integer array regardless. Pre-existing, and reachable from the shipped panel by typing
    `[4 -3 1]` for zirconium. The row now reports the family form, which is what is known.

**Two loose scientific claims tightened.** The bcc entry named "110-type" reflections, but the zone
law `h + k = 0` excludes 110 from the [110] zone entirely — what is present is 11̄0. And the fcc
entry called 220 "the next ring in" when it is further out than 200. All three entries' claims are
now *measured on the plates the gallery produces*: equal lengths at 90° with a √2 diagonal at 45°
for fcc, a perpendicular √2 rectangle for bcc, √3·a/c for hcp, and 0001 absent while 0002 is
present.

### Verification

- `ruff check .`, `mypy src` (136 files) and the full `pytest tests/unit` lane are green.
- `python -m sphinx -b html docs/site` exits 0; `workflows/tem_pattern_indexing` and
  `examples/generated/saed_practice_patterns` both render. The seven non-autodoc warnings are
  pre-existing docstring formatting in `cbed`, `holz`, `models` and `plotting.runtime`; none comes
  from the new modules or pages. The remaining 595 are the known autodoc duplicate-object noise on
  `api/full_reference`.
- Three new worked examples recompute on every run against analytic values: `L·λ/d` for the
  camera-constant identity, `√3·a/c` for the hcp prism-zone aspect ratio, and the exactly-90°
  basal-to-prism angle.
- Class-model atlas figures and the count on its page were regenerated: two new modules moved it to
  272 public classes and 255 dataclasses, with still only 6 inheritance relations.
- `pytest --cov=pytex --cov-fail-under=87` passes at **91% overall**, with the new surfaces at 93%
  (`tem/synthetic.py`), 95% (`tem/atlas.py`), 97% (`app/tem_gallery.py`) and 88%
  (`app/services/tem.py`). Coverage did not decrease.

### What was deliberately not done

- **No dynamical intensities.** The plates are kinematic, so relative brightness is indicative and
  double diffraction is absent — a forbidden reflection a real plate shows through double
  diffraction is missing here. `pytex.diffraction.kinematic` can add and flag those; wiring it into
  the gallery would change what the plates teach and belongs in its own change.
- **No automatic spot detection on an uploaded micrograph.** Auto-pick works from the simulated
  truth, not from image analysis. Detecting spots in a real plate is a different problem — local
  maxima, background, a beam stop — and pretending otherwise would be worse than not offering it.
- **No fourth practice plate.** The gallery is read from the manifest, so adding one is a Python
  edit with no frontend change; three cover fcc, bcc and hcp, which is what was asked for.
- **The atlas takes about two seconds** for twelve rows, because reachability is computed by the
  real planner rather than a cheaper approximation. That was the right trade — a row marked
  reachable here is reachable there, and a test pins the two to agree — but it is the obvious place
  to optimise if the row count grows.

### Next task

None claimed. This goal is complete.

## Workbench Visual Modernization And GUI Completeness — COMPLETE (2026-08-14)

**Objective.** Improve the shared desktop/web workbench with a professional, responsive modern
interface; manifest-backed inline help linked to the Sphinx documentation; user-editable
visualization styling (beginning with diffraction spot shape, scale and colour); canonical
ready-to-run data for every GUI feature; a scientifically explicit Burgers beta-Zr (bcc) to
alpha-Zr (hcp) example; and human-style end-to-end testing of both shells with every discovered
defect fixed.

**Architecture kept.** Both shells continue to use the byte-identical frontend and the existing
JSON-in/JSON-out service layer in `pytex.app`. Scientific calculations remain in Python; browser
code controls presentation and camera state only. The zero-build, zero-third-party-JavaScript
deployment rule remains in force.

### Initial audit

- Worktree is clean on `main` at `9a916a9`, matching `origin/main`.
- Six panels exist (calculator, crystal, TEM, diffraction, variants, texture), each already has
  at least three runnable manifest examples and tests execute every example.
- Existing responsiveness tests cover 390, 768 and 1440 px widths; prior work fixed hidden tabs,
  overflowing plot toolbars, desktop export, drag lifetime, opening defaults and legend focus.
- Diffraction spots are still hard-coded as circles with one fourth-root radius mapping and
  generated variant colours; users cannot select marker shape, scale, palette or parent colour.
- Operation help is manifest-backed but has no first-class link to the matching Sphinx page.
- Burgers examples currently use the Fe-bcc phase as a symmetry-compatible stand-in for beta Zr,
  while labelling the case as zirconium. A real built-in beta-Zr phase is therefore required.

### Planned increments

| # | Increment | Status | Commit |
| --- | --- | --- | --- |
| 1 | Shared visual design refresh, explicit responsive/theme/accessibility checks | done | (this commit) |
| 2 | Manifest-backed documentation links and richer inline help | done | (this commit) |
| 3 | Reusable visualization-style controls; diffraction spot shape/scale/colour first | done | (this commit) |
| 4 | Canonical example audit and beta-Zr → alpha-Zr Burgers scenario | done | (this commit) |
| 5 | Full powder-XRD workbench with rigorous service, plot, help and examples | done | (this commit) |
| 6 | Browser + real desktop human-style feature pass and defect fixes | done | (this commit) |
| 7 | User guide/architecture/testing docs, full verification, closeout | done | (this commit) |

### Increment 1 result

- Reworked the shared visual tokens around the canonical PyTex ink/blue/violet language, with
  softer hierarchy, modern card/control geometry, disciplined shadows and a subtle stage glow.
- Added a three-state Auto/Light/Dark theme control. The choice persists in local storage, Auto
  continues to follow the operating system, and the button retains a useful accessible name when
  its visible label collapses on a phone.
- Repaired the crystal toolbar's desktop composition: the global full-width `select` rule had
  made the figure-format picker consume a complete flex row, stacking a six-control toolbar into
  three lines at 1440 px. Toolbar selects are now intrinsically sized while the toolbar still
  wraps at narrow widths.
- At 390×844, all six tabs and all three masthead actions are visible, the actions collapse to
  icons with titles/ARIA labels, and document width equals viewport width. At 768×1024 the rail
  moves below the stage and both are exactly viewport width; at 1440×900 the rail remains beside
  the plot and the toolbar stays on one line.
- Playwright verified theme cycling, persistence across reload, light and dark rendering, and a
  clean browser console. Source-level tests pin the three theme states, shared-shell placement,
  persistence, and the toolbar width exception. `tests/unit/test_app_server.py` (35 tests) and
  targeted Ruff checks are green.

### Increment 2 result

- Every operation now carries a typed `DocumentationLink` in the manifest, selected by panel from
  the closest canonical Sphinx concept, workflow or theory page.
- Manifest tests prove that every target exists below `docs/site/`, rejects traversal, and checks
  the emitted source URL. A renamed or deleted page now fails the base lane instead of leaving a
  silently stale Help button.
- The shared help drawer renders a prominent Sphinx documentation card before its input reference,
  so concise inline guidance remains close to the control and the deeper scientific explanation is
  one deliberate click away. The target opens separately and cannot replace a user's working GUI
  state.
- Because PyTex has no published Sphinx host yet, the link opens GitHub's rendered view of the
  canonical MyST source. This preserves the single-source documentation doctrine; the architecture
  note records how an intranet deployment may mirror the same source.
- The user guide now documents the field help, operation drawer, documentation link, theme states,
  and narrow-screen accessible names. Browser driving verified the new card in dark mode and a
  clean console. The combined manifest/server target (306 tests) is green.

**Exact next action.** Build a reusable presentation-control surface and apply it to composite
SAED: marker shape, spot-size scale, intensity sizing mode, parent colour and variant palette,
with immediate redraw from one simulation and no scientific recomputation.

### Increment 3 result

- Added a reusable, dependency-free marker-style module with circle, square, triangle, diamond,
  star and cross shapes; filled and unfilled rendering; a 0.5–2.5× visual scale; perceptual,
  intensity-proportional-area and constant sizing;
  editable parent/product colours; distinct, colourblind-safe and single-product palettes; and a
  complete reset action.
- Applied it to composite diffraction without changing the service contract. Appearance changes
  redraw the existing rows immediately, while coordinates, indices, intensity values, hover
  records and exports remain unchanged. Legend swatches follow palette, marker shape and fill.
  Double-diffraction reflections remain dashed when the global display is unfilled.
- Made composite visibility deliberate rather than implicit: individual source chips still toggle
  one lattice and preserve keyboard focus, while Show all, Parent only and Focus a variant provide
  one-action recovery, decluttering and parent-plus-variant isolation for a 24-variant pattern.
- Added concise in-place guidance explaining the scientific/display boundary and the trade-off of
  each sizing and palette choice. The workflow guide now documents those controls and the
  architecture note pins presentation choices to the frontend side of the service boundary.
- Human browser driving changed circle markers to diamonds, changed to the colourblind-safe cycle,
  switched intensity sizing, and reset every option. It isolated Variant 7 plus the parent (32 of
  192 spots), reduced to the parent alone (24), restored all sources, then hid Variant 3 alone
  while focus remained on its chip. Layout had zero horizontal overflow and the console stayed
  clean. Targeted JavaScript syntax, Ruff and 59 server/diffraction tests are green.

**Exact next action.** Add a first-class beta-zirconium bcc phase with a cited high-temperature
lattice parameter, replace the Fe-bcc stand-in in each Burgers example, and pin the 12-variant
beta-Zr to alpha-Zr case in manifest and GUI tests.

### Increment 4 result

- Added `zr_bcc_beta` as a first-class catalogue phase: monatomic bcc Zr, `Im-3m` (229), two-Zr
  conventional basis, and `a = 3.6090 Å` at 863 °C. The display name includes the temperature so a
  high-temperature allotrope is not presented as a room-temperature constant. Provenance names the
  Zuzek phase assessment and the IUCr journal table with DOI `10.1107/S1600577515009054`.
- Replaced the Fe-bcc symmetry stand-in in all three Burgers GUI examples. Calculator,
  Diffraction and Variants now share the explicit `zr_bcc_beta` parent and `zr_hcp` child; manifest
  tests pin that exact pair and catalogue tests pin the cell, space group and Zr basis.
- The focused application lane is green (397 passed, one optional-stack skip). Human browser
  driving ran all three examples: Calculator reported the named zirconium pair and 12 variants;
  Diffraction rendered 426 spots from 12 variant sources; Variants rendered 12 basal poles in six
  packets of two. Both phase pickers showed the expected catalogue identifiers, all layouts had no
  horizontal overflow and the browser console remained clean.
- The architecture record and workbench guide now distinguish beta-bcc and alpha-hcp zirconium,
  document the temperature-specific cell parameter, and explain that the examples exercise real Zr
  spacings, scattering and provenance rather than iron with compatible symmetry.

**Exact next action.** Perform the complete human-style web pass across every example, control,
export and responsive width, then read the computer-use confirmation/API guidance and repeat the
critical workflows in the real desktop shell. Fix every reproducible GUI defect before closeout.

### Increment 5 result

- Added a seventh shared workspace, **XRD**, backed directly by
  `pytex.diffraction.xrd.generate_xrd_pattern`. The service exposes the canonical phase catalogue,
  Cu/Mo/Co single-line or doublet radiation, tabulated or teaching scattering models, scan range,
  Gaussian/pseudo-Voigt profile, width, mixing, angular sampling and reflection-index limit.
- The result contract carries a normalized sampled diffractogram and one exportable/hoverable row
  per primary reflection family: canonical label, $2\theta$, $d$, relative integrated intensity,
  multiplicity, $|F|$ and Lorentz–polarization factor. Its prose explicitly bounds the result as
  kinematic, background-free and unsuitable for quantitative phase or Rietveld analysis.
- Added live presentation-only profile/stick colours, line width, area fill, reflection sticks,
  peak labels and threshold, plus linear, square-root and log-like display transforms. These redraw
  the existing arrays without a second service request. Four canonical examples cover the fcc
  sequence and resolved Cu doublet, silicon extinctions, Mo wavelength shift and alpha-Zr metrics.
- Scientific tests pin nickel 111 at 44.495° and the conventional 111/200/220/311 sequence, prove
  single/doublet and Gaussian/pseudo-Voigt distinctions, require shared hover/export columns, and
  exercise every example. Human browser driving verified the indexed chart, log transform, live
  colour/width/fill/stick changes, dark-theme presentation and a clean console.
- Extended composite-SAED appearance with triangle and star symbols plus a global filled/unfilled
  mode. Variant identity can use any selectable combination of shape, size and colour; the default
  three-channel encoding gives all 24 product variants distinct symbols and the interactive legend
  mirrors them. Unfilled double-diffraction spots retain their dashed edge, while the transmitted
  beam has an independent double-ring marker and explicit `(000)` plot/legend labels. Human browser
  driving rendered a 192-spot pattern as enlarged outline stars, changed its parent colour live,
  and found no console warning or error.

### Increments 6–7 result and closeout

- Playwright exercised the final shared frontend at desktop and 390-pixel phone widths. The XRD
  workspace produced an indexed nickel pattern, responded live to log-like scale, colour, width,
  fill and stick controls, retained all seven tabs without horizontal overflow, and emitted no
  console warning or error. Composite SAED proved 24 distinct shape/size signatures even under a
  single colour, 24 distinct full signatures in the default mode, all seven selectable channel
  combinations, and a real SVG `(000)` transmitted label.
- The restarted native desktop shell loaded the same final assets. It rendered the XRD chart and
  log transform, filled and unfilled star spots, the default 24-way shape/size/colour pattern, the
  independent double-ring transmitted marker and the complete variant-encoding control. No
  desktop-only interaction or layout defect was found.
- Final verification is green: JavaScript syntax checks; focused server, diffraction and XRD tests;
  Ruff over `src` and `tests`; mypy over all 133 source files; Sphinx HTML with warnings as errors;
  and the complete unit suite (one expected optional-stack skip). The unit lane took 722.5 seconds,
  including live worked-example and teaching-notebook checks.
- The objective is complete. The shared workbench now has seven responsive, documented workspaces;
  canonical runnable examples for every manifest operation; explicit beta-Zr to alpha-Zr Burgers
  examples; presentation-only diffraction decluttering and styling; and a rigorous powder-XRD
  surface with clearly stated scientific limitations.

**Exact next action.** None for this objective. Resume future application work from a new ledger
entry after reviewing this closeout and the commits listed above.

## Application Platform: Desktop + Intranet Workbench — COMPLETE (2026-08-12)

**Objective.** Build one interactive application over the PyTex library that ships both as a
desktop app and as an intranet web app, sharing all scientific code and the entire user interface,
diverging only in how a window opens and where files are written. Flagship surfaces: an
interactive 3D crystal viewer with arbitrary superimposed planes/directions/annotations; an
interactive TEM diffraction-pattern solver (upload → calibrate → pick spots → index → plan the
tilt to the next zone axis); and a crystallographic calculator for interplanar angles, symmetric
families, and cross-phase geometry. Every result exports as a publication figure *and* as
re-plottable numbers (CSV/XLSX/JSON).

**Design record.** `docs/architecture/application_platform.md` — read it first. It fixes six
decisions: (1) the shared layer is a JSON-in/JSON-out **service layer**, not a widget library;
(2) the UI is **generated from a self-describing operation manifest**, so help text cannot drift;
(3) **zero mandatory runtime dependencies** — stdlib `http.server`, hand-written ES modules, no
bundler, no third-party JS, because the deployment target is an offline intranet host;
(4) interactive 3D **projects in the browser, publishes from Python** (scene built once in Python,
camera lives in the browser, export replays the camera through `pytex.plotting`);
(5) the desktop shell is the web shell in a window (`pywebview` if present, else default browser);
(6) every result is exportable in a re-plottable form; (8) **every plot carries a live cursor
readout in its own physical units and hover detail on every drawn entity** — a spot shows its hkl,
d, |g|, intensity, phase and variant, taken from the same row the CSV export writes; (7) **every tab ships with 3–4 canonical
worked examples** so a user with no data of their own can still exercise every feature — the shared
material set is NaCl, Fe-fcc (austenite), Fe-bcc (ferrite), and Zr-hcp, all present in the built-in
phase catalogue with cited parameters and full atomic bases.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 1 | Architecture and decision record | done | (this commit) |
| 2 | `pytex.app` foundation: errors, JSON envelope, operation manifest/registry | done | (this commit) |
| 3 | Phase specification + built-in phase catalog (no pymatgen required) | done | (this commit) |
| 4 | Calculator service (`calc.*`) + tests | done | (this commit) |
| 4b | Example-scenario registry + named-OR operation; 7 calculator examples | done | (this commit) |
| 5 | Stdlib HTTP server + routing + in-process loopback tests | done | (this commit) |
| 6 | Frontend shell: tabs, manifest-driven controls, palette, calculator tab | done | (this commit) |
| 7 | Crystal viewer service + panel (browser projection, Python export) | done | (this commit) |
| 8 | Export surface: PNG/SVG figures, CSV/XLSX/JSON data, stdlib xlsx writer | done | (this commit) |
| 9 | TEM solver service + panel (calibrate, pick, index, tilt plan) | done | (this commit) |
| 10 | Desktop shell (`pywebview`, browser fallback) + `python -m pytex.app` verbs | done | (this commit) |
| 11 | Diffraction tab: composite SAED of parent + product variants | done | (this commit) |
| 11d | Drive the running app in both shells; fix what only running it reveals | done | d81b26b, e1b627e, 2b45b7e |
| 11e | Second driving pass: every operation from its own defaults | done | (this commit) |
| 11b | Orientation-relationship visualization (variants, pole figures, packets) | done | (this commit) |
| 11c | Texture tab: pole figures, IPF, ODF | done | (this commit) |
| 12 | User guide, worked example, docs index wiring | done | (this commit) |

### Current worktree state

Steps 1–11 and 11d landed. Every panel carries at least three runnable examples; the manifest
test enforces it and executes each one, so an example cannot rot.

**Step 11d — what driving the application found.** The panels were tested by clicking them, in
the browser and in the real pywebview window, with a synthetic Ni [001] and [011] SAED pattern
generated to a known camera constant so the indexer could be checked against an answer known in
advance (it recovers both zone axes, 100% of spots, residual < 0.002 Å⁻¹). Eight defects that
no unit test was positioned to see:

1. **The crystal viewer's Figure button was dead** — `crystal.render` rejected three parameters
   `crystal.scene` accepts, so the request failed the moment a user touched the unit-cell
   outline, the atom labels, or the bond tolerance. Fixed by giving the renderer the parameters;
   a test now asserts every scene parameter is one the figure accepts.
2. **Mathtext leaked to the browser and to prose** — plane overlays were labelled `$(100)$` on
   screen, and the tilt planner told an operator to drive to `$[0\bar{1}1]$`. The repository
   convention was already plain style for `describe()`; `pytex.tem` had not followed it.
3. **The closed overlays covered the whole application** — `.palette, .drawer { display: flex }`
   outranks the user agent's `[hidden] { display: none }`, so both full-viewport layers sat open
   over the page and swallowed every click. Nothing in the JavaScript was wrong, which is why no
   panel test could see it.
4. **A drag turned the crystal once and stopped** — turning redraws the scene, a redraw builds a
   new SVG, and the drag handlers belonged to the SVG, so every movement after the first reached a
   node with no drag in progress. The handlers now live on the frame, which outlives the picture.
5. **Two-phase operations failed on their opening press** — neither the named orientation
   relationship nor the composite SAED declared a default phase, so the picker started parent and
   child on the same catalogue entry and the relationship refused them. Every operation is now run
   from its declared defaults by a test, substituting for an absent phase exactly what the picker
   sends.
6. **An edited phase kept the last one's identity** — editing nickel into a tetragonal cell left
   the name "Nickel (fcc)" on every title and export, and left Fm-3m applying an F centring that
   deleted two thirds of the reflection families. Edited phases are now named "(edited)", the
   space group clears when the crystal system changes, and `PhaseSpec` refuses a space group from
   another system outright.
7. **The desktop shell could export nothing** — an anchor with `download` is accepted and
   silently dropped by the embedded web view, so every export in the desktop shell produced no
   file and no error. Fixed with a `SaveBridge` js_api that writes through a native save dialog,
   `GET /api/shell` so the page is told which shell it is in rather than sniffing, and a
   `pytex:saved` toast in both shells so an export never again passes in silence.

8. **The viewer could not publish an SVG, and its PNG could not be saved on the desktop** — the
   Figure button hard-coded PNG, so the format the help explains when to choose was unreachable;
   and the PNG went out through a private copy of the anchor-download trick, so it failed silently
   in the desktop shell while the SVG beside it saved. A format control is on the toolbar, and
   every file now leaves through `saveBlob`, which a test enforces by refusing a second
   `createObjectURL` anywhere in the frontend.

Also: the page names its own favicon, so a load no longer ends in a 404 for `/favicon.ico` in the
log of a server that otherwise logs only real requests.

**Verification.** Both shells were driven to the end. In the browser: every panel, all ten
calculator operations, the palette, the help drawer, the examples, custom phases, the error paths,
and the exports in all three formats. In the real pywebview window, through `evaluate_js`: all
four panels compute, a drag turns the crystal, a CSV and a 481 kB PNG figure reach disk through
the save bridge with the toast naming the path, and the console is clean.

**Step 11e — a second driving pass, this time over opening states.** The application was driven
again in the browser, with every one of the ten calculator operations selected in turn and its
button pressed *without touching a single control*, which is what a user does first and what step
11d had not systematically covered. All ten produce a result and no console error. Two defects:

9. **A machine slug reached a heading.** The orientation-relationship result was titled
   `kurdjumov-sachs: Austenite (fcc Fe) to Ferrite (bcc Fe)` — a Python identifier with its
   underscore swapped, printed where two surnames belong. Three call sites reconstructed the
   display name from the identifier by hand, and none of them could reach the proper name the
   choice list already carried. A single `relationship_name()` now derives it from the choice
   labels, and titles, prose, and error messages all take it from there.
10. **The cross-phase angle operation opened on a question with no content.** Its defaults left
   the second phase unset, so the phase picker sent the same phase twice and the opening press
   compared nickel with nickel through a null rotation — a table of angles a single-phase
   operation already answers. The defaults are now austenite against ferrite under the first
   Kurdjumov-Sachs variant, so the first press produces the parallelism the help text promises:
   austenite (111) against ferrite (011) at 0°, with 45°, 54.74° and 80.41° beneath it. The
   rotation is stored as a literal because a manifest is static, and a test recomputes it from
   `OrientationRelationship` and fails if the two ever part company.

The general lesson, now a rule for new panels: **an operation's declared defaults are a user-facing
surface and must be exercised as one.** `test_an_operation_runs_from_its_own_defaults` proved they
*run*; neither it nor anything else asked whether what they produce is worth looking at.

**Step 11f — the layout at widths nobody had opened it at.** The page was measured in a browser at
390 × 844, 768 × 1024 and 1440 × 900, walking all four panels at each width and collecting every
element whose box crossed the viewport edge. Two defects, both invisible at the desktop width the
application had always been developed at:

11. **Three of the four workspaces were unreachable on a phone.** `.tabs` was a horizontal scroll
    container with `scrollbar-width: none`. At 390 px it was 65 px wide holding 395 px of tabs, so
    it showed *Crystal Viewer* and nothing else — no cut-off edge, no scrollbar, no hint that a
    sideways drag would reveal the rest. Navigation is the one thing that must never be hidden,
    and hiding a scrollbar hides the only sign that content is off-screen. The bar now wraps.
12. **The figure toolbar left the window.** `.plot__header` has visible overflow by design — the
    cursor readout and the detail popover sit outside the flow — so an over-wide toolbar gets no
    scrollbar and no clip: at 390 px the preset buttons, the format select and the export button
    rendered out to x = 497 on a 390 px screen. Header and toolbar now both wrap.

Also in this pass, deliberately rather than as a fix: the masthead reflows to two rows below 48 rem
with the tab bar full-width beneath the brand; controls take a 2.75 rem minimum target under
`(pointer: coarse)`, keyed off the pointer rather than the width so a narrow desktop window keeps
its compact controls; and `prefers-reduced-motion` is honoured.

After the change, at all three widths and on all four panels: no horizontal overflow anywhere, all
four tabs on screen, and result tables still scrolling inside their own `.table-wrap` box as
intended. Two source-level tests pin the wrap rules, because a stylesheet regression is exactly the
kind a Python suite cannot see.

### Step 11b — the variants panel

A fifth workspace, `pytex.app.services.variants` plus `js/panels/variants.js`, answering the two
questions a variant set is actually asked.

**`variants.pole_figure`** projects the chosen child plane family of every variant into the parent
frame, coloured by packet. The counts are the published ones and are pinned as such: Kurdjumov-Sachs
gives 24 variants in 4 packets of 6 on the parent {111}, Nishiyama-Wassermann 12 in 4 packets of 3,
Burgers 12 in 6 packets of 2 on the parent {110}, Bain 3.

**`variants.intervariant_misorientations`** gives every variant pair its disorientation, with the
nearest low-index axis and how far that label is from the exact one. The Kurdjumov-Sachs spectrum
comes out as the ten angles of Morito et al., Table 2 — 10.53, 14.88, 20.61, 21.06, 47.11, 49.47,
50.51, 51.73, 57.21, 60.00 — and the within-packet subset as the three exact ones, 10.53 about
⟨111⟩ and ⟨110⟩, 49.47, and 60 about ⟨111⟩ and ⟨011⟩. Nothing in the test file is compared against
a previous output of this code.

**`variants.render`** publishes the same poles through `pytex.plotting.spherical`, called *through*
the pole-figure handler rather than recomputing, so the figure and the CSV cannot describe different
poles.

Four things the implementation had to decide, recorded because none is obvious from the code:

- **Equal-area is not normalised by the library.** `project_directions` returns each projection at
  its natural radius, and equal-area reaches √2 where stereographic reaches 1. Drawing both in the
  same circle puts 41% of the equal-area figure outside the rim. The service divides each by its own
  equatorial radius, so the rim means 90° in both and the exported x and y are comparable.
- **The packet colours exist twice**, as hex for matplotlib and as HSL for CSS, because neither
  renderer can read the other's constant. A test converts and compares them; it caught a real
  disagreement (`#2f7fd4` on screen against `#2c7fdd` in the figure) the first time it ran.
- **A disorientation axis is a line, not a ray**, and the sign out of a symmetry reduction is
  arbitrary, so the label is canonicalised to a positive leading index — otherwise ⟨1 -1 1⟩ and
  ⟨-1 1 -1⟩ appear as two rows of the same table.
- **The axis index limit is 8, not 6.** At 6 the 14.88° pair was labelled ⟨2 5 0⟩, nearly four
  degrees away — a label worse than none, since the literature quotes these up to ⟨5 5 7⟩.

**Step 11b — what driving it found.** Two defects, both in code the tests passed over:

13. **Toggling a legend threw away the keyboard user's place.** Both plotting panels rebuilt the
    whole legend as part of the redraw a legend click causes, so the button just pressed left the
    document and the browser moved focus to `body`. Measured: focus went from "Packet 2 (6
    variants)" to BODY on one click. The legend is now built once per result and updated in place.
    The defect was older than this panel — the diffraction legend had it too, and was fixed with it.
14. **The Figure button disabled the wrong control**, because the format select carries `.button`
    for its styling and `.plot__toolbar .button` matches the select first. Held by name now. In the
    same place: the format list opened on PNG while the operation's own default and its help text
    both say SVG is the right artefact for line art. SVG is first now, and the three agree.

### Step 11c — the texture panel

A sixth workspace, `pytex.app.services.texture` plus `js/panels/texture.js`, showing one model
polycrystal the three ways texture is read: `texture.pole_figure`, `texture.inverse_pole_figure`
and `texture.odf_sections`.

**Why a model texture rather than a file import.** The application has no measurement to load, and
a texture panel that can only say "import an EBSD map first" teaches nothing. The decisive argument
is not convenience: a model built from the components the literature names has a *known answer*,
and a data set does not. A random texture must read 1 m.r.d. everywhere; Goss is {011}⟨100⟩, so its
(011) poles must sit at the centre of the figure. Those are the tests, and neither is available
from a data set. The panel offers cube, Goss, brass, copper, S, the standard fcc and bcc rolling
mixtures, and random, with controls for component spread, grain count, kernel halfwidth and seed.

**Two defects found by checking numbers against what they claimed to be.** Both are the same kind
— a quantity wearing a unit it was not in — and neither could have been caught by a test that
compared a figure with itself.

15. **A public library method crashed on two independently built phases.** `Phase` is a dataclass
    with array-valued fields, so `phase != other` raises `ValueError: The truth value of an array
    with more than one element is ambiguous`. Four guards in `pytex/texture/models.py` used `!=`,
    including `ODF.evaluate_pole_density`, which therefore failed for any caller who built its pole
    from its own copy of the phase rather than passing the identical object. `phases_semantically_match`
    already exists for exactly this and documents itself as the thing to prefer "whenever the two
    objects may have been constructed independently". Fixed at all four sites, with a regression
    test in `tests/unit/test_texture.py` that uses two separately built phases — with one shared
    object the bug is invisible — and a second test that the guard still refuses a genuinely
    different phase.
16. **The ODF sections were 24 times too large, in a column labelled "m.r.d."**
    `ODF.evaluate(normalized=True)` normalises the kernel over the whole of SO(3), but a
    symmetry-aware evaluation folds every query into the fundamental zone, which is 1/|G| of it. A
    random texture therefore read |G| rather than 1: measured at 23.9 for m-3m and 11.9 for 6/mmm,
    against proper-operator counts of 24 and 12. A sharp cube component came out at 737 m.r.d.,
    which is not a texture strength any material has. Dividing by the operator count puts it on the
    same scale as the pole figures beside it: cube now reads 19.6 m.r.d. and the fcc rolling texture
    4.9, which are the published ranges. Both symmetries are checked, because one cannot distinguish
    a correct normalisation from a constant fudge factor.

Also in this pass, from driving it: **the on-screen table is now a capped preview.** A texture ODF
is 1083 rows, and past a couple of hundred the scroll box has no search, no sort and no way to
reach row 900 except dragging, while the export buttons directly above it carry every row. The
first 200 are shown with the caption saying so — a table that quietly shows a subset is worse than
one that is too long — and the export still sends the whole result: measured at 200 rows on screen
and 863 lines in the CSV of an 862-row figure. DOM nodes on the ODF view fell from 6789 to 2360.

### Step 12 — documentation

- **`docs/site/workflows/workbench_application.md`** — the user guide: running both shells, what
  the six workspaces answer, how to read a result, the two kinds of export and when each format is
  right, scripting the service layer directly, and what degrades without `pywebview` or
  `matplotlib`. It ends with two panels a reader can check against answers fixed before the
  calculation runs.
- **`worked_examples/examples/workbench_service_layer.py`** — four executable examples over the
  service layer, rendered into `docs/site/examples/generated/workbench-service-layer.md`: the
  Kurdjumov-Sachs packet structure and the ten intervariant angles against Morito Table 2, the
  m.r.d. mean identity checked on a random, a single-component and a five-component texture at
  once, and the Goss (011)-on-ND claim that a Miller label makes.
- **Index wiring** — the guide is in the Sphinx workflows toctree; the application-platform design
  record now has a site stub (it had none, so `{doc}` links into it did not resolve) and is in the
  architecture toctree; `docs/README.md` points at the guide from the design record's entry.
- **`test_workbench_guide_quotes_numbers_the_code_actually_produces`** in
  `tests/unit/test_documentation_policy.py` — the guide's three quantitative claims are recomputed
  from the service layer and matched against the text, per the `AGENTS.md` rule that documentation
  numbers must not be hand-transcribed.

### Where the application stands

All twelve steps are done. Six workspaces, 21 operations, 30 runnable examples, both shells driven
end to end, and the layout measured at 390, 768 and 1440 px on every panel.

**Resume point.** No step is outstanding. The application platform task is **complete**; if it is
picked up again, the natural next increments are an EBSD-map import path for the texture panel
(which today builds a model texture only, deliberately — see step 11c), and a `variants.render`
equivalent for the texture figures, which currently publish from the browser SVG rather than
through `pytex.plotting`.

Note for whoever continues: every defect recorded above was invisible to a passing test suite and
obvious within a minute of using the application or of reading its numbers against what they
claimed to be. New panels are driven the same way before they are called done, and the driving now
covers four things in particular — the opening press with nothing touched, the layout at 390 px,
where focus goes after every control that redraws, and whether a reported quantity is really in the
unit its column header names.

### Decisions taken during implementation

- Phases entering the app are described by a **`PhaseSpec`** (six cell parameters + point group +
  optional space group and atomic sites), which is JSON-round-trippable and needs no optional
  dependency. CIF import stays an optional path through `Phase.from_cif`, so the app runs fully
  without pymatgen.
- The built-in phase catalog is defined by literal parameters in Python with citations, not by
  reading the CIF fixtures, so the app's starting catalog cannot be broken by an optional
  dependency being absent.

## Deepening The Theory Layer: Closed Forms, Known Answers, Non-Obvious Cases — COMPLETE (2026-08-11)

**Objective.** Increase the breadth and depth of the theory and algorithm notes across the repo.
For every surface treated: state the closed-form solution where one exists, show how the code
evaluates it, and demonstrate agreement with an independently known answer. Prioritise the
**non-obvious, non-textbook** decisions — the ones a reader cannot reconstruct from a textbook
because they are implementation conventions rather than physics.

**The governing test for what deserves a note.** Not "is this important?" but *"would a competent
crystallographer reading the code be unable to tell why it does that?"* IPF colouring is the
canonical example: every textbook shows the coloured triangle, none states that the colour is a
barycentric coordinate, that it is raised to a power to control saturation, or that it is then
rescaled so the largest channel saturates. Those three choices determine every pixel of an IPF map
and appear in no reference.

### Gap survey (2026-08-11)

Searched all 37 notes for each concept. Findings:

| Concept | Code | Notes mentioning | Verdict |
| --- | --- | ---: | --- |
| IPF colour keys | `plotting/ipf.py` (331) | 0 with the algorithm | **gap** — workflow page is 61 lines, zero mathematics |
| Taylor factor, Schmid | `properties/taylor.py` (141), `slip.py` (252) | **0** | **gap** |
| Elastic homogenisation, Voigt/Reuss/Hill | `properties/tensors.py` (689) | **0** | **gap** |
| Misorientation distribution (Mackenzie) | `core/misorientation_distribution.py` (248) | **0** | **gap** |
| Sphere sampling / equal-area grids | `core/sphere.py` (793) | partial | thin |
| Texture kernels | `texture/kernels.py` (484) | partial | thin |
| Texture components and fibres | `texture/components.py`, `fibres.py` | partial | thin |

`src/pytex/properties/` — roughly 1,080 lines covering slip systems, Schmid tensors, Taylor
factors, elastic tensors, Voigt/Reuss/Hill homogenisation, and directional modulus surfaces — has
**no theory note at all**, despite being the richest source of closed-form results and classic
published reference values in the repository.

### Why these four first

Each has an exact closed form *and* a published number to check against, so each note can end in a
verified agreement rather than an assertion:

1. **IPF colour keys** — named by the user; the barycentric/gamma/renormalisation chain is pure
   convention and is invisible in the literature.
2. **Mackenzie distribution** — an exact analytic density with a published mean disorientation
   (42.03° for random cubic), so the sampler can be checked against theory rather than itself.
3. **Elastic anisotropy** — closed-form directional Young's modulus for cubic symmetry, the
   Voigt/Reuss bounds that must bracket any real aggregate, and the Zener ratio.
4. **Taylor/Schmid** — maximum Schmid factor exactly 0.5 at 45°/45°, and the Taylor factor ≈3.06
   for a randomly textured fcc polycrystal (Taylor 1938).

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 1 | Gap survey across all 37 notes and the source tree | done | (this commit) |
| 2 | IPF colour keys note + worked examples | done | (this commit) |
| 3 | Misorientation distribution / Mackenzie note + worked examples | done | (this commit) |
| 4 | Elastic anisotropy and homogenisation note + worked examples | done | (this commit) |
| 5 | Taylor factor and Schmid analysis note + worked examples | done | (this commit) |
| 6 | Wire new notes into `theory/index.md`, `docs/README.md`, cross-links | done | 0f223d2 |
| 7 | Pole-figure arithmetic and the m.r.d. scale note + worked examples | done | 35809d9 |
| 8 | Directional statistics / mean axes note + worked examples | done | f3bd460 |
| 9 | Ghost problem / odd harmonics note + worked examples | done | (this commit) |

### Notes added so far

| Note | The non-obvious thing it fixes |
| --- | --- |
| `ipf_color_keys` | Colour is a barycentric coordinate raised to $1/\gamma_{s}$ and renormalized; every IPF colour is fully saturated, so brightness cannot carry a second variable. |
| `random_disorientation_baseline` | The random cubic mean is 40.7° and the *median* is 42.3°; 2.2% of random boundaries are low-angle by chance. |
| `elastic_anisotropy_and_homogenization` | Compliance carries Voigt factors that stiffness does not; cubic $K_V = K_R$ exactly, so the whole bound gap is in shear. |
| `schmid_and_taylor_plasticity` | Schmid $\le 1/2$ in two lines; the Taylor LP replaces a 792-way enumeration; $M$ is unique but the slip pattern is not. |
| `pole_figure_arithmetic_and_mrd` | An unweighted raster mean is biased by exactly +50% and refinement does not fix it; skipping m.r.d. normalization is a 34× error. |
| `directional_statistics_and_mean_axes` | The vector mean of axes cancels identically; the orientation tensor is the fix, and its eigenvalue gap decides whether a "mean axis" means anything. |
| `ghost_problem_and_odd_harmonics` | Pole figures determine only even harmonic degrees — ~47% of the basis is undetermined by *any* amount of data — and EBSD has no ghost problem at all. |

### Method that is working

Each note follows the same shape and it should be kept: state the closed form, show the code
agreeing with it to a stated precision, then draw out the consequence that constrains use. Every
expected value in a worked example is derived by hand or published — never a recorded program
output — and where no closed form exists (the cubic disorientation mean) an *independently written*
implementation is used instead, because a check sharing code with the thing it checks is not a
check.

Probing the code before writing has repeatedly found things worth documenting that were not
visible in the source: the 7.7 GiB allocation in the disorientation reduction, the exact
$[110] = [112]$ stiffness coincidence, and the fact that the barycentric sum-normalization in the
IPF path is arithmetically inert.

**Resume point.** Remaining candidates: texture kernels and the halfwidth convention
(`texture/kernels.py`, 484 lines); texture components and fibres; equal-area vs latitude-longitude
S2 grid construction (partly covered now by the raster-bias treatment).

## Retiring `docs/tex/`: LaTeX Notes Become Rendered MyST Pages — COMPLETE (2026-08-11)

**Objective.** Make MyST Markdown the single canonical source for the scientific notes that
previously lived as LaTeX under `docs/tex/`. The notes must render as first-class pages inside the
Sphinx site instead of being emitted as raw `.tex` download links, and PDF output must come from
Sphinx's own `latexpdf` builder rather than a separate `latexmk` toolchain.

**Why.** The Sphinx site is already the primary browsable and searchable documentation surface, but
Sphinx has no LaTeX-parsing extension configured, so every `docs/tex/*.tex` link was copied verbatim
into `_downloads/<hash>/` and surfaced as a download link. The derivations — the actual scientific
content of these notes — were therefore invisible on the site. Keeping LaTeX canonical while also
rendering it in Sphinx would require a conversion pipeline and would leave two representations of
the same content free to drift. Collapsing to one source removes the drift risk and the extra
toolchain at the same time.

**Why the conversion is safe.** A survey of all 37 notes found a deliberately narrow LaTeX
vocabulary and, critically, **no** `\cite`, `\ref`, `\includegraphics`, `\newcommand`, `\input`, or
TikZ anywhere. Environments in use are only: `itemize`, `enumerate`, `equation`, `align`, `align*`,
`cases`, `tabular`, `center`, `quote`, `thebibliography`, `document`. All of these have direct MyST
equivalents, and `myst_enable_extensions` already carries `amsmath` and `dollarmath`. The only
cross-reference machinery present is `\label`/`\eqref` on equations, which maps onto MyST's
labelled-equation syntax. Nothing in the corpus requires a LaTeX-only feature.

**Second finding: the LaTeX was already not what the standard claimed.** 26 of the 37 files have no
`\documentclass` at all — they are section fragments that were never standalone-compilable. The
"canonical LaTeX compiled with latexmk" posture in `AGENTS.md` and
`docs/standards/scientific_notes_and_figures.md` described a build that could not have run for most of the
corpus. This migration makes the documented posture and the actual state agree.

**Scope decision — format migration only, not content dedup.** Four notes overlap existing
`docs/site/algorithms/` pages, which are strict supersets of them. Those pages are a *different
documentation layer* by design (`algorithms/index.md` lists "theory notes" as its own row), so all
37 notes convert and the layering is preserved. Deduplicating the derivation layer against the
algorithms layer is a separate editorial question and is explicitly **not** part of this task.

**Target layout.** All 37 notes become pages under `docs/site/theory/`, flat (all basenames are
unique), grouped by `toctree` into Theory / Algorithms / Validation / Foundations. `docs/tex/` is
deleted. PDF comes from `sphinx -b latexpdf`.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 1 | Survey LaTeX vocabulary; confirm clean convertibility | done | 474a071 |
| 2 | Record objective, rationale, and target layout in ledger | done | 474a071 |
| 3 | Convert 37 `.tex` notes to MyST pages under `docs/site/theory/` | done | 474a071 |
| 4 | Rewrite `docs/site/theory/index.md` with grouped toctrees | done | 474a071 |
| 5 | Rewire repo-wide references, `src/pytex/cli.py`, `scripts/check_repo_integrity.py` | done | (this commit) |
| 6 | Make MyST canonical in `AGENTS.md`, `mission.md`, `specifications.md`, both standards | done | (this commit) |
| 7 | Delete `docs/tex/`; confirm build clean and `pytest`/`ruff`/`mypy` green | done | (this commit) |

### Verification

- `pytest`: 5022 passed. `ruff`, `mypy` (114 files), and `check_repo_integrity.py`: clean.
- Sphinx HTML build: 602 warnings, down from 604, and all of them pre-existing (autodoc
  duplicate-object and docstring definition-list warnings that predate this task). The two removed
  were the `eq-wedge` dangling references described below. **No new warning was introduced**, and
  no `.tex` file is emitted into `_downloads/` any more.
- Spot-checked rendering: display and inline mathematics emit MathJax `\[…\]`, `align*` renders
  through the amsmath container, both converted tables are real `<table>` elements, and a labelled
  equation renders as numbered `(1)` with a working `#equation-eq-cbed-s` anchor.

### Decisions worth carrying forward

- **The notes are hand-maintained sources now, not generated output.** The migration script lived in
  the session scratchpad and was deliberately never committed; committing a generator would wrongly
  imply the `.md` files are regenerable and would invite someone to re-run it over edited prose.
- **`docs/standards/latex_and_figures.md` was renamed** to `scientific_notes_and_figures.md`. A file
  by the old name asserting that MyST is canonical is exactly the incoherence this task removed.
  16 references plus two tests were repointed.
- **Two tests encoded the retired policy.**
  `test_foundational_docs_agree_on_hybrid_documentation_policy` asserted the literal string `latex`
  in every foundational document. It is now `..._layered_documentation_policy` (Sphinx / notes /
  SVG), joined by `test_no_foundational_doc_still_calls_latex_canonical`, which fails if the old
  posture reappears.
- **One conversion defect was found by the build, not by inspection.** A `\label` inside an `align`
  environment was dropped, dangling `\eqref{eq:wedge}`. An audit of the deleted sources confirmed
  `eq:wedge` was the only such label in the corpus; it is now a `{math}` directive with `:label:`.
- **Historical records were repointed, not rewritten.** CHANGELOG and older ledger entries had their
  paths updated so they still resolve, but their claims were left alone — except one phrase ("Two
  canonical LaTeX theory notes") whose wording contradicted the path beside it.

### Deliberately not done

Four notes (`orientation_relationship_determination`, `orientation_relationship_index_correspondence`,
`saed_ratio_angle_indexing`, `tem_specimen_tilt_navigation`) substantially overlap their
`docs/site/algorithms/` counterparts, which are strict supersets of them. Both layers were kept,
because `algorithms/index.md` declares them distinct layers, and `theory/index.md` now carries a
table cross-linking each pair. **Whether the derivation layer should be merged into the algorithms
layer is an open editorial question**, excluded here on purpose so that a format migration could not
quietly become a content deletion.

## Double Diffraction In The Kinematic SAED Engine — CAPABILITY COMPLETE (2026-08-11)

**Objective.** Let the kinematic zone-axis engine optionally emit the reflections that appear in
a real TEM pattern only because a diffracted beam re-diffracts: a reflection whose structure
factor vanishes still appears when its indices are the algebraic sum of two reflections that are
themselves excited. Silicon 002 along [110] and hcp 0001-type absences are the canonical cases.
Such spots must be **labelled as kinematically forbidden and present only through double
diffraction**, never silently mixed in with genuine reflections.

**Why the integer-sum rule is the right model here.** The engine is kinematic; it cannot solve
the dynamical coupling that actually produces these spots. What it *can* do exactly is the
selection rule: the doubly-diffracted beam leaves the crystal along `g1 + g2`, so the set of
reachable spots is the set of pairwise integer sums of the excited reflections. That set is a
geometric fact independent of the dynamical theory, which is why the sum rule is what textbooks
state. Only the *intensity* is a model choice, and it is declared as indicative.

**Key structural fact discovered while designing this.** Lattice-centring conditions define a
sublattice of reciprocal space, and a sublattice is closed under addition: the sum of two
centring-allowed reflections is always centring-allowed. So double diffraction can *never*
revive a centring absence (F-centred 100 stays dark), only a **basis** absence from a glide
plane, a screw axis, or a motif — which is exactly what the literature reports. The
implementation therefore needs no special centring handling; the closure does the work, and a
test pins it.

**Where the code goes.** `src/pytex/diffraction/kinematic.py` — the vectorized engine and
`SpotTable`, which is what `pytex.diffraction.composite` and the TEM indexing surfaces consume.
The legacy loop-based `saed.generate_saed_pattern` is not extended; its docstring is corrected
to point at the engine that models this.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 0 | Survey the SAED surfaces; open this ledger entry | done | (step 0) |
| 1 | Engine: config flags, sum-rule enumeration, `SpotTable` designation arrays, `describe()` | done | (step 1) |
| 2 | Export/reflection-table column and JSON contract | done | (step 2) |
| 3 | Plotting designation (distinct marker for forbidden spots) | done | (step 3) |
| 4 | Docs, worked example, validation matrix, CHANGELOG | done | (step 3) |

### Current worktree state

Step 1 landed. `pytex.diffraction.kinematic` carries
`KinematicSimulationConfig.include_double_diffraction` / `.double_diffraction_coupling`, the
public `double_diffraction_sums` selection rule, `SpotTable.is_double_diffraction` /
`.double_diffraction_parents` / `.forbidden_mask()` /
`.double_diffraction_origin_label(row)`, and a `describe()` that states the designation.
Verified on silicon [110]: 002, 222 and 442 appear, each flagged, each with a parent pair of
{111}-type reflections, and all weaker than every genuine reflection. The class-model
diffraction SVG was regenerated for the two new `SpotTable` fields.

Step 2 landed. `ReflectionTableRow.double_diffraction_origin` (with the derived
`is_double_diffraction`), two appended `REFLECTION_TABLE_COLUMNS`, a `(dd)` marker plus footnote
in the Markdown rendering, a `describe()` sentence, the two config knobs in the manifest, and
the matching required properties in `schemas/composite_saed_manifest.schema.json`.

One thing to know when writing further tests: the bcc/hcp Burgers composite along the **[110]**
parent zone yields no double-diffraction spots at all, because the mapped child zone puts every
basis-absent hcp reflection parallel to the beam. The [001] parent zone yields six.

Steps 3 and 4 landed; the capability is complete. Rendering splits forbidden spots into their
own hollow collection (same marker shape, own gid and legend entry). Documentation: a theory
subsection in `docs/site/theory/reciprocal_space_and_kinematic_spots.md`, a section in
`docs/site/workflows/saed_generation.md`, a validation-matrix row, registry entries for
`I_dd` and `c`, corrected limitation claims in tutorials 12 and 21, corrected legacy
`generate_saed_pattern` docstrings, a CHANGELOG entry, and the worked example
`kinematic-silicon-double-diffraction-002`.

**On the worked example's expected value.** The check is `r(002)/r(004) = 0.5` on silicon
[110], taken only if the engine also flags the row as forbidden. It is the sharpest single
number available here: the ratio follows from `|g_00l| = l/a` for a cubic cell, so it is
independent of the lattice parameter and the camera constant, and it verifies the selection
rule, the flagging, and the detector projection at once. No spot count was used as an expected
value, because a count depends on the enumeration cut-offs rather than on physics.

Verified: `ruff` and `mypy src` clean, full `tests/unit` suite green.

### Next task

None claimed. Two follow-ons a later session could pick up, in rough priority order:

1. The legacy `saed.generate_saed_pattern` still has no double-diffraction path and now
   documents that it delegates the capability. If that surface is meant to stay, it should
   either grow the option or be formally deprecated in favour of the engine.
2. `double_diffraction_coupling` is a single global constant. A thickness-aware estimate would
   be the natural next increment, but it needs a defensible source before it is worth having.

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
`docs/site/theory/stereographic_kikuchi_maps.md`, the worked example
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
`docs/site/theory/dynamical_cbed_and_symmetry_determination.md`.

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
`docs/site/theory/orientation_representations.md`.

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
`docs/site/theory/convergent_beam_electron_diffraction.md`.

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

## Kearns f Parameter: Four Estimation Routes, Theory Note, Tutorial — COMPLETE (2026-08-13)

**Objective.** Make the Kearns orientation parameter `f` a first-class, explainable surface of
`pytex.texture`, estimated by all three routes the Zr literature uses — the original
diffractogram/inverse-pole-figure route, the basal pole figure route, and the ODF (or discrete
orientation) route — with the corrections each one needs, a canonical theory note, executable
worked examples with cited expected values, and a tutorial notebook that both demonstrates the
methods on real XRDML data and cross-checks them against simulated textures whose `f` is known
in closed form.

**Normative sources read (PDFs local-only, `kearns_parameter_data_references/`, untracked).**

- J. J. Kearns, *Thermal expansion and preferred orientation in Zircaloy*, WAPD-TM-472,
  Bettis Atomic Power Laboratory (November 1965). Defines `f` (its Eqs. 1-7) and the
  diffractogram route; Table 3 is the pinned numerical baseline.
- J. L. Baron *et al.*, *Interlaboratories tests of textures of Zircaloy-4 tubes. Part 1*,
  Textures and Microstructures **12** (1990) 125-140, doi:10.1155/TSM.12.125. Its Eqs. 4-5 are
  the pole-figure route, with the incomplete-pole-figure pseudo-norm of Kern and Bergmann.
- R. A. Holt and S. A. Aldridge, *J. Nucl. Mater.* **135** (1985) 246-259 — `F_d = sum V(theta) cos^2(theta)`,
  the resolved-basal-pole form used throughout the CANDU literature.
- K. V. Mani Krishna *et al.*, *J. Nucl. Mater.* **414** (2011) 492-497,
  doi:10.1016/j.jnucmat.2011.04.065 — comparison of the IPF, PF, ODF and EBSD routes, their
  cross-section dependence, and the normalization that the IPF route needs.

**The unifying formulation adopted.** All four routes estimate the same object: the second-moment
(orientation) tensor of the basal-pole direction distribution,
`A = <c c^T>`, whence `f(d) = d^T A d` for any specimen direction `d`. The literature's
sum rule `f_RD + f_TD + f_ND = 1` is then `tr(A) = 1`, which holds identically rather than
approximately, and the isotropic value `1/3` is `A = I/3`. Each experimental route is a different
estimator of `A` (or, for the diffractogram route, of one diagonal element of it).

**Finding to carry into the docs.** The reference 2-theta XRDML scans in
`reference_exp_data/*/2theta*.xrdml` are *fixed-omega detector scans* (`scanAxis="2Theta"`,
`Omega` a `commonPosition` of 15 deg), not coupled Bragg-Brentano scans (which appear in the same
corpus as `scanAxis="Gonio"` with both axes ranged, e.g. `stress01-13.xrdml`). The diffraction
vector therefore sits at `theta - omega` from the specimen normal, reaching 45 deg at the top of
the range, so the Kearns assignment of each reflection's intensity to the section normal does not
hold as written. The module carries a per-reflection specimen direction rather than assuming the
surface normal, and reports the spread as a diagnostic.

**Second finding.** Kearns 1965 Table 3 contains an arithmetic slip in the transverse-section
block: the `70-80` row's `V cos^2 phi` cell reads `0.0214` where `0.353 * cos^2(75 deg) = 0.0237`,
so the quoted `f = 0.0508` should be `0.0529` from his own volume fractions, or `0.0526`
recomputed from his intensities. The longitudinal block reproduces exactly
(`f = 0.4879` against the quoted `0.488`), so the longitudinal block is what gets pinned as the
worked example, with the discrepancy documented rather than silently averaged away.

### Step ledger

| # | Step | Status | Commit |
| --- | --- | --- | --- |
| 1 | Reference analysis, formulation, ledger | done | (this commit) |
| 2 | `src/pytex/texture/kearns.py` + tests (70 tests) | done | (this commit) |
| 3 | Theory note, symbol registry, docs index, parity matrix | done | (this commit) |
| 4 | Worked examples + gallery regeneration | done | (this commit) |
| 5 | Tutorial notebook 31 (55 cells) | done | (this commit) |

**Worktree state.** `kearns_parameter_data_references/` is untracked and must stay that way: it
holds copyrighted PDFs and several MB of instrument data. A `.gitignore` entry lands with step 2.

**Third finding, fixed in this commit.** `raster_solid_angle_weights` extends the outermost ring's
band outwards by its own half step, which for a hemispherical raster ending exactly at 90 degrees
pushes that band past the equator and gives the equatorial ring close to twice the solid angle it
owns. Because `cos^2` is zero there, the excess pulls every Kearns integral down: the spherical
mean of `cos^2` over a 5 degree raster comes out at 0.3196 against the exact 1/3, a -4.1 percent
error that the existing `pole-figure-raster-weighted-mean-converges` worked example documents as
"converging from below". A new optional `polar_max_deg` keyword bounds the outermost band; passing
90.0 reduces the error to -0.06 percent. The default is unchanged so the pinned worked example is
untouched, and `kearns_from_pole_figure` passes 90.0 for antipodal figures. **Whether the default
should change is a separate question left open**; it would require regenerating that example.

**Also noted, not fixed.** `KernelSpec.evaluate` silently degenerates for de la Vallee Poussin
halfwidths below about 0.5 degrees: `np.isclose(cos(halfwidth/2), 1.0)` becomes true at the default
tolerance and the exponent collapses to 1, giving a near-uniform kernel instead of a very sharp one.
Out of scope here; no realistic halfwidth reaches it.

**Fourth finding, from writing the tutorial.** The `normalization="harris"` option originally
allowed the diffractogram route to run with no reference intensities at all, treating raw peak
areas as pole densities. That is wrong — in a *random* powder the alpha-Zr reflections span a
factor of twenty from structure factor and multiplicity alone — and it silently produced a
plausible number. Both normalizations now require reference intensities (calculated ones are fine,
since only ratios enter). The tutorial also makes explicit what the Harris rescaling does *not*
buy: it changes absolute pole densities by about the 23 percent Kearns measured, and changes `f`
by exactly nothing, because `f` is a ratio and a common scale cancels. A regression test pins the
two normalizations to agree to 1e-14.

**Reference-data findings recorded for the notebook.** On the three principal sections of the
rolled Zircaloy-2 plate, the pole-figure route gives f = 0.658, 0.523, 0.251 along the three
section normals, summing to 1.432 instead of 1 — a 43 percent one-sided excess from truncation at
85 degrees plus uncorrected defocusing, both of which suppress the high-tilt rings where cos^2 is
smallest. The ODF route on four incomplete pole figures per section reproduces the same values to
within 0.03 and the same 1.44 sum, confirming that inversion inherits the bias rather than
repairing it. The theta-2theta route on the same plate's scans gives 0.472, 0.301, 0.140, summing
to 0.913; after normalising both triads the two independent routes agree to a few hundredths,
which is about the interlaboratory spread Baron et al. (1990) report. The dataset contains no
random-powder standard, so the defocusing curve cannot be measured from it — which is itself the
lesson the notebook draws.

### Verification

- `ruff check`, `ruff format --check`, `mypy` (132 source files) and the full `pytest tests/unit`
  lane are green. 76 new tests in `tests/unit/test_kearns_parameter.py`.
- `python -m sphinx -b html docs/site` exits 0. `theory/kearns_parameter_and_basal_pole_texture`
  and `tutorials/notebooks/31_kearns_parameter` both render, the notebook executes under
  `nb_execution_mode = "cache"` with live outputs, and no build warning mentions either page. The
  603 warnings in the log are the pre-existing autodoc duplicate-object noise on
  `api/full_reference`.
- The notebook was executed both with the reference data present and with
  `PYTEX_KEARNS_DATA` pointing at a missing directory; it completes either way, which is what the
  docs build needs since the data is untracked.
- Worked-example gallery regenerated; `tests/unit/test_worked_examples.py` recomputes all four
  Kearns examples against their cited values on every run.

### What was deliberately not done

- **The `raster_solid_angle_weights` default was not changed.** Bounding the outermost band at the
  measured edge is the better quadrature everywhere, not only for Kearns integrals, but flipping
  the default would move the pinned `pole-figure-raster-weighted-mean-converges` worked example
  (0.31960 -> 0.33354 at a 5 degree step) and any downstream m.r.d. normalization. That is a
  separate, deliberate change with its own commit; the keyword exists so it can be made
  incrementally.
- **`KernelSpec.evaluate` degenerates below about a 0.5 degree de la Vallee Poussin halfwidth**
  (`np.isclose(cos(halfwidth/2), 1.0)` becomes true and the exponent collapses to 1, giving a
  near-uniform kernel where a very sharp one was asked for). Pre-existing, unreachable at any
  realistic halfwidth, and out of scope here.
- **No defocusing model.** The reference corpus contains no random-powder standard, so the
  defocusing curve cannot be measured from it; `PoleFigureCorrectionSpec.defocus_factors` already
  accepts one when a standard exists. A measured-standard workflow remains the gap.
- **No EBSD-specific Kearns surface.** `kearns_from_orientations` already takes grain weights, so
  a `CrystalMap` convenience wrapper would add API without adding capability.

### Next task

None claimed. This goal is complete.

## Follow-up: Table 3 reproduced column by column (2026-08-13)

Tutorial 31's Kearns-Table-3 section was expanded from 3 cells to 18: both blocks of his table
are now reproduced one column at a time, with his printed values beside the recomputed ones at
every step, so a disagreement localises to a cell rather than to a final number.

Two numbers in the earlier commits were wrong and are corrected here. Recomputing the transverse
block from Kearns' own volume-fraction column gives **0.0529**, not the 0.0532 stated in the
theory note, CHANGELOG, ledger and worked example; from his intensity column it gives 0.0526. And
the notebook's closing sum-rule cell paired the wrong value with the rod axis. The correct
pairing follows his Figure 6 captions: the *longitudinal section* measures a direction
perpendicular to the rod axis, the *transverse section* measures the axis itself. A swaged rod is
axially symmetric, so the triad is (f_perp, f_perp, f_axial) = (0.488, 0.488, 0.053), closing to
**1.028** — a 2.8 percent closure error, which validates his measurement rather than impugning
it. The earlier cell computed f_long + 2 f_trans = 0.593 and drew the opposite conclusion.

Worth recording as a methodological point the section now makes: correcting the bad cell moves
the triad sum by 0.002, well inside its own scatter, so the sum rule could not have found this
error. The row arithmetic did. The two checks are complementary, and both need the intermediate
columns to be published.
