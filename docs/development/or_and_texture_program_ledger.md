# Program Ledger: OR Visualization, Interfaces, And Texture Quantification

The resumable state of the long-horizon program defined in
[Vision And Plan: Interface Crystallography, Composite Visualization, And Texture Quantification](../roadmap/vision_interface_crystallography_and_texture_quantification.md).

**How to resume.** Read §1 for where the program stands, then §2 for the milestone currently in
progress and its next concrete step. Everything below §3 is append-only history. This file and the
git log are the only things a later session needs; nothing about this program lives only in an
agent's head.

**Scope note.** MTEX parity is explicitly **deferred by the user** (2026-08-29) and is not a gate on
any milestone here. It returns as a gate only when the paper is being prepared for submission. No
document may claim MTEX parity in the meantime — deferring the campaign does not license the claim.

---

## 1. Program Status

| Milestone | Content | State |
| --- | --- | --- |
| **M1** | Kearns parameter in the GUI (T1) | **Complete** (2026-08-29) |
| **M2** | F15 variant-aware composite scenes, F18 OR stereogram | **Complete** (2026-08-29) |
| **M3** | F19 composite crystal viewer in the workbench, F17 OR dossier | **Complete** (2026-08-29) |
| **M4** | F21–F23 measured-pair OR workbench | Not started |
| **M5** | T3 axial specimen symmetry, T2 ghost correction | Not started |
| **M6** | F16 interface crystallography, Program D contracts + CLI, T5 uncertainty | Not started |
| **M7** | F20 PTMC / habit-plane prediction | Not started — user has committed to this as the long-horizon goal |

The user's stated order: **Kearns GUI first, PTMC/habit-plane last as the long-horizon goal**, with
everything else in between. The vision document's M4 go/no-go on PTMC is therefore resolved: it is a
**go**, scheduled last as M7.

---

## 2. Current Milestone — M4: The Measured-Pair OR Workbench

**Goal.** F21–F23: take two *measured* orientations — typed in, or picked off an EBSD map — and
answer "what relationship is this, stated in integers, and what does it look like". Everything M2
and M3 built is the answer half; M4 is the input half.

**Why now.** M3 is complete, so a relationship declared by name already produces a composite scene,
a stereogram, a contact sheet and an exportable dossier. What no route yet does is start from a
measurement. The existing `characterize_orientation_relationship` and
`fit_orientation_relationship` are the numerical core and need no new crystallography — this is a
plumbing and presentation milestone, mostly in `app/`.

### Sub-milestones

| Step | Content | State |
| --- | --- | --- |
| **M4a** | F21: a panel taking two orientations (Euler, matrix, or axis/angle) and reporting the relationship — catalog deviations, the fitted rotation, the integer parallelism statement | Not started |
| **M4b** | F22: the OR stereogram beside the 3-D view, sharing the camera | Not started |
| **M4c** | F23: pick two grains on an EBSD map and flow their mean orientations into M4a with no retyping | Not started |

### Next concrete step

**M4a — the measured-pair panel.** Add a service operation taking two parent/child orientations and
returning: the deviation from every catalog relationship (`characterize_orientation_relationship`),
the fitted relationship and its residuals (`fit_orientation_relationship`), and the integer
parallelism statement of the fit (`describe_orientation_relationship`, plus the F4 search). Then a
panel over it.

The output of that operation is the natural input to everything M2/M3 landed: hand its fitted
relationship to `or_dossier`, `plot_or_stereogram` or `variants.composite_scene` and the whole
answer half already works. Build the operation so that handing over is one call, not a re-derivation.

Two traps to carry in, unchanged since M2:

1. Name parallelisms from the **variant's own** `parallel_planes` / `parallel_directions`, and spell
   planes through `pytex.core.miller.canonicalize_sign`, so one plane has one spelling everywhere.
2. Say what each deviation measures. A fit residual, a catalog separation and a rationalization
   residual are three different quantities, and a panel that labelled them all "deviation" would be
   stating something false about two of them.

**Known and deliberately untouched.**---

## 3. History

### 2026-08-29 — M3c complete, and M3 with it: the OR dossier (F17)

**What shipped.** `src/pytex/core/or_dossier.py` — `or_dossier(...)`, `ORDossier` and its four block
types, `or_dossier_schema_path()` — with `schemas/or_dossier.schema.json`, 29 tests in
`tests/unit/test_or_dossier.py`, a worked example, the workflow page
`docs/site/workflows/or_dossier.md` with its toctree and API-index entries, and
`Lattice.volume_angstrom3()`.

**Where it lives, and why not where the plan said.** The vision document names
`pytex.core.transformation.or_dossier`. That module is already ~4100 lines, and the dossier must
import `find_parallel_planes`, `variant_close_packed_groups` and (lazily) the plotting stack, so
putting it there would have meant either more growth or an import cycle. It is a sibling module
instead, exported from `pytex` and `pytex.core` under the documented names. The deviation is stated
here rather than left for a later reader to discover.

**The rule, and how it is checked.** The dossier calls the existing functions and reimplements none
of them. The tests do not compare against recorded output: each value is asserted equal to the
function a reader would check it against — the volume against `Lattice.volume_angstrom3`, the
correspondence matrices against `OrientationRelationship`, the spectrum against
`intervariant_misorientation_angles_deg` — and those functions are in turn pinned to the published
figures (24 variants, 4 packets, Morito's ten angles, 42.85 deg).

**Two presentation decisions worth carrying forward.**

1. *The `origin` column.* A `defining` parallelism's deviation is zero by construction and measures
   nothing; a `discovered` one's is a rationalization residual. Listing them in one table without
   saying which is which would repeat M1's mistake of printing a number that is true and not the
   number the reader assumes. The column, the schema description and the prose all say it.

2. *The interface block is `null`, not absent.* F16 is not implemented, and a missing section reads
   as an oversight where an explicit null with a sentence reads as a scope statement.

**Verification.** `ruff check .` clean; `mypy src` clean over 153 files; `pytest tests/unit` green
(full suite); the worked-example gallery and the class-model atlas regenerated (six new public
dataclasses moved the atlas counts from 284/267 to 290/273).

**Not done.** The figure bundle omits the variant pole figure and the per-variant SAED patterns:
they need a measured parent orientation and a diffraction geometry respectively, neither implied by
a relationship, and the docstring says so rather than shipping an empty figure. No application
surface over the dossier — it is a library object for now.

### 2026-08-29 — M3b complete: the composite crystal viewer (F19)

**What shipped.** Two new views in the Variants panel — *Both crystals of one variant* and *Every
variant at once* — built on a new `core/compositescene.js`, a `setTitle` on the shared plot frame,
contact-sheet styles in `app.css`, and three `export` keywords on `panels/crystal.js`. One Playwright
test (the browser suite is 53 green), and rewritten prose in
`docs/site/workflows/workbench_application.md`, which previously said these operations had no user
interface.

**No second renderer.** The panels compose scene payloads into the shape the crystal viewer's
renderer already takes, which was the whole design of M3a's payloads. That buys more than reuse: one
renderer means one global depth sort, so parent and child atoms occlude each other correctly instead
of one crystal being drawn wholly in front of the other.

**Three findings.**

1. *Colour had to stop meaning the element.* Both phases of an orientation relationship are usually
   the same element — austenite and ferrite are both iron — so the first working composite was a
   single orange blob in which no parallelism could be read. Colour now carries the **phase**, and
   the legend says so in those words, because a viewer that quietly changed what a colour meant
   would be worse than one that could not tell the crystals apart.

2. *A stale response could overwrite a new one.* Switching view started a request without cancelling
   the previous one, and the pole figure already running could land after the composite scene,
   overwrite `state.result`, and hand the new view's drawing code the old view's data. It threw here
   because the shapes disagree; where two views' shapes happen to agree it would have drawn stale
   numbers silently, which is worse. Runs now carry a token. This was a pre-existing defect in the
   panel that the new views exposed.

3. *One plane, two spellings.* The variant table named variant 10's parent plane `(1 1 -1)` while
   the 3-D overlay beside it named the same plane `(-1 -1 1)`. A plane has no sign, and the sign a
   symmetry image comes back with is an artefact. Every producer of a plane label in this path now
   goes through `pytex.core.miller.canonicalize_sign` — the rule that already existed — rather than
   through a private copy: the first fix added a fourth implementation to `core.notation` and was
   replaced by the existing one. Directions deliberately keep their sign, and a test asserts it.

**Verification.** `ruff check .` clean; `mypy src` clean over 152 files; `pytest tests/unit` green;
53 Playwright tests green against a loopback server, the new one run three times in a row after the
race was fixed; both views driven by hand in the browser with an empty console.

**Not done.** No publication renderer for the composite views: the *Figure* button and its format
picker are disabled outside the pole figure rather than silently producing the wrong figure. A
`variants.composite_render` belongs with the dossier's figure bundle in M3c.

### 2026-08-29 — M3a complete: the composite-scene service layer

**What shipped.** `variants.composite_scene` and `variants.contact_sheet` in
`src/pytex/app/services/variants.py`, three example scenarios, 24 tests in
`tests/unit/test_app_variants.py`, a section in `docs/site/workflows/workbench_application.md`
saying plainly that these two operations have no panel yet, and a changelog entry.

**The payload decision worth recording.** `composite_scene` returns both crystals **already placed
in one world frame** (the parent crystal frame), because that is what makes a shared camera free
and undriftable. `contact_sheet` deliberately does the opposite: the two structures go once, in
their own crystal frames, with a 3x3 placement matrix per variant. Twenty-four placed copies of
both crystals would be tens of megabytes to say what a matrix multiply says exactly, and applying
one of those matrices is the same arithmetic the camera already does — so this still honours the
vision document's Decision 4 (all crystallography in Python; the browser multiplies matrices only).
The two shapes are different, so a test asserts they agree: the matrix the sheet gives for variant
k is the matrix `composite_scene` places variant k by. A grid and a detail view that disagreed
would be worse than either alone.

**Two smaller findings.** The side-by-side offset is measured from the two *unplaced* scenes, so it
does not change with the variant — measured per variant, the crystals would jump as a user stepped
through the family, and a test pins that. And plane labels are sign-canonicalized: a plane has no
sign, and leaving `(111)` and `(-1 -1 -1)` as separate entries made the 24 Kurdjumov-Sachs variants
name eight parent planes in a table whose packet column says four.

**Verification.** `ruff check .` clean; `mypy src` clean over 152 files; `pytest tests/unit` green
(full suite); `test_app_variants.py`, `test_app_manifest.py` and `test_app_export.py` green, which
covers the runs-from-defaults gate and the example executor.

**Not done.** No panel and no JavaScript — that is M3b, and the documentation says so rather than
describing a UI that does not exist. No worked example: the workbench service group covers this
surface class already, and both operations return scenes rather than a numerical result with an
independent reference value.

### 2026-08-29 — M2b complete, and M2 with it: the OR stereogram (F18)

**What shipped.** In `src/pytex/plotting/spherical.py`: `ORStereogramPair`,
`or_stereogram_pairs(...)`, `build_or_stereogram_figure_spec(...)` and `plot_or_stereogram(...)`,
exported from `pytex.plotting` and (the plot function) from `pytex`. 17 tests in the new
`tests/unit/test_or_stereogram.py`, one worked example, a section in
`docs/site/workflows/stereographic_projections.md`, and a changelog entry.

The figure carries what F18 asked for: parent poles open and child poles filled, both in the parent
crystal frame; a tie-line per pair labelled with its deviation; and, for plane pairs, the great
circles of both planes, parent dashed and child solid, so a plane parallelism reads as two
coincident circles rather than two coincident points.

**Two findings worth keeping.**

1. *The deviation label means something narrower than it looks.* `find_parallel_planes` maps a
   parent plane to its **exact** child image — which is parallel by construction — and then
   rationalizes it, so the angle it reports is the **rationalization residual**, not a departure
   from parallelism. The figure draws the parent pole against the *rationalized* child pole, so the
   gap on the net and the printed number are the same quantity and agree; but the docstrings and
   the workflow page now say what that quantity is, because "0.33 deg" invites the wrong reading.
   This is the same class of defect as M1's triad closure: a number that is true but not the number
   the reader assumes.

2. *An equatorial pole can split a zero-deviation tie-line across the whole net.* Kurdjumov-Sachs
   variants 7 and 9 have a defining direction lying in the equatorial plane, and the variant
   rotation returns the child copy at `z = -8e-16`. The antipodal fold breaks ties on the equator,
   so the two ends folded to **opposite rims** and the figure drew a diameter where the
   crystallography has a single point. The pair is now folded once, from the parent pole's
   decision, with the child's sub-noise dip flattened onto it; a pair that genuinely straddles the
   equator is left alone and its tie-line is split rather than drawn as a chord. A test walks all
   24 variants and pins the worst pole gap and the worst great-circle gap below 1e-9, and the
   worked example states the same identity.

Separately, pair angles are computed as `2 atan2(|a-b|, |a+b|)` rather than `arccos(a·b)`: arccos
loses half its significant digits exactly where these angles live, and reported 8.5e-7 deg for a
parallelism exact to machine precision.

**Verification.** `ruff check .` clean; `mypy src` clean over 152 files; `pytest tests/unit` green
(full suite); worked-example, documentation-policy, repo-integrity and public-API-docstring suites
green; the gallery regenerated with `python scripts/generate_worked_examples.py`.

**Repository mechanics worth knowing before the next increment.** Adding one public dataclass
(`ORStereogramPair`) broke three tests in `tests/unit/test_class_model_atlas.py`: the committed
`docs/figures/class_hierarchy.svg` and `class_model_architecture.svg` are byte-compared against
`python scripts/generate_class_model_figures.py`, and `docs/site/architecture/class_model_atlas.md`
states the class and dataclass counts in prose, which a test parses. Any new public class therefore
costs a figure regeneration and a two-number edit on that page (283/266 became 284/267). Nothing is
wrong with that — it is the no-hand-transcribed-numbers rule doing its job — but it is not obvious
from the failure message.

**Not done, deliberately.** No theory note: the figure states an existing relationship rather than
introducing new theory, and every number it prints is owned by a function that already has one. The
shared Wulff-net grid's own rim jumps were left alone — see the note in §2. The Sphinx warning count
was not re-checked against a clean-worktree build.

### 2026-08-29 — M2a complete: variant-aware composite scenes (F15)

**What shipped.** `WorldScene3D.from_orientation_relationship(..., variant=...)` accepting a
`TransformationVariant` or a one-based index; `WorldScene3D.variant_scenes(...)`;
`plotting.scene3d.render_variant_contact_sheet(...)`; and, in `core/transformation.py`, the
properties `TransformationVariant.parallel_planes`, `.parallel_directions`,
`.parent_symmetry_operator`, `.child_symmetry_operator`. Exported from `pytex.plotting` and, for
the renderer, from `pytex`. 13 new tests across `tests/unit/test_scene3d_composition.py` and
`tests/unit/test_transformation.py`, one worked example, and prose in
`docs/site/concepts/visualization_primitives.md`.

**The trap, and where it is now closed.** A variant is `V = S_c R S_p^T`, so the objects actually
parallel under it are the defining pair carried by those operators —
`(S_p n_parent) || (S_c n_child)` — not the nominal pair the relationship was written with. The
fix deliberately lives on `TransformationVariant`, not in the scene builder: any future consumer
(the M3 viewer, the F17 dossier, the F18 stereogram) that asks a variant for its parallelism gets
the right pair without knowing the derivation. Concretely, over the 24 Kurdjumov-Sachs variants
the property names **four** distinct parent {111} members — the Morito packet planes — where the
nominal pair would name one; a test asserts that, and another asserts that the nominal pair opens
a visible angle on variant 17 while the variant's own pair closes to 1e-12.

The arrows and patches are now labelled with the indices they draw, e.g. `(1 -1 1) ∥ (0 1 1)`,
through `pytex.core.notation.format_miller_indices` rather than a generic "∥ plane". A figure that
states which plane it is drawing cannot quietly be the wrong one.

**Verification.** `ruff check .` clean; `mypy src` clean over 152 files; `pytest tests/unit` green
(exit 0, full suite, before the docs and example additions; the worked-example, documentation-policy
and repo-integrity suites re-run green after). The worked-example gallery was regenerated with
`python scripts/generate_worked_examples.py`.

**Not done, deliberately.** No theory-note change: the derivation `V (S_p n) = S_c n'` is one line
and sits in the property's docstring and the worked example's reference, and
`docs/site/theory/phase_transformation_relationship_construction.md` already carries the variant
algebra. The Sphinx warning count was not re-checked against a clean-worktree build; the prose
added is one section in an existing concepts page.

### 2026-08-29 — M1b complete, and a defect in M1a corrected

**The upload routes.** `kearns.from_pole_figure` (one XRDML figure, Baron Eq. (5)) and
`kearns.from_odf` (several figures inverted to an ODF first, with the density/support tensor choice
exposed rather than defaulted). Both reuse `uploaded_file` and `read_xrdml_pole_figure`; both
declare a required `files` parameter with no default and are therefore exempt from the
runs-from-defaults gate, as `texture.measured_pole_figures` already is. The panel gained a file
picker shown only for those routes.

**The defect M1a shipped.** The panel presented `f_RD + f_TD + f_ND = 1` as a green *passed* check.
It is not a check at all for any route that builds a single pole orientation tensor: the tensor
averages `c cᵀ` over unit vectors, so its trace is 1 whatever the data were. That covers three of
the five routes, including both new ones.

`fixtures/xrdml/synthetic_random_standard.xrdml` proved it concretely and is now the regression. It
is a random standard truncated at 60 deg of tilt, so its true `f` is 1/3 in every direction. The
pole-figure route reports **0.518** — wrong by more than half — while the triad sums to **1.0000**.
The panel now writes the verdict three ways: "closes by construction" (styled as neither pass nor
failure) for single-tensor routes, a real pass/fail where sections were measured independently, and
nothing at all for a single-section route. Beside it sits the **coverage** note, which is the
diagnostic that does test these routes.

This is worth remembering beyond this milestone: a quantity that is identically true is not
evidence, and presenting one as a passed test is worse than printing nothing.

**Verification.** `ruff check .` clean; `mypy src` clean over 152 files; `pytest tests/unit` green;
52 Playwright tests green, including a new one asserting that the panel shows no green pass on the
truncated fixture. 45 tests in `tests/unit/test_app_kearns.py`.

**Deliberately not done (M1c).** No worked example and no theory-note change. The
executable-examples standard covers stable public *numerical* surfaces; this milestone added an
application surface over `pytex.texture.kearns`, which already carries both a theory note
(`docs/site/theory/kearns_parameter_and_basal_pole_texture.md`) and its own tests. Revisit only if
the panel grows a numerical surface of its own.

### 2026-08-29 — M1a complete: the Kearns tab

Landed the Kearns panel with three routes, green on the base lane and the browser lane.

**What shipped.** `src/pytex/app/services/kearns.py` (three operations, three examples);
`src/pytex/app/static/js/panels/kearns.js`; the `kearns` documentation link in
`registry.py`; the Texture workspace regrouped in `main.js`; the triad and closure styles in
`app.css`; `tests/unit/test_app_kearns.py` (32 tests); one Playwright test; the workbench guide and
changelog.

**The design that made it worth doing.** All three routes' defaults describe *one* synthetic
specimen — a basal fibre of 30 deg spread about ND, true `f_ND` = 0.6423, generated by the
`synthetic_reflections` recipe of `tests/unit/test_kearns_parameter.py`. Measured on those defaults:
exact 0.6412, diffractogram 0.6421, tilt profile 0.6416. Pressing the button demonstrates a route
recovering a known answer; running all three demonstrates the agreement the panel claims. Both the
service suite and the browser test assert it, so the lesson cannot rot into a false statement.

**Two things found and fixed while building.**

1. The model fibre's help text claimed a wide spread "approaches random, where f goes to 1/3". It
   does not: a Gaussian truncated to the quadrant still leans towards the axis, and at the 90 deg
   maximum `f` is 0.374. The text now says so and a test pins the ceiling, because it would
   otherwise read as a defect to the next person.
2. `_report_payload` initially overwrote the report's own documented `phase` key (a plain name) with
   the application's `PhaseSpec` JSON. That silently changed a published JSON contract; the
   specification is now added beside it as `phase_spec` and a test guards the shape.

**Verification.** `ruff check .` clean; `mypy src` clean over 152 files; `pytest tests/unit` green
(exit 0, full suite); 51 Playwright tests green against a loopback server; the panel driven by hand
in the browser with no console errors. The plot viewBox was widened to 380x140 after measuring the
drawing filling only 40 percent of the frame's canvas width; it now fills 62 percent, the rest being
the frame's ordinary letterboxing.

**Not done, deliberately.** No worked example and no theory-note change: the executable-examples
standard covers stable public *numerical* surfaces, and this increment adds an application surface
over `pytex.texture.kearns`, which already has both. The Sphinx warning count was not re-checked
against a clean-worktree build; the only prose added is in an existing workflow page.

### 2026-08-29 — program opened

The vision and plan document landed on `main` as `41d8518`, with a Sphinx include stub and a toctree
entry. The user then set the order of work: Kearns GUI first, PTMC/habit-plane prediction as the
long-horizon goal, MTEX parity deferred. This ledger opened at that point.
