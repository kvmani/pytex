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
| **M2** | F15 variant-aware composite scenes, F18 OR stereogram | Not started |
| **M3** | F19 composite crystal viewer in the workbench, F17 OR dossier | Not started |
| **M4** | F21–F23 measured-pair OR workbench | Not started |
| **M5** | T3 axial specimen symmetry, T2 ghost correction | Not started |
| **M6** | F16 interface crystallography, Program D contracts + CLI, T5 uncertainty | Not started |
| **M7** | F20 PTMC / habit-plane prediction | Not started — user has committed to this as the long-horizon goal |

The user's stated order: **Kearns GUI first, PTMC/habit-plane last as the long-horizon goal**, with
everything else in between. The vision document's M4 go/no-go on PTMC is therefore resolved: it is a
**go**, scheduled last as M7.

---

## 2. Current Milestone — M2: Variant-Aware Composite Scenes And The OR Stereogram

**Goal.** Render the parent and product crystals of *any* variant of an orientation relationship,
with that variant's own parallel planes and directions drawn on them, and put the OR's parallelism
statement on a stereogram with tie-lines. Pure Python; the workbench viewer that consumes it is M3.

**Why now.** M1 is complete, so the texture half of the user's request is delivered. M2 is the
first half of the OR-visualization request, and it carries no browser risk: it extends
`plotting/scene3d.py`, which already builds the two-crystal composite for variant 1 and is covered
by `tests/unit/test_scene3d_composition.py`.

### Sub-milestones

| Step | Content | State |
| --- | --- | --- |
| **M2a** | `WorldScene3D.from_orientation_relationship(..., variant=...)`, `variant_scenes(...)`, and the contact-sheet renderer | **Complete** (2026-08-29) |
| **M2b** | The OR stereogram (F18): tie-lines between OR-parallel pairs, great circles of the parallel planes, deviations annotated | Not started |

### M1 in one paragraph, for context

The Kearns parameter is now reachable by five routes in a Texture-workspace sub-tab. The design
worth carrying forward: all self-contained defaults describe *one* synthetic specimen with a known
answer, so the opening press demonstrates a route recovering a truth rather than producing a number;
and the panel refuses to present the triad sum as a passed check where it closes by construction.
Both decisions are recorded in the history below and enforced by tests.

### Next concrete step

**M2b — F18, the OR stereogram.** M2a is complete (see the history entry below), so the variant
half of M2 is delivered and M3's viewer has something to consume.

What is missing is the *pairing* on the net. `variant_pole_figure` and `plot_variant_pole_figure`
already put child variant poles and optional parent poles on one stereogram. F18 adds:

- **tie-lines** joining each OR-parallel pair, labelled with the pair's deviation;
- **great circles** of the parallel planes;
- the deviations annotated rather than left implicit.

`find_parallel_planes` / `find_parallel_directions` supply the pairs and the deviations, and
`TransformationVariant.parallel_planes` / `.parallel_directions` (landed in M2a) supply the
per-variant statement a tie-line must be drawn from — do not re-derive the nominal pair here, for
the same reason M2a did not.

Surfaces verified present and now also usable: `generate_variants`, `variant_close_packed_groups`,
`find_parallel_planes`, `find_parallel_directions`, `variant_pole_figure`,
`plot_variant_pole_figure`, `WorldScene3D.variant_scenes`, `render_variant_contact_sheet`,
`TransformationVariant.parallel_planes` / `.parallel_directions`.

---

## 3. History

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
