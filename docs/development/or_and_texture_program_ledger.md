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
| **M2a** | `WorldScene3D.from_orientation_relationship(..., variant=...)`, `variant_scenes(...)`, and the contact-sheet renderer | Not started |
| **M2b** | The OR stereogram (F18): tie-lines between OR-parallel pairs, great circles of the parallel planes, deviations annotated | Not started |

### M1 in one paragraph, for context

The Kearns parameter is now reachable by five routes in a Texture-workspace sub-tab. The design
worth carrying forward: all self-contained defaults describe *one* synthetic specimen with a known
answer, so the opening press demonstrates a route recovering a truth rather than producing a number;
and the panel refuses to present the triad sum as a passed check where it closes by construction.
Both decisions are recorded in the history below and enforced by tests.

### Next concrete step

**M2 — F15 variant-aware composite scenes, and F18 the OR stereogram.** Pure Python, no GUI, so it
carries no browser risk and unblocks M3.

F15, in `src/pytex/plotting/scene3d.py`:

- `WorldScene3D.from_orientation_relationship` currently places the child by
  `relationship.parent_to_child_rotation.inverse()` and has no `variant` parameter — variant 1 only.
  Add `variant: int | TransformationVariant | None = None`, placing the child by
  `variant.parent_to_child_rotation.inverse()`, consistent with the regression-pinned composition
  `g_child = g_parent ∘ Vᵀ`.
- **The trap:** `_orientation_relationship_primitives` draws `relationship.parallel_directions` and
  `parallel_planes` — the *nominal* pair. Under variant k the parent-side objects are the symmetry
  images under that variant's operator, so drawing the nominal pair on variant 17 produces a figure
  that looks right and is wrong. Re-derive them per variant.
- Add `WorldScene3D.variant_scenes(...)` and a contact-sheet renderer.
- Validation: for every variant, the world-frame images of the defining parallel plane normals of
  parent and child must coincide to 1e-12; the child orientations across scenes must equal
  `generate_variants()` as a set under child symmetry.

F18, the OR stereogram: `variant_pole_figure` and `plot_variant_pole_figure` already put child
variant poles and optional parent poles on one net. What is missing is the *pairing* — tie-lines
joining OR-parallel pairs labelled with their deviation, and the great circles of the parallel
planes. `find_parallel_planes` / `find_parallel_directions` supply the pairs and the deviations.

Existing surfaces to build on, all verified present: `generate_variants`,
`variant_close_packed_groups`, `find_parallel_planes`, `find_parallel_directions`,
`variant_pole_figure`, `plot_variant_pole_figure`, `plotting/primitives.py` (`crystal_plane_patch`,
`Arrow3D`, `PlanePatch3D`), `tests/unit/test_scene3d_composition.py`.

---

## 3. History

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
