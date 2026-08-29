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
| **M1** | Kearns parameter in the GUI (T1) | **In progress** |
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

## 2. Current Milestone — M1: The Kearns Tab

**Goal.** Make the Kearns parameter reachable from the workbench by every route
`pytex.texture.kearns` implements, with the triad, its closure check, and the evidence behind each
number on screen.

**Why first.** All five routes are already implemented and validated in `texture/kearns.py` (1667
lines, `tests/unit/test_kearns_parameter.py`). Nothing scientific is missing; only the application
surface is. Highest impact per hour in the whole program.

### Sub-milestones

| Step | Content | State |
| --- | --- | --- |
| **M1a** | The self-contained routes: `kearns.from_orientations` (model texture), `kearns.from_tilt_profile`, `kearns.from_diffractogram`; the panel, the triad readout, the tilt-profile plot | **Complete** (2026-08-29) |
| **M1b** | The upload routes: `kearns.from_pole_figure`, `kearns.from_odf`, reusing the XRDML machinery already behind `texture.measured_pole_figures` | Not started |
| **M1c** | Worked example, theory-note cross-links, Playwright coverage, changelog | Not started |

### Design decisions taken

1. **A new panel `kearns`, grouped into the texture workspace.** `main.js` already supports grouped
   workspaces (`TEM_ANALYSIS`, `EBSD_ANALYSIS`); the texture workspace becomes the third, with
   `texture` and `kearns` as sub-tabs. This is what the user asked for — "one more tab in the
   texture module" — and it avoids a fourth view crammed into the existing texture view picker.
2. **The panel is manifest-driven, in the `calculator.js` style**, with one panel-specific view: the
   triad readout, the orientation tensor, and the tilt profile. No crystallography in JavaScript.
3. **Every operation must run from its own defaults** (enforced by
   `tests/unit/test_app_manifest.py::test_an_operation_runs_from_its_own_defaults`). The
   self-contained routes therefore ship a realistic default dataset; the upload routes declare a
   required parameter with no default and are exempt, as `texture.measured_pole_figures` already is.
4. **The default diffractogram is synthetic and says so.** It is generated from a known basal-fibre
   texture by the recipe in `tests/unit/test_kearns_parameter.py::synthetic_reflections`, so the
   answer is known before the calculation runs — the texture panel's own stated design principle.
   A service test pins the recovered `f` against the orientations route.
5. **The triad closure check is the panel's headline diagnostic.** `f_RD + f_TD + f_ND = 1` is exact
   for any texture, so a departure measures the systematic error of the measurement. The panel shows
   the sum whenever the directions form a triad, and says what a departure means.

### Next concrete step

**M1b — the upload routes.** Add `kearns.from_pole_figure` and `kearns.from_odf`, reusing the XRDML
machinery already behind `texture.measured_pole_figures`:

- `_measured_pole_figures` in `src/pytex/app/services/texture.py` shows the pattern: an
  `ObjectParameter(name="files", required=True)` with no default (which exempts it from
  `test_an_operation_runs_from_its_own_defaults`), then `uploaded_file(...)` and
  `read_xrdml_pole_figure(...)`, then `measurement.to_pole_figure(plane, specimen_frame=...)`.
- The resulting `PoleFigure` goes straight into `pytex.texture.kearns.kearns_from_pole_figure`.
  Note its `sampling` attribute decides the quadrature weights, and getting it wrong over-counts
  the pole of a tilt raster by up to 50 percent — the function reads it from the figure rather than
  assuming, so nothing extra is needed here, but the panel must report which was used.
- For `kearns.from_odf`, `_measured_odf` in the same module already inverts a set of pole figures
  into an `ODF`; `kearns_from_odf` takes it directly. Its `deconvolve_kernel` flag is a real
  scientific choice (support tensor versus density tensor) and must be a declared parameter, not a
  hidden default.
- Real fixture data exists at `kearns_parameter_data_references/reference_exp_data/` — three
  specimens with `.xrdml` files, including a PHWR clad tube with 002/011/012/013 figures. **These
  are untracked reference data, not repository assets**; check the cardinal repository-content rule
  before adding any of them to git.

The pole-figure route is *not* usable in a section where the basal peak is negligible (the ND-TD
section of strongly basal-textured zirconium), where the normalisation divides by noise. The panel
must say so where the user is, not only in the docstring — that failure is exactly what Mani Krishna
et al. (2011) traced the inconsistent literature f_RD values to.

---

## 3. History

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
