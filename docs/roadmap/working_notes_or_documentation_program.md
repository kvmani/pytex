# Working Notes: OR Scientific Documentation Program (2026-07-18)

> **Policy note (2026-08-08).** This note predates the change to notebook handling.
> Notebooks are now committed *without* outputs and executed by the Sphinx build
> (`nb_execution_mode = "cache"`); `scripts/execute_notebooks.py` has been removed.
> Instructions below to run it are historical record, not current procedure.

Purpose: resumable progress ledger for the orientation-relationship (OR)
documentation program. The goal is documentation that works simultaneously as
scientific reference (equations, conventions, algorithms), user tutorial
(runnable notebooks with rich rendered outputs), and developer documentation
(algorithm internals, validation pointers) — a standard-setting layer over the
implemented F1–F14 surface.

If interrupted, resume by reading this file top to bottom; each phase lists
status, deliverables, and exactly what remains. Master context:
`docs/development/active_task_progress.md` (phases 0–28 landed the OR feature
surface this program documents).

## Program Design

Three pillars, delivered phase by phase (one verified commit + push each):

1. **Scientific SVG diagram assets** (`docs/site/_static/or/`), hand-authored,
   colorful, consistent palette (parent = indigo `#3f51b5`, child = orange
   `#f57c00`, math/annotation = teal `#00897b`, neutral greys), reusable from
   notebooks, concept pages, and README-class material:
   - `or_doctrine_pipeline.svg` — rotation vs correspondence vs deformation
     doctrine with the report objects each produces.
   - `variant_generation_algorithm.svg` — flowchart: defining parallelism →
     rotation → parent-symmetry orbit → child-symmetry dedup → variants →
     consumers (intervariant table, packets, pole figures).
   - `bain_correspondence_cells.svg` — the classic two-fcc-cells /
     inscribed-bct-cell construction with contraction/expansion arrows.
   - `f7_identification_refinement_pipeline.svg` — EBSD map → grain graph →
     boundary misorientations → fingerprint scoring → refinement loop →
     parent reconstruction.
   - `burgers_family_basal_alignment.svg` — hcp basal hexagon on bcc {110}
     with the shared close-packed direction and the Burgers / Pitsch-Schrader
     (5.26 deg) / Potter (~1.4 deg) rotations about it.

2. **Executed tutorial notebooks** (at the time these were produced by a
   notebook generator, since removed — they are now hand-authored `.ipynb`
   files — then executed with `scripts/execute_notebooks.py` so the Sphinx
   site, which runs with `nb_execution_mode = "off"`, renders the committed
   rich outputs):
   - `18_orientation_relationships_fundamentals` — conventions and math
     (crystal→specimen orientations, $g_{\text{child}} = g_{\text{parent}} \circ \mathbf{V}^{\mathsf{T}}$, coset
     misorientation representative), KS as flagship: 24 variants, intervariant
     angle spectrum (Morito), packets 4×6, variant pole figure, deviation +
     fitting (GT recovered from KS nominal), index correspondence.
   - `19_lattice_correspondence_and_transformation_strain` — correspondence
     matrices `M = A_c^{-1} R A_p`, deformation gradients
     $\mathbf{F} = \mathbf{R}^{\mathsf{T}} \mathbf{A}_c \mathbf{M}_{\mathrm{int}} \mathbf{A}_p^{-1}$, Bain principal stretches analytically and
     live, KS/NW/GT polar rotations, 3D lattice-point renderings of the Bain
     cell-in-cell construction and the Burgers basal alignment.
   - `20_or_catalogs_identification_and_reconstruction` — live tour of all
     five standard catalogs with computed variant counts and separations
     (PS–Burgers 5.26 deg, Potter–Burgers ~1.37 deg, Bagaryatsky–Isaichev
     3.586 deg), then the experimental F7/F8 pipeline end to end on a
     synthetic lath-martensite microstructure: identification → rotation
     refinement → parent-grain reconstruction, with honesty statements.

3. **Cross-linking** — notebook atlas (`docs/site/tutorials/notebooks.md`),
   concept page (`docs/site/concepts/orientation_relationships.md`), ledger
   and CHANGELOG entries.

Conventions used throughout (must match `docs/site/concepts/*` and code):
orientations are crystal→specimen rotations; the canonical composition is
$g_{\text{child}} = g_{\text{parent}} \circ \mathbf{V}^{\mathsf{T}}$ (variant rotations map parent-crystal to
child-crystal frames); all reported inter-OR distances are symmetry-reduced.

## Phase Ledger

| Phase | Deliverable | Status | Commit |
| --- | --- | --- | --- |
| 29 | Program plan (this file) + 5 SVG diagram assets | done | `b125d40` |
| 30 | Notebook 18 (fundamentals, executed) + executor script | done | `5a308f5` |
| 31 | Notebook 19 (correspondence/strain, executed) | done | `89ae8ad` |
| 32 | Notebook 20 (catalogs/pipeline, executed) + final wiring | done, committing | — |

Program v1 complete: all three notebooks executed and committed, five SVG
assets in use, concept page cross-linked. Natural follow-ons (not scheduled):
a LaTeX theory note consolidating the notebook math, an interactive variant
explorer, and measured-data notebook sections once the reconstruction
fixtures land.

Verification gates per phase: `ruff check .`, `mypy src`,
`python scripts/check_repo_integrity.py`, full pytest, and
`python -m sphinx -b html docs/site docs/_build/html` when site content
changes. Notebook phases additionally re-run
`python scripts/execute_notebooks.py --only <new>` to prove execution is
clean. (Historical note: these phases also ran a notebook generator, since
removed — see the Resume Instructions below.)

## Resume Instructions

- SVGs are hand-authored files; edit directly, keep the palette above, keep
  text as `<text>` elements (searchable, translatable), 12–16 px labels.
- **Superseded (2026-07-20):** notebook content used to live in
  `scripts/generate_tutorial_notebooks.py`, under a policy that notebooks were
  generated and never hand-edited. That generator has been **removed**. The
  `.ipynb` files are now the source of truth and are edited directly; after
  editing, run `python scripts/execute_notebooks.py --only <prefix>` and commit
  the executed notebook. See
  `docs/roadmap/working_notes_composite_saed_program.md` for why.
- The smoke-test policy (`tests/unit/test_notebooks.py`) requires every
  notebook to appear in `notebooks.md` and contain markdown + code cells.
