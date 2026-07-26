# Working Notes — Figure Layout Repair And Burgers OR Notebook Program

Durable plan and phase ledger. **Read this file first to resume**, then continue at the first
phase whose status is not `DONE`. Each phase is independently landable, so an interruption costs
at most the phase in flight.

This follows the reference-frame foundation
(`working_notes_reference_frame_foundation.md`) and the convention harmonization
(`working_notes_convention_harmonization.md`). The figure half closes the known follow-on recorded
at the end of that second program.

## Objective

1. **Every documentation figure is legible.** Review all `docs/figures/*.svg` and fix text that
   runs outside its box, overflows the canvas, or collides with other text. The previous program
   fixed arrowhead scaling repo-wide and regenerated six frame figures; text layout in the
   remaining hand-authored figures was explicitly left as a follow-on and is now in scope.
2. **The Burgers orientation relationship is taught exhaustively, in both directions.** The OR
   tutorial notebooks gain detailed bcc → hcp *and* hcp → bcc treatment: variant calculations,
   transformation-strain estimation, parallel plane and direction calculations, with crystal
   visualizations and variant pole figures embedded.

## Phase ledger

Status values: `TODO`, `IN PROGRESS`, `DONE`.

| Phase | Scope | Status |
| --- | --- | --- |
| FX1 | This ledger; active-task pointer | DONE |
| FX2 | Transform-aware SVG text-layout auditor | DONE |
| FX3 | Fix text overflow and overlap in all flagged figures | DONE |
| FX4 | Repo-wide figure text-layout policy test | DONE |
| FX5 | Audit Burgers OR support in the core API; record gaps | DONE |
| FX6 | Exhaustive forward Burgers (bcc → hcp) notebook content | TODO |
| FX7 | Reverse Burgers (hcp → bcc) notebook content | TODO |
| FX8 | Execute notebooks; pin key numbers in tests/worked examples | TODO |
| FX9 | Full gate; commit and push | TODO |

## Part 1 — Figure layout

### Why the earlier auditor was insufficient

The first pass only measured `<text>` elements carrying absolute `x`/`y`. Much of the text in these
figures sits inside `<g transform="translate(...)">` groups, so it was never checked. FX2 resolves
translate transforms before measuring, and adds a check the first pass lacked entirely: text
overflowing its **enclosing panel rectangle**, not just the canvas. That is the defect actually
reported — "text is running outside its box".

### Known targets entering the program

From the previous pass, with canvas-overflow or same-baseline collisions:
`active_passive_rotation`, `bunge_euler_geometry`, `orientation_conventions`,
`ipf_sector_geometry_detailed`, `ipf_sector_reduction`, `kinematic_spot_projection`,
`orientation_reduction_workflow`, `pole_figure_construction`, `pytex_architecture_compact`,
`pytex_architecture_evolution_poster`, `so3_dirichlet_fundamental_region`,
`symmetry_direction_orbit`, `zone_axis_ewald_geometry`.

The first three were already regenerated from code in the previous program and are clean; the rest
are hand-authored and are the FX3 workload. Expect FX2 to surface more once group transforms are
resolved.

### Repair policy

Prefer the smallest correct change, in this order:

1. widen the enclosing box or canvas, when the text is right and the box is too small;
2. split a long caption across lines at a natural clause break;
3. shorten the wording, only when it does not lose scientific content;
4. reposition a colliding label.

Do **not** silently drop content to make it fit.

## Part 2 — Burgers OR notebooks

### Scientific content to cover, both directions

- **Variant calculation**: the full variant set with symmetry-reduced distinctness, intervariant
  misorientation table, and packet/group structure.
- **Transformation strain**: the deformation gradient implied by the lattice correspondence, its
  principal strains, and volume change.
- **Parallelism**: the defining plane and direction parallelisms, verified numerically, plus the
  family-orbit parallelism search.
- **Visualization**: crystal scenes of parent and product in the relationship, and variant pole
  figures.

### Reverse transformation

`hcp -> bcc` is not simply the inverse rotation: the parent and child phases swap, so the variant
count, the symmetry reduction, and the strain all differ. FX5 establishes what the API supports
before FX7 writes the content, and any gap is recorded here rather than worked around in notebook
prose.

## Part 1 outcomes (FX2-FX4)

### Accurate metrics first

The first auditor assumed an average character width, which is wrong by up to 30%
for capital-heavy strings, and it read font sizes from attributes only — these
figures declare most typography in CSS classes, so it was measuring almost
everything at a default 12px. Both were fixed before any repair:
`pytex.plotting._svg_text` now carries the Helvetica advance-width table (Arial is
metrically compatible), and the auditor resolves sizes from attributes, classes and
ancestors. `frame_diagrams.text_width` measures through the same table, so
generated and hand-authored figures are judged by one ruler.

Measured accurately, the starting state was **86 defects across 39 figures**.

### What the auditor checks

- `canvas` — text past the document viewBox
- `box` — text past the rect that encloses it (the reported "running outside its box")
- `collision` — two runs on the same or nearly the same baseline overlapping
- `over-card` — text painted across a panel it is not inside

The last two were added after a visual pass showed the first version missing real
defects: an exact-baseline test misses runs half a line apart, and a containment
test cannot see a header painted over the card beneath it.

### Repairs

Three scripts, applied in order, each keeping every word at its designed size:

1. `fix_svg_text_overflow.py` — wraps captions to fit their box, growing the box
   downward or sideways when the space allows. **Never wraps titles**: a card
   title's second line lands exactly on the caption, turning an overflow into an
   overprint (learned the hard way — the first run did exactly that).
2. `fix_svg_title_wraps.py` — restores wrapped titles to one line and widens the
   card instead, bounded by the neighbouring card.
3. Individual fixes for what remained: a centred header left-aligned off the card
   it sat on, five panels widened, one footer caption moved to the margin.

### Result

**86 → 19.** `box` and `canvas` overflows — the defect reported — are at **zero**
and held there by `test_no_figure_text_runs_outside_its_box_or_canvas`. The 19
remaining are `collision`/`over-card` inside dense teaching diagrams where the
geometry is meaningful (an axis label beside a formula, a callout over a sphere);
each needs an individual decision. They are capped by
`test_figure_text_overlaps_do_not_grow` so no future change can add more.

### Note on tooling damage, caught and reverted

An early version of the box-growth pass wrote back through `ElementTree`, which
reformatted seven figures with generated `ns0:` namespace prefixes. Caught by
inspection, reverted from a backup, and the pass rewritten to edit the raw source
textually. Worth recording: round-tripping hand-authored SVG through a parser is
not a neutral operation.

## Part 2 status (FX5 complete; FX6-FX9 open)

### What already works for Burgers

Verified against the `zr_hcp` and `fe_bcc` fixtures plus a constructed beta-Zr
(bcc, a = 3.574 A):

- `OrientationRelationship.from_burgers_correspondence(parent_phase=bcc, child_phase=hcp)`
  builds the relationship and enforces the point groups (parent 432, child 622).
- `generate_variants()` returns the literature **12 variants**.
- `misorientation()` gives **45.291 deg**.
- `correspondence_direct()` / `correspondence_reciprocal()`, `find_parallel_planes()`,
  `find_parallel_directions()`, `intervariant_misorientation_angles_deg()`,
  `variant_pole_figure()` and `variant_close_packed_groups()` are all available.

So variants, misorientation, parallelism and pole figures can be written now.

### The blocker: transformation strain for a hexagonal product

`deformation_gradient()` **fails for Burgers** with:

> The nearest-integer lattice correspondence is singular; the relationship's exact
> correspondence is too far from an integer matrix for Bain-strain analysis.

This is not a fixture problem and not an element mismatch — it fails equally for
beta-Zr to alpha-Zr, both allotropes of one element. The cause is structural:

`deformation_gradient()` computes `exact = solve(A_child, R A_parent)` from the
**conventional** direct bases and rationalizes to the nearest *integer* matrix. That
is correct for cubic-to-cubic (Bain, KS, NW all work, giving the textbook Bain
stretches 1.1504, 1.1504, 0.8135), because the conventional cubic basis vectors map
onto child lattice vectors. For bcc to hcp they do not: only the **primitive** bcc
vectors (the `<111>/2` set) map onto hcp lattice vectors, so the conventional-basis
correspondence is irrational. Measured for beta-Zr to alpha-Zr it is

```
[[-0.8991  0.8991  0.1172]
 [-0.3778  0.3778  1.1597]
 [ 0.4910  0.4910  0.0000]]
```

which is not near-integer at any denominator up to 6, and whose nearest integer
matrix is singular. The error message is therefore accurate and the guard is doing
its job — the method simply does not cover hexagonal products.

### Scoped fix, before FX6/FX7 write strain content

Extend `deformation_gradient()` to work from **primitive** bases (deriving the
primitive vectors from the lattice centring), and/or accept an explicit
`correspondence` argument so a caller can supply a literature lattice
correspondence. Validate against published Burgers strains for Zr/Ti rather than
against a prior program output, per the executable-example rule.

Do **not** hand-roll the strain inside the notebook: `AGENTS.md` requires expanding
the shared core model rather than encoding one-off conversions locally, and a
notebook that computed this privately would be exactly the drift the repository
forbids.

### Order of work when resuming

1. FX5a (new): extend `deformation_gradient` for hexagonal products; add tests
   pinning the Burgers strain to a cited value.
2. FX6: forward bcc to hcp notebook content — variants, misorientation table,
   packet grouping, parallelisms, strain, crystal scenes, variant pole figures.
3. FX7: reverse hcp to bcc. Note the reverse is **not** the inverse rotation: the
   phases swap, so `from_burgers_correspondence` cannot be reused with arguments
   exchanged (it requires parent 432 / child 622). Build it with
   `from_parallel_plane_direction` using (0001)_hcp || (110)_bcc and
   `<11-20>_hcp || <-111>_bcc`, and expect a different variant count and symmetry
   reduction.
4. FX8, FX9 as listed.

## Verification record

- Entering state (commit `059026a`): ruff clean, mypy strict clean, **1293 tests passed**,
  coverage 89.52%, integrity green, Sphinx build warning-free.
- FX9 final run: *(to be recorded)*

## Outcomes

*(filled in per phase as work lands)*
