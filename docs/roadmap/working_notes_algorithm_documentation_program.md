# Working Notes: Algorithm Documentation Program (TD)

Running ledger for the **TD** program: making the algorithms, mathematics and constraints
behind every delivered feature readable in the Sphinx site, with publication-grade SVG
diagrams.

This file exists so an interrupted session can resume without reconstructing context from
chat history (AGENTS.md, "Durable progress and resumability"). Keep it current *before*
every long verification run and every commit.

## Objective

For every feature of the TX program (and the OR/diffraction surfaces it builds on):

1. **The algorithm**, stated as steps a reader could reimplement.
2. **The mathematics**, rendered in the Sphinx site rather than hidden in `.tex` sources —
   MyST has `dollarmath` and `amsmath` enabled, so equations belong on the page.
3. **The constraints**: assumptions, tolerances and what each one is calibrated against,
   failure modes, complexity, and what the surface deliberately does *not* do.
4. **A publication-grade SVG** per algorithm — flow sheet or geometry — following
   `docs/standards/visualization_style_guide.md`.
5. **Two worked systems throughout**: one **cubic** (Kurdjumov-Sachs, fcc→bcc) and one
   **hexagonal** (Burgers, bcc→hcp), so hexagonal-specific behaviour (four-index
   Miller-Bravais labelling, the non-cubic metric) is never left implicit.

## Ground rules

- Everything on `main`; commit and **push** at each satisfactory phase (user instruction).
- Figures are **generated**, not hand-authored, so a diagram cannot drift from the model it
  illustrates — the precedent set by `scripts/generate_reference_frame_figures.py`.
- Numbers quoted in the algorithm pages must be computed by a test, a worked example, or the
  figure generator. No hand-transcribed values.
- Reuse before invention: the SVG primitives already exist in
  `pytex.plotting.frame_diagrams`; extract rather than duplicate them.

## Phase status

| Phase | Scope | Status | Commit |
| --- | --- | --- | --- |
| TD0 | Audit, ledger, page plan | DONE | (TD1 commit) |
| TD1 | Shared SVG primitives + `algorithm_diagrams` module + generator | DONE | (this commit) |
| TD2 | OR determination from measured orientations (TX1) | DONE | (this commit) |
| TD3 | Variant correspondence tables (TX2) | DONE | (this commit) |
| TD4 | Composite SAED assembly and anchoring (TX3, TX4) | DONE | (this commit) |
| TD5 | SAED pattern indexing (TX5) | DONE | (this commit) |
| TD6 | LaTeX notes, index wiring, validation matrix, closure | DONE | (this commit) |

## Audit at TD0 (verified against the live repository, 2026-08-04)

**What already documents these features.** Prose exists and is good, but it explains *what*
the surfaces do and *when* to use them; it does not state the mathematics or the complexity.

- `docs/site/concepts/orientation_relationships.md` — concept-level, with the seeding
  subtlety and the preference rationale for TX1/TX2.
- `docs/site/workflows/composite_or_diffraction.md` — TX3/TX4 usage, the centring trap, the
  spot-order note.
- `docs/site/workflows/saed_pattern_solving.md` — TX5 usage, the honesty properties, limits.
- Docstrings carry purpose/when/inputs/outputs per the documentation contract.
- `docs/site/examples/generated/` — verified numbers, but one identity per example.

**What is missing**, and is therefore this program's scope:

- No rendered equations anywhere for the TX surfaces: the double-coset seed, the eigen-mean,
  the excitation-error selection, the ratio/angle admissibility test, the crystal-to-pattern
  triad construction.
- No complexity or vectorization statements, though the code is written for them.
- No single place that lists each algorithm's constraints and failure modes together.
- No SVG for any TX algorithm. The existing figure set covers frames, orientation reduction,
  IPF sectors, diffraction geometry and Ewald construction — nothing for OR determination,
  variant correspondence, composite assembly, or pattern indexing.

**Infrastructure that exists and will be reused.**

- `pytex.plotting.frame_diagrams` — canonical tokens, `text_width` estimation, `_card`,
  `_text`, `_relationship_arrow`, `_callout`, `_document` (emits `<title>`/`<desc>`),
  `_plain_marker` (absolute marker units). All private; TD1 promotes them.
- `tests/unit/test_figure_markers.py` — guards marker units, arrowhead scale, well-formed
  XML, and, for the *generated* figure list, title/desc, canonical font, and text
  overflow/collision. New generated figures must join that list.
- `scripts/audit_figure_text_layout.py` — the overflow/collision checker behind those tests.

## Page plan

New Sphinx section `docs/site/algorithms/`, one page per algorithm family, each with the
five elements above:

| Page | Covers | Figure |
| --- | --- | --- |
| `orientation_relationship_determination.md` | TX1: double-coset seed, symmetry-aware eigen-mean, catalog ranking, parallelism extraction | flow sheet |
| `variant_correspondence.md` | TX2: index maps, rationalization, family grouping | flow sheet + index-map geometry |
| `composite_saed_assembly.md` | TX3/TX4: excitation-error selection, shared detector basis, child anchoring, centring | flow sheet |
| `saed_pattern_indexing.md` | TX5: calibration, ratio/angle seeding, triad construction, verification, ambiguity | flow sheet |

## Verification command set

```
python -m pytest
python scripts/check_repo_integrity.py
python scripts/generate_algorithm_figures.py
python -m ruff check .
python -m mypy src
python -m sphinx -b html docs/site docs/_build/html
```

## Ledger

### TD0 (2026-08-04) — audit and plan

- Audited the delivered TX surfaces against the existing documentation; the gap is
  mathematics, constraints and complexity, not usage prose.
- Chose **generated** figures over hand-authored ones, following the reference-frame
  precedent, so the diagrams cannot drift and the existing layout guards apply.
- Fixed the two teaching systems for the whole program: Kurdjumov-Sachs (cubic) and Burgers
  (hexagonal product).

### TD1 (2026-08-04) — figure infrastructure

- **`pytex.plotting.svg_primitives`** now holds the shared elements: the style-guide tokens,
  `text_width`, a new `wrap_text`, `card`, `relationship_arrow`, `callout`, `header_width`,
  `arrow_marker` (always `userSpaceOnUse`) and `document` (always emits `<title>`/`<desc>`).
  `frame_diagrams` imports them instead of defining its own copies; its generated figures are
  byte-identical afterwards, which is the evidence the extraction is behaviour-preserving.
- **`pytex.plotting.algorithm_diagrams`** adds `AlgorithmStage`, `SideNote` and
  `algorithm_flow_svg`: a lane-based flow sheet, one row per phase, with constraint notes
  attached beside the stage they govern rather than collected in a legend. Role colour never
  carries information the label does not also carry.
- **`scripts/generate_algorithm_figures.py`** writes the four algorithm figures.

**A reproducibility defect found and fixed.** `pytex.plotting.frames` derived its arrowhead
marker ids from `builtins.hash(frame.name)`, and Python randomizes string hashing per
process — so regenerating a committed figure changed its bytes every run and `git diff` was
permanently dirty. A generated asset whose bytes move for no reason cannot be checked for
drift, which is the entire reason these figures are generated rather than drawn. Now a
`zlib.crc32` digest. `test_algorithm_figures_are_deterministic` pins the property for the new
figures.

**A layout defect caught by the existing guard**, on the first generation: a formula string
overflowed its card by 7 units. Rather than shortening that one string, `_stage_card` now
wraps formulas through `wrap_text` like every other label, so the class of defect cannot
recur. The figure guard tests gained `ALGORITHM_FIGURES`, so the new figures are held to the
strict overflow/collision contract rather than only the marker check.

### TD2a (2026-08-04) — a defect found while gathering the numbers to document

Documenting a quantity means first checking it. Computing the fingerprint sizes the
reconstruction robustness study quotes showed **both were wrong**, and the reason was a real
defect rather than a typo.

- `intervariant_boundary_fingerprint` deduplicated on **quaternions**, which need a sign
  convention. The convention was "make the largest-magnitude component positive", and two
  components tie in magnitude for the 90 and 180 degree elements of a crystal point group, so
  `argmax` broke the tie arbitrarily: numerically identical rotations canonicalized to `q` and
  `-q` and were counted twice — 81 of them for Kurdjumov-Sachs. Rounding the keys added an
  independent failure at rounding boundaries, which made the count of a group-theoretic set
  depend on the lattice parameters that entered the rotation (10664 vs 10665).
- Deduplication now runs on the **matrices**, which have no sign ambiguity, via a SciPy
  spatial query. A lexicographic sort was tried first and is wrong here: duplicates are not
  reliably adjacent, because a distinct element can agree with them in the leading entries and
  sort between them. Measured before believing it.
- **Not a results change** — the distance kernel takes a maximum over the set, so duplicates
  never altered a distance, a grouping, or an identification. It mattered because the size is
  quoted as science.
- True counts: **10 584** (Kurdjumov-Sachs), **684** (Burgers) — a factor of 15, which is the
  mechanism behind Burgers reconstructing more robustly. Three tests now pin the counts, their
  independence from lattice parameters and c/a, and the absence of duplicates.
- Corrected in `docs/testing/reconstruction_robustness_study.md`, the archived reconstruction
  ledger, and the CHANGELOG.

### TD2-TD5 (2026-08-04) — the four algorithm pages

New Sphinx section `docs/site/algorithms/`, wired into the site toctree with its own
caption, holding one page per algorithm family. Each page carries the same five elements:
what is computed, the mathematics rendered on the page, the algorithm as reimplementable
steps, the constraints and failure modes as admonitions beside the step they govern, the
cost, and a verification table naming where each claim is checked.

| Page | Covers |
| --- | --- |
| `orientation_relationship_determination` | double-coset seed, symmetry-aware eigen-mean, catalog ranking, parallelism extraction |
| `variant_correspondence` | the two index maps and why they differ, rationalization, family grouping |
| `composite_saed_assembly` | Ewald/excitation error, shared detector basis, child anchoring, centring |
| `saed_pattern_indexing` | calibration, ratio/angle seeding, triad construction, symmetry deduplication, ambiguity |

**Every number on these pages was computed, not recalled**, and several are new to the
documentation: the symmetry-reduced catalog separation matrix (5x5, fcc->bcc), the
equivalence-group breakdown of the KS `(111)` table (4 groups of 6) and the Burgers
`(011)beta` table (2+2+4+4, uneven because the hexagonal child group is smaller), the
per-variant child-zone deviations for Burgers down beta [110] (two variants exactly on
zone at 54 reflections, the rest at 1.0-1.3 deg), and the solver's noise envelope
(100% at 0.5 mm, no solution at 2.0 mm).

**Cubic and hexagonal side by side on every page**, as required: Kurdjumov-Sachs and
Burgers, with the hexagonal-specific behaviour — four-index Miller-Bravais labelling, the
non-cubic metric making the two index maps genuinely differ, the smaller child point group
— stated rather than left implicit.

Sphinx builds zero-warning with all four pages and their figures.

### TD6 (2026-08-04) — theory notes, wiring, closure

- **Two canonical LaTeX notes**, both promised by the TX specification and outstanding until
  now: `docs/site/theory/orientation_relationship_determination.md` and
  `docs/site/theory/saed_ratio_angle_indexing.md`. They carry the derivations and the
  convention fixing; the Sphinx pages carry the same algorithms with worked numbers and
  figures, and each side points at the other.
- **Wiring**: the theory index, `docs/README.md` (LaTeX tree, the four new figures, and the
  TD ledger), and three rows in the plotting validation matrix covering the generated
  algorithm figures, generation determinism, and the shared SVG primitives.
- **CHANGELOG** entries under `Added` (the algorithms section, the theory notes, the figure
  modules) and `Fixed` (the non-reproducible marker ids).

## Program outcome

Every delivered TX feature now has its algorithm, mathematics and constraints readable in the
Sphinx site, with a generated publication-grade figure, and cubic and hexagonal examples side
by side.

| Feature | Algorithm page | Figure | LaTeX note |
| --- | --- | --- | --- |
| OR determination (TX1) | `orientation_relationship_determination` | yes | new |
| Variant correspondence (TX2) | `variant_correspondence` | yes | index-correspondence note (existing) |
| Composite SAED (TX3, TX4) | `composite_saed_assembly` | yes | kinematic-spot note (existing) |
| Pattern indexing (TX5) | `saed_pattern_indexing` | yes | new |

**Two defects were found by doing this work**, both because documenting a quantity means first
computing it:

1. **Generated figures were not reproducible** — arrowhead marker ids came from
   `builtins.hash`, which Python randomizes per process, so a committed figure's bytes moved
   on every regeneration and drift could not be detected. Fixed with a stable digest and
   pinned by a test.
2. **The intervariant fingerprint contained duplicates and its size moved with lattice
   parameters** — quaternion sign canonicalization broke on magnitude ties. Fixed by
   deduplicating on matrices; two published numbers corrected; three tests added. (TD2a.)

Both are recorded in the CHANGELOG.

### Open follow-ons

- **Geometry figures**, as opposed to flow sheets: what an excitation error *is* on the Ewald
  sphere, what the crystal-to-pattern triad construction looks like. The existing
  `zone_axis_ewald_geometry.svg` covers part of the first. These would suit a dedicated
  renderer rather than `algorithm_flow_svg`.
- **Algorithm pages for the surfaces outside the TX program** — parent-grain reconstruction,
  OR identification from boundaries, harmonic ODF inversion, EBSD grain segmentation. The
  section and its conventions are now in place for them.
- The pre-existing ~16 ruff findings (RUF022/RUF043/RUF059) in files this program did not
  touch, and the two mypy `to_hex` stub-drift errors in `plotting/crystal3d.py`.
