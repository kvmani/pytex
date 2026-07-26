# Working Notes — Crystallographic Convention Harmonization

Durable plan and phase ledger for making PyTex's human-facing notation fully consistent with
established international crystallographic convention, and for fixing those conventions centrally
so they cannot drift. Read this file first to resume; continue at the first phase not marked
`DONE`.

This program follows directly from the reference-frame foundation
(`working_notes_reference_frame_foundation.md`), which introduced starred reciprocal axis labels
and exposed that PyTex had no central rule for them.

## Objective

1. **Reciprocal quantities carry the star, everywhere.** Reciprocal basis vectors and reciprocal
   axis labels render as `a*, b*, c*` in every surface — frames, bases, plots, prose, docs — and
   the rule is a foundational principle, not a local habit.
2. **Codify the other international conventions PyTex already relies on but never fixed
   centrally**, each anchored to a normative source per `docs/standards/reference_canon.md`.
3. **Harmonize the repository** so code, docstrings, docs, figures, notebooks, and worked examples
   all speak the same notation.
4. **Use the new reference-frame visualizations consistently** wherever a frame is at stake in the
   documentation.

## What "starred" does and does not mean

This is the part most easily got wrong, so it is fixed explicitly. Following IUCr convention
(International Tables Volume A; Volume C §1.1 for reciprocal-space definitions):

| Quantity | Correct form | Starred? |
| --- | --- | --- |
| direct basis vectors | **a**, **b**, **c** | no |
| **reciprocal basis vectors** | **a\***, **b\***, **c\*** | **yes** |
| reciprocal-space axis labels of a frame | `a*`, `b*`, `c*` | **yes** |
| reciprocal lattice parameters | `a*`, `b*`, `c*`, `alpha*`, `beta*`, `gamma*` | **yes** |
| Miller plane indices | `(hkl)` | **no** |
| plane family | `{hkl}` | **no** |
| direction indices | `[uvw]` | **no** |
| direction family | `<uvw>` | **no** |
| reciprocal lattice vector | **g**\_hkl = h**a\*** + k**b\*** + l**c\*** | bold g, starred basis |

The star marks the **basis**, not the indices. `(hkl)` are already reciprocal-basis components by
definition, so starring them would be wrong and would read as a different quantity. A stable
surface that quotes reciprocal-basis components must instead *say* the basis is reciprocal — which
is what the starred axis labels and `BasisKind.RECIPROCAL` accomplish.

## Additional conventions to fix centrally (HC4)

Each of these is something the repository already assumes in prose or code but has never fixed in
one authoritative place.

1. **Bracket families.** `(hkl)` a specific plane; `{hkl}` a symmetry-related family of planes;
   `[uvw]` a specific direction; `<uvw>` a symmetry-related family of directions. PyTex had no
   formatter for the family forms even though `properties/slip.py` and the multiplicity surfaces
   discuss them.
2. **Negative indices use an overbar**, not a leading minus, in publication-facing output
   (`1bar` rendered as an overbar over the digit).
3. **Four-index Miller-Bravais forms** `(hkil)`, `{hkil}`, `[uvtw]`, `<uvtw>` with the redundancy
   constraints `h + k + i = 0` and `u + v + t = 0`.
4. **Units at stable boundaries**: angstrom for direct lengths, inverse angstrom for reciprocal
   lengths, degrees for angles on the public API, radians only where a function name says `_rad`.
5. **Zone law** `h u + k v + l w = 0` as the statement relating a plane and a direction lying in
   it.
6. **Hermann-Mauguin symbols** as the primary point-group and space-group naming form.

## Phase ledger

Status values: `TODO`, `IN PROGRESS`, `DONE`.

| Phase | Scope | Status |
| --- | --- | --- |
| HC1 | This plan; convention inventory | DONE |
| HC2 | Starred-reciprocal + family-bracket notation surface in `core/notation.py`, used everywhere | DONE |
| HC3 | Codify the reciprocal-star rule in the foundational docs, with an enforcing policy test | DONE |
| HC4 | Codify the further international conventions listed above | DONE |
| HC5 | Roll the new reference-frame visualizations through the docs consistently | DONE |
| HC6 | Full verification, then commit and push to `main` | DONE |

## Verification record

- Entering state (reference-frame foundation complete, uncommitted): ruff clean, mypy strict
  clean, **1191 tests passed**, coverage 89.31%, integrity and Sphinx build green.
- HC6 final run (all green):
  - `python -m ruff check .` - All checks passed
  - `python -m mypy src` - no issues in 86 source files (strict)
  - `python scripts/check_repo_integrity.py` - passed
  - `python -m sphinx -b html docs/site docs/_build/html` - build succeeded with **zero warnings**
    (the 5 pre-existing `reference_canon.md` cross-reference warnings were fixed as part of the
    harmonization)
  - `python -m pytest` - **1236 passed** (1191 entering this program; +45), no warnings
  - coverage **89.37%**, above the 87% CI gate
  - worked-example gallery regenerated; all 21 notebooks re-executed, and the ones whose
    diff was only an execution timestamp were reverted so the commit shows real changes
  - committed and pushed to `main` as `73f4dc8`

## Outcomes

### New notation surface (`pytex.core.notation`)

- **Reciprocal stars**: `RECIPROCAL_STAR`, `format_reciprocal_axis_label`,
  `format_reciprocal_axis_labels`, `format_reciprocal_lattice_vector`,
  `is_reciprocal_axis_label`, `strip_reciprocal_star`. Starring is **idempotent**, so a label
  passing through two layers cannot become `a**`.
- **Bracket families**: `format_miller_indices` gained a `scope` parameter, plus
  `format_plane_family_indices` / `format_direction_family_indices` for `{hkl}` and `<uvw>`.
  Mathtext emits the escaped `\{...\}` and `\langle ... \rangle` forms, verified to parse in
  matplotlib.

### Ambiguity fix found by an existing test

Replacing an inline formatter with the central one made a pinned test fail, which exposed a real
defect: `format_miller_indices` concatenated components unconditionally, so `[1-10]` could be read
as `[1, -1, 0]` **or** `[1, -10]`, and `(1210)` as `(1, 2, 1, 0)` or `(12, 1, 0)`. The formatter now
inserts a separator whenever a component is negative in plain style or any component has more than
one digit; single-digit non-negative indices keep the classical concatenated form `(110)`. The
composite-SAED module had already solved this locally — the central formatter now matches it.

### Scientific corrections, not just cosmetics

- **Pole figures and powder reflections are families.** A pole figure plots the whole symmetry
  orbit and a powder reflection *is* its multiplicity, so both now read `{hkl}` rather than
  `(hkl)`. Because a `PoleFigure` can be built without family expansion, the object now records
  `includes_symmetry_family` and the notation follows the record rather than an assumption.
  Contracts round-trip the field with a backward-compatible default.
- **Reciprocal frames** built anywhere (including `Lattice.reciprocal_basis`) carry `a*, b*, c*`.
- Five modules that formatted indices inline now route through `pytex.core.notation`.

### Enforcement

`tests/unit/test_notation_conventions.py` (43 tests) pins the conventions *and* enforces them:
one test fails if any module reintroduces inline index formatting, another if any module appends a
reciprocal star by hand, and a third checks the rules are actually documented in `AGENTS.md`, the
notation standard, and the terminology registry. Every mathtext form is rendered through matplotlib
so an unparseable label fails as a test rather than as a broken figure.

### Documentation

- The starred-reciprocal rule is a non-negotiable in `AGENTS.md`, a section of
  `notation_and_conventions.md` (with the explicit starred/not-starred table), entries in the
  terminology registry, a topic in `reference_canon.md` classified
  `normative from IUCr/International Tables`, and a rule in `canonical_data_model.md`.
- The generated frame figures now appear on the site index, `core_model`, `core_foundation`, the
  technical glossary (which gained an index-notation table and the sample frame), and the
  diffraction-geometry, SAED, composite-OR, crystal-visualization, plotting-primitives, and
  EBSD-import workflow pages — each documenting the capability that page actually needs rather
  than repeating one picture.
- Notebook 01 gained sections on the reciprocal star and the bracket families; all 21 notebooks
  were re-executed so no stored output shows the old notation.
- Two new worked examples pin the convention with computed values: the reciprocal frame has
  exactly three starred axes and the crystal frame none; every Miller bracket form contains zero
  stars.
