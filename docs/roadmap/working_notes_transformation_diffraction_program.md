# Working Notes: Transformation Crystallography And Diffraction Program (TX)

Running ledger for the **TX** program. The normative specification is
[Transformation Crystallography And Composite Diffraction Program](../architecture/transformation_crystallography_and_diffraction_program.md).

This file exists so an interrupted session can resume without reconstructing context from chat
history (AGENTS.md, "Durable progress and resumability"). Keep it current *before* every long
verification run and every commit.

## Objective

Deliver five user-facing answers on top of existing PyTex primitives:

- **(a) TX1** — measured parent/child Euler angles ⇒ *what is the OR?*
- **(b) TX2** — an OR + an arbitrary parent plane/direction ⇒ the parallel planes/directions in
  every product variant.
- **(c) TX3** — a parent zone axis ⇒ a robust composite kinematic SAED, exportable as graphics
  **and** reflection tables.
- **(d) TX4** — a *product-variant* zone axis ⇒ the same composite, with matrix and sibling
  variants generated around it.
- **(e) TX5** — a measured SAED pattern ⇒ solved, from spots picked interactively or listed in a
  YAML file.

Closed out by **TX6**: a Burgers β↔α notebook demonstrating (a)–(e) end to end.

## Ground rules for this program

- Everything on `main`; no feature branches (user instruction).
- Commit after each phase's gates pass, so no work is ever at risk.
- Reuse before invention: no TX surface may re-derive a rotation convention, a symmetry reduction,
  a rationalization, or a detector basis that already exists in the core.
- Every new report object gets `describe()` and a JSON contract in lockstep.
- Every new numerical surface gets an executable worked example with **independent** provenance.

## Phase status

| Phase | Scope | Status | Commit |
| --- | --- | --- | --- |
| TX0 | Specification + ledger | DONE | (this commit) |
| TX1 | OR characterization from measured orientations | TODO | |
| TX2 | Variant correspondence tables | TODO | |
| TX3 | Composite SAED robustness + export layer | TODO | |
| TX4 | Child-zone-anchored composite patterns | TODO | |
| TX5a | Measured-pattern YAML + calibration + solver core | TODO | |
| TX5b | Variant assignment + interactive picker | TODO | |
| TX6 | Burgers notebook, figures, theory notes, ledger closure | TODO | |

## Baseline established at TX0 (verified against live code, 2026-08-03)

What already exists, so later phases do not rebuild it:

- `pytex.core.transformation`: `OrientationRelationship` with eleven correspondence constructors
  (Bain, NW, KS, GT, Pitsch, Shoji–Nishiyama, Burgers, Pitsch–Schrader, Potter, Bagaryatsky,
  Isaichev), `generate_variants()`, `map_{plane,direction}_to_{child,parent}`,
  `map_{plane,direction}_across_variants`, `find_parallel_{planes,directions}`,
  `fit_orientation_relationship`, `or_deviation`, `intervariant_boundary_fingerprint`,
  `boundary_fingerprint_distances_deg`, `deformation_gradient`, `variant_pole_figure`.
- `pytex.core.parent_reconstruction`: `OrientationRelationshipCatalog` and the
  `standard_{fcc_bcc,bcc_hcp,fcc_hcp,hcp_bcc,ferrite_cementite}_relationships` builders.
- `pytex.experimental`: `identify_orientation_relationship` (child–child boundaries),
  `refine_orientation_relationship_from_boundaries`, `reconstruct_parent_grains`.
- `pytex.diffraction.kinematic`: `simulate_zone_axis_spots`, `zone_basis_from_axis`, `SpotTable`,
  `KinematicSimulationConfig`, `centering_allowed_mask`, `electron_structure_factors`.
- `pytex.diffraction.composite`: `simulate_composite_saed`, `CompositeSAEDPattern`,
  `VariantZonePattern`, `rationalize_zone_axis`, `find_spot_coincidences`, `sweep_parent_zone_axes`.
- `pytex.diffraction.models`: `DiffractionGeometry`, `KinematicSimulation`, `index_saed_pattern`,
  `estimate_zone_axis` — a **detector-geometry-driven** indexing path, distinct from TX5's
  calibrated-spot-list path; TX5 must state the difference and reuse what it can.
- `pytex.plotting.composite_saed`: renderer, `CompositeSAEDPlotConfig`, annotation engine.

Identified gaps are enumerated per feature in the specification §3.2, §4.2, §5.2, §6, §7.

## Known pre-existing issues on this machine (inherited, not caused by TX)

Carried forward from the reconstruction-stabilization ledger; do **not** fix inside a TX
scientific-behavior commit:

1. ~~Six phase-fixture SHA-256 mismatches~~ — **FIXED in TX0.** The cause was confirmed, not
   guessed: `fixtures/phases/fe_bcc/phase.cif` held 19 CRLF pairs on disk and hashed to
   `e512334e…`, while the LF form hashes to `8afe4f95…` — exactly the digest pinned in
   `fixtures/phases/catalog.json`. Git's Windows default `core.autocrlf=true` was rewriting
   checksum-pinned artifacts on checkout. A `.gitattributes` marking
   `fixtures/phases/**`, `fixtures/mtex_parity/**` and `*.ipynb` as `-text` disables the
   conversion; `scripts/check_repo_integrity.py` now passes. **The full test suite no longer
   needs the two deselects.**
2. ~20 ruff findings from newer rule versions (RUF022/RUF059/RUF043) in untouched files.
3. Two mypy `to_hex` arg-type errors in `plotting/crystal3d.py` (matplotlib stub drift).
4. No MATLAB/MTEX on this machine, so no new MTEX parity claim can be executed here.

## Verification command set

```
python -m pytest
python scripts/check_repo_integrity.py
python -m ruff check .
python -m mypy src
python -m sphinx -b html docs/site docs/_build/html
python scripts/generate_worked_examples.py
```

## Ledger

### TX0 (2026-08-03) — specification and ledger

- Read the governing documents and audited the live code surface for all five asks; the baseline
  above is what that audit found, not an assumption.
- Wrote the normative specification with per-feature API signatures, algorithms, and validation
  plans, including the honest-limits requirements (`is_conclusive` semantics for TX1 and TX5, the
  `z` vs `-z` SAED ambiguity, kinematic-only intensities).
- Recorded the phase order and its dependencies: TX4 needs TX3's geometry, TX5b needs TX4.
- **Fixed the inherited fixture-hash gate failure** (see item 1 above) by adding `.gitattributes`.
  This was worth doing first: the integrity check is a gate for every TX phase, and leaving it red
  would mean every later phase runs with two deselects and cannot tell a new breakage from the old
  one.
- Artifact hygiene: `.gitignore` now covers `docs/site/_build/`, the regenerated
  `fixtures/mtex_parity/results/pytex/*/` directories, and the stray root `package.json`, all of
  which were sitting untracked in the worktree.

### Next action

Start **TX1**: refactor `fit_orientation_relationship`'s align/average loop into a shared
`_fit_from_seed`, add the seedless double-coset eigen-mean start, then build
`ORCharacterizationReport` with catalog ranking and parallelism extraction.
