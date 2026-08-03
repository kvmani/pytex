# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history.

## Current Goal — COMPLETE (2026-08-04)

**Transformation Crystallography And Composite Diffraction Program (TX)** — five user-facing
phase-transformation answers built on the existing OR and diffraction primitives. All phases
TX0–TX6 are delivered and committed to `main`.

| Ask | Delivered surface |
| --- | --- |
| **(a)** measured parent/child Euler angles ⇒ *what is the orientation relationship?* | `characterize_orientation_relationship`, `orientation_relationship_from_euler`, `describe_orientation_relationship` |
| **(b)** an OR + an arbitrary parent plane/direction ⇒ the parallel planes and directions in every product variant | `variant_correspondence_table` |
| **(c)** a parent zone axis ⇒ a composite kinematic SAED exportable as graphics **and** reflection tables | `composite_reflection_table`, `export_composite_saed`, `CompositeSAEDPattern.centering_audit` |
| **(d)** a product-variant zone axis ⇒ the same composite, matrix and siblings around it | `simulate_composite_saed_from_child_zone` |
| **(e)** a measured SAED pattern ⇒ solved, from spots picked interactively or listed in YAML | `solve_saed_pattern`, `MeasuredSAEDPattern`, `SAEDSpotPicker`, `assign_transformation_variant` |

Demonstrated end to end on Burgers β→α in the committed-executed notebook
`docs/site/tutorials/notebooks/23_transformation_crystallography_end_to_end.ipynb`.

### Where the durable records are

- **Specification (normative, now marked delivered):**
  [`docs/architecture/transformation_crystallography_and_diffraction_program.md`](../architecture/transformation_crystallography_and_diffraction_program.md)
- **Phase ledger with the full outcome, the defects found, and the open follow-ons:**
  [`docs/roadmap/working_notes_transformation_diffraction_program.md`](../roadmap/working_notes_transformation_diffraction_program.md)

### Four defects found and fixed during the program

Three were pre-existing and are the reason the ledger is worth reading before starting
anything adjacent:

1. Checksum-pinned fixtures failed on every Windows clone (`core.autocrlf` rewriting
   hash-pinned artifacts). The suite now runs green with no deselects.
2. The seedless OR fit averaged tied double-coset representatives, turning planted Bain into
   a meaningless 26.9 deg that read as Kurdjumov-Sachs.
3. The shared Burgers worked-example setup declared no space groups, so it had been
   simulating β-titanium without body-centring absences and listing forbidden reflections.
4. Kinematic spot ordering was decided by floating-point noise at symmetry-equivalent ties,
   so the same pattern reached two ways came out permuted.

### Open follow-ons from this program

Measured-EBSD fixtures for the OR determination (validation is synthetic and says so); JSON
round-trip contracts for the new report objects; canonical SVG figures and LaTeX theory notes
for the OR-statement extraction and the ratio/angle algorithm; structure-aware cubic-to-cubic
catalog dispatch. None blocks use of the delivered surfaces. Details in the ledger.

## Previous Goal — archived

The reconstruction-stabilization program (Phases 1–5, commits `7dd77d7b` … `c8c6eb1b`) is
complete and its handoff record is archived at
[`docs/development/archive/reconstruction_stabilization_2026_07.md`](archive/reconstruction_stabilization_2026_07.md).
Its three open follow-ons — running the MTEX side of `or_transformation_v1` on a machine with
MATLAB, a measured-data reconstruction fixture, and irregular grain geometry in the map sweep
— are **still open** and remain blockers to moving parent-grain reconstruction out of
`experimental`. They were not part of the TX program.
