# Active Task Progress

This file is the durable handoff record for the current substantial repository task. Keep it
current enough that work can resume after an interrupted agent session without relying on chat
history.

## Current Goal (2026-08-03)

**Transformation Crystallography And Composite Diffraction Program (TX)** — deliver five
user-facing phase-transformation answers on top of the existing OR and diffraction primitives:

- **(a)** measured parent/child Euler angles ⇒ *what is the orientation relationship?*
- **(b)** an OR + an arbitrary parent plane/direction ⇒ the parallel planes and directions in
  every product variant
- **(c)** a parent zone axis ⇒ a robust composite kinematic SAED, exportable as graphics **and**
  reflection tables
- **(d)** a product-variant zone axis ⇒ the same composite, with the matrix and sibling variants
  generated around it
- **(e)** a measured SAED pattern ⇒ solved, from spots picked interactively or listed in YAML

Closed out by a Burgers β↔α notebook demonstrating (a)–(e) end to end.

### Where the live records are

- **Specification (normative):**
  [`docs/architecture/transformation_crystallography_and_diffraction_program.md`](../architecture/transformation_crystallography_and_diffraction_program.md)
- **Running phase ledger (read this to resume):**
  [`docs/roadmap/working_notes_transformation_diffraction_program.md`](../roadmap/working_notes_transformation_diffraction_program.md)

The ledger carries phase status, the verified code baseline, the inherited pre-existing failures
on this machine, the verification command set, and the exact next action.

## Previous Goal — archived

The reconstruction-stabilization program (Phases 1–5, commits `7dd77d7b` … `c8c6eb1b`) is complete
and its handoff record is archived at
[`docs/development/archive/reconstruction_stabilization_2026_07.md`](archive/reconstruction_stabilization_2026_07.md).
Its three open follow-ons — running the MTEX side of `or_transformation_v1` on a machine with
MATLAB, a measured-data reconstruction fixture, and irregular grain geometry in the map sweep — are
**still open** and remain blockers to moving parent-grain reconstruction out of `experimental`.
They are not part of the TX program.
