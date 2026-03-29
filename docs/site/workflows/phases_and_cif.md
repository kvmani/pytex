# Core Model: Phases, Unit Cells, And CIF Import

PyTex now supports creating canonical phase primitives directly from crystallographic structures and CIF files.

This support lives at the core-model layer, not as an afterthought in a downstream importer, because phase construction fixes semantics that the rest of the library depends on:

- lattice parameters
- crystal frame ownership
- point-group symmetry used for orientation reduction
- space-group identity used for structure definition
- unit-cell atomic basis
- provenance of the imported structure

## Scope

- `Lattice.from_pymatgen_lattice(...)`
- `UnitCell.from_pymatgen_structure(...)`
- `Phase.from_pymatgen_structure(...)`
- `Phase.from_cif(...)`
- `Phase.from_cif_string(...)`
- `phase_fixture_catalog_path()`
- `list_phase_fixtures()`
- `get_phase_fixture(...)`

These constructors use `pymatgen` as the optional crystallographic parser and structure source, but the returned objects are pure PyTex primitives.

## Built-In Phase Fixtures

PyTex also ships a small pinned in-repo phase-fixture corpus under `fixtures/phases/README.md`.
These fixtures exist for:

- structure-import validation
- benchmark reproducibility
- small canonical demos
- teaching workflows that need stable crystallographic examples

The public fixture helpers make that corpus discoverable without hard-coding file paths:

```python
from pytex import (
    FrameDomain,
    Handedness,
    ReferenceFrame,
    get_phase_fixture,
    list_phase_fixtures,
)

print([record.fixture_id for record in list_phase_fixtures()])

crystal = ReferenceFrame(
    "crystal",
    FrameDomain.CRYSTAL,
    ("a", "b", "c"),
    Handedness.RIGHT,
)

zr = get_phase_fixture("zr_hcp").load_phase(crystal_frame=crystal)
print(zr.name, zr.space_group_symbol, zr.space_group_number)
```

The fixture catalog is hash-pinned, so the built-in corpus participates in the reproducibility and integrity surface rather than behaving like an informal sample-data directory.

## Fixture-Corpus Audit Surface

The built-in fixture corpus is also consumed by the structure-import audit workflow so every pinned
fixture participates in the same reproducibility surface:

- `fe_bcc` for BCC cubic structure-import and diffraction smoke cases
- `zr_hcp` for hexagonal-axis conventions, structure import, and crystal-visualization teaching
- `ni_fcc` for FCC structure import and diffraction external-baseline coverage
- `nicl` for multi-species unit-cell normalization and manifest-backed validation
- `diamond` for covalent cubic structure import and teaching-oriented crystal examples

The pinned audit summary for this corpus lives in
`benchmarks/structure_import/phase_fixture_audit_summary.json`.

## Why This Lives In The Core Model

PyTex does not want CIF import to become a sidecar convenience that bypasses the canonical data model. If a CIF file defines a phase, that phase should enter the library as:

- a `Lattice`
- a `UnitCell`
- a `Phase`
- a `SymmetrySpec`
- a `SpaceGroupSpec`

with explicit frame ownership and stable invariants from the start.

## Installation

The standard contributor and documentation installs now include the CIF-backed structure-import
stack used by this workflow:

```bash
python -m pip install -e '.[dev]'
```

For the docs build path:

```bash
python -m pip install -e '.[dev,docs]'
```

The broader `adapters` extra remains optional for heavier interoperability work outside the normal
structure-import baseline.

## Minimal Example: From CIF Text

```python
from pytex import FrameDomain, Handedness, Phase, ReferenceFrame

crystal = ReferenceFrame(
    "crystal",
    FrameDomain.CRYSTAL,
    ("a", "b", "c"),
    Handedness.RIGHT,
)

phase = Phase.from_cif_string(cif_text, crystal_frame=crystal)

print(phase.name)
print(phase.space_group_symbol, phase.space_group_number)
print(phase.space_group.symbol, phase.space_group.number)
print(phase.symmetry.point_group)
print(phase.chemical_formula)
print(len(phase.unit_cell.sites))
```

## Example: From A CIF File On Disk

```python
from pytex import FrameDomain, Handedness, Phase, ReferenceFrame

crystal = ReferenceFrame(
    "crystal",
    FrameDomain.CRYSTAL,
    ("a", "b", "c"),
    Handedness.RIGHT,
)

phase = Phase.from_cif(
    "NaCl.cif",
    crystal_frame=crystal,
    phase_name="rocksalt",
)
```

## Example: From An Existing `pymatgen` Structure

```python
from pymatgen.core import Structure

from pytex import FrameDomain, Handedness, Phase, ReferenceFrame

crystal = ReferenceFrame(
    "crystal",
    FrameDomain.CRYSTAL,
    ("a", "b", "c"),
    Handedness.RIGHT,
)

structure = Structure.from_file("quartz.cif")
phase = Phase.from_pymatgen_structure(structure, crystal_frame=crystal)
```

## What PyTex Derives Automatically

When a phase is created from a CIF-backed structure, PyTex derives:

- lattice parameters from the crystallographic lattice
- point-group symmetry from the structure’s space-group analysis
- `UnitCell` atomic sites from the crystallographic basis
- reduced chemical formula
- space-group symbol and number
- a first-class `SpaceGroupSpec`
- provenance indicating CIF-backed creation

The critical point is that PyTex uses the structure source only to populate stable canonical objects. Downstream texture, EBSD, and diffraction code then works on PyTex types instead of on `pymatgen` objects.

## Conventional Versus Primitive Cells

The CIF constructors accept `primitive=True` when you want PyTex to reduce the imported structure to a primitive cell before building the canonical objects.

- default behavior keeps the parsed structure as provided by the crystallographic reader
- `primitive=True` requests primitive-cell reduction before `Lattice`, `UnitCell`, and `Phase` are created

## Disordered Sites And Occupancy

PyTex stores `AtomicSite` as one species plus one occupancy per site record. If the source structure contains a disordered site with multiple species, the CIF-backed construction expands that source site into multiple PyTex `AtomicSite` entries with the same fractional coordinates and separate occupancy values.

This keeps the public `AtomicSite` type explicit and avoids smuggling site disorder through ambiguous dictionaries.

## Interpretation Notes

- PyTex currently derives point-group symmetry from space-group analysis, then normalizes that to the supported point-group surface already implemented in `SymmetrySpec`.
- `SpaceGroupSpec` stores the structure-facing space-group identity, while `SymmetrySpec` stores the point-group-facing reduction surface used in orientation workflows.
- The constructors create canonical PyTex objects; they do not expose `pymatgen` structures as part of the stable PyTex public API.
- CIF-backed construction is phase-centric and unit-cell-centric. It is not yet a full crystallographic database or structure-editing subsystem.

## Current Limits

- CIF import currently depends on `pymatgen` rather than a native PyTex parser
- space-group information is retained on `Phase`, but PyTex symmetry algorithms still operate on the proper point-group layer implemented in `SymmetrySpec`
- magnetic symmetry, modulated structures, and richer crystallographic metadata are still ahead
- broader literature-backed structure-import benchmarking is still ahead of the current invariant-driven baseline

## Related Material

- {doc}`../concepts/core_model`
- {doc}`../concepts/how_pytex_differs`
- {doc}`../tutorials/notebooks/04_phases_lattices_space_groups_and_cif`
- [../../tex/theory/canonical_data_model.tex](../../tex/theory/canonical_data_model.tex)
- [../../tex/theory/crystal_structures_and_cif_import.tex](../../tex/theory/crystal_structures_and_cif_import.tex)
