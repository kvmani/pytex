# Core Model: Phases, Unit Cells, And CIF Import

PyTex now supports creating canonical phase primitives directly from crystallographic structures and CIF files.

This support lives at the core-model layer, not as an afterthought in a downstream importer, because phase construction fixes semantics that the rest of the library depends on:

- lattice parameters
- crystal frame ownership
- point-group symmetry used for orientation reduction
- unit-cell atomic basis
- provenance of the imported structure

## Scope

- `Lattice.from_pymatgen_lattice(...)`
- `UnitCell.from_pymatgen_structure(...)`
- `Phase.from_pymatgen_structure(...)`
- `Phase.from_cif(...)`
- `Phase.from_cif_string(...)`

These constructors use `pymatgen` as the optional crystallographic parser and structure source, but the returned objects are pure PyTex primitives.

## Why This Lives In The Core Model

PyTex does not want CIF import to become a sidecar convenience that bypasses the canonical data model. If a CIF file defines a phase, that phase should enter the library as:

- a `Lattice`
- a `UnitCell`
- a `Phase`
- a `SymmetrySpec`

with explicit frame ownership and stable invariants from the start.

## Installation

The CIF constructors require the optional adapters stack:

```bash
python -m pip install -e '.[adapters]'
```

If `pymatgen` is not installed, the constructors fail with an explicit error explaining that the adapters extra is required.

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
- The constructors create canonical PyTex objects; they do not expose `pymatgen` structures as part of the stable PyTex public API.
- CIF-backed construction is phase-centric and unit-cell-centric. It is not yet a full crystallographic database or structure-editing subsystem.

## Current Limits

- CIF import currently depends on `pymatgen` rather than a native PyTex parser
- space-group information is retained on `Phase`, but PyTex symmetry algorithms still operate on the proper point-group layer implemented in `SymmetrySpec`
- magnetic symmetry, modulated structures, and richer crystallographic metadata are still ahead

## Related Material

- {doc}`../concepts/core_model`
- {doc}`../concepts/how_pytex_differs`
- [../../tex/theory/canonical_data_model.tex](../../tex/theory/canonical_data_model.tex)
- [../../tex/theory/crystal_structures_and_cif_import.tex](../../tex/theory/crystal_structures_and_cif_import.tex)
