# Crystal Structures And CIF Import

PyTex now supports phase construction from crystallographic structures and CIF data at the canonical data-model layer rather than only through external adapters.

## Core Rule

If a crystallographic source defines a phase, PyTex should construct canonical primitives directly:

- a lattice,
- a unit cell,
- a phase,
- a point-group symmetry specification derived from the source structure.

This ensures that downstream orientation, EBSD, and diffraction workflows operate on PyTex primitives rather than on source-specific structure objects.

## Current Implementation

The present implementation uses `pymatgen` as an optional crystallographic parser and structure provider. PyTex exposes:

- `Lattice.from_pymatgen_lattice(...)`,
- `UnitCell.from_pymatgen_structure(...)`,
- `Phase.from_pymatgen_structure(...)`,
- `Phase.from_cif(...)`,
- `Phase.from_cif_string(...)`.

The Workbench's shared phase selector uses this same `Phase.from_cif_string(...)` boundary for
user-loaded `.cif` files. It converts the resulting canonical `Phase` to the application's
JSON-ready phase specification, so Crystal Viewer, XRD, TEM, EBSD forward simulation, texture and
orientation-relationship tools do not acquire parser-specific semantics.

The imported structure is normalized into PyTex objects, including:

- lattice parameters,
- atomic basis,
- reduced chemical formula,
- space-group symbol and number,
- crystal point group reduced to the supported `SymmetrySpec` surface.

## Disordered Sites

PyTex stores one species and one occupancy per `AtomicSite`. When a source structure contains a disordered crystallographic site with multiple species, the current import path expands that source site into multiple PyTex `AtomicSite` records sharing the same fractional coordinates and carrying explicit occupancy values.

## Current Limits

- CIF support currently depends on `pymatgen` rather than on a native PyTex parser.
- Space-group information is retained on the phase object, but symmetry algorithms currently operate at the point-group level implemented in `SymmetrySpec`.
- Magnetic symmetry, superspace descriptions, and richer crystallographic metadata remain future work.

## Normative And Informative References

- International Union of Crystallography, *Crystallographic Information Framework (CIF)*. <https://www.iucr.org/resources/cif>.
- S. P. Ong et al., *Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis*, Computational Materials Science 68 (2013) 314–319. DOI: <https://doi.org/10.1016/j.commatsci.2012.10.028>.
