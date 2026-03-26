# LaTeX And Figures Standard

LaTeX remains the canonical source for major scientific notes, while SVG remains the canonical source for scientific figures.

This standard works together with `documentation_architecture.md`, where Sphinx is defined as the primary browsable and searchable documentation surface.

## Canonical Sources

- `docs/tex/` contains the canonical scientific notes.
- `docs/figures/` contains the canonical SVG figure sources.

## Required Scientific Note Set Per Major Stable Feature

- mathematical or theory note
- algorithm and implementation note
- validation and limitations note

## Figure Requirements

Where geometry or conventions matter, figures must:

- be maintained canonically as SVG
- identify frames and axes explicitly
- label vectors, planes, poles, angles, and units
- be suitable for papers, lectures, and documentation reuse
- cite the scientific convention or source they are illustrating when the figure fixes a standard

## Required Figure Topics

At minimum, the repository must maintain canonical figures for:

- shared reference-frame domains and transforms
- Euler and quaternion storage conventions
- diffraction geometry
- hexagonal and trigonal indexing conventions, including HCP 3-index and 4-index notation

## Completion Rule

A major scientifically substantial stable feature is incomplete until its LaTeX note, SVG figures, tests, and validation note all exist.

## References

### Normative

- `documentation_architecture.md`
- `reference_canon.md`

### Informative

- `../tex/README.md`
