# LaTeX And Figures Standard

LaTeX remains the canonical source for major scientific notes, while SVG remains the canonical source for scientific figures.

This standard works together with `documentation_architecture.md`, where Sphinx is defined as the primary browsable and searchable documentation surface.

## Canonical Sources

- `docs/tex/` contains the canonical scientific notes.
- `docs/figures/` contains the canonical SVG figure sources.
- Runtime plotting APIs may return ordinary Matplotlib figures for user code; the SVG rule applies to canonical repository-tracked documentation figures.
- Runtime plotting regression tests must not introduce repo-tracked SVG byte baselines. Publication-facing runtime surfaces should be validated through structural and semantic assertions unless a figure is intentionally promoted into `docs/figures/` as a canonical documentation asset.

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

This SVG rule is for repository assets and documentation figures. It does not require every user-generated runtime plot to be exported as SVG.

When a plotting surface is also part of the runtime API, its visual defaults should be controlled through the shared YAML style system rather than duplicated inside workflow examples or documentation scripts.

## Required Figure Topics

At minimum, the repository must maintain canonical figures for:

- shared reference-frame domains and transforms
- Euler and quaternion storage conventions
- diffraction geometry
- powder XRD and SAED geometry where those figures fix conventions or detector mappings
- crystal-structure and plane-geometry schematics for the 3D visualization subsystem
- hexagonal and trigonal indexing conventions, including HCP 3-index and 4-index notation

## Completion Rule

A major scientifically substantial stable feature is incomplete until its LaTeX note, SVG figures, tests, and validation note all exist.

## References

### Normative

- `documentation_architecture.md`
- `reference_canon.md`

### Informative

- `../tex/README.md`
