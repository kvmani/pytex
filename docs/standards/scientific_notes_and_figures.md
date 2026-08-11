# Scientific Notes And Figures Standard

MyST Markdown is the canonical source for major scientific notes, while SVG remains the canonical source for scientific figures.

This standard works together with `documentation_architecture.md`, where Sphinx is defined as the primary browsable and searchable documentation surface.
Canonical SVG assets that explain architecture, process flow, governance, validation, workflow, or
teaching material must also follow `visualization_style_guide.md`.

## Canonical Sources

- `docs/site/theory/` contains the canonical scientific notes, authored as MyST Markdown.
- `docs/figures/` contains the canonical SVG figure sources.

## Why The Notes Are Not LaTeX

Until 2026-08-11 these notes were authored as LaTeX under `docs/tex/` and the Sphinx site linked
them as raw `.tex` downloads. Sphinx has no LaTeX-parsing extension configured, so the derivations
never rendered: a reader following a link got a source file instead of the mathematics. Keeping
LaTeX canonical *and* rendering it on the site would have required a conversion pipeline and left
two representations free to drift.

One source removes both problems. The notes render natively because `myst_enable_extensions`
carries `amsmath` and `dollarmath`, and a typeset PDF comes from Sphinx's own builder:

```bash
python -m sphinx -b latexpdf docs/site docs/_build/latex
```

Nothing was lost in the move. The corpus used no `\cite`, `\ref`, `\includegraphics`,
`\newcommand`, `\input`, or TikZ, and 26 of its 37 files had no `\documentclass` at all — they were
fragments that had never been standalone-compilable, so the documented `latexmk` build could not
have run for most of them.

## Authoring Rules For Notes

- Write notes as MyST Markdown in `docs/site/theory/`, one H1 per page, and add each new note to a
  `toctree` group in `docs/site/theory/index.md`.
- Use `$…$` for inline mathematics and `$$…$$` for display mathematics. Number an equation that is
  referred to later by appending a label — ``$$ … $$ (eq-name)`` — and cite it with
  ``{eq}`eq-name` ``.
- `align`, `gather`, `cases`, and the other amsmath environments may be written directly, without
  `$$` fences; the `amsmath` extension renders them as block mathematics.
- Cross-reference other notes and pages with the `{doc}` role rather than a file path, so a moved
  page fails the build instead of rotting into a dead link.
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
- follow the central visualization style guide when the figure is an architecture, process-flow,
  workflow, validation, or teaching illustration

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

A major scientifically substantial stable feature is incomplete until its theory note, SVG figures, tests, and validation note all exist.

## References

### Normative

- `documentation_architecture.md`
- `reference_canon.md`
- `visualization_style_guide.md`

### Informative

- `../site/theory/index.md`
