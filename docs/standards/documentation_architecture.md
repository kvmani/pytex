# Documentation Architecture

PyTex uses a hybrid documentation architecture so the project can be both scientifically rigorous and easy to use.

## Core Policy

- Sphinx is the primary browsable and searchable documentation surface.
- LaTeX is the canonical source for major scientific notes.
- SVG is the canonical source for scientific figures and geometry schematics.

The public documentation entry point is therefore the HTML docs layer, not the PDF layer alone.

## Content Types

Use the Sphinx layer for:

- concept pages
- tutorials
- workflow guides
- installation and contributor guidance
- curated API reference
- links to deeper scientific notes

Use the LaTeX layer for:

- formal theory notes
- algorithm notes
- convention-fixing documents
- validation and limitations notes

## Cross-Link Rule

- Every major concept or workflow page should link to the deeper theory or algorithm note where one exists.
- Major scientific notes should point back to the relevant user-facing concept or workflow page when practical.
- API docs should not stand alone; they should connect back to concepts or workflows where scientific meaning matters.

## Explicitness Requirement

For major conventions, frame relationships, symmetry actions, and core algorithms, the documentation must expose all three of the following surfaces together:

- prose that states the scientific meaning in plain language
- mathematics that fixes the convention, mapping, or reduction unambiguously
- annotated figures that make the geometry or transformation legible at a glance

If one of these surfaces is missing for a stable foundational feature, the documentation is incomplete even if the code and tests exist.

## Deliverables

The intended docs surfaces are:

- HTML documentation built from Sphinx/MyST
- downloadable PDF scientific notes built from LaTeX
- reusable SVG figures shared across both surfaces
- mathematics-and-graphics-backed convention pages for foundational semantics

## Planned Site Shape

The Sphinx layer should eventually expose top-level sections for:

- concepts
- tutorials
- workflows
- API reference
- theory and standards index
- validation and benchmarking index

## Tooling Expectations

- Sphinx and MyST are the default tools for the HTML layer.
- Math-heavy web pages should use MathJax-backed equations.
- Citations in the HTML layer should use a BibTeX-backed workflow.
- The LaTeX layer should remain suitable for publication-grade mathematical exposition.

## Cross-Platform Requirement

- Ordinary users should be able to browse and search the built HTML docs uniformly across Linux, Windows, and macOS.
- The documentation toolchain should remain practical across Linux, Windows, and macOS for contributors.
- Full PDF builds may require a documented TeX toolchain, but the public documentation system must not depend on users having LaTeX installed locally.

## References

### Normative

- `reference_canon.md`
- `scientific_citation_policy.md`

### Informative

- `../site/README.md`
