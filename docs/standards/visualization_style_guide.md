# Visualization Style Guide

This document defines the shared visual language for PyTex documentation illustrations.

It governs canonical SVG assets in `docs/figures/`, especially architecture maps, workflow
flowsheets, process diagrams, frame schematics, and teaching illustrations. Runtime plots may use
the YAML plotting theme system, but promoted documentation figures should follow this guide unless
a scientific convention requires a different visual treatment.

## Design Goals

- Make scientific structure legible at a glance without weakening mathematical precision.
- Use one recognizable PyTex visual language across Sphinx pages, LaTeX-linked figures, notebooks,
  talks, and README summaries.
- Prefer authored SVG for complex or publication-facing diagrams instead of cramped Mermaid blocks.
- Keep visual styling subordinate to scientific meaning: colors group concepts, arrows show
  dependency or data flow, and labels use repository terminology.

## Canonical Tokens

| Token | Value | Use |
| --- | --- | --- |
| ink | `#07122f` | titles, primary labels, arrow strokes |
| muted text | `#40506f` | subtitles, captions, secondary labels |
| paper | `#fbfdff` | figure background |
| core blue | `#2563eb` | canonical core, primary data path |
| teal | `#0f9f9f` | validation, measurement, acquisition, accepted flow |
| violet | `#7c3aed` | documentation, standards, teaching, semantic governance |
| amber | `#f59e0b` | adapters, boundary layers, caution or decision nodes |
| rose | `#e11d48` | errors, rejected branches, risk, incomplete status |
| green | `#16a34a` | implemented, verified, successful outputs |

SVG assets should use soft gradient fills derived from these tokens, restrained shadows, and
consistent stroke colors. Avoid one-note palettes that make every subsystem look equivalent.

## Layout Rules

- Use authored SVG for diagrams that explain architecture, process flow, governance, validation
  posture, or multi-step scientific workflows.
- Avoid long single-row or single-column flowsheets when the content naturally has phases. Prefer
  two-row, lane-based, radial, or grouped layouts that use the available screen area.
- Keep repeated node cards near 8px corner radius.
- Use `Arial, Helvetica, sans-serif` for all canonical SVG text. Do not use serif fallback fonts for
  titles, labels, symbols, or callouts unless a figure explicitly embeds a publication source
  facsimile.
- Use 18px or larger node labels and at least 15px body text inside SVG assets.
- Use 44px icon medallions or equivalent visual anchors for major nodes when icons clarify the
  diagram.
- Keep arrows orthogonal or gently curved; avoid dense crossing lines.
- Put legends inside the SVG only when color, stroke style, or line style carries meaning.
- Every canonical SVG must include a `<title>` and `<desc>` element.

### Class And Object Model Diagrams

Class-model diagrams are the one registered exception to the 18px-node-label rule above. A UML
card is a table of declared field names, not a poster node: at poster type sizes a fifteen-class
domain view would not fit any screen, and truncating the field lists would remove the content the
diagram exists to show. They use their own fixed scale:

| Element | Size |
| --- | --- |
| class name | 18px bold |
| attribute lines | 15px |
| module path, stereotype, relation labels | 13px |

Everything else in this guide applies unchanged: canonical tokens, Arial, 8px card radius,
mandatory `<title>` and `<desc>`, and arrowheads in `userSpaceOnUse`. These figures are generated
by `pytex.plotting.class_diagrams` from the introspected model in `scripts/class_model.py`; do not
hand-edit them. Because they hold their drawn size rather than scaling to the content column, the
Sphinx pages that embed them must use the `class-atlas-figure` scroll container.

## Diagram Categories

- `architecture maps`: layered or lane-based diagrams that show modules, governance, boundaries,
  and future foundations.
- `process flows`: staged diagrams that show input, validation, computation, failure branches, and
  outputs.
- `scientific geometry`: diagrams that fix frames, bases, vectors, planes, angles, detector axes,
  or projection conventions.
- `validation and evidence`: diagrams that show how tests, fixtures, manifests, benchmarks, and
  claims relate.
- `teaching summaries`: simplified diagrams allowed only when labels make the simplification clear
  and do not contradict the canonical model.

## Mermaid Policy

Mermaid remains acceptable for small local diagrams, drafts, or pages where the diagram is not a
canonical visual reference. Complex diagrams that are reused across the public documentation,
architecture atlas, teaching pages, or governance material should be promoted to SVG in
`docs/figures/`.

When replacing Mermaid with SVG:

- preserve the scientific relationships from the original diagram
- improve layout density instead of copying the same cramped row or column
- use repository terminology from the terminology and symbol registry
- add meaningful alt text where the figure is referenced
- add or update integrity checks when the figure becomes part of the canonical documentation set

## References

### Normative

- [Documentation Architecture](documentation_architecture.md)
- [LaTeX And Figures](latex_and_figures.md)
- [Terminology And Symbol Registry](terminology_and_symbol_registry.md)

### Informative

- <a href="../validation/plotting_validation_matrix.html">Plotting Validation Matrix</a>
