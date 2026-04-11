# Sphinx Site

This directory contains the primary Sphinx documentation surface for PyTex and is the public documentation entry point.

## Build

```bash
python -m pip install -e '.[docs]'
sphinx-build -b html docs/site docs/_build/html
```

For full platform-specific install, notebook, and PDF-build guidance, see
{doc}`tutorials/installation_and_build`.

## Structure

- `index.md`: public documentation entry point
- `concepts/`: conceptual documentation for the canonical model
- `tutorials/`: runnable onboarding material
- `tutorials/installation_and_build.md`: full install and build guidance for Windows, macOS, and Linux
- `tutorials/notebooks/`: executable deep-dive tutorials in `.ipynb` form
- `workflows/`: multi-step domain workflows
- `api/`: curated API navigation
- `theory/`: pointers into the canonical LaTeX notes
- `validation/`: testing strategy and parity posture
- [Implementation roadmap](../roadmap/implementation_roadmap.md): active development guide for the next phase

## Policy

- Sphinx is the public documentation entry point.
- Concept pages, tutorials, workflows, and curated API reference live here.
- Canonical scientific notes remain under `../tex/`.
- Canonical SVG figures remain under `../figures/`.

## Architecture

The architecture atlas is the best starting point for understanding the library as a whole:

- {doc}`concepts/library_structure`
- `../figures/pytex_architecture_poster.svg`
- `../figures/pytex_architecture_evolution_poster.svg`
- `../figures/pytex_architecture_compact.svg`
