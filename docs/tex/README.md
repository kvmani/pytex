# LaTeX Documentation

This directory contains the canonical scientific documentation for PyTex.

## Policy

- Major scientific notes are authored here first.
- The Sphinx-facing docs layer summarizes, indexes, and links to these documents.
- SVG figures under `../figures/` are the canonical editable figure sources.

## Current Documents

- `foundations/project_philosophy.tex`
- `theory/canonical_data_model.tex`
- `theory/crystal_visualization_geometry.tex`
- `theory/hexagonal_conventions.tex`
- `theory/orientation_space_and_disorientation.tex`
- `theory/reference_frames.tex`
- `algorithms/discrete_odf_and_pole_figures.tex`
- `algorithms/powder_xrd_and_saed.tex`
- `validation/validation_program.tex`

## Build Notes

These documents are intended to be built with a modern LaTeX distribution. If `svg` package support is used for direct SVG inclusion, the build will typically require `latexmk -pdf -shell-escape`.
