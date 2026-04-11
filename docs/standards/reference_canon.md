# Reference Canon

This document defines the allowed source corpus for fixing PyTex scientific semantics.

PyTex is not allowed to adopt conventions opportunistically from whatever library or notebook already implements a feature. Stable scientific behavior must be anchored to an explicit reference class.

## Core Rule

Every scientific rule in PyTex must be classified as one of:

- `normative from IUCr/International Tables`
- `normative by PyTex adoption of a named literature or tool convention`
- `informative only`

If a rule is informative only, it must not silently fix stable API semantics.

## Repository Working Corpus

PyTex maintains a working reference corpus under `../../references/` to make authoritative formulas and examples discoverable during implementation.

- [Reference Index](../../references/reference_index.md) is the page-level discovery index.
- [Formulation Summary](../../references/formulation_summary.md) is the first-stop summary for formulas, notation crosswalks, and ready-to-use examples.
- [Feature Opportunities](../../references/feature_opportunities.md) records source-grounded implementation ideas and gap analysis.

These files are repository guidance tools, not a replacement for the source hierarchy below. When they summarize a source, contributors should still verify the cited pages before fixing stable semantics.

## Normative Corpus By Topic

### Symmetry And Crystallographic Structure

- IUCr / *International Tables for Crystallography, Volume A*
- IUCr / *International Tables for Crystallography, Volume G*

These sources define the normative backbone for point-group, space-group, and CIF/data-exchange semantics.

### Orientation And Texture Conventions

- Bunge, *Texture Analysis in Materials Science*
- IUCr and International Tables material where it fixes crystallographic convention directly

PyTex may additionally adopt named MTEX-facing or literature-facing conventions where the docs explicitly say so.

### Frame Discipline Across Texture And Neutron Workflows

- Nolze et al. (2023), *Journal of Applied Crystallography*, DOI: <https://doi.org/10.1107/S1600576723009275>

This paper is especially important when PyTex fixes coordinate-system doctrine for texture and diffraction workflows spanning specimen and laboratory reasoning.

### EBSD Geometry And Interpretation

- Schwartz, Kumar, Adams, Field, *Electron Backscatter Diffraction in Materials Science*
- relevant maintained tool documentation only when PyTex intentionally adopts a tool-specific interpretation

### Diffraction Geometry And TEM Reasoning

- IUCr and *International Tables for Crystallography, Volume C*
- De Graef, *Introduction to Conventional Transmission Electron Microscopy*
- Fultz and Howe, *Transmission Electron Microscopy and Diffractometry of Materials*

### Phase Transformation

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*
- peer-reviewed orientation-relationship and variant-generation literature named explicitly in the relevant theory notes

## Source Hierarchy

Use this hierarchy when a new rule must be fixed:

1. IUCr and International Tables
2. formal standards and dictionaries such as CIF-related standards
3. canonical textbooks used in the field
4. peer-reviewed papers
5. maintained project documentation for tools such as MTEX or ORIX
6. vendor or educational notes

## Adoption Rule For Tool Conventions

PyTex may adopt a tool convention from MTEX, ORIX, KikuchiPy, diffsims, or vendor software only when:

- the adoption is explicit
- the scientific meaning is described in PyTex’s own docs
- the convention does not conflict with a higher-ranked normative source
- the validation note explains why PyTex follows that convention

## Disallowed Shortcuts

- Do not define stable semantics purely from existing code behavior.
- Do not treat a parser library’s default behavior as normative without checking the scientific source it is meant to implement.
- Do not make a maintained-tool behavior normative unless the docs state that PyTex intentionally adopts it.

## References

### Normative

- Hahn, Th. (ed.), *International Tables for Crystallography, Volume A: Space-Group Symmetry*
- Hall, S. R. and McMahon, B. (eds.), *International Tables for Crystallography, Volume G: Definition and Exchange of Crystallographic Data*
- Bunge, H.-J., *Texture Analysis in Materials Science: Mathematical Methods*
- Nolze et al. (2023), DOI: <https://doi.org/10.1107/S1600576723009275>

### Informative

- `scientific_citation_policy.md`
- [Reference Index](../../references/reference_index.md)
- [Formulation Summary](../../references/formulation_summary.md)
- MTEX documentation: <https://mtex-toolbox.github.io/>
