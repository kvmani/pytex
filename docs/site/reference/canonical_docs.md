# Canonical Documentation Map

PyTex uses the Sphinx site as the main browsing surface, but the full documentation program also
includes architecture notes, standards, validation ledgers, and canonical LaTeX research notes
outside `docs/site/`.

This page pulls those sources into one browsable map so the generated HTML site works as a teaching
surface, a library guide, and an entry point into the deeper research documentation.

## Architecture

- {doc}`../architecture/overview`
- {doc}`../architecture/canonical_data_model`
- {doc}`../architecture/orientation_and_texture_foundation`
- {doc}`../architecture/ebsd_foundation`
- {doc}`../architecture/diffraction_foundation`
- {doc}`../architecture/multimodal_characterization_foundation`

## Standards

- {doc}`../standards/engineering_governance`
- {doc}`../standards/notation_and_conventions`
- {doc}`../standards/hexagonal_and_trigonal_conventions`
- {doc}`../standards/documentation_architecture`
- {doc}`../standards/development_principles`
- {doc}`../standards/data_contracts_and_manifests`
- {doc}`../standards/terminology_and_symbol_registry`

## Validation And Testing

- {doc}`../validation/testing_strategy`
- {doc}`../validation/mtex_parity_matrix`
- {doc}`../validation/diffraction_validation_matrix`
- {doc}`../validation/plotting_validation_matrix`
- {doc}`../validation/structure_validation_matrix`

## Canonical LaTeX Notes

The rendered site links these through {doc}`../theory/index`, but the full source corpus lives
under `docs/tex/`:

- `../../tex/README.md`
- <a href="../../tex/README.md">LaTeX README</a>
- `../../tex/theory/`
- `../../tex/algorithms/`
- `../../tex/validation/`

## How To Use The Site

- Start with `Concepts` if you need scientific semantics and vocabulary.
- Use `Tutorials` and `Workflows` for executable, end-to-end teaching paths.
- Use `API` for the curated stable surface.
- Use `Theory` for linked canonical LaTeX notes and algorithm references.
- Use this page when you need the complete repository documentation map, including architecture and
  standards that anchor the public API.

## References

### Normative

- {doc}`../standards/documentation_architecture`
- {doc}`../standards/engineering_governance`

### Informative

- <a href="../../README.md">Repository README</a>
