# Canonical Documentation Map

PyTex uses the Sphinx site as the main browsing surface, but the full documentation program also
includes architecture notes, standards, validation ledgers, and canonical LaTeX research notes
outside `docs/site/`.

This page pulls those sources into one browsable map so the generated HTML site works as a teaching
surface, a library guide, and an entry point into the deeper research documentation.

## Architecture

- `../../architecture/overview.md`
- `../../architecture/canonical_data_model.md`
- `../../architecture/orientation_and_texture_foundation.md`
- `../../architecture/ebsd_foundation.md`
- `../../architecture/diffraction_foundation.md`
- `../../architecture/multimodal_characterization_foundation.md`

## Standards

- `../../standards/engineering_governance.md`
- `../../standards/notation_and_conventions.md`
- `../../standards/hexagonal_and_trigonal_conventions.md`
- `../../standards/documentation_architecture.md`
- `../../standards/development_principles.md`
- `../../standards/data_contracts_and_manifests.md`
- `../../standards/terminology_and_symbol_registry.md`

## Validation And Testing

- `../../testing/strategy.md`
- `../../testing/mtex_parity_matrix.md`
- `../../testing/diffraction_validation_matrix.md`
- `../../testing/plotting_validation_matrix.md`
- `../../testing/structure_validation_matrix.md`

## Canonical LaTeX Notes

The rendered site links these through {doc}`../theory/index`, but the full source corpus lives
under `docs/tex/`:

- `../../tex/README.md`
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

- `../../standards/documentation_architecture.md`
- `../../standards/engineering_governance.md`

### Informative

- `../../README.md`
