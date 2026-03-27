# PyTex

PyTex is a pure-Python-first crystallographic texture and diffraction library built around an explicit canonical data model for frames, symmetry, orientations, EBSD maps, and diffraction geometry.

PyTex is intentionally opinionated about scientific semantics. Stable APIs are expected to make frame, basis, symmetry, and provenance meaning explicit instead of hiding that meaning in arrays or import-time conventions.

```{toctree}
:maxdepth: 2
:caption: Concepts

concepts/core_model
concepts/core_foundation
concepts/library_structure
concepts/how_pytex_differs
concepts/reference_frames_and_conventions
concepts/orientation_texture
concepts/texture_foundation
concepts/ebsd_foundation
concepts/diffraction_foundation
concepts/symmetry_and_fundamental_regions
```

```{toctree}
:maxdepth: 2
:caption: Tutorials

tutorials/quickstart
tutorials/notebooks
```

```{toctree}
:maxdepth: 2
:caption: Workflows

workflows/ebsd_kam
workflows/ebsd_grains
workflows/ebsd_import_normalization
workflows/phases_and_cif
workflows/texture_odf_inversion
workflows/plotting_primitives
workflows/diffraction_geometry
workflows/diffraction_spots
workflows/ipf_colors
```

```{toctree}
:maxdepth: 2
:caption: Reference

api/index
theory/index
validation/index
```

## Why PyTex Exists

Many texture and diffraction workflows are scientifically fragile at tool boundaries. A rotation matrix, a quaternion, or an EBSD map is not self-describing unless the associated frame, symmetry, and convention choices remain attached to it. PyTex treats those choices as part of the stable product surface.

This leads to three concrete design commitments:

- semantic types before convenience arrays
- one shared frame and symmetry model across texture, EBSD, and diffraction
- documentation and validation artifacts treated as part of feature completion

## How PyTex Differs

If you are evaluating PyTex against MTEX, orix, LaboTex, or vendor EBSD environments, read {doc}`concepts/how_pytex_differs`. That page is the dedicated comparison surface and explains the specific design choices that make PyTex a different kind of project rather than just another feature checklist.

```{note}
PyTex is designed to serve both research and teaching. The Sphinx pages explain concepts and workflows; the LaTeX notes under `docs/tex/` remain the canonical deep scientific source for theory, algorithms, and validation posture.
```

## Current Scope

- Stable core primitives for frames, symmetry, lattice semantics, provenance, batch semantics, rotations, orientations, and transformation records.
- Foundational texture-domain support for pole figures, inverse pole figures, IPF color keys, and discrete kernel ODF evaluation.
- Runtime semantic plotting support for vectors, symmetry elements and orbits, rotations, orientations, PF/IPF objects, and discrete ODF views.
- Real EBSD workflow support for regular-grid neighbor relationships, KAM, grain segmentation, GROD, boundary extraction, cleanup, and import-manifest normalization contracts.
- Core-model phase creation from crystallographic structures and CIF files, with explicit point-group and space-group semantics.
- Shared multimodal acquisition primitives for geometry, calibration, quality, scattering, and experiment-manifest emission.
- Stable manifest families for import, experiment, benchmark, validation, and workflow-result interchange.
- Diffraction-domain containers for geometry and pattern data, with foundational experiment-context integration.
- Executable notebook tutorials that mirror the implemented features and their mathematical contracts.

## Recommended Reading Path

### If You Are New To Texture Or EBSD

1. Start with {doc}`tutorials/quickstart`.
2. Read {doc}`concepts/core_model` to understand why frames and symmetry are explicit.
3. Continue to {doc}`concepts/orientation_texture` for Euler, quaternion, and symmetry-reduction semantics.
4. Then move to {doc}`workflows/ebsd_kam` and {doc}`workflows/ebsd_grains`.

### If You Are Evaluating Scientific Rigor

1. Read {doc}`concepts/core_model`.
2. Check {doc}`validation/index`.
3. Follow the linked architecture notes and LaTeX theory or algorithm notes.
4. Review the MTEX parity ledger in `docs/testing/mtex_parity_matrix.md`.

## Canonical Assets

- Architecture overview: `docs/architecture/overview.md`
- Canonical data model: `docs/architecture/canonical_data_model.md`
- Orientation and texture foundation: `docs/architecture/orientation_and_texture_foundation.md`
- EBSD foundation: `docs/architecture/ebsd_foundation.md`
- Diffraction foundation: `docs/architecture/diffraction_foundation.md`
- Multimodal characterization foundation: `docs/architecture/multimodal_characterization_foundation.md`
- Validation program: `docs/testing/strategy.md`

![Reference Frames](../figures/reference_frames_vectors.svg)

## Scientific Posture

PyTex does not treat external-tool agreement as a vague aspiration. For currently covered areas, MTEX is the validation floor, not the ceiling. Fixture-backed parity tests, unit tests, and theory notes are all part of the definition of a stable feature.

See {doc}`validation/index` for the current validation surface and explicit current limits.

## Normative And Informative References

### Normative For PyTex Conventions

- `docs/standards/notation_and_conventions.md`
- `docs/standards/hexagonal_and_trigonal_conventions.md`
- `docs/architecture/canonical_data_model.md`

### Informative External References

- Hielscher, Schaeben and collaborators, [MTEX documentation](https://mtex-toolbox.github.io/).
- Nolze et al., [Texture measurements on quartz single crystals to validate coordinate systems for neutron time-of-flight texture analysis](https://journals.iucr.org/j/issues/2023/06/00/xx5024/), *Journal of Applied Crystallography*.
- Bunge, *Texture Analysis in Materials Science* (1982).
