# PyTex

PyTex is a pure-Python-first crystallographic texture and diffraction library built around an
explicit canonical data model for frames, symmetry, orientations, EBSD maps, diffraction geometry,
and research-grade provenance.

PyTex is intentionally opinionated about scientific semantics. Stable APIs are expected to make
frame, basis, symmetry, and provenance meaning explicit instead of hiding that meaning inside
arrays or import-time conventions.

```{toctree}
:maxdepth: 2
:caption: Concepts

concepts/core_model
concepts/core_foundation
concepts/library_structure
concepts/how_pytex_differs
concepts/reference_frames_and_conventions
concepts/miller_planes_directions
concepts/technical_glossary_and_symbols
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
tutorials/installation_and_build
tutorials/notebooks
```

```{toctree}
:maxdepth: 2
:caption: Workflows

workflows/ebsd_import_normalization
workflows/ebsd_kam
workflows/ebsd_grains
workflows/ebsd_to_texture_outputs
workflows/phase_transformation_manifests_and_scoring
workflows/vectorized_miller_workflows
workflows/orix_kikuchipy_interop
workflows/phases_and_cif
workflows/texture_odf_inversion
workflows/plotting_primitives
workflows/stereographic_projections
workflows/diffraction_geometry
workflows/diffraction_spots
workflows/xrd_generation
workflows/saed_generation
workflows/crystal_visualization
workflows/style_customization
workflows/combined_structure_diffraction_visualization
workflows/ipf_colors
workflows/xrdml_texture_import
workflows/labotex_open_pole_figures
```

```{toctree}
:maxdepth: 2
:caption: Architecture

architecture/index
```

```{toctree}
:maxdepth: 2
:caption: Standards

standards/index
```

```{toctree}
:maxdepth: 2
:caption: Validation And Research

validation/index
benchmarks/index
theory/index
reference/canonical_docs
api/index
```

## Why PyTex Exists

Many texture and diffraction workflows are scientifically fragile at tool boundaries. A rotation
matrix, a quaternion, or an EBSD map is not self-describing unless the associated frame, symmetry,
phase, and convention choices remain attached to it. PyTex treats those choices as part of the
stable product surface.

This leads to four concrete design commitments:

- semantic types before convenience arrays
- one shared frame and symmetry model across texture, EBSD, diffraction, and transformation work
- validation and benchmark manifests treated as part of feature completion
- teaching-grade and research-grade documentation served from the same browsable site

## Current Scope

- Stable core primitives for frames, transforms, symmetry, lattice semantics, provenance, batch
  semantics, rotations, orientations, Miller objects, and transformation records.
- Foundational texture-domain support for pole figures, inverse pole figures, IPF color keys, and
  discrete kernel ODF evaluation.
- Runtime semantic plotting support for vectors, symmetry elements, EBSD maps, PF/IPF objects, and
  discrete ODF views.
- Multiphase and graph-backed EBSD workflows with phase-aware normalization, KAM, grain
  segmentation, boundary extraction, and direct ODF/PF/IPF outputs from normalized crystal maps.
- Dedicated manifest families for import, experiment, benchmark, validation, workflow result, and
  transformation interchange.
- Bounded experimental parent-candidate scoring staged under `pytex.experimental` rather than
  overstated as a stable reconstruction API.
- Diffraction-domain containers for geometry and pattern data, with powder XRD generation, SAED
  spot simulation, and pinned external-baseline starter cases.

## Recommended Reading Path

### If You Are New To Texture Or EBSD

1. Start with {doc}`tutorials/quickstart`.
2. Read {doc}`concepts/core_model` and {doc}`concepts/technical_glossary_and_symbols`.
3. Continue to {doc}`concepts/orientation_texture` and {doc}`concepts/ebsd_foundation`.
4. Then move to {doc}`workflows/ebsd_import_normalization`,
   {doc}`workflows/ebsd_kam`, and {doc}`workflows/ebsd_to_texture_outputs`.

### If You Are Evaluating Scientific Rigor

1. Read {doc}`architecture/index`.
2. Read {doc}`standards/index`.
3. Check {doc}`validation/index` and {doc}`benchmarks/index`.
4. Follow the linked theory notes in {doc}`theory/index`.

```{note}
PyTex is designed to serve both research and teaching. The Sphinx pages explain concepts and
workflows; the LaTeX notes under `docs/tex/` remain the canonical deep scientific source for
theory, algorithms, and validation posture.
```

![Reference Frames](../figures/reference_frames_vectors.svg)

## Scientific Posture

PyTex does not treat external-tool agreement as a vague aspiration. For covered areas, MTEX is the
validation floor, not the ceiling. Fixture-backed parity tests, benchmark manifests, validation
ledgers, and theory notes are all part of the definition of a stable feature.
