# Full API Reference

This page is the exhaustive Sphinx-rendered API surface for the current PyTex package. It is
generated from the Python modules with `sphinx.ext.autodoc`, so class methods, functions,
properties, dataclass fields, signatures, type hints, inheritance, and currently undocumented
members remain visible in the browsable documentation.

Use {doc}`index` first for the curated scientific entry points. Use this page when you need the
complete callable surface of a module, including lower-level helpers that support stable workflows.

```{note}
Private implementation helpers whose names begin with an underscore are intentionally not expanded
as public API entries. Modules with leading underscores are listed only where they contain reusable
support objects that are imported elsewhere in the library.
```

## Package Surface

```{eval-rst}
.. automodule:: pytex
   :no-members:
```

```{eval-rst}
.. automodule:: pytex.contracts
```

```{eval-rst}
.. automodule:: pytex.cli
```

## Core

```{eval-rst}
.. automodule:: pytex.core
   :no-members:
```

```{eval-rst}
.. automodule:: pytex.core.acquisition
```

```{eval-rst}
.. automodule:: pytex.core.batches
```

```{eval-rst}
.. automodule:: pytex.core.conventions
```

```{eval-rst}
.. automodule:: pytex.core.fixtures
```

```{eval-rst}
.. automodule:: pytex.core.frames
```

```{eval-rst}
.. automodule:: pytex.core.hexagonal
```

```{eval-rst}
.. automodule:: pytex.core.lattice
```

```{eval-rst}
.. automodule:: pytex.core.miller
```

```{eval-rst}
.. automodule:: pytex.core.notation
```

```{eval-rst}
.. automodule:: pytex.core.orientation
```

```{eval-rst}
.. automodule:: pytex.core.orientation_geometry
```

```{eval-rst}
.. automodule:: pytex.core.parent_reconstruction
```

```{eval-rst}
.. automodule:: pytex.core.provenance
```

```{eval-rst}
.. automodule:: pytex.core.representations
```

```{eval-rst}
.. automodule:: pytex.core.symmetry
```

```{eval-rst}
.. automodule:: pytex.core.transformation
```

### Core Support Modules

```{eval-rst}
.. automodule:: pytex.core._arrays
```

```{eval-rst}
.. automodule:: pytex.core._chemistry
```

## Texture

```{eval-rst}
.. automodule:: pytex.texture
   :no-members:
```

```{eval-rst}
.. automodule:: pytex.texture.models
```

```{eval-rst}
.. automodule:: pytex.texture.projections
```

```{eval-rst}
.. automodule:: pytex.texture.reconstruction
```

```{eval-rst}
.. automodule:: pytex.texture.harmonics
```

## EBSD

```{eval-rst}
.. automodule:: pytex.ebsd
   :no-members:
```

```{eval-rst}
.. automodule:: pytex.ebsd.gnd
```

```{eval-rst}
.. automodule:: pytex.ebsd.models
```

```{eval-rst}
.. automodule:: pytex.ebsd.texture_workflow
```

## Diffraction

```{eval-rst}
.. automodule:: pytex.diffraction
   :no-members:
```

```{eval-rst}
.. automodule:: pytex.diffraction.cbed
```

```{eval-rst}
.. automodule:: pytex.diffraction.diffraction_groups
```

```{eval-rst}
.. automodule:: pytex.diffraction.dynamical
```

```{eval-rst}
.. automodule:: pytex.diffraction.holz
```

```{eval-rst}
.. automodule:: pytex.diffraction.kikuchi
```

```{eval-rst}
.. automodule:: pytex.diffraction.models
```

```{eval-rst}
.. automodule:: pytex.diffraction.physics
```

```{eval-rst}
.. automodule:: pytex.diffraction.preferred_orientation
```

```{eval-rst}
.. automodule:: pytex.diffraction.saed
```

```{eval-rst}
.. automodule:: pytex.diffraction.stereonets
```

```{eval-rst}
.. automodule:: pytex.diffraction.xrd
```

## Plotting

```{eval-rst}
.. automodule:: pytex.plotting
   :no-members:
```

```{eval-rst}
.. automodule:: pytex.plotting.builders
```

```{eval-rst}
.. automodule:: pytex.plotting.crystal3d
```

```{eval-rst}
.. automodule:: pytex.plotting.diffraction
```

```{eval-rst}
.. automodule:: pytex.plotting.ebsd
```

```{eval-rst}
.. automodule:: pytex.plotting.ipf
```

```{eval-rst}
.. automodule:: pytex.plotting.runtime
```

```{eval-rst}
.. automodule:: pytex.plotting.spherical
```

```{eval-rst}
.. automodule:: pytex.plotting.styles
```

### Plotting Support Modules

```{eval-rst}
.. automodule:: pytex.plotting._render
```

```{eval-rst}
.. automodule:: pytex.plotting._plotting_validation_cases
```

## Adapters

```{eval-rst}
.. automodule:: pytex.adapters
   :no-members:
```

```{eval-rst}
.. automodule:: pytex.adapters.ebsd
```

```{eval-rst}
.. automodule:: pytex.adapters.kikuchipy_workflows
```

```{eval-rst}
.. automodule:: pytex.adapters.labotex
```

```{eval-rst}
.. automodule:: pytex.adapters.manifests
```

```{eval-rst}
.. automodule:: pytex.adapters.orix
```

```{eval-rst}
.. automodule:: pytex.adapters.orix_miller
```

```{eval-rst}
.. automodule:: pytex.adapters.xrdml
```

## Experimental

```{eval-rst}
.. automodule:: pytex.experimental
   :no-members:
```

```{eval-rst}
.. automodule:: pytex.experimental.phase_transformation
```
