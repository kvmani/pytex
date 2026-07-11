# Library Structure

This page is the architecture atlas for PyTex.

:::{raw} html
<p class="architecture-lead">
It is meant to give a reader a fast, accurate mental model of the library before they dive into implementation details. The diagrams below do not merely mirror the folder tree. They explain how PyTex is organized semantically, how scientific data moves through the system, how documentation and validation govern the codebase, and how the current implementation relates to the next planned foundations.
</p>
:::

## At A Glance

PyTex is organized around four non-negotiable ideas:

- one canonical crystallographic core
- multiple domain modules that reuse that core
- adapters that normalize external data into PyTex semantics
- documentation and validation as first-class product surfaces

That structure is what lets the project remain texture-led while expanding toward EBSD, diffraction, and future phase-transformation and multimodal workflows.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Semantic Nucleus
:class-card: sd-shadow-sm
PyTex is built around one canonical crystallographic core for frames, symmetry, lattice, orientation, reciprocal-space, and provenance semantics.
:::

:::{grid-item-card} Domain Expansion
:class-card: sd-shadow-sm
Texture, EBSD, diffraction, and plotting are downstream scientific layers that must reuse the same core model rather than inventing private conventions.
:::

:::{grid-item-card} Boundary Discipline
:class-card: sd-shadow-sm
Adapters live at the normalization boundary. They are responsible for translation into PyTex semantics, not for defining the stable product surface.
:::

:::{grid-item-card} Completion Doctrine
:class-card: sd-shadow-sm
Documentation, figures, tests, ledgers, and benchmark placeholders are architectural surfaces. They are part of feature completion, not release polish.
:::
::::

:::{figure} ../../figures/pytex_architecture_poster.svg
:alt: PyTex architecture poster
:class: architecture-poster-figure
:::

:::{figure} ../../figures/pytex_architecture_evolution_poster.svg
:alt: PyTex architecture evolution poster
:class: architecture-poster-figure
:::

## Poster Series

The poster below is designed as the primary reusable figure for talks, README summaries, and high-level architecture discussions.

It is intentionally poster-like: it compresses the whole project into a single view while preserving the semantic hierarchy of governance, core model, domain modules, boundary layers, and next foundations.

## Visual Legend

::::{grid} 1 1 4 4
:gutter: 2

:::{grid-item-card} Governance
:class-card: sd-shadow-sm
Mission, standards, documentation, and validation layers that control what the code is allowed to claim.
:::

:::{grid-item-card} Core
:class-card: sd-shadow-sm
Canonical objects and semantics shared across all scientific workflows.
:::

:::{grid-item-card} Domain
:class-card: sd-shadow-sm
Texture, EBSD, diffraction, and presentation-facing modules built on the core.
:::

:::{grid-item-card} Growth
:class-card: sd-shadow-sm
Boundary layers and future foundations that extend the system without weakening the existing contracts.
:::
::::

## 1. System Structure

:::{raw} html
<p class="architecture-section-intro">
This first diagram is the top-level structural map. It shows the code structure as a semantic system rather than as a bare directory tree.
</p>
:::

:::{figure} ../../figures/pytex_system_structure.svg
:alt: Layered PyTex architecture map showing governance, canonical core, domain modules, boundary layers, and next foundations.
:class: architecture-poster-figure
:::

### Reading This View

- The canonical core is the semantic nucleus of the library.
- Texture, EBSD, diffraction, and plotting are downstream scientific layers, not separate semantic islands.
- Adapters are boundary modules. They translate into PyTex semantics rather than owning them.
- Documentation and validation sit above the codebase because in PyTex they define what counts as a finished feature.
- The future-foundations block is not speculative fluff. It is the documented expansion path for the next stable API families.

## 2. Scientific Data Flow

:::{raw} html
<p class="architecture-section-intro">
This second diagram shows how scientific information is intended to move through the library: normalization at the boundary, canonical objects in the middle, domain workflows downstream, and documentation plus validation surrounding the flow.
</p>
:::

:::{figure} ../../figures/pytex_scientific_data_flow.svg
:alt: Scientific data-flow diagram showing external sources normalized into canonical PyTex objects, domain workflows, outputs, documentation, and validation.
:class: architecture-poster-figure
:::

### Reading This View

- PyTex is designed around one-way normalization into canonical objects.
- Internal algorithms should operate on PyTex domain types, not on source-native objects.
- Documentation and validation are attached to the flow, not bolted on afterward.
- Outputs are expected to remain traceable back to conventions, inputs, and normalization boundaries.
- The data-flow model is intentionally conservative because scientific clarity matters more than API cleverness.

## 3. Governance And Completion Model

:::{raw} html
<p class="architecture-section-intro">
This third diagram explains why PyTex development is deliberately constrained by standards, architecture docs, and validation policy.
</p>
:::

:::{figure} ../../figures/pytex_governance_completion_model.svg
:alt: Governance and completion model showing how scientific capabilities pass through standards, architecture, implementation, documentation, validation, and claim gates.
:class: architecture-poster-figure
:::

### Reading This View

- PyTex does not treat implementation alone as completion.
- Standards and architecture documents are upstream decision surfaces.
- A stable feature claim is only valid after docs and validation exist alongside the code.
- This is the mechanism that keeps the project from drifting as it grows beyond texture into broader crystallography.

## 4. Current State Versus Planned Expansion

:::{raw} html
<p class="architecture-section-intro">
This final diagram separates the current implemented foundation from the next layers that are implemented foundationally today versus still ahead as broader expansion work.
</p>
:::

:::{figure} ../../figures/pytex_current_state_expansion.svg
:alt: Current state versus planned expansion diagram separating implemented foundations from next validation and workflow hardening targets.
:class: architecture-poster-figure
:::

### Reading This View

- PyTex already has a substantial implemented foundation.
- The next major work is not “add more features blindly.” It is to deepen validation and workflow breadth on top of the now-stronger semantic foundation.
- The gap between “implemented today” and “planned next” is intentionally documented so future contributors do not guess.

## Summary

The architecture of PyTex can be reduced to four ideas:

- one canonical crystallographic core
- multiple downstream scientific domains
- adapters at the boundary, not in the center
- documentation and validation as first-class architectural surfaces

That combination is what makes PyTex a scientific infrastructure project rather than just a collection of analysis routines.

## Related Material

- {doc}`core_model`
- {doc}`core_foundation`
- {doc}`orientation_texture`
- {doc}`texture_foundation`
- {doc}`ebsd_foundation`
- {doc}`diffraction_foundation`
- {doc}`../architecture/overview`
- {doc}`../architecture/multimodal_characterization_foundation`
- {doc}`../architecture/phase_transformation_foundation`

## References

### Normative

- {doc}`../architecture/overview`
- {doc}`../architecture/canonical_data_model`
- {doc}`../standards/reference_canon`

### Informative

- {doc}`../architecture/repo_review_2026_foundation_audit`
