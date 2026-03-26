# Mission

## Vision

Build a world-class Python library for texture-led, crystallography-centered, multimodal materials characterization.

PyTex should give materials scientists, educators, and method developers one coherent environment for orientation analysis, texture analysis, EBSD post-processing, diffraction geometry, CIF-backed phase semantics, and future phase-transformation workflows across EBSD, XRD, neutron diffraction, and TEM.

## Product Position

PyTex is not intended to be a thin Python clone of any single external tool. It should combine:

- the scientific ambition and validation discipline associated with MTEX
- the Python accessibility associated with ORIX
- the EBSD ecosystem connectivity associated with KikuchiPy and PyEBSDIndex
- the crystallographic structure interoperability associated with pymatgen
- the geometric and diffraction breadth needed for XRD, neutron, and TEM reasoning

The differentiator is not breadth alone. It is semantic coherence across those domains.

## Problem Statement

Current crystallographic characterization workflows are fractured across tools that disagree or stay underspecified about:

- reference-frame chains between crystal, specimen, map, detector, laboratory, and reciprocal spaces
- Euler-angle conventions and orientation mapping direction
- reciprocal-space normalization and scattering geometry
- point-group versus space-group semantics
- vendor-specific EBSD import and normalization behavior
- provenance, calibration state, uncertainty, and reproducibility
- the relationship between teaching diagrams and production code

Researchers and students spend avoidable time verifying whether two orientations, pole figures, diffraction setups, or map coordinates mean the same thing before they can do science.

## Product Mission

Create a scientific library that enables users to:

- represent crystallographic objects with explicit frame, symmetry, basis, calibration, and provenance semantics
- normalize imported data once into canonical PyTex semantics without discarding source-system meaning
- express core texture workflows through named domain primitives rather than naked arrays
- connect EBSD, diffraction, and future phase-transformation workflows through one shared crystallographic model
- generate publication-quality and teaching-quality figures from the same scientifically explicit contracts
- trace every major result back to the conventions, transforms, references, and inputs used to produce it

## Core Outcomes

- Eliminate avoidable frame and convention ambiguity from texture and diffraction workflows.
- Establish a canonical internal data model that can support EBSD, XRD, neutron, TEM, and transformation workflows without subsystem-specific semantic forks.
- Keep scientific documentation, figures, and validation artifacts on equal footing with code.
- Build a validation culture where parity and literature alignment are explicit, versioned, and extendable.
- Support both advanced research and PhD-level teaching without maintaining separate conceptual models.

## Product Principles

1. Explicit semantics over convenience
   No important frame, symmetry, basis, calibration, or provenance assumption should be implicit.

2. Texture-led core model first
   Orientation and texture semantics are the nucleus of the product, but they must scale to multimodal characterization.

3. Documentation as scientific infrastructure
   Sphinx concepts and workflows, LaTeX scientific notes, and SVG figures are required deliverables.

4. Validation before claims
   PyTex should not claim scientific equivalence without traceable parity, benchmark, or literature-backed validation.

5. Adapter-friendly, not adapter-owned
   External libraries may strengthen PyTex, but they must not define the stable PyTex public surface.

6. Teaching and research together
   The same library should support derivation-level explanation, classroom use, and production analysis.

7. Publication quality by design
   Figures, notation, validation notes, and public outputs must be suitable for papers, theses, and lectures.

8. Provenance, calibration, and uncertainty matter
   Scientific meaning is incomplete unless the acquisition and conversion context remain inspectable.

## Non-Negotiable Requirements

- Pure-Python-first project structure
- GPL-compatible licensing posture
- Canonical internal convention set
- Canonical data structures for frames, symmetry, structure, orientations, maps, acquisition geometry, and diffraction geometry
- Sphinx as the primary browsable and searchable documentation surface
- LaTeX as the canonical source for major scientific notes
- SVG-first scientific diagrams where geometry or conventions matter
- MTEX parity as a floor for texture and EBSD validation, with additional PyTex-specific validation beyond it
- Separation between stable scientific APIs and experimental research tracks

## Success Criteria

- A user can import or construct crystallographic data without ambiguity about frames, conventions, or core symmetry meaning.
- Core texture and diffraction workflows can be expressed through PyTex public primitives rather than ad hoc arrays or tool-specific objects.
- Documentation explains not only how to run a workflow, but what mathematical object is being transformed and why.
- Validation artifacts make clear what has been checked against MTEX, literature, or future diffraction baselines.
- Future EBSD, neutron, XRD, TEM, and phase-transformation algorithms can be added without inventing local frame or symmetry semantics.

## References

### Normative

- `docs/standards/reference_canon.md`
- `docs/standards/notation_and_conventions.md`
- `docs/architecture/canonical_data_model.md`

### Informative

- Bunge, H.-J., *Texture Analysis in Materials Science: Mathematical Methods*
- MTEX documentation: <https://mtex-toolbox.github.io/>
