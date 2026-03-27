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
- represent point-group and space-group meaning explicitly, without collapsing structure semantics into orientation reduction semantics
- normalize imported data once into canonical PyTex semantics without discarding source-system meaning
- express core texture workflows through named domain primitives rather than naked arrays
- express large vectorized workflows through first-class batch primitives rather than anonymous `(n, ...)` arrays
- connect EBSD, diffraction, and future phase-transformation workflows through one shared crystallographic model
- carry experiment, benchmark, validation, and workflow-result context through stable machine-readable manifests rather than ad hoc side channels
- generate publication-quality and teaching-quality figures from the same scientifically explicit contracts
- inspect texture results through first-class plotting surfaces, including contour pole figures and ODF section views, without leaving the PyTex semantic model
- generate diffraction and structure visuals through first-class plotting surfaces, including powder XRD plots, SAED spot maps, and VESTA-like 3D crystal scenes
- control visual output through reusable YAML themes and overrides so notebooks, workflow pages, and user scripts share one plotting style system
- learn the library through executable, mathematics-backed notebooks that match the formal documentation instead of drifting away from it
- trace every major result back to the conventions, transforms, references, and inputs used to produce it
- rely on one stable terminology and symbol vocabulary across documentation, theory notes, notebooks, and code-facing explanations rather than redefining notation per page

## Core Outcomes

- Eliminate avoidable frame and convention ambiguity from texture and diffraction workflows.
- Treat vectorized operations on scientific primitives as a first-class capability rather than a convenience layer.
- Establish a canonical internal data model that can support EBSD, XRD, neutron, TEM, and transformation workflows without subsystem-specific semantic forks.
- Keep scientific documentation, figures, and validation artifacts on equal footing with code.
- Build a validation culture where parity and literature alignment are explicit, versioned, and extendable.
- Support both advanced research and PhD-level teaching without maintaining separate conceptual models.

## Product Principles

1. Explicit semantics over convenience
   No important frame, symmetry, basis, calibration, or provenance assumption should be implicit.

1. Batch semantics are part of the product
   Vectorized operations must preserve frame, convention, symmetry, and provenance meaning through first-class batch primitives.

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
- Canonical data structures for frames, symmetry, structure, orientations, vectorized batch primitives, maps, acquisition geometry, and diffraction geometry
- Sphinx as the primary browsable and searchable documentation surface
- LaTeX as the canonical source for major scientific notes
- SVG-first scientific diagrams where geometry or conventions matter
- MTEX parity as a floor for texture and EBSD validation, with additional PyTex-specific validation beyond it
- Separation between stable scientific APIs and experimental research tracks

## Success Criteria

- A user can import or construct crystallographic data without ambiguity about frames, conventions, or core symmetry meaning.
- Core texture and diffraction workflows can be expressed through PyTex public primitives rather than ad hoc arrays or tool-specific objects.
- A user can perform large batched operations on vectors, Euler angles, quaternions, rotations, and orientations without dropping semantic metadata.
- Documentation explains not only how to run a workflow, but what mathematical object is being transformed and why.
- Notebook tutorials provide executable concept explanations for the implemented surface without becoming an informal competing documentation layer.
- Plotting documentation and runtime plotting APIs stay aligned, so the same semantic objects can drive notebooks, workflow pages, and user figures.
- Runtime plotting, documentation figure generation, and YAML theme governance remain aligned, so XRD, SAED, PF/IPF/ODF, and crystal-structure visuals all reuse one styling and rendering doctrine.
- Validation artifacts make clear what has been checked against MTEX, literature, or future diffraction baselines.
- Stable workflows can emit or consume manifests for experiment context, benchmark identity, validation evidence, and workflow outputs.
- Future EBSD, neutron, XRD, TEM, and phase-transformation algorithms can be added without inventing local frame or symmetry semantics.

## References

### Normative

- `docs/standards/reference_canon.md`
- `docs/standards/notation_and_conventions.md`
- `docs/architecture/canonical_data_model.md`

### Informative

- Bunge, H.-J., *Texture Analysis in Materials Science: Mathematical Methods*
- MTEX documentation: <https://mtex-toolbox.github.io/>
