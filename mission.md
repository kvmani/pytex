# Mission

## Vision

Build a world-class Python library for crystallographic texture and diffraction that gives materials scientists and educators a coherent, rigorously documented, publication-ready environment for orientation analysis, texture analysis, EBSD post-processing, stereographic reasoning, and diffraction workflows.

PyTex should combine the scientific ambition associated with MTEX, the Python accessibility of ORIX, the EBSD ecosystem connectivity of KikuchiPy and PyEBSDIndex, and the geometric breadth seen in diffraction-focused tools, while correcting the fragmentation, convention drift, and documentation gaps that still slow real scientific work.

## Problem Statement

Current texture and diffraction workflows are often split across multiple tools with different assumptions about:

- reference frames
- Euler angle conventions
- reciprocal-space normalization
- symmetry representations
- vendor-specific EBSD import semantics
- plotting standards
- provenance and reproducibility

Researchers and students routinely spend avoidable time translating between these systems, verifying whether two orientations mean the same thing, and rebuilding figures or derivations for publications and teaching.

## Product Mission

Create a scientific library that enables users to:

- represent crystallographic objects with explicit frame, symmetry, basis, and provenance semantics
- convert safely between vendor and software conventions without hidden ambiguity
- compute and visualize pole figures, inverse pole figures, and ODFs in a reproducible Python workflow
- analyze EBSD-derived orientation maps with shared core semantics
- model diffraction geometry and connect texture analysis to diffraction reasoning
- produce publication-quality figures and teaching-quality explanatory material from the same workflows
- trace every major scientific result back to the conventions, inputs, and transforms used to produce it

## Core Outcomes

- Eliminate avoidable reference-frame and convention confusion from scientific workflows.
- Provide a canonical internal data model for texture and diffraction work in Python.
- Make scientific documentation equal in quality to the code.
- Build a validation culture where parity with established tools is explicit and extendable.
- Support both advanced research and classroom explanation without forking the conceptual model.

## Product Principles

1. Explicit semantics over convenience
   No important frame, symmetry, or basis assumption should be implicit.

2. Core model first
   Domain primitives are part of the product, not implementation detail.

3. Documentation as scientific infrastructure
   Sphinx concepts and tutorials, LaTeX scientific notes, and SVG figures are required deliverables.

4. Validation before claims
   PyTex should not claim equivalence to established tools without traceable evidence.

5. Adapter-friendly, not adapter-owned
   External libraries should strengthen PyTex, but not define PyTex public semantics.

6. Teaching and research together
   The same library should support both derivation-level explanation and production analysis.

7. Publication quality by design
   Figures, notation, and outputs must be fit for papers, theses, and lectures.

8. Scientific provenance matters
   Inputs, conventions, transforms, and software context must remain inspectable.

## Non-Negotiable Requirements

- Pure-Python-first project structure
- GPL-compatible licensing posture
- Canonical internal convention set
- Canonical data structures for frames, symmetry, orientations, maps, and diffraction geometry
- Sphinx as the primary browsable/searchable documentation surface
- LaTeX as the canonical source for major scientific notes
- SVG-first scientific diagrams where geometry matters
- MTEX parity matrix and stronger PyTex-specific validation
- Separation between stable scientific APIs and experimental research tracks

## Success Criteria

- A user can import or construct crystallographic data without ambiguity about frames or conventions.
- Core texture workflows can be expressed through PyTex public primitives rather than ad hoc arrays.
- Documentation explains not only how to run a workflow, but what the mathematics and conventions mean.
- Validation artifacts make it clear what has been checked against MTEX, ORIX, or literature references.
- Future algorithms can be added without reinventing frame or symmetry semantics.
