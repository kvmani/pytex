# How PyTex Differs

This page is deliberately prominent because PyTex is easy to misunderstand if it is compared only by checklist count.

PyTex is not trying to be:

- a clone of MTEX
- a thin Python wrapper around someone else’s data model
- a vendor-tied EBSD analysis shell
- a plotting-first texture package with hidden conventions

PyTex is trying to build a single explicit scientific model that can be shared across orientations, texture, EBSD, diffraction, documentation, and validation.

## The Short Version

The distinctive PyTex emphasis is:

- one canonical frame and symmetry model across texture, EBSD, and diffraction
- stable scientific APIs centered on named domain primitives rather than unnamed arrays
- machine-readable manifests and provenance as part of workflow stability
- research-grade and teaching-grade documentation treated as part of feature completion
- pure-Python-first implementation with optional adapters at the boundaries

## Feature Table

The table below is intentionally conservative. It compares documented emphasis and public product surface, not marketing claims about who is “best.”

| Capability or emphasis | PyTex | MTEX | orix | LaboTex | TSL / EDAX OIM Analysis |
| --- | --- | --- | --- | --- | --- |
| Primary environment | Python | MATLAB toolbox | Python | Commercial desktop software | Commercial desktop software |
| Explicit phase creation from CIF / crystal structures in the public workflow | Yes, through canonical `Phase` / `UnitCell` constructors | Documented crystal symmetry and structure import workflows | Documented `Phase` plus structure / space-group integration | Texture-analysis workflows are documented; structure import is not the public centerpiece in the same way | Phase handling is part of EBSD workflows, but not exposed as a Pythonic scientific core model |
| EBSD map workflows | Yes | Yes | Yes | Documented EBSD support | Yes |
| Texture-domain ODF / PF focus | Yes | Yes, major focus | More orientation / EBSD focused in public docs | Yes, major focus | Not the primary public identity |
| Diffraction-domain primitives in the same library surface | Yes | Limited compared with PyTex’s explicit diffraction-core direction | Limited compared with PyTex’s diffraction-core direction | Not the public focus | EBSD indexing and analysis are the public focus |
| One shared canonical frame/symmetry model across texture, EBSD, and diffraction as a primary design rule | Core design commitment | Not the central public framing | Not the central public framing | Not the central public framing | Not the central public framing |
| Machine-readable manifests and schemas as part of stable workflow interchange | Core design commitment | Not the public centerpiece | Not the public centerpiece | Not the public centerpiece | Vendor workflow files, but not an open manifest-first scientific contract |
| Sphinx + LaTeX + SVG teaching/research docs treated as part of feature completion | Core design commitment | Strong documentation, different documentation architecture | Strong Python docs | Product documentation and manuals | Product documentation and manuals |
| Open-source, GPL-compatible pure-Python-first core | Yes | Open source, MATLAB-based | Open source, Python-based | Commercial | Commercial |

## What This Means In Practice

### Compared With MTEX

MTEX is the most important parity baseline for many PyTex texture and EBSD semantics, and PyTex explicitly treats MTEX as a validation floor where applicable.

But PyTex is architecturally different in two ways:

- PyTex is Python-native rather than MATLAB-native
- PyTex is built around a single canonical data model and machine-readable workflow contracts, not only around analysis routines

So PyTex should be read as “Python scientific infrastructure with MTEX-grade parity targets where relevant,” not as “MTEX rewritten line-for-line.”

### Compared With orix

orix is the closest open Python relative in spirit. It documents crystal maps, phases, orientations, and structure-aware phase metadata very well.

PyTex differs mainly in scope and product posture:

- PyTex is explicitly trying to unify texture, EBSD, and diffraction under one canonical model
- PyTex treats manifests, provenance, and theory/figure assets as part of the stable surface
- PyTex is stricter about making frames, bases, and symmetry contracts visible in public APIs

### Compared With LaboTex

LaboTex is a long-standing texture-analysis environment with strong practical workflows around pole figures and ODF-centered texture analysis.

PyTex differs by being:

- Python-native and scriptable as a library-first surface
- explicit about typed frames, symmetry, and provenance in public APIs
- designed to support diffraction and EBSD canonical modeling in the same codebase rather than focusing only on texture workflows

### Compared With TSL / EDAX OIM Analysis

OIM Analysis is a mature commercial EBSD analysis environment centered on indexing, cleanup, grain analysis, and industrial or lab workflows.

PyTex differs by being:

- an open scientific library rather than a closed desktop analysis product
- designed around transparent scientific primitives and documentation assets
- intended to make the canonical model reusable across research code, reproducible workflows, and teaching materials

## The Most Important Unique Aspect

The most important PyTex claim is not “more features than everyone else.” That would be false today and not useful anyway.

The real claim is narrower and stronger:

PyTex makes semantic explicitness itself part of the product surface.

That means:

- reference-frame geometry is documented, not implied
- symmetry reduction rules are documented, not hidden
- workflow interchange is meant to become schema-backed, not notebook folklore
- scientific docs, figures, tests, and parity notes are part of completion, not post hoc extras

## Sources Used For This Comparison

- MTEX documentation: [mtex-toolbox.github.io](https://mtex-toolbox.github.io/)
- orix documentation: [orix.readthedocs.io](https://orix.readthedocs.io/)
- LaboTex material from LaboSoft: [labosoft.com.pl](https://labosoft.com.pl/)
- EDAX OIM Analysis product material: [edax.com](https://www.edax.com/)
- pymatgen documentation for CIF and structure handling: [pymatgen.org](https://pymatgen.org/)
