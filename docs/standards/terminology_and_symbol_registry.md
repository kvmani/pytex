# Terminology And Symbol Registry

This document fixes the repository-wide vocabulary and symbol policy for PyTex.

## Purpose

PyTex is large enough that symbol drift and term drift become real scientific risks. This registry exists so the same term and the same symbol keep the same meaning across:

- Sphinx concept pages
- workflow guides
- notebook tutorials
- LaTeX theory and algorithm notes
- code explanations and docstrings where symbols are discussed

## Policy

- Stable scientific terms should be defined once here and reused elsewhere.
- Stable mathematical symbols should be introduced here or in the closest canonical theory note, then reused consistently.
- If a page needs a local symbol extension, it should state that extension explicitly and keep the core registry symbols unchanged.
- Pages that rely on registry terms should link back to the user-facing glossary page and, when needed, to this standards document.
- When a canonical SVG figure labels one of these terms or symbols, it should use the same wording and symbol form as this registry unless the figure explicitly documents a local teaching simplification.

## Core Terms

| Term | Fixed meaning |
| --- | --- |
| reference frame | A named, domain-typed coordinate frame such as crystal, specimen, map, detector, laboratory, or reciprocal. |
| orientation | A crystal-to-specimen mapping carried by an explicit `Orientation` object. |
| rotation | A geometric active rotation that does not by itself define crystallographic source and target meaning. |
| symmetry | Point-group-facing operator set used for orientation and direction reduction. |
| space group | Structure-facing crystallographic identity used for phases and CIF-backed construction. |
| pole figure | Distribution of crystal directions or plane normals expressed relative to specimen directions. |
| inverse pole figure | Distribution of specimen directions expressed in crystal coordinates and reduced by symmetry where appropriate. |
| ODF | Orientation distribution function over orientation space. |
| zone axis | Direct-space crystallographic direction defining an electron-diffraction viewing or incidence condition. |
| powder pattern | Grid-sampled XRD spectrum built from discrete reflections and an optional broadening model. |
| crystal scene | Reusable geometry bundle for 3D crystal rendering. |

## Core Symbols

| Symbol | Meaning |
| --- | --- |
| \(\mathbf{v}\) | Generic vector in an explicitly named frame. |
| \(\mathbf{R}\) | Rotation matrix acting actively on vectors. |
| \(q\) | Unit quaternion in `w, x, y, z` storage order. |
| \((\phi_1, \Phi, \phi_2)\) | Bunge Euler angles. |
| \(\mathbf{a}, \mathbf{b}, \mathbf{c}\) | Direct-lattice basis vectors. |
| \(\mathbf{a}^{*}, \mathbf{b}^{*}, \mathbf{c}^{*}\) | Reciprocal-lattice basis vectors under the PyTex normalization rule. |
| \(\mathbf{g}_{hkl}\) | Reciprocal-lattice vector associated with Miller indices \((hkl)\). |
| \(d_{hkl}\) | Interplanar spacing for the \((hkl)\) family. |
| \(\theta\) | Bragg half-angle. |
| \(2\theta\) | Powder-diffraction scattering angle reported in XRD plots. |
| \(F_{hkl}\) | Reflection structure-factor quantity or current PyTex proxy where explicitly stated. |
| \(\hat{\mathbf{z}}\) | Unit zone-axis direction in direct space. |
| \(u, v\) | Detector-plane plotting coordinates in SAED or detector geometry contexts. |

## References

### Normative

- `notation_and_conventions.md`
- `reference_canon.md`

### Informative

- `../site/concepts/technical_glossary_and_symbols.md`
