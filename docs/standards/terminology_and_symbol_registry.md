# Terminology And Symbol Registry

This document fixes the repository-wide vocabulary and symbol policy for PyTex. It
is the single source of truth for nomenclature: the same term and the same symbol
keep the same meaning across every documentation surface and in code explanations.

## Purpose

PyTex is large enough that symbol drift and term drift become real scientific
risks. This registry exists so the same term and the same symbol keep the same
meaning across:

- Sphinx concept pages
- workflow guides
- notebook tutorials
- executable worked examples
- LaTeX theory and algorithm notes
- code explanations and docstrings where symbols are discussed
- canonical SVG figures

## Policy

- Stable scientific terms should be defined once here and reused elsewhere.
- Stable mathematical symbols should be introduced here or in the closest
  canonical theory note, then reused consistently.
- If a page needs a local symbol extension, it should state that extension
  explicitly and keep the core registry symbols unchanged.
- Pages that rely on registry terms should link back to the user-facing glossary
  page and, when needed, to this standards document.
- When a canonical SVG figure labels one of these terms or symbols, it should use
  the same wording and symbol form as this registry unless the figure explicitly
  documents a local teaching simplification.
- Executable worked examples must reference symbols from this registry; a symbol
  used in a worked example that is not yet registered must be added here first.
  See [Executable Worked Examples](executable_examples.md).

## Symbol Registration Policy

New notation is introduced through this registry, not ad hoc in a single page:

1. Prefer an existing registered symbol.
2. If a new symbol is genuinely needed, add it to the appropriate table below with
   a one-line fixed meaning before using it in prose, math, figures, notebooks,
   worked examples, or docstrings.
3. Do not reuse a registered symbol for a different meaning; choose a distinct
   symbol or an explicitly scoped local extension.
4. Keep storage-order and normalization conventions (quaternion order, reciprocal
   normalization, Euler labeling) identical to those fixed in
   [Notation and Conventions](notation_and_conventions.md).

## Core Terms

| Term | Fixed meaning |
| --- | --- |
| reference frame | A named, domain-typed coordinate frame such as crystal, specimen, map, detector, laboratory, or reciprocal. |
| orientation | A crystal-to-specimen mapping carried by an explicit `Orientation` object. |
| rotation | A geometric active rotation that does not by itself define crystallographic source and target meaning. |
| misorientation | The orientation mapping between two orientations, before symmetry reduction. |
| disorientation | The symmetry-reduced misorientation of minimal rotation angle in the fundamental zone. |
| symmetry | Point-group-facing operator set used for orientation and direction reduction. |
| space group | Structure-facing crystallographic identity used for phases and CIF-backed construction. |
| pole figure | Distribution of crystal directions or plane normals expressed relative to specimen directions. |
| inverse pole figure | Distribution of specimen directions expressed in crystal coordinates and reduced by symmetry where appropriate. |
| ODF | Orientation distribution function over orientation space. |
| fundamental zone | The symmetry-reduced subset of orientation space used for canonical orientation keys. |
| zone axis | Direct-space crystallographic direction defining an electron-diffraction viewing or incidence condition. |
| powder pattern | Grid-sampled XRD spectrum built from discrete reflections and an optional broadening model. |
| multiplicity | Number of symmetry-equivalent members of a plane or direction family under the phase point group. |
| crystal scene | Reusable geometry bundle for 3D crystal rendering. |

## Core Symbols

### Frames, vectors, and rotations

| Symbol | Meaning |
| --- | --- |
| \(\mathbf{v}\) | Generic vector in an explicitly named frame. |
| \(\hat{\mathbf{v}}\) | Unit vector (normalized) in an explicitly named frame. |
| \(\mathbf{R}\) | Rotation matrix acting actively on vectors. |
| \(\mathbf{T}\) | Rigid placement of geometry into a shared world frame: \(\mathbf{T}(\mathbf{x}) = \mathbf{R}\,\mathbf{x} + \mathbf{t}\) (the `Transform3D` visualization primitive). |
| \(q\) | Unit quaternion in `w, x, y, z` storage order. |
| \((\phi_1, \Phi, \phi_2)\) | Bunge Euler angles. |
| \((\mathbf{n}, \omega)\) | Axis-angle pair: rotation axis \(\mathbf{n}\) and angle \(\omega\). |

### Lattice and reciprocal lattice

| Symbol | Meaning |
| --- | --- |
| \(\mathbf{a}, \mathbf{b}, \mathbf{c}\) | Direct-lattice basis vectors. |
| \(a, b, c, \alpha, \beta, \gamma\) | Lattice parameters (edge lengths and angles). |
| \(\mathbf{a}^{*}, \mathbf{b}^{*}, \mathbf{c}^{*}\) | Reciprocal-lattice basis vectors under the PyTex normalization rule. |
| \(\mathbf{g}_{hkl}\) | Reciprocal-lattice vector associated with Miller indices \((hkl)\). |
| \(\mathbf{G}\) | Direct-space metric tensor. |
| \(d_{hkl}\) | Interplanar spacing for the \((hkl)\) family. |

### Miller indices and crystallographic geometry

| Symbol | Meaning |
| --- | --- |
| \((hkl)\) | Miller plane indices; \((hkil)\) for the four-index hexagonal form. |
| \([uvw]\) | Miller direction indices; \([uvtw]\) for the four-index hexagonal form. |
| \(\mathbf{n}\) | Plane normal direction. |
| \(\angle(\mathbf{n}_1, \mathbf{n}_2)\) | Angle between two plane normals (interplanar angle). |
| \(\angle(\mathbf{d}_1, \mathbf{d}_2)\) | Angle between two lattice directions. |
| \(m_{\{hkl\}}\) | Symmetry multiplicity of a plane family under the phase point group. |

### Orientation and misorientation

| Symbol | Meaning |
| --- | --- |
| \(g\) | An orientation (crystal-to-specimen mapping). |
| \(\Delta g\) | A misorientation between two orientations. |
| \(\omega\) | Disorientation angle: minimal misorientation angle over the symmetry group. |
| \(\Sigma\) | Coincidence-site-lattice index of a boundary (for example \(\Sigma 3\)). |

### Diffraction

| Symbol | Meaning |
| --- | --- |
| \(\lambda\) | Radiation wavelength. |
| \(\theta\) | Bragg half-angle. |
| \(2\theta\) | Powder-diffraction scattering angle reported in XRD plots. |
| \(F_{hkl}\) | Reflection structure-factor quantity or current PyTex proxy where explicitly stated. |
| \(\hat{\mathbf{z}}\) | Unit zone-axis direction in direct space. |
| \(u, v\) | Detector-plane plotting coordinates in SAED or detector geometry contexts. |

## References

### Normative

- [Notation and Conventions](notation_and_conventions.md)
- [Executable Worked Examples](executable_examples.md)
- [Reference Canon](reference_canon.md)

### Informative

- <a href="../site/concepts/technical_glossary_and_symbols.md">Technical Glossary and Symbols</a>
