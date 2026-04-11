# Reference Index

This index is the discovery surface for the PDF corpus in `references/`.
Future tasks should use it together with `formulation_summary.md` before mining the PDFs directly.

| Aspect | Reference PDF | Page number(s) | Remarks |
| --- | --- | --- | --- |
| Rotation representations overview and convention risks | `MathsOfrotations_RolletDegraef.pdf` | pp. 1-3 | Why sign, frame, and representation conventions drift; good basis for PyTex orientation docs. |
| Right-handed Cartesian frames | `MathsOfrotations_RolletDegraef.pdf` | p. 3 | Explicitly fixes right-handed frame usage. |
| Passive rotation interpretation | `MathsOfrotations_RolletDegraef.pdf` | pp. 4-5 | Useful when documenting orientation matrix meaning and frame mapping direction. |
| Bunge Euler convention and angle ranges | `MathsOfrotations_RolletDegraef.pdf` | pp. 5, 21 | Gives the `zxz` Bunge convention and explicit matrix form. |
| Axis-angle and quaternion conventions | `MathsOfrotations_RolletDegraef.pdf` | pp. 6-7, 23-25 | Good source for quaternion storage/order and axis-angle formulas. |
| Quaternion rotation operator and composition | `MathsOfrotations_RolletDegraef.pdf` | pp. 7-8, 14-17 | Useful for rotation implementation audits and docs. |
| Orientation descriptors overview | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | p. 34 | Overview of rotation matrix, Miller/Miller-Bravais, Euler, angle-axis, Rodrigues. |
| Specimen and crystal coordinate systems | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 35-37 | Strong source for frame-definition figures and docs. |
| Orientation matrix construction | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 37-40 | Good for orientation-matrix documentation and worked examples. |
| Miller or Miller-Bravais ideal-orientation notation | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 41-43 | Useful when relating orientation labels to matrix form. |
| Pole figure construction and specimen-frame angles | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 43-46 | Basis for PF/IPF docs and redraws. |
| Inverse pole figure meaning | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 45-46 | Key source for IPF figures and reduction docs. |
| Bunge Euler angle definition and matrix relation | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 47-49 | Useful for PyTex Euler docs and parity explanation. |
| Euler space, asymmetric unit, and symmetry reduction | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 49-50 | Basis for Euler-space and fundamental-region documentation. |
| Harmonic ODF series-expansion method | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 105-107 | Primary source for symmetrized generalized spherical harmonics, PF/ODF coefficient linkage, truncation, and odd-order limitations in diffraction PFs. |
| Harmonic ghost-error discussion | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 107-109 | Useful when documenting why antipodal PFs motivate conservative even-degree defaults. |
| Fixed-sample inverse pole figure measurement idea | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | p. 103 | Useful for diffraction-texture workflow ideas. |
| Kikuchi pattern overview | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 148-151 | Good for EBSD/TEM Kikuchi documentation. |
| Kikuchi bands, zone axes, and pattern content | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 150-151 | Grounds Kikuchi-band and indexing figures. |
| Kikuchi pattern as gnomonic projection | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | p. 153 | Important for gnomonic projection and pattern-center figures. |
| Qualitative Kikuchi interpretation | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 153-155 | Useful for pattern-quality and boundary-analysis docs. |
| Orientation determination from Kikuchi patterns | `Introduction_to_Texture_Analysis__Macrotexture_Microtexture_and_Orientation_Mapping.pdf` | pp. 157-159 | Strong source for orientation-indexing workflow docs. |
| Hexagonal plane 4-index rule `i = -(h + k)` | `hexagnoal 4index mathematics.pdf` | pp. 1-2 | Clear derivation for plane notation. |
| Hexagonal direction 4-index rule and inverse transform | `hexagnoal 4index mathematics.pdf` | pp. 2-3 | Direct source for `uvw <-> UVTW` conversions. |
| Four-index zone equation | `hexagnoal 4index mathematics.pdf` | pp. 2-4 | Useful for future `zone_law` helpers in Miller-Bravais notation. |
| Worked four-index examples and addition/subtraction rules | `hexagnoal 4index mathematics.pdf` | pp. 3-4 | Ready-made test and documentation examples. |
| Gnomonic construction for hexagonal zone symbols | `hexagnoal 4index mathematics.pdf` | pp. 5-7 | Useful for crystallographic projection figures. |
| Miller planes from intercepts | `crystallographY_calcualtions.pdf` | book pp. 9-10 (PDF pp. 19-20) | Clear derivation of Miller indices from intercepts. |
| Reciprocal basis definition `a_i . a_j* = delta_ij` | `crystallographY_calcualtions.pdf` | book pp. 10-11 (PDF pp. 20-21) | Good source for basis and reciprocal-lattice docs. |
| Reciprocal vector normal to `(hkl)` | `crystallographY_calcualtions.pdf` | book pp. 11-13 (PDF pp. 21-22) | Grounds plane-normal and reciprocal-vector semantics. |
| `|g_hkl| = 1 / d_hkl` | `crystallographY_calcualtions.pdf` | book pp. 13-14 (PDF pp. 22-23) | Essential for d-spacing and diffraction docs/tests. |
| Reciprocal metric tensor and direct/reciprocal transforms | `crystallographY_calcualtions.pdf` | book pp. 14-18 (PDF pp. 23-27) | Strong basis for metric-tensor feature ideas. |
| Non-Cartesian cross product in crystallographic bases | `crystallographY_calcualtions.pdf` | book pp. 18-20 (PDF pp. 27-29) | Useful for advanced basis and reciprocal-space work. |
| Lattice geometry, directions, and Weiss zone law | `Kelly & Groves.pdf` | pp. 3-11 | Good general crystallography reference for directions/planes/zones. |
| Hexagonal point-group/system overview | `Kelly & Groves.pdf` | pp. 56-59 | Useful for symmetry-facing hexagonal documentation. |
| Crystallographic calculations, reciprocal lattice, and matrices | `Kelly & Groves.pdf` | Appendix 1, pp. 435-449 | Broad reference for vector algebra, reciprocal lattice, rotations, quaternions. |
| Stereographic projection constructions | `Kelly & Groves.pdf` | Appendix 2, pp. 451-468 | Strong source for PF/IPF/reference-frame figures. |
| Interplanar spacings, including hexagonal formula | `Kelly & Groves.pdf` | Appendix 3, pp. 469-472 | Useful for exact d-spacing formulas in docs/tests. |
| Hexagonal index transformations after unit-cell changes | `Kelly & Groves.pdf` | Appendix 4, pp. 473-477 | Useful for future lattice-setting conversion helpers. |
| Basis definition, reciprocal basis, and metrics in worked-example style | `Bhadesia_crystallography.pdf` | pp. 2-6 | Practical source for basis-transform and reciprocal-basis docs. |
| Orientation relations between crystals | `Bhadesia_crystallography.pdf` | pp. 12-24 | Useful for future orientation-relationship and variant tools. |
| Homogeneous deformation, stretch and rotation, invariant-plane strain | `Bhadesia_crystallography.pdf` | pp. 25-49 | Good source for transformation and interface features. |
| Bragg and Laue diffraction equations | `williamsandcarter.pdf` | pp. 78-79 | Good concise source for diffraction-angle documentation. |
| Angle conventions in TEM diffraction geometry | `williamsandcarter.pdf` | pp. 64-65 | Useful for detector and scattering-angle docs. |
| Specimen orientation from diffraction patterns | `williamsandcarter.pdf` | pp. 317-318 | Grounds future zone-axis and indexing helpers. |
| Kikuchi maps for fcc, bcc, and c.p.h. crystals | `kikuchi maps of cubic and hexgonal crystals.pdf` | Appendix 2, pp. 4-11 | Reference corpus for Kikuchi-map redraws and future indexing tools. |

## Notes

- `page number(s)` refers to the printed page where the OCR or outline makes that reliable.
- Where the PDF is only a partial extract, the remark includes both the printed page and the PDF page when needed.

## References

### Normative

- [Reference Canon](../docs/standards/reference_canon.md)
- [Testing Strategy](../docs/testing/strategy.md)

### Informative

- [Formulation Summary](formulation_summary.md)
- [Feature Opportunities](feature_opportunities.md)
