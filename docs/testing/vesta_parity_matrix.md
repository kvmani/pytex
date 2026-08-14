# VESTA Parity Matrix

This document is the authoritative ledger for PyTex crystal-structure visualization parity against
VESTA (Visualization for Electronic and STructural Analysis), the reference desktop application for
publication-grade crystal rendering.

Reference baseline:

- VESTA 3.x feature set as documented in the VESTA manual (Momma & Izumi, J. Appl. Cryst. 44 (2011)
  1272-1276)
- PyTex surface: `pytex.plotting.crystal3d` (single crystal), `pytex.plotting.scene3d` (composite
  world scenes), `pytex.plotting.primitives` (geometric primitives)

The goal is stated in the module contract: *VESTA-class* rendering quality, plus capabilities a
scriptable library can offer that a GUI cannot (composition, placement transforms, regression
testing, theme systems). PyTex does not aim to reproduce VESTA's GUI or its electronic-structure
(volumetric-data) features; those rows are marked `n/a` with reasons.

## Status Keys

- `implemented`: PyTex has the capability with automated regression coverage
- `partial`: a usable subset exists; the gap is described in Notes
- `planned`: accepted target, not yet implemented
- `exceeded`: PyTex offers a strictly stronger capability than the VESTA equivalent
- `n/a`: out of PyTex scope, with explanation

## Matrix

| Feature | VESTA behavior | PyTex status | Notes |
| --- | --- | --- | --- |
| Ball-and-stick rendering | lit spheres + cylinders | implemented | `render_style="ball_and_stick"` (default). Blinn-Phong-style lighting and one globally depth-sorted mesh give publication figures correct atom/bond occlusion (`plot_crystal_structure_3d`); the shared desktop/web workbench uses responsive radial sphere shading and layered cylindrical bonds for the same visual reading while the camera moves. |
| Space-filling rendering | full-radius spheres | implemented | `render_style="space_filling"`: Slater atomic radii (`atomic_radius_angstrom`), bonds suppressed. Ionic-radius selection per oxidation state is planned. |
| Polyhedral rendering | coordination polyhedra | implemented | `render_style="polyhedral"` auto-selects every ≥4-coordinated species; `polyhedra_species` narrows it. Convex-hull faces with outward normals ride in the depth-sorted mesh. Polyhedron volume/distortion metrics are planned. |
| Stick rendering | uniform thin cylinders | implemented | `render_style="stick"`: bond cylinders and atom caps share one radius. |
| Wireframe rendering | line bonds | implemented | `render_style="wireframe"`: the bond network as lines only, atom bodies hidden (`atom_render_mode="none"`), no lit mesh — the VESTA wireframe convention. |
| Thermal displacement ellipsoids | anisotropic ADP ellipsoids | planned | The data model carries isotropic `b_iso` only; anisotropic ADP support in `AtomicSite` must land first. |
| Partial occupancy display | pie-sliced spheres | implemented | Sites sharing one position (mixed species) render as azimuthal sectors of one shared-radius sphere; occupancy below one leaves a vacancy sector in the theme `vacancy_color` (VESTA's white). Automatic from `AtomicSite.occupancy`. |
| Boundary atoms | translated copies complete the cell | implemented | `include_boundary_atoms=True` (default), with coincident-duplicate removal. |
| Multiple unit cells | supercell blocks | implemented | `repeats=(nx, ny, nz)` with per-cell overlays (`show_unit_cells`, `CrystalCellOverlay`). |
| Hexagonal cell outline | hexagonal prism | implemented | `CrystalCellOverlay(kind="hexagonal_prism")` with lattice-consistency validation. |
| Lattice (hkl) cut planes | translucent plane sections | implemented | `plane_hkls` / `CrystalPlaneOverlay` with offset control, polygon clipping to the cell block, and `(hkl)` labels. |
| Crystallographic direction arrows | vectors through the cell | implemented | `CrystalDirectionOverlay` with `[uvw]` labels, plus the composable `direction_arrow` primitive. |
| Vectors on atoms (moments, displacements) | per-site arrows | implemented | `site_vectors={label: vector}` draws the arrow on every periodic copy of the site (crystal-Cartesian angstrom). |
| Atom labels | per-atom text | implemented | `atom_label_mode="species"` or `"site"`. |
| Two-tone bonds | each half in its atom color | implemented | Default `bond_color_mode="two_tone"`; uniform mode available. |
| Dashed / hydrogen-bond styles | per-bond-type line styles | planned | Bonds currently share one solid style per scene. |
| Depth cueing (fog) | distance fade | implemented | `depth_cue_strength` theme key; fades mesh faces toward the background along the view direction. The workbench exposes the strength and recomputes the cue on every camera change; publication export recomputes it for the exported view. |
| Lighting model | ambient/diffuse/specular | implemented | Theme-controlled ambient, diffuse, specular strengths, shininess, and light direction. The workbench exposes those controls with glossy, matte, and flat presets, and transforms its screen-space light into the crystal frame for camera-aware export. |
| Perspective / orthographic projection | both | implemented | `projection="persp"` / `"ortho"`. |
| View along crystal axes / directions | align view to [uvw] or a,b,c | implemented | `view_preset="a"/"b"/"c"`, `view_direction=` vector or `CrystalDirection`. |
| Distance measurement | interactive click readout | exceeded | Programmatic: `CrystalScene.bond_lengths_angstrom()` and `bond_length_summary()` (per species-pair count/min/mean/max) are scriptable and regression-tested rather than click-driven. |
| Angle / dihedral measurement | interactive readout | planned | Bond-angle and dihedral helpers around a coordination center. |
| Symmetry expansion from asymmetric unit | space-group generation | partial | CIF import via pymatgen expands symmetry-equivalent sites (`Phase.from_cif`); native space-group generation without pymatgen is not implemented. |
| Volumetric data (electron density, isosurfaces) | isosurfaces and 2D slices | n/a | Electronic-structure data display is outside the texture/diffraction library scope. |
| Interactive GUI editing | click-driven scene edits | n/a | PyTex is a scriptable library; matplotlib's 3D interactivity applies, and reproducible scene *code* replaces GUI state. |
| Image export | raster screenshots + vector | implemented | `export_figure` writes SVG/PDF/PNG at publication DPI under YAML themes; SVG text stays editable. |
| Multiple structures in one figure | one structure per window | exceeded | `WorldScene3D` places any number of crystals with `Transform3D` placements in one globally depth-sorted scene — e.g. `from_orientation_relationship` renders two phases in a KS/NW/GT/Pitsch/Burgers OR, which VESTA cannot express. |
| Styling system | per-file GUI settings | exceeded | Central YAML theme system (`journal`, `presentation`, `dark`) with per-call overrides; every rendering knob is themed and testable. |

## Verification

- Automated coverage: `tests/unit/test_vesta_parity_features.py`, `tests/unit/test_crystal3d_rendering.py`, `tests/unit/test_scene3d_composition.py`.
- Executable worked examples: `docs/site/examples/generated/visualization.md` (placement identities,
  measured bond length against the exact NaCl-type a/2 geometry).

## References

### Normative

- Momma, K. & Izumi, F. (2011). "VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data." J. Appl. Cryst. 44, 1272-1276.

### Informative

- [Visualization Style Guide](../standards/visualization_style_guide.md)
- [MTEX Parity Matrix](mtex_parity_matrix.md)
