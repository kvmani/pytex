# Changelog

All notable changes to PyTex are recorded here, newest first. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions
follow the pre-1.0 semantics of
[API Stability And Deprecation](docs/standards/api_stability_and_deprecation.md)
(minor versions may break with deprecation warnings; patch versions do not).
Every release entry must state scientific behavior changes explicitly —
"fixed" for correctness, "changed" for convention or semantics — because
downstream analyses depend on them.

## [Unreleased]

### Added

- **The same-parent boundary fingerprint is now a public core surface.**
  `intervariant_boundary_fingerprint(relationship)` returns the deduplicated
  set $G_c \left(R\,G_p\,R^{\mathsf{T}}\right) G_c$ of misorientations that two
  child grains of one parent can exhibit, and
  `boundary_fingerprint_distances_deg(relative_matrices, fingerprint)` scores
  measured boundaries against it with a memory-bounded blocked kernel. Both are
  exported from `pytex.core` and the top level.

  This is the quantity that answers "could these two product grains share a
  parent?", and it was previously an undocumented private helper inside
  `pytex.experimental.or_identification`, duplicated in weakened
  (angle-only) form by parent-grain reconstruction. Reconstruction and OR
  identification now share the one definition, per the repository's
  one-shared-helper rule. New worked example
  `or-ks-same-parent-boundary-fingerprint` pins two identities: the Sigma3 twin
  relation (60 deg about $\langle 111 \rangle$) is an admissible
  Kurdjumov-Sachs same-parent boundary — it is Morito's published V1-V20
  intervariant pair — and all 276 variant-pair boundaries of a common parent
  sit at zero distance from the set they generate.

- **Reference frames are now a first-class shared foundation.** Frames were
  previously a thin label plus three axis names, and each module built the ones
  it needed inline. The foundation replaces that with one model used everywhere.
  See [Reference Frame Foundation](docs/architecture/reference_frame_foundation.md).

  - `ReferenceFrame` now carries **axis geometry**: `axis_vectors` (the
    components of its three labelled axes in the canonical right-handed
    Cartesian reference `X, Y, Z`), optional `axis_descriptions` long names, a
    `basis_matrix` property, `axis_index` / `axis_vector` / `unit_axis_matrix`
    accessors, `is_orthonormal` / `is_right_handed` / `determinant` reporting,
    `with_axis_vectors` / `renamed` / `rotated` derivation, and `describe()`.
    Construction now rejects linearly dependent axes and a declared handedness
    that contradicts the axis-vector determinant. The geometry is stored as a
    hashable tuple of float triples, so frames stay comparable — frame equality
    gates `VectorSet`, `FrameTransform`, `Orientation`, and `SymmetrySpec`
    consistency checks.
  - `FrameTransform` gained `from_rotation`, `from_bunge_euler`,
    `from_axis_angle`, `from_axis_correspondence` (state a vendor axis
    convention in words instead of hand-writing a permutation matrix),
    `between_frames`, `as_rotation`, `rotation_angle_deg`, `rotation_axis`,
    `is_identity`, `source_axes_in_target`, and `describe()`. **New:**
    `apply_to_directions` applies the rotation only — directions, plane normals,
    and poles are translation-invariant, so an origin offset must not move them;
    `apply_to_vectors` keeps applying rotation *and* translation for positions.
  - `FrameGraph` registers frames and declared transforms and resolves the
    transform between any two connected frames by composing the **shortest**
    declared chain (fewest matrix products, least accumulated error). Edges are
    usable in both directions.
  - `pytex.core.frame_catalog` builds the standard frames once:
    `CARTESIAN_FRAME`, `SPECIMEN_FRAME`, `SAMPLE_RD_TD_ND_FRAME` (`RD/TD/ND`),
    `CRYSTAL_FRAME`, `MAP_FRAME`, `DETECTOR_FRAME`, `LABORATORY_FRAME`, with
    matching builders, `reciprocal_frame_for`, `rolling_frame_graph`,
    `get_standard_frame`, and `list_standard_frames`. Catalog defaults are
    pinned to the field values the repository's modules already used, so
    adopting the catalog is identity-preserving — asserted directly in
    `tests/unit/test_frame_catalog.py`.
  - `pytex.plotting.frames` renders the same frame three ways from one geometry
    computation (`FrameTriad`): `frame_triad` / `frame_triad_primitives` for 3D
    scenes, `add_frame_indicator` as an **embeddable corner gizmo** for any 2D
    figure (SAED diffractograms, pole figures, IPF maps, crystal-viewer panels;
    works on polar axes), and `reference_frame_svg` / `frame_catalog_svg` as
    standalone documentation SVG generated in pure Python with **no matplotlib
    dependency**. `project_orthographic` and `TRIAD_AXIS_COLORS` are public.
  - Three renderers accept a frame gizmo directly, all **opt-in** so existing
    figures are unchanged: `plot_saed_pattern(show_frame_indicator=True)` shows
    the detector `u`/`v` axes; `render_composite_saed` with
    `CompositeSAEDPlotConfig(show_frame_indicator=True)` shows the *parent
    crystal* axes as they land on the detector; and
    `plot_crystal_structure_3d(show_frame_indicator=True)` shows the phase's
    `a`/`b`/`c` axes from the lattice basis at the figure's own view angles.
  - New generated canonical figures `docs/figures/reference_frame_catalog.svg`
    and `docs/figures/sample_frame_rd_td_nd.svg`, produced by
    `scripts/generate_reference_frame_figures.py` from the same public code path
    users call, so a documentation figure cannot drift from the model.
  - New executable worked examples (`reference_frames` group) checking the
    rotation implied by a declared axis correspondence, the resulting
    components, multi-hop graph composition, exact round-trip invertibility, and
    the right-handed determinant convention.

- **Crystallographic notation is now fixed centrally and enforced.**
  `pytex.core.notation` is the single place PyTex turns crystallographic
  quantities into text, and the conventions it implements are anchored to the
  IUCr *International Tables* in `docs/standards/notation_and_conventions.md`.

  - **The reciprocal star** marks the *basis*, never the indices:
    `format_reciprocal_axis_label(s)` produce `a*, b*, c*` for reciprocal basis
    vectors and reciprocal-frame axes, while Miller indices stay unstarred
    because `(hkl)` are already reciprocal-basis components. Starring is
    idempotent, so a label passing through two layers cannot become `a**`.
    `format_reciprocal_lattice_vector` renders $\mathbf{g}_{hkl}$.
  - **Bracket families** are now expressible: `format_plane_family_indices` and
    `format_direction_family_indices` give $\{hkl\}$ and $\langle uvw \rangle$, alongside the
    existing `(hkl)` and `[uvw]`. `format_miller_indices` gained a `scope`
    parameter.
  - The rule is a non-negotiable in `AGENTS.md` and is **enforced** by
    `tests/unit/test_notation_conventions.py`, which fails if any module
    reintroduces inline index formatting or hand-rolled starring, and which
    renders every mathtext form through matplotlib so an unparseable label fails
    as a test rather than as a broken figure.

### Fixed

- **Parent-grain reconstruction linked grains on the misorientation angle
  alone, merging unrelated parents.** This is a scientific behavior change:
  reconstruction groupings and parent counts change, and previously reported
  results on real microstructures should be regenerated.

  `reconstruct_parent_grains` (and therefore
  `reconstruct_parent_grains_from_graph`) decided whether two neighbouring
  child grains descend from a common parent by reducing the intervariant table
  to its **distinct angles** and asking whether the boundary disorientation
  angle fell within `tolerance_deg` of any of them. The misorientation **axis
  was discarded**. For a cubic-cubic relationship those angles are spread
  densely enough over the accessible range that the test was far too
  permissive: against 20 000 uniformly random, entirely unrelated boundaries it
  accepted **52.8%** of them at the default 3 deg tolerance (28.6% at 1 deg,
  62.6% at 5 deg) for Kurdjumov-Sachs, and 39.1% at 3 deg for
  Nishiyama-Wassermann.

  The consequence at map scale was silent merging of distinct parent grains.
  On a 12-parent synthetic microstructure with one contact edge between
  consecutive parents, the angle-only rule linked 5 to 8 of the 11 cross-parent
  boundaries and recovered only **4 to 7 of the 12** planted parents.

  The edge test now matches the **full rotation**, against the admissible
  same-parent set $G_c \left(R\,G_p\,R^{\mathsf{T}}\right) G_c$ — exact, because
  two children of one parent satisfy
  $\mathbf{C}_i^{\mathsf{T}}\mathbf{C}_j = \mathbf{V}_i\mathbf{V}_j^{\mathsf{T}}$
  with $\mathbf{V}_i = R\,S_{p,i}$. False acceptance of unrelated boundaries
  drops to 7.1% at 3 deg and 0.26% at 1 deg, and the same 12-parent fixture now
  recovers 10 to 12 of 12. **No sensitivity is lost:** true same-parent
  boundaries score zero against the fingerprint to 1.2e-6 deg (the
  quaternion/matrix round-trip floor), and no true edge was missed in any
  measurement above.

- **The fingerprint distance kernel allocated gigabytes at map scale.** The
  comparison in `identify_orientation_relationship` was written as
  `einsum("eij,kij->ek", ...)`, which materializes one float per (edge,
  fingerprint element) pair — **4.3 GB for 50 000 edges** against a cubic-cubic
  fingerprint. It is now a blocked `(512, 9) @ (9, k)` GEMM in the shared
  `boundary_fingerprint_distances_deg`: numerically identical to 4e-13, 2.5x
  faster, and bounded at ~22 MB regardless of edge count.

- **spglib 2.7 broke the test suite under the warnings-as-errors policy.** It
  announces its own error-handling migration from inside the library on every
  call through the legacy path, which no caller-side filter or adapter shim can
  suppress at the point of emission. Added as a fourth narrow, commented
  exemption in `pyproject.toml` alongside the existing pymatgen ones.

- **Documentation figures rendered with runaway arrowheads.** SVG markers
  default to `markerUnits="strokeWidth"`, which multiplies the arrowhead by the
  stroke width of the line it terminates, so a figure declaring a 12-unit head
  and drawing a `stroke-width="4"` line rendered a **48-unit** head. Across
  `docs/figures/` this left arrowheads occupying 11% to 125% of the lines they
  annotated — in the worst case the head was longer than the whole arrow, and
  the reference-frame triads were unreadable.

  Six frame and orientation-convention figures are now **generated from the
  model** by `scripts/generate_reference_frame_figures.py`:
  `reference_frames.svg` (the canonical chain), `reference_frames_vectors.svg`
  and `orientation_mapping_semantics.svg` (the crystal-to-specimen mapping, the
  second showing the inverse as a separate relationship),
  `active_passive_rotation.svg`, `bunge_euler_geometry.svg` (one computed panel
  per Euler step), and `hcp_reference_frame.svg` (basal axes read from
  `Lattice.direct_basis()`). Their axis directions are therefore the modelled
  axis directions, and their layout is computed so text cannot overflow or
  collide.

  The remaining 30 hand-authored figures were corrected in place by the new
  `scripts/fix_svg_marker_units.py`, which switches each marker to absolute
  units and pre-scales its geometry to preserve the figure's intended visual
  weight while bounding the head against the lines it terminates. Median
  head-to-line ratios dropped from as high as 1.11 to at most 0.25.
  `tests/unit/test_figure_markers.py` fails if the defect reappears.

  While regenerating, one scientific error in the old chain figure was also
  fixed: it drew the reciprocal frame as a link in the linear chain, implying a
  `laboratory -> reciprocal` step. Duality relates the reciprocal frame to the
  **crystal** frame, so it is now drawn off the chain, matching the canonical
  frame chain in the notation standard.

- **Index formatting was ambiguous for negative and multi-digit components.**
  `format_miller_indices` concatenated components unconditionally, so `[1-10]`
  could be read as `[1, -1, 0]` *or* `[1, -10]`, and `(1210)` as `(1, 2, 1, 0)`
  or `(12, 1, 0)`. A separator is now inserted whenever a component is negative
  in plain style or any component has more than one digit; single-digit
  non-negative indices keep the classical concatenated form `(110)`. This
  changes user-visible label text — `describe()` output and figure labels for
  such indices now read `[1 -1 0]` where they previously read `[1-10]`.

### Changed

- **Pole figures and powder reflections are labelled as families.** A pole
  figure plots the whole symmetry-related orbit of its pole, and a powder
  reflection *is* its multiplicity, so both now read $\{hkl\}$ rather than
  `(hkl)`; writing a single member misstated the quantity. Because a
  `PoleFigure` can be built with `include_symmetry_family=False`, the object now
  records `includes_symmetry_family` and titles follow the record rather than an
  assumption. JSON contracts round-trip the new field, defaulting to `True` for
  payloads written before it existed.
- Five modules that formatted indices inline (`plotting/diffraction.py`,
  `plotting/builders.py`, `plotting/composite_saed.py`,
  `diffraction/composite.py`, `diffraction/kinematic.py`) now route through
  `pytex.core.notation`.
- Composite SAED reciprocal-space axis labels now write the scattering vector in
  bold per IUCr vector convention.
- Fixed five long-standing Sphinx cross-reference warnings in
  `docs/standards/reference_canon.md`; the documentation build is now
  warning-free.

- **Every module now builds frames through the shared catalog.**
  `adapters/scan_files.default_ebsd_frames`, `diffraction/saed`,
  `core/lattice.Lattice.reciprocal_basis`, the CLI core demo, and the plotting
  validation cases no longer construct `ReferenceFrame` inline. This is
  behaviour-preserving: the catalog defaults reproduce the previous field values
  exactly, so the frames compare equal and no downstream consistency check
  changes. The one visible difference is that a reciprocal frame's axis labels
  are now starred (`a, b, c` becomes `a*, b*, c*`), which makes a
  reciprocal-space vector impossible to mistake for a direct-space one.
- `plotting.primitives.reference_frame_triad` now honours a frame's own
  `axis_vectors` instead of always drawing the canonical Cartesian triad, so a
  frame recorded as rotated draws rotated. An explicit `basis` argument still
  wins.
- `pytex.contracts` serializes `axis_vectors` and `axis_descriptions` for
  reference frames. Deserialization is backward compatible: payloads written
  before these fields existed get the identity triad and no long names, which
  reproduces exactly the frame they described, so older files still round-trip
  to equal objects.

- **Symmetry-reduced disorientation is now a single dense product.** The
  reduction $\min_{S_l, S_r} \angle\!\left(S_l \mathbf{M} S_r^{\mathsf{T}}\right)$ previously expanded an
  `(n, |S_l|, |S_r|, 3, 3)` candidate array through a chain of einsums. Because
  the disorientation angle depends only on the *scalar* part of the
  symmetry-conjugated relative quaternion, and that scalar part is linear in
  the quaternion, the whole reduction collapses to one precomputed matrix of
  linear functionals and a single `(n, 4) @ (4, k)` product per memory-bounded
  block. Rows that agree up to sign are redundant under `|.|`, so for
  same-phase cubic symmetry the 24 x 24 operator pairs deduplicate to 24
  functionals. Results are unchanged to 5e-14 rad across twelve point groups
  and four cross-symmetry pairs, pinned by
  `test_reduced_disorientation_kernel_matches_trace_reference`.

  Measured on this machine (fcc nickel scan, `m-3m`, 4-connected):

  | operation | before | after | speedup |
  | --- | --- | --- | --- |
  | reduction kernel, 200 000 pairs | 2.504 s | 0.110 s | 22.8x |
  | KAM, 13 000 points | 0.625 s | 0.012 s | 52.4x |
  | KAM, 61 600 points | 6.491 s | 0.142 s | 45.8x |
  | `segment_grains`, 2 080 points | 5.73 s | 0.06 s | 95x |
  | `segment_grains`, 13 000 points | 114.4 s | 1.02 s | 112x |
  | `or_deviation`, 5 000 pairs | 1.591 s | 0.137 s | 11.6x |
  | `reconstruct_parent_grains`, 400 grains | 3.946 s | 0.428 s | 9.2x |

- `CrystalMap.segment_grains` no longer allocates an `(n, n)` angle matrix per
  grain, and `pytex.ebsd.models` no longer expands an unbounded
  `(pairs, |S|, |S|, 3, 3)` candidate array for neighbour misorientations —
  that array reached 1.07 GB for a 13 000-point scan and 5.09 GB for a
  61 600-point one. Both paths now accumulate through the shared,
  block-bounded kernel, so memory is flat in the number of pairs.

- **Grain representative orientations now resolve exact ties deterministically.**
  A grain whose members are symmetric about its centre has no unique medoid;
  under bare `argmin` the representative — and therefore the reference
  orientation GROD is measured from — could be decided by summation order, the
  BLAS build or the machine. Members within a relative `1e-9` of the minimum
  total disorientation are now treated as tied and the lowest index wins. On
  the reference scans this changes the representative for grains that had two
  candidates agreeing to ~1e-10 rad, with the next distinct candidate 1e-5 to
  3e-4 rad away; GROD maps are unchanged.

### Fixed

- **Small-angle misorientation accuracy (scientific).** Neighbour
  misorientations were computed as $\arccos\!\left((\operatorname{tr} - 1)/2\right)$ on a triple matrix
  product, which is ill-conditioned exactly where EBSD measures — KAM, GROD and
  low-angle boundaries all live below 1 degree. Against a well-conditioned
  `atan2` reference the old path erred by up to 3.5e-8 rad; the quaternion path
  errs by 4.5e-13 rad, roughly five orders of magnitude better. Reported KAM
  and GROD values shift in the eighth decimal.

- Notebook 07 (`07_ebsd_regular_grid_workflows`) built a `CrystalMap` whose
  orientations lived in the specimen frame while its grid lived in the map
  frame, then called `to_experiment_manifest()`, which raises because the
  specimen-to-map relationship is undefined. The notebook now supplies an
  explicit `AcquisitionGeometry` with a `specimen_to_map` `FrameTransform`.
  The defect was latent because the notebook had never been executed.

### Added

- Tutorial notebooks are now hand-authored `.ipynb` files edited directly.
  `scripts/generate_tutorial_notebooks.py` has been **removed**: it constrained
  how notebooks could be written and rewrote every notebook with empty outputs
  on each run. Removing it exposed that notebooks 01-17 had never been executed
  and were publishing as bare code listings; all 21 notebooks are now committed
  executed (59 images render across the site, up from 13). Two guard tests
  (`test_every_notebook_is_committed_executed`,
  `test_no_notebook_contains_error_output`) now enforce this.
- Burgers beta->alpha (bcc -> hcp) is now a canonical case alongside
  Kurdjumov-Sachs across the composite-diffraction tests, examples and
  documentation. Hexagonal phases are labelled in four-index Miller-Bravais
  notation throughout (`is_hexagonal_phase`,
  `RationalizedZoneAxis.indices_bravais`, bravais flags on `SpotCoincidence`,
  `format_hkl(..., bravais=True)`); cubic phases keep three-index labels.
- Executed tutorial notebook
  `21_composite_or_diffraction_patterns` documenting the composite OR
  diffraction surface with Ewald/excitation-error theory, original diagrams,
  and both canonical cases (KS and Burgers).
- Composite orientation-relationship SAED simulation and rendering
  (kinematic only). New `pytex.diffraction.kinematic` provides a fully
  vectorized zone-axis engine (`simulate_zone_axis_spots`, `SpotTable`,
  `KinematicSimulationConfig`, `zone_basis_from_axis`,
  `electron_wavelength_angstrom` — relativistic, pinned to standard
  values; excitation-error reflection selection handling irrational zones;
  vectorized centering absences and electron structure factors). New
  `pytex.diffraction.composite` assembles a parent phase plus any subset of
  OR variants on one shared parent-anchored detector for an arbitrary parent
  zone axis (`simulate_composite_saed`, `CompositeSAEDPattern`,
  `VariantZonePattern`, `rationalize_zone_axis` for nearest-rational child
  zone labels), and quantifies which reflections superimpose
  (`find_spot_coincidences`, `SpotCoincidenceReport`) plus a
  `sweep_parent_zone_axes` survey iterator. New
  `pytex.plotting.composite_saed` renders it with a typed, publication-grade
  configuration (`render_composite_saed`, `CompositeSAEDPlotConfig`,
  `SpotStyle`, `SpotAnnotationConfig`): per-variant marker/color/size styling,
  variant subsetting, in-plane rotation, mm / $\text{\AA}^{-1}$ axes, and coincidence-merging,
  crowding-aware spot annotation. Report objects carry `describe()`; two
  worked examples and a workflow page document the surface. See
  `docs/roadmap/working_notes_composite_saed_program.md`.
- Orientation-relationship analysis flagship (development-guide Cycles A-B and
  follow-ons): index correspondence with rationalization and angular residuals
  (`correspondence_direct`/`correspondence_reciprocal`,
  `map_plane_to_child`/`map_direction_to_child` and inverses, across-variant
  tables); the misorientation representation (`misorientation()`) and
  deviation metric (`or_deviation`); parallelism finders
  (`find_parallel_planes`/`find_parallel_directions`); OR fitting
  (`fit_orientation_relationship`); variant packet classification
  (`variant_close_packed_groups`); variant pole figures
  (`variant_pole_figure`, `plot_variant_pole_figure`); named KS, GT, Pitsch,
  Burgers, Shoji-Nishiyama, Pitsch-Schrader, Potter, Bagaryatsky, and
  Isaichev constructors with standard catalogs; intervariant
  misorientation tables.
- Experimental OR identification from child-child boundaries
  (`pytex.experimental.identify_orientation_relationship`): ranks candidate
  relationships by their double-coset intervariant fingerprint, no parent
  orientations required.
- Orientation-relationship documentation program: executed tutorial notebooks
  18-20 (fundamentals; lattice correspondence and transformation strain;
  catalogs, identification, and reconstruction) with equations, rendered
  figures, and five reusable scientific SVG diagrams under
  `docs/site/_static/or/`; `scripts/execute_notebooks.py` executes notebooks
  in place so the site renders their outputs.
- Experimental boundary-based OR rotation refinement
  (`pytex.experimental.refine_orientation_relationship_from_boundaries`):
  recovers the operative rotation from child-child boundary misorientations
  alone by alternating coset-element assignment with least-squares updates.
- Experimental map-scale parent-grain reconstruction
  (`pytex.experimental.reconstruct_parent_grains`, `..._from_graph`) with
  intervariant-fingerprint edge testing, union-find clustering,
  quaternion-averaged parent refinement, and EBSD grain-graph wiring.
- Explainable-results doctrine: `describe()` prose on every stable
  transformation report, substring-validated in tests.
- Transformation deformation gradients (`deformation_gradient()`,
  `DeformationGradientReport`): nearest-integer lattice correspondence, polar
  decomposition, textbook Bain stretches and the literature KS/NW rigid
  rotations pinned.
- Texture kernel breadth: `GaussianSO3Kernel` (Gauss-Weierstrass spectrum)
  and `AbelPoissonKernel` beside de la Vallee Poussin, with closed-form
  Chebyshev coefficients and halfwidth-defined construction;
  `KernelSpec.as_so3_kernel()` routes all three.
- Engineering: warnings-as-errors test policy with zero-warning suite;
  coverage ratchet (87%) and ubuntu+macos x Python 3.11-3.13 CI matrix;
  Hypothesis property suites (rotation algebra, Miller-Bravais round trips,
  correspondence invariants, hexagonal metrics); runnable transformation
  performance benchmark lane; `OrientationSet` slicing;
  `CrystalDirection.from_cartesian`; public `phases_semantically_match`.

### Fixed

- **Orientation-convention bug (scientific):** the transformation stack
  composed predicted children as `V @ P`, contradicting the normative
  crystal-to-specimen orientation convention; all prediction, deviation,
  fitting, scoring, and reconstruction surfaces now compose
  $g_{\text{child}} = g_{\text{parent}} \circ \mathbf{V}^{\mathsf{T}}$, pinned by a specimen-space parallelism
  regression test. Synthetic data was internally consistent either way; real
  measured orientations would have received wrong variant assignments.
- `SymmetrySpec` equality raised `ValueError` on distinct-but-equal instances
  and the class was unhashable; equality is now explicit semantic identity.
- Matplotlib figure leaks in the test suite; spglib's force-enabled
  dict-interface deprecation warning is intercepted at the adapter boundary.

### Changed

- `OrientationRelationship.parallel_directions` stores typed
  `CrystalDirection` pairs (index meaning preserved); the JSON contract emits
  typed payloads and still reads legacy float triples.

## [0.1.0.dev0]

Baseline development snapshot predating this changelog: canonical core model
(frames, symmetry, lattice, orientations, batches), texture (PF/IPF/ODF,
harmonics), EBSD (crystal maps, grains, KAM/GROD, CSL), diffraction (powder
XRD, SAED, scattering factors), plotting (pole figures, IPF maps, ODF
sections, VESTA-class crystal viewer, visualization primitives), adapters,
manifests, and the executable worked-example documentation system.
