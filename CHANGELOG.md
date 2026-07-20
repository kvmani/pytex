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
  variant subsetting, in-plane rotation, mm/Å⁻¹ axes, and coincidence-merging,
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
  `g_child = g_parent o V^T`, pinned by a specimen-space parallelism
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
