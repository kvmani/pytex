# Plotting Validation Matrix

PyTex treats plotting as part of the stable scientific surface, not as an optional decoration layer.
This matrix records what is currently validated for the runtime plotting subsystem.

Status terms:

- `implemented`: automated coverage and documentation exist for the stable plotted surface
- `foundational`: the implementation exists and is scientifically structured, but the external or parity-backed validation surface is not yet complete
- `planned`: the category is accepted but not yet validated adequately

## Current Matrix

| Category | Validation basis | Status | Notes |
| --- | --- | --- | --- |
| Semantic input validation | unit tests over frame, symmetry, and convention mismatches | implemented | Plot builders reject invalid frame or symmetry combinations before rendering. |
| Runtime figure generation | Matplotlib figure creation and export tests | implemented | Public plotters return ordinary Matplotlib figures and can be exported to raster or vector formats by users. |
| Documentation SVG policy | documentation policy tests and `save_documentation_figure_svg(...)` coverage | implemented | Repo-tracked docs figures remain SVG while runtime plotting stays backend-native. |
| Pole-figure contour rendering | builder and runtime tests for contour layers | implemented | Contour pole figures are rendered from smoothed projected density grids with deterministic builder behavior. |
| ODF contour rendering | builder and runtime tests for Euler-space contour layers | implemented | The current contour path is a discrete support inspection surface in Euler space. |
| ODF Bunge-section rendering | builder and runtime tests for multi-panel section figures | implemented | Section plots now support both the discrete ODF inspection path and the harmonic ODF evaluation path on explicit Bunge-section grids. |
| YAML theme loading and merging | unit tests over theme catalogs, overrides, and invalid payloads | implemented | Runtime plotting style is centralized rather than distributed across plotting routines. |
| Powder XRD plotting | runtime plotting tests and style-config coverage | implemented | XRD figures are deterministic Matplotlib surfaces built on shared style resolution. |
| SAED plotting | runtime plotting tests and style-config coverage | implemented | Detector-space spot plots are validated for figure creation and semantic labeling. |
| 3D crystal visualization | scene-builder tests and runtime figure generation | implemented | Unit-cell atoms, heuristic bonds, repeated cell overlays, hexagonal-prism overlays, bounded plane overlays, direction overlays, and scientific Miller annotations are covered by current tests. |
| Stereographic direction and plane plotting | geometry tests, runtime plotting tests, and structural figure-property assertions | implemented | Wulff-net overlays, projected crystal directions, plane poles, and great-circle traces are validated without repo-tracked SVG baselines. |
| Rotational symmetry-element plotting | order-specific symbol mapping tests plus structural figure-property assertions | implemented | Proper rotational axes are rendered with semantic dyad, triangle, square, and hexagon symbols together with shared Miller-style annotations. |
| IPF plotting | workflow docs, notebook examples, and runtime figure-structure checks | implemented | Stable plotting coverage checks the structural output surface rather than byte-stable SVG artifacts. |
| Visual parity against external tools | MTEX-style rendered comparisons and pinned image baselines | planned | Useful for future visual regression and semantics checks, but not yet a stable requirement. |
| Publication presets and house styles | deterministic YAML themes plus structural plotting assertions | implemented | Journal, presentation, and dark presets are validated through shared style resolution and figure-surface checks rather than tracked runtime SVG baselines. |

## Interpretation Notes

- Plot validation is not only about pixels. The higher priority is that plotted geometry respects the same frame, symmetry, and reduction contracts as the computational API.
- The current contour and section plots are scientifically meaningful inspection surfaces for the implemented discrete and harmonic texture models, but they are still review plots rather than claims of full external visual parity.
- PyTex does not yet claim visual parity with MTEX or any other plotting package across every rendered detail.

## References

### Normative

- `strategy.md`
- `../standards/documentation_architecture.md`
- `../standards/latex_and_figures.md`

### Informative

- `../site/workflows/plotting_primitives.md`
- `../site/workflows/texture_odf_inversion.md`
- `../site/workflows/xrd_generation.md`
- `../site/workflows/saed_generation.md`
- `../site/workflows/crystal_visualization.md`
- `../site/workflows/stereographic_projections.md`
