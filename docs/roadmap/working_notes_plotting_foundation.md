# Working Notes: Publication Plotting Foundation (2026-07-13)

Purpose: resumable ledger for the plotting-focused development phase. Goal:
foundational, reusable plotting classes/modules that make every PyTex figure
publication-ready and rich enough to exceed MTEX (texture) and OVITO (atomic
visualization) output quality. If the session is interrupted, resume by
reading this file top to bottom.

Gates per task (unchanged): `python -m ruff check .`, `python -m mypy src`,
`python -m pytest -q --deselect tests/unit/test_phase_fixtures.py --deselect
tests/unit/test_repo_integrity.py`. Commit + push per task.

## Foundation Ledger

| # | Task | Status | Commit |
| --- | --- | --- | --- |
| A | `plotting/colormaps.py`: scientific colormap + palette foundation | done, tests green | 6a8eada |
| B | `plotting/figure.py`: rc bridge, `PanelGrid`, panel labels, scale bar, multi-format export | done, tests green | 4d4cfef |
| C | Wire foundations into ODF sections + EBSD maps | done, tests green | this commit |

### What the foundation provides (reuse surface)

- `pytex.plotting.colormaps`: `ColormapSpec`, `register_pytex_colormaps()`
  (idempotent), `get_pytex_colormap`, `categorical_colors` (fixed-order
  Okabe-Ito), `srgb_to_lightness` (CIELAB L* check used by tests).
  Registered maps: `pytex.texture` (white-anchored monotone-lightness m.r.d.
  ramp), `pytex.misorientation` (white-to-deep-blue single hue),
  `pytex.diverging` (blue/red, neutral midpoint).
- `pytex.plotting.figure`: `rc_params_from_style` (theme -> full rcParams,
  including `svg.fonttype: none` for editable SVG text and TrueType PDF/PS
  fonts), `publication_style` context manager, `PanelGrid` (uniform panel
  sizing, `label()` for (a),(b),(c), `shared_colorbar`, `hide_unused`,
  `export`), `label_panels`, `add_scale_bar` (anchored micrograph bar),
  `export_figure` (SVG/PDF/PNG, stem- or suffix-addressed).
- Consumers wired so far: `plot_odf_phi2_sections` (PanelGrid layout,
  Bunge-Euler math labels, m.r.d. colorbar, `pytex.texture` default),
  `plot_phase_map` (Okabe-Ito identity palette), `plot_kam_map`
  (`pytex.misorientation` default), and `scale_bar=` on IPF/KAM/property/
  phase maps.

## Next Candidates (rough priority)

1. Pole-figure / IPF publication upgrade: route `builders.py` contour paths
   through `pytex.texture`, add m.r.d. colorbars, RD/TD specimen-axis
   annotations on the projection rim, and named-component markers.
2. ODF sigma sections (`ODF.sigma_sections` + plotter) to pair with phi2
   sections (roadmap 4.1 "section plotting parity").
3. IPF color-key legend panel: render the `IPFColorKey` sector triangle as an
   inset/companion axis for IPF maps (MTEX-style key-next-to-map).
4. Crystal3d OVITO-grade pass: theme-driven specular/ambient defaults exist;
   add depth-sorted bond/atom composition checks, orthographic scale bar, and
   `export_figure` wiring for ray-clean SVG/PNG.
5. Map annotation kit: north/RD-TD arrows, grain-boundary legend entries,
   inset zoom boxes (reusing `PanelGrid`).
6. Gallery documentation page generating all foundation figures via
   `export_figure` into `docs/figures/` with integrity checks.

## Known Constraints

- Matplotlib is an optional extra: every foundation module keeps imports
  inside functions and raises the standard install hint.
- Existing tests pin plot titles ("Phase Map", "Kernel Average
  Misorientation", "IPF Map (z)") and the phi2 3-sections + colorbar axes
  count; upgrades must keep those contracts or update tests deliberately.
- The structural validation matrix (`_plotting_validation_cases.py`) is the
  place to register new canonical figures.
