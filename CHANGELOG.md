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

## [0.2.0] - 2026-08-23

This production release adds the configurable feedback and onboarding experience and makes
Miller-index input unambiguous throughout the workbench. It changes no crystallographic
convention or numerical algorithm. The compatibility action for service callers is to replace
an unseparated value such as `"111"` with `"1 1 1"` or `[1, 1, 1]`.

Release packaging now uses the SPDX license expression supported by current setuptools, removing
the deprecated license table/classifier metadata without changing PyTex's GPL-3.0-or-later terms.

The exhaustive autodoc reference remains fully rendered without duplicating the curated API
guide's object-index entries, restoring a useful Sphinx warning signal for release builds.

### Added

- **A feedback and feature-request form, in both shells.** The masthead carries a **Feedback**
  button that opens a short form: what this is about, the note itself, and — all optional — who
  sent it. The page attaches which workspace and panel were open. Every submission is appended to
  a JSON file *before* anything is sent, so a relay outage costs a notification rather than a
  note, and the receipt on screen distinguishes "filed here" from "filed and e-mailed". A
  deployment that configures an internal SMTP relay also mails it to the maintainer. New modules
  `pytex.app.feedback` and routes `GET /api/experience`, `POST /api/feedback`.
- **A deployment configuration file.** `pytex.app.config` reads three sections — `feedback`,
  `relay` and `tour` — from `PYTEX_APP_CONFIG`, `./pytex_app.yml` or `~/.pytex/pytex_app.yml`, and
  works with no file at all. `config/pytex_app.example.yml` is the annotated template, loaded by a
  test through the same loader so a renamed key cannot stay documented under its old name. Nothing
  in it changes a number PyTex computes; an unknown key is an error at startup; a relay password
  is named as an environment variable rather than written in the file.
- **A welcome message and a skippable tour**, on by default and switchable from the same file. Seven
  steps over the workspaces, the generated controls, the feature search, the message log and the
  feedback form, each highlighting the element it is about. Skip is on every step, skipping is
  remembered per browser and per version, and the tour reopens from the Help panel.

### Changed

- **Miller indices are entered one index per box, and a run of digits is now refused.**
  Users were typing `110` for what they meant as `1 1 0`, and the service layer guessed it
  correctly — which was the problem. The guess only works while every index is a single digit:
  `10 10 0` typed as `10100` is `(10, 10, 0)`, `(1, 0, 100)` and `(101, 0, 0)` with nothing to
  choose between them, and the wrong reading returns a plausible number rather than an error. The
  workbench now renders every `indices` and `indices-list` parameter as one small box per index,
  named `h`, `k`, `l` — or `u`, `v`, `w` where the parameter's label says so — so there is no
  place to put a run of digits; a list parameter is a stack of those rows with add and remove
  buttons. `IndicesParameter` and `IndicesListParameter` reject an unseparated run whatever its
  length, naming the separated form of what was typed. **API note:** a caller that was passing
  `"111"` to the service layer must pass `"1 1 1"` or `[1, 1, 1]`.

### Fixed

- **Only the first index box of a row was reachable.** The boxes were styled by class, and the
  stylesheet's base rule for a text input is `width: 100%` on an attribute selector, which
  outweighs a bare class: every box stretched to the full width of the control rail, so the first
  one filled the panel and the second and third were pushed off the edge of the screen with no way
  to reach or see them. A zone axis could therefore only be given its `u`. The width now carries
  the same attribute selector and holds four monospaced characters — three digits and a sign,
  which covers every index anyone types — so all three boxes sit inside the column the fields
  above and below them occupy. A browser test measures the boxes against the rail so the rule
  cannot be outweighed again unnoticed.

## [0.1.1] - 2026-08-21

The first cut release. `0.1.0.dev0` was a development snapshot with no tag; this section is
everything that has landed since, cut into a version so that a result can name the code that
produced it.

### Release Notes

**What this release is.** PyTex became an application as well as a library. The workbench went
from nine flat tabs to seven workspaces holding sixteen panels, and the transmission-electron and
EBSD surfaces are each one subject seen several ways rather than several subjects. Every panel
ships runnable examples that the test suite executes, so nothing on offer can rot quietly.

**The scientific additions.** Kikuchi geometry is now a first-class surface in both settings it
occurs in: the crystal's whole band network on a stereographic map, the bands of a solved TEM
plate superimposed on the plate itself, and — new here — the pattern an EBSD camera would record,
configured in the terms a microscope is configured in. Alongside them, SAED simulation, dynamical
CBED with diffraction-group determination, the Kearns parameter, pole-figure arithmetic on the
m.r.d. scale, EBSD scan readers for `.ang`, `.ctf` and OIM HDF5, and the orientation-relationship
programme's variant, correspondence and reconstruction machinery.

**What to check if you are upgrading.** One output changed shape rather than value: zone axes
reported by `simulate_kikuchi_pattern` are now reduced to coprime indices, because a zone axis is a
direction and `[002]` is `[001]`. Code that counted the axes of a pattern will get a smaller and
more meaningful number; code that looked one up by a multiple of its indices will not find it.
Everything else in this release is additive.

**Where to start.** `docs/site/tutorials/` — thirty-one executable notebooks, each of which now
carries a short "Good to know" note — and `docs/site/workflows/workbench_application.md` for the
application.

### Added

- **`ebsd.simulate_kikuchi_pattern`, and a *Kikuchi simulator* in the EBSD
  workspace.** The forward problem behind every indexed point of a scan:
  a phase, an orientation and a camera in, the pattern that camera would record
  out. It is configured the way the microscope is — stage tilt, camera elevation
  and azimuth, pattern centre, camera distance as the fraction of the screen
  width a calibration reports — and it takes no scan, because it is the problem
  the other six views of that workspace live downstream of. Bands are drawn as
  the gap between their two Kossel-cone edges with the plane trace dashed
  between them, named along themselves, with the zone axes and the pattern
  centre marked.

- **`DiffractionGeometry.for_ebsd(...)`** turns those EBSD terms into the one
  canonical diffraction geometry, and owns the frame convention so that no
  caller composes the stage rotation and the camera placement by hand. The beam
  is the laboratory z axis, the stage tilt axis is x, and the camera sits on -y;
  the specimen normal then projects at gnomonic radius
  `tan(90° - tilt + elevation)`, which is the arithmetic that checks the whole
  convention and is what the new worked example
  `diffraction-ebsd-specimen-normal-radius` computes.

- **The workbench is organised into workspaces, and every transmission-electron
  and EBSD surface is grouped under one.** Nine flat tabs became seven
  workspaces holding fifteen panels. *TEM Analysis* carries **SAED Simulator**,
  **TEM Solver**, **CBED** and **Composite SAED**; *EBSD* carries **IPF map**,
  **GROD**, **KAM**, **Scan summary**, **Distributions** and **Pole figures**.
  A sub-tab strip below the masthead names the view; the tab bar names the
  subject. The five other workspaces are single panels and are unchanged. An
  EBSD scan opened in any of that workspace's views is the scan every view
  analyses, because it belongs to the session rather than to a panel.

- **`tem.simulate_saed` — the forward diffraction problem.** A phase and a
  direction in, the reflections that land on the detector out, with the
  calibration attached so a simulated plate can be carried to the solver without
  a number being retyped. The crystal is pointed either by zone axis and roll or
  by Bunge Euler angles; in the second mode the operation resolves the
  orientation to the nearest low-index zone axis and **states the deviation**,
  because a spot pattern exists only on a zone axis and an orientation between
  two of them has no pattern of its own to draw. Its Kikuchi overlay is
  requested from the existing `tem.kikuchi_overlay` with the matrix the
  simulation reports, so the bands and the spots cannot disagree.

- **`SyntheticSAEDImage.crystal_to_pattern()`** states the orientation a
  simulated plate was built from, in the same convention the indexer reports and
  the Kikuchi overlay consumes: `R @ g` gives a reflection's detector
  coordinates and its excitation error in the pattern frame.

- **`crystal.kikuchi_map` — the crystal's band network about a chosen zone
  axis**, drawn as a third figure in the Crystal Viewer's orientation dock. It
  does not turn with the camera: the pole figure and the inverse pole figure are
  the same view of the same crystal as the structure, while a map is fixed to
  the crystal and what moves on it is the marker showing which direction the
  current view has on the beam. The centre defaults to `[001]` and is editable.

- **`pytex.diffraction.projected_trace_runs`**, and
  **`KikuchiMapBand.centre_directions` / `.edge_directions`.** A one-hemisphere
  projection folds the far half onto the near one, which is right for a point
  and wrong for a curve: a band edge straddling the equator gains a spurious
  chord across the map. The splitting and the unprojected trace samples a
  renderer needs were private to the plotting layer; they are public now,
  because the browser draws the same map and must break the same curves in the
  same places. The plotting layer delegates to both.

- **A pole-figure fly-by** in the crystal viewer, beside the inverse pole
  figure's, drawn as a fading cloud because a pole family turns rigidly and
  every member of it moves at once.

- **`crystal.orientation` reports Bunge angles unconditionally** (`euler_bunge`),
  whichever convention was requested, and the dock's readout leads with them. It
  is the convention every EBSD file, every ODF section and every published
  orientation is written in.

- **`ebsd.scan_summary`, `ebsd.distribution` and `ebsd.discrete_figure`.** A
  scan summary in four sections — acquisition, indexing quality, phases,
  microstructure — with every microstructural number quoted against the grain
  threshold that produced it and every quality channel given its median and 5th
  and 95th percentiles rather than a mean that hides whether a scan is uniformly
  fair or half excellent and half unindexed. A histogram of grain size, area,
  aspect ratio, boundary misorientation angle, KAM, GROD or a measured channel;
  the misorientation distribution carries a random reference computed from
  randomly paired points **of the same scan**, so a textured material is not
  accused of a boundary preference it does not have. And the discrete pole and
  inverse pole figures — the scatter a contour is an estimate of — with the
  subsampling stated and seeded.

- **`fixtures/tem/zr_hcp_basal_saed.png`**, a 17 kB simulated zirconium `[0001]`
  plate with a JSON sidecar carrying its answer, written by
  `scripts/generate_tem_test_pattern.py`. It exists so that opening a pattern
  from disk can be tested the way a user does it; the browser test opens it,
  picks the three seed reflections and requires zirconium down `[0001]` with
  every spot indexed.

- **EBSD scans read from EDAX OIM HDF5 files (`.oh5` and `.h5`).**
  `pytex.adapters.read_oh5` reads an OIM Analysis HDF5 scan into the same
  `CrystalMap` and `EBSDImportManifest` that `read_ang` and `read_ctf` produce.
  `.oh5` and `.h5` are one container under two extensions, and one reader
  serves both. Checked against a real OIM 8.6 scan and the `.ang` OIM exported
  from it: the orientation matrices agree to 8e-6, the coordinates to 4e-7, and
  the confidence-index and fit channels to 5e-4 — the residual being the five
  decimal places the text export rounds to. Both EDAX formats resolve symmetry
  through the same TSL code table (`LGsymID` in HDF5, `# Symmetry` in text),
  carry Bunge angles in radians, and keep or drop the same points, so a scan
  imports identically whichever way it was saved. The HDF5 export carries more:
  every per-point scalar channel in the file is read, not only the columns an
  `.ang` row has room for, with the four channels `.ang` shares keeping
  `read_ang`'s names (`image_quality`, `confidence_index`, `detector_signal`,
  `fit`). `h5py` is an optional dependency (extra `hdf5`), imported lazily.

- **The workbench opens HDF5 scans too.** *Open a scan* in the EBSD panel now
  accepts `.oh5` and `.h5` beside `.ang` and `.ctf`, and dispatches through
  `read_scan`, so a format the readers gain is a format the panel accepts.
  Binary files ride the same JSON request as text ones, base64-encoded in the
  same field (`pytex.app.uploads.uploaded_bytes`), rather than through a second
  transport. A scan file that cannot be read now says so under the file button
  and in a toast: its control is hidden, so the error had nowhere to land and
  the panel failed silently while going on drawing the practice dataset.

- **`pytex.adapters.read_scan` picks the reader from the path.** Extension
  dispatch across `.ang`, `.ctf`, `.oh5` and `.h5` is written once, with
  `SCAN_FILE_SUFFIXES` naming what is accepted and `scan_reader_for` exposing
  the choice on its own.

### Changed

- **A Kikuchi band's indices are written along the band, wherever bands are
  drawn**, including `plot_kikuchi_pattern(..., label_bands=True)` in the
  plotting layer, which could not name its bands at all before. A band is identified by which line it is, so a horizontal caption
  beside a steeply running band belongs — visually — to whichever line is
  nearest the text; on a zone-axis plate a dozen bands cross within a few tens
  of pixels. The convention now lives in one place and is shared by the
  simulated-plate overlay, the crystal's stereographic map, and the new EBSD
  simulator.

- **The Crystal Viewer's Kikuchi map is magnified with the wheel**, moved with a
  drag and fitted again with a double-click, and it names its bands. At the size
  a dock figure occupies, a band that is a fraction of a degree wide on a map
  spanning sixty is a picture of a network rather than something to read indices
  off; magnification is what turns it back into an atlas, and the number of
  named bands rises with it.

- **Zone axes reported by `simulate_kikuchi_pattern` are reduced to coprime
  indices.** A zone axis is a *direction*: `[002]` is the same axis as `[001]`,
  satisfies the same zone law and projects to the same point, so listing both
  counted one hub twice and made the number of axes on a pattern meaningless.

- **The TEM solver's measurements moved out from under the pattern.** The
  measured-pick table and the live cursor readout were drawn inside the drawing
  area, which hid the measurement behind the very spots it measures on the one
  panel where the figure is worked on rather than looked at. `plotFrame` gained
  a **readout bar** under the drawing — capped at a quarter of the card and
  scrolled internally, so a long pick list cannot starve the figure — and the
  measured table gained a **|g| column** beside R and d. The pattern and the
  stereogram now stack below 86 rem rather than 72 rem, because at 1280 px the
  split wrapped the plate's eight tools into a 172 px header over a 169 px
  drawing.

- **The `diffraction` panel is titled *Composite SAED***, which is what it
  draws. Operation ids are unchanged.

- **EBSD maps draw at interactive speed.** The `ebsd.map` panel operation is 8x to 106x faster
  with an identical result — the base64 image, the IPF colour checksum, the grain and
  boundary-segment counts and every grain's reference-orientation index match the previous
  implementation value for value on all three practice datasets. On a 200-point-a-side map
  (40 000 measurements): `sigma3_twin` 66.9 s → 0.63 s, `equiaxed_polycrystal` 39.7 s → 1.98 s,
  `bicrystal_gradient` 143.1 s → 3.51 s. Three kernels changed. The grain reference orientation,
  which was 95% of the wall time, no longer searches the symmetry group on every pair: the
  rotation angle is a bi-invariant metric, so an angle taken on a single branch and found below
  half the group's smallest non-identity rotation *is* the disorientation, and a grain whose
  members all lie within a quarter of that angle of their mean therefore reduces to one dense Gram
  matrix — with a fallback to the full group search whenever that certificate does not hold.
  `segment_grains` now floods the neighbour graph with `scipy.sparse.csgraph.connected_components`
  instead of a per-point Python union-find, and splits grain members with one sort rather than one
  full scan per grain. The symmetry-orbit reduction behind IPF colouring forms the orbit as a
  single dense product and no longer renormalizes it member by member. Derived in
  [EBSD Grain Segmentation And GROD Foundations](docs/site/theory/ebsd_grain_segmentation_and_grod.md).

- **A map's grain boundaries are one artist, not one per face.** `plot_ipf_map`, `plot_kam_map`
  and the other `pytex.plotting.ebsd` maps draw `boundary_overlay` as a single `LineCollection`
  rather than an `ax.plot` call per boundary face. On a 20 000-face network the figure builds in
  0.379 s instead of 7.057 s and redraws in 0.166 s instead of 1.730 s. The drawing is unchanged
  to the pixel, cap style included.

- **A figure's own readout no longer covers the figure.** The panel readout pinned to the top-left
  of a plot — the TEM measured-picks table above all — is now painted *under* the drawing, so the
  opaque part of the picture masks it and the data behind that corner is visible again. Bringing
  the pointer onto its rectangle raises it in full; leaving restores the clear view. Clicking is
  unchanged: the card has always been `pointer-events: none` and a click in that corner has always
  reached the figure.

### Fixed

- **A diffraction pattern opened from disk was drawn at the previous pattern's
  zoom.** The TEM solver preserved the camera across every redraw, which is
  right while picking and wrong on an open: at 448 % on a 1024 px practice
  plate, a freshly opened 400 px micrograph arrived as an 89 px crop of itself.
  Nothing errored, so it read as the panel refusing to display the file. The
  camera is now keyed to the pattern rather than to the frame.

- **A distribution of a quantity with no spread crashed the binning.** A
  coherent twin's boundaries are every one of them at 60 degrees, which is a
  distribution — and a very informative one — but numpy refuses a zero-width
  range. It comes back as a single bin holding the whole population. The
  degeneracy test is relative: sixty degrees computed two ways differ in the
  last bit, which is enough to pass a test against zero and then produce
  duplicate bin edges.

- **A uniform grain's reference orientation no longer depends on the machine.** When every member
  of a grain carries the same orientation, every candidate has the same total disorientation to
  the others — zero — and the previous relative tie tolerance had nothing to be relative to, so
  the representative was decided by floating-point residue (`arccos` loses half its significant
  digits at its endpoint, leaving about 2e-8 rad of noise per pair angle by any route). Such a
  cluster is now recognised before the sum is taken, below a millionth of a radian, and answered
  with the lowest member index by definition. Grains with real orientation spread are unaffected:
  the tie band remains purely relative and resolves exactly as before.

### Added

- **Kikuchi bands on the solved TEM pattern.** `tem.kikuchi_overlay` draws the bands the accepted
  solution predicts over the plate it was solved from, behind a **Kikuchi** toggle that appears
  beside *Lattice* and *Calculated* as soon as the pattern is indexed, drawing from the accepted
  solution when there is one and otherwise from the candidate currently selected. A detector records directions of
  outgoing electrons, and the spots and the bands are placed in that one angular space by the same
  reciprocal lattice and the same orientation, so a plane and its normal are drawn as what they
  are — pole and polar of one another. The metrics agree too: **a band is exactly as wide as the
  000→g distance of its own plane and perpendicular to it**, which is a check the user can make by
  eye on the pattern in front of them, and the reason the overlay needs no calibration beyond the
  pixel scale that already indexed the pattern — not the diffraction rotation, not the parity, not
  λ or L separately. Naming a target zone axis draws the connecting band distinctly and labels it
  `follow (200) toward [011]` with the low-index waypoints named, from
  `pytex.tem.path.connecting_band`, which says *no single band connects these zones* rather than
  inventing one. Band edges are the exact Kossel-cone conics, not the straight-line small-angle
  approximation. The result and the status line state what is *not* modelled: excess and deficient
  contrast, relative darkness and HOLZ lines are dynamical; bands move rigidly with the crystal
  while spots do not, which is why bands navigate and spots identify; a thin foil may show strong
  spots and no bands at all; and the overlay is a prediction from the accepted orientation rather
  than independent evidence for it.
- Each candidate returned by `tem.solve_pattern` now carries its own `crystal_to_pattern`, not only
  the best-ranked one, because accepting a candidate is a judgement and everything drawn from the
  accepted solution needs *its* orientation.
- **A stereogram beside the TEM pattern.** `tem.stereogram` projects every zone axis of the phase
  onto the upper hemisphere in *holder* coordinates, so a pole's position on the drawing is the tilt
  that reaches it. It carries the holder envelope as the region of poles the stage can actually
  reach, marks the axis on the beam where the stage puts it, and — when a target is named — draws
  the geodesic to it, which is the connecting Kikuchi band, with the low-index zones along the way
  ringed as re-indexing waypoints. Hovering anywhere reports the α and β that would bring that point
  onto the beam, together with the pole under the cursor. The panel now lays the pattern and the
  stereogram side by side, which is where the previously empty right-hand half of the stage went.
- **The workbench reads user data files.** `pytex.app.uploads` is the one place the browser-to-reader
  gap is crossed: a file's text travels as an ordinary request parameter and is materialised as a
  temporary path for the library's own importer, so no format gains a second implementation.
  - **EBSD `.ang` and `.ctf` scans.** *Open a scan* replaces the practice dataset, and every control
    — colouring, scalar modulation, boundaries, grain threshold — acts on the imported map exactly as
    it does on a constructed one. Hexagonal scans, which EDAX writes by default, are drawn on a
    half-step raster: each measurement gets its own cell, the cells between them are drawn empty
    rather than interpolated. An imported map reports **no known answer**, where a practice dataset
    states one.
  - **Texture `.xrdml` pole figures.** `texture.measured_pole_figures` opens a set of Panalytical
    files and draws them in tabs, full size, on **one shared intensity scale** — two figures on two
    scales cannot be compared, and comparing them is the reason to measure more than one. Contour
    levels are set explicitly (`1, 2, 4, 7, 10`) or spaced automatically, and apply to every figure
    in the set. m.r.d., peak-normalised and as-recorded intensities are all available; only m.r.d.
    makes two instruments comparable. **Reconstruct the ODF** inverts the opened set into an
    orientation distribution and adds it as a further tab, sliced at the three phi-2 sections
    texture papers print — with its residual reported beside it, because the inversion is
    ill-posed and a peak without a residual is not a result.
- **Calibration from the image itself, in the TEM workspace.** The camera equation uses one number —
  how much reciprocal space a pixel spans — and an image that arrives without its recorded camera
  length usually still carries something whose length is known. **Calibrate** measures a line across
  it: a reciprocal length (Å⁻¹ or nm⁻¹, from a scale bar or a known reflection) sets that scale
  directly through the new *pixels with a measured scale* coordinate mode, and a real length on the
  plate (cm or mm) sets the pixel size instead. The scale can also simply be typed. Indexing through
  a measured scale gives the same answer as indexing through a camera constant with that quotient,
  which is checked on a practice plate both ways.
- **An About panel** (`pytex.app.about`): version, description, author and licence, served on the
  application manifest so the version shown is the version that answered the request. The licence it
  displays is checked against `pyproject.toml` and `LICENSE` by test.

- **An EBSD workspace in the workbench.** `ebsd.map` draws one orientation map from four
  independent choices — what the colour means, whether a scalar channel modulates it, whether the
  grain boundaries are superimposed, and what counts as a grain — so combinations nobody would
  enumerate as separate buttons stay reachable. Colourings are IPF along a chosen specimen
  direction, grain identity, GROD, KAM, and the measured confidence-index, fit and image-quality
  channels. Any of them can be greyed by any scalar channel, which is how an IPF map is made to show
  where the indexing should be believed without giving up the orientation; fit is inverted, since it
  is an error rather than a quality. The map travels as a raster at its native grid resolution — one
  pixel is one measurement, so nothing is interpolated — with the boundary network as line geometry
  on top, reconstructed onto the pixel faces so it stays sharp under zoom. Three practice datasets
  (`pytex.app.ebsd_gallery`) are constructions with answers fixed before the calculation runs: a
  bicrystal at exactly 40° with a linear deformation gradient, coherent Σ3 twins at exactly 60°, and
  a twelve-grain equiaxed polycrystal. Each result carries its dataset's known answer, so a user can
  check the numbers rather than trust them.

- **A CBED workspace in the workbench.** `cbed.pattern` simulates a zone-axis convergent-beam
  pattern over `pytex.diffraction.cbed` by either the two-beam closed form or the coupled
  Bloch-wave solution, and reports the overlap regime explicitly — the discs being separated is
  what makes each one an independent rocking curve, and therefore what decides whether a thickness
  can be measured at all. The discs arrive as one rasterised intensity image with a stated extent
  in millimetres, drawn under a vector overlay of outlines, labels and HOLZ rings, so the fringes
  survive while the geometry stays measurable under the cursor and hoverable per disc.
  `cbed.thickness_from_fringes` inverts the two-beam fringe minima for a foil thickness *and* an
  extinction distance, plotted as the straight line it is so a misassigned fringe order shows as a
  point off the line. `cbed.holz_rings` reports the higher-order Laue-zone radii and the
  reciprocal-lattice repeat along the beam — the one dimension a zone-axis pattern is blind to.
  Symmetry determination is offered only on a Bloch-wave pattern and refused on a two-beam one,
  because a two-beam pattern is symmetric in the excitation error by construction and its symmetry
  is the method's rather than the crystal's.

- **One centralized message log, in both shells.** The workbench now carries a severity-graded
  console pinned to the bottom of the window, and every module reports into it continuously —
  critical, error, warning, important, success, info, and progress ticks that carry a percentage
  and an ETA (`50% progress. ETA: 2 min 30 sec`). `pytex.app.logbook` defines the record, the
  severity vocabulary and a bounded thread-safe buffer; `pytex.app.contracts.execute` narrates
  every operation of every service without any service opting in; each call's records ride back on
  its own envelope, and the events that belong to no call are polled from `GET /api/log?since=`.
  Standard-library logging is bridged into the same stream, so what a desktop user could previously
  only have read in a terminal they do not have is now on screen. The console filters by severity
  and by text, groups a long task's ticks onto one line that counts up, and copies the whole
  session as plain text. It replaces the former activity strip, which counted calculations rather
  than describing them.

### Changed

- **One stereographic projection, in one documented public place.** `project_directions` and
  `fold_upper_hemisphere` now live in `pytex.core.sphere` and are exported from `pytex.core`, with
  the full docstring contract: both radial laws stated (`r = tan(rho/2)` stereographic,
  `r = 2 sin(rho/2)` equal area), when each is the right map, and what the antipodal folding means.
  `pytex.texture.projections` re-exports them, so no caller changed. The two private copies that
  had grown beside it — the publication tilt stereogram's `_project` and the TEM workbench
  stereogram's `_stereographic` — are gone; both call the shared helper. No drawn output changes,
  and the tangent-half-angle law is pinned directly on the helper in
  `tests/unit/test_sphere.py`.

### Fixed

- **The server could die mid-session under two concurrent operations.** `ThreadingHTTPServer` ran
  every request on its own thread, including the calculations, over a scientific stack that is not
  thread-safe — pyplot's state is global by construction. Two operations arriving together produced
  a native access violation that took the process down with no Python traceback, which became
  routine the moment a panel drew two figures on mount. Operations and exports now run one at a
  time; static files, the manifest and the log poll still answer concurrently, which is what the
  threading was for.
- **The fitted lattice was drawn outside the pattern.** A viewBox is a coordinate system, not a
  boundary: the grid is generated outwards until it certainly covers the image diagonal, and the
  overshoot was painted across the blank margins beside the picture, asserting lattice where there
  is no image. Every overlay is now clipped to the image rectangle, and a click in the margin is
  refused instead of placing a spot on no pixels.
- **A number field holding non-numeric text was reported as empty.** An `<input type="number">`
  whose content is not a number exposes an empty `value`, so text the user could see on screen
  reached the server as a missing parameter and came back as "this required field is required".
  The control now detects the browser's own `badInput` state and answers with what actually
  happened — `Invalid format of the input: only integers are allowed!` — beside the control and in
  the message log.

- **A lattice fitted to two or three picks was silently wrong.** `fit_planar_lattice` refined the
  beam centre whenever it was asked to, including where the picks cannot support it: the model has
  six free parameters and each spot gives two equations, so three spots make the system exactly
  determined and two make it rank-deficient. `lstsq` answers a rank-deficient system with the
  minimum-norm solution rather than failing, which shrinks the basis and displaces the origin while
  still passing through every pick. Two orthogonal picks 208 units from the beam produced basis
  vectors of 25 and 54 units, nodes (3, 4) and (−8, 1), and a centre 12.7 units from the click —
  with small residuals, so nothing looked wrong. The centre is now refined only from four spots up
  and held with a stated reason below that, every design matrix is rank-checked before it is
  solved, and two picks take an exact path: they *are* the basis, about the clicked origin,
  unreduced so that `a` still points at the spot picked first. `PlanarLatticeFit.centre_refined`
  reports which happened. **Any two- or three-pick lattice fit, and any beam centre taken from
  one, should be recomputed.**

- **A TEM pick made after zooming or panning landed in the wrong place.** The solver panel
  converted pointer events against the element's bounding box, which ignores the `viewBox` camera
  introduced with the zoom and pan controls, so every pick after a camera move was displaced by the
  camera offset. Picking now goes through the plot frame that owns the camera, and the view is
  preserved across the redraw that follows a pick, so zooming in to place a spot precisely works.

### Added

- **The TEM picks are editable coordinates, not only clicks.** Every pick — the transmitted beam
  included — appears in the rail as an editable `x` and `y` with its radius, d-spacing and angle
  from the first spot beside it. Type a coordinate, or select a row and nudge it with the pad or
  the arrow keys down to 0.1 px, and the lattice re-fits live. Rows can be removed or promoted to
  the transmitted beam, which swaps with the current beam rather than discarding it, and the whole
  set reads and writes as text so a measurement survives the session. A click is snapped to the
  spot it landed on — the intensity-weighted centroid on an uploaded micrograph, the known centre
  on a practice plate — while a click on background stays where it was put.

- **The pattern overlay says which marks are measured and which are fitted.** The fitted centre is
  drawn as its own mark joined to the clicked beam by a dashed line, and only when the fit actually
  solved for it, instead of the grid being drawn from one origin and the crosshair from another
  with nothing to say so. The picks, the fitted grid, the basis vectors and the calculated pattern
  now carry four distinct colours, and a basis arrow that lands on no pick is dashed rather than
  drawn as if it were a picked spot.

- **The figure fits the window, and its controls travel with it.** Every plot card is now sized by
  the visible stage rather than by the drawing inside it, so the complete figure is on screen when a
  panel opens instead of running past the bottom of the window. A legend that toggles what is drawn
  moved inside the card, under the figure, where it stays visible with the plot it controls. Zoom
  now runs below 100% as well as above it, and every plot — including the Crystal Viewer's 3-D
  camera — carries a pan tool alongside the existing Shift-drag and middle-drag. Presentation-only
  control groups start collapsed behind a settings mark, and the TEM solver's four workflow steps
  are one accordion that advances with the work.

- **Measured picks on the TEM pattern.** The solver draws the picked spots' radius, d-spacing, ratio
  to the first pick and angle from it in the top-left of the pattern, taken from the clicked
  coordinates and the calibration with no solution involved — the numbers a zone axis is identified
  from, beside the pattern while picking rather than in a table under it.

- **Named-component ODF fitting.** `fit_odf_components` fits declared texture-component kernels
  plus an optional random term with non-negative fractions constrained to sum to one.
  `ODFComponentFit` carries observed/predicted normalized densities, RMS/maximum residuals, $R^2$,
  `describe()`, and a JSON contract; rank-deficient supports raise. An executable analytic 70/30
  cube-Goss mixture pins the fractions and zero residual.

- **Finite-thickness SAED shape factors.** `FiniteThicknessShapeFactor` implements the normalized
  plane-parallel amplitude `sinc(t s_g)` and intensity `sinc^2(t s_g)`, with its first zero,
  `describe()`, and JSON round trip. Both vectorized SAED engines accept a physical foil thickness;
  the legacy Lorentzian remains available only when thickness is absent. The analytic 100 angstrom
  slab landmarks are pinned by tests and an executable worked example.

- **The beam centre is solved for, not clicked.** `pytex.diffraction.lattice_fit`
  fits the plane lattice that a zone-axis pattern's spots must lie on, which
  over-determines the transmitted beam with four or more picks and names any pick
  the lattice cannot explain. This matters because the camera equation measures
  every radius *from the beam*: an error there biases every spacing in the same
  direction and yields a self-consistent answer for the wrong material rather
  than an obvious failure. Exposed as `tem.fit_lattice`, and drawn live over the
  pattern in the workbench, with the two basis vectors as labelled arrows from
  the beam to the picks that generate them.

  - Seeding uses *differences between spots*, never offsets to the picked centre:
    a difference is a lattice vector however badly the centre was picked, while
    an offset is one only if the pick was already right. The alternative fails in
    exactly the case the method exists to repair.
  - Candidate lattices are ranked by evidence rather than by inlier count, since
    halving a cell explains every spot it explained before plus the mis-picked
    one. A cell of area `A` puts `π t² / A` of the plane within tolerance of a
    node, so each inlier carries `log(A / π t²)` nats and a denser lattice pays
    `log 4` per inlier for its extra nodes.
  - The reported basis is Gauss-reduced, so its lengths and included angle are
    lattice invariants. Before that, a square lattice could be reported as two
    vectors 135° apart, which is correct and useless.
  - Two limits are documented and tested rather than papered over: a centre wrong
    by an exact lattice vector is undetectable from geometry, and refinement is
    leashed at half a spacing because beyond that it relabels which node the
    origin is.

- **Candidate solutions carry deviations, a configurable score, and their own
  calculated pattern.** `pytex.diffraction.solution_scoring` reports measured
  against calculated d-spacings per spot and angles per pair — the evidence — and
  fuses them through an explicit, documented `ScoringWeights` policy that travels
  with every number it produces. `tem.solve_pattern` now ranks candidates by that
  score rather than by the solver's sort key, which is deliberately not a
  quality, and says so when the two orders disagree.

  - Angles are weighted above lengths because a wrong camera constant scales
    every length and leaves every angle untouched: an angular disagreement is
    evidence about the crystallography, a length disagreement may only be
    evidence about the instrument. A 5 percent calibration error moves the length
    deviation to 5.01 percent and the angle deviation not at all.
  - Every candidate returns the pattern it predicts, in picking coordinates, so a
    solution can be accepted by looking at whether it lands on the measured
    spots. Accepting is explicit, and is what carries the phase and axis into
    zone-axis listing and tilt planning.
  - The per-spot table gained a Δd column. The *same* deviation on every spot is
    the signature of a calibration error; a scatter of them is an indexing error.

- **A human-readable export, on every result in every panel.** `Report` writes a
  Markdown page — the answer in prose, the caveats, the data, the exact inputs,
  the citations — which is the one thing CSV, XLSX and JSON are each unsuited to.
  `EXPORT_FORMATS` is now published in the manifest and the browser generates its
  export buttons from it, so a format added in Python reaches every panel without
  a matching edit in JavaScript.

- **The TEM panel can be used without a micrograph, and its answers can be
  checked.** Three practice SAED plates ship with the workbench — aluminium fcc
  down [001], ferrite bcc down [110], zirconium hcp down [2̄110] — each a real
  kinematic calculation projected onto a stated detector raster, so the whole
  workflow (calibrate, pick the beam, pick the spots, index) runs on a pattern
  whose zone axis is known by construction. `tem.solve_pattern` gained an
  optional `expected_zone_axis`, and the verdict it returns compares the two
  directions **up to symmetry**: a bcc [110] pattern is indistinguishable from a
  [101] one, and calling that a mismatch would be wrong.

  - `pytex.tem.synthetic` projects a simulated zone-axis pattern onto a
    `DetectorRaster` and returns pixel coordinates, relative intensities and
    display radii. The raster is the recorded image — no handedness flip — so a
    picked `(column, row)` is the detector `(X, Y)` divided by the pixel pitch.
    The beam may sit off the geometric centre, the pattern carries a roll about
    the beam, and a seeded sub-pixel scatter makes indexing residuals realistic.
    Positions are exact; intensities are kinematic and double diffraction is not
    modelled.
  - `SyntheticSAEDImage.independent_seed_spots` picks strong reflections whose
    directions from the beam all differ. Selecting the top *n* by brightness
    alone returns Friedel pairs, which are collinear through the beam and cannot
    seed an index.
  - The gallery's camera constant is *computed* as `L·λ` from the camera length
    and accelerating voltage, not typed, so the calibration field and the
    geometry cannot drift apart.

- **A zone-axis atlas: which axis to tilt to, not just whether you can reach the
  one you named.** `pytex.tem.atlas.zone_axis_atlas` enumerates the
  symmetry-distinct zone-axis families of a phase with the angle from the axis on
  the beam, the family size, the number of reflections inside a stated cut-off,
  and the pattern's apparent *n*-fold symmetry — the last measured on the
  simulated spot set rather than deduced from the point group, so it reports what
  an operator will actually recognise on arrival. Exposed as
  `tem.zone_axis_atlas`, whose reachability column is computed by the same
  planner `tem.plan_tilt` uses, against the same envelope and roll.

### Fixed

- **⟨110⟩ could vanish from a 45° zone-axis search.** The angle between two zone
  axes is an arccos through a basis product and lands a few ulps either side of
  exactly 45°, so a bare comparison against the search limit dropped the single
  most-wanted cubic target about half the time.
- **The axis already on the beam was offered as somewhere to tilt to.** The same
  arccos is square-root-behaved near 1, turning 1e-16 of cosine error into ~1e-6
  of a degree — enough to pass a `> 1e-6` test. For ferrite [110] the atlas
  summary read "the nearest reachable one is [110] at 0.00°".
- **A catalogue phase carried into the TEM panel was renamed "(edited)".** The
  gallery sent an expanded `PhaseSpec`, and the phase picker treats a full
  description as a user-edited phase; the indexed result was then titled
  "Aluminium (fcc) (edited)" on a phase nobody had touched. The gallery now sends
  a catalogue reference, and the panel carries forward the value the picker
  emitted rather than the description the result echoes back.
- **A tilt plan counted "members of [012]".** A specific direction has no
  members; the sentence that counts symmetry-equivalent members now writes the
  family form ⟨012⟩, per the notation standard, while the title keeps the
  direction the user typed.
- **`tem.plan_tilt` failed with a 500 on a high-index hexagonal target.**
  `TiltSolution.orbit_member_indices` is `None` when the direction placed on the
  beam has no low-index integer form within the navigation module's bound, which
  happens routinely for a hexagonal family such as ⟨4 3̄ 1⟩; the service read it
  as an integer array regardless. The plan itself was correct — only the name of
  the member was unavailable — so the row now reports the *family* form, which is
  what is actually known.

### Changed

- **The canonical scientific notes are MyST Markdown, and their mathematics now
  renders on the site.** The 37 theory, algorithm, and validation notes moved
  from `docs/tex/*.tex` to `docs/site/theory/*.md`. Sphinx had no LaTeX-parsing
  extension, so every note had been copied verbatim into `_downloads/<hash>/`
  and surfaced as a download link: the derivations — the actual content — were
  invisible to a reader of the documentation. They are now ordinary pages,
  grouped in `theory/index.md` and reachable from search.

  - No new dependency: `myst_enable_extensions` already carried `amsmath` and
    `dollarmath`. Equation labels became MyST labelled equations, so
    cross-referenced equations render numbered and link to their anchors.
  - PDF output comes from Sphinx itself, `sphinx -b latexpdf docs/site`,
    replacing a per-document `latexmk` step. Because it renders the same
    sources as the HTML, the print and web forms cannot drift apart.
  - The conversion was safe to automate: the corpus used no `\cite`, `\ref`,
    `\includegraphics`, `\newcommand`, `\input`, or TikZ, and 26 of the 37
    files had no `\documentclass` at all — they were fragments, so the
    documented `latexmk` build could not have run for most of them.
  - Content is unchanged. Four notes overlap `docs/site/algorithms/` pages,
    which are a separate documentation layer by design; both are kept and
    `theory/index.md` cross-links each pair.
  - `docs/standards/latex_and_figures.md` is renamed
    `docs/standards/scientific_notes_and_figures.md`, and `pytex docs inventory`
    now lists the notes rather than `.tex` files.

### Added

- **The Kearns parameter, computed four ways from one tensor.** `f` is the scalar
  texture index zirconium components are specified against, and the four
  techniques that measure it do not always agree. `pytex.texture.kearns`
  implements all four as estimators of a single object — the second-moment
  tensor of the basal-pole distribution, `A = <c c^T>`, with
  `f(d) = d^T A d` — which turns two facts the literature states empirically
  into identities: the values along an orthonormal triad sum to `tr(A) = 1` for
  every texture, and a random texture gives `1/3` in every direction. A measured
  triad that misses 1 is therefore reporting the systematic error of the
  measurement, and needs no reference specimen to say so.

  - `kearns_from_orientations` (EBSD or simulation, exact),
    `kearns_from_odf`, `kearns_from_pole_figure` (Baron *et al.* 1990 Eq. 5),
    `kearns_from_tilt_profile` (Kearns 1965 Eq. 5) and
    `kearns_from_diffractogram` (the original theta-2theta route), plus
    `pole_orientation_tensor`, `basal_tilt_angle_deg` — computed from the phase's
    reciprocal metric rather than a transcribed table — and
    `harris_texture_coefficients`.
  - `KearnsReport.describe()` states the route, the values against the `1/3`
    reference, the triad sum, the measured tilt coverage, and the caveats each
    route carries; `to_json()` is kept in lockstep with it.
  - `kernel_axis_shrinkage` gives the closed-form factor
    `beta = (3 rho - 1)/2` by which a kernel-density ODF's own smoothing pulls
    every departure from `1/3` toward isotropy — 0.94 at a 10 degree halfwidth,
    0.78 at 20, so an `f` of 0.70 reads as 0.62 there. Derived from Rodrigues'
    formula and the moments of a uniform rotation axis, and checked against a
    Monte-Carlo integral over SO(3). Report the halfwidth with any `f` taken
    from an ODF.
  - Pinned against the sources of record: computed basal tilt angles reproduce
    Kearns' Table 2 to 0.2 degrees, and his Table 3 longitudinal-section
    intensities reproduce his tabulated `f = 0.488`. His transverse-section
    block does not, and the theory note says why: one cell reads 0.0214 where
    `0.353 cos^2(75 deg) = 0.0237`. Recomputing his transverse block gives 0.0529
    from his own volume fractions and 0.0526 from his intensities, against the
    quoted 0.0508.
  - Theory note `kearns_parameter_and_basal_pole_texture`, tutorial notebook 31,
    and four worked examples whose expected values are elementary consequences
    of the definition rather than prior outputs of this code.

- **Double diffraction in the kinematic SAED engine.** A kinematic simulation
  shows a real pattern's spots minus the ones dynamical scattering puts there,
  and the largest such class is double diffraction: a beam diffracted by `g1`
  is itself an incident beam inside the crystal, so diffracting it again by
  `g2` sends it out along `g1 + g2`, and a reflection whose structure factor
  vanishes appears anyway. Silicon 002 along [110] is the standard case.

  - `KinematicSimulationConfig(include_double_diffraction=True)` adds those
    reflections. The selection rule — the pairwise integer sums of the excited
    reflections, exposed as `pytex.diffraction.kinematic.double_diffraction_sums`
    — is exact, because it follows from the additivity of scattering vectors
    and not from the dynamical theory the engine cannot solve. Off by default;
    every existing pattern is unchanged.
  - They are never mixed in unlabelled. `SpotTable.is_double_diffraction` marks
    them, `double_diffraction_parents` records the strongest contributing pair,
    `double_diffraction_origin_label(row)` names the path as `g1 + g2 = g`,
    `describe()` states the designation, the exported reflection table gains
    `double_diffraction` and `double_diffraction_origin` columns (appended, so
    the declared `REFLECTION_TABLE_COLUMNS` contract is preserved), and
    `render_composite_saed` draws them hollow in a separate collection with its
    own legend entry.
  - The intensity is `c * sum over paths of I1 * I2`, an observability estimate
    scaled by `double_diffraction_coupling`, documented as such everywhere it
    appears. The kinematic intensity of a forbidden reflection is exactly zero,
    so there is nothing else it could be.
  - Worth knowing: this can never revive a **centring** absence. Centring
    conditions define a sublattice of reciprocal space and a sublattice is
    closed under addition, so only basis absences — glide plane, screw axis,
    motif — can be revived. That is exactly what is observed, and it is pinned
    by a test for both I- and F-centred lattices.
  - Validated on silicon [110], where the forbidden 002 appears at exactly half
    the detector radius of 004 (worked example
    `kinematic-silicon-double-diffraction-002`), with theory in
    `docs/site/theory/reciprocal_space_and_kinematic_spots.md`.

- **Dynamical CBED: many-beam coupling, absorption, HOLZ lines, and the point
  group.** Tutorial-28-era CBED computed each disc as an independent two-beam
  calculation, and said so: no coupling, no absorption, no HOLZ lines, and no
  diffraction-group symmetry determination. All four now exist, and they turn
  out to be one capability rather than four — admitting higher-order Laue zone
  beams is what produces the lines *and* what breaks the projection symmetry
  that would otherwise make every pattern look centrosymmetric.

  - `pytex.diffraction.dynamical` solves the coupled many-beam problem by the
    Bloch-wave method. The absolute scale is inherited rather than re-asserted:
    the coupling coefficient is `nu_g = lambda F_g / (pi V_c cos theta_g)`, so
    `|nu_g| = 1/xi_g` for the extinction distance already validated against
    Williams and Carter Table 23.1, and the two-beam limit reproduces
    `two_beam_rocking_curve` to 2e-15 — which pins the diagonal convention, the
    off-diagonal scale and the `i*pi` in the propagator at once. Without
    absorption the propagator is unitary and the beams sum to one to 1e-12,
    which is the check that catches obtaining the Bloch-wave excitation
    amplitudes by projection when the eigenvectors are not orthogonal.
  - `AbsorptionModel` adds the imaginary optical potential. The *structure* is
    exact — anomalous absorption emerges from the eigenvector structure rather
    than being applied to the output, and is tested by the
    Hashimoto-Howie-Whelan theorem that the bright-field rocking curve becomes
    asymmetric while the dark-field one stays symmetric. The *magnitudes* are
    the customary phenomenological ratios of Hirsch et al., and `describe()`
    says so; absorptive form factors are not computed. Normal absorption is
    proved to factor out exactly as `exp(-2 pi t / xi'_0)`, so the
    phenomenological number cannot contaminate any conclusion about shape,
    position or symmetry.
  - `pytex.diffraction.holz` gives the HOLZ line loci in closed form, with
    their chords in both discs, their angular width for a given thickness, their
    intersections, and their sensitivity to strain and to voltage. Positions are
    checked against `pytex.diffraction.dynamical`, which derives the excitation
    error independently, to 1e-15. The strain/wavelength degeneracy is asserted
    rather than described: a fractional lattice strain and a fractional
    wavelength change of the same size cancel to 1e-16 at every reflection
    simultaneously, which is why quantitative HOLZ metrology calibrates the
    accelerating voltage first.
  - `pytex.diffraction.diffraction_groups` **derives** the 31 diffraction
    groups from PyTex's own operator tables rather than storing Buxton's table,
    and obtains the point-group correspondence by inversion. Bright-field and
    whole-pattern symmetry are derived too. The centrosymmetry statement is an
    exact correspondence — `2_R` is present at every beam direction of a
    centrosymmetric crystal and none of an acentric one, asserted over all 32
    point groups — so the `+-g` observation alone splits them into exactly 21
    and 11.
  - `ConvergentBeamConfig` gains `method`, `absorption`, `laue_zones` and the
    HOLZ search bounds; `CBEDPattern` gains `beam_set`, `holz_lines`,
    `predicted_diffraction_group()`, `symmetry_observations()` and
    `determine_point_group()`. Down `[001]`, zincblende GaAs and diamond
    silicon are separated by whole-pattern symmetry, `2mm` against `4mm`, and
    the determination returns `{-42m, -43m}` and "not centrosymmetric" for the
    polar one. Confine the beam set to the zeroth Laue zone and the same
    crystal looks centrosymmetric, so `symmetry_observations()` refuses a
    projection calculation unless asked twice and a two-beam pattern outright.
  - Tutorial `29_dynamical_cbed_and_point_groups`, the theory note
    `docs/site/theory/dynamical_cbed_and_symmetry_determination.md`, and six
    worked examples whose expected values are analytic identities or published
    counts.

  Not implemented, and stated in `describe()` rather than left to be
  discovered: Bethe perturbation of weak beams, computed absorptive form
  factors, Buxton's dark-field and `+-g` observations for reflections on
  symmetry lines, and specimen realism (wedge, bending, strain gradient, probe
  aberration, inelastic background).

- **Pole-figure arithmetic.** Two pole figures could not previously be combined
  at all: `PoleFigure` holds *scattered* specimen directions, so two figures
  generally share no sampling direction, and the arithmetic was blocked
  structurally rather than merely unwritten. The chain is now complete.

  - `PoleFigure.sampling` records whether the intensities are per-pole weights
    of a pole cloud (`"scattered_poles"`) or densities evaluated at the given
    directions (`"sampled_density"`). The two demand different resampling
    estimators — a weighted sum versus a weighted mean — and using the sum on a
    latitude-longitude raster biases the result towards the poles, where such a
    raster oversamples. Round-trips through the JSON contract; payloads
    predating the tag read as pole clouds, which is what `from_orientations`
    produced.
  - `PoleFigure.on_grid` resamples onto any `S2Grid` by kernel smoothing,
    choosing the estimator from `sampling`. This is what gives two figures a
    common support.
  - `PoleFigure.normalize_to_mrd` and `PoleFigure.spherical_mean` put
    intensities on the multiples-of-random scale, backed by the new
    `pytex.core.sphere.raster_solid_angle_weights` for measured rasters. The
    XRDML and LaboTex readers gain an `"mrd"` normalization mode.
  - `PoleFigure.__add__`, `__mul__`, `__truediv__`, `difference`/`__sub__`,
    plus `rotate`, `symmetrize` and `restrict_polar_range`. Every way two
    figures can differ while still looking combinable — pole, specimen frame,
    antipodal convention, family flag, support — raises with the reason.
  - `PoleFigureDifference` holds a signed residual. Subtraction does not return
    a `PoleFigure` because a pole density is non-negative by invariant while a
    difference is signed, and the sign is its entire content.

- **Residual pole figures for ODF inversion.** `PoleFigureResidualReport` gains
  `difference_figure()` and a `describe()`; `pytex.plotting.plot_pole_figure_difference`
  draws it on the diverging colormap with limits symmetric about zero.
  A residual norm says how badly an ODF misses its own input data; the residual
  *figure* says where, which is what distinguishes an unmodelled component from
  counting noise. `ScatterLayer2D` gains `vmin`/`vmax` to support this.

### Fixed

- **`raster_solid_angle_weights` over-weighted the equatorial ring of a
  hemispherical raster.** The outermost ring's band was extended outwards by its
  own half step, which for a raster ending exactly at 90 degrees pushed it past
  the equator and gave that ring close to twice the solid angle it owns. Because
  `cos^2` vanishes there, it dragged every Kearns integral down: the spherical
  mean of `cos^2` over a 5 degree raster came out at 0.3196 against the exact
  `1/3`, a -4.1 percent error. A new optional `polar_max_deg` bounds the band and
  takes the error to -0.06 percent. **The default is unchanged**, so existing
  results and the pinned `pole-figure-raster-weighted-mean-converges` worked
  example are unaffected; `kearns_from_pole_figure` passes `polar_max_deg=90.0`
  for antipodal figures.

- **The harmonic ODF could not use a bandwidth above 6.** `HarmonicODF` builds
  its basis from Wigner small-`d` functions, whose coefficient is a ratio of
  factorials. At degree 7 the numerator `factorial(2*l)` already reaches 8.7e18
  and the formula multiplies four such terms, so the product became a Python big
  integer that NumPy could hold only as an object — and `np.sqrt` on an object
  array raises `TypeError`. Any `degree_bandlimit >= 7` therefore crashed. Since
  the default is 6, this surfaced only when a user raised the bandwidth to
  resolve a sharp texture, which is exactly when it is needed. The coefficient
  is now evaluated in log-gamma, which agrees with the exact integer form to
  2.7e-14 over every term of degrees 0-6 and stays finite at any degree.

- **Pole-figure residuals compared two different scales.**
  `PoleFigureResidualReport.from_odf` subtracted a discrete `ODF`'s pole density
  straight from the measured figure, but `ODF.evaluate_pole_density` returns a
  kernel-weighted *response*, not multiples of random — a random texture returns
  the kernel's spherical mean, about 0.016 at a 12 degree halfwidth. A perfect
  fit therefore reported a relative residual near 1.0, condemning every sound
  inversion. The prediction is now converted with the new
  `pytex.texture.models.random_pole_density`, which also replaces the private
  duplicate in `pytex.diffraction.preferred_orientation`. `HarmonicODF` is
  unaffected: its densities are already in m.r.d., so that path is left alone.

- **Comparing two Miller indices raised instead of answering.** `MillerIndex`,
  `CrystalDirection`, `ZoneAxis` and `ReciprocalLatticeVector` inherited the
  generated dataclass `__eq__`, which compares their index arrays with `==` and
  raises `ValueError: The truth value of an array with more than one element is
  ambiguous` for every distinct-but-equal pair. `a == b` on two separately
  constructed but identical indices was therefore impossible, which blocked any
  operation that must first check two objects describe the same quantity.
  `SymmetrySpec` already carried the custom `__eq__`/`__hash__` for this
  reason; these four now do too.

- **`import pytex` failed on a clean install of the declared dependencies.**
  `pytex.diffraction.solving` and `pytex.plotting.styles` import `yaml` at module
  level, but `pyyaml` was not declared, so a fresh `pip install pytex` raised
  `ModuleNotFoundError` on the very first import. `pyyaml` is now a declared
  runtime dependency — it backs two data contracts on the core import path, the
  plotting style themes and the measured-SAED pattern format. Verified by
  building the wheel and importing it in a clean virtual environment.

- **Matplotlib was silently mandatory.** `pytex.plotting.crystal3d` and
  `pytex.plotting.scene3d` imported `matplotlib.colors` at module level, which
  contradicted the `plotting` extra and the `_require_matplotlib()` guards used
  everywhere else. Those imports are now lazy, so the library imports and
  computes with numpy, scipy and pyyaml alone; matplotlib is required only when a
  plot is actually drawn.

- **The version string existed in four files.** `pyproject.toml`,
  `pytex.__version__`, and both manifest writers each carried their own literal,
  so a release bump would have left exported manifests stamping the previous
  version. `src/pytex/_version.py` is now the single source: the packaging
  metadata reads it statically, and `CITATION.cff` states the same version.

- **`nye_dislocation_density_tensor` destroyed its own measurable output.**
  Subtracting the trace as `trace * identity` propagated the `NaN` of the
  unmeasurable out-of-plane curvature into every off-diagonal component, because
  `NaN * 0` is `NaN`. The trace is now applied to the diagonal alone.

- **ODF-weighted preferred-orientation factors were mis-scaled by ~40x.**
  `ODF.evaluate_pole_density` returns a kernel-weighted response, not a value in
  multiples of random — the smoothing kernel peaks at 1 rather than integrating
  to 1 — so a uniform texture returns the kernel's spherical mean. The
  correction now divides by that mean, computed by Gauss-Legendre quadrature, and
  a uniform ODF gives a factor of 1 as it must.

- **`GnomonicProjection.contains` rejected exact detector edges**, so
  `detector_corner_coordinates()` reported its own corners as off-detector. Now
  takes an explicit pixel tolerance.

- **Five `pytest.raises(match=...)` patterns contained unescaped regex
  metacharacters**, silently matching more loosely than they read.

### Added

- **Kikuchi band geometry and the gnomonic projection**
  (`pytex.diffraction.kikuchi`): `GnomonicProjection`, `KikuchiBand`,
  `KikuchiZoneAxis`, `KikuchiPattern`, `simulate_kikuchi_pattern`, and
  `plot_kikuchi_pattern`. The geometric layer shared by EBSD and TEM. Band
  centre lines are exactly straight in gnomonic coordinates at any detector tilt;
  band edges are computed on the **exact Kossel cones**, so they are the conics
  they physically are rather than the usual small-angle straight-line
  approximation. Validated against closed-form anchors: gnomonic radius
  `tan(psi)`, Bragg's law by hand for the 2.42-degree Ni{111} band at 20 kV, and
  `[011]` projecting to gnomonic radius exactly 1 at the cube orientation, which
  pins the whole crystal-to-detector frame chain.

- **Preferred-orientation corrections for powder intensities**
  (`pytex.diffraction.preferred_orientation`): `march_dollase_factors`,
  `MarchDollaseModel`, `ODFPreferredOrientationModel`, the
  `PreferredOrientationModel` protocol, `apply_preferred_orientation`, and a
  `preferred_orientation=` argument on `generate_xrd_pattern`. The ODF-weighted
  model drives powder intensities directly from a measured texture with no fitted
  parameter — the texture core feeding diffraction. March-Dollase factors are
  averaged over the full symmetry family, so they cannot depend on which family
  representative the enumeration emitted, and the distribution's exact spherical
  normalization is pinned as a test.

- **Lattice curvature and GND density** (`pytex.ebsd.gnd`):
  `lattice_curvature_tensor`, `nye_dislocation_density_tensor`,
  `geometrically_necessary_dislocation_density` (Nye and KAM routes), and
  `plot_gnd_density_map`. Completes the KAM/GROD/GOS/GAM family with the
  dislocation content those gradients imply, in m^-2. Unmeasurable components are
  reported as `NaN` rather than zero: a 2-D map determines six of nine curvature
  and five of nine Nye components. Densities are documented as lower bounds and
  as resolution dependent, with the step-size scaling pinned as a test.

- **`phase_fixtures_available()`**, so code can test for the checksum-pinned
  phase-fixture corpus rather than discover its absence by exception. The corpus
  is a repository asset and is not shipped in the wheel; the loader now says so
  and names the alternatives.

- **Complete public-API docstring coverage.** Every one of the 411 names in
  `pytex.__all__` and every public member of every exported class now carries a
  docstring — 574 members plus 118 class docstrings that previously had only
  Python's auto-generated dataclass signature. Repo-wide public docstring
  coverage rose from 33.1% to 90.2%.

- **Three new test ratchets** that would have caught the defects above:
  `tests/unit/test_public_api_docstrings.py` (every export documented, with a
  real summary line), and `tests/unit/test_release_metadata.py` (no undeclared
  module-level imports, one version literal, citation metadata in step).

- **New Sphinx section `docs/site/algorithms/`** documenting how each scientific
  surface computes what it computes. Four pages — OR determination from measured
  orientations, variant-resolved plane and direction correspondence, composite
  SAED assembly, and SAED pattern indexing — each carrying the mathematics
  rendered on the page, the algorithm as reimplementable steps, the constraints
  and failure modes beside the step they govern, the cost, and a verification
  table naming where every claim is checked. Cubic (Kurdjumov-Sachs) and
  hexagonal (Burgers) examples run side by side throughout, so hexagonal-specific
  behaviour is stated rather than left implicit.

  Every number on these pages is computed rather than recalled, and several are
  new to the documentation: the symmetry-reduced catalog separation matrix for
  the fcc→bcc family, the equivalence-group breakdown of the KS `(111)` and
  Burgers `(011)β` correspondence tables, the per-variant child-zone deviations
  for Burgers down β[110], and the solver's noise envelope.

- **Two canonical theory notes** the transformation program had promised:
  `docs/site/theory/orientation_relationship_determination.md` (double-coset
  seeding, symmetry-aware rotation averaging, catalog distance, the
  non-uniqueness of a parallelism statement) and
  `docs/site/theory/saed_ratio_angle_indexing.md` (calibration, the
  admissibility test, triad construction, and the intrinsic zone-sense
  ambiguity).

- **`pytex.plotting.algorithm_diagrams` and `pytex.plotting.svg_primitives`**,
  with `scripts/generate_algorithm_figures.py`. Algorithm flow sheets are
  *generated* rather than drawn, so a figure cannot drift from the algorithm it
  illustrates, and they are held to the repository's figure layout guards —
  title/desc, canonical font, absolute marker units, no text overflow or
  collision, and byte-for-byte reproducibility. The shared primitives were
  extracted from `frame_diagrams`, whose figures are byte-identical afterwards.

- **One call now answers "I measured two phases by EBSD — what is the
  orientation relationship?"** `characterize_orientation_relationship(parents,
  children)` fits the operative rotation, names it against the standard catalog
  for the two crystal systems, states it as parallel planes and directions, and
  judges whether the identification can be trusted — returning an
  `ORCharacterizationReport` with `describe()` and a JSON payload.
  `orientation_relationship_from_euler(...)` is the same surface taking two
  columns of Bunge Euler angles in degrees, which is how measurements arrive.

  The fit needs no nominal relationship. The starting estimate comes from the
  data: the first pair is reduced to its minimum-angle representative in the
  double coset $G_c V_0 G_p$, which absorbs the parent symmetry operation that
  distinguishes one variant from another, so pairs drawn from different
  variants still align to it. **Only one pair is reduced, deliberately** —
  reducing every pair and averaging looks more robust and is not, because the
  maximum-trace element is not unique when the relationship's own rotation is
  symmetric. Bain is the concrete failure: 45 deg about <100> with three
  variants averages to a meaningless 26.9 deg, which then reads as
  Kurdjumov-Sachs. A regression test pins the correct behavior.

  The verdict is deliberately conservative. `is_conclusive` requires the winner
  to fit within tolerance *and* to lead the runner-up by more than both the
  measurement scatter and its own misfit. Measured behavior on planted KS data:
  identified conclusively up to 2 deg of added scatter; at 5 deg — comparable
  to the 2.40 deg that separates KS from Greninger-Troiano — the report
  degrades to an explicit "NOT conclusively identified" rather than to a
  confident wrong answer.

- **`describe_orientation_relationship(relationship)` reads a rotation back as
  crystallography**, recovering the parallel-plane and parallel-direction
  clauses that define it: `(111) || (011)` with `[10-1] || [11-1]` for
  Kurdjumov-Sachs, `(011) || (0001)` with `[-111] || [-12-10]` for Burgers,
  with hexagonal phases labeled in four-index Miller-Bravais form. A rotation
  typically satisfies several exact low-index parallelisms at once, all true;
  which one the literature quotes depends on the structures rather than on the
  rotation, so the search takes a preference — by default the relationship's own
  recorded defining families, and for a fitted relationship those of the
  catalog member it matched. The reported deviations then verify the statement
  against the fitted rotation instead of asserting it.

- **Tutorial notebook 23, `23_transformation_crystallography_end_to_end`**,
  committed executed. Answers the five questions a transformation study actually
  asks, on Burgers β→α in one pass: the orientation relationship determined from
  measured Euler angles with no nominal supplied, the variant correspondence
  table for an arbitrary parent plane, the composite SAED down [110]β with its
  reflection table and manifest, the same composite re-anchored on [0001]α of one
  variant, and a measured pattern solved back to phase, zone axis, spot indices
  and transformation variant. Each section states the crystallography first and
  then computes it, and closes with what was *not* shown — kinematic only, no
  HOLZ, no dynamical intensities, synthetic validation, no MTEX claim.

- **Measured SAED patterns can now be solved.** New module
  `pytex.diffraction.solving`: `solve_saed_pattern(pattern, phases)` determines
  the phase, the zone axis, the crystal orientation in the pattern frame, and the
  Miller indices of every spot, from picked spot positions plus a camera
  constant. `MeasuredSAEDPattern` reads and writes a documented YAML contract
  (`schemas/measured_saed_pattern.schema.json`), `solve_saed_pattern_file` does
  both in one call, and `assign_transformation_variant` names which variant a
  solved product pattern belongs to when the parent's orientation is known.

  The algorithm is classical ratio/angle indexing: two non-collinear reflections
  fix the zone, so the two shortest measured vectors seed the solution, and every
  remaining spot is indexed by projection and scored. **Intensities are never
  used** — a kinematic intensity model is not reliable enough to index against,
  and a printed pattern rarely carries calibrated intensities at all.

  This is a *different* surface from `diffraction.models.index_saed_pattern`,
  which starts from a calibrated `DiffractionGeometry` and works in detector
  pixels. Use that one when the detector model is known; use this one when only
  the pattern is. Both docstrings say so.

  Three honesty properties are built in rather than bolted on. A single SAED
  pattern cannot distinguish a zone axis from its reverse for a centrosymmetric
  reflection set, and the report says so instead of presenting one sense as the
  answer. Symmetry-equivalent descriptions are deduplicated and rewritten into
  the conventional one, so a cubic cube-axis pattern reports `[001]` rather than
  the equally valid `[0-10]` the seed search happened to find first. And *not*
  solving is a legitimate outcome: `best()` raises rather than guessing, and
  `is_conclusive` is `False` whenever spots are left unindexed or a genuinely
  different candidate explains the pattern equally well.

- **`pytex.plotting.saed_picker` picks spots interactively.** `SAEDSpotPicker`
  displays a pattern and collects clicks (left add, right remove nearest, middle
  set beam centre, `u` undo, `c` clear), producing a `MeasuredSAEDPattern` or a
  YAML file. The picking *logic* is `SpotPickerState`, a plain object with no
  Matplotlib dependency, so it is tested headlessly and the GUI is not on the
  critical path — an interactive tool that cannot be tested is a liability in a
  scientific library.

- **Composite patterns can be anchored on a product-variant zone axis.**
  `simulate_composite_saed_from_child_zone(relationship, child_zone_axis,
  anchor_variant_index=k, ...)` matches how the microscope is used: the operator
  tilts to a low-index zone of the *product* and wants the matrix and the sibling
  variants around it. The anchor variant's rotation `R_k` carries parent Cartesian
  vectors into that child's frame, so the requested child zone corresponds to the
  parent direction `R_k^T z_c` — generally irrational, and reported exactly
  alongside its nearest rational label, the same honesty child zone axes already
  got. `align_child_g` works in the child's own indices.

  Because the geometry then goes through the same `zone_basis_from_axis` call the
  parent-anchored path uses, there is one detector-geometry definition and a
  testable identity: anchoring on variant `k`'s image of a parent zone reproduces
  the parent-anchored pattern for that zone exactly (verified to 1e-13 mm, and
  pinned as a worked example). `simulate_composite_saed` accordingly accepts an
  irrational `CrystalDirection` as well as a `ZoneAxis`, and
  `CompositeSAEDPattern` records `anchor_variant_index` and
  `nearest_parent_zone_axis` so every export states which crystal defined the
  geometry.

- **Composite SAED patterns now export.** New module `pytex.diffraction.export`:
  `composite_reflection_table(pattern)` produces one row per rendered spot —
  source (parent or variant `k`), phase, Miller indices and label, `d`, `|g|`,
  detector position and radius, excitation error, `|F|`, relative intensity —
  with `to_csv`, `to_markdown`, `to_records`, `to_json_dict` and `describe()`.
  `export_composite_saed(pattern, directory)` writes the table, the rendered
  figure in the requested formats, the parent/child coincidence table, and a
  JSON manifest validated by the new
  `schemas/composite_saed_manifest.schema.json`. Every table value is read from
  the `SpotTable` the engine produced, so the table and the figure cannot
  disagree, and figures are closed after writing.

- **`CompositeSAEDPattern.centering_audit()` and
  `phase_centering_is_declared(phase)` expose a silent failure mode.**
  `ReflectionCondition.from_phase` reads the lattice centering from the first
  letter of a phase's space-group symbol and falls back to primitive when the
  phase carries none — so a body-centred phase supplied without that metadata is
  simulated as primitive and keeps reflections its real structure forbids, with
  nothing in the spot list to say so. The audit reports, per phase, the centering
  applied and whether it was *declared* or *assumed*; `describe()`, the
  reflection table and the manifest all carry the statement, and an assumed
  centering produces an explicit warning.

  This was not hypothetical: the shared Burgers worked-example setup declared no
  space groups, so it had been simulating beta-titanium without body-centring
  absences and listing forbidden reflections. Both phases now declare theirs, and
  a worked example pins that no `h + k + l` odd beta reflection survives.

- **`variant_correspondence_table(relationship, objects)` answers "what does this
  plane (or direction) become in every product variant?"** as a table rather than
  a loop: one row per (source, variant) carrying the exact image, its nearest
  integer indices, the angular residual between them, and an equivalence-group
  label that collapses variants giving crystallographically equivalent images.
  It takes one object or a list, runs either mapping sense, labels hexagonal
  phases in four-index Miller-Bravais form, and exports through `to_csv`,
  `to_markdown`, `to_records`, `to_json_dict` and `describe()`.

  The grouping turns 24 rows into a readable answer. Under Kurdjumov-Sachs the
  austenite `(111)` has exactly four distinct images across the 24 variants, six
  variants each; the six giving a `{011}` ferrite plane at zero residual are
  exactly the packet `variant_close_packed_groups` returns, which the tests
  assert against each other rather than against a stored number. The reverse map
  is not selective in the same way: the child `(011)` maps back onto a `{111}`
  parent plane in all 24. `exact_rows()` isolates the correspondences that are
  genuinely exact, which — unlike the grouping of the irrational images — do not
  depend on the rationalization bound.

- **`default_relationship_catalog(parent_phase, child_phase)`** resolves the
  standard named catalog from the two crystal systems through one auditable
  dispatch table, and returns `None` rather than forcing an inapplicable list
  when no standard catalog covers the pair.

- **`is_hexagonal_phase` moved into `pytex.core.hexagonal`** and is exported
  from `pytex.core`; `pytex.diffraction.composite` re-exports it. Core
  notation needs the same test the diffraction labels use, and the library
  keeps one definition of it.

- **Program specification** for transformation crystallography and composite
  diffraction:
  [`docs/architecture/transformation_crystallography_and_diffraction_program.md`](docs/architecture/transformation_crystallography_and_diffraction_program.md),
  with the running phase ledger in
  [`docs/roadmap/working_notes_transformation_diffraction_program.md`](docs/roadmap/working_notes_transformation_diffraction_program.md).

- **Parent-grain reconstruction no longer treats connectivity as proof of a
  shared parent.** This is a scientific behavior change: groupings and parent
  counts change again, and strictly for the better.

  Union-find over linked edges assumed every link was certain. It is not: two
  unrelated parents can share a boundary that genuinely lies *inside* the
  same-parent fingerprint, and no edge test can reject those because they are
  indistinguishable from real same-parent boundaries. On a dense grain graph
  each parent pair shares several boundaries, so one coincidence anywhere along
  a shared boundary merged the pair irreversibly. At map scale (100 tiled
  parents, 900 grains, ~1740 edges) that cost **30% of all parents** under
  Kurdjumov-Sachs at the default tolerance.

  Each connected cluster is now split by agreement instead: every member
  proposes the parent it implies, each proposal is scored by how many members
  it explains (``C_j^T P`` near the variant-description set), and the
  best-supported proposal claims its supporters, with unexplained members
  repeating the vote. A cluster spanning two parents separates because no
  single orientation explains all of it; a genuine single-parent cluster is
  returned whole, so the sparse-adjacency behavior is unchanged.

  On the identical map-scale sweep this turns **69.7 recovered parents into
  99.7 of 100**, with 95% of grains landing in single-parent clusters at the
  default 3 deg tolerance. At 1 deg and below the partition is recovered
  *exactly* — every grain, every parent — for both Kurdjumov-Sachs and Burgers,
  including with 0.25 deg added scatter. Pinned by a regression fixture built
  around a searched worst case: two parents 49.18 deg apart whose children
  share a boundary only 0.36 deg from the fingerprint, which the edge test
  correctly cannot reject and which consistency nonetheless separates.

- **Reconstruction results now report how reliable their own clustering is.**
  `ParentGrainReconstructionResult.chance_link_probability` is the probability
  that two *unrelated* grains would be linked at the tolerance in use — the
  fraction of uniformly random misorientations lying within `tolerance_deg` of
  the same-parent fingerprint. It depends only on the relationship and the
  tolerance, not on the data, so it is available on real maps where no ground
  truth exists. `describe()` warns when the expected number of chance links
  across the tested edges reaches one, because a single chance link merges two
  parents irreversibly:

  > WARNING: at this tolerance 7.3% of unrelated grain pairs fall within the
  > same-parent fingerprint, so about 2 of the 29 tested edge(s) are expected to
  > link by chance alone.

  This exists because the map-scale sweep showed the failure is **not** a false
  positive in the edge test — that rate stayed at zero throughout — but genuine
  physical coincidence, which no binary edge test on orientations can avoid.

- **Parent-grain reconstruction has a measured operating envelope.**
  `scripts/study_reconstruction_robustness.py` sweeps orientation noise, edge
  tolerance, and grain count against planted ground truth (48 cells, 25 seeds
  each), and
  [Reconstruction Robustness Study](docs/testing/reconstruction_robustness_study.md)
  records the result. Headline: **the false-link rate is exactly zero in every
  cell** — the edge test never merged two separable parents anywhere in the
  sweep. The remaining failure mode is the opposite one, splitting a parent by
  rejecting a noisy same-parent boundary, which yields the practical rule
  **set `tolerance_deg` to at least four times the per-grain orientation
  scatter** (at 2x the partition collapses; at 4x it is essentially always
  exact). Parent-orientation error tracks $\sigma/\sqrt{n}$ in the grain count,
  confirming the quaternion-eigen-mean refinement averages noise rather than
  inheriting it.

  A second sweep at **map scale** — 100 tiled parents, 900 grains, dense
  four-connected adjacency — showed the tolerance also has an *upper* bound,
  from a different mechanism. Because one chance link anywhere along a shared
  boundary merges two parents irreversibly, and a dense graph offers many such
  boundaries, connectivity-only clustering recovered **only 69.7 of 100 parents
  under Kurdjumov-Sachs** at the default tolerance. The false-link rate among
  separable boundaries stayed exactly zero throughout, so every one of those
  merges was a genuine physical coincidence rather than an edge-test error.
  The effect is strongly relationship-dependent: Burgers kept 97 of 100 even at
  3.0 deg, because 12 variants give a far smaller admissible set than the 24 of
  Kurdjumov-Sachs.

  That upper bound is **no longer the binding constraint** — the single-parent
  consistency split described above recovers most of the loss (69.7 to 99.7 of
  100 at the same tolerance). The default is unchanged, and the
  `chance_link_probability` diagnostic reports the residual risk.

- **An MTEX parity campaign for the orientation-relationship stack.**
  `fixtures/mtex_parity/campaigns/or_transformation_cases.json` defines shared
  cases for the OR-as-misorientation representative, the variant count, and
  recovery of the operative parent-to-child rotation from measured pairs (the
  `calcParent2Child` comparison), with a MATLAB handler
  (`scripts/mtex_generators/mtex_parity_transformation.m`) and PyTex-side
  generation. The PyTex results are generated and reproduce the literature
  values (Kurdjumov-Sachs 42.85 deg, Nishiyama-Wassermann 45.99 deg, and an
  exact Greninger-Troiano recovery from a Kurdjumov-Sachs nominal reporting the
  2.40 deg separation).

  **The MTEX side is ungenerated and the MATLAB handler is unrun**, because no
  MATLAB/MTEX installation was available where the campaign was authored. No
  PyTex document claims MTEX parity for this stack, and the validation ledger
  and generator README both state the limitation explicitly.

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

- **Generated figures were not reproducible.** `pytex.plotting.frames` derived
  its arrowhead marker ids from `builtins.hash(frame.name)`, and Python
  randomizes string hashing per process, so regenerating a committed figure
  changed its bytes on every run and `git diff` was permanently dirty. A
  generated asset whose bytes move for no reason cannot be checked for drift,
  which is the entire reason these figures are generated rather than drawn. Now
  a `zlib.crc32` digest, with a test pinning byte-for-byte reproducibility.

- **The intervariant boundary fingerprint contained duplicate elements, and its
  size moved with lattice parameters.** `intervariant_boundary_fingerprint`
  returns the deduplicated double coset `G_c (R G_p R^T) G_c` — the set of
  misorientations two children of one parent can show. Deduplication ran on
  quaternions, which need a sign convention, and the convention was "make the
  largest-magnitude component positive". Two components tie in magnitude for the
  90 and 180 degree elements of a crystal point group, so `argmax` broke the tie
  arbitrarily: numerically identical rotations canonicalized to `q` and `-q`,
  landed far apart, and were counted twice. Rounding the keys to a fixed number
  of decimals added a second, independent failure — a value near a rounding
  boundary rounds either way depending on floating-point noise.

  Together these made the size of a purely group-theoretic set depend on the
  lattice parameters that entered the rotation: the Kurdjumov-Sachs fingerprint
  returned 10664 elements for one cubic pair and 10665 for another. The true
  counts are **10 584** for Kurdjumov-Sachs and **684** for Burgers.

  Deduplication now runs on the rotation matrices, which carry no sign
  ambiguity, using a SciPy spatial query rather than a sort — duplicates are not
  reliably adjacent in lexicographic order, because a distinct element can agree
  with them in the leading entries and sort between them.

  **This is not a results change.** `boundary_fingerprint_distances_deg` takes a
  maximum over the set, so 81 duplicate elements never altered a distance,
  a reconstruction grouping, or an OR identification; they wasted 0.8% of each
  evaluation. It mattered because the set's size is quoted as a scientific
  quantity — the reconstruction robustness study reasons from it, and had
  published "about 2 800" for Burgers against "about 10 700" for
  Kurdjumov-Sachs. Both figures are corrected, and both counts plus their
  independence from lattice parameters are now pinned by tests.

- **Kinematic spot ordering was not stable against floating-point ties.** The
  sort keys are decreasing intensity, then detector radius, then lexicographic
  `hkl`. Symmetry-equivalent reflections have mathematically equal intensity and
  radius that differ in the last few ULPs depending on how the detector basis was
  built, so the *noise* decided the order before the exact `hkl` tie-break was
  ever reached. Reaching the same pattern by two equivalent routes — a parent
  zone axis, or that same axis recovered from a child variant's zone — therefore
  produced correctly-positioned but **permuted** spot tables, and any exported
  table or pinned figure inherited the permutation.

  Both continuous keys are now quantized before sorting: 1 pm of detector radius
  and 1e-12 of full-scale relative intensity, far below anything physical and far
  above the ~1e-14 noise they suppress. Ties fall through to the exact index
  comparison as intended. Row order in exported reflection tables may change for
  reflections that were previously ordered by noise.

- **Checksum-pinned fixtures no longer fail integrity checks on Windows.**
  `fixtures/phases/fe_bcc/phase.cif` held CRLF on disk and hashed to
  `e512334e…`, while its LF form hashes to `8afe4f95…` — the digest pinned in
  `fixtures/phases/catalog.json`. Git's Windows default `core.autocrlf=true`
  was rewriting checksum-pinned artifacts on checkout, so six phase fixtures
  failed `scripts/check_repo_integrity.py` and `tests/unit/test_phase_fixtures.py`
  on any such clone. A new `.gitattributes` marks `fixtures/phases/**`,
  `fixtures/mtex_parity/**` and `*.ipynb` as `-text`, disabling the conversion.

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
