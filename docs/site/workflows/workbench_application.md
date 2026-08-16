# The PyTex Workbench: Desktop And Intranet Application

The workbench is an interactive application over the PyTex library. It ships as a desktop app and
as an intranet web app from **one codebase**: the scientific code, the service layer, and the entire
user interface are shared, and the two diverge only in how a window opens and where a saved file is
written.

The design record — including the decisions that fix the architecture — is
{doc}`../architecture/application_platform`. This page is the user guide.

## Running It

Both verbs come from the same module, and neither needs anything beyond the library's own
dependencies: the server is `http.server` from the standard library, and the frontend is
hand-written ES modules with no bundler and no third-party JavaScript.

```{code-block} console
$ python -m pytex.app desktop
```

Opens a maximized native window through `pywebview` if it is installed, and the default browser
otherwise. Maximizing by default keeps the calculation-activity bar inside the operating system's
work area; the window remains resizable down to 960 × 640. Use `--browser` to force the browser, and
`--port` to pin the loopback port.

```{code-block} console
$ python -m pytex.app serve
```

Serves the application over HTTP, on loopback by default. Pass `--host 0.0.0.0` to serve
colleagues on an intranet.

```{warning}
The application has **no authentication**. Binding it to anything other than loopback exposes every
operation to everyone who can reach the port, which is why it is not the default and has to be asked
for explicitly. Serve it only on a trusted network.
```

A third verb prints the operation manifest as JSON, which is how to script against the service layer
without a browser:

```{code-block} console
$ python -m pytex.app operations
```

## What Is In It

Nine workspaces. Every one of them ships with runnable examples, so a user with no data of their own
can still exercise every feature — the manifest test executes each example, so an example cannot rot.

| Workspace | What it answers |
| --- | --- |
| **Crystal Viewer** | What does this structure look like, with these planes and directions drawn on it? |
| **TEM Solver** | Which zone axis is this pattern, where should I go next, and can the holder get there? |
| **CBED** | What would a convergent-beam exposure of this zone show, how thick is the foil, and what is the point group? |
| **Diffraction** | What does a two-phase SAED pattern contain, and which variant is that spot? |
| **XRD** | Which powder peaks should this structure produce, and how do radiation and profile choices change the diffractogram? |
| **EBSD** | What does this orientation map show — orientation, grains, local or grain-referenced misorientation — and where should it be believed? |
| **Variants** | Where do the child orientations of one parent grain point, and how do they meet? |
| **Texture** | Where does a crystal plane point across a whole polycrystal? |
| **Calculator** | Interplanar angles, symmetry families, d-spacings, orientation relationships. |

Every panel is generated from a **self-describing operation manifest**. A capability registered in
Python appears in the interface — with its controls, its units, its help text, and its citations —
without anyone editing a menu, which is what keeps the documentation from drifting away from the
behaviour.

Press <kbd>Ctrl</kbd>+<kbd>K</kbd> to search every operation, example and material at once.

## Help, Documentation, And Appearance

The **?** beside a field explains that input in place. The **Help** button explains the active
operation, every input, the returned result and the scientific sources. Each operation also names
the closest canonical Sphinx page; **Read the full guide** opens that page for the theory,
conventions, worked examples and related APIs. The manifest test verifies that every operation has
a link and that its MyST source exists under `docs/site/`, so documentation moves cannot leave a
stale button behind.

The **About** button opens the identity panel: the running version, what the program is, the author
and their institute, contact addresses, and the licence with its warranty disclaimer. Every fact in
it is served by {mod}`pytex.app.about` on the application manifest, so the version shown is the
version that answered the request rather than a number typed into a page — and the licence shown is
the one `pyproject.toml` and `LICENSE` declare, which `tests/unit/test_app_about.py` checks against
both files.

The colour-theme control cycles through **Auto**, **Light**, and **Dark**. Auto follows the operating
system; an explicit choice is remembered by the shared frontend, so it behaves identically in the
browser and desktop window. On a narrow screen the masthead actions collapse to icons while
retaining their full titles and accessible names, and the nine workspace tabs wrap rather than
disappearing into an unmarked horizontal scroller.

The Diffraction workspace adds an **Appearance** group below **Simulate pattern**. Filled or
unfilled circle, square, triangle, diamond, star and cross symbols, overall spot-size scale,
intensity-to-size mapping, parent colour and variant palette redraw the current result immediately.
Variant encoding can use shape, size, colour or any combination of those channels. The default
shape + size + colour mode gives every one of the 24 product variants a distinct symbol, and each
legend chip mirrors the exact marker shown in the pattern.
These are presentation choices: changing them does not run the simulation again and cannot alter
detector coordinates, indices, intensities, hover data or exports. **Perceptual intensity** keeps
weak reflections legible, **Area follows intensity** makes marker area quantitative, and
**Constant size** is useful when comparing geometry alone. The colourblind-safe cycle is the
preferred palette when colour must distinguish variants. Double-diffraction spots keep their
dashed edge when all spots are unfilled, preserving their scientific designation. The transmitted
beam is always rendered independently as a double-ring centre mark, labelled **(000) transmitted**
in the plot and **Transmitted beam (000)** in the legend.

The legend is also a visibility controller. Click any source chip to hide or restore that source;
its pressed state remains visible and keyboard focus is preserved. **Show all** restores the full
composite, **Parent only** removes every product lattice, and **Focus a variantâ€¦** retains the
parent reference plus one selected variant. Visibility is presentation-only too: hidden rows remain
in the result table and every numerical export.

The **XRD** workspace uses the same canonical phase picker to generate structure-aware powder
patterns. It exposes Cu, Mo and Co radiation, single-line or Kα1/Kα2 doublets, tabulated X-ray form
factors, Gaussian or pseudo-Voigt broadening, scan range, angular sampling and index limit. The
result retains both the sampled profile and an indexed Kα1 reflection table with $2\theta$, $d$,
multiplicity, $|F|$ and Lorentz–polarization factor. Hovering a peak reads that exported row.
Profile and reflection-stick colours, line width, area fill, labels and linear, square-root or
log-like display scaling update immediately without recomputing or altering the scientific result.
Four built-in standards demonstrate fcc indexing and Cu doublet splitting, silicon extinctions,
the wavelength shift under Mo radiation and the hcp zirconium metric.

## Reading A Result

Every result carries four things, and the fourth is the one that makes the application usable as an
instrument rather than a picture generator.

1. **A prose summary** that states what was computed, with its conventions explicit. This is the
   explainable-results doctrine of `AGENTS.md` at the application surface: a result that cannot say
   in words what it means is not finished.
2. **A table**, one row per reported entity, at full precision.
3. **Notes and citations**, including the caveats that apply to the specific settings used — the
   kinematic-intensity warning appears on a monatomic phase because that is where it matters most.
4. **A live cursor readout and hover detail on every drawn entity.** Hovering a diffraction spot
   gives its indices, spacing, intensity, phase and variant; hovering a pole gives its variant,
   packet and Miller indices. The payload is the same row the CSV export writes, so what the screen
   says and what the file contains cannot disagree.

Every two-dimensional graphic uses the same inspection gestures: scroll over the graphic to zoom
about the pointer, Shift-drag (or middle-drag) to pan, and press **Fit** to restore the full figure.
The percentage in the plot toolbar is the current magnification. The Crystal Viewer uses scroll to
zoom and ordinary drag to rotate its three-dimensional camera, so it keeps its own Reset control.

The Crystal Viewer's **Object properties** section separates structure from presentation. It can
hide object classes; scale and fade atoms; recolour each species; set bond and cell edge weight;
adjust plane and direction colour/opacity; scale annotations; show or hide the axis gizmo; and tune
the studio lighting and depth cue. Atoms are shaded as spheres and bonds as cylinders by default;
**Matte** softens the highlight and **Flat** gives an intentionally diagrammatic disc-and-line view.
Azimuth and elevation move the light in screen space, while ambient, diffuse, specular, highlight
sharpness, and depth cue controls expose the professional rendering model without mixing it into
the crystallographic inputs. A change redraws the existing scene immediately—no crystallographic
service is called, and no atom, plane or index moves. Publishing the figure sends the validated
appearance object with the camera; the browser converts its screen-space light into the crystal
frame so the Python renderer uses the same camera-relative illumination, visible object classes,
colours, opacities, and relative sizes.

```{note}
The on-screen table is a **preview**: the first 200 rows, with the caption saying so. Every export
carries all of them. A texture ODF is over a thousand rows, and a scroll box with no search is not
how anyone reads a thousand numbers.
```

## The Message Log

One console is docked at the bottom of the window, and every part of the application reports into
it: the panels when a spot is picked or a control is rejected, the service layer when a calculation
starts and finishes, and Python's own `logging` — the last of which a desktop user could otherwise
never see, since a native window has no terminal behind it.

It exists because the same session's story used to be split four ways: a toast that vanished after
nine seconds, a strip that counted calculations without describing them, an error line beside a
control, and a terminal. There was no way to answer "what did the application just do, and in what
order".

Each record carries a severity, a timestamp, the surface that reported it, and structured detail:

| Level | For |
| --- | --- |
| **Progress** | A long task's percentage and estimated time remaining, e.g. `50% progress. ETA: 2 min 30 sec`. |
| **Debug** | Detail only a maintainer wants. |
| **Info** | Ordinary narration: a workspace opened, a calculation started. |
| **Important** | A result worth noticing — the number of grains found, the disc regime, where a file was saved. |
| **Success** | Something completed, with how long it took. |
| **Warning** | Something suspect that did not stop the work, such as a malformed entry in a control. |
| **Error** | Work that failed, quoting the same sentence shown beside the offending control. |
| **Critical** | A failure that leaves the application unusable, such as the server becoming unreachable. |

Collapsed, the bar shows the newest message and a count of unseen warnings and errors. Opened, the
stream filters by severity and by text, and copies as plain text for pasting into an issue. A long
task's progress ticks collapse onto a single line that counts up, rather than burying everything
else under a thousand entries.

Two delivery paths, deliberately. A call's own records ride back **on its response envelope**, so a
message about a calculation can never arrive before the result it describes; everything else —
server start-up, background work — is polled from `GET /api/log?since=<sequence>`. Merging on the
sequence number is what keeps a record delivered by both paths from appearing twice.

The rule for what belongs here: a record is a finished sentence, addressed to the researcher, with
the numbers in it. `Spot 2 is selected: coordinates are (452.48, 709.30) px.` — not
`handleClick fired`.

## Exporting

Two kinds of export, and the distinction is deliberate.

**A publication figure.** The Figure button on a plotting panel renders the same scene through
`pytex.plotting` in the journal style. The interactive view and the figure are the same geometry,
because the camera and the settings travel with the request.

Choose the format for the drawing, not by habit: a lit sphere is a mesh and a vector format writes
every facet of it, so PNG is right for the crystal viewer; a pole figure is line art, so SVG is
right there and is roughly twenty times smaller.

**Re-plottable numbers, and a readable account of them.** Four formats, on every result in every
panel. The CSV is one row per entity at full precision; the XLSX adds a sheet recording the inputs;
the JSON is the complete result, reloadable. The **Report** is a Markdown page written for a person
rather than a program: the answer in prose, the caveats, the data, the exact inputs that produced
it, and the citations — the thing to paste into a notebook entry, which none of the other three is.
No result in the application is exportable only as a picture, and a result with no table at all
still exports as a report, because the prose and the provenance are the point.

The buttons are generated from the manifest, like everything else here, so a format added in Python
appears on every result at once.

```{note}
In the desktop shell every export goes through a native save dialog, because an embedded web view
silently discards a browser download. Both shells announce where the file went. An export that
passes in silence is the failure this path exists to prevent.
```

## Worked Reading: Two Panels With Known Answers

Two of the panels can be checked against answers that are fixed before the calculation runs, which
makes them the right places to start.

### The variants panel, against Morito's table

Open **Variants** and run the first example. One austenite grain under Kurdjumov-Sachs produces
**24 variants**, and colouring by packet shows them falling into **4 groups of 6** — one per member
of the parent {111} family, because each variant carries exactly one {111} into exact parallelism
with a child {110}. That grouping is the packet a lath martensite micrograph shows as a block.

Switch the view to the misorientation spectrum. The 276 variant pairs collapse onto **ten discrete
disorientations**: 10.53°, 14.88°, 20.61°, 21.06°, 47.11°, 49.47°, 50.51°, 51.73°, 57.21° and
60.00°. Those are the values of Morito *et al.*, Acta Mater. **51** (2003) 1789, Table 2. Of the ten,
only three occur *within* a packet — 10.53°, 49.47° and 60° — and those three are the ones with an
exactly rational axis, ⟨111⟩ and ⟨110⟩. The 60°/⟨111⟩ pair is the Σ3 twin relation.

This is the panel's whole argument: a measured misorientation histogram from prior-austenite grains
should show peaks at those angles and nowhere else, and a peak elsewhere is a boundary between two
*different* parent grains.

### The texture panel, against the m.r.d. scale

Open **Texture** and run the random baseline first. It reads **1 m.r.d. everywhere**, and the status
line reports the area-weighted mean as 1.000. That is not a coincidence to admire but the definition
of the scale: multiples of a random distribution means exactly that a texture-free material is 1, and
the area-weighted mean of any correctly normalised pole figure is 1 whatever the texture. A figure
whose mean is not 1 has not been normalised, and its numbers mean nothing outside itself.

Then choose the Goss component and plot (011). Goss is written {011}⟨100⟩ — the {011} plane lies in
the sheet plane — so its (011) pole must sit at the **centre** of the figure, which is ND. It does.
The notation is a testable claim about where the poles go, and checking it needs no reference
figure.

### Contour presentation

Pole figures and ODF sections can be shown as **contour lines**, **filled contours**, or both. The
**Contour properties** section controls the automatic number of levels or an exact positive level
list (for example `0.5, 1, 2, 4` m.r.d.), the upper colour limit, palette, isoline colour/weight,
fill opacity, and pole-figure display-grid resolution. The default m.r.d. palette keeps 1 m.r.d.
as the neutral random baseline; Viridis and Turbo are available for presentation needs.

Pole-figure values arrive on a spherical sampling support, so the browser inverse-distance
interpolates them onto a clipped projection grid for display. Marching squares draws isolines and
the same grid is quantised into filled bands, guaranteeing that line and fill boundaries use the
same declared levels. This interpolation is presentation-only: hover targets remain the computed
samples, and CSV/XLSX/JSON/Report exports retain the service values. **SVG** in the plot toolbar
saves the current line/filled appearance through the same saver used by both app shells.

### Calculation progress and history

Every Python service call narrates itself into the message console described above: one record when
it starts, one when it finishes with its elapsed time, and an error record quoting the same sentence
shown beside the offending control when it fails. While work is in flight the collapsed bar names
the operation and shows an indeterminate indicator.

The console is an interface aid, not a scientific record. Reproducible parameters and provenance
remain in result exports and reports; the log is bounded and is discarded when the process ends.

### Browser verification

The repository exercises the shared frontend in Chromium with `npm run test:browser`. The critical
suite opens every workspace, completes each panel's default scientific calculation (including the
TEM auto-pick/index path), verifies the visible error envelope and the message console, rejects page
and console errors, and checks that all nine workspaces remain reachable at 390 × 844. These are
real requests to a loopback PyTex server. Playwright exists only in the development lockfile; the
application delivered to users remains hand-written ES modules with no third-party browser code.

```{warning}
The kernel halfwidth is a **smoothing choice, not a property of the material**. Too small and the
figure shows the individual grains that were sampled; too large and real detail is washed out. It is
the setting most often left unreported in the literature, which is why the application puts it in the
control rail rather than burying it in a default.
```

## The TEM Solver, Step By Step

This panel is laid out as the microscope session is, and its four numbered sections are the four
questions in order.

**1 · Open a pattern.** Three practice plates ship with the application — aluminium down [001],
ferrite down [110], zirconium down [2̄110] — each a real crystallographic calculation with the
answer attached, so the indexing workflow can be tried and *checked* without a micrograph. Or open
your own image: it never leaves the machine, and only the coordinates you click are sent. The
accelerating voltage and camera length are microscope settings shared by all three plates, and the
camera constant is computed from them rather than typed.

**2 · Calibrate and index.** Click the transmitted beam first — it is the origin every spot is
measured from, not a reflection — then the spots. **Auto-pick** places the beam and six strong,
mutually non-collinear reflections, which is what indexing needs: the brightest spots of a pattern
almost always include a Friedel pair, and a pair collinear through the beam cannot seed anything.
**Show answer** labels every simulated spot with its indices, for checking your reading against the
construction. When a practice plate is loaded, the indexed result is compared with the axis it was
built from — up to symmetry, because a bcc [110] pattern is indistinguishable from a [101] one.

As soon as two spots are picked, a **lattice fitted to them** is drawn over the pattern, with the
two basis vectors as labelled arrows from the beam to the picks that generate them. Move a spot an
arrow points at and every line in the grid turns with it, so the picks worth being careful about
are obvious. The **beam-centre tool** beside the form reports where the beam is and where the spots
say it should be; nudge it with the pad and watch the grid, or press *Refine from the spots* and
let least squares place it. That number matters more than it looks: the camera equation measures
every radius from the beam, so an error there biases every spacing at once and yields a
self-consistent answer for the wrong material.

Indexing then lists **candidates ranked by an accuracy score**, each with three bars — how well
d-spacings agree, how well angles agree, and how many spots were explained. Selecting one draws
**the pattern it predicts** over the one you measured, which is the honest way to choose between
them; *Accept this solution* is a separate act, and it is what carries the phase and axis into the
two steps below. The *Where it disagrees* card gives the evidence: the same deviation on every spot
is a camera constant, a scatter of them is an indexing error, and angles do not depend on the
calibration at all.

**3 · Where to go next.** Every zone-axis family of the phase near the one on the beam, with its
angle, its number of symmetry-equivalent members, how many reflections its pattern shows, the
rotational symmetry you will recognise on arrival, and whether the holder can reach it. The answer
is usually not the nearest axis: a two-spot zone 11° away buys less than a thirty-six-spot six-fold
zone at 35°. Choosing a row sets it as the tilt target below, so no indices are retyped.

**4 · Tilt to the target.** The plan, the alpha and beta angles, the margin against the tilt
envelope, and the tilt map showing at a glance which candidates sit comfortably inside the holder's
range and which are pressed against a stop.

### The stereogram beside the pattern

The stage carries two figures. The pattern says what is on the beam; the **stereogram** to its right
says what else is within reach and in which direction it lies. It follows the tilt step's inputs
directly — change the axis, the target, the stage reading or the holder limits and it redraws — so
it is a view of that form rather than a second set of controls.

It is drawn in **holder coordinates**, not crystal coordinates: the centre is the holder's zero-tilt
axis, α increases upwards and β to the left, so a pole's position on the drawing *is* the tilt that
reaches it. Four things are on it.

- **Every zone axis of the phase** up to the index limit, projected onto the upper hemisphere.
  Poles within the labelling limit are drawn large and named; the rest are small ticks. A pole
  filled with the accent colour can be brought onto the beam within the holder's range; a faint one
  cannot.
- **The holder envelope**, drawn as the region of poles the stage can actually reach — the image of
  the α/β range itself, not a circle approximating it.
- **The axis on the beam**, marked with a crosshair where the stage puts it, which is the centre
  only at zero tilt. After *Accept this solution* that is the indexed axis, so the drawing shows the
  solved orientation rather than an assumed one.
- **The route to the target**, as a dotted geodesic — the same great circle as the connecting
  Kikuchi band — with the low-index zones lying along it ringed as waypoints. Re-indexing at each
  waypoint is what keeps a long tilt from accumulating rotation error, so the intermediate zones are
  part of the plan rather than scenery.

Hovering anywhere on the drawing reports the stage reading that would bring **that point** onto the
beam: α, β, the angle from the holder axis, and the name of the pole if the cursor is on one.
Hovering a pole itself gives its full row, including Δα and Δβ from where the stage is now.

```{note}
The stage reading beside a pole is the *principal branch*: three other branches reach the same pole
and a real holder usually cannot set them. And where a pole sits on the drawing depends on the
rotation about the beam, which one indexed pattern does not determine — the angles between poles do
not depend on it at all.
```

## The CBED Workspace

Three views of one technique, chosen with the **View** picker.

**CBED pattern** simulates a zone-axis exposure. Every zeroth-Laue-zone reflection becomes a disc
whose angular radius is the convergence semi-angle, and the disc is filled with the diffracted
intensity at the excitation error of every incident direction in the illumination cone. Position
inside a disc is therefore *incident-beam direction*, which is what makes its fringes readable.

The convergence semi-angle is the parameter that matters, and the panel reports what it produced:

- **Kossel–Moellenstedt**, discs separated. Each disc is an independent rocking curve, and its
  fringes measure the foil thickness. This is the regime a thickness measurement requires.
- **Kossel**, discs overlapping. A disc is no longer an independent rocking curve. The overlaps are
  drawn as the sum of the contributing discs, and the result says so — the *interference* between
  the overlapping beams, which is what an experiment records there, is not modelled.

Whether the discs separate depends on the material and the zone, not on the instrument alone: the
disc diameter `2α/λ` has to fall below the closest reciprocal-lattice spacing of that zone. Silicon
down [001] separates at 3 mrad and overlaps at 6.

**Two-beam** computes each disc independently — cheap, and exactly the model the thickness
measurement inverts. It is also symmetric in the excitation error by construction, so a two-beam
pattern displays a symmetry belonging to the *method* rather than to the crystal; symmetry
determination is refused on one, rather than answered with a caveat. **Bloch wave** solves the
coupled many-beam problem, and is the only method whose relative intensities and symmetry mean
anything — and the only one that can decide whether the crystal has a centre of symmetry, which
kinematic diffraction cannot determine at all.

**Thickness from CBED fringes** inverts the two-beam minima. `(s_n/n)²` is linear in `1/n²` with
intercept `t⁻²` and slope `−ξ_g⁻²`, so one least-squares fit returns the thickness *and* the
extinction distance. The extinction distance is the check: compare it with the value the structure
predicts, because a large disagreement means the fringe orders were misassigned, which is this
measurement's usual failure. The panel draws the fit as the straight line it is, so a misassignment
shows as a point off the line rather than only as a wrong number.

**HOLZ ring radii** reports where the higher-order Laue zones fall, and the reciprocal-lattice layer
spacing `H` behind them. This is the one dimension a zone-axis pattern is blind to — every
zeroth-zone reflection is perpendicular to the zone axis — and a change in `H` is a change in the
lattice parameter along the beam, which is how CBED measures local strain and composition.

## The EBSD Workspace

One map, from four independent choices. Declaring them as parameters rather than as a dozen separate
buttons is what makes a combination like *grain boundaries superimposed on a GROD map, greyed by
confidence index* reachable — nobody would have enumerated it, and it is what a real analysis asks
for.

**Colour by** decides what a pixel's colour means. *IPF* is the standard orientation map: each pixel
is coloured by which crystal direction lies along the specimen axis you choose, folded into the
symmetry fundamental sector so symmetrically equivalent orientations share a colour. The colour key
belongs to the point group, so two maps of different symmetries are not colour-comparable, and one
direction alone does not fix an orientation — X, Y and Z together do. *Grain identity* colours each
segmented grain arbitrarily. *GROD* and *KAM* render the two misorientation fields, and the
remaining choices render the measured channels.

**Modulate by** darkens any colouring by a scalar channel while keeping its hue. This is how an IPF
map is made to show *where the indexing should be believed* without giving up the orientation it is
showing: at a boundary the diffraction pattern overlaps two lattices, so the confidence index falls
and those pixels recede. Fit is an error rather than a quality, so its scale is inverted — without
that, modulating by fit would darken exactly the pixels indexed best.

**Grain boundaries** superimposes the network on whatever is underneath, classified into low- and
high-angle. They are drawn as line geometry rather than as pixels, so they stay sharp under zoom;
each segment is one pixel face, reconstructed onto the face itself rather than smoothed through the
midpoints, which would round off the corners where three grains meet.

**Grain threshold** is what counts as one grain, and it changes the grain table, GROD and the
boundary network together. The KAM threshold is separate and does a different job: it excludes
neighbour pairs above it from the average, which is the standard way to keep grain boundaries out of
an intragranular KAM.

### GROD and KAM are not the same picture

The `bicrystal_gradient` dataset exists to make this concrete. Its right-hand grain carries a linear
orientation gradient, and the same microstructure looks completely different in the two maps:

- **KAM is flat** across the whole gradient, at *half* the per-step rotation. Half, because the
  four-neighbour kernel averages two neighbours a full step along the gradient with two across it
  that are identical to the centre point. A KAM is an average over a kernel, not a gradient
  magnitude.
- **GROD is strongly graded**, and is a *deviation* rather than a ramp: it falls to zero at the
  grain's own reference orientation and rises on both sides of it.

### The practice datasets

Each is a construction whose answer is fixed before the calculation runs, and the answer travels
with the result rather than living only in the documentation:

| Dataset | Known by construction |
| --- | --- |
| **Bicrystal with a deformation gradient** | Boundary misorientation exactly 40° about [001]; KAM exactly half the per-step rotation; GROD reaching 7°. |
| **Annealing twins in a cubic metal** | Every boundary segment a Σ3 twin at exactly 60° about [111], so the misorientation histogram is a single spike. |
| **Equiaxed polycrystal** | Twelve grains, so a segmentation below the smallest boundary misorientation must find twelve. |

They are constructions, not measurements: no detector geometry, no indexing step, and no noise model
beyond the stated spread.

## Using Your Own Material

Every phase control offers the built-in catalogue — NaCl, austenite, ferrite, beta-bcc and
alpha-hcp zirconium, nickel and the rest, all with cited parameters and full atomic bases — or six
cell parameters and a point group typed in directly. A CIF can be loaded through `Phase.from_cif`
in the library, but the
application never requires it: the catalogue is defined by literal parameters in Python, so the
starting materials cannot be broken by an optional dependency being absent.

For the Burgers orientation relationship, choose the built-in example in Calculator, Diffraction
or Variants. All three use the same explicit `Zirconium (bcc, beta at 863 °C)` parent and
`Zirconium (hcp, alpha)` child. The beta phase uses `Im-3m`, `a = 3.6090 Å`, and a two-Zr-atom
conventional basis; the examples therefore exercise zirconium scattering, spacings and provenance,
not merely the cubic symmetry of an iron stand-in.

```{note}
An edited phase is renamed "(edited)" and its space group is cleared when the crystal system changes.
Both matter: an edited cell that kept the original name would put the wrong material on every title
and every export, and a space group from another crystal system would apply centring conditions that
delete reflection families the edited cell actually has.
```

## Scripting The Same Operations

The service layer is JSON-in, JSON-out and knows nothing about HTTP, so every operation the interface
offers is callable directly:

```{code-block} python
from pytex.app import REGISTRY

result = REGISTRY.call(
    "variants.intervariant_misorientations",
    {
        "phase": {"builtin": "austenite_fcc"},
        "child_phase": {"builtin": "fe_bcc"},
        "relationship": "kurdjumov_sachs",
        "packet_plane": [1, 1, 1],
        "merge_equal_angles": True,
    },
)
print(result["summary"])
for row in result["table"]["rows"]:
    print(row["angle_deg"], row["pairs"])
```

`REGISTRY.manifest()` describes every operation, parameter, default and example, which is the same
description the interface builds itself from.

## Where It Runs

The application is designed for an offline intranet host and has **zero mandatory runtime
dependencies** beyond the library: standard-library `http.server`, hand-written ES modules, no
bundler, no CDN. Two things are optional and degrade gracefully:

- **`pywebview`** — without it, `desktop` opens the default browser instead of a native window.
- **`matplotlib`** — without it, every calculation still runs and every data export still works;
  only the publication-figure buttons report a missing dependency, naming it.

## See Also

- {doc}`../architecture/application_platform` — the design record and the six decisions
- {doc}`composite_or_diffraction` — the composite two-phase pattern the Diffraction panel simulates
- {doc}`tem_pattern_indexing` — the full TEM Solver workflow, from a plate to the next zone axis
- {doc}`saed_pattern_solving` — the indexing the TEM Solver performs
- {doc}`crystal_visualization` — the renderer behind the Crystal Viewer's figures
- {doc}`stereographic_projections` — the projections the Variants and Texture panels use
- {doc}`../theory/pole_figure_arithmetic_and_mrd` — why the mean of a pole figure is 1 m.r.d.
