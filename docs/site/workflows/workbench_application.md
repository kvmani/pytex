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

Seven workspaces holding ten panels. Every one of them ships with runnable examples, so a user with
no data of their own can still exercise every feature — the manifest test executes each example, so
an example cannot rot.

| Workspace | What it answers |
| --- | --- |
| **Crystal Viewer** | What does this structure look like, with these planes and directions drawn on it? |
| **TEM Analysis** | Everything transmission-electron, in four sub-tabs (below). |
| **XRD** | Which powder peaks should this structure produce, and how do radiation and profile choices change the diffractogram? |
| **EBSD** | What does this orientation map show — orientation, grains, local or grain-referenced misorientation — and where should it be believed? |
| **Variants** | Where do the child orientations of one parent grain point, and how do they meet? |
| **Texture** | Where does a crystal plane point across a whole polycrystal? |
| **Calculator** | Interplanar angles, symmetry families, d-spacings, orientation relationships. |

The **TEM Analysis** workspace carries four sub-tabs, because they are four views of one technique
rather than four subjects:

| Sub-tab | What it answers |
| --- | --- |
| **SAED Simulator** | What pattern would this crystal give down this axis, and where do its Kikuchi bands fall? |
| **TEM Solver** | Which zone axis is this pattern, where should I go next, and can the holder get there? |
| **CBED** | What would a convergent-beam exposure of this zone show, how thick is the foil, and what is the point group? |
| **Composite SAED** | What does a two-phase SAED pattern contain, and which variant is that spot? |

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
retaining their full titles and accessible names, and the workspace tabs wrap rather than
disappearing into an unmarked horizontal scroller. A grouped workspace's sub-tabs wrap the same
way, in their own strip below the masthead.

The Composite SAED panel adds an **Appearance** group below **Simulate pattern**. Filled or
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

## The Crystal Viewer's Orientation Dock

Turning a structure by hand answers *what does it look like from here*. It used to leave *where is
here* unanswered, which is the question a texture or a slip analysis actually needs. The dock beside
the structure answers it three ways at once, and all three describe the view you are looking at
rather than a view you have to set up separately. A fourth figure sits with them and deliberately
does *not* move with the view: the Kikuchi map, which is the atlas the other three are read
against.

**The frame, stated first, because everything below depends on it.** The screen is the specimen
frame: `RD` points right, `TD` points up, and `ND` points out of the screen towards you. That is a
right-handed triad, and it is
`pytex.core.frame_catalog.SAMPLE_RD_TD_ND_FRAME` with its default identity axes, so the
camera the viewer accumulates from your drag *is* an orientation matrix $g$ in PyTex's
crystal-to-specimen convention, $\mathbf{v}_{\text{specimen}} = g\,\mathbf{v}_{\text{crystal}}$.
Nothing new is defined for the viewer.

```{note}
Pole figures in the literature are usually drawn with RD at the *top*. This one has RD to the
right, which is the same ninety-degree difference the texture panel already carries, and it is
deliberate: the pole figure here is the same view of the same crystal as the structure beside it,
and rotating the projection would break exactly the correspondence the figure exists to show.
```

**The pole figure** plots one pole family — `{100}`, `{110}`, `{111}` on a cubic phase, the basal,
prismatic and pyramidal sets on a hexagonal one — expanded over the point group and projected
stereographically onto the upper hemisphere. Press a chip to add or remove a family; the chip
carries the colour its poles are drawn in, so the legend is the control. Any planes you superimposed
on the structure appear as labelled diamonds, because the figure should show the features the
picture shows.

**The inverse pole figure** asks the opposite question: which crystal direction currently lies along
a specimen axis. The chips choose the axis — ND by default, RD and TD available — and the point is
folded into the fundamental sector of the phase's point group, which is what makes the standard
triangle standard. Its corners are labelled with the low-index direction each is, computed rather
than assumed, so a hexagonal phase gets `[0001]`, `[2-1-10]` and `[10-10]` without anyone hard-coding
the cubic triangle.

**The fly-by** is the trail of dots behind the moving points, on *both* figures. In the inverse
pole figure it is the path one specimen axis takes; in the pole figure it is the path of a whole
family, drawn as a fading cloud rather than a single track, since a family turns rigidly and every
member of it is going somewhere at once. On the pole figure it is the plainest demonstration that
poles move on circles about the rotation axis; on the inverse figure it is the fastest way to see
that the standard triangle is a fundamental region.

The inverse figure's trail Turn the crystal and the trail draws the path
the specimen axis takes through the triangle, which is the fastest way to *see* that the standard
triangle is a fundamental region: the path jumps when the direction crosses a symmetry boundary,
because the two sides of that boundary are the same direction. Both trails are dots rather than
lines for the same reason — a line would draw a chord the crystal never took, across the triangle
in one figure and between two different poles in the other. They clear on a reset, on an axis
button, and whenever an orientation is set outright.

**The Kikuchi map** is the crystal's whole band network, projected stereographically about a zone
axis you choose — `[001]` to begin with, and editable in the figure's own caption. Each band is
drawn as its two edges, which are what a plate records, with the plane's trace dashed between them,
which is a construction rather than something visible; a band's angular width is twice its Bragg
angle, so the widest bands belong to the *finest*-spaced planes. Where bands cross, a zone axis is
marked and labelled, sized by how many bands meet there — the count that is the n-fold symmetry of
the pattern you will see on arriving.

Unlike the two figures above it, this one does not turn with the camera, and that is the point. The
pole figure and the inverse pole figure are the same view of the same crystal as the structure. The
map is the atlas: fixed to the crystal, centred where you put it, with a cross marking which
direction the current view has on the beam. Turning the crystal moves the cross across a stationary
map, which is how a map is read. It comes from `crystal.kikuchi_map`, once per phase and centre
rather than per frame, because the crystal's band network does not depend on where you happen to be
looking from.

**Setting an orientation.** The **Orientation** group in the control rail holds the angle triple,
in Bunge $(\phi_1, \Phi, \phi_2)$ or Matthies $(\alpha, \beta, \gamma)$. Type three angles and press
**Set view** and the structure turns to them. Drag the structure and the same three fields report
where you have got to, along with the rotation as an axis and angle and the crystal direction along
each specimen axis. The Bunge triple $(\phi_1, \Phi, \phi_2)$ heads the readout **whichever
convention the picker holds**, because it is the convention every EBSD file, every ODF section and
every published orientation is written in, and a reader working in Matthies still needs it in front
of them to compare with anything outside this program. The named ideal orientations below — cube, Goss, brass, copper, S and the two
rotated variants — are one press each; they are offered only for cubic phases, because that
catalogue is the rolling-texture catalogue of cubic metals and "Goss" on a hexagonal phase would
name a relationship it does not describe.

**Why the pictures are live and the numbers are not, quite.** The figures redraw at drag rate,
because the whole crystallographic content of both — the symmetry operators, the expanded families,
the fundamental sector with its bounding-plane normals — travels with the scene in the Cartesian
crystal frame, and the browser only multiplies by the camera and projects. Euler angles are a
different matter: they depend on a convention, and a convention implemented twice will eventually
disagree with itself. So the readout is a `crystal.orientation` call that fires about a fifth of a
second after the drag settles. The pictures never lag; the numbers arrive a moment later, and they
are the only Euler angles in the application.

**The layout.** On a wide window the dock is a column beside the structure, so the crystal and both
its figures are on screen together — that simultaneity is the point, and a figure that has to be
scrolled to cannot be watched. On a narrower window it becomes a strip beneath the structure and the
structure gives up part of its height to it. On a phone it starts folded away. It is a disclosure
group everywhere, so a reader who came to look at bonds can put it away entirely.

`crystal.orientation` is a service operation like any other, so the same numbers are available
without the interface:

```python
from pytex.app import REGISTRY

result = REGISTRY.call(
    "crystal.orientation",
    {
        "phase": {"builtin": "cu_fcc"},
        "euler_convention": "bunge",
        "angle1": 0.0,
        "angle2": 45.0,
        "angle3": 0.0,
    },
)
result["data"]["camera_matrix"]  # nine numbers, row-major, crystal to specimen
result["table"]["rows"]  # the crystal direction along RD, TD and ND
```

Passing `camera_matrix` instead of the angles runs the conversion the other way, which is what the
readout does. See {doc}`../examples/generated/workbench-service-layer` for the checked worked
examples: that Goss puts ND forty-five degrees from `[001]`, and that the two conversions are exact
inverses.

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

```{seealso}
{mod}`pytex.plotting.tilt_stereogram` draws the *publication* form of this figure: a two-panel
matplotlib stereogram in the **crystal** frame, with the trajectory as spaced dots, the reachable
region as exact circular arcs, and every symmetry-equivalent target marked reachable or not. The two
are complementary rather than duplicates — that one is the figure for a paper or a logbook, drawn
from a complete `TiltPlanReport`; this one is the interactive map beside the pattern, drawn in the
holder frame so that position *is* tilt, and drawn whether or not a target has been chosen yet. Both
take their geometry from the same `pytex.tem.stage` forward model.
```

### Kikuchi bands on the solved pattern

Once the pattern is **indexed**, a **Kikuchi** toggle joins *Lattice* and *Calculated* on the
pattern toolbar, and draws the bands the solution predicts as fine dotted lines under the spots. It
draws from the accepted solution when one has been accepted, and otherwise from the candidate
currently selected in the list — deciding between candidates is done by looking, and the bands are
one more thing to look at, exactly as the calculated spots are.

It is the same angular space. A detector records the directions of the outgoing electrons, and both
the spots and the bands are placed in that space by the same reciprocal lattice and the same
orientation, so superimposing them mixes nothing: a plane $(hkl)$ and its normal $\mathbf{g}$ are
one crystallographic object, and the spot at $\mathbf{g}$ and the band centre line for $(hkl)$ are
pole and polar of one another.

The metrics agree too, and that is what makes the overlay checkable. A band's width is
$L\,2\theta_B \approx \lambda L / d = r_g$, so

> the band for $(hkl)$ is exactly as wide as the $000 \rightarrow \mathbf{g}$ spot distance, and
> perpendicular to it.

At an exact zone axis the band edges therefore bisect the $000 \rightarrow \mathbf{g}$ and
$000 \rightarrow -\mathbf{g}$ segments, and every band of the zone runs through the transmitted
beam. That is also why the bands are labelled *out* where they separate and never at the pole: at
000 every band of the zone crosses, which is the most crowded and least informative point of the
figure, and where the beam marker and the picks live.

**Nothing new has to be calibrated.** The overlay lives entirely in the pattern frame and is driven
by the accepted solution's `crystal_to_pattern` plus the pixel scale that already indexed the
pattern — not the diffraction rotation $\varphi_D$, not the parity, not $\lambda$ or $L$
separately. Picking two non-collinear spots measures their azimuths on the recorded image, so
`crystal_to_pattern` is fully determined by the picks, *including* the roll. What one pattern does
not determine is pattern-to-holder, and nothing here needs it.

That is the real value of the feature rather than a technicality. "Keep the (200) band aligned and
travel along it" is an instruction in the pattern frame, so it survives exactly the calibration
nobody has, where "tilt α by +12.3°" does not. This is why microscopists navigate by bands, and why
naming a **target zone axis** in the tilt form draws the connecting band distinctly and labels it
`follow (200) toward [011]`, with the low-index waypoint zones named beside it. The band comes from
{func}`pytex.tem.path.connecting_band`, which reports *no single band connects these zones* rather
than inventing one when the two axes span no low-index plane.

```{warning}
The positions and widths are exact geometry; the **contrast is not modelled at all**. Which side of
a band is excess and which deficient, how dark one band is against another, and the HOLZ lines
crossing a zone axis are dynamical. A thin foil may show strong spots and no visible bands, because
the diffuse internal source needs thickness. And the overlay is a *prediction from the accepted
solution*, not independent evidence for it — though it is checkable, since the pattern it is drawn
on was recorded before the prediction was made.
```

### Calibrating from the image

The camera equation uses one number: how much reciprocal space one pixel spans. An image that
arrives without its recorded camera length — a printed plate, a figure from a colleague — usually
still carries something whose length is known, and **Calibrate** measures against it.

Press it, click the two ends of a length you know, and say what that length is. Two kinds of answer
are accepted, and they calibrate different things:

- **A reciprocal length** (Å⁻¹ or nm⁻¹) — a scale bar, or a reflection whose spacing you already
  know. This fixes the scale directly: the coordinate units switch to *pixels with a measured
  scale*, and the **Scale** field holds the result as `1 px = … Å⁻¹`. Type into that field instead
  if the value is already known.
- **A real length on the plate** (cm or mm). This fixes the **pixel size**, and the camera constant
  goes on doing its usual job.

Either way the answer lands in the fields the indexing already reads, so nothing downstream learns a
new concept, and the cursor readout, the fitted lattice and every d-spacing change together.

```{note}
A measured scale is not an approximation of the camera-constant route: the ratio it states is the
only quantity the camera equation ever uses, so indexing through it gives the same answer.
`tests/unit/test_app_tem_calibration.py` runs one practice plate both ways and compares them.
```

## The SAED Simulator

The forward problem, where the solver does the inverse one. Choose a phase and point the crystal —
either by **zone axis** with a roll about it, which is how one thinks at the column, or by **Bunge
Euler angles**, which is how a measured orientation arrives from EBSD or from an indexed pattern —
and the panel draws the pattern that would land on the plate at your camera length and detector
pitch.

Beside the drawing it *states* the orientation rather than implying it: the zone axis on the beam,
the same orientation in Bunge angles, the roll about the beam, and — in orientation mode — how many
degrees the orientation asked for is from the zone actually drawn. That last number is never
silent, because a spot pattern exists only on a zone axis: an orientation five degrees off [011]
has no pattern of its own, and what is drawn is [011].

**Kikuchi** superimposes the bands the same orientation predicts, fetched through the same overlay
operation the solver uses. The relation to check by eye is that the band belonging to (hkl) runs
perpendicular to the line from the transmitted beam to its own spot and is exactly as wide as that
spot is far out — both consequences of the band and the spot being pole and polar of one reciprocal
vector. Seeing it hold on a pattern whose answer is known is how one learns to check it on a plate
whose answer is not.

The plate is painted by the same drawing that paints the solver's practice patterns, so the two
look identical; they are the same calculation.

## The CBED Panel

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

For the Burgers orientation relationship, choose the built-in example in Calculator, Composite SAED
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

## Using Your Own Data Files

Three workspaces read measurements as well as models. In every case the file goes through the same
library importer a script would call — {func}`pytex.adapters.read_scan`, which dispatches to
{func}`pytex.adapters.read_ang`, {func}`pytex.adapters.read_ctf` or
{func}`pytex.adapters.read_oh5`, and {func}`pytex.adapters.read_xrdml_pole_figure` — so phases,
symmetry, grid topology and quality channels come from the file's own header rather than from
anything the application assumes. The contents travel as an ordinary request parameter — as text
for a text format, base64-encoded for an HDF5 one — and are written to a temporary path only for as
long as the reader has it open; nothing is retained between requests. See {mod}`pytex.app.uploads`.

**EBSD — `.ang`, `.ctf`, `.oh5` and `.h5`.** *Open a scan* in the EBSD rail takes an EDAX/TSL
`.ang`, an Oxford/HKL `.ctf`, or an EDAX OIM HDF5 scan under either of its two extensions, `.oh5`
and `.h5` (one format, two names). While one is open it replaces the practice dataset, and every
control means exactly what it means for a practice map: the colouring, the scalar modulation, the
boundary overlay and the grain threshold all act on the imported map. Hexagonal scans — which EDAX
writes by default — have no rectangular shape, so they are drawn on a half-step raster: every
measurement lands on its own cell, the cells between them are drawn empty rather than filled in, and
nothing is interpolated.

The HDF5 formats carry more than their text export does: every per-point scalar channel in the file
is read, not only the columns an `.ang` row has room for, so a channel a processing tool wrote back
into the scan is available to colour or modulate the map. Two practical notes. Reading them needs
`h5py`, which is optional — install PyTex with the `hdf5` extra — and the panel says so rather than
failing obscurely if it is missing. And a scan saved *with its diffraction patterns* is far larger
than a request can carry, which is a property of the file and not of the transport: export it
without the patterns, or read it with {func}`pytex.adapters.read_scan` in a script.

```{warning}
An imported map has **no known answer**, and the panel says so where a practice dataset states one.
The practice datasets are checkable because they were constructed; a measurement is only as good as
the scan behind it. Read the confidence-index or fit channel beside the orientation map before
believing a feature in it.
```

**Texture — `.xrdml` pole figures.** The **Measured pole figures** view opens one or more Panalytical
XRDML files, one per measured reflection, and draws them in **tabs** — full size, one at a time,
rather than shrunk into a row. Give the plane of each file in the order the files were opened: XRDML
records the diffraction angle, not the reflection, so the assignment cannot be read from the file.

Two controls decide whether the set can be read as a set:

- **One scale for every figure** puts them all on the same intensity range. Two pole figures of one
  specimen drawn on separate scales cannot be compared, and comparing them is the reason to measure
  more than one. Turn it off only to inspect a weak figure alone.
- **Contour levels** are a reading decision, not a property of the data. Type them — `1, 2, 4, 7, 10`
  is the sequence most of the texture literature uses — or leave the field empty for evenly spaced
  ones. Whatever is chosen applies to every figure in the set.

**Reconstruct the ODF** inverts the opened set into an orientation distribution and adds it as a
further tab, sliced at φ₂ = 0°, 45° and 65°. This is the classical inverse problem of quantitative
texture analysis and it is **ill-posed**: pole figures are projections and lose the odd-order
information, so the answer depends on the dictionary, the kernel and the regularization. One pole
figure cannot constrain it at all and three from different planes is the usual minimum. The residual
is printed beside the sections; read it before the peaks.

Normalisation decides what the numbers mean at all. A measured figure arrives in detector counts,
which depend on the counting time and the instrument; **m.r.d.** rescales it so a texture-free
specimen reads 1 everywhere, which is the only form in which two instruments' figures mean the same
thing. Defocusing and absorption corrections are *not* applied, so the outer rim of an uncorrected
figure is the part to distrust.

**TEM — any image.** Covered above: the micrograph never leaves the browser, and **Calibrate**
measures its scale from a length you know.

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
- {doc}`composite_or_diffraction` — the composite two-phase pattern the Composite SAED panel simulates
- {doc}`tem_pattern_indexing` — the full TEM Solver workflow, from a plate to the next zone axis
- {doc}`saed_pattern_solving` — the indexing the TEM Solver performs
- {doc}`crystal_visualization` — the renderer behind the Crystal Viewer's figures
- {doc}`stereographic_projections` — the projections the Variants and Texture panels use
- {doc}`../theory/pole_figure_arithmetic_and_mrd` — why the mean of a pole figure is 1 m.r.d.
