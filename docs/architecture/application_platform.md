# Application Platform: One Codebase, Two Shells

PyTex is a library. This document specifies the *application* built on top of it: an interactive
workbench that a researcher runs on their own machine as a desktop app, and that the same group
runs on one intranet host so colleagues can use it in a browser with nothing installed.

The governing requirement is that these are **not two applications**. They are one application with
two launchers. Everything a user can do — every calculation, every plot, every export — is
implemented exactly once, in Python, in the library-facing service layer described below. The
desktop shell and the web shell differ only in how a window appears and where files are written.

## The Layering

```text
┌──────────────────────────────────────────────────────────────┐
│  Shells (thin, the only place the two targets differ)        │
│    pytex.app.desktop   native window, local filesystem        │
│    pytex.app.server    intranet HTTP, download responses      │
├──────────────────────────────────────────────────────────────┤
│  Frontend (one copy, byte-identical in both shells)           │
│    pytex/app/static/   hand-written ES modules, no build step │
├──────────────────────────────────────────────────────────────┤
│  Operation manifest + JSON envelope                           │
│    pytex.app.registry  self-describing operations             │
│    pytex.app.contracts request/response encoding              │
├──────────────────────────────────────────────────────────────┤
│  Services (pure Python, no web or GUI imports)                │
│    pytex.app.services.*  every user-visible capability        │
├──────────────────────────────────────────────────────────────┤
│  PyTex library (core, diffraction, tem, texture, plotting)    │
└──────────────────────────────────────────────────────────────┘
```

The service layer is the product. The shells are plumbing. A capability that exists only inside a
shell is a defect, because it is a capability the other shell does not have.

## Decision 0 — Standalone first, portal integration second

The PyTex desktop and intranet web workbench are complete standalone deployment modes. Integration
with `ml_server` is optional and must use a stable service boundary; it must not require the portal
or any unrelated companion repository to be installed. A colleague who needs only PyTex must be
able to install, test, run, and update PyTex independently.

## Decision 1 — The Shared Layer Is A Service Layer, Not A Widget Library

The obvious way to share code between a desktop and a web app is to share *widgets*, through a
cross-platform GUI toolkit. We do not, because the sharing boundary would then run through the
presentation layer, where the two targets genuinely differ (window management, file dialogs, DPI,
event loops), and would not run through the scientific layer, where they do not differ at all.

Instead the boundary is drawn at **JSON in, JSON out**. A service takes a plain dictionary, does
crystallography, and returns a plain dictionary plus, where relevant, an SVG string. It imports
nothing from any GUI or web framework. This has three consequences worth stating explicitly:

1. Services are testable in the base lane with `pytest` alone — no browser, no display, no event
   loop, no screenshot comparison. The application's correctness is checked the same way the
   library's correctness is checked.
2. Services are usable from a notebook, a script, or the CLI, not only from the UI. The UI is one
   client of the service layer, not its owner.
3. Adding a third shell later (a JupyterLab panel, a batch runner, a REST client in another
   language) costs a launcher, not a rewrite.

## Decision 2 — The UI Is Generated From A Self-Describing Manifest

The requirement to "expose every functionality in an intuitive, non-cluttered manner with
sufficient inline help and systematic feature discovery" is not achievable by hand-maintaining
forms and help text alongside the code they describe: they drift apart within weeks, and the drift
is invisible until a user is misled by a stale tooltip.

So each operation declares itself. `OperationSpec` carries the operation's title, one-line summary,
long-form help, its parameters (each with a label, units, default, constraints, and its own help
string), what it returns, and its literature citations. The frontend requests
`GET /api/manifest` and *builds the controls from that*: label, input widget, unit suffix, default
value, validation, and the inline help popover all come from the same declaration the server
validates against.

The properties this buys are exactly the ones the requirement asks for:

- **No drift.** Help text lives beside the implementation and is served from it.
- **A path to depth.** Every operation also carries the closest canonical Sphinx document. The
  help drawer keeps the concise answer in place and offers the full theory, conventions, worked
  examples, and related APIs without making the control rail read like a manual.
- **Feature discovery.** A global command palette and a searchable capability index are generated
  from the manifest, so a new operation is discoverable the moment it is registered — nobody has to
  remember to add it to a menu.
- **Uniformity.** Every numeric field behaves the same way because there is one control renderer.

Hand-written panels remain for the three flagship surfaces (crystal viewer, TEM solver,
calculator), because their interactions are genuinely bespoke. Even there, the *parameters* come
from the manifest; only the layout and the canvas are custom.

Documentation targets are stored as source-relative Sphinx document names and checked against
`docs/site/` in the manifest tests, so a renamed page fails the base lane instead of leaving a
quietly broken Help link. Until PyTex has a dedicated hosted Sphinx URL, the manifest opens the
GitHub-rendered canonical MyST source; an intranet deployment may mirror that source alongside the
application without changing the operation declarations.

## Decision 3 — Zero Mandatory Runtime Dependencies For The Shells

The server is `http.server` from the standard library, not FastAPI or Flask. The frontend is
hand-written ES modules served as static files, with no bundler, no `npm install`, and no
third-party JavaScript.

This is a deliberate trade of developer convenience for deployment reality:

- The stated deployment is an **intranet host, used by colleagues without installation**. Such
  hosts are frequently offline, air-gapped, or behind a proxy that blocks PyPI and npm. `python -m
  pytex.app serve` must work on a machine that has PyTex and nothing else.
- The API surface is a dozen JSON POST endpoints and a static file tree. FastAPI would contribute
  routing sugar and OpenAPI generation; we already generate a richer, domain-specific manifest
  (Decision 2), and the routing table fits on one screen.
- No third-party JavaScript means no supply-chain surface in a tool that researchers will point at
  unpublished data, and no build step that can rot between sessions.
- `AGENTS.md` forbids silent import-time coupling to heavy stacks. A web framework in the import
  path of a crystallography library is exactly that.

The cost is that we write our own request routing (~200 lines), our own reactive rendering
(~150 lines), and our own `.xlsx` writer (~150 lines, since `.xlsx` is a zip of XML). Each is small,
auditable, and stable. Performance is a non-issue: `ThreadingHTTPServer` with a handful of
concurrent researchers is bounded by NumPy, not by the socket layer.

**Where this decision would be revisited.** If the app ever needs authentication, per-user
persistence, or WebSocket streaming of long computations, the stdlib server stops being the right
size and a real ASGI framework earns its dependency. That threshold is recorded here so the choice
is re-examined on evidence rather than defended out of habit.

## Decision 4 — Interactive Geometry Is Projected In The Browser, Published From Python

The crystal viewer must rotate smoothly under the mouse *and* export a publication-quality figure.
These pull in opposite directions: smooth rotation wants the projection in the browser, and
publication quality wants matplotlib and `pytex.plotting` in Python.

The resolution is to split by *what changes*:

- The **scene** — atoms, bonds, cell edges, plane polygons, direction arrows, labels, all in
  three-dimensional crystal coordinates — is built once in Python by
  `pytex.plotting.crystal3d.build_crystal_scene` and sent to the browser as JSON.
- The **camera** — a rotation matrix, a projection, a zoom — lives in the browser. Dragging applies
  a rotation to vertices already in hand and redraws SVG. No round trip, no server load, no lag.
- **Presentation** — sphere and cylinder shading, surface finish, light direction and strengths,
  and depth cue — is also evaluated in the browser for immediate feedback. These settings change
  only how the projected geometry is painted. Reusable SVG gradients keep the interaction light,
  and an explicit flat mode preserves a diagrammatic alternative.
- **Export** sends the camera back to Python, which renders the identical scene through the
  existing publication renderer at the requested size, DPI, and format. The screen-space light is
  transformed through the camera rotation before export so highlights remain camera-aware rather
  than silently reverting to a fixed laboratory axis.

The duplicated mathematics is a 3×3 rotation and an orthographic divide — under fifty lines of
JavaScript — and it is duplication of *viewing*, not of crystallography. Nothing about lattices,
symmetry, or Miller indices is ever computed in the browser. That line is absolute: the browser
does arithmetic on numbers Python has already made meaningful.

## Decision 5 — The Desktop Shell Is The Web Shell In A Window

`pytex.app.desktop` starts the same server bound to loopback on an ephemeral port, then opens it in
a window. It prefers `pywebview` (an optional extra) for a native, chrome-free window; when
`pywebview` is absent it opens the user's default browser at the loopback URL instead. Either way
the running application is byte-identical, and the desktop path also has **no mandatory
dependency** — `python -m pytex.app desktop` works on a bare install.

Two things genuinely differ, and are the only permitted divergences:

| Concern | Desktop | Web |
| --- | --- | --- |
| Export destination | writes to a chosen local path | streams a download response |
| Data source | may read local files by path | uploads only |

Both differences are confined to `pytex.app.shell_capabilities`, which reports to the frontend
which mode it is in. Feature code never branches on it; only file-handling code does.

## Decision 6 — Every Result Is Exportable In A Re-Plottable Form

A figure a researcher cannot re-plot is a dead end. Every operation that produces a visual also
returns the numbers behind it, and the export surface offers, per result:

| Format | Purpose |
| --- | --- |
| SVG, PDF | publication figure, vector, embedded fonts per the style guide |
| PNG | at a requested DPI, for slides and drafts |
| CSV | the tabular result, one row per plotted entity |
| XLSX | the same table plus a metadata sheet recording inputs and provenance |
| JSON | the full result object, schema-tagged, round-trippable back into the app |

The JSON export is the contract that matters most: it carries the `ProvenanceRecord` and the exact
input parameters, so a figure in a paper can be regenerated from the file that produced it. This is
the `describe()`/explainable-results doctrine of `AGENTS.md` applied to the application surface.

Presentation-only controls are deliberately outside that scientific contract. Marker shape,
visual scale and display palette operate on the rows already returned by the service and trigger a
frontend redraw, not a second scientific request. Their shared implementation lives in
`static/js/core/visualstyle.js`; operation inputs and exported numerical values remain manifest- and
service-owned. This boundary lets users prepare legible displays without implying that cosmetic
choices changed a calculation.

## Decision 7 — Every Tab Ships With Canonical Worked Examples

A researcher opening the TEM solver for the first time may have no pattern to hand, and a
researcher opening the crystal viewer may not yet know which phase makes the feature obvious. An
empty panel with a file-upload box teaches nothing.

So **every tab ships with three or four ready-to-run canonical examples**, reachable from a "Try an
example" control that is present before any input is given. An example is not a screenshot: it is a
complete, named set of inputs that runs the real service and lands the user in a working state they
can then modify.

The canonical materials are shared across tabs so a user builds familiarity rather than meeting a
new structure in every panel:

| Material | Why it is in the set |
| --- | --- |
| NaCl (halite) | Two species, F-centred, the textbook case for structure-factor absences |
| Fe-fcc (austenite) | The parent phase of every orientation-relationship example |
| Fe-bcc (ferrite) | The product phase; together with austenite it demonstrates cross-phase work |
| Zr-bcc (beta, 863 °C) | The physically correct parent for the canonical Burgers beta-to-alpha example |
| Zr-hcp | Hexagonal, so four-index notation and c/a effects appear rather than being described |

Each of these is in the built-in phase catalogue (`pytex.app.phases.BUILTIN_PHASES`) with literal,
cited parameters and a full atomic basis, so every tab can reach them with no data files and no
optional dependency. On top of the materials, each panel registers its own **example scenarios** —
for the TEM solver, a simulated pattern down a named zone axis, already calibrated; for the crystal
viewer, a structure with a slip plane and Burgers direction superimposed; for the calculator, a
question with a known textbook answer. Examples are registered beside the operations they exercise
and appear in the manifest, so the same "no drift, discoverable" property of Decision 2 applies to
them.

The beta-zirconium entry is `Im-3m` with `a = 3.6090 Å` at 863 °C, rather than ferrite used as a
cubic symmetry stand-in. The value is tabulated by Maimaitiyili *et al.*, *J. Synchrotron Rad.* 22
(2015) 995–1000, [doi:10.1107/S1600577515009054](https://doi.org/10.1107/S1600577515009054),
from the Zuzek *et al.* zirconium phase assessment. Temperature belongs in the display name because
beta Zr is the high-temperature allotrope and its lattice parameter is not a room-temperature
material constant.

## Decision 8 — Every Plot Answers "What Is Under My Cursor?" And Lets It Be Inspected

Two interactions are mandatory on every plot in the application, not optional per panel, because
they are how a researcher interrogates a figure rather than merely looking at one:

1. **A live cursor readout.** The bottom corner of every plot shows the current pointer position in
   the plot's own physical units — Å⁻¹ for a diffraction pattern, degrees for a stereogram, Å for a
   crystal view — updating as the pointer moves. A plot whose axes carry units but whose cursor
   does not is a picture, not an instrument.
2. **Hover detail on every drawn entity.** Pointing at a diffraction spot shows its indices, its
   d-spacing, its |g|, its relative intensity, and — in a composite pattern — which phase and which
   variant produced it. The detail comes from the same row the CSV export writes, so what a user
   reads on screen and what they get in the file cannot disagree.

For two-dimensional figures, the same frame also owns one viewport language: the mouse wheel zooms
about the pointer, Shift-drag, middle-drag or the toolbar's **pan tool** moves the camera, the
toolbar reports the magnification, and **Fit** restores the complete figure. Zoom runs below 100% as
well as above it — a figure seen whole with room around it is as ordinary a request as a spot seen
closely — and the camera is bounded by its centre rather than its edges so that panning behaves the
same at every magnification. Zoom changes the SVG view box rather than applying a CSS transform,
so the coordinate readout and hover hit regions remain in the plotted coordinate system at every
magnification. A three-dimensional panel supplies its own camera and gets the same four controls
over it; the Crystal Viewer does, because dragging there rotates unless the pan tool or Shift says
otherwise.

Two further slots belong to the frame rather than to the stage below it, because both are read
against the drawing:

- **The control strip**, under the figure and inside the card, for anything that changes what is
  drawn — the packet legend of a variant pole figure, the source chips of a composite pattern. As a
  sibling below the card these were pushed under the fold by the result tables, so toggling a source
  and seeing the effect needed a scroll in each direction.
- **The overlay**, in the top-left of the drawing, for a panel's own live readout. The TEM solver
  puts the measurements taken off the picks there — d-spacing, ratio to the first pick, and the
  angle between them — because those are the three numbers a pattern is actually identified from,
  and they are worth having beside the pattern while picking rather than in a table under it.

The card is sized by the window, not by its content: it fills the visible stage bar a sliver, so the
complete figure and its controls are on screen when a panel opens and the result tables announce
themselves below. An SVG left in flow sizes itself from its own aspect ratio, which on a wide stage
made a square figure taller than the window; the drawing is therefore positioned out of flow inside
the card.

All of these behaviours are provided by one shared frontend module (`core/plotframe.js`) that every
panel mounts its SVG into, so a new plot gets them by construction and no panel re-implements them.
The readout formatter is chosen per plot from a declared unit; the hover payload is whatever row the
service attached to the entity.

This is why services attach the full row to each drawn entity rather than only its coordinates: the
extra bytes are the difference between a figure and a tool.

## Decision 9 — An Ambiguous Input Is Refused, Not Guessed

Users type `110` where they mean `1 1 0`, and for a long time the service layer read it — a
single digit per index, so the split is unambiguous. That reading is the defect, not the
convenience.

It only holds while every index is one digit. `10 10 0` typed the same way is `10100`, which is
`(10, 10, 0)`, `(1, 0, 100)` and `(101, 0, 0)` with nothing in the string to choose between them,
and a hexagonal four-index row makes it worse. The failure mode is the one this repository takes
most seriously: not an error, but a plausible number computed for indices nobody entered.

So the rule is stated once and applies to both halves of the application:

- **The workbench does not offer a control that can express the ambiguity.** Every `indices`
  parameter renders as one narrow box per index, and every `indices-list` parameter as a stack of
  those rows with add and remove buttons. There is no text field to run digits together in.
- **Each box is named**, because three identical fields in a row is the same ambiguity wearing a
  different hat. The letters come from the parameter's own label where the label states them —
  "Current zone axis [uvw]" names its boxes `u`, `v`, `w` — and otherwise from the width: `hkl`
  for three, the Bravais-Miller `hkil` for four. `pytex.app.registry._index_symbols` decides this
  in Python and publishes it as `symbols` in the manifest, so the browser never parses a label.
- **The service layer refuses the run outright**, including the cases it used to read correctly.
  Accepting the easy case is what teaches the habit that produces the ambiguous one. The message
  shows the separated form of what was typed.

The one check that stayed in the browser is the half-filled row. An empty box arrives at the
server as `''`, and "`''` is not a whole number" names a box the user cannot see a name for; by
the time the request is built, three boxes are one list with a hole in it.

## Decision 10 — What The Deployment Decides Lives In One YAML File, And Never Touches A Number

Three things are genuinely a site decision rather than a code decision: whether the feedback form
is offered, whether feedback is forwarded through an internal SMTP relay, and whether a first-time
visitor is greeted with the tour. `pytex.app.config` reads them from
`PYTEX_APP_CONFIG`, `./pytex_app.yml` or `~/.pytex/pytex_app.yml`, in that order, and works with
no file at all. `config/pytex_app.example.yml` is the annotated template, and a unit test loads it
through the same loader so a renamed key cannot stay documented under its old name.

The boundary is stated in the module's own docstring and is the test for whether anything may be
added to it: **nothing in this file changes a number PyTex computes.** A scientific result must
not depend on a file somebody edited.

Two consequences that are easy to get wrong:

- **An unknown key is an error, not a warning.** The usual way to get one is a misspelling, and a
  misspelled `smtp_host` is a relay that silently never delivers.
- **A configuration file that cannot be read does not take the workbench down.** It is logged and
  the defaults apply. The science does not depend on this file, so refusing to draw a pole figure
  because a relay host is mistyped would be a poor trade.

## Decision 11 — Feedback Is Stored Before It Is Sent

`POST /api/feedback` appends the submission to a JSON store — atomically, through a temporary file
and a rename — and only then attempts the relay. The ordering is the design. The alternative,
mailing it and storing it only if the mail failed, loses submissions exactly when the site is
having a bad day, which is when people have the most to say.

Three further properties follow from the same reasoning:

- **The store is not optional.** A deployment with no relay configured is not a deployment where
  feedback is lost; it is one where the person who runs the server reads a file.
- **The clock is the server's**, and so is which shell answered. A workstation with a wrong clock
  cannot file a note under next year, and a submission cannot claim to have come from a shell it
  did not. What the page reports — the open workspace and panel — is recorded as *claimed*
  context, separately from what the server knows.
- **`GET /api/experience` publishes whether a note will be mailed, never where or as whom.** The
  workbench has no authentication, so everything that route returns is readable by anyone who can
  reach the port; a relay host and an envelope sender are the administrator's business. A unit
  test asserts the exact key set for that reason.

The welcome tour is governed by the same route and the same file. It is skippable from every step,
the skip is remembered per browser and per version, and it is reachable again from the Help panel
— which is what makes remembering it safe rather than final.

## Frontend Architecture

No framework, but not ad hoc either. The frontend is four layers:

- `core/state.js` — a small observable store; panels subscribe to the slices they render.
- `core/api.js` — one `call(operation, params)` function over the JSON envelope, with error
  surfacing that shows the server's user-facing message and hint rather than a stack trace, and
  start/finish activity events consumed by the shared progress/history bar.
- `core/controls.js` — the manifest-driven control renderer (Decision 2), including the
  one-box-per-index Miller control of Decision 9.
- `core/feedback.js` — the feedback and feature-request drawer, built from the invitation the
  server publishes rather than from text in the page (Decision 11).
- `core/tour.js` — the welcome and the skippable tour (Decision 11).
- `panels/*.js` — one module per tab.

Layout follows the "visualisation gets the room" rule: a persistent tab bar, a single large canvas
region, a width-bounded control rail, and a compact activity bar that expands over the stage only
when its history is requested. At tablet width the control rail becomes a vertically bounded sheet
below the figure. Mid-size mastheads progressively omit duplicated tagline/action text before any
workspace navigation wraps, preserving vertical plot space.

## Workspace Inventory

Seven workspaces hold sixteen panels. A workspace is a *subject*; a panel inside it is a *view* of
that subject, which is why the two grouped workspaces have sub-tabs and the five single ones read
exactly as flat tabs did.

| Workspace | Panels | Services |
| --- | --- | --- |
| Crystal Viewer | the structure with superimposed planes, directions and annotations, and an orientation dock carrying a pole figure, an inverse pole figure and the crystal's Kikuchi map | `crystal.*` |
| TEM Analysis | SAED Simulator, TEM Solver, CBED, Composite SAED | `tem.*`, `cbed.*`, `diffraction.*` |
| XRD | powder peaks, broadened profiles, radiation and profile choices, indexed inspection | `xrd.*` |
| EBSD | IPF map, GROD, KAM, Scan summary, Distributions, Pole figures, Kikuchi simulator | `ebsd.*` |
| Variants | variant pole figures, packets, and the intervariant misorientation spectrum | `variants.*` |
| Texture | pole figures, inverse pole figures, ODF sections | `texture.*` |
| Calculator | interplanar angles, d-spacings, symmetry families, zone axes, cross-phase angles | `calc.*` |

Six of the seven EBSD panels read one scan that belongs to the session rather than to a panel; the
seventh, the Kikuchi simulator, reads none, because it is the forward problem the other six live
downstream of.

## Testing Obligations

- Every service operation has unit tests against known crystallographic answers, in the base lane.
- The manifest is validated structurally: every operation reachable, every parameter documented,
  every declared default accepted by its own validator. A missing help string fails a test.
- The HTTP layer is tested against a live loopback server started in-process, covering routing,
  error envelopes, content types, and expected client disconnects.
- Frontend logic that is worth testing (the projection maths) is mirrored by a Python
  implementation and checked against it, so the browser cannot silently disagree with the exporter.
- A dedicated Playwright/Chromium lane opens the real loopback application and checks all seven
  workspaces, the critical default calculation path in every panel, deliberate service-error
  surfacing, browser console/page errors, and the 390 × 844 responsive layout. Playwright is a
  test-only development dependency: no third-party JavaScript is shipped to the browser.
- The welcome tour is dismissed by the browser lane's own `openWorkbench` helper rather than
  disabled in the fixture, so every test in that lane re-proves the property that matters about a
  greeting: one click removes it and the application underneath is immediately usable.
- The feedback path is tested at both ends — the store, the relay and the route in the base lane
  against a fake SMTP server, and the form itself in the browser lane, including that a
  half-written note survives closing the drawer.

## See Also

- `docs/architecture/overview.md`
- `docs/standards/visualization_style_guide.md`
- `docs/architecture/tem_tilt_navigation_foundation.md`
- `docs/standards/data_contracts_and_manifests.md`
