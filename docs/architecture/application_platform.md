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
- **Export** sends the camera back to Python, which renders the identical scene through the
  existing publication renderer at the requested size, DPI, and format.

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
| Zr-hcp | Hexagonal, so four-index notation and c/a effects appear rather than being described |

Each of these is in the built-in phase catalogue (`pytex.app.phases.BUILTIN_PHASES`) with literal,
cited parameters and a full atomic basis, so every tab can reach them with no data files and no
optional dependency. On top of the materials, each panel registers its own **example scenarios** —
for the TEM solver, a simulated pattern down a named zone axis, already calibrated; for the crystal
viewer, a structure with a slip plane and Burgers direction superimposed; for the calculator, a
question with a known textbook answer. Examples are registered beside the operations they exercise
and appear in the manifest, so the same "no drift, discoverable" property of Decision 2 applies to
them.

## Decision 8 — Every Plot Answers "What Is Under My Cursor?"

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

Both are provided by one shared frontend module (`core/plotframe.js`) that every panel mounts its
SVG into, so a new plot gets both behaviours by construction and no panel re-implements them. The
readout formatter is chosen per plot from a declared unit; the hover payload is whatever row the
service attached to the entity.

This is why services attach the full row to each drawn entity rather than only its coordinates: the
extra bytes are the difference between a figure and a tool.

## Frontend Architecture

No framework, but not ad hoc either. The frontend is four layers:

- `core/state.js` — a small observable store; panels subscribe to the slices they render.
- `core/api.js` — one `call(operation, params)` function over the JSON envelope, with error
  surfacing that shows the server's user-facing message and hint rather than a stack trace.
- `core/controls.js` — the manifest-driven control renderer (Decision 2).
- `panels/*.js` — one module per tab.

Layout follows the "visualisation gets the room" rule: a persistent tab bar, a single large canvas
region, and controls in a collapsible side rail plus a context strip under the canvas. Panels are
responsive down to tablet width, at which point the side rail becomes a sheet.

## Tab Inventory

| Tab | Status | Services |
| --- | --- | --- |
| Crystal Viewer | 3D structure, arbitrary superimposed planes/directions/annotations | `crystal.*` |
| TEM Pattern Solver | upload, calibrate, pick spots, index, plan the tilt to the next zone axis | `tem.*` |
| Calculator | interplanar angles, d-spacings, symmetry families, zone axes, cross-phase angles | `calc.*` |
| Diffraction | simulated SAED, Kikuchi maps, powder XRD | `diffraction.*` |
| Texture | pole figures, IPF, ODF | `texture.*` |

The first three are the starting scope; the last two follow, and both are already backed by
library code.

## Testing Obligations

- Every service operation has unit tests against known crystallographic answers, in the base lane.
- The manifest is validated structurally: every operation reachable, every parameter documented,
  every declared default accepted by its own validator. A missing help string fails a test.
- The HTTP layer is tested against a live loopback server started in-process, covering routing,
  error envelopes, and content types — with no browser.
- Frontend logic that is worth testing (the projection maths) is mirrored by a Python
  implementation and checked against it, so the browser cannot silently disagree with the exporter.

## See Also

- `docs/architecture/overview.md`
- `docs/standards/visualization_style_guide.md`
- `docs/architecture/tem_tilt_navigation_foundation.md`
- `docs/standards/data_contracts_and_manifests.md`
