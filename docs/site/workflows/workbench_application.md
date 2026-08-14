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

Opens a native window through `pywebview` if it is installed, and the default browser otherwise.
Use `--browser` to force the browser, and `--port` to pin the loopback port.

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

Seven workspaces. Every one of them ships with runnable examples, so a user with no data of their own
can still exercise every feature — the manifest test executes each example, so an example cannot rot.

| Workspace | What it answers |
| --- | --- |
| **Crystal Viewer** | What does this structure look like, with these planes and directions drawn on it? |
| **TEM Solver** | Which zone axis is this pattern, and how do I tilt to the next one? |
| **Diffraction** | What does a two-phase SAED pattern contain, and which variant is that spot? |
| **XRD** | Which powder peaks should this structure produce, and how do radiation and profile choices change the diffractogram? |
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

The colour-theme control cycles through **Auto**, **Light**, and **Dark**. Auto follows the operating
system; an explicit choice is remembered by the shared frontend, so it behaves identically in the
browser and desktop window. On a narrow screen the three masthead actions collapse to icons while
retaining their full titles and accessible names, and the seven workspace tabs wrap rather than
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

```{note}
The on-screen table is a **preview**: the first 200 rows, with the caption saying so. Every export
carries all of them. A texture ODF is over a thousand rows, and a scroll box with no search is not
how anyone reads a thousand numbers.
```

## Exporting

Two kinds of export, and the distinction is deliberate.

**A publication figure.** The Figure button on a plotting panel renders the same scene through
`pytex.plotting` in the journal style. The interactive view and the figure are the same geometry,
because the camera and the settings travel with the request.

Choose the format for the drawing, not by habit: a lit sphere is a mesh and a vector format writes
every facet of it, so PNG is right for the crystal viewer; a pole figure is line art, so SVG is
right there and is roughly twenty times smaller.

**Re-plottable numbers.** CSV, XLSX and JSON. The CSV is one row per entity at full precision; the
XLSX adds a sheet recording the inputs; the JSON is the complete result, reloadable. No result in
the application is exportable only as a picture.

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

```{warning}
The kernel halfwidth is a **smoothing choice, not a property of the material**. Too small and the
figure shows the individual grains that were sampled; too large and real detail is washed out. It is
the setting most often left unreported in the literature, which is why the application puts it in the
control rail rather than burying it in a default.
```

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
- {doc}`saed_pattern_solving` — the indexing the TEM Solver performs
- {doc}`crystal_visualization` — the renderer behind the Crystal Viewer's figures
- {doc}`stereographic_projections` — the projections the Variants and Texture panels use
- {doc}`../theory/pole_figure_arithmetic_and_mrd` — why the mean of a pole figure is 1 m.r.d.
