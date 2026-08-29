# Vision And Plan: Interface Crystallography, Composite Visualization, And Texture Quantification

**Status:** proposed program for the next development phase (drafted 2026-08-29).
**Precedence:** subordinate to `mission.md`, `specifications.md`, and the
[Critical Review And Development Guide](critical_review_and_development_guide.md); normative over
the visualization and texture-quantification portions of the older roadmaps once adopted. It
extends the [Orientation Relationship Analysis Foundation](../architecture/orientation_relationship_analysis_foundation.md)
(the OR flagship program) rather than replacing it, and continues that document's feature
numbering at **F15**.

**Purpose:** define the work that turns PyTex's existing orientation-relationship engine into an
integrated, publishable instrument for phase-transformation crystallography — one that renders the
parent and product crystals of any variant of any OR in three dimensions, draws the interface
between them, produces the whole crystallographic dossier of an OR in a single call, determines the
OR from two measured EBSD grains and shows it as a locked composite view — and, in the texture half
of the library, closes the three quantification gaps (Kearns factors in the GUI, ghost correction,
axial specimen symmetry) that currently bound what the tool may honestly claim.

The terminal deliverable is a peer-reviewed methods paper describing the algorithms and the
integrated tool. Section 9 states the paper plan and the gates that must close before submission.

---

## 1. The Thesis

A researcher studying a phase transformation asks a small number of questions over and over:

1. What is the orientation relationship, stated as parallel planes and directions in integer
   indices, and how far is my measurement from it?
2. Which variants does it produce, and what does each one *look like* sitting on the parent?
3. What is the interface between them — which plane, how well do the two lattices match on it, and
   what does the matching look like?
4. What are the numbers behind all of that: metric tensors, structure matrices, the correspondence
   matrix, the deformation gradient, the misorientation?
5. What will each variant give me in the microscope — SAED patterns, pole figures, stereograms?

Today these questions are answered by five different tools with five different conventions, and the
conventions are usually implicit. PyTex already answers most of them individually, with the
conventions explicit and pinned. **The thesis of this phase is that answering them together, from
one declaration of an OR, in one reproducible artifact, is itself the scientific contribution** —
because the errors that ruin transformation crystallography are boundary errors between tools, not
arithmetic errors inside them.

That is the claim the paper will make, and this program is what must exist for the claim to be true.

---

## 2. Verified Starting State

Checked directly against the working tree on 2026-08-29, not inherited from prior notes.

### 2.1 What is already built, and strong

| Capability | Where | Note |
| --- | --- | --- |
| `OrientationRelationship` with 12 named catalog constructors plus the generic parallel-plane/direction constructor | `core/transformation.py` (4278 lines) | Bain, KS, NW, GT, Pitsch, Burgers, Shoji-Nishiyama, Pitsch-Schrader, Potter, Bagaryatsky, Isaichev |
| Variant generation, packets, intervariant misorientation tables | `core/transformation.py` | Morito-validated; literature-correct counts |
| Index correspondence (`correspondence_direct` / `correspondence_reciprocal`), plane and direction mapping in both senses, rationalization with angular residual | `core/transformation.py` | F1–F3 |
| Parallelism finders over symmetry families | `find_parallel_planes`, `find_parallel_directions` | F4 |
| Variant correspondence tables with CSV/Markdown/JSON export | `variant_correspondence_table` | F2b |
| OR deviation, OR fitting, OR **determination** from measured pairs and from Euler angles | `or_deviation`, `fit_orientation_relationship`, `characterize_orientation_relationship`, `orientation_relationship_from_euler` | F5, F6, F6b — the double-coset seed is regression-pinned |
| OR determination with no parent, from child–child boundaries | `experimental/or_identification.py`, `experimental/or_refinement.py` | F7, experimental |
| Map-scale parent-grain reconstruction | `core/parent_reconstruction.py`, `experimental/parent_grain_reconstruction.py` | F8, experimental |
| Variant pole figures | `variant_pole_figure`, `plot_variant_pole_figure` | F10 |
| Deformation gradients with polar decomposition | `deformation_gradient()` | F12; Bain stretches $(1.127, 1.127, 0.797)$ pinned |
| Single-crystal 3D rendering with plane, direction, and cell overlays | `plotting/crystal3d.py` (2217 lines) | atoms, bonds, polyhedra, lit meshes, depth cue |
| **Composite two-crystal world scene from an OR** | `plotting/scene3d.py` — `WorldScene3D.from_orientation_relationship` | parallel-direction arrows and parallel-plane patches already drawn |
| Kearns factor by four independent routes | `texture/kearns.py` (1667 lines) | `kearns_from_orientations`, `kearns_from_pole_figure`, `kearns_from_odf`, `kearns_from_diffractogram`, plus `basal_tilt_profile` and `harris_texture_coefficients` |
| Harmonic ODF, PF→ODF inversion, discrete ODF, named-component fitting | `texture/harmonics.py`, `texture/models.py` | |
| SAED, CBED, Kikuchi, HOLZ simulation and composite SAED | `diffraction/`, `tem/` | the "SAED for each variant" half of the vision already exists |
| Self-describing GUI: 41 operations across 13 service modules, manifest-driven forms, browser-side live rotation | `app/` (6.7k Python lines, 15.4k JS lines) | Decision 4: Python decides the geometry, the browser only multiplies matrices |

### 2.2 The gaps this program closes

Each is a *verified absence*, not a suspicion:

1. **`WorldScene3D.from_orientation_relationship` has no `variant` parameter.** It places the child
   by `relationship.parent_to_child_rotation.inverse()` — variant 1 only. There is no way to render
   the 8, 12, or 24 variants a relationship produces, and no N-up contact sheet.
2. **The composite scene is not exposed in the GUI at all.** A search for `WorldScene` across
   `src/pytex/app` returns nothing. Its only consumers are one unit test and one tutorial notebook.
   The GUI's `variants` panel produces a stereographic pole figure and a misorientation table, and
   nothing three-dimensional.
3. **There is no interface-crystallography surface.** `TransformationVariant.habit_plane_pairs`
   exists and its own docstring says it is *descriptive only*; nothing computes it, nothing renders
   it, and there is no object representing a user-declared interface plane, its lattice misfit, or
   its in-plane matching.
4. **There is no aggregate OR report.** Metric tensors, structure matrices, correspondence
   matrices, deformation gradient, misorientation, variant table, parallelisms, and figures are all
   reachable — through eight separate calls returning eight separate types. Nothing assembles them.
5. **There is no OR stereogram.** `variant_pole_figure` plots child variant poles and optional
   parent poles on one net, but does not pair them, draw the great circles of the parallel planes,
   or annotate the parallelism deviations — the standard figure of every OR paper.
6. **Two measured grains cannot be taken to a picture.** `orientation_relationship_from_euler`
   answers the numerical question in Python. There is no GUI path, no reduction of the fitted
   statement to a *nearest-integer relationship object*, and no locked composite viewer.
7. **The Kearns parameter is absent from the GUI.** Searching `src/pytex/app` for `kearns` returns
   nothing. Four validated Python routes exist and none is reachable from the workbench.
8. **Ghost correction is not implemented.** `HarmonicODF.invert_pole_figures` says so in its own
   docstring, and `docs/site/theory/ghost_problem_and_odd_harmonics.md` states it plainly. This
   bounds every quantitative texture-strength claim the library makes.
9. **Axial (fibre) specimen symmetry does not exist.** `_SPECIMEN_SYMMETRY_POINT_GROUPS` in
   `core/symmetry.py` offers `triclinic`, `monoclinic`, `orthorhombic`, and `orthotropic` — nothing
   for wires and rods, which is the geometry the Kearns parameter is most often measured on.

---

## 3. Program A — The OR Dossier And Composite Crystal Visualization

Continues the OR foundation's feature numbering. Every feature lands with a theory note, a worked
example, a validation row, and `describe()` support, per the standing doctrine.

### F15 — Variant-aware composite scenes (**M**)

Make the composite two-crystal figure a function of the *variant*, not of the relationship alone.

```python
WorldScene3D.from_orientation_relationship(
    relationship,
    *,
    variant: int | TransformationVariant | None = None,   # NEW; None keeps today's behaviour
    ...
)

WorldScene3D.variant_scenes(relationship, *, variants=None, **kwargs) -> tuple[WorldScene3D, ...]
plotting.scene3d.render_variant_contact_sheet(scenes, *, columns=4, ...) -> Figure
```

The child placement rotation must be `variant.parent_to_child_rotation.inverse()`, consistent with
the regression-pinned composition
$g_{\text{child}} = g_{\text{parent}} \circ \mathbf{V}^{\mathsf{T}}$. The parallel-plane and
parallel-direction primitives must be re-derived **per variant** — the parent-side objects are the
symmetry images under that variant's operator, not the nominal pair, and drawing the nominal pair on
variant 17 is a silently wrong figure.

*Deliverable:* eight scenes for a Shoji-Nishiyama fcc↔hcp relationship, twelve for Burgers,
twenty-four for KS, each individually renderable and available as one contact sheet.

*Validation:* for every variant, the world-frame images of the defining parallel plane normals of
parent and child must coincide to $10^{-12}$; the set of child orientations across scenes must equal
`generate_variants()` as a set under child symmetry.

### F16 — Interface specification, quantification, and rendering (**L**)

A new module `pytex.core.interface`, treating the interface as a first-class crystallographic object
rather than a drawing.

```python
@dataclass(frozen=True)
class InterfaceSpec:
    relationship: OrientationRelationship
    variant_index: int
    parent_plane: CrystalPlane                    # the terrace plane, parent side
    child_plane: CrystalPlane | None = None       # None -> mapped through the correspondence
    parent_in_plane_direction: CrystalDirection | None = None   # reference direction in the interface
    child_in_plane_direction: CrystalDirection | None = None

@dataclass(frozen=True)
class InterfaceReport:
    """What the interface *is*, numerically."""
    plane_deviation_deg: float           # angle between the two plane normals in the world frame
    d_spacing_parent: float
    d_spacing_child: float
    d_spacing_mismatch: float
    in_plane_misfit: np.ndarray          # 2x2 planar misfit strain in the interface basis
    principal_misfit: tuple[float, float]
    invariant_line_in_plane: CrystalDirection | None   # a direction of zero misfit, if one exists
    coincidence: InterfaceCoincidence | None           # near-coincident interfacial site lattice
    step_character: StepCharacter | None               # terrace/ledge/riser, when a step is declared
    def describe(self) -> str: ...
```

Rendering, in `plotting/scene3d.py`:

- the interface as a shared translucent patch in the world frame, with both crystals optionally
  **cut** by its half-space, so the joint is visible rather than buried inside two interpenetrating
  boxes;
- a companion 2-D **interfacial matching diagram**: both lattices projected into the interface
  plane and superimposed, with the near-coincidence cell outlined and the misfit annotated — the
  figure that makes "good matching" a measurement instead of an adjective;
- terrace/ledge decomposition drawn when a step direction is declared.

*Honest scope statement, to be repeated in the docstrings:* F16 analyses an interface **the user
declares**. It does not *predict* the habit plane. Prediction is F20.

*Validation:* a coherent $\{111\}_{\gamma} \parallel \{110\}_{\alpha}$ terrace under KS must report
zero plane deviation and the correct planar misfit from the lattice parameters alone; the
invariant-line search must reproduce the classic $\langle 110 \rangle$-family results for the Bain
strain; the misfit tensor must be invariant to the choice of in-plane basis.

### F17 — The OR dossier: one declaration, one artifact (**M**)

```python
pytex.core.transformation.or_dossier(
    relationship, *, variant=None, poles=(...), directions=(...), interface=None,
) -> ORDossier
```

`ORDossier` is a typed, explainable, serializable aggregate with five blocks:

1. **Lattice block** — both phases' cell parameters, direct and reciprocal metric tensors, structure
   matrices $\mathbf{A}_p$ and $\mathbf{A}_c$, cell volumes, point groups.
2. **Transformation block** — the rotation $\mathbf{R}$, the direction correspondence $\mathbf{C}$,
   the plane correspondence $\mathbf{C}^{-\mathsf{T}}$, the deformation gradient $\mathbf{F}$ with
   its polar decomposition, principal strains, and volume change.
3. **Misorientation block** — the symmetry-reduced axis/angle representative, the variant count and
   grouping, the intervariant misorientation spectrum.
4. **Parallelism block** — defining parallelisms plus discovered near-parallelisms from F4, with
   deviations, in publication notation via `pytex.core.notation`.
5. **Figure block** — the composite 3-D scene per variant (F15), the variant pole figure (F10), the
   OR stereogram (F18), the interface figures (F16), and the per-variant SAED patterns from the
   existing diffraction stack.

`ORDossier.describe()` produces the prose; `to_json()` closes part of F11; `export(directory)` writes
the whole bundle — figures as SVG, tables as CSV and Markdown, numbers as JSON against a schema under
`schemas/`. **That bundle is the paper's reproducibility artifact.**

*Rule:* the dossier must call the existing functions, never reimplement them. A dossier number that
disagrees with the function it came from is the exact class of defect this repository exists to
prevent.

### F18 — The OR stereogram (**S–M**)

One stereographic net carrying:

- parent poles (open symbols) and child poles (filled) for the nominated families, in the parent
  frame, for one variant or all;
- **tie-lines** joining OR-parallel pairs, labelled with the deviation in degrees;
- great circles of the parallel planes, so a plane parallelism reads as two coincident circles rather
  than two coincident points;
- zone-axis and trace annotation consistent with the notation standard.

This is the figure by which ORs are read in the literature, and it is the natural teaching object for
the parallelism statement.

### F19 — The composite crystal viewer in the workbench (**L**)

The GUI deliverable, and the one a user will actually touch. A new workspace (or a promoted
`variants` workspace) with these operations:

| Operation | Returns |
| --- | --- |
| `variants.composite_scene` | both crystals' `scene_payload`s in a common world frame, the OR primitives, and the variant's rotation matrix |
| `variants.contact_sheet` | the same for every variant, for the N-up grid |
| `variants.interface` | the `InterfaceReport` plus the interface patch and the 2-D matching diagram |
| `variants.dossier` | the `ORDossier` JSON plus its figure bundle |

Browser behaviour, following the existing Decision 4 split — **all crystallography stays in Python;
the browser multiplies matrices only**:

- one camera rotation $\mathbf{R}_{\text{cam}}$ drives *both* crystals; the child's scene is
  pre-placed in the world frame by Python, so the lock is free and cannot drift;
- a variant selector, and a contact-sheet mode showing every variant at the current camera;
- toggles for parallel planes, parallel directions, the interface patch, the cut-away, cells, bonds,
  and labels — reusing `crystal.js`'s existing appearance controls and depth sort;
- side-by-side and interpenetrating placement modes;
- a live readout of what is currently parallel to the screen normal.

*Reuse, not rewrite:* `scene_payload` (`app/services/crystal.py`), `rotation3.js`, `plotframe.js`,
and the painter's-algorithm depth sort in `crystal.js` all carry over. The new JavaScript is a
composition layer, not a new renderer.

### F20 — PTMC and habit-plane prediction (**L**, stretch, explicit go/no-go)

Invariant-plane-strain analysis: lattice-invariant shear options, predicted habit-plane normals,
predicted OR, and the shape-strain magnitude — populating the long-dead `habit_plane_pairs` slot with
computed, provenance-carrying results. This is OR-foundation **F13**, open since the July review.

*Recommendation:* schedule it, but **gate it behind F15–F19 and Program C**. It is the single item
most likely to elevate the paper from "an integrated tool" to "an integrated tool that predicts", and
also the single item most likely to consume the whole phase. Decide at milestone M4 (§10) with
F15–F19 already landed, so that a no-go costs nothing.

---

## 4. Program B — From Two Measured Grains To A Picture

### F21 — Measured-pair OR determination, reduced to integers (**M**)

The numerical core exists: `orientation_relationship_from_euler` and
`characterize_orientation_relationship` return an `ORCharacterizationReport` with the fitted rotation,
the catalog ranking, the recovered parallelism statement, and an `is_conclusive` verdict. Three
things are missing.

1. **A rational reduction that produces an object, not prose.**

   ```python
   ORCharacterizationReport.as_rational_relationship(
       *, max_index: int = 4, tolerance_deg: float = 3.0,
   ) -> RationalizedORResult
   ```

   returning a genuine `OrientationRelationship` built from the nearest-integer parallel plane and
   direction pair, **plus** the angular cost of the idealization: the deviation of each rationalized
   index from the exact one, the residual rotation between the measured and idealized relationships,
   and provenance stating in plain language that this is an idealization and how large it is. A
   rational OR handed back without its cost is exactly the kind of silent boundary error the library
   is built to prevent.

2. **Multi-pair aggregation.** A list of grain pairs feeds `fit_orientation_relationship` (exists),
   with residual statistics, outlier flags, and a per-pair deviation table.

3. **The GUI operation** `variants.or_from_grains`, taking two Euler triples with an explicit
   convention selector (Bunge by default; the convention machinery already exists in
   `core/conventions.py` and `app/services/crystal.py`), two phases, and an optional tolerance —
   returning the report, the catalog ranking, the rationalized statement, and the deviation.

### F22 — The locked composite viewer for measured grains (**M**)

The F19 viewer, seeded from the measured pair instead of a catalog OR:

- both crystals placed by their *measured* orientations, so the view is the specimen frame the EBSD
  data came in;
- rotating one rotates the other, because the relative placement is fixed by the measurement;
- the **measured** parallel planes and directions drawn as overlays, found by `find_parallel_planes`
  and `find_parallel_directions` at the fitted rotation and labelled with their actual deviations — a
  visual guide that is honest about being approximate;
- a toggle between the measured relationship and its rationalized idealization (F21), so the user
  sees what the idealization costs geometrically;
- the OR stereogram (F18) beside the 3-D view, sharing the camera.

### F23 — Closing the loop from the map (**S–M**)

The `ebsd.map` panel already loads and segments scans. Add: pick two grains on the map, and their
mean orientations flow into F21 with no retyping. This makes the whole path — scan, grains, OR,
integer statement, composite figure — a workbench workflow rather than a scripting exercise.

---

## 5. Program C — Texture Quantification

### T1 — The Kearns tab (**M**)

A new subtab in the texture workspace exposing all four existing routes, so a user can compute the
same specimen's Kearns factors three independent ways and compare them — precisely the comparison the
literature never makes conveniently, and a natural paper figure.

| Operation | Python surface (exists) | Input |
| --- | --- | --- |
| `texture.kearns_from_diffractogram` | `kearns_from_diffractogram`, `harris_texture_coefficients`, `basal_tilt_profile` | reflection list with $2\theta$ and integrated intensity, plus a random-powder reference |
| `texture.kearns_from_pole_figure` | `kearns_from_pole_figure` | a measured or imported pole figure |
| `texture.kearns_from_odf` | `kearns_from_odf` | a reconstructed or modelled ODF |
| `texture.kearns_from_orientations` | `kearns_from_orientations` | EBSD orientations |
| `texture.kearns_tilt_profile` | `kearns_from_tilt_profile`, `basal_tilt_angle_deg` | a $\chi$-tilt intensity profile |

The panel must show the **triad** $(f_1, f_2, f_3)$ with its closure check $\sum_i f_i = 1$, the pole
orientation tensor, the specimen-frame axis labels, the assumed specimen symmetry, and the
`describe()` prose. Where the input is a diffractogram, the Harris coefficients and the basal tilt
profile are shown as intermediate evidence, because a Kearns factor quoted without them is a number
no reviewer can check.

*Cost:* this is almost entirely plumbing — registry parameter declarations, service handlers, one JS
panel, example scenarios, and Playwright coverage. The science is already validated. It is the
highest impact-per-hour item in the whole program.

### T2 — Ghost correction (**L**)

The largest scientific-credibility gap in the texture half of the library, named as such by the
2026-08 capability review. Pole figures determine only the even-order coefficients under Friedel's
law; the odd part is currently left unconstrained, and the reconstruction must be read as the even
part alone.

```python
HarmonicODF.invert_pole_figures(
    ..., ghost_correction: GhostCorrectionSpec | None = None,
) -> HarmonicODFReconstructionReport   # gains a GhostCorrectionReport
```

Two methods, both classical, both declared explicitly in the output:

1. **Zero-range method (Bunge–Esling).** Identify the zero range from the even-part reconstruction,
   then determine the odd coefficients by requiring $f \geq 0$ everywhere and $f = 0$ in the zero
   range.
2. **Positivity / iterative method.** A non-negative or WIMV-type iterative solution as the
   alternative, for cases where the zero range is small or ill-defined.

`GhostCorrectionReport` must carry the method, the zero-range volume fraction, the recovered odd
part's contribution, the change in texture index and entropy, and an honest statement of what remains
unconstrained. `docs/site/theory/ghost_problem_and_odd_harmonics.md` — which currently states the
absence — is updated in the same change.

*Validation:* synthetic ODFs with a known odd part must be recovered to a stated tolerance; the
texture index must move toward truth, not merely change; an MTEX parity row is required before any
document claims parity.

### T3 — Axial (fibre) specimen symmetry (**M**)

Wires, rods, and drawn products are the geometry the Kearns parameter is most often reported on, and
PyTex cannot express their sample symmetry. Add axial symmetry as a first-class specimen symmetry:

```python
SymmetrySpec.specimen("axial", axis=..., reference_frame=...)   # infinity-fold about the specimen axis
```

Two implementation paths, both needed, both explicit:

- **Harmonic path (exact).** In the Bunge expansion, fibre symmetry about the specimen axis
  annihilates every specimen index except $n = 0$. This is an exact projection, not an
  approximation, and it is where the axial assumption should live for ODF work.
- **Discrete path (approximate, honest).** For discrete ODFs, pole figures, and orientation sets,
  approximate the continuous group by an $N$-fold rotation set with $N$ declared in the output and a
  stated error bound. The API must never silently choose $N$.

Downstream: pole-figure and IPF plotting must honour it, and the Kearns triad must reduce to
$(f_{\perp}, f_{\perp}, f_{\parallel})$. The α-Zr case already recorded in the ledger,
$(0.488, 0.488, 0.053)$, is the pinned regression.

### T4 — Assumption transparency in the measured-ODF service (**S**)

With T2 and T3 landed, the `texture.measured_pole_figures` and `texture.odf_sections` operations must
state, in every result: the defocus correction applied, whether ghost correction ran and by which
method, and the imposed specimen symmetry. An assumption the GUI does not display is an assumption
the user did not make.

### T5 — Uncertainty on the reported numbers (**M**)

Bootstrap confidence intervals on Kearns factors, volume fractions, and texture index. The paper will
quote these numbers; quoting them without intervals invites the reviewer question that has no good
answer.

---

## 6. Program D — Cross-Cutting Engineering

- **JSON contracts (F11 closure).** Schemas under `schemas/` for `OrientationRelationship`,
  `InterfaceSpec` / `InterfaceReport`, `ORDossier`, `ORCharacterizationReport`, and `KearnsReport` —
  reconstruction-grade, not lossy summaries.
- **A dossier CLI.** `python -m pytex.dossier --relationship kurdjumov_sachs --parent austenite_fcc
  --child fe_bcc --variants all --out figures/`, producing the complete figure and table bundle from
  one command. This is what the paper's data-availability statement points at.
- **Vectorization.** `intervariant_misorientations` and the per-element `Rotation` list building in
  `PhaseTransformationRecord` remain Python loops — already flagged P2, and now on the hot path of
  the contact sheet and the dossier.
- **Browser coverage.** Playwright tests for every new panel, in the existing critical lane.
- **Parity ledger rows.** MTEX rows for ghost correction and OR fitting; VESTA rows for the composite
  and interface renderings.

---

## 7. High-Impact Opportunities Found In The Repository Scan

Outside this program's scope, ranked by impact, recorded so the next planning pass does not have to
rediscover them:

1. **The MTEX parity campaign is defined but never executed.** `fixtures/mtex_parity/campaigns/` has
   the MATLAB handler; no MATLAB/MTEX installation exists in the development environment. This is the
   single largest credibility risk for a paper that positions itself against MTEX. Resolving it needs
   a machine, not code — treat it as a logistics task and start it now, because it has the longest
   lead time of anything here.
2. **Out-of-core / chunked EBSD backing.** The EDAX HDF5 reader loads a scan whole, which bounds the
   size of real dataset PyTex can ingest — the most common practical reason a user would leave.
3. **Measured powder-XRD I/O and multiphase fitting.** Scattering, profiles, and preferred
   orientation exist; measured-pattern import, background and multiphase fitting, and refinement do
   not. The Kearns $2\theta$ route in T1 makes measured-pattern import newly urgent.
4. **The stable→experimental import boundary defect** — stable `core.parent_reconstruction` imports an
   experimental scoring primitive. A stability-policy violation in the flagship subsystem.
5. **The 558-symbol flat `__all__`.** Growing without a namespace plan; every feature in this program
   adds to it.
6. **No performance regression lane in CI.** Benchmarks exist and run; nothing guards against
   regression.
7. **Ring and polycrystalline SAED, and specimen-thickness distributions** — the remaining breadth gap
   in the diffraction stack, and a natural companion to the per-variant SAED figures.

---

## 8. Validation Program

Per repository doctrine, features are incomplete until all four rows exist.

**Analytic.** Composite scene round-trips (world-frame parallel normals coincide to $10^{-12}$ for
every variant); correspondence composed with its inverse is the identity; misfit tensors invariant to
in-plane basis choice; the fibre-symmetry harmonic projection idempotent.

**Literature.** KS $42.85^{\circ}$ about $\langle 0.968\;0.178\;0.178 \rangle$; Bain principal
stretches $(1.127, 1.127, 0.797)$; Morito intervariant and packet tables; Burgers 12-variant tables;
the α-Zr Kearns triad $(0.488, 0.488, 0.053)$; classical habit-plane results if F20 proceeds.

**Parity.** MTEX for OR fitting, parent reconstruction, and ghost-corrected ODFs; VESTA for the
composite and interface renderings. **No document may claim parity before the comparison has actually
run** — the standing rule, and the one most at risk while the MATLAB machine is missing.

**Property-based.** Random lattices and ORs: a plane and its symmetric equivalents must map to the
same child family; rationalization residuals invariant to index scaling; the interface report
invariant to relabelling the two phases where the geometry is symmetric.

**Measured fixtures.** The outstanding item across F6b, F8, and now F21–F23: real EBSD data with a
known relationship. Acquiring or licensing one dataset is a prerequisite for the paper's measured case
study.

---

## 9. The Paper

**Working title.** *Semantically explicit orientation-relationship and interface crystallography: an
integrated open-source tool for phase-transformation analysis.*

**Target.** *Journal of Applied Crystallography* (IUCr), Computer Programs section — the natural home
for a crystallographic tool whose contribution is convention explicitness and integration, and whose
audience already accepts IUCr notation as normative. Alternatives, in order: *Computer Physics
Communications*, *Materials Characterization*, *Journal of Microscopy*, *SoftwareX* (lowest weight;
short format).

**The claim.** Not "we implemented known crystallography" — everyone has. The claim is:

1. an explicit separation of the **three objects** (rotation, index correspondence, deformation) that
   transformation-crystallography tools routinely conflate, carried through a whole toolchain;
2. an OR **determination** algorithm from measured pairs that is honest about identifiability — the
   double-coset seed and the `is_conclusive` verdict, with the measured discrimination limits
   (conclusive to $2^{\circ}$ of scatter, correctly inconclusive at $5^{\circ}$, against the
   $2.40^{\circ}$ that separates KS from GT);
3. a **rationalization with a stated cost** — integer OR statements that carry the angular price of
   being integers;
4. **integration as a scientific result**: one OR declaration producing variants, interfaces,
   stereograms, pole figures, deformation gradients, and SAED patterns in one reproducible bundle,
   with every convention stated and machine-checked.

**Structure.**

1. Introduction — the boundary-error problem in transformation crystallography.
2. Conventions and theory — the three-objects doctrine; variant identity as convention; interface
   geometry; the notation standard.
3. Algorithms — OR determination and its identifiability limits; rationalization and its cost;
   variant-aware composite placement; interface misfit and coincidence; ghost correction; the Kearns
   triad by four routes.
4. Implementation — the self-describing operation registry, the Python/browser split, the
   explainable-results doctrine, the parity ledgers.
5. Validation — literature pins, MTEX and VESTA parity, property-based invariants.
6. Case studies — (a) KS fcc→bcc: 24 variants, dossier, interface, per-variant SAED; (b) Burgers
   β→α Zr: the same specimen's Kearns factor by diffractogram, pole figure, and ODF, under axial
   specimen symmetry; (c) a measured two-grain OR from EBSD, determined, rationalized, and rendered.
7. Availability and reproducibility — the dossier CLI, the version DOI, the fixture set.

**Figures (target 8–10), all produced by the tool:** the three-objects schematic; the 24-variant
contact sheet; one composite variant with its parallel planes and directions; the OR stereogram with
tie-lines; the interface patch plus the 2-D matching diagram; the variant pole figure; the
identifiability curve (deviation versus scatter, with the KS/GT separation marked); the three-route
Kearns comparison; ghost correction before and after; the measured two-grain case.

**Gates that must close before submission.**

| Gate | Status today | Action |
| --- | --- | --- |
| MTEX parity numbers actually generated | **Blocked** — no MATLAB/MTEX machine | Logistics; start immediately, longest lead time |
| A measured EBSD dataset with a known OR | **Missing** | Acquire or license; needed for case study (c) |
| Ghost correction validated | **Not implemented** | T2 |
| Kearns reachable in the tool being described | **Not in the GUI** | T1 |
| Axial symmetry for the wire case | **Not implemented** | T3 |
| Version DOI and `CITATION.cff` pinned to the release described | `CITATION.cff` exists; needs a Zenodo DOI and a tagged release | Program D |

The paper must describe the tool **as it will be at submission**, and every claim in it must be
traceable to a passing test. That is not a stylistic preference here; it is the repository's own
standing rule applied to the manuscript.

---

## 10. Sequencing

Each milestone ends on `main`, green, with its ledger entry — per the cardinal resumability rule.

| Milestone | Content | Size | Why here |
| --- | --- | --- | --- |
| **M1** | T1 Kearns tab | M | Highest impact per hour; the science is already validated, only the surface is missing. Ships user-visible value in one focused session. |
| **M2** | F15 variant-aware scenes, F18 OR stereogram | M | Pure Python, no GUI risk; unblocks everything visual. |
| **M3** | F19 composite viewer in the workbench, F17 dossier | L | The centrepiece. The dossier is assembled from F15/F18 and existing calls, so it comes cheap once M2 lands. |
| **M4** | F21–F23 measured-pair workbench; **F20 go/no-go decision** | L | Closes the EBSD-to-figure loop. Decide PTMC here, with M1–M3 already banked. |
| **M5** | T3 axial symmetry, then T2 ghost correction | L | T3 first: smaller, and the Kearns wire case needs it. T2 is the long pole of Program C. |
| **M6** | F16 interface crystallography; Program D contracts, CLI, parity rows; T5 uncertainty | L | Interface work benefits from every earlier surface; the CLI is the paper's artifact and must come last so it captures the finished tool. |

Parallel, starting **now** and not on this critical path: the MATLAB/MTEX machine, and the measured
EBSD dataset. Both have lead times measured in weeks and gate the paper, not the code.

---

## 11. Risks And Honest Limits

- **Scope.** This is a large program. F20 (PTMC) is the item most likely to swallow it; the M4 gate
  exists for that reason. Program C alone is a phase's worth of work if ghost correction proves
  stubborn.
- **The parity claim.** Without a MATLAB machine, no MTEX parity claim may appear in the code, the
  docs, or the manuscript. This is the standing rule, and it will feel expensive at submission time.
- **Rationalization is an idealization.** F21's integer statements are useful and lossy. Every surface
  that produces one must carry its cost, or the feature actively misleads.
- **Axial symmetry is an assumption about the process, not the specimen.** Imposing it on a specimen
  that is not axisymmetric fabricates symmetry. The API must make it opt-in and must say so.
- **Browser complexity.** The composite viewer roughly doubles the geometric state the frontend holds.
  Decision 4 must not erode: if crystallography starts migrating into JavaScript to make the viewer
  feel faster, the coherence claim the paper rests on quietly dies.
- **Figure count versus depth.** An "everything in one go" dossier can become an unreadable dump. Each
  block must stand alone as a publishable figure or table, or it does not belong in the dossier.

---

## References

### Normative

- [Critical Review And Development Guide](critical_review_and_development_guide.md)
- [Orientation Relationship Analysis Foundation](../architecture/orientation_relationship_analysis_foundation.md)
- [Phase Transformation Foundation](../architecture/phase_transformation_foundation.md)
- [Application Platform](../architecture/application_platform.md) — Decision 4, the Python/browser split
- [Notation And Conventions](../standards/notation_and_conventions.md)
- [Data Contracts And Manifests](../standards/data_contracts_and_manifests.md)

### Informative

- Bunge, H. J., *Texture Analysis in Materials Science* — harmonic method, ghost problem, specimen symmetry
- Esling, C., Bunge, H. J., et al. — the zero-range method for odd-order coefficients
- Matthies, S., Vinel, G. W., Helming, K. — WIMV and positivity-based ODF reconstruction
- Kearns, J. J. — the basal-pole texture parameter for zirconium alloys
- Morito, S., et al. — lath-martensite variants, packets, and the V1–V24 convention
- Burgers, W. G. — the bcc-to-hcp relationship
- Bhadeshia, H. K. D. H., *Geometry of Crystals, Polycrystals, and Phase Transformations* — PTMC, invariant-line analysis
- Bollmann, W., *Crystal Defects and Crystalline Interfaces* — the O-lattice and interfacial coincidence
- Sutton, A. P., Balluffi, R. W., *Interfaces in Crystalline Materials* — interfacial misfit and step character
- MTEX documentation — `calcParent2Child`, parent-grain reconstruction, ODF estimation
