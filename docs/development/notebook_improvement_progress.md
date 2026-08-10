# Notebook Improvement Progress

Durable handoff note for the tutorial-notebook overhaul under
`docs/site/tutorials/notebooks/`.

---

# Round 2 (opened 2026-08-10) — raise notebooks 01–25 to the tutorial-26–29 standard

## Objective

Round 1 (below, complete 2026-08-08) removed the machine-generated boilerplate and made
every notebook run, assert, and plot. It succeeded at that and stopped there. Tutorials
26–29 — orientation representations, TEM indexing, CBED, dynamical CBED — were then written
to a visibly higher standard, and the gap is now the defect. Round 2 closes it.

The instruction, in the requester's words: the earlier notebooks must

- meet the standard, rigour and depth of the recent ones;
- demonstrate not only PyTex's capabilities but **the teaching-worthy content of the subject
  matter itself** — what is interesting about rotations, about ODF inversion, about the
  reciprocal lattice;
- increase the focus on **the mathematics** behind the concepts and on **the algorithms** the
  implementation actually runs;
- carry auxiliary notes that produce an "awe" reaction in students and researchers alike.

## The rubric, distilled from tutorials 26–29

This is the contract. A notebook is done when every row holds. Tutorial 26
(`26_orientation_representations.ipynb`) and tutorial 28 are the reference exemplars; when
in doubt, open 26 and match it.

| # | Requirement | What it looks like in 26/28 |
| --- | --- | --- |
| R1 | **Open with the idea, not the API.** The first cell says what the subject *is* and why anyone should care, in the voice of someone who finds it interesting. No "this notebook introduces `Rotation`". | 26: "A crystal orientation is one geometric object. The numbers that denote it are many." 28: "One CBED exposure contains the rocking curve." |
| R2 | **Learning goals as questions**, five or so, that a reader could be examined on. | 26: "Why does uniform sampling of Euler angles produce a badly non-uniform set of orientations, and what fixes it?" |
| R3 | **Focused imports** in a `## 0. Setup` cell — only what the notebook uses. | |
| R4 | **Numbered sections**, each *theory → code → interpretation*: display-math derivation, then a cell that computes it live, then a markdown cell that says what the numbers mean. The interpretation cell is not optional; a printed table with no reading is half a section. | 28 §3 derives `s_g(θ)`, computes the disc fields, then reads the fringe orientation off the linearity |
| R5 | **At least one exact analytic value** verified live — an identity, a closed form, a group order, a conservation law — not a tolerance against a prior run. | 26: ⟨ω⟩ = π/2 + 2/π = 126.4756°; cubochoric cube volume π² |
| R6 | **At least one boxed algorithm note** (`> **Algorithm — …**` blockquote) giving the steps the implementation actually runs, including its guards. | 28 §5: the four-step two-beam thickness fit |
| R7 | **The awe note.** At least one auxiliary aside whose content is a genuine surprise: a measure that is not what intuition says, a degeneracy, an invariant, a number that comes out exact. Marked so the reader can find it. | 26 §6, the measure problem: large-angle rotations are *far* more numerous, so the obvious sampler is wrong |
| R8 | **A failure-mode section, deliberately triggered.** Run the wrong thing, show the wrong answer, show what the library does about it. | 28 §6 forces the wrong fringe order and gets an unphysical fit; 26 §13 mixes Rodrigues for homochoric |
| R9 | **Honest limits.** A short section stating what the implementation does *not* do, matching what `describe()` says. | 28 §8 |
| R10 | **Figures that teach**, each with a caption sentence saying what to look at. Multi-panel where a comparison is the point. | 26 §8.3: the cube face and its curvilinear image |
| R11 | **Close with takeaways** (imperative, memorable) **and Further reading**: the `docs/tex/` note, sibling tutorials by number and name, and literature with volume and page. | 28 §9 |
| R12 | **Notation per the registry.** `docs/standards/notation_and_conventions.md` — starred reciprocal *basis vectors*, unstarred Miller indices, `{hkl}`/`⟨uvw⟩` for families, overbars for negatives. | |
| R13 | Committed with **outputs and execution counts cleared**; runs error-free under the Sphinx build; `tests/unit/test_notebooks.py` green. | |

Two anti-requirements, also learned from 26–29:

- **Do not duplicate a sibling.** 26 owns the representation zoo; 02 must therefore own the
  *metric and the group quotient*, and point at 26 for the charts. Where two notebooks touch
  one subject, each states which half is its own and links the other.
- **No `plt.show()`.** Under pytest's warnings-as-errors it raises on the Agg backend.

## Scale reference

| | round 1 typical | 26–29 typical | round 2 target |
| --- | --- | --- | --- |
| cells | 14–16 | 29–54 | 30–50 |
| markdown characters | 4–6 k | 15–20 k | ≥ 12 k |
| figures | 1 | 4–8 | ≥ 3 |
| sections | 5 | 9–14 | ≥ 8 |

Character counts are a symptom, not the goal; a notebook that hits 12 k characters of padding
has failed the rubric. They are recorded because the round-1 notebooks fail them by a factor
of three, which is the honest measure of how much was missing.

## Target list and status

Order is by mathematical richness and by what the requester named explicitly (rotations, ODF
inversion) rather than by notebook number.

| Notebook | Subject it must own | Status |
| --- | --- | --- |
| 02 rotations & batch primitives | the metric on SO(3), the symmetry quotient, disorientation, quaternion averaging, the Mackenzie distribution | **done** — 34 cells, 25.6 k md chars, 3 figures |
| 06 texture/ODF/PF inversion | kernel normalization, the ghost problem, conditioning of the forward operator, regularization trade-off | pending |
| 03 symmetry & fundamental regions | point groups as groups, orbit–stabilizer, the fundamental sector's area, Laue vs proper | pending |
| 01 reference frames | the metric tensor as the frame, active vs passive, the crystal→Cartesian convention choice | pending |
| 04 lattices, space groups, CIF | metric tensor arithmetic, reciprocal duality, systematic absences from the space group | pending |
| 17 Miller vectorized | index arithmetic as linear algebra in dual bases, Miller–Bravais as a projection | pending |
| 08 diffraction geometry | Ewald construction, the excitation error, why zone-axis patterns look the way they do | pending |
| 11 powder XRD | structure factor as a lattice sum, multiplicity, the Lorentz–polarization factor | pending |
| 12 SAED | camera constant, the zone law, double diffraction vs true extinction | pending |
| 07 EBSD grid | KAM/GROD as discrete differential geometry on the grid, segmentation as a graph problem | pending |
| 16 EBSD → texture | from discrete orientations to a density: the estimator and its bandwidth | pending |
| 09 phase transformations | the OR as a rotation, variant counting from coset decomposition | pending |
| 10 plotting primitives | stereographic vs equal-area: what each preserves and the distortion cost | pending |
| 13 crystal visualization | projection geometry, and what a zone-axis view actually shows | pending |
| 05 acquisition manifests | provenance as a contract; the schema as the interface | pending |
| 14 YAML style | precedence resolution as a deterministic algorithm | pending |
| 15 structure → diffraction pipeline | the chain of contracts from one phase to three outputs | pending |
| 18–25 (OR + PF arithmetic track) | already hand-authored and strong; audit against R1–R13 and top up | pending |

## The companion capability: stereographic Kikuchi maps

Requested in the same instruction and tracked in `docs/development/active_task_progress.md`:
Kikuchi maps in the stereographic plane that guide a TEM operator in planning the path to the
next zone axis, plus tutorial 30 with cubic and hexagonal inline graphics. That work adds a
notebook rather than improving one, so it is ledgered there, not here.

## Working method

One notebook per commit, ledger row flipped in the same commit, pushed. Each notebook is
smoke-executed before commit (`tests/unit/test_notebooks.py` plus a direct cell execution of
the notebook under change), because the Sphinx build executes them and a notebook that raises
breaks the docs build for everyone.

---

# Round 1 (complete 2026-08-08) — remove the boilerplate, make every notebook run

## Objective

Rewrite the tutorial notebooks so each is a genuinely pedagogical, scientifically
rigorous, executable demonstration of the library: theory + conventions + citations,
focused imports (no kitchen-sink boilerplate), live-computed results verified against
analytic/literature values, and at least one figure. Notebooks must execute error-free
(the Sphinx build runs them) and pass `tests/unit/test_notebooks.py`.

## Motivation / baseline problem

Notebooks 01–16 were thin, machine-generated stubs carrying identical dead boilerplate
(an ~80-name kitchen-sink `from pytex import (...)` plus five unused helper functions
`make_context`, `describe_phase_fixture`, `load_zr_hcp_phase`, `load_diamond_phase`,
`publication_crystal_style`). Notebooks 17–21 (OR track) were already hand-authored and
richer — used as the quality template.

## Build tooling

- **Superseded 2026-08-08.** Notebooks are authored and edited by hand as `.ipynb` files,
  with no generator and no execution step. They are committed **without outputs**: the Sphinx
  site sets `nb_execution_mode = "cache"` and executes them at build time, and
  `nb_execution_raise_on_error = True` fails the build on a notebook that no longer runs.
  `scripts/execute_notebooks.py` existed only to bake outputs into the committed file and has
  been removed.
- Check: `tests/unit/test_notebooks.py` rejects any committed output, execution count, or
  run-specific metadata, and smoke-executes the priority notebooks' code cells.

## Status (rewritten + executed error-free, with verified assertions + figure)

- [x] 01 reference frames — frames/domains, FrameTransform apply/inverse/compose, VectorSet guard
- [x] 02 rotations/orientations — parameterizations, geodesic distance, misorientation/disorientation, OrientationSet batch
- [x] 03 symmetry — point/Laue/proper groups, operators, orbits/multiplicity, fundamental sector
- [x] 04 phases/CIF — metric tensor, cell volume, reciprocal d-spacings, atomic basis vs fixture site counts, CIF round-trip
- [x] 05 acquisition/manifests — AcquisitionGeometry, ExperimentManifest JSON round-trip + schema validation
- [x] 06 texture/ODF/PF inversion — kernel ODF, evaluate/volume-fraction, PF/IPF, invert_pole_figures recovers weights
- [x] 07 EBSD grid — CrystalMap, KAM edge detection, grain segmentation, GROD/GOS, IPF+KAM maps
- [x] 08 diffraction geometry — |g|=1/d, electron wavelength/Ewald, Bragg 2theta, kinematic spots + zone law
- [x] 09 phase-transformation foundations — K-S/N-W ORs, describe(), literature misorientation angles, variants, intervariant spectrum, PhaseTransformationRecord
- [x] 10 plotting semantic primitives — Wulff net, directions/planes/symmetry elements, PF/IPF
- [x] 11 powder XRD (priority) — Bragg + FCC selection rules + multiplicity/structure factor + anode shift
- [x] 12 SAED (priority) — zone law + camera constant R=C|g| + forbidden-spot extinction
- [x] 13 crystal visualization (priority) — CrystalScene, plane/direction overlays, Miller-Bravais, zone-axis view
- [x] 14 YAML style customization — theme→file→override precedence, two-theme render
- [x] 15 structure->diffraction pipeline (priority) — one phase → scene+XRD+SAED + workflow/validation manifests
- [x] 16 EBSD -> texture outputs — plane-direction seeding, texture_report ODF/PF/IPF, IPF+KAM maps
- [x] 17 miller vectorized — fully rewritten: d-spacings, cubic angles, families, Miller-Bravais, projection
- [x] 18-21 OR track — already excellent hand-authored; verified green; added house-style cross-links

## Constraints learned (API gotchas)

- `Misorientation.disorientation()` and `.angle_deg`/`.angle_rad` — disorientation is a METHOD.
- `misorientation_to(other, reduce_by_symmetry=False)` to show raw vs reduced.
- `Rotation.from_bunge_euler(phi1,Phi,phi2, degrees=True)` (degrees default True).
- `CrystalPlane.d_spacing_angstrom` and `.reciprocal_lattice_vector` are PROPERTIES.
- `Lattice.direct_basis().matrix` (returns Basis object).
- `CrystalMapPhase.name` must equal `phase.name` (e.g. "nickel-fcc"), not an arbitrary label.
- Segmentation `grain_sizes()`/`grain_orientation_spread_deg()` return dicts; `grod_map_deg()` is a method.
- pymatgen emits CIF warnings on the minimal fixtures; suppress with warnings.filterwarnings.
- Priority notebooks (04,11,12,13,15) must contain `get_phase_fixture`/`list_phase_fixtures`.
- notebooks index `docs/site/tutorials/notebooks.md` must mention each notebook stem.
