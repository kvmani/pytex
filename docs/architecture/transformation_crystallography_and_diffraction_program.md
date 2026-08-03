# Transformation Crystallography And Composite Diffraction Program

**Status:** specification (2026-08-03). Program identifier: **TX**.
**Owner document:** this file is the normative specification. The running phase ledger is the
repository-local working note `docs/roadmap/working_notes_transformation_diffraction_program.md`,
which is a development record and deliberately not part of the rendered documentation site. The
repository documentation index links it.

---

## 1. Purpose

PyTex already has strong, separately-verified pieces of transformation crystallography: an
`OrientationRelationship` model with correspondence constructors, variant enumeration, fitting,
parallelism finders, intervariant fingerprints, parent-grain reconstruction, and a composite
kinematic SAED engine. What it does **not** have is a coherent *user-facing answer path* for the
five questions a phase-transformation researcher or student actually asks:

| Ask | Question | Program feature |
| --- | --- | --- |
| (a) | I measured Euler angles for a parent grain and a child grain by EBSD. **What is the orientation relationship?** | **TX1** |
| (b) | Given an OR, for an arbitrary parent direction or plane, **what are the parallel directions and planes in every product variant?** | **TX2** |
| (c) | Given a parent zone axis, give me a **composite kinematic SAED** containing all (or selected) variants — robust, and exportable as both graphics and reflection tables. | **TX3** |
| (d) | Let me instead **choose a zone axis of a product variant**, and generate the matrix and the other variants from that. | **TX4** |
| (e) | Here is a **measured SAED pattern**. Let me pick spots (by hand, or from a YAML file) and **solve it**. | **TX5** |

The program's job is to turn existing primitives into these five answers, with the repository's
standing quality contract: explicit conventions, `describe()` explainability, executable worked
examples with independent provenance, publication-grade figures, and tests that bite.

### Non-goals (explicitly out of scope for TX v1)

- Dynamical (Bloch-wave / multi-beam) diffraction intensities. Kinematic only, stated everywhere.
- HOLZ / higher-order Laue zone rings and Kikuchi-line overlays (tracked separately).
- A desktop GUI. TX5's "interactive" surface is a Matplotlib picker plus a YAML file contract,
  which is scriptable, testable and headless-safe.
- Phenomenological theory of martensite crystallography (PTMC) habit-plane prediction — that is
  OR-foundation F13 and has its own program.
- Automated spot *detection* from raw image data beyond what `DiffractionPattern.cluster_observations`
  already offers. TX5 consumes picked/listed spot coordinates.

---

## 2. Conventions pinned for the whole program

These are inherited, not invented here. Restated because every TX surface depends on them.

1. **Orientation convention.** Bunge ZXZ, passive, crystal→specimen matrices stored such that a
   crystal-frame vector `v_c` maps to specimen as `v_s = M v_c`, where `M = Orientation.as_matrix()`.
2. **Parent→child rotation.** `OrientationRelationship.parent_to_child_rotation` is the rotation
   `R` for which a parent crystal vector expressed in parent Cartesian coordinates becomes
   `R v_p` in child Cartesian coordinates. For measured pairs the operative rotation of a pair is
   `V = C^T P` (see `fit_orientation_relationship`), with `P`, `C` the parent/child crystal→specimen
   matrices. **All TX code uses this single definition; no local re-derivation.**
3. **Variant indexing** is 1-based and follows `OrientationRelationship.generate_variants()` order
   (parent-symmetry operators applied on the parent side, reduced by the child-symmetry orbit).
4. **Planes vs directions.** Plane normals transform with the reciprocal basis, directions with the
   direct basis. Every TX surface that maps indices routes through
   `OrientationRelationship.map_plane_to_child` / `map_direction_to_child`, never through a bare
   matrix product.
5. **Notation** is produced only by `pytex.core.notation`: `(hkl)` / `[uvw]` for specific,
   `{hkl}` / `<uvw>` for families, overbars for negatives, stars on reciprocal *axes* only.
6. **Diffraction geometry.** The composite engine's shared detector basis is parent-anchored
   (`zone_basis_from_axis`), with column 2 the zone axis; excitation error drives the intensity
   weighting; camera constant `L·λ` relates `|g|` to detector radius.
7. **Units.** Angles in degrees on all public surfaces (`_deg` suffix mandatory), lengths in mm on
   the detector, reciprocal lengths in Å⁻¹, wavelengths in Å.

---

## 3. TX1 — Orientation relationship from measured orientations

### 3.1 The scientific question

EBSD gives Euler angle triples for a parent grain and one or more child grains. The user wants:

1. the operative parent→child **rotation** (axis/angle representative in the disorientation sense),
2. **which named OR it is** (or that it is none of them), with an honest margin statement,
3. the **crystallographic statement** — which parent planes are parallel to which child planes, and
   which parent directions to which child directions — because that, not a matrix, is how the
   literature reports an OR,
4. a quantified **scatter/quality** figure so the answer can be trusted or rejected.

### 3.2 What exists and what is missing

Existing: `fit_orientation_relationship` (symmetry-aware eigen-mean refinement from a nominal),
`or_deviation`, `OrientationRelationshipCatalog` + `standard_*_relationships` builders,
`find_parallel_planes` / `find_parallel_directions`, `identify_orientation_relationship`
(child–child boundaries only, requires candidates).

Missing:
- a **catalog-ranking** entry that does not require the user to already guess the nominal;
- the **parallelism extraction** from a fitted rotation into a human-readable OR statement;
- an ergonomic **Euler-angle entry point** (the user has degrees in a spreadsheet, not
  `OrientationSet` objects);
- one **report object** that carries all four answers with `describe()` and a JSON contract.

### 3.3 API (new, in `pytex.core.transformation` unless noted)

```python
@dataclass(frozen=True, slots=True)
class ORParallelismStatement:
    """One '(hkl)_p || (hkl)_c' or '[uvw]_p || [uvw]_c' clause of an OR statement."""
    kind: str                      # "plane" | "direction"
    parent_indices: np.ndarray     # (3,) int
    child_indices: np.ndarray      # (3,) int
    deviation_deg: float
    parent_label: str              # from pytex.core.notation
    child_label: str

@dataclass(frozen=True, slots=True)
class ORCharacterizationReport:
    relationship: OrientationRelationship        # the fitted OR
    pair_count: int
    residuals_deg: np.ndarray                    # per-pair, symmetry-aligned
    catalog_names: tuple[str, ...]
    catalog_deviations_deg: np.ndarray           # symmetry-reduced distance to each catalog OR
    best_catalog_name: str | None
    best_catalog_deviation_deg: float
    margin_deg: float                            # runner-up minus winner
    plane_statements: tuple[ORParallelismStatement, ...]
    direction_statements: tuple[ORParallelismStatement, ...]
    provenance: ProvenanceRecord | None = None

    @property
    def mean_residual_deg(self) -> float: ...
    @property
    def is_conclusive(self) -> bool: ...          # margin > residual scatter and winner < tolerance
    def describe(self) -> str: ...
    def to_json_dict(self) -> dict[str, Any]: ...


def characterize_orientation_relationship(
    parent_orientations: OrientationSet,
    child_orientations: OrientationSet,
    *,
    catalog: OrientationRelationshipCatalog | tuple[OrientationRelationship, ...] | None = None,
    nominal: OrientationRelationship | None = None,
    parallelism_tolerance_deg: float = 3.0,
    max_index: int = 3,
    max_statements: int = 4,
    provenance: ProvenanceRecord | None = None,
) -> ORCharacterizationReport: ...


def orientation_relationship_from_euler(
    parent_euler_deg: ArrayLike,          # (n, 3) Bunge degrees
    child_euler_deg: ArrayLike,           # (n, 3) Bunge degrees
    parent_phase: Phase,
    child_phase: Phase,
    **kwargs,
) -> ORCharacterizationReport: ...


def describe_orientation_relationship(
    relationship: OrientationRelationship,
    *,
    tolerance_deg: float = 3.0,
    max_index: int = 3,
    max_statements: int = 4,
) -> tuple[ORParallelismStatement, ...]: ...
```

### 3.4 Algorithm

**Step 1 — starting estimate without a nominal.** With `n ≥ 1` pairs compute
`V_i = C_i^T P_i`. Reduce every `V_i` into the double-coset fundamental zone (choose, over
`S_c V_i S_p`, the representative with maximum trace). Take the eigen-mean of the reduced set as
the seed. This removes the requirement that the caller already knows the answer; when `nominal`
*is* supplied it is used as the seed instead (reproducing today's behavior exactly).

**Step 2 — refine.** Delegate to the existing `fit_orientation_relationship` machinery with the
seed relationship. Do not duplicate the align/average loop — refactor it into a private
`_fit_from_seed` used by both entry points (one-shared-helper rule).

**Step 3 — catalog ranking.** For each catalog member compute the symmetry-reduced angle between
the fitted rotation and the member (existing `_symmetry_reduced_angle_between_deg`). Rank; record
the winner and the margin. The default catalog is chosen by the parent/child crystal systems
(cubic→cubic ⇒ fcc/bcc set; cubic→hexagonal ⇒ bcc/hcp Burgers set; etc.), and this dispatch is a
single documented table, not scattered `if`s.

**Step 4 — parallelism extraction.** For the fitted rotation, run the existing parallelism finder
over low-index parent planes and directions (`max_index`), keep matches within
`parallelism_tolerance_deg`, sort by deviation then by index magnitude, and return the top
`max_statements` **non-redundant** clauses. Redundancy rule: a plane clause and a direction clause
are independent only if the direction is not the cross product of two already-reported plane
normals (and vice versa); rank preference goes to close-packed families. This yields the familiar
form `(111)_γ ‖ (011)_α`, `[10̄1]_γ ‖ [11̄1]_α`.

**Step 5 — honesty.** `is_conclusive` is False whenever the winner's margin over the runner-up is
smaller than the pair scatter, or the winner deviation exceeds `parallelism_tolerance_deg`.
`describe()` must **say so in words**, per the explainable-results doctrine. A single measured pair
can never be conclusive about scatter; `describe()` states that too.

### 3.5 Validation

- **Round trip:** synthesize `n` child orientations from a known parent through known variants of
  KS/NW/GT/Pitsch/Burgers, add controlled scatter, recover; assert the fitted rotation is within
  the analytic bound and the catalog winner is the planted OR.
- **Literature identity:** KS recovers `(111)_γ ‖ (011)_α` and `[10̄1]_γ ‖ [11̄1]_α` at 0 deg;
  NW recovers `(111)_γ ‖ (011)_α` and `[11̄2]_γ ‖ [01̄1]_α` at 0 deg;
  Burgers recovers `(0001)_α ‖ (011)_β` and `[112̄0]_α ‖ [11̄1]_β` at 0 deg.
  These are definitional identities, not copied program output.
- **Discrimination:** KS vs NW are 5.26 deg apart; the report must separate them from data with
  scatter well below that and must declare *inconclusive* when scatter approaches it.
- **Non-degeneracy:** a random rotation is reported as matching no catalog member.

---

## 4. TX2 — Variant-resolved parallel directions and planes

### 4.1 The scientific question

"For this OR, if I have `[uvw]` (or `(hkl)`) in the parent, what is the parallel direction (plane)
in each of the product variants?" — the everyday tool for interpreting trace analysis, habit
directions, and diffraction spot correspondence.

### 4.2 What exists and what is missing

Existing: `map_direction_to_child` / `map_plane_to_child` (single variant),
`map_direction_across_variants` / `map_plane_across_variants` (tuples of correspondences),
`find_parallel_planes` / `find_parallel_directions` (search over families).

Missing: a single **table** object over one *or many* input indices, both mapping senses, with
grouping of variants that give crystallographically equivalent answers, exact-vs-rationalized
residual reporting, `describe()`, and CSV/JSON/Markdown export. Users currently have to loop and
format by hand, which is exactly what a library should absorb.

### 4.3 API (new, in `pytex.core.transformation`)

```python
@dataclass(frozen=True, slots=True)
class VariantCorrespondenceRow:
    variant_index: int
    source_indices: np.ndarray        # (3,) int, parent (or child) input
    exact_components: np.ndarray      # (3,) float, unrationalized image
    indices: np.ndarray               # (3,) int, rationalized image
    residual_deg: float               # angle between exact image and rationalized indices
    source_label: str
    image_label: str
    equivalence_group: int            # variants sharing a symmetry-equivalent image

@dataclass(frozen=True, slots=True)
class VariantCorrespondenceTable:
    relationship_name: str
    kind: str                          # "plane" | "direction"
    direction: str                     # "parent_to_child" | "child_to_parent"
    rows: tuple[VariantCorrespondenceRow, ...]
    max_index: int
    provenance: ProvenanceRecord | None = None

    def rows_for(self, source_indices) -> tuple[VariantCorrespondenceRow, ...]: ...
    def rows_for_variant(self, variant_index: int) -> tuple[VariantCorrespondenceRow, ...]: ...
    def distinct_image_count(self) -> int: ...
    def to_records(self) -> list[dict[str, Any]]: ...
    def to_csv(self, path) -> Path: ...
    def to_markdown(self) -> str: ...
    def to_json_dict(self) -> dict[str, Any]: ...
    def describe(self) -> str: ...


def variant_correspondence_table(
    relationship: OrientationRelationship,
    *,
    planes: Sequence[CrystalPlane] | None = None,
    directions: Sequence[CrystalDirection] | None = None,
    variants: tuple[TransformationVariant, ...] | None = None,
    sense: str = "parent_to_child",
    max_index: int = 6,
    provenance: ProvenanceRecord | None = None,
) -> VariantCorrespondenceTable | tuple[VariantCorrespondenceTable, ...]: ...
```

`equivalence_group` is computed by reducing each rationalized image to its symmetry-canonical
family representative under the *image phase* point group; variants sharing a representative share
a group id. This directly answers "how many crystallographically distinct answers are there really"
— for KS, a `{111}_γ` plane gives one distinct `{011}_α` answer across all 24 variants, which the
table must show as one equivalence group.

### 4.4 Validation

- KS: `(111)_γ` maps to `(011)_α` with 0 deg residual in every variant (definitional).
- Burgers: `(0001)_α` ↔ `(011)_β` and `[112̄0]_α` ↔ `[11̄1]_β` at 0 deg (definitional).
- Round trip: `child_to_parent` of the `parent_to_child` image returns the input within the
  rationalization residual.
- Rationalization sanity: raising `max_index` never increases the residual.
- Export round trip: `to_records` → CSV → re-read reproduces indices exactly.

---

## 5. TX3 — Robust, exportable composite SAED from a parent zone axis

### 5.1 What exists

`simulate_composite_saed` (parent-anchored, all/selected variants, exact irrational child zones with
nearest-rational labels), `CompositeSAEDPattern` with `describe()`, `find_spot_coincidences`,
`sweep_parent_zone_axes`, and the renderer/annotation engine in `pytex.plotting.composite_saed`.

### 5.2 Gaps to close (robustness)

**Corrected 2026-08-04 after auditing the live engine.** Three of the five gaps this section
originally listed were already closed, and the specification was wrong to claim otherwise:

- *Spot ordering* is already fully deterministic — `simulate_zone_axis_spots` ends in
  `np.lexsort((l, k, h, radius, -intensity))`, and `SpotTable`'s docstring states it.
- *Guard rails* already exist: `KinematicSimulationConfig.__post_init__` validates every
  numeric field, `simulate_zone_axis_spots` rejects a phase mismatch and validates the
  detector basis, and `CompositeSAEDPattern` rejects duplicate variant indices.
- *Intensity normalization* is already an explicit, documented decision rather than an
  accident: each sub-pattern is max-normalized, and `describe()` states why — kinematic
  cross-phase intensity ratios are not defined at this level of theory. Adding a "shared"
  normalization option would manufacture a number the theory does not support, so it is
  **rejected**, not deferred.

The genuine gaps are:

1. **Systematic absence audit.** `ReflectionCondition.from_phase` falls back to primitive
   (`P`) when a phase carries no space-group symbol, so a phase lacking that metadata
   silently keeps reflections that its real centering would forbid. The composite report and
   the export manifest must state, per phase, whether the centering was *declared* or
   *assumed*.
2. **`describe()` completeness** — should name the centering actually applied per phase
   alongside the existing camera constant and kinematic-only caveat.

### 5.3 Export layer (new module `pytex.diffraction.export`)

```python
@dataclass(frozen=True, slots=True)
class CompositeSAEDExport:
    directory: Path
    reflection_table_path: Path
    coincidence_table_path: Path | None
    manifest_path: Path
    figure_paths: tuple[Path, ...]
    def describe(self) -> str: ...


def composite_reflection_table(
    pattern: CompositeSAEDPattern,
    *,
    intensity_threshold: float = 0.0,
) -> ReflectionTable: ...                 # tabular object with to_csv/to_markdown/to_records


def export_composite_saed(
    pattern: CompositeSAEDPattern,
    directory,
    *,
    stem: str = "composite_saed",
    figure_formats: tuple[str, ...] = ("svg",),
    include_coincidences: bool = True,
    coincidence_tolerance_mm: float = 0.05,
    plot_config: CompositeSAEDPlotConfig | None = None,
) -> CompositeSAEDExport: ...
```

The reflection table has one row per rendered spot with: source (`parent` or `variant N`), phase
name, `h k l`, formatted label, `d` (Å), `|g|` (Å⁻¹), detector `x`, `y`, `r` (mm), excitation error,
`|F|`, relative intensity, and the zone axis it belongs to. The manifest is JSON against a new
`schemas/composite_saed_manifest.schema.json` recording relationship, phases, zone axes (exact and
rationalized), configuration, variant selection, and file inventory — the data-contract rule for
anything crossing a tool boundary.

### 5.4 Validation

- Table row count equals rendered spot count for a pinned case; every row's `d` matches
  `1/|g|` to 1e-12.
- Friedel symmetry: for a centrosymmetric phase the table is symmetric under `g → -g` in both
  position and intensity.
- Manifest validates against its schema; CSV re-read reproduces detector coordinates to 1e-9 mm.
- Figure export produces non-empty SVG containing the expected variant labels (structural
  assertion, not a byte baseline — per the repository's SVG policy).

---

## 6. TX4 — Composite patterns anchored on a product-variant zone axis

### 6.1 The scientific question

In practice one tilts to a **low-index zone of the product** (e.g. `[0001]_α` of a particular
martensite/α variant), and then wants to know what the parent and the *other* variants contribute
to that same pattern. TX3 only supports anchoring on the parent zone.

### 6.2 API (new, in `pytex.diffraction.composite`)

```python
def simulate_composite_saed_from_child_zone(
    relationship: OrientationRelationship,
    child_zone_axis: ZoneAxis,
    *,
    anchor_variant_index: int = 1,
    variant_indices: tuple[int, ...] | list[int] | None = None,
    include_parent: bool = True,
    config: KinematicSimulationConfig | None = None,
    child_config: KinematicSimulationConfig | None = None,
    align_child_g: MillerIndex | None = None,
    in_plane_rotation_deg: float = 0.0,
    rationalize_max_index: int = 6,
    provenance: ProvenanceRecord | None = None,
) -> CompositeSAEDPattern: ...
```

**Algorithm.** The anchor variant's rotation `R_k` maps parent Cartesian to child Cartesian. The
requested child zone `z_c` therefore corresponds to the exact parent direction
`z_p = R_k^T z_c` (generally irrational). Build the shared detector basis from `z_p` — reusing the
*same* `zone_basis_from_axis` path as TX3 so there is exactly one geometry definition — then
delegate to `simulate_composite_saed`. The in-plane alignment reference `align_child_g` is mapped
to the parent frame before basis construction so that "put this child reflection along +u" works
in the child's own indices.

`CompositeSAEDPattern` gains an `anchor` field (`"parent"` or `("variant", k)`) so `describe()`
and every export state which crystal defined the geometry. The parent zone axis is reported both
exactly and as its nearest rational label, mirroring how TX3 reports child zones.

### 6.3 Validation

- **Consistency identity:** anchoring on variant `k`'s image of a parent zone `z_p` must reproduce
  the TX3 pattern for `z_p` exactly (same detector coordinates to 1e-12 mm). This is the strongest
  available test and is definitional, not empirical.
- Anchoring on `[0001]_α` for Burgers must place the parent `{011}_β` reflections at the radii
  predicted by the β lattice parameter (independent analytic check).
- The anchor variant's own spot table must have an exactly rational zone axis with zero
  rationalization residual.

---

## 7. TX5 — Solving a measured SAED pattern from picked spots

### 7.1 The scientific question

Given a measured pattern — a set of spot positions relative to the transmitted beam, in pixels or
mm or Å⁻¹ — determine: the phase, the zone axis, the indexing of every spot, the crystal
orientation, and (when an OR is supplied) which variant. Report residuals and alternatives.

### 7.2 Input contract: the measured-pattern YAML

A single documented file format, validated against
`schemas/measured_saed_pattern.schema.json`:

```yaml
schema: pytex.measured_saed_pattern/1
name: burgers_alpha_zone_01
calibration:
  # exactly one of these three blocks
  camera_constant_mm_angstrom: 20.5        # L*lambda; r_mm = camera_constant * |g|
  # or: {camera_length_mm: 800, beam_energy_kev: 200}
  # or: units: reciprocal_angstrom          (coordinates already in 1/A)
  pixel_size_mm: 0.014                      # required when coordinates are in px
  centre_px: [512.0, 512.0]                 # transmitted beam position
spots:
  - {x: 612.0, y: 498.0, intensity: 1.0, label: A}
  - {x: 545.0, y: 631.0, intensity: 0.8, label: B}
candidates:
  phases: [alpha_zr, beta_zr]               # named fixtures or inline phase definitions
  max_index: 4
tolerances:
  length_relative: 0.03                     # |g| relative tolerance
  angle_deg: 2.0
relationship:                               # optional; enables variant assignment
  name: burgers_bcc_hcp
```

### 7.3 API (new module `pytex.diffraction.solving`)

```python
@dataclass(frozen=True, slots=True)
class MeasuredSpot:
    position_mm: np.ndarray        # (2,) relative to the transmitted beam
    intensity: float | None
    label: str | None
    @property
    def g_magnitude_inv_angstrom(self) -> float: ...

@dataclass(frozen=True, slots=True)
class MeasuredSAEDPattern:
    name: str
    spots: tuple[MeasuredSpot, ...]
    calibration: PatternCalibration
    provenance: ProvenanceRecord | None = None
    @classmethod
    def from_yaml(cls, path) -> MeasuredSAEDPattern: ...
    def to_yaml(self, path) -> Path: ...
    def describe(self) -> str: ...

@dataclass(frozen=True, slots=True)
class SolvedSpot:
    measured_index: int
    indices: np.ndarray            # (3,) int hkl
    label: str
    predicted_position_mm: np.ndarray
    residual_mm: float

@dataclass(frozen=True, slots=True)
class PatternSolution:
    phase_name: str
    zone_axis: ZoneAxis
    zone_axis_label: str
    orientation_in_pattern_frame: Rotation      # crystal -> detector frame
    solved_spots: tuple[SolvedSpot, ...]
    unindexed_indices: tuple[int, ...]
    mean_residual_mm: float
    max_residual_mm: float
    matched_fraction: float
    variant_index: int | None                   # when a relationship was supplied
    score: float
    def describe(self) -> str: ...
    def to_json_dict(self) -> dict[str, Any]: ...

@dataclass(frozen=True, slots=True)
class PatternSolutionReport:
    pattern_name: str
    solutions: tuple[PatternSolution, ...]      # ranked, best first
    considered_phase_names: tuple[str, ...]
    def best(self) -> PatternSolution: ...
    def is_conclusive(self) -> bool: ...
    def describe(self) -> str: ...
    def to_json_dict(self) -> dict[str, Any]: ...


def solve_saed_pattern(
    pattern: MeasuredSAEDPattern,
    phases: Sequence[Phase],
    *,
    max_index: int = 4,
    length_tolerance_relative: float = 0.03,
    angle_tolerance_deg: float = 2.0,
    relationship: OrientationRelationship | None = None,
    max_solutions: int = 5,
    provenance: ProvenanceRecord | None = None,
) -> PatternSolutionReport: ...


def solve_saed_pattern_file(path, **kwargs) -> PatternSolutionReport: ...
```

### 7.4 Algorithm (two-spot seed, all-spot verification)

1. **Calibrate.** Convert every picked spot to a reciprocal-space vector magnitude
   `|g| = r_mm / (L·λ)` and an in-plane azimuth. The detector plane is the `(u, v)` plane; the
   beam direction is `+w`.
2. **Seed pairs.** Choose the two shortest non-collinear spots (and, for robustness, the best
   `k` such pairs). For each candidate phase enumerate allowed reflections up to `max_index`,
   filtered by the phase's reflection conditions.
3. **Match.** A candidate `(g_1, g_2)` assignment is admissible if
   `| |g_i^calc| - |g_i^obs| | / |g_i^obs| ≤ length_tolerance_relative` for both, and the
   interplanar angle agrees within `angle_tolerance_deg`. This is the classical ratio/angle
   indexing test; it is vectorized as a single pairwise comparison over the reflection list.
4. **Zone axis.** `z = g_1 × g_2` in reciprocal space → direct-lattice zone axis via the metric
   tensor; rationalize with the existing `rationalize_zone_axis`.
5. **Orientation.** Build the crystal→detector rotation that carries `g_1`, `g_2` onto their
   observed in-plane directions (right-handed Gram–Schmidt on the calculated pair, matched to the
   observed pair). The handedness ambiguity (`z` vs `-z`) is resolved by which sign indexes more
   spots, and when both index equally the report says the pattern **cannot** distinguish them —
   a genuine and well-known property of a single SAED pattern that the tool must not paper over.
6. **Verify.** Project every allowed zone reflection to the detector and assign each measured spot
   to its nearest prediction within tolerance. Score by matched fraction first, then mean residual.
7. **Variant assignment.** When a relationship is supplied, compare the solved crystal→detector
   rotation against each variant's prediction (given the parent orientation implied by the
   solution, or given a supplied parent) and report the best-fitting variant index with its
   deviation.
8. **Rank and report.** Return up to `max_solutions` ranked solutions; `is_conclusive` requires
   the best score to lead the runner-up by a stated margin *and* the runner-up not to be a mere
   symmetry equivalent of the best.

### 7.5 Interactive picking (new, in `pytex.plotting.saed_picker`)

```python
class SAEDSpotPicker:
    """Matplotlib click-to-pick front end that produces a MeasuredSAEDPattern."""
    def __init__(self, image=None, *, calibration=None, extent=None, existing=None): ...
    def connect(self) -> SAEDSpotPicker: ...       # left click add, right click remove
    def pattern(self, *, name: str) -> MeasuredSAEDPattern: ...
    def save_yaml(self, path) -> Path: ...
```

The picker's *state machine* (add/remove/centre-set/undo) is a plain, fully-tested object;
Matplotlib event callbacks are a thin adapter over it, so the logic is testable headlessly and the
GUI is not on the critical path. `pattern()` works identically whether spots were clicked or read
from YAML — the file contract is the boundary, so a solved pattern is fully reproducible from a
committed text file.

### 7.6 Validation

- **Synthetic closure:** take a `simulate_zone_axis_spots` pattern for a known phase/zone,
  export its detector coordinates as a `MeasuredSAEDPattern`, solve, and recover the zone axis
  exactly and every spot's indices up to the symmetry/handedness ambiguity.
- **Noise robustness:** the same test with Gaussian position noise; recovery must hold to a stated
  noise level and must degrade to *inconclusive*, not to a confidently wrong answer.
- **Discrimination:** an fcc pattern must not be solved as bcc (systematic absences differ), and a
  solve run with both phases as candidates must rank the true one first.
- **Analytic pin:** for cubic `[001]`, the `{200}` and `{220}` spot radius ratio is exactly `√2`
  and the inter-spot angle exactly 45 deg — the solver's output is checked against these
  identities, not against prior program output.
- **Composite case:** a composite pattern from TX3 with parent + one variant, solved with the OR
  supplied, returns both phases and the correct variant index.
- **Picker:** the state machine's add/remove/undo/centre behavior is unit-tested; a YAML
  round trip reproduces every coordinate bit-for-bit.

---

## 8. Cross-cutting deliverables

| Deliverable | Requirement |
| --- | --- |
| Documentation | A concept page per feature under `docs/site/concepts/`, one workflow page tying (a)–(e) into an end-to-end story, API pages, and updates to this file plus the OR and phase-transformation foundations. |
| Theory | LaTeX notes in `docs/tex/` for the OR-statement extraction (TX1) and the ratio/angle indexing algorithm (TX5). |
| Figures | Canonical SVGs: OR statement geometry, variant correspondence fan, child-anchored composite geometry, and the SAED solving flow — all per the visualization style guide. |
| Worked examples | At least one executable worked example per feature, each with independent provenance (analytic identity or cited literature value). |
| Notebook | `23_burgers_transformation_crystallography.ipynb` — Burgers β↔α as the end-to-end teaching case covering (a)–(e), committed executed. |
| Schemas | `composite_saed_manifest.schema.json`, `measured_saed_pattern.schema.json`, `pattern_solution.schema.json`. |
| Validation ledger | New rows in `docs/testing/phase_transformation_validation_matrix.md` and `docs/testing/diffraction_validation_matrix.md`. |
| Terminology | New symbols registered in `docs/standards/terminology_and_symbol_registry.md` before use. |

## 9. Phasing

The ordering is chosen so each phase is independently useful and testable, and so later phases
consume earlier ones rather than anticipating them.

| Phase | Scope | Depends on |
| --- | --- | --- |
| TX0 | This specification + ledger + schema stubs | — |
| TX1 | OR characterization from measured orientations | — |
| TX2 | Variant correspondence tables | — |
| TX3 | Composite SAED robustness + export layer | — |
| TX4 | Child-zone-anchored composite patterns | TX3 |
| TX5a | Measured-pattern YAML contract + calibration + solver core | TX3 |
| TX5b | Variant assignment + interactive picker | TX5a, TX4 |
| TX6 | Burgers notebook, figures, theory notes, ledger closure | all |

Every phase ends in a verified commit on `main` with the standing gates green:
`pytest`, `check_repo_integrity.py`, `ruff`, `mypy src`, and the Sphinx build.

## References

### Normative

- [Canonical Data Model](canonical_data_model.md)
- [Orientation Relationship Analysis Foundation](orientation_relationship_analysis_foundation.md)
- [Phase Transformation Foundation](phase_transformation_foundation.md)
- [Diffraction Foundation](diffraction_foundation.md)
- [Notation And Conventions](../standards/notation_and_conventions.md)
- [Data Contracts And Manifests](../standards/data_contracts_and_manifests.md)
- [Executable Examples](../standards/executable_examples.md)

### Informative

- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*, 3rd ed.
- Williams and Carter, *Transmission Electron Microscopy*, 2nd ed. — SAED indexing and camera
  constant calibration.
- Edington, *Practical Electron Microscopy in Materials Science*, Monograph 2 — ratio/angle
  indexing of single-crystal patterns.
- Burgers, *Physica* 1 (1934) 561 — the bcc↔hcp orientation relationship.
- Morito et al., *Acta Materialia* 51 (2003) 1789 — KS variant numbering and packet structure.
