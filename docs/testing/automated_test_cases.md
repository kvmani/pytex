# Automated Test Cases

PyTex uses automated tests as a scientific audit surface, not merely as a CI gate.

This page is written for domain experts who want to inspect whether a method, formula, or
convention is scientifically credible without first reading Python implementation details.

## What This Page Must Do

For every major crystallographic conversion, orientation convention, or diffraction calculation,
the validation surface should make five things visible in one place:

- the authoritative reference and exact page anchor
- the source notation and governing formula
- the algorithm currently implemented by PyTex
- the source-derived expected value for a worked example
- the actual value returned by the current code

The goal is not to decorate the test suite. The goal is to make the numerical pathway auditable by
a crystallographer, diffraction specialist, or EBSD expert who wants to challenge assumptions,
notation choices, or interpretation directly.

## Review Pattern

Each documented case on this page follows the same audit structure:

1. `Code surface`
   The method, class, and automated tests that define the behavior.
2. `Reference basis`
   Exact source and page numbers.
3. `Source formulation`
   The formula or convention written as close to the source notation as practical.
4. `Algorithm used by PyTex`
   A prose summary of the implementation path, without requiring code reading.
5. `Audit example`
   A worked example with source-derived expected values, current code outputs, and an
   interpretation column.

## Existing Cases

The current source-controlled audit cases remain:

- Euler conversion audit
- Hexagonal 4-index audit
- Reciprocal vector and interplanar spacing audit
- Powder XRD and Bragg-law audit
- Harmonic ODF reconstruction audit
- Foundation feature priority audit

## Case 6: Structure-Import Fixture Audit

### Scope

This case documents how PyTex validates the bundled CIF-backed starter corpus against pinned
space-group, point-group, lattice, and site-count expectations.

### Code surface

- `tests/unit/test_phase_fixtures.py`
- `benchmarks/structure_import/phase_fixture_audit_summary.json`
- `fixtures/phases/catalog.json`

### Reference basis

- IUCr crystallographic conventions
- pinned fixture metadata under `fixtures/phases/*/metadata.json`

### Algorithm used by PyTex

PyTex currently:

1. verifies the hash-pinned catalog and metadata
2. loads each fixture through the optional `pymatgen`-backed CIF path in the full scientific lane
3. records conventional and primitive expectations in a pinned audit summary
4. checks that the loaded phases match the expected space-group, point-group, lattice, and site-count semantics

### Audit examples

Representative pinned expectations from the current audit summary:

| Fixture | Expected conventional semantics | Expected primitive semantics | Interpretation |
| --- | --- | --- | --- |
| `fe_bcc` | `Im-3m`, point group `m-3m`, site count `2`, `a=b=c=2.8665 A` | `Im-3m`, site count `1`, rhombohedral primitive cell angles `109.4712206345 deg` | Confirms that PyTex preserves both the conventional and primitive-cell interpretation for BCC iron. |
| `zr_hcp` | `P6_3/mmc`, point group `6/mmm`, site count `2`, `gamma=120 deg` | same as conventional | Confirms the pinned HCP starter case used for hexagonal-axis teaching and validation. |
| `ni_fcc` | `Fm-3m`, point group `m-3m`, site count `4`, `a=b=c=3.52387 A` | `Fm-3m`, site count `1`, primitive angles `60 deg` | Confirms the pinned FCC starter case that also feeds the diffraction baseline corpus. |

### Current claim boundary

- The full scientific lane justifies that the bundled fixture corpus reconstructs into canonical
  PyTex objects with pinned semantics.
- It does not yet justify broad parity against large external CIF corpora.

## Case 7: Transformation Starter-Family Audit

### Scope

This case documents the current literature-tracked starter transformation helpers and their tested
correspondence semantics.

### Code surface

- `tests/unit/test_transformation.py`
- `tests/unit/test_experimental_phase_transformation.py`
- `src/pytex/core/transformation.py`

### Reference basis

- Bain correspondence
- Nishiyama-Wassermann correspondence
- Porter, Easterling, Sherif, *Phase Transformations in Metals and Alloys*

### Algorithm used by PyTex

PyTex currently:

1. constructs named orientation relationships from explicit plane-direction correspondences
2. enforces cubic-phase guards for the named Bain and Nishiyama-Wassermann helpers
3. generates symmetry-derived variants and tests uniqueness
4. uses those variants in stable child-orientation prediction and bounded experimental parent scoring

### Audit examples

| Helper | Source-facing correspondence | Current tested behavior | Interpretation |
| --- | --- | --- | --- |
| `OrientationRelationship.from_bain_correspondence(...)` | `(001)_p || (001)_c`, `[110]_p || [100]_c` | the mapped parent `[110]` direction normalizes onto child `[100]`; non-cubic parents are rejected | This is a literature-tracked starter family, not yet a broad family-validation program. |
| `OrientationRelationship.from_nishiyama_wassermann_correspondence(...)` | `(111)_p || (011)_c`, `[1-10]_p || [100]_c` | the mapped parent `[1-10]` direction normalizes onto child `[100]`; non-cubic children are rejected | PyTex now records this family explicitly rather than leaving it as an ad hoc custom rotation. |

### Current claim boundary

- PyTex now has explicit, tested starter transformation families and variant semantics.
- It does not yet claim broad literature parity across transformation-family catalogs or full parent reconstruction workflows.

## Case 8: Orientation-Construction And Rodrigues Audit

### Scope

This case documents the parity-oriented constructor surface for quaternions, matrices, axis-angle,
 Euler angles, and Rodrigues / Rodrigues-Frank coordinates.

### Code surface

- `tests/unit/test_orientation.py`
- `src/pytex/core/orientation.py`
- `src/pytex/core/batches.py`

### Reference basis

- the quaternion, axis-angle, and Rodrigues relations documented in
  `../tex/algorithms/orientation_representations_and_plane_direction_construction.tex`

### Algorithm used by PyTex

PyTex currently:

1. normalizes every quaternion-backed constructor onto the canonical unit-quaternion surface
2. maps Rodrigues coordinates through the explicit tangent-half-angle relation
3. preserves the Rodrigues-Frank singularity at `omega = pi` as `+inf` rather than clipping it
4. checks that scalar and batch constructors recover the same rotation matrices

### Current claim boundary

- The stable constructor surface is now explicit and test-backed for Euler, matrix, quaternion,
  axis-angle, Rodrigues, and Rodrigues-Frank inputs.
- It does not yet claim exhaustive external parity for every MTEX constructor idiom or sampling helper.

## Case 9: Miller Family And Hexagonal-Index Audit

### Scope

This case documents the current Miller plane and direction family behavior, including four-index
hexagonal and trigonal conversions.

### Code surface

- `tests/unit/test_miller_objects.py`
- `tests/unit/test_orientation.py`
- `src/pytex/core/miller.py`
- `src/pytex/core/orientation.py`

### Algorithm used by PyTex

PyTex currently:

1. reduces Miller indices to canonical integer representatives
2. converts `hkl <-> hkil` and `uvw <-> UVTW` through the documented four-index relations
3. expands families by applying crystal symmetry in cartesian space and recovering integer indices
4. constructs orientations from Miller scalar/set objects and from raw three- or four-index arrays

### Current claim boundary

- The Miller surface now supports vectorized family, angle, projection, JSON-contract, orix-adapter,
  and orientation-construction workflows without ambiguous array-only helpers.
- It does not yet claim full parity for every MTEX plotting or crystal-class convenience workflow.

## Case 10: Fundamental-Region Reduction Audit

### Scope

This case documents class-by-class exact orbit reduction across the supported proper point groups.

### Code surface

- `tests/unit/test_orientation_utilities.py`
- `src/pytex/core/orientation.py`
- `src/pytex/core/symmetry.py`

### Algorithm used by PyTex

PyTex currently:

1. enumerates the `s g c` orbit under specimen and crystal symmetry
2. canonicalizes orbit candidates in the unit-quaternion hemisphere
3. selects the minimum-angle representative relative to the identity
4. resolves boundary ties deterministically through lexicographic quaternion keys

### Current claim boundary

- Exact orbit reduction is now explicitly regression-tested across the supported proper point groups
  already generated in the codebase.
- Closed-form external boundary catalogs and broader external parity fixtures still remain ahead.

## Case 11: Harmonic Reconstruction Residual Audit

### Scope

This case documents the retained-basis and residual-diagnostic surface for harmonic PF-to-ODF
reconstruction.

### Code surface

- `tests/unit/test_harmonic_odf.py`
- `src/pytex/texture/harmonics.py`
- `src/pytex/texture/reconstruction.py`

### Algorithm used by PyTex

PyTex currently:

1. builds a symmetry-projected raw harmonic basis on a Bunge quadrature grid
2. orthonormalizes the retained basis after dropping numerically rank-deficient modes
3. solves the regularized least-squares reconstruction problem
4. reports residual norms, matrix rank, retained-basis size, symmetry orders, and density diagnostics

### Current claim boundary

- The harmonic reconstruction surface now exposes auditable residual and retained-basis diagnostics
  in addition to reconstruction outputs.
- It does not yet claim full experimental correction doctrine or broad external harmonic parity.

## Current Audit Findings

The documented cases above already show issues that are scientifically important even when the
current tests pass:

- The Bunge matrix in the reference texts and the matrix surfaced by PyTex differ by a transpose
  because they describe different mapping conventions. This is not a bug by itself, but it must be
  documented explicitly.
- The current `abg` or `matthies` export agrees numerically with the Roe-style relation for the
  worked orientation example. It should not be described as Kocks unless PyTex adds a distinct
  Kocks conversion and validates it.
- The powder-XRD angle pathway is source-auditable and strong. The intensity pathway is
  intentionally pedagogical and should stay labeled as such.
- The bundled structure-fixture corpus now supports stronger reproducibility claims than
  parser-success testing alone, but it is still a starter corpus rather than a broad external
  parity suite.
- The transformation layer now supports literature-tracked starter correspondences, but it remains
  foundational until broader family breadth and curated datasets land.
- The expanded orientation constructor surface is scientifically explicit and internally consistent,
  but parity status should remain `foundational` until broader external comparison cases land.

## Maintenance Rule

When a future task changes a formula, convention, or public numerical pathway:

- update the relevant automated tests
- update the corresponding audit case on this page
- refresh the source-derived expected values if the reference basis changes
- refresh the current code output block from a real local run
- describe any mismatch as a scientific interpretation issue, not only as a pass or fail test
  event

## References

### Normative

- <a href="testing_strategy.html">Testing Strategy</a>
- <a href="../standards/documentation_architecture.html">Documentation Architecture</a>
- <a href="../standards/development_principles.html">Development Principles</a>

### Informative

- `references/formulation_summary.md`
- `references/reference_index.md`
