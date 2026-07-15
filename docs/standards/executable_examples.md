# Executable Worked Examples

This standard defines PyTex's **documentation-as-test** discipline: the rule that
the numerical examples in the documentation are computed live from the code and
checked against independently known reference values. It is normative for every
stable public numerical surface.

The goal is a documentation system that teaches, that demonstrates the API, and
that fails loudly when the code and the docs disagree — all from a single source
of truth.

## Motivation

A scientific library earns trust when a reader can see a claim, see the exact
code that produces it, and know that the printed number came from that code and
not from a stale copy-paste. Hand-written example outputs rot: a refactor changes
a result and the documentation silently lies. PyTex removes that failure mode by
making every tabulated value the output of a real call into the public API,
regenerated on every documentation build and asserted on every test run.

The canonical illustration is the angle between two crystallographic vectors. The
documentation does not merely *assert* that the angle between `(100)` and `(110)`
in a cubic crystal is 45 degrees; it *computes* it with
`angle_plane_plane_rad(...)`, prints the result, and compares it to the analytic
identity `arccos(1/sqrt(2)) = 45` — so the same artifact teaches the geometry,
demonstrates the function, and guards the implementation.

## The Worked-Example Contract

A worked example is a `WorkedExample` object (see the `worked_examples/` package)
that carries all of the following, together, in one place:

- **id** — a stable, unique identifier used in tests and generated pages.
- **scenario** — prose stating *when and why* a user computes this quantity, and
  *where in a real workflow* the corresponding module, class, or method is used.
- **setup + code** — one Python source that is *both* rendered verbatim in the
  documentation *and* executed to bind a `result` variable. The displayed code is
  the executed code; there is no separate, drift-prone transcript.
- **expected** — an independently known reference value: an analytic identity, an
  International-Tables convention, or a cited textbook/standard value. It must not
  be copied from a previous program run.
- **unit** and **tolerance** — the physical unit and the absolute tolerance of the
  computed-versus-expected comparison.
- **reference** — the derivation or provenance of the expected value.
- **citation** — a normative or informative source, per the
  [Scientific Citation Policy](scientific_citation_policy.md).
- **symbols** — the registry symbols the example uses, kept consistent with the
  [Terminology and Symbol Registry](terminology_and_symbol_registry.md).
- **see_also** — navigable cross-links to the theory note, API surface, and
  concept page, per the [Documentation Architecture](documentation_architecture.md).

## The Three Surfaces Rule

Each worked example must drive all three surfaces from the same object:

1. **Test** — `tests/unit/test_worked_examples.py` runs every example and asserts
   `computed` matches `expected` within `tolerance`. A divergence is a test
   failure, not a documentation nit.
2. **Documentation** — `scripts/generate_worked_examples.py` renders every example
   into the Sphinx gallery under `docs/site/examples/`, showing the computed value
   and the expected value side by side, with the scenario, code, provenance, and
   citation. Generated pages carry a "do not edit by hand" banner and are
   byte-stable, so a staleness test can detect un-regenerated docs.
3. **Notation** — every example references symbols fixed in the terminology and
   symbol registry, which keeps one nomenclature across prose, mathematics,
   figures, notebooks, and code.

## When A Worked Example Is Required

- Every **stable public numerical function or method** whose result a user could
  reasonably want to verify (angles, spacings, multiplicities, misorientation
  angles, structure factors, diffraction angles, tensor invariants, texture index
  values, and similar) must have at least one worked example.
- Any change that alters such a numerical result must update the worked example in
  the same change, and the new expected value must retain independent provenance.
- Experimental surfaces (`src/pytex/experimental/`) may add worked examples but are
  not obligated to; when they do, the scenario must state the experimental status.

## Choosing Expected Values

Expected values must be defensible without running PyTex:

- **Analytic identities** are preferred where a closed form exists for the crystal
  system (for example, cubic interplanar angles, or the perpendicularity of the
  basal and prismatic normals in hexagonal metals for any `c/a`).
- **Standard reference values** (International Tables multiplicities, ICDD powder
  positions, established coincidence-site boundary angles) are used where an
  analytic identity is unwieldy.
- A value copied from a prior PyTex run is **not** an acceptable expected value; it
  would only test that the code agrees with itself.

## Tolerances

- Tolerances are absolute and explicit per example.
- Exact integer or exact-identity results use a tolerance at or near floating-point
  round-off (for example `1e-9`), or `0.0` for integer-valued results.
- Reference values with rounding or physical-constant uncertainty (for example a
  tabulated powder angle) use a stated, justified tolerance.
- Tolerance policy is consistent with the
  [Benchmark and Tolerance Governance](benchmark_and_tolerance_governance.md).

## Relationship To Other Documentation Surfaces

Worked examples do not replace the deeper surfaces; they complete them:

- **LaTeX theory and algorithm notes** remain the canonical derivations.
- **Notebook tutorials** remain the staged, interactive narratives.
- **SVG figures** remain the canonical geometry illustrations.
- Worked examples are the **verifiable numerical bridge** between the theory and
  the running code, and they are part of the definition of done for stable
  numerical features (see [Development Principles](development_principles.md)).

## Authoring Checklist

Before a numerical feature is considered documented:

- [ ] A `WorkedExample` exists in the appropriate `worked_examples/examples/` module.
- [ ] The `code` binds `result` and the `setup` uses only the public `pytex` API.
- [ ] The `expected` value has independent provenance and a citation.
- [ ] The scenario states the real use case and where the code is applied.
- [ ] The symbols used are present in the terminology and symbol registry.
- [ ] `python scripts/generate_worked_examples.py` was run and the generated pages
      committed.
- [ ] `python -m pytest tests/unit/test_worked_examples.py` passes.

## References

### Normative

- [Documentation Architecture](documentation_architecture.md)
- [Terminology and Symbol Registry](terminology_and_symbol_registry.md)
- [Development Principles](development_principles.md)
- [Scientific Citation Policy](scientific_citation_policy.md)
- [Benchmark and Tolerance Governance](benchmark_and_tolerance_governance.md)

### Informative

- The `worked_examples/` package (framework, registry, and example modules).
- `scripts/generate_worked_examples.py` (documentation rendering).
- `tests/unit/test_worked_examples.py` (agreement and staleness tests).
