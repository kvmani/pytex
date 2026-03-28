# Development Principles

This document makes the repository's engineering posture explicit so future stages do not drift in style or semantics.

## Core Principles

1. Semantic correctness before feature breadth
   New capability must extend the canonical data model when scientific meaning would otherwise be implicit.

1. Vectorized semantics are stable API surface
   If a large-data workflow depends on shared frame, convention, symmetry, or provenance meaning, the public API should use a semantic batch primitive rather than a raw ndarray alone.

2. Normalize once at the boundary
   Vendor, adapter, and file-format conventions should be converted into PyTex canonical semantics once and then carried explicitly.

3. Fail fast on inconsistent metadata
   Mismatched frames, phases, symmetries, units, and notation should be rejected at construction time wherever practical.

4. Centralize conventions
   Frame rules, indexing conventions, symmetry choices, and tolerance policy belong in shared standards documents and core utilities.

5. Separate stable and experimental work
   `src/pytex/experimental/` may explore new methods, but it must not redefine stable semantics or weaken validation requirements.

6. Document with the implementation
   Stable scientific behavior changes require code, tests, Markdown guidance, LaTeX theory or algorithm notes, explicit mathematical exposition where conventions or algorithms are involved, and figures when geometry matters.

7. Validate before claiming equivalence
   PyTex should not claim parity with MTEX or other tools until the relevant behavior is covered by automated tests or benchmark notes.

8. Optimize after semantics are fixed
   Performance work is welcome, but only after the canonical representation, failure modes, and validation criteria are clear.

9. Treat machine-readable interchange as product surface
   Major stable outputs should not remain trapped in Python memory only. When a result is important for reproducibility, external tooling, validation, or workflow chaining, it should gain a versioned JSON contract that preserves enough scientific meaning for reconstruction.

## Definition Of Done For Stable Features

A stable feature is not complete until all of the following exist:

- explicit domain types or clearly justified API boundaries
- semantic batch types for vectorized workflows where shared metadata affects interpretation
- construction-time invariant checks
- deterministic tests
- a benchmark or parity note
- an example workflow
- an executable notebook tutorial when the feature benefits from staged interactive exposition
- a reusable plotting path when the feature naturally produces scientific figures
- citations and scientific documentation
- a canonical JSON contract for major outputs that are intended to cross workflow or tool boundaries
- explicit mathematical definitions for major conventions, mappings, reductions, and algorithms
- SVG figures when frames, vectors, geometry, or reduction logic matter

## Documentation Escalation Rule

The repository default is no longer “text first, figures later” for foundational scientific behavior.

- Major conventions must be defined with explicit mathematical meaning, not just prose summaries.
- Frame mappings, symmetry actions, reduction rules, and projection conventions must be illustrated with annotated graphics suitable for both teaching and review.
- If a change alters a major algorithm, the corresponding docs update must explain the mathematical object being transformed, the reduction or selection rule, and the current implementation boundary.
- Notebook tutorials should agree with the formal docs and must not introduce private conventions, hidden helper assumptions, or unexplained magic numbers.
- Runtime plotting should reuse shared semantic builders rather than inventing ad hoc Matplotlib logic per workflow.
- A code change that improves scientific capability but leaves the foundational documentation less clear than the implementation is incomplete.
- Stable terms and symbols should be added to the shared terminology registry and linked from affected docs rather than redefined independently.

## Review Questions

Before merging a substantial change, answer these questions explicitly:

- Does it introduce a new scientific convention?
- Does it bypass or duplicate an existing core primitive?
- Does it preserve provenance and source semantics?
- Does it overstate validation or maturity?
- Does it make the next stage easier rather than more ad hoc?

## References

### Normative

- `engineering_governance.md`
- `reference_canon.md`

### Informative

- `../roadmap/implementation_roadmap.md`
