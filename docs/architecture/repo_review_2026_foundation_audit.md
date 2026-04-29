# Repo Review: 2026 Foundation Audit

This memo records the current state of the repository after the first repo-wide foundation
hardening pass and now serves as the handoff into the validation-and-interoperability hardening
phase.

## Overall Assessment

PyTex is materially ahead of a normal pre-alpha repository.

Strengths:

- clear canonical frame, symmetry, and orientation semantics
- unusually mature documentation and governance posture for this stage
- real test coverage across core, texture, EBSD, and diffraction foundations
- explicit parity doctrine instead of vague compatibility claims

Primary risks:

- validation breadth is uneven across subsystems
- optional scientific dependency expectations are still too implicit unless the contributor follows
  the right install path
- multimodal delivery depth is still shallower than the long-horizon mission breadth
- some future algorithm families remain only foundational rather than externally benchmarked

## Review Buckets

### Doc Drift

- README previously understated the implemented EBSD and diffraction breadth.
- Specifications previously described the Sphinx surface as planned even though it is live.
- Roadmap prose previously mixed chronological phases with current-state claims in a way that
  obscured what is already implemented.
- Contributor setup docs previously blurred the lightweight base environment and the heavier full
  scientific lane.

### Foundation Gaps

- validation and benchmark evidence still lag the now-stable structure, manifest, acquisition, and
  transformation semantics
- multimodal acquisition semantics are implemented, but broader modality-specific workflow adoption
  remains ahead
- phase transformation now has stable primitives, but not yet a strong literature-backed validation
  program

### Validation Gaps

- texture and EBSD have an explicit MTEX-backed ledger
- diffraction has strong unit coverage but a weaker external-baseline doctrine
- structure import still needs a more explicit validation story than parser-success tests
- phase transformation needs a broader literature-tracked evidence surface before stronger family
  claims are made

## Subsystem Scorecard

| Subsystem | Semantics | Validation | Documentation | Readiness | Notes |
| --- | --- | --- | --- | --- | --- |
| Core model | strong | strong | strong | foundationally ready | Best current asset in the repo. |
| Orientation and texture | strong | strong | strong | foundationally ready | Needs later harmonic and exact-boundary expansion, not semantic rescue. |
| EBSD | strong | strong | strong | foundationally ready | Regular-grid workflow base is real; richer modality scope is still ahead. |
| Diffraction | strong foundation | moderate | strong | foundational | Geometry, powder XRD, and SAED workflows are real, but external-baseline doctrine is thinner than the implemented surface. |
| Scientific visualization | strong foundation | moderate | strong | foundational | Shared plotting, YAML styles, and 3D crystal viewing now exist, but publication presets and external visual regression are still thinner than the core semantics. |
| CIF and structure import | strong | moderate | strong | foundational | Space-group semantics and a dedicated validation ledger now exist; broader external validation is still ahead. |
| Multimodal acquisition | strong | moderate | strong | foundational | Core primitives, experiment manifests, and workflow entry points now exist; broader modality depth is still ahead. |
| Phase transformation | foundational | weak | strong | foundational | Primitive family now exists; validation and algorithm breadth are still ahead. |

## Learnings

- The repository's strongest differentiator is semantic explicitness, not feature count.
- The current codebase already proves that PyTex can carry texture, EBSD, and diffraction logic on
  one shared model.
- The next major risk is adding code faster than foundational documents define the scientific
  contracts.
- "Validation-first" needs to mean more than MTEX parity once the repo moves deeper into
  diffraction and transformation work.

## Immediate Recommendations

1. Keep the base lane and full scientific lane both explicit and green in CI.
2. Grow the manifest family into richer experiment, benchmark, validation, and workflow-result use
   cases.
3. Build structure-import, diffraction, and transformation validation doctrine before broadening
   public claims in those areas.
4. Treat transformation algorithms as validation-first follow-on work on top of the new primitive
   family.

## References

### Normative

- [Architecture Overview](overview.md)
- [Canonical Data Model](canonical_data_model.md)
- [Reference Canon](../standards/reference_canon.md)

### Informative

- [Orientation And Texture Foundation](orientation_and_texture_foundation.md)
- [EBSD Foundation](ebsd_foundation.md)
- [Diffraction Foundation](diffraction_foundation.md)
