# Repo Review: 2026 Foundation Audit

This memo records the current state of the repository after the first repo-wide foundation hardening pass.

## Overall Assessment

PyTex is materially ahead of a normal pre-alpha repository.

Strengths:

- clear canonical frame, symmetry, and orientation semantics
- unusually mature documentation and governance posture for this stage
- real test coverage across core, texture, EBSD, and diffraction foundations
- explicit parity doctrine instead of vague compatibility claims

Primary risks:

- the multimodal ambition is broader than the currently explicit foundational doctrine
- validation breadth is uneven across subsystems
- some necessary future primitive families are defined only implicitly today

## Review Buckets

### Doc Drift

- README previously understated the implemented EBSD and diffraction breadth.
- Specifications previously described the Sphinx surface as planned even though it is live.
- Roadmap prose mixed chronological phases with current-state claims in a way that obscured what is already implemented.

### Foundation Gaps

- no stable `SpaceGroupSpec` layer distinct from current point-group-heavy semantics
- no shared acquisition, calibration, or uncertainty object family
- no stable phase-transformation architecture yet
- no cross-domain manifest family beyond the EBSD import manifest
- no single explicit multimodal frame-chain doctrine until this hardening pass

### Validation Gaps

- texture and EBSD have an explicit MTEX-backed ledger
- diffraction has strong unit coverage but a weaker external-baseline doctrine
- structure import still needs a more explicit validation story than parser-success tests

## Subsystem Scorecard

| Subsystem | Semantics | Validation | Documentation | Readiness | Notes |
| --- | --- | --- | --- | --- | --- |
| Core model | strong | strong | strong | foundationally ready | Best current asset in the repo. |
| Orientation and texture | strong | strong | strong | foundationally ready | Needs later harmonic and exact-boundary expansion, not semantic rescue. |
| EBSD | strong | strong | strong | foundationally ready | Regular-grid workflow base is real; richer modality scope is still ahead. |
| Diffraction | strong foundation | moderate | strong | foundational | Geometry and indexing scaffolding are real, but external-baseline doctrine is thinner. |
| CIF and structure import | moderate to strong | moderate | moderate to strong | foundational | Good architectural shape; still needs broader validation and space-group hardening. |
| Multimodal acquisition | weak | weak | now foundationally defined | planned | Docs now define the target, but stable APIs are ahead. |
| Phase transformation | weak | weak | now foundationally defined | planned | Keep out of stable APIs until the primitive family is implemented and validated. |

## Learnings

- The repository’s strongest differentiator is semantic explicitness, not feature count.
- The current codebase already proves that PyTex can carry texture, EBSD, and diffraction logic on one shared model.
- The next major risk is adding code faster than foundational documents define the scientific contracts.
- “Validation-first” needs to mean more than MTEX parity once the repo moves deeper into diffraction and transformation work.

## Immediate Recommendations

1. Keep lint, type checking, tests, and docs builds green in CI.
2. Treat multimodal acquisition, calibration, and uncertainty as first-class foundational work.
3. Treat phase-transformation semantics as an architecture-first area, not an opportunistic feature area.
4. Grow diffraction validation doctrine before broadening diffraction-facing claims.

## References

### Normative

- `overview.md`
- `canonical_data_model.md`
- `../standards/reference_canon.md`

### Informative

- `orientation_and_texture_foundation.md`
- `ebsd_foundation.md`
- `diffraction_foundation.md`
