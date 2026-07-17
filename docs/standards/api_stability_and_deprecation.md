# API Stability And Deprecation

This standard defines what "public API" means for PyTex, which stability
guarantees each surface carries, and the only sanctioned way to remove or
change public behavior.

## Public Surface Definition

- **Stable public API**: every symbol exported from the top-level `pytex`
  namespace (`pytex.__all__`), plus the documented submodule surfaces
  re-exported there (`pytex.core`, `pytex.texture`, `pytex.ebsd`,
  `pytex.diffraction`, `pytex.properties`, `pytex.plotting`,
  `pytex.adapters`, `pytex.contracts`).
- **Experimental API**: everything under `pytex.experimental`. Experimental
  surfaces may change or disappear in any release without deprecation. A
  symbol moves from experimental to stable only with tests, documentation,
  and an entry in the relevant validation ledger.
- **Private**: any module or name with a leading underscore, and everything
  not importable from the surfaces above. No guarantees.
- **Data contracts**: versioned JSON schemas under `schemas/` follow the
  data-contract policy in
  [Data Contracts And Manifests](data_contracts_and_manifests.md); schema
  changes bump the schema version rather than using runtime deprecation.

## Stability Guarantees (pre-1.0)

PyTex is pre-1.0; the contract below is what "pre-alpha" means here in
practice:

1. Stable public symbols are never removed or behavior-changed silently.
   Removal or incompatible change requires a deprecation period of **at least
   two minor releases** (e.g. deprecated in 0.3.0, removable no earlier than
   0.5.0).
2. During the deprecation period the old surface keeps working and emits a
   `DeprecationWarning` created through `pytex._deprecation` (below), naming
   the replacement when one exists and the first release that may remove it.
3. Additive changes (new functions, new keyword-only parameters with
   defaults, new fields at the end of frozen dataclasses with defaults) are
   allowed in any release.
4. Scientific-behavior fixes (a computed value was wrong) are bug fixes, not
   breaking changes; they land with a regression test and a note in the
   release notes, without deprecation.

## The Deprecation Helper

All deprecations go through `pytex._deprecation` so warnings are uniform,
grep-able, and testable:

```python
from pytex._deprecation import deprecated, warn_deprecated

@deprecated(since="0.2.0", removal="0.4.0", replacement="pytex.new_name")
def old_name(...): ...

def old_parameter_path(...):
    warn_deprecated(
        "old_parameter_path(legacy=True)",
        since="0.2.0",
        removal="0.4.0",
        replacement="mode='modern'",
    )
```

Rules:

- `since` and `removal` are release strings; `removal` is the first release
  allowed to delete the surface, not a promise that it will.
- The warning category is always `DeprecationWarning` with `stacklevel`
  pointing at the caller, so user code (not PyTex internals) is implicated.
- Every deprecation lands in the same commit as its replacement and is listed
  in the commit message.
- A deprecated symbol keeps its tests until removal; the tests assert both
  the warning and the preserved behavior.

## Removal Checklist

1. Confirm the deprecation has been released for two minor versions.
2. Delete the surface, its tests, and its `__all__` entries in one commit.
3. Note the removal in the commit message ("removes X, deprecated since Y").

## Release And Changelog Policy

- `CHANGELOG.md` (repository root, Keep-a-Changelog format) is the durable
  user-facing record of behavior changes; the `Unreleased` section accumulates
  entries as features land and is cut into a version section at release time.
- Scientific behavior changes must be stated explicitly and categorized
  honestly: correctness fixes under **Fixed** (even when embarrassing),
  convention or semantics changes under **Changed** — downstream analyses
  depend on them, and silence is a correctness bug of the documentation.
- A release consists of: version bump in `pyproject.toml`, changelog section
  cut, a git tag `vX.Y.Z`, and green CI on the tagged commit. Pre-1.0, minor
  versions may break with deprecation warnings per the stability guarantees
  above; patch versions must not.

## References

### Normative

- [Engineering Governance](engineering_governance.md)
- [Data Contracts And Manifests](data_contracts_and_manifests.md)

### Informative

- `docs/roadmap/world_class_feature_roadmap.md` (item I3)
