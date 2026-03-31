# Benchmarks And Manifests

PyTex uses benchmark and validation manifests as part of its scientific quality system, not as
optional bookkeeping.

## Benchmark Families

- `benchmarks/ebsd/`
- `benchmarks/diffraction/`
- `benchmarks/texture/`
- `benchmarks/plotting/`
- `benchmarks/structure_import/`
- `benchmarks/transformation/`
- `benchmarks/validation/`

## High-Value Current Benchmark Surfaces

- graph-backed and multiphase EBSD map workflows
- pinned `pymatgen` powder XRD starter baselines
- pinned `diffsims` SAED shell-geometry starter baselines
- structure-import fixture audit and validation manifests
- texture import and inversion fixtures for XRDML and LaboTex-compatible workflows
- transformation variant prediction and experimental parent-candidate scoring

## CLI Support

Use the lightweight CLI inventory when you want a quick audit of the manifest surface:

```bash
pytex benchmarks inventory
pytex validate manifests
```

## Related Material

- {doc}`../validation/index`
- {doc}`../architecture/phase_transformation_foundation`
- {doc}`../standards/data_contracts_and_manifests`
