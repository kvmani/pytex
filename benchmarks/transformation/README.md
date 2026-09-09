# Transformation Benchmarks

This directory holds benchmark manifests for implemented transformation foundations and
experiment-facing research staging.

Current scope:

- orientation-relationship and variant-generation identity
- variant-indexed child-orientation prediction
- experimental candidate-parent scoring for reconstruction studies

These artifacts define benchmark identity and validation links without overstating full parent
reconstruction maturity.

Run `python scripts/benchmark_transformation_performance.py --quick` for the smoke sizes,
or omit `--quick` for 5,000 paired orientations. The weighted-OR case uses planted KS
variants with nonuniform weights, including an excluded pair; it asserts recovery of the
known relationship and records timing plus peak traced allocation. Tracing covers the fit
call, not input construction, native BLAS workspaces or whole-process resident memory.
Results remain in ignored local output, not in the benchmark manifest directory.
