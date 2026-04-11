# Orientation And Texture Foundation

This document records the current implementation posture for Phase 2.

## Implemented Core Behavior

- proper rotational point-group generation for common crystallographic groups
- alias handling for Laue-style names such as `m-3m` and `6/mmm`
- symmetry-aware disorientation by minimum-angle reduction over crystal symmetry operators
- deterministic symmetry-reduced orientation canonicalization
- Bunge Euler import and export on top of quaternion-backed rotations
- convention-aware Euler conversion support for Bunge and Matthies/ABG parity workflows
- mapping between crystal vectors and specimen vectors for individual orientations and orientation sets
- pole-figure synthesis from orientation sets
- inverse-pole-figure synthesis from orientation sets
- class-specific IPF sector reduction for supported proper point groups
- explicit orientation projection to a symmetry-reduced representative
- kernel-backed ODF evaluation and simple volume-fraction queries
- discrete pole-figure inversion over an explicit orientation dictionary with a regularized non-negative solver
- band-limited harmonic ODF reconstruction with explicit crystal and specimen symmetry handling
- construction-time validation of frame, phase, and symmetry consistency across the orientation and texture domain models
- cached proper point-group operator generation to keep repeated symmetry construction cheap

## Deliberate Current Limits

- exact polyhedral fundamental-region boundaries for every crystal class are not yet implemented
- broad experimentally calibrated PF inversion doctrine beyond the current kernel-regularized harmonic model is still ahead
- exact orientation-space polyhedral regions for all crystal classes are not yet implemented

## Why This Still Moves The Project Forward

The current implementation now covers both explicit discrete inversion and a first harmonic reconstruction path. The remaining work is therefore no longer semantic groundwork alone; it is broader validation, benchmark depth, and higher-fidelity experimental doctrine.

## References

### Normative

- [Canonical Data Model](canonical_data_model.md)
- [Reference Canon](../standards/reference_canon.md)

### Informative

- <a href="../tex/theory/orientation_space_and_disorientation.tex">Orientation Space And Disorientation</a>
- <a href="../tex/algorithms/discrete_odf_and_pole_figures.tex">Discrete ODF And Pole Figures</a>
- <a href="../tex/algorithms/harmonic_odf_reconstruction.tex">Harmonic ODF Reconstruction</a>
