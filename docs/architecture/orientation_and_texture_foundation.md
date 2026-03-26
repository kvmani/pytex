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
- construction-time validation of frame, phase, and symmetry consistency across the orientation and texture domain models
- cached proper point-group operator generation to keep repeated symmetry construction cheap

## Deliberate Current Limits

- exact polyhedral fundamental-region boundaries for every crystal class are not yet implemented
- harmonic ODF expansion and full inversion from experimental PF data are not yet implemented
- exact orientation-space polyhedral regions for all crystal classes are not yet implemented

## Why This Still Moves The Project Forward

The current implementation fixes the shared semantics that later harmonic, inversion, EBSD, and plotting work will depend on. It is intentionally more concerned with explicit meaning and reproducibility than with fully optimized texture algorithms at this stage.
