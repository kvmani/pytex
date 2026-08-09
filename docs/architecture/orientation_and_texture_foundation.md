# Orientation And Texture Foundation

This document records the current implementation posture for Phase 2.

## Implemented Core Behavior

- proper rotational point-group generation for common crystallographic groups
- alias handling for Laue-style names such as `m-3m` and `6/mmm`
- symmetry-aware disorientation by minimum-angle reduction over crystal symmetry operators
- deterministic symmetry-reduced orientation canonicalization
- Bunge Euler import and export on top of quaternion-backed rotations
- convention-aware Euler conversion support for Bunge and Matthies/ABG parity workflows
- the two **equal-volume** charts of SO(3), homochoric and cubochoric, with vectorized
  conversions to and from every other representation, so uniform sampling and dictionary
  grids need no chart-specific correction
- a single-call all-representations report (`orientation_representations`) that converts one
  orientation into all ten numerical forms with `describe()` and a JSON contract
- recovery of the ideal `(hkl)[uvw]` texture-component name from an orientation, reported
  with its residual angles and with the four-index form for hexagonal and trigonal phases
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
- recorded pole-figure sampling semantics: whether intensities are per-pole weights of a cloud or densities evaluated at given directions, which decides the correct resampling estimator
- spherical resampling of a pole figure onto an arbitrary `S2Grid`, giving two figures a common support
- multiples-of-random normalization from solid-angle integration weights, including a raster weighting for measured tilt/rotation grids
- pole-figure arithmetic on a shared support, with a signed `PoleFigureDifference` for subtraction
- residual pole figures from an ODF reconstruction, as a plottable diagnosis rather than a scalar norm

### Why pole-figure arithmetic needed three steps, not one

The operators were not merely unwritten; they were undefined. A `PoleFigure`
carries *scattered* specimen directions, so two figures generally share no
sampling direction and there is nothing to combine pointwise. Measured
intensities additionally arrive in detector counts or scaled by their own
maximum, so even on a shared support their magnitudes are not comparable.
Resampling supplies the shared support, m.r.d. supplies the shared scale, and
only then does arithmetic mean anything. The sequence is a dependency chain,
not a convenience ordering.

Subtraction is the one operator that cannot return a `PoleFigure`. A pole
density is non-negative and the type enforces it; a difference is signed, and
its sign — which regions a model over-predicts and which it under-predicts — is
the entire content of the result. Hence `PoleFigureDifference`.

## Deliberate Current Limits

- exact polyhedral fundamental-region boundaries for every crystal class are not yet implemented
- broad experimentally calibrated PF inversion doctrine beyond the current kernel-regularized harmonic model is still ahead
- exact orientation-space polyhedral regions for all crystal classes are not yet implemented
- m.r.d. normalization over a partial pole figure averages over the *measured*
  cap; it equals the true spherical mean only if the unmeasured region has the
  same mean. Defocusing limits the reachable tilt, so this assumption is real
  and is stated at the call site rather than hidden
- the resampling kernel is a fixed von Mises-Fisher shape; a kernel library on
  S2 matching the existing SO(3) one is not yet implemented
- ghost correction and zero-range methods are still absent, so the odd part of
  a reconstructed ODF remains unconstrained

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
