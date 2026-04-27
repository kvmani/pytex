# Foundation Feature Priorities: Implemented Surfaces

This page records the first implementation slice for the next foundation cycle. The goal is not
to claim final breadth in any one scientific area. The goal is to make the next layers of texture,
symmetry, diffraction, EBSD, and phase-transformation work explicit, testable, and serializable.

## Texture Reconstruction Suite V2

PyTex now separates pole-figure correction, reconstruction configuration, and residual auditing:

- `PoleFigureCorrectionSpec`
- `ODFReconstructionConfig`
- `PoleFigureResidualReport`
- `residual_reports_for_pole_figures(...)`

The correction layer applies a deterministic intensity transform:

$$
I_i^{\mathrm{corr}} = s \left(\frac{I_i}{d_i} - b\right),
$$

where \(s\) is the scale, \(d_i\) is an optional defocussing factor, and \(b\) is a constant
background. Negative corrected intensities are either clipped to zero or rejected, depending on
the configured policy.

The residual report compares a measured pole figure against a fitted ODF:

$$
r_i = I_i^{\mathrm{fit}} - I_i^{\mathrm{obs}},
\qquad
\rho = \frac{\lVert r\rVert_2}{\max(\lVert I^{\mathrm{obs}}\rVert_2,\epsilon)}.
$$

This keeps PF correction and PF refitting reviewable by domain experts without requiring them to
inspect the solver internals.

```python
config = ODFReconstructionConfig(
    algorithm="discrete",
    correction=PoleFigureCorrectionSpec(scale=2.0, background=0.25),
)
report = config.reconstruct(pole_figures, orientation_dictionary=dictionary)
residual = PoleFigureResidualReport.from_odf(report.odf, corrected_pole_figure)
```

## Exact Symmetry And Orientation-Space Geometry Kernel

The new geometry helpers expose named convention and boundary objects:

- `EulerConventionTransform`
- `IPFSectorBoundary`
- `OrientationFundamentalRegion`

For the Randle and Engler Bunge example used in the validation page, PyTex now exposes the
Roe-style and Kocks-style angle transforms as named operations:

$$
(\Psi,\Theta,\Phi)_{\mathrm{Roe}}
=
(\phi_1 - 90^\circ,\Phi,\phi_2 + 90^\circ),
$$

$$
(\Psi,\Theta,\Phi)_{\mathrm{Kocks}}
=
(\phi_1 - 90^\circ,\Phi,90^\circ-\phi_2).
$$

`IPFSectorBoundary.from_symmetry(...)` exposes the currently implemented fundamental-sector
vertices and human-readable boundary equations. For cubic \(432\), the antipodal sector is
represented by inequalities such as:

$$
z \ge 0,\qquad y \ge 0,\qquad x \ge y,\qquad z \ge x.
$$

These are intentionally surfaced as reviewable scientific statements rather than only as hidden
branch logic.

## Structure-To-Diffraction Physics Layer

The diffraction layer now has public physics primitives:

- `ScatteringFactorTable`
- `StructureFactor`
- `ReflectionCondition`
- `DiffractionIntensityModel`

The structure factor is evaluated as:

$$
F_{hkl}
=
\sum_j o_j f_j(|\mathbf{g}_{hkl}|)
\exp\!\left(-\frac{B_j|\mathbf{g}_{hkl}|^2}{16\pi^2}\right)
\exp\!\left(2\pi i\,\mathbf{h}\cdot\mathbf{x}_j\right),
$$

where \(o_j\), \(f_j\), \(B_j\), and \(\mathbf{x}_j\) are occupancy, scattering factor,
isotropic displacement parameter, and fractional coordinates for site \(j\). The first stable
scattering table is deliberately conservative: it supports unit scattering, atomic-number
scattering, and a smooth X-ray proxy. It is not a tabulated Cromer-Mann implementation.

Reflection conditions support common centered lattices. For example, face-centered reflections
are allowed only when \(h\), \(k\), and \(l\) have the same parity.

$$
I_{hkl}
=
m_{hkl}|F_{hkl}|^2 L_p(2\theta),
$$

with the current powder correction

$$
L_p(2\theta)
=
\frac{1+\cos^2(2\theta)}{\sin^2\theta\cos\theta}.
$$

`PowderReflection` now stores optional structure-factor real/imaginary components, the
Lorentz-polarization factor, and the intensity-model name.

## Workflow-Grade EBSD To Texture Pipeline

`EBSDTextureWorkflow` turns the existing `CrystalMap` methods into one reproducible entry point:

1. resolve an optional phase selection
2. normalize orientation quality weights
3. optionally segment grains
4. compute ODF, PF, and IPF outputs
5. emit an experiment manifest and a typed `EBSDTextureWorkflowResult`

Weights are normalized after validity masking:

$$
w_i^{\mathrm{norm}} =
\frac{m_i w_i}{\sum_j m_j w_j}.
$$

This makes EBSD-derived texture workflows ready for later confidence, calibration, and detector
quality handling without changing the public workflow shape.

```python
workflow = EBSDTextureWorkflow(
    poles=(plane_111,),
    weights=OrientationQualityWeights([0.2, 0.3, 0.5]),
)
result = workflow.run(crystal_map)
```

## Stable Parent-Reconstruction Track

The experimental scoring kernel now has a stable workflow wrapper:

- `ParentReconstructionConfig`
- `ParentReconstructionReport`
- `OrientationRelationshipCatalog`
- `VariantSelectionReport`
- `reconstruct_parent_orientation(...)`

For each candidate parent orientation \(p_m\), PyTex predicts the child orientations through the
selected orientation relationship and variants, computes angular residuals against observed child
orientations, and reduces them by a configured statistic:

$$
s_m = R\left(\Delta\theta_{m1},\Delta\theta_{m2},\ldots,\Delta\theta_{mn}\right),
\qquad
R \in \{\mathrm{mean},\mathrm{median},\max\}.
$$

The report exposes the best candidate, its score in degrees, and all candidates within the
configured ambiguity tolerance. This is a stable staged workflow, not yet a full parent
reconstruction ecosystem.

## JSON Contract Coverage

The JSON contract layer now includes the previously missing texture and EBSD outputs:

- `KernelSpec`
- `PoleFigure`
- `InversePoleFigure`
- `ODF`
- `ODFInversionReport`
- `HarmonicODF`
- `CrystalMap`
- `TextureReport`

This closes the immediate interchange gap for texture and EBSD results. Figures are not serialized
inside `TextureReport`; the contract preserves the scientific data needed to reconstruct or review
the result.

## Automated Coverage

The implementation is covered by `tests/unit/test_foundation_feature_priorities.py`.

The tests assert:

- pole-figure correction, reconstruction configuration, residual reports, and JSON round trips
- Roe and Kocks convention transforms from the reference Bunge example
- IPF sector boundary membership and exposed equations
- structure-factor amplitude, face-centered reflection conditions, and intensity calculation
- workflow-grade EBSD-to-texture execution and JSON round trips
- stable parent candidate scoring and ambiguity reporting

## Current Limits

- `ScatteringFactorTable(model="xray_gaussian_proxy")` is a smooth proxy, not a tabulated
  physical scattering-factor database.
- `IPFSectorBoundary` exposes current exact vector-sector inequalities; full class-by-class
  orientation-space polyhedral catalogs remain a later expansion.
- `EBSDTextureWorkflow` carries quality weights into texture outputs but does not yet propagate
  uncertainty into every downstream statistic.
- Parent reconstruction is still candidate-scoring based. Grain-graph reconstruction and
  literature-wide OR catalogs remain future work.

## Related Material

- {doc}`texture_odf_inversion`
- {doc}`harmonic_odf_reconstruction`
- {doc}`diffraction_spots`
- {doc}`xrd_generation`
- {doc}`ebsd_to_texture_outputs`
- {doc}`phase_transformation_manifests_and_scoring`
- [../../tex/algorithms/foundation_feature_priorities.tex](../../tex/algorithms/foundation_feature_priorities.tex)
