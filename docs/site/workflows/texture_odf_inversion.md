# Texture: ODF Evaluation And Pole-Figure Inversion

PyTex now supports a richer texture workflow surface than the earlier ODF-evaluation-only foundation.

Current scope:

- `PoleFigure.from_orientations(...)`
- `InversePoleFigure.from_orientations(...)`
- `ODF.from_orientations(...)`
- `ODF.evaluate(...)`
- `ODF.reconstruct_pole_figure(...)`
- `ODF.reconstruct_pole_figures(...)`
- `ODF.invert_pole_figures(...)`
- `ODFInversionReport`
- `plot_pole_figure(..., kind="scatter" | "histogram" | "contour")`
- `plot_odf(..., kind="scatter" | "contour" | "sections")`

## Discrete Inversion Model

The current inversion path is intentionally explicit and dictionary-based.

Given an orientation dictionary `g_j`, measurement directions `s_m`, and a pole family `H`, PyTex builds the forward model

$$
A_{mj} = \frac{1}{|H|}\sum_{h \in H} K\bigl(\angle(s_m, g_j h)\bigr),
$$

then solves a regularized nonnegative least-squares problem for the dictionary weights.

This is not yet a harmonic inversion framework or a claim of full experimental PF inversion breadth. It is a scientifically explicit discrete inversion surface that fits the current PyTex data model.

## Estimation And Inspection Surface

PyTex currently exposes two closely related texture-estimation surfaces:

- a discrete ODF represented by an explicit orientation support together with non-negative weights
- a dictionary-based inversion path that estimates those weights from one or more pole figures

Once that discrete model exists, the same support can be used to:

- evaluate the ODF at query orientations
- reconstruct pole figures from the estimated support
- inspect the support in Euler space through scatter or contour views
- inspect the estimated texture through kernel-smoothed Bunge-section plots

That means plotting is part of the texture interpretation workflow, not an unrelated presentation layer.

## Example

```python
import numpy as np

from pytex import (
    CrystalPlane,
    FrameDomain,
    Handedness,
    KernelSpec,
    Lattice,
    MillerIndex,
    ODF,
    OrientationSet,
    Phase,
    ReferenceFrame,
    SymmetrySpec,
)

crystal = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
specimen = ReferenceFrame("specimen", FrameDomain.SPECIMEN, ("x", "y", "z"), Handedness.RIGHT)
symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
phase = Phase("fcc_demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)

dictionary = OrientationSet.from_euler_angles(
    np.array(
        [
            [0.0, 0.0, 0.0],
            [90.0, 0.0, 0.0],
            [35.0, 25.0, 10.0],
        ]
    ),
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
    phase=phase,
)

true_odf = ODF.from_orientations(
    dictionary,
    weights=[4.0, 2.0, 1.0],
    kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=10.0),
)

poles = (
    CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase),
    CrystalPlane(miller=MillerIndex([1, 1, 1], phase=phase), phase=phase),
)
pole_figures = true_odf.reconstruct_pole_figures(
    poles,
    include_symmetry_family=False,
    antipodal=False,
)

report = ODF.invert_pole_figures(
    pole_figures,
    orientation_dictionary=dictionary,
    kernel=true_odf.kernel,
    regularization=1e-8,
    include_symmetry_family=False,
)

print(report.converged)
print(report.odf.normalized_weights)
```

## Plotting The Estimated Texture

```python
from pytex import plot_odf, plot_pole_figure

pole_figure_plot = plot_pole_figure(
    pole_figures[0],
    kind="contour",
    bins=81,
    sigma_bins=1.5,
    levels=12,
)

odf_sections = plot_odf(
    report.odf,
    kind="sections",
    section_phi2_deg=(0.0, 15.0, 30.0, 45.0, 60.0, 75.0),
    section_phi1_steps=121,
    section_big_phi_steps=61,
    levels=10,
)
```

## Interpretation Notes

- The dictionary is explicit. PyTex does not hide the orientation support over which the inversion is performed.
- The inversion report is explicit. Residual norm, iteration count, and convergence state are part of the returned surface.
- The current method is appropriate for the present discrete kernel foundation. It should not yet be read as a replacement for broader harmonic or experimentally calibrated inversion frameworks.
- Pole-figure contours are built from a smoothed projected density grid over the discrete pole data.
- ODF section plots are kernel-smoothed Bunge-section inspection views over the current discrete support.

## Current Limits

- the inversion path is dictionary-based, not harmonic
- sample-symmetry handling is currently conservative and explicit rather than aggressively implicit
- ODF section plots are inspection surfaces for the current discrete model, not a harmonic ODF implementation
- richer experimental PF inversion doctrine and benchmarking are still ahead

## Related Material

- {doc}`../concepts/orientation_texture`
- {doc}`../tutorials/notebooks`
- [../../tex/algorithms/discrete_odf_and_pole_figures.tex](../../tex/algorithms/discrete_odf_and_pole_figures.tex)
