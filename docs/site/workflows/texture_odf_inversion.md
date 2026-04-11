# Texture: Discrete And Harmonic ODF Reconstruction

PyTex now supports two explicit PF-to-ODF reconstruction surfaces:

- `ODF.invert_pole_figures(...)` for dictionary-based non-negative inversion over an
  explicit orientation support
- `HarmonicODF.invert_pole_figures(...)` for band-limited harmonic reconstruction with
  explicit crystal and specimen symmetry handling
- `ODF.evaluate_pole_density(...)` and `HarmonicODF.evaluate_pole_density(...)` for
  comparing reconstructed ODFs against explicit sample-direction measurement grids

Both methods use the same `PoleFigure` object, the same phase and symmetry semantics, and
the same plotting and validation surfaces. The distinction is algorithmic, not semantic.

## Choosing Between The Two

Use the discrete dictionary path when:

- the orientation support should remain fully explicit and inspectable
- a teaching workflow needs direct support-point reasoning
- the inversion is naturally posed on a known orientation dictionary

Use the harmonic path when:

- a band-limited ODF field is the desired output
- crystal and specimen symmetry should be enforced directly in the reconstruction basis
- Bunge-section inspection and PF refitting are the primary review surfaces

## Discrete Dictionary Inversion

The discrete path builds a forward matrix over an orientation dictionary `g_j`:

```text
A_mj = (1 / |H|) sum_{u in H} K(angle(y_m, g_j u))
```

and solves for non-negative weights over that explicit support. This remains the most
inspectable reconstruction route in the library.

```python
report = ODF.invert_pole_figures(
    pole_figures,
    orientation_dictionary=dictionary,
    kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=10.0),
    regularization=1e-8,
    include_symmetry_family=True,
)
```

## Harmonic Reconstruction

The harmonic path expands the ODF in a truncated symmetry-projected harmonic basis and
solves a regularized least-squares system for the retained coefficients.

```python
harmonic_report = HarmonicODF.invert_pole_figures(
    pole_figures,
    degree_bandlimit=6,
    regularization=1e-6,
    include_symmetry_family=True,
    specimen_symmetry=sample_symmetry,
    pole_kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=7.5),
    phi1_step_deg=30.0,
    big_phi_step_deg=20.0,
    phi2_step_deg=30.0,
)
```

See {doc}`harmonic_odf_reconstruction` for the full mathematical explanation.

## Common Inspection Surface

Both inversion reports are intended to be review objects, not opaque success tokens.

Discrete report:

- residual norm
- relative residual
- mean and maximum absolute error
- predicted intensities
- observation-to-dictionary coverage ratio

Harmonic report:

- residual norm
- relative residual
- mean and maximum absolute error
- matrix rank
- condition number
- retained basis size and raw basis size
- quadrature size
- mean, minimum, and maximum reconstructed density

## Plotting The Result

Discrete ODF:

```python
discrete_sections = plot_odf(
    report.odf,
    kind="sections",
    section_phi2_deg=(0.0, 15.0, 30.0, 45.0, 60.0, 75.0),
)
```

Harmonic ODF:

```python
harmonic_sections = plot_odf(
    harmonic_report.odf,
    kind="sections",
    section_phi2_deg=(0.0, 15.0, 30.0, 45.0, 60.0, 75.0),
    section_phi1_steps=121,
    section_big_phi_steps=61,
)
```

For harmonic reconstructions, Bunge sections are the preferred first review surface.

## Example Workflow

```python
from pytex import (
    CrystalPlane,
    HarmonicODF,
    KernelSpec,
    MillerIndex,
    ODF,
    OrientationSet,
    plot_odf,
)

dictionary = OrientationSet.from_bunge_grid(
    crystal_frame=crystal,
    specimen_frame=specimen,
    symmetry=symmetry,
    phase=phase,
    phi1_step_deg=15.0,
    big_phi_step_deg=15.0,
    phi2_step_deg=15.0,
)

discrete_report = ODF.invert_pole_figures(
    pole_figures,
    orientation_dictionary=dictionary,
    kernel=KernelSpec(name="von_mises_fisher", halfwidth_deg=10.0),
)

harmonic_report = HarmonicODF.invert_pole_figures(
    pole_figures,
    degree_bandlimit=6,
    specimen_symmetry=sample_symmetry,
    pole_kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=7.5),
)

print(discrete_report.relative_residual_norm)
print(harmonic_report.relative_residual_norm)

figure = plot_odf(harmonic_report.odf, kind="sections")
```

## Scientific Review Guidance

- Do not compare the discrete and harmonic routes only by residual norm. Compare also the
  conditioning, basis size, PF refit behavior, and section smoothness.
- For antipodal diffraction pole figures, remember that the harmonic path defaults to
  even degrees only.
- Treat kernel halfwidth, quadrature density, and degree bandlimit as declared scientific
  choices in reports, notebooks, and publications.

## Related Material

- {doc}`harmonic_odf_reconstruction`
- {doc}`xrdml_texture_import`
- {doc}`../validation/automated_test_cases`
- [../../tex/algorithms/discrete_odf_and_pole_figures.tex](../../tex/algorithms/discrete_odf_and_pole_figures.tex)
- [../../tex/algorithms/harmonic_odf_reconstruction.tex](../../tex/algorithms/harmonic_odf_reconstruction.tex)

## References

### Normative

- {doc}`../concepts/orientation_texture`
- {doc}`../validation/mtex_parity_matrix`

### Informative

- `references/reference_index.md`
- `references/formulation_summary.md`
