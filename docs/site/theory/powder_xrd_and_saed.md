# Powder XRD And SAED Foundations

## Powder XRD

PyTex models a powder reflection from an integer Miller triplet $(hkl)$ through the reciprocal-lattice vector

$$
\mathbf{g}_{hkl} = h \mathbf{a}^{*} + k \mathbf{b}^{*} + l \mathbf{c}^{*}
$$

with

$$
d_{hkl} = \frac{1}{\lVert \mathbf{g}_{hkl} \rVert}
$$

For a wavelength $\lambda$, the current XRD surface applies Bragg's law

$$
2 d_{hkl} \sin \theta = \lambda
$$

and reports the observable angle $2\theta$.

The powder intensity model is a deliberately bounded kinematic approximation:

$$
I_{hkl} \propto m_{hkl} \, |F_{hkl}|^2 \, L_{p}(2\theta)
$$

where $m_{hkl}$ is the multiplicity inferred from the phase point-group symmetry, $F_{hkl}$
comes from the selected tabulated-X-ray, constant-atomic-number, or unit-amplitude model, and
$L_p$ is a Lorentz-polarization term.

Enumeration returns one deterministic representative for each symmetry orbit $\{hkl\}$, while
$m_{hkl}$ records the size of that orbit. Emitting every equivalent index as a separate reflection
would apply multiplicity twice when constructing a powder spectrum and is therefore excluded by the
family-uniqueness invariant.

This is suitable for phase-identification workflows, controlled comparison, and teaching. It is
not a calibrated instrument model or a Rietveld refinement.

## Spectrum Construction

The broadened powder spectrum is constructed by depositing each reflection onto a sampled
$2\theta$ grid. Gaussian and pseudo-Voigt profiles are available with constant width or a Caglioti
angle-dependent width. For the Gaussian case centered at $2\theta_{hkl}$ with width $\sigma$,

$$
I(2\theta) = \sum_{hkl} I_{hkl} \exp\left[-\frac{(2\theta - 2\theta_{hkl})^2}{2\sigma^2}\right]
$$

## Measured Profiles And Whole-Profile Comparison

`MeasuredPowderPattern` is the experiment-facing object. Its angular support must be finite and
strictly increasing; intensity must be finite and non-negative; optional standard uncertainties
must be positive. It records whether intensity means counts, counts per second, or arbitrary units,
and it carries radiation, metadata, provenance, and an explicit `synthetic` label. The whitespace or
CSV reader accepts comment metadata without silently normalizing, subtracting background, or
resampling the observation.

Comparison interpolates the simulated profile onto measured angles only inside their shared range.
It then solves

$$
I_{\mathrm{calc},i} = a I_{\mathrm{sim},i} + b
$$

for non-negative scale $a$ and optional constant background $b$ by weighted least squares. If the
measurement supplies standard uncertainties $\sigma_i$, $w_i=1/\sigma_i^2$; otherwise $w_i=1$.
PyTex reports the IUCr pdCIF whole-profile agreement factors

$$
R_p = \frac{\sum_i |I_{\mathrm{obs},i}-I_{\mathrm{calc},i}|}
           {\sum_i I_{\mathrm{obs},i}},
\qquad
R_{wp} = \left[
\frac{\sum_i w_i(I_{\mathrm{obs},i}-I_{\mathrm{calc},i})^2}
     {\sum_i w_i I_{\mathrm{obs},i}^2}
\right]^{1/2}.
$$

These numbers describe pointwise profile agreement after only scale/background alignment. They are
not an expected R factor, reduced $\chi^2$, structure refinement, phase quantification, or permission
to compare patterns whose peaks are substantially displaced. The residual array remains attached to
the result because the difference profile is more diagnostic than a scalar alone.

## SAED

PyTex treats the SAED zone axis as a direct-space direction while the diffracted reflections live in reciprocal space. This distinction is explicit in the public data model.

Given a zone axis unit vector $\hat{\mathbf{z}}$, a reflection is accepted into the current in-zone SAED construction when

$$
|\mathbf{g}_{hkl} \cdot \hat{\mathbf{z}}| \le \varepsilon
$$

for a configurable tolerance $\varepsilon$ in inverse angstroms.

An orthonormal detector basis $(\hat{\mathbf{u}}, \hat{\mathbf{v}}, \hat{\mathbf{z}})$ is then constructed, and the reciprocal vector is projected into detector coordinates:

$$
u = C \, (\mathbf{g}_{hkl} \cdot \hat{\mathbf{u}}), \qquad
v = C \, (\mathbf{g}_{hkl} \cdot \hat{\mathbf{v}})
$$

where $C$ is the camera-constant-style scale factor in millimeter-angstrom units.

The current spot intensity is again a proxy quantity intended for ranking and visualization, not a dynamical diffraction model.

## Current Limits

- Powder comparison does not refine peak shift, lattice parameters, profile width, structure,
  specimen displacement, or an instrument response.
- Powder XRD remains kinematic; absorption, fluorescence, extinction and detector response are not
  inferred from a measured profile.
- SAED spot intensity remains a proxy in this module; dynamical electron diffraction is provided by
  the separate CBED/Bloch-wave surfaces rather than by the kinematic SAED generator.

## References

- IUCr, [Powder CIF dictionary: profile R-factor definitions](https://www.iucr.org/resources/cif/dictionaries/browse/cif_pd).
- McCusker *et al.*, “Rietveld refinement guidelines,” *J. Appl. Cryst.* 32 (1999),
  [doi:10.1107/S0021889898009856](https://doi.org/10.1107/S0021889898009856).
- Young (ed.), *The Rietveld Method*, IUCr Monographs on Crystallography 5 (1993).
