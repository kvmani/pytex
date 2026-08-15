# Foundation Feature Priorities

This note records the mathematical contract for the first implementation slice of the next PyTex
foundation cycle. The aim is to make reconstruction, symmetry, diffraction physics, EBSD texture
workflows, and parent reconstruction inspectable and serializable before adding broader algorithmic
breadth.

## Pole-Figure Correction And Residual Auditing

For a pole-figure observation $I_i$, scale $s$, background $b$, and optional defocusing factor
$d_i$, PyTex applies

$$
I_i^{\mathrm{corr}} = s\frac{I_i-b}{d_i}.
$$

Background is subtracted before division: background is not specimen diffraction signal and must
not be amplified by the defocusing correction. Negative values are either clipped or rejected
according to the correction policy.

For an untextured reference measured on polar-angle rings, let $I_{ij}^{\mathrm{rand}}$ be the
intensity at tilt ring $i$ and azimuth $j$. After subtracting the standard's declared background,
PyTex reduces each ring by a mean or median, $\bar I_i^{\mathrm{rand}}$, and normalizes to the
lowest measured tilt $i=0$:

$$
d_i = \frac{\bar I_i^{\mathrm{rand}}}{\bar I_0^{\mathrm{rand}}}.
$$

The result retains ring counts and the within-ring relative standard deviation. A true random
standard should be azimuthally uniform; large scatter is therefore evidence against the calibration
assumption rather than texture in the specimen. The curve is tied to the same reflection and
instrument configuration and is interpolated only inside its measured tilt interval. PyTex refuses
extrapolation because no data constrain the loss beyond that interval.

Given a fitted ODF-derived pole density $I_i^{\mathrm{fit}}$, the residual report records

$$
r_i = I_i^{\mathrm{fit}} - I_i^{\mathrm{obs}},
\qquad
\rho = \frac{\|r\|_2}{\max(\|I^{\mathrm{obs}}\|_2,\epsilon)}
$$

### References

- MTEX, [ODF reconstruction from X-ray diffraction data of an Al alloy rolled sheet](https://mtex-toolbox.github.io/ExAlODFReconstruction.html), “Background and Defocusing Correction.”
- MTEX, [`PoleFigure.correct`](https://mtex-toolbox.github.io/PoleFigure.correct.html).
- Welzel and Leoni, “Use of polycapillary X-ray lenses in the X-ray diffraction measurement of
  texture,” *J. Appl. Cryst.* 35 (2002),
  [doi:10.1107/S0021889802000481](https://doi.org/10.1107/S0021889802000481).

## Euler Convention Transforms

For Bunge angles $(\phi_1,\Phi,\phi_2)$, the implemented Roe-style transform is

$$
(\Psi,\Theta,\Phi)_{\mathrm{Roe}}
=
(\phi_1 - 90^\circ,\Phi,\phi_2 + 90^\circ)
$$

and the implemented Kocks-style transform is

$$
(\Psi,\Theta,\Phi)_{\mathrm{Kocks}}
=
(\phi_1 - 90^\circ,\Phi,90^\circ-\phi_2)
$$

These formulas are exposed through a named public transform so validation pages can distinguish
Roe/ABG-style behavior from Kocks notation.

## Structure Factors And Powder Intensity

For reciprocal vector $\mathbf{g}_{hkl}$ and fractional site coordinates $\mathbf{x}_j$, the
implemented structure-factor surface evaluates

$$
F_{hkl}
=
\sum_j o_j f_j(|\mathbf{g}_{hkl}|)
\exp\!\left(-\frac{B_j|\mathbf{g}_{hkl}|^2}{16\pi^2}\right)
\exp\!\left(2\pi i\,\mathbf{h}\cdot\mathbf{x}_j\right)
$$

The initial scattering-factor tables are deliberately limited to unit, atomic-number, and smooth
proxy models. The powder intensity model records

$$
I_{hkl} = m_{hkl}|F_{hkl}|^2 L_p(2\theta),
\qquad
L_p(2\theta)=
\frac{1+\cos^2(2\theta)}{\sin^2\theta\cos\theta}
$$

## EBSD Texture Weights

For raw orientation weights $w_i$ and validity mask $m_i$, the workflow normalizes

$$
w_i^{\mathrm{norm}} =
\frac{m_i w_i}{\sum_j m_j w_j}
$$

These normalized weights are then passed to ODF, pole-figure, and inverse-pole-figure construction.

## Parent Candidate Scoring

For candidate parent orientation $p_m$, observed child orientations, and a selected orientation
relationship or variant family, PyTex computes child residual angles $\Delta\theta_{mi}$ and
reports

$$
s_m = R(\Delta\theta_{m1},\Delta\theta_{m2},…,\Delta\theta_{mn}),
\qquad
R \in \{\mathrm{mean},\mathrm{median},\max\}
$$

Candidates within the configured ambiguity tolerance of the best score are retained in the
`ParentReconstructionReport`.
