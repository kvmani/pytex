# Texture Analysis In The Workbench: Measured Pole Figures And Kearns Sections

The **Texture** workspace of the workbench has three sub-panels, and this page documents the two
that read measurements:

- **Measured texture** — pole figures, the crystal and the sample symmetry are stated once, and
  one analysis returns the measured, symmetrized, recalculated and difference pole figures, the
  ODF sections and the volume fractions of ideal orientations, each in its own tab.
- **Kearns parameter**, route *Kearns f from three measured scans* — three symmetric
  $\theta$–$2\theta$ scans, one per principal section, give $f_a$, $f_r$ and $f_t$ in one
  calculation, with every peak, intensity and quadrature step on screen.

The third, **Texture**, draws model textures built from named components; see
[the workbench page](workbench_application.md).

Both analyses follow one rule: a number is shown with the intermediate results it came from, so it
can be checked without trusting the program. Each also runs, with no file open, on demonstration
data generated from a texture whose answer is known, so the method can be read against the truth
before it is applied to a specimen.

## The Inputs, Stated Once

| Input | What it decides |
| --- | --- |
| Pole-figure files (XRDML) | The measurement. One file per reflection. |
| Plane of each file $\{hkl\}$ | Which reflection each file is, in the order the files were opened: XRDML records the diffraction angle, not the plane. |
| Phase | The crystal symmetry: it folds the ODF, fixes the Euler range of the sections and chooses the ideal-orientation catalogue. |
| Sample symmetry | Triclinic (none), monoclinic, orthorhombic or axial; imposed on the figures before inversion. |

The files are held by the workspace rather than by a panel (`static/js/core/texturefiles.js`), so a
set opened in *Measured texture* is the set the Kearns pole-figure and ODF routes read. Changing
tab never re-enters an input. The service caches the inversion on everything that changes it —
files, phase, planes, sample symmetry and inversion settings — so choosing other sections or
another tolerance returns without solving again.

The operation behind the panel is `texture.analysis` in `pytex/app/services/texture_analysis.py`.

## Sample Symmetry, Including Axial

A sample symmetry is an assumption about the process: rolling is conventionally orthorhombic about
RD, TD and ND; drawing, extrusion and a tube read along its axis are **axial** (fibre-symmetric)
about their axis (Bunge 1982). Imposing it averages the measured figure over the
group's operations, which improves the statistics when the specimen has the symmetry and fabricates
it when it does not.

{func}`~pytex.texture.sample_symmetry.impose_sample_symmetry` applies it on the figure's own
measured directions, so nothing outside the measured cap is extrapolated.

- **Finite groups.** Every direction's orbit under the proper operators is formed, and each image
  is read from the nearest measured direction. The two-fold axes of the monoclinic and orthorhombic
  groups map a pole at polar angle $\theta$ to $180^{\circ}-\theta$, which the antipodal fold
  returns to $\theta$, so every image stays on the same tilt ring.
- **Axial symmetry** is the continuous Curie group $\infty/mm$, whose proper part is every rotation
  about the fibre axis with the two-folds perpendicular to it. No finite operator list represents
  it, so it is not approximated by one. An axially symmetric figure depends on the polar angle
  alone, and the group average is the azimuthal mean on each ring,

$$
\bar{P}(\theta) = \frac{1}{2\pi}\oint P(\theta, \psi)\,d\psi
\approx \frac{\sum_{j} w_{j} P(\theta, \psi_{j})}{\sum_{j} w_{j}},
\qquad w_{j} = \tfrac{1}{2}\left(\psi_{j+1} - \psi_{j-1}\right),
$$ (eq-axial-ring-average)

  taken with circular trapezoid weights so that an unevenly sampled ring is integrated rather than
  averaged. {meth}`SymmetrySpec.specimen("axial") <pytex.core.symmetry.SymmetrySpec.specimen>`
  returns the group (aliases `fibre`, `fiber`, `cylindrical`); where an operator array is
  unavoidable it holds the closed dihedral group $D_{72}$, rotations every $5^{\circ}$.

**Reading the assumption.** Stage 2 of the result reports the RMS change the symmetry made to each
measured figure. A change comparable with the counting noise means the specimen has the symmetry;
a change comparable with the texture means it does not. The *Symmetrized* column of the
pole-figure plate shows what was averaged away.

## Measured, Recalculated And Difference Pole Figures

The inversion is ill-posed — pole figures are projections and lose the odd-order part of the ODF
(Bunge 1982) — so an ODF is not accepted on its own residual. It is projected back onto
every measured direction with
{meth}`~pytex.texture.reconstruction.PoleFigureResidualReport.from_odf`, and the *Pole figures* tab
draws, for every figure, four readings on one intensity scale:

| Reading | Definition |
| --- | --- |
| Measured | As recorded, normalised to multiples of a random distribution. |
| Symmetrized | With the sample symmetry imposed: the data actually inverted. |
| Recalculated | $P_{\mathrm{calc}}$, the pole density the ODF predicts at each measured direction. |
| Difference | $P_{\mathrm{calc}} - P_{\mathrm{obs}}$, on its own diverging scale centred on zero. |

The difference has its own scale because it answers a different question: $+0.3$ m.r.d. there is
not a weak pole but the ODF over-predicting by $0.3$ at that direction. Noise spread over the whole
figure is expected; a coherent lobe is a component the ODF missed, or a systematic error in the
measurement.

Its size is summarised by the **RP factor** (Matthies, Wenk and Vinel 1988),

$$
RP = \frac{100}{N}\sum_{P_{\mathrm{obs}} > \varepsilon}
\frac{\left|P_{\mathrm{calc}} - P_{\mathrm{obs}}\right|}{P_{\mathrm{obs}}},
\qquad \varepsilon = 0.5\ \text{m.r.d.},
$$ (eq-rp-factor)

reported twice: against the symmetrized figures that were inverted (the quality of the fit alone)
and against the measured figures (the fit plus the sample-symmetry assumption). **The gap between
the two is what the assumption costs.** Below about $10\%$ the service calls the fit good.

On the demonstration data the two cases read as they should. A rolled fcc sheet that genuinely is
orthorhombic, inverted with orthorhombic symmetry over 800 dictionary orientations, gives
$RP = 6.1\%$ against the inverted figures and $6.7\%$ against the measured ones. A split-basal
zirconium sheet — basal poles $30^{\circ}$ either side of ND towards TD, which is orthorhombic but
not axial — inverted with axial symmetry gives $7.8\%$ and $59.6\%$, and a symmetry change of
$3.1$ m.r.d. against $0.11$ with orthorhombic symmetry: the ODF fits what it was given, and the
measurement says it was given the wrong thing.

## ODF Sections: Constant φ₂, φ₁ Or σ, And The LaboTex Plate

{func}`~pytex.texture.sections.odf_sections` slices either ODF representation at constant
$\varphi_2$ (the default and the convention of most texture papers), constant $\varphi_1$, or
constant $\sigma = \varphi_1 + \varphi_2$ (the view of the bcc $\gamma$ fibre), in m.r.d. for both.
The box each section covers is not fixed at the cubic $[0, 90^{\circ}]^3$; it is derived from the
symmetries by {func}`~pytex.texture.sections.euler_section_ranges`
(Bunge 1982; Randle and Engler 2009):

| Angle | Range | Set by |
| --- | --- | --- |
| $\varphi_2$ | $360^{\circ}/n$ | $n$, the order of the crystal's rotation axis along its third axis: $90^{\circ}$ cubic and tetragonal, $60^{\circ}$ hexagonal, $120^{\circ}$ trigonal. |
| $\Phi$ | $90^{\circ}$ or $180^{\circ}$ | $90^{\circ}$ when a two-fold axis lies perpendicular to the third axis of the crystal or of the specimen. |
| $\varphi_1$ | $360^{\circ}$, $180^{\circ}$ or $90^{\circ}$ | The specimen symmetry: triclinic, monoclinic, orthorhombic. Axial: the density does not depend on $\varphi_1$ at all; $90^{\circ}$ is drawn. |

A section drawn over the cubic box for a hexagonal phase repeats a third of itself, and one drawn
over $[0, 90^{\circ}]$ of $\varphi_1$ for a triclinic specimen hides three quarters of the texture;
both were possible before the range was derived.

Three presets choose the sections:

- **Standard sections** — $\varphi_2 = 0, 45, 65^{\circ}$ for cubic crystals, which between them
  carry cube, Goss, brass, copper and S; $\varphi_2 = 0, 30^{\circ}$ for hexagonal ones.
- **Every section (LaboTex plate)** — every section at a fixed step over the whole range, on one
  shared scale, as LaboTex presents an ODF (Pawlik and Ozga 1999). Both edges of the box are
  drawn, so a component on the boundary is read without wrapping round.
- **Chosen values** — any list; a value outside the range the symmetries require is refused.

The kind and preset can be changed from the *ODF sections* tab itself, without returning to the
rail.

## Volume Fractions Of Ideal Orientations

The *Volume fractions* tab integrates the ODF over a misorientation ball of radius $w$ about each
ideal orientation of the crystal system, and sets the fraction beside what a random texture holds in
the same ball. Because the density is invariant under the $|G|$ crystal-symmetry operators, the
integral over the $|G|$ equivalent balls is

$$
V_{\mathrm{component}} = |G|\,\frac{w - \sin w}{\pi}\,\langle f \rangle_{\mathrm{ball}},
\qquad V_{\mathrm{random}} = |G|\,\frac{w - \sin w}{\pi},
$$ (eq-component-volume-fraction)

where $(w - \sin w)/\pi$ is the Haar measure of a rotation ball and $\langle f \rangle$ is the mean
density inside it, estimated from orientations sampled uniformly within the ball. The ratio is the
*times random* column. The expression holds while neighbouring equivalent balls do not overlap:
$w < 45^{\circ}$ for cubic and $w < 30^{\circ}$ for hexagonal crystals. Components are independent,
so overlapping balls of different components may both claim the same volume and the fractions need
not sum to one. The API is {func}`~pytex.texture.components.odf_component_volume_fractions` and
{func}`~pytex.texture.components.random_component_fraction`.

| Crystal | Catalogue |
| --- | --- |
| Cubic | Cube, Goss, brass, copper, S, rotated cube, rotated Goss |
| Hexagonal | Basal ($c \parallel$ ND); basal $30^{\circ}$ from ND towards TD; towards RD; $c \parallel$ TD; $c \parallel$ RD |

The hexagonal components are named by where the basal pole lies. With Bunge angles the crystal
$[0001]$ appears in specimen axes at $(\sin\varphi_1\sin\Phi, -\cos\varphi_1\sin\Phi, \cos\Phi)$, so
$\Phi$ is the basal tilt from ND and $\varphi_1$ the direction it tilts towards; each position is
pinned by a test.

## Kearns Parameters From Three Measured Scans

The route *Kearns f from three measured scans* (`kearns.from_three_sections`, in
`pytex/app/services/kearns_sections.py`) applies Kearns' diffractogram method
(Kearns 1965) to one scan of each principal section and returns the triad. The theory is
in [the Kearns note](../theory/kearns_parameter_and_basal_pole_texture.md); this is what the panel
does with it.

**Which scan goes in which slot.** $f$ along a direction is measured on the section whose *surface
normal* is that direction, because a symmetric scan diffracts only from planes parallel to the
surface.

| Slot | Tube | Plate | Surface prepared |
| --- | --- | --- | --- |
| Axial | $f_a$ | $f_{\mathrm{RD}}$ | A cross-section: surface normal along the tube axis. |
| Radial | $f_r$ | $f_{\mathrm{ND}}$ | The outer surface, flattened: surface normal radial. |
| Transverse | $f_t$ | $f_{\mathrm{TD}}$ | A longitudinal section tangent to the wall: normal along the hoop. |

**What is done to each scan**, with each step's numbers in its own result stage:

1. Peaks are detected against the counting noise and fitted with pseudo-Voigt profiles, the
   K-$\alpha_2$ partner modelled at its Bragg position.
2. Every reflection the phase predicts in the angular range is matched to the nearest fitted peak
   within the tolerance. Its fate is recorded in the reflection table and drawn on the scan:
   *used*; *not detected*, counted as zero intensity because a missing peak is texture too;
   *overlaps another reflection* or *too weak in a random powder*, excluded.
3. The fitted area (or height) is divided by the same reflection's random-powder intensity —
   calculated from multiplicity, $|F|^2$ and the Lorentz-polarisation factor, or measured on an
   opened random standard. The ratio is the basal-pole density at the reflection's tilt to
   $[0001]$.
4. The densities are interpolated onto tilt bins and integrated by Kearns' Eq. (5); the quadrature
   table lists every node, and its contributions sum to $f$.

**The closure check is genuine here.** The three values come from three surfaces measured
independently, so their sum is not 1 by construction. It departs from 1 because the reflections
sample the tilt range unevenly and each section carries its own background and absorption; Kearns
found sums between $0.94$ and $1.06$. The panel reads the sum against that range, and reports the
values normalised to sum to one beside the raw ones, as Mani Krishna *et al.* (2011) does.

**Against a known answer.** With no scans open, three scans are generated from a pilgered-tube
texture — basal poles $30^{\circ}$ either side of radial towards transverse, scattered $15^{\circ}$
in rotation space — with Poisson counting noise, and its exact Kearns parameters, the mean of
$\cos^2$ over its basal poles, are shown beside the measured ones. The route returns
$f_a = 0.0268$, $f_r = 0.7256$, $f_t = 0.2561$ against exact values of $0.0212$, $0.7232$ and
$0.2556$, summing to $1.0085$.

An earlier version of that demonstration added Gaussian noise to each Euler angle, which piles
basal poles up at $\Phi = 0$, where the coordinates are singular. The $(0002)$ density then came out
eighteen times the true peak and $f_r$ was $0.855$ against $0.716$. It is recorded here because it
is exactly the failure Kearns' interpolation is sensitive to: a sharp density spike at a tilt with
no neighbouring reflection is carried across a whole tilt bin.

## Verification

| Claim | Where it is checked |
| --- | --- |
| Axial and orthorhombic averaging remove exactly the terms the group forbids; unevenly sampled rings are trapezoid-integrated | `tests/unit/test_texture_sections_and_sample_symmetry.py` |
| Euler ranges from the operators; random texture 1 m.r.d. on every section kind; $\varphi_1$ and $\varphi_2$ sections agree where they cross | same |
| Hexagonal component basal-pole positions; Haar random fraction; a random ODF is $1\times$ random in every component | same |
| One analysis fills every view; difference = recalculated − measured; the ODF reproduces its data; views reuse the inversion; axial rings are constant | `tests/unit/test_app_texture_analysis.py` |
| The three-section route recovers the model's exact triad within 0.02; quadrature contributions sum to $f$; uploaded scans reproduce the demonstration | `tests/unit/test_app_kearns_sections.py` |
| The panels, end to end in a browser | `tests/browser/workbench.spec.js`: *measured texture: one set of inputs, every reading in tabs*; *Kearns from three section scans shows the triad and every scan* |

## Related Material

- [The Kearns parameter and basal-pole texture](../theory/kearns_parameter_and_basal_pole_texture.md)
- [Texture: discrete and harmonic ODF reconstruction](texture_odf_inversion.md)
- [Harmonic ODF reconstruction](harmonic_odf_reconstruction.md)
- [XRDML texture import](xrdml_texture_import.md)
- [Pole-figure presentation](pole_figure_presentation.md)
- [The PyTex workbench](workbench_application.md)

## References

### Normative

- H.-J. Bunge, *Texture Analysis in Materials Science*, Butterworths (1982), sections 2.3 and 4.2:
  Euler-space ranges and sample symmetry.
- J. J. Kearns, *Thermal Expansion and Preferred Orientation in Zircaloy*, WAPD-TM-472, Bettis
  Atomic Power Laboratory (1965), Eqs. (1)–(7) and Tables 2–3.

### Informative

- V. Randle and O. Engler, *Introduction to Texture Analysis*, 2nd ed., CRC Press (2009), chapter 5.
- S. Matthies, H.-R. Wenk and G. W. Vinel, *J. Appl. Cryst.* 21 (1988) 285–304,
  doi:10.1107/S0021889888000184: the RP factor of recalculated pole figures.
- K. Pawlik and P. Ozga, LaboTex: The Texture Analysis Software, *Göttinger Arbeiten zur Geologie
  und Paläontologie* SB4 (1999).
- K. V. Mani Krishna *et al.*, *J. Nucl. Mater.* 414 (2011) 492–497,
  doi:10.1016/j.jnucmat.2011.04.065: comparison of the Kearns routes and normalisation of the
  triad.
- U. F. Kocks, C. N. Tomé and H.-R. Wenk, *Texture and Anisotropy*, Cambridge University Press
  (1998), chapter 1.
