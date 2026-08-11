# Convergent-Beam Electron Diffraction

Selected-area diffraction illuminates the specimen with a parallel beam, so each reflection
is a spot and the recorded pattern carries the projected symmetry of the zone and nothing
more. Converging the beam to a probe of semi-angle $\alpha$ turns every spot into a
*disc* in which each point corresponds to a different incident direction — that is,
to a different deviation from the Bragg condition. One CBED exposure therefore records the
rocking curve that a parallel beam could obtain only by tilting through a series of
exposures.

This note derives the disc geometry, the two-beam rocking curve, the thickness
determination it supports, and the HOLZ ring radii. The implementation is
`pytex.diffraction.cbed`; the parallel-beam companion is
{doc}`reciprocal_space_and_kinematic_spots`.

## Disc Geometry

Let the zone axis point toward the gun, so the beam propagates along $-\hat{z}$, and let
the incident direction be tilted by $(\theta_x, \theta_y)$ inside the convergence cone
$|\boldsymbol{\theta}| \le \alpha$. The incident wavevector is
$\mathbf{k} = \hat{n}/\lambda$ with
$\hat{n} \simeq (\theta_x, \theta_y, -1)$ to first order in the tilt.

The Ewald condition for reflection $\mathbf{g}$ is
$|\mathbf{k}+\mathbf{g}| = |\mathbf{k}|$, that is
$2\mathbf{k}\cdot\mathbf{g} + g^{2} = 0$. Defining the deviation parameter along the zone
axis in the usual way,

$$
s_{g} = -\frac{2\mathbf{k}\cdot\mathbf{g} + g^{2}}{2|\mathbf{k}|}
        = g_{z} - \theta_x g_u - \theta_y g_v - \tfrac{1}{2}\lambda g^{2}
$$ (eq-cbed-s)

where $(g_u, g_v, g_z)$ are the components of $\mathbf{g}$ in the zone basis. At zero
tilt this reduces to $s_g = g_z - \lambda g^{2}/2$, the parallel-beam expression used by
`pytex.diffraction.kinematic`, so the two engines share one convention.

Two consequences follow directly from {eq}`eq-cbed-s`.

**The fringes are straight and perpendicular to $\mathbf{g}$.** $s_g$ is
*linear* in the tilt, with gradient $-\mathbf{g}_\perp$. Lines of constant $s_g$
are therefore straight lines perpendicular to $\mathbf{g}_\perp$, which is exactly what a
Kossel–Möllenstedt disc looks like.

**The disc radius.** A tilt $\boldsymbol{\theta}$ displaces the diffracted beam
by $\boldsymbol{\theta}/\lambda$ in reciprocal space, so the disc radius is
$\alpha/\lambda$, and on a detector calibrated by the camera constant $C = L\lambda$,

$$
R_{\text{disc}} = C\,\frac{\alpha}{\lambda} = L\alpha .
$$

The wavelength cancels: the disc size on the plate is set by the camera length and the
convergence angle alone.

**The two regimes.** Neighbouring discs touch when
$2\alpha/\lambda = |\mathbf{g}|_{\min}$. Below that the discs are separated
(*Kossel–Möllenstedt*) and each is an independent rocking curve; above it they
overlap (*Kossel*) and the overlap regions carry interference between beams. Only the
first regime supports fringe counting, which is why
`CBEDPattern.is_kossel_moellenstedt` exists and why its
`describe()` states the regime before anything else.

**Where the exact Bragg condition sits.** A disc is centred at
$s_g = -\lambda g^{2}/2$, and spans $\pm\alpha|\mathbf{g}_\perp|$. The exact Bragg
condition $s_g = 0$ therefore lies inside the disc only when
$\alpha > \lambda|\mathbf{g}|/2 = \theta_{B}$: the convergence angle must exceed the Bragg
angle. Below that the disc shows only one wing of the rocking curve, and the two branches
of fringes are unequal — a fact that matters when choosing which side to measure.

## The Two-Beam Rocking Curve

With only the transmitted and one diffracted beam coupled, the Howie–Whelan equations have
the closed solution (no absorption)

$$
I_{g}(s) = \frac{\sin^{2}\!\left(\pi t s_{\mathrm{eff}}\right)}
                  {\left(\xi_{g}s_{\mathrm{eff}}\right)^{2}},
  \qquad
  s_{\mathrm{eff}} = \sqrt{s^{2} + \xi_{g}^{-2}},
  \qquad
  I_{0} = 1 - I_{g}
$$ (eq-two-beam)

with $t$ the foil thickness and $\xi_{g}$ the extinction distance. At exact Bragg
incidence $I_{g} = \sin^{2}(\pi t/\xi_{g})$: the beams exchange intensity completely as
the thickness increases, the Pendellösung oscillation, and that is what makes
$\xi_{g}$ a measurable length rather than a bookkeeping constant. The kinematic limit is
$\xi_{g}\to\infty$, where $s_{\mathrm{eff}}\to s$ and the prefactor becomes
$(\pi t/\xi_{g})^{2}$; the kinematic approximation is already poor once
$t \gtrsim \xi_{g}/3$.

## The Extinction Distance And Its Absolute Scale

$$
\xi_{g} = \frac{\pi V_{c}\cos\theta_{B}}{\lambda\,|F_{g}|}
$$ (eq-xi)

with $V_{c}$ the unit-cell volume and $F_{g}$ the electron structure factor *in
units of length*. Getting $F_g$ onto an absolute scale is the step at which a plausible
but wrong answer is easiest to produce, so it is worth stating in full.

Electrons are scattered by the electrostatic potential, not by the charge density, so the
electron scattering factor follows from the X-ray form factor by the Mott–Bethe relation

$$
f_{e}(s) = \frac{Z - f_{x}(s)}{8\pi^{2}a_{0}\,s^{2}}
           = \frac{Z - f_{x}(s)}{41.78214\,s^{2}},
  \qquad s = \frac{\sin\theta}{\lambda} = \frac{|g|}{2}
$$ (eq-mott-bethe)

in ångström. The numerator $Z - f_{x}$ is the nuclear charge screened by the electron
cloud, which is what the incident electron sees. PyTex stores X-ray form factors in the
De Graef–McHenry form $f_{x} = Z - 41.78214\,s^{2}\sum_i a_i e^{-b_i s^{2}}$; substituting
into {eq}`eq-mott-bethe` returns exactly $\sum_i a_i e^{-b_i s^{2}}$, so the same fitted
coefficients *are* the electron scattering-factor parameters and the constant
$41.78214$ is precisely $8\pi^{2}a_{0}$. No further constant is introduced.

The structure factor then sums over the cell, with the relativistic factor
$\gamma = 1 + E/m_{0}c^{2}$ applied because the incident electron is fast — 1.20 at
100 kV, 1.39 at 200 kV, so omitting it would lengthen every extinction distance by tens of
percent:

$$
F_{g} = \gamma \sum_{j} o_{j}\,f_{e}^{j}(s)\,e^{-B_{j}s^{2}}\,
          e^{2\pi i\,\mathbf{g}\cdot\mathbf{r}_{j}} .
$$

**Validation.** Against Williams and Carter (Table 23.1), aluminium at 100 kV:

| reflection | $\{111\}$ | $\{200\}$ | $\{220\}$ |
| --- | ---: | ---: | ---: |
| tabulated (Å) | 556 | 673 | 1057 |
| PyTex (Å) | 555 | 664 | 1063 |

within 1.4 percent throughout. The fitted parametrization is least accurate at small $s$
for heavy elements, where the agreement degrades to roughly ten percent — which is one
reason CBED practice *measures* $\xi_{g}$ rather than tabulating it, as the next
section does.

## Thickness Determination

The minima of {eq}`eq-two-beam` occur where $t\,s_{\mathrm{eff}} = n$ for integer
$n$. Substituting $s_{\mathrm{eff}}^{2} = s^{2} + \xi_{g}^{-2}$ and dividing by
$n^{2}$,

$$
\left(\frac{s_{n}}{n}\right)^{2}
    = \frac{1}{t^{2}} - \frac{1}{\xi_{g}^{2}}\,\frac{1}{n^{2}} .
$$ (eq-kelly)

Plotting $(s_n/n)^{2}$ against $1/n^{2}$ therefore gives a *straight line* whose
intercept is $t^{-2}$ and whose slope is $-\xi_{g}^{-2}$. One least-squares fit returns
both unknowns, so the thickness does not inherit the error of a tabulated extinction
distance. This is the linearization of Kelly *et al.* (1975), implemented as
`thickness_from_fringe_minima`.

**Algorithm.**

1. Read the dark-fringe positions $s_n$ from one disc, on the branch that shows the more fringes.
2. Assign consecutive integers $n_0, n_0+1, …$ to them, innermost first.
3. Fit {eq}`eq-kelly` by least squares; require intercept $> 0$ and slope $< 0$, which is what makes the result physical.
4. If $n_0$ is unknown, repeat over candidate $n_0$ and take the assignment with the best $R^{2}$ among the physical fits.

**The known failure.** The order of the innermost *visible* minimum is rarely
1: the first few minima fall outside the disc when the convergence angle is small or the
foil is thin. Assuming $n_0 = 1$ then biases the thickness, and often makes the fit
outright unphysical (a positive slope, i.e.\ a negative $1/\xi_{g}^{2}$) — which is why
the routine raises rather than returning a number in that case. With three or more minima
the correct $n_0$ is identifiable from the data, because a wrong assignment curves the
plot; with exactly two, every assignment fits exactly and $n_0$ must be supplied from the
physics.

## Higher-Order Laue Zones

A zone-axis pattern is blind to the lattice repeat *along the beam*, because every
zeroth-Laue-zone reflection is perpendicular to the zone axis. HOLZ reflections restore it.
Writing $H = 1/|\mathbf{r}_{uvw}|$ for the reciprocal-lattice layer spacing along the zone
axis, the Ewald sphere of radius $1/\lambda$ intersects the $n$-th layer on a circle of
projected radius

$$
G_{n} \simeq \sqrt{\frac{2nH}{\lambda}}
$$ (eq-holz)

the small-angle approximation, accurate to better than a percent for the first few zones at
ordinary TEM voltages. Since $G_n \propto \sqrt{H}$, a HOLZ ring radius converts directly
into the lattice repeat along the beam, which is the basis of CBED lattice-parameter and
strain metrology.

Lattice centering can extinguish an entire Laue zone, so PyTex reports a ring only for those
orders that carry at least one allowed reflection; reporting a ring that cannot appear would
be worse than reporting none.

## What Is Not Implemented

- **Many-beam dynamical coupling.** Each disc is computed in its own two-beam approximation, so the discs of one simulated pattern are not mutually consistent and their relative intensities are not meaningful. At a zone axis — where many reflections are simultaneously excited — this is the two-beam picture at its worst.
- **Absorption.** The fringes therefore do not decay with thickness as they do in practice.
- **HOLZ lines** inside the bright-field disc, the sharp features used for lattice-parameter metrology. Only the ring radii are given.
- **Diffraction-group symmetry determination.** The symmetry within and between discs determines the diffraction group and hence the point group *including* the presence of a centre of symmetry, which Friedel's law hides from kinematic SAED. That determination needs full dynamical intensities and is not attempted; `CBEDPattern.describe()` says so explicitly.

## Normative And Informative Sources

- P. M. Kelly, A. Jostsons, R. G. Blake and J. G. Napier, “The determination of foil thickness by scanning transmission electron microscopy”, *Physica Status Solidi (a)* **31** (1975) 771–780 — the linearization {eq}`eq-kelly`.
- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy*, 2nd ed., Springer, 2009 — two-beam theory, extinction distances (Table 23.1), and CBED practice.
- M. De Graef, *Introduction to Conventional Transmission Electron Microscopy*, Cambridge University Press, 2003 — the Howie–Whelan equations and the scattering-factor conventions.
- J. W. Steeds, “Convergent beam electron diffraction”, in *Introduction to Analytical Electron Microscopy*, Plenum, 1979 — the diffraction-group method not implemented here.
- N. F. Mott and H. S. W. Massey, *The Theory of Atomic Collisions*, 3rd ed., Oxford, 1965 — the Mott–Bethe relation {eq}`eq-mott-bethe`.
- P. A. Doyle and P. S. Turner, “Relativistic Hartree-Fock X-ray and electron scattering factors”, *Acta Crystallographica A* **24** (1968) 390–397 — the parametrization family behind the tabulated coefficients.
