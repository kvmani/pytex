# Dynamical CBED: Many-Beam Coupling, Absorption, HOLZ Lines, and Point-Group Determination

The companion note {doc}`convergent_beam_electron_diffraction`
derives the geometry of a convergent-beam pattern and the two-beam rocking curve inside one
disc. That treatment is deliberately incomplete in four connected ways, and this note closes
them: the discs of one pattern are computed independently and so are not mutually consistent;
nothing is absorbed, so the fringes never decay; the higher-order Laue zones appear only as
ring radii and not as the sharp lines used for metrology; and the symmetry analysis that
yields the crystal point group — including the presence or absence of a centre of
symmetry — is absent altogether.

Those four are one capability rather than four. Many-beam coupling makes the pattern a
single object; absorption makes it a physically realizable one; admitting higher-order Laue
zone beams into the same beam set both produces the HOLZ lines and destroys the projection
symmetry that would otherwise make every pattern look centrosymmetric; and the symmetry that
survives that destruction *is* the diffraction group. The implementation is
`pytex.diffraction.dynamical`, `pytex.diffraction.holz` and
`pytex.diffraction.diffraction_groups`.

## The Coupled Beam Equations

Write the wavefield inside the crystal as a sum over reciprocal-lattice beams,
$\Psi(\mathbf{r}) = \sum_{g}\psi_{g}(z)\exp[i(\mathbf{k}+\mathbf{g})\cdot\mathbf{r}]$.
Substituting into the Schrödinger equation for a fast electron in the crystal potential and
neglecting backscattering (the column approximation) gives, for a parallel-sided slab,

$$
\frac{\mathrm{d}\psi}{\mathrm{d}z} = i\pi\,\mathbf{A}\,\psi,
  \qquad
  A_{gg} = 2 s_{g},
  \qquad
  A_{gh} = \nu_{g-h} \quad (g \neq h)
$$ (eq-dyn-system)

with the excitation errors $s_{g}(\boldsymbol{\theta})$ of Eq. (1) of the companion note
and the Fourier coefficients of the scaled potential

$$
\nu_{g} = \frac{\lambda F_{g}}{\pi V_{c}\cos\theta_{g}},
  \qquad |\nu_{g}| = \frac{1}{\xi_{g}} .
$$ (eq-dyn-nu)

The scaling in {eq}`eq-dyn-nu` is chosen so that the modulus of the coupling coefficient is
exactly the reciprocal of the two-beam extinction distance already validated against Williams
and Carter's Table 23.1. The absolute scale of the many-beam calculation is therefore not a
second, independent claim: it is the same claim, and it is pinned by the two-beam limit below.

Because $\mathbf{A}$ does not depend on $z$, {eq}`eq-dyn-system` integrates in closed
form,

$$
\psi(t) = \exp(i\pi\mathbf{A}t)\,\mathbf{e}_{0}
$$ (eq-dyn-propagator)

the incident boundary condition being unit amplitude in the transmitted beam. Writing
$\mathbf{A} = \mathbf{C}\,\mathrm{diag}(\gamma_{j})\,\mathbf{C}^{-1}$,

$$
\psi_{g}(t) = \sum_{j} C_{gj}\,\alpha_{j}\,e^{i\pi\gamma_{j}t},
  \qquad
  \boldsymbol{\alpha} = \mathbf{C}^{-1}\mathbf{e}_{0},
  \qquad
  I_{g} = |\psi_{g}|^{2}.
$$ (eq-dyn-bloch)

The eigenvectors are the *Bloch waves*: wavefields that propagate unchanged in shape,
each attenuated at its own rate $\mathrm{Im}\,\gamma_{j}$. The real parts of the
$\gamma_{j}$ trace the dispersion surface.

**Why the excitation amplitudes are solved for and not projected.** $\mathbf{A}$ is
Hermitian without absorption and neither Hermitian nor normal with it, so $\mathbf{C}$ is
not orthogonal in general. Obtaining $\boldsymbol{\alpha}$ by projection, $\alpha_{j} = C_{0j}^{*}$, is the classic silent error in a Bloch-wave implementation: it produces
rocking curves of the right shape with the wrong contrast, and it breaks intensity
conservation in a way that is only visible if one checks for it.

## Three Exact Properties

**The two-beam limit.** For a single reflection, $\mathbf{A} = s\mathbf{I} + \mathbf{B}$, where $\mathbf{B}$ is the traceless matrix with diagonal $(-s, +s)$ and
both off-diagonal entries equal to $\nu = 1/\xi_{g}$. It has eigenvalues
$\pm s_{\mathrm{eff}}$, $s_{\mathrm{eff}} = \sqrt{s^{2}+\xi_{g}^{-2}}$, and the exponential of a traceless $2\times 2$ matrix is
$\cos(\pi s_{\mathrm{eff}}t)\mathbf{I} + i\sin(\pi s_{\mathrm{eff}}t)\mathbf{B}/ s_{\mathrm{eff}}$, so

$$
I_{g} = \frac{\sin^{2}(\pi t s_{\mathrm{eff}})}{(\xi_{g}s_{\mathrm{eff}})^{2}}
$$

which is the Howie–Whelan expression the companion note derives. This is asserted to
machine precision, and it pins the diagonal convention, the off-diagonal scale and the
factor $i\pi$ simultaneously.

**Unitarity.** Without absorption $\mathbf{A}$ is Hermitian, so
$\exp(i\pi\mathbf{A}t)$ is unitary and $\sum_{g}I_{g} = 1$ at every thickness and every
incident direction. This is the only exact global check available on a many-beam
calculation, and it is lost the moment absorption is switched on — which is a reason to
run it once with absorption off.

**Normal absorption is a scalar.** The mean absorptive coefficient
$i/\xi'_{0}$ sits on *every* diagonal element, so it commutes out of the exponential:

$$
\exp\left[i\pi\left(\mathbf{A} + \tfrac{i}{\xi'_{0}}\mathbf{I}\right)t\right]
    = e^{-\pi t/\xi'_{0}}\exp(i\pi\mathbf{A}t),
  \qquad
  I_{g} \to e^{-2\pi t/\xi'_{0}}\,I_{g}.
$$ (eq-dyn-normal)

It can therefore change no relative intensity, no fringe position and no symmetry. Every
qualitative effect of absorption comes from the *off-diagonal* absorptive terms.

## Absorption as an Imaginary Optical Potential

Electrons leave the coherent elastic wavefield by processes the calculation does not follow:
thermal diffuse scattering above all, then plasmon and core losses. The standard
representation adds an imaginary part to each Fourier coefficient,

$$
\frac{1}{\xi_{g}} \;\longrightarrow\; \frac{1}{\xi_{g}} + \frac{i}{\xi'_{g}}
$$

so that {eq}`eq-dyn-system` becomes $A_{gg} = 2s_{g} + i/\xi'_{0}$ and $A_{gh} = \nu_{g-h}(1 + i r)$ with $r = \xi_{g}/\xi'_{g}$.

**What is derived and what is assumed.** The *structure* is not an
approximation: an imaginary optical potential is the correct representation of loss from a
coherent wavefield, and the observable consequences below follow from the eigenvector
structure rather than being applied afterwards. The *magnitudes* are phenomenological.
`pytex` carries them as the two ratios $\xi_{0}/\xi'_{0}$ and $\xi_{g}/\xi'_{g}$,
whose customary working value for a metal near 100–200 kV is about
$0.1$ (Hirsch *et al.*, *Electron Microscopy of Thin Crystals*, Ch. 12). A
first-principles absorptive form factor — the Einstein-model thermal-diffuse integral of
Hall and Hirsch, as parametrized by Bird and King — is not implemented, and the
documentation says so rather than implying a computed constant.

**Positivity.** The absorptive matrix has eigenvalues $\xi_{0}'^{-1} \pm \xi_{g}'^{-1}$ in the two-beam case. Since $\xi_{0} < \xi_{g}$ always, requiring the
reflection ratio not to exceed the mean ratio is sufficient for both to be positive. A
larger reflection ratio would produce a Bloch wave that *gains* intensity with depth,
which no absorption process can do; the constructor refuses it.

**Anomalous absorption.** The two Bloch waves of a two-beam calculation have their
intensity maxima on and between the atomic planes respectively. The first is absorbed
strongly, the second weakly, and which is preferentially excited depends on the sign of
$s$. The consequence is the Hashimoto–Howie–Whelan result: with absorption, the
bright-field rocking curve becomes asymmetric about $s = 0$ while the dark-field one
remains exactly symmetric. `pytex` asserts both halves of that statement numerically;
they are a stronger test of the absorption implementation than any single number, because
they are properties no incorrect implementation reproduces by accident.

## Friedel's Law as a Propagator Theorem

For a real potential $\nu_{-g} = \nu_{g}^{*}$, so $\mathbf{A}$ is always Hermitian. It is
*symmetric* only when every included $\nu_{g}$ is real, which happens exactly when the
structure sampled by the beam set has a centre of symmetry with the origin on it.

Relabelling the beams by $g \mapsto -g$ and the incident tilt by $\boldsymbol{\theta} \mapsto -\boldsymbol{\theta}$ maps the diagonal $2s_{g}(\boldsymbol{\theta}) \mapsto 2s_{-g}(-\boldsymbol{\theta})$ and the off-diagonal $\nu_{g-h} \mapsto \nu_{h-g} = \nu_{g-h}^{*}$. In the zeroth Laue zone $s_{-g}(-\boldsymbol{\theta}) = s_{g}(\boldsymbol{\theta})$ exactly, so the relabelled matrix is
$\mathbf{A}^{\mathsf{T}}$, the propagator becomes $M^{\mathsf{T}}$, and

$$
I_{g}(\boldsymbol{\theta}) = I_{-g}(-\boldsymbol{\theta})
  \iff |M_{g0}| = |M_{0g}|
  \iff \mathbf{M} \text{ symmetric}
  \iff \text{the sampled structure is centrosymmetric.}
$$ (eq-dyn-friedel)

Friedel's law is thus recovered as a theorem about the propagator rather than as a kinematic
accident, and {eq}`eq-dyn-friedel` is the mechanism behind the point-group determination.

**The word “sampled” carries the whole difficulty.** A beam set confined to the
zeroth Laue zone samples the potential *projected* along the beam, and that projection
is frequently centrosymmetric when the crystal is not. Zincblende down $[111]$ is the
standard example: every ZOLZ coefficient is real, so a projection calculation reports
Friedel's law to machine precision and cannot see the polarity of the crystal. Admitting the
first-order Laue zone gives the coefficients phases, breaks the symmetry of the propagator,
and separates the two discs of a $\pm\mathbf{g}$ pair by tens of percent. `pytex`
pins this with a controlled pair: zincblende GaAs and a rocksalt structure built on the same
lattice from the same two species, differing only by the position of the second sublattice
and therefore only by the centre of symmetry. The Friedel violation separates them by more
than three orders of magnitude.

*Higher-order Laue zone interaction is not a refinement here; it is the entire
mechanism.* A symmetry conclusion drawn from a zeroth-Laue-zone calculation is worthless, and
`BeamSet.holz_mask` exists so that this can be checked rather than assumed.

## Beam Selection

A many-beam calculation is defined only once the beams are named. `pytex` selects a
reflection when it comes within a stated window of the Bragg condition *for some
incident direction in the illumination cone*. Because $s_{g}$ is affine in the tilt,

$$
\min_{|\boldsymbol{\theta}| \le \alpha} |s_{g}(\boldsymbol{\theta})|
    = \max\bigl(0,\; |s_{g}(\mathbf{0})| - \alpha|\mathbf{g}_{\perp}|\bigr)
$$ (eq-dyn-selection)

which is evaluated in closed form. The distinction matters: a HOLZ reflection is far from
Bragg at the centre of the pattern and exactly at Bragg somewhere inside the bright-field
disc — that locus is a HOLZ line — so selecting on the zero-tilt excitation error alone
would discard every HOLZ reflection and, with it, the symmetry-breaking mechanism of the
previous subsection.

**Cost.** Every selected beam is solved exactly; there is no Bethe perturbation of
weak beams, so the cost is $O(m n^{3})$ for $m$ incident directions and $n$ beams. The
economy that works is a tighter window in {eq}`eq-dyn-selection`, which removes beams that
were barely coupled, rather than coarser tilt sampling, which degrades the fringe positions
the calculation exists to predict.

## HOLZ Lines

A higher-order Laue zone reflection is far from the Bragg condition at the centre of the
pattern and exactly at it somewhere inside the bright-field disc. Because $s_g$ is affine in
the incident tilt, that locus is a straight line,

$$
\boldsymbol{\theta}\cdot\hat{\mathbf{g}}_{\perp} = d_{g},
  \qquad
  d_{g} = \frac{g_{z} - \tfrac{1}{2}\lambda|\mathbf{g}|^{2}}{|\mathbf{g}_{\perp}|}
$$ (eq-holz-line)

with unit normal $\hat{\mathbf{g}}_{\perp}$ and signed distance $d_g$ from the pattern
centre; it crosses the bright-field disc when $|d_g| \le \alpha$. The line appears dark
(*deficiency*) in the direct disc and bright (*excess*) at the same incident tilts in
the disc of $\mathbf{g}$ itself. Nothing beyond {eq}`eq-dyn-selection` is approximated, so
line positions are exact to the accuracy of $s_g$; `pytex` checks them by asking the
dynamical module for $s_g$ at points on the line and requiring zero to $10^{-15}$.

**Sharpness.** The rocking curve has its first zero at $|s| = 1/t$, and $s$ varies
across the disc at the rate $|\mathbf{g}_{\perp}|$ per radian, so the line's angular
half-width is $\Delta\theta = 1/(t|\mathbf{g}_{\perp}|)$. It narrows in proportion to the foil
thickness — HOLZ metrology wants a *thick* specimen, the opposite of the usual thin-foil
instinct.

**The metrology, and the trap in it.** Scaling the lattice by $1+\varepsilon$ shrinks
every $\mathbf{g}$ by the same factor. Substituting into {eq}`eq-holz-line`,

$$
d_{g}(\varepsilon, \lambda) = \frac{g_{z}}{|\mathbf{g}_{\perp}|}
    - \frac{\lambda|\mathbf{g}|^{2}}{2(1+\varepsilon)|\mathbf{g}_{\perp}|},
  \qquad
  \left.\frac{\partial d_{g}}{\partial\varepsilon}\right|_{0}
    = \frac{\lambda|\mathbf{g}|^{2}}{2|\mathbf{g}_{\perp}|}
$$

and the wavelength enters the same term with the opposite sign. Setting $\lambda \to \lambda(1+\varepsilon)$ therefore cancels a lattice strain $\varepsilon$ exactly, at every
reflection simultaneously. *A fractional change in lattice parameter and a fractional
change in wavelength are indistinguishable from HOLZ line positions.* This is not a limitation
of the model: it is why quantitative HOLZ metrology begins by calibrating the accelerating
voltage against a standard of known lattice parameter, and why an uncalibrated measurement of a
lattice parameter is a measurement of the high-tension supply. `pytex` asserts the
degeneracy to $10^{-16}$ rather than describing it.

**Why intersections are measured, not lines.** For nickel down $[001]$ at 200 kV in a
$1000$ Å foil, the best single line moves $0.059$ mrad per unit strain against a
half-width of $0.21$ mrad: a strain of $3.6\times10^{-3}$ shifts it by its own width, far
short of the $10^{-4}$ the technique is known for. The crossing of two lines meeting at angle
$\phi$, however, moves as $1/\sin\phi$ times faster, and the best near-parallel pair
resolves $6.3\times10^{-5}$. That is where the sensitivity comes from, and
`HOLZLinePattern.intersections` reports the amplification explicitly.

## Diffraction Groups

Fix a beam direction $\mathbf{b}$. Every operator $S$ of the crystal point group falls into
one of three cases: it fixes $\mathbf{b}$, it reverses it, or it does neither and relates
different patterns. The first two contribute the restriction $T = S|_{\perp}$ to the plane
normal to the beam, the second tagged with the reciprocity flag $R$ because such an operator
is a symmetry only in combination with the reciprocity theorem. The map $S \mapsto (S|_\perp, \text{tag})$ is a homomorphism onto a subgroup of $G_2 \times \mathbb{Z}_2$ with $G_2$ one
of the ten plane point groups, and those subgroups are of three kinds: no tagged element (10 of
them), the full direct product (10, written with the suffix $1_R$), and the graph of a
surjection onto $\mathbb{Z}_2$ (11, with $_R$ on individual generators). That is Buxton's
31, obtained by construction rather than transcription, and `pytex` verifies the count
and the membership against the published set.

**The two observables.** A tagged element requires reciprocity, which relates a point in
one disc to a point in another at an incident direction outside the illumination cone. Only
untagged elements are rigid symmetries of the recorded pattern:

$$
\text{WP} = \{\,T : (T, \text{untagged}) \in D\,\}.
$$

Inside the direct disc the reciprocity displacement is proportional to $\mathbf{g}_\perp$ and
therefore vanishes, leaving only reciprocity's own inversion of the incident direction, so a
tagged element acts there as $-T$:

$$
\text{BF} = \varphi(D),
  \qquad \varphi(T, \text{untagged}) = T,
  \qquad \varphi(T, \text{tagged}) = -T .
$$

$\varphi$ is a homomorphism because $-1$ is central in two dimensions. The familiar
consequences follow: a two-fold axis perpendicular to the beam gives a bright-field mirror the
whole pattern does not have ($m_R$); a mirror perpendicular to the beam gives a bright-field
two-fold alone ($1_R$); and $\bar{6}m2$ down its three-fold shows a *six*-fold direct
disc over a $3m$ whole pattern ($3m1_R$).

**The centre of symmetry.** The element $2_R$ requires an operator acting as $-1$ on
the beam direction *and* as $-1$ on the transverse plane, which is the inversion and
nothing else. Hence $2_R \in D$ at every beam direction of a centrosymmetric crystal and at
none of an acentric one — an exact correspondence, verified over all 32 point groups. Knowing
only whether $2_R$ is present therefore partitions the 32 point groups into exactly 21 and 11,
which is the arithmetic of the whole technique.

**What $\pm\mathbf{G}$ is, and is not.** Buxton's $\pm\mathbf{G}$ observation
compares the $+\mathbf{g}$ and $-\mathbf{g}$ *dark-field* discs, each recorded with its
own reflection at the Bragg condition — two exposures at different specimen tilts, related by
reciprocity. It is *not* a two-fold rotation of a single zone-axis pattern. Taking it to be
one gives a test that fails, and the reason is visible in the excitation errors:
$s_{-g}(-\boldsymbol{\theta}) - s_{g}(\boldsymbol{\theta}) = -2g_{z}$, which vanishes only in
the zeroth Laue zone. Once higher-order beams are admitted the two-fold is broken for a
centrosymmetric crystal too, and numerically the residual *grows* with the beam set for
centric and acentric structures alike — so it is physics, not truncation. `pytex`
therefore does not measure $\pm\mathbf{G}$ from a zone-axis simulation, and says so instead of
reporting a number that would sometimes be wrong. The determination does not need it: at any
zone whose candidate diffraction groups differ in bright-field or whole-pattern symmetry, those
two settle the question.

**The worked case.** Down $[001]$, zincblende $\bar{4}3m$ has diffraction group
$4_Rmm_R$ — a four-fold bright-field disc over a merely two-fold whole pattern — while a
centrosymmetric structure on the same lattice has $4mm1_R$, four-fold in both. Simulating both
with the coupled method and reading the symmetry back separates them with residuals of $0.00$
against $0.32$, and inverting the observation returns $\{\bar{4}2m, \bar{4}3m\}$ for the
polar structure with the centrosymmetry verdict *false*. Confine the beam set to the
zeroth Laue zone and the same crystal reports the four-fold whole-pattern symmetry of a
centrosymmetric one: the projected potential of zincblende down $[001]$ *is*
centrosymmetric. That is the single most important caveat in CBED symmetry work, and it is
enforced rather than described — `symmetry_observations` refuses a projection
calculation unless asked a second time.

**Not implemented.** Buxton's table also lists dark-field and $\pm\mathbf{G}$
symmetries for reflections lying on symmetry lines, recorded at their own Bragg conditions.
These sharpen determinations that the two implemented observations leave open — $4mm1_R$
against $4mm$, for instance. They are absent here, and the report says so and recommends a
second zone axis instead, naming the tool that finds one.
