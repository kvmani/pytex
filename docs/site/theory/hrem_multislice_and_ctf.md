# High-Resolution TEM Multislice Simulation and Contrast Transfer Function

High-Resolution Transmission Electron Microscopy (HRTEM or HREM) forms images by phase
interference between the forward-transmitted beam and diffracted electron waves. Unlike
selected-area electron diffraction (SAED) or convergent-beam electron diffraction (CBED), which
record the intensity in reciprocal space, HRTEM records the interference pattern in real
space on the detector plane.

Because electrons interact strongly with atomic Coulomb potentials via Rutherford and dynamical
scattering, the exit wave emerging from the specimen carries complex phase modulations. The
objective lens transfer function then modulates and dampens these spatial frequencies through
geometric aberrations, aperture cutoffs, and partial coherence damping envelopes.

This note develops the electron optics contrast transfer function (CTF), higher-order aberration
corrections, chromatic and spatial coherence damping envelopes, the Cowley–Moodie multislice
propagation algorithm, and Thon-ring power spectra. The canonical Python implementation is
`pytex.diffraction.hrem` and the multislice bridge is `pytex.adapters.abtem`.

---

## Relativistic Electron Kinematics and Interaction

Electrons accelerated through high potential $V$ (typically $60\text{--}300\text{ kV}$) attain
relativistic velocities. The relativistic kinetic energy is $E = e V$. The de Broglie wavelength
$\lambda$ in vacuum is:

$$
\lambda = \frac{h}{\sqrt{2 m_0 e V \left(1 + \frac{e V}{2 m_0 c^2}\right)}}
$$ (eq-hrem-lambda)

where $h$ is Planck's constant, $m_0$ is the electron rest mass, and $c$ is the speed of light.
At $200\text{ kV}$, $\lambda \approx 0.02508\text{ \AA}$; at $300\text{ kV}$,
$\lambda \approx 0.01969\text{ \AA}$.

The relativistic interaction parameter $\sigma$ links the projected Coulomb potential
$V_{\mathrm{proj}}(x, y) = \int V(x, y, z)\,dz$ to the phase shift experienced by the electron
wave:

$$
\sigma = \frac{2\pi}{\lambda V} \left(\frac{m_0 c^2 + e V}{2 m_0 c^2 + e V}\right)
       = \frac{2\pi m_0 e \lambda}{h^2} \left(1 + \frac{e V}{m_0 c^2}\right)
$$ (eq-hrem-sigma)

At $200\text{ kV}$, $\sigma \approx 0.000729\text{ V}^{-1}\text{\AA}^{-1}$
($0.729\text{ V}^{-1}\text{nm}^{-1}$).

---

## The Objective Lens Contrast Transfer Function (CTF)

In the linear imaging approximation for a weak phase object, the electron wave entering the
objective lens is modulated by the lens phase aberration function $\chi(\mathbf{q})$. For an
azimuthally symmetric lens, the wave aberration $\chi(q)$ as a function of spatial frequency
$q = 1/d$ is:

$$
\chi(q) = \pi \Delta f \lambda q^2 + \tfrac{1}{2} \pi C_s \lambda^3 q^4 + \tfrac{1}{3} \pi C_5 \lambda^5 q^6
$$ (eq-hrem-chi)

where:
- $\Delta f$ is the defocus (underfocus is conventionally negative, $\Delta f < 0$ in the Kirkland
  convention, while Scherzer overfocus is positive; PyTex adheres to the Kirkland standard where
  underfocus is negative),
- $C_s \equiv C_{30}$ is the third-order spherical aberration coefficient (typically
  $0.5\text{--}2.0\text{ mm}$ for uncorrected instruments; $< 10\text{ }\mu\text{m}$ for
  $C_s$-corrected instruments),
- $C_5 \equiv C_{50}$ is the fifth-order spherical aberration coefficient.

### Non-round aberrations, and why a single radial cut is not the lens

Equation {eq}`eq-hrem-chi` describes a *round* lens. Its three terms depend on $q$ alone, so a
single radial cut of $\chi$ describes the instrument completely and every direction in the image is
resolved equally. Real objective lenses are not round. In the Krivanek numbering $C_{nm}$, where
$n$ is the order in $q$ and $m$ the azimuthal multiplicity, the terms with $m > 0$ carry an
azimuthal dependence:

$$
\chi(q, \theta) = \underbrace{\pi \Delta f \lambda q^{2}
  + \tfrac{1}{2}\pi C_{s}\lambda^{3}q^{4}
  + \tfrac{1}{3}\pi C_{5}\lambda^{5}q^{6}}_{\text{round}}
  + \underbrace{\pi C_{12}\lambda q^{2}\cos 2(\theta - \varphi_{12})
  + \tfrac{2}{3}\pi C_{21}\lambda^{2}q^{3}\cos(\theta - \varphi_{21})
  + \tfrac{2}{3}\pi C_{23}\lambda^{2}q^{3}\cos 3(\theta - \varphi_{23})}_{\text{non-round}}
$$ (eq-hrem-chi-azimuthal)

with $C_{12}$ the two-fold astigmatism, $C_{21}$ the axial coma and $C_{23}$ the three-fold
astigmatism (trefoil), each at its own azimuth $\varphi_{nm}$.

Three consequences follow, and the third is why the workbench reports what it does.

**Two-fold astigmatism is a directional defocus.** Setting $\theta = \varphi_{12}$ in
{eq}`eq-hrem-chi-azimuthal` collapses its cosine to $+1$, and the $C_{12}$ term becomes
$\pi C_{12} \lambda q^{2}$ — algebraically indistinguishable from adding $C_{12}$ to the defocus.
At $\theta = \varphi_{12} + 90^{\circ}$ the cosine is $-1$ and the same term subtracts it. A lens
carrying $C_{12}$ therefore behaves as two different round lenses, at defocus $\Delta f + C_{12}$
and $\Delta f - C_{12}$, along two orthogonal directions. This identity is what
`tests/unit/test_hrem_azimuthal_ctf.py` pins, to a relative tolerance of $10^{-12}$, and it is the
provenance of the worked example: no recorded output is involved.

**Each term has a period fixed by its multiplicity.** $\chi$ repeats every $180^{\circ}$ in
$\theta$ for $C_{12}$, every $360^{\circ}$ for $C_{21}$ and every $120^{\circ}$ for $C_{23}$.
The single-fold period of coma is the reason it transfers differently in opposite directions,
displacing image detail asymmetrically rather than merely blurring it.

**Correction moves the problem rather than removing it.** A corrector drives $C_{s}$ toward zero,
which is what extends the passband past the Scherzer boundary
({eq}`eq-hrem-scherzer-res` and the NCSI discussion below). What then limits the instrument
is the residual it leaves: the non-round terms, and $C_{5}$. Quoting a point resolution from one
radial cut of such a lens states the resolution of one direction and says nothing about the others,
which is why PyTex reports the **range of point resolution over azimuth** and its spread — the
resolution anisotropy — alongside any cut, and why both `describe()` surfaces name the residual
terms that are present.

The coherence envelopes of the next section are Frank's isotropic forms, derived for a round lens:
the spatial envelope uses the radial derivative of the round part of $\chi$. The anisotropy PyTex
reports is therefore that of the transfer oscillation $\sin\chi$, not of the damping. This is a
stated limitation, not an approximation whose error is quantified here.

### Classical Scherzer Imaging

In an uncorrected transmission electron microscope ($C_s > 0$), Scherzer (1949) showed that an
optimal passband of uniform negative phase shift $\sin\chi(q) \approx -1$ is achieved at the
Scherzer defocus:

$$
\Delta f_{\mathrm{Sch}} = -1.2 \sqrt{C_s \lambda}
$$ (eq-hrem-scherzer-defocus)

The corresponding interpretable point resolution $d_{\mathrm{Sch}}$ corresponds to the first zero
crossing of the transfer function:

$$
d_{\mathrm{Sch}} \approx 0.64 \left(C_s \lambda^3\right)^{1/4}
$$ (eq-hrem-scherzer-res)

For example, at $200\text{ kV}$ ($\lambda = 0.02508\text{ \AA}$) and $C_s = 1.0\text{ mm}$,
$\Delta f_{\mathrm{Sch}} \approx -601\text{ \AA}$ and $d_{\mathrm{Sch}} \approx 2.44\text{ \AA}$.

### Double-Corrected TEM and Negative $C_s$ Imaging (NCSI)

Modern aberration-corrected microscopes incorporate multipole hexapole or quadrupole-octupole
correctors (Haider et al., 1998; Kabius et al., 2009) that cancel both $C_s$ and chromatic
aberration $C_c$.

1. **$C_s$-Corrected Regime:** By tuning the hexapole corrector, $C_s$ is reduced to near zero
   ($|C_s| < 5\text{ }\mu\text{m}$). The passband extends to much higher spatial frequencies,
   and the point resolution approaches the information limit.
2. **Negative $C_s$ Imaging (NCSI):** Urban et al. (2009) demonstrated that deliberately setting
   a small negative spherical aberration ($C_s \approx -10\text{ to }-30\text{ }\mu\text{m}$)
   combined with a small positive overfocus ($\Delta f \approx +30\text{ to }+80\text{ \AA}$)
   produces maximum phase contrast where atomic columns appear bright on a dark background.
   This regime yields high contrast for light elements (such as oxygen or carbon) adjacent to heavy
   metal columns.

---

## Partial Coherence Damping Envelopes

Real electron beams are neither perfectly monochromatic nor perfectly parallel. Temporal and
spatial incoherence act as multiplicative damping envelopes on the transfer function (Frank, 1973;
Kirkland, 2010):

$$
T(q) = E_c(q) \cdot E_s(q) \cdot A(q) \cdot \sin\chi(q)
$$ (eq-hrem-transfer)

### Temporal Coherence Envelope $E_c(q)$

Variations in accelerating voltage $\Delta V$, objective lens current fluctuations $\Delta I$,
and the natural energy spread of the field-emission electron gun $\Delta E$ produce a focal spread
$\Delta$:

$$
\Delta = C_c \sqrt{\left(\frac{\Delta V}{V}\right)^2 + 4\left(\frac{\Delta I}{I}\right)^2 + \left(\frac{\Delta E}{E}\right)^2}
$$ (eq-hrem-focal-spread)

The temporal envelope is Gaussian in $q^2$:

$$
E_c(q) = \exp\left(-\tfrac{1}{2} \pi^2 \lambda^2 \Delta^2 q^4\right)
$$ (eq-hrem-ec)

Note that $E_c(q)$ depends on $q^4$, imposing an absolute information limit beyond which no phase
information can be transferred, regardless of defocus. In double-corrected microscopes,
$C_c$-correction reduces $\Delta$ from $\sim 30\text{ \AA}$ to $< 10\text{ \AA}$, dramatically
extending the temporal envelope.

### Spatial Coherence Envelope $E_s(q)$

Due to the finite source size, the illumination cone has an effective convergence semi-angle
$\alpha_s$ (typically $0.1\text{--}0.5\text{ mrad}$). The spatial coherence envelope dampens
spatial frequencies according to the slope of the wave aberration $\nabla_q \chi$:

$$
E_s(q) = \exp\left(-\pi^2 \alpha_s^2 \left[\Delta f q + C_s \lambda^2 q^3\right]^2\right)
$$ (eq-hrem-es)

Near Scherzer defocus or when $C_s \approx 0$ and $\Delta f \approx 0$, the bracket
$[\Delta f q + C_s \lambda^2 q^3]$ remains small over a wide frequency band, minimizing spatial
damping.

### Objective Aperture Mask $A(q)$

An objective aperture inserted in the back focal plane blocks electrons scattered beyond a
cutoff angle $\alpha_{\mathrm{obj}}$:

$$
A(q) = \begin{cases}
1, & q \le \frac{\alpha_{\mathrm{obj}}}{\lambda} \\
0, & q > \frac{\alpha_{\mathrm{obj}}}{\lambda}
\end{cases}
$$ (eq-hrem-aperture)

---

## The Multislice Algorithm

For specimens thicker than a few nanometres, multiple scattering and dynamical interactions
invalidate the single-phase-object approximation. Cowley & Moodie (1957) formulated the multislice
method, which slices the atomic potential along the beam propagation direction $z$ into thin
layers of thickness $\Delta z \approx 1\text{--}2\text{ \AA}$.

Each step consists of two operations:

1. **Phase Transmission:** The electron wave $\psi_n(x, y)$ interacting with slice $n$ undergoes
   a transmission phase shift:
   $$
   t_n(x, y) = \exp\left(i \sigma v_n(x, y)\right)
   $$ (eq-hrem-trans)
   where $v_n(x, y) = \int_{z_n}^{z_n+\Delta z} V(x, y, z) \, dz$ is the projected potential of
   atoms within the slice, computed via parameterized Kirkland atomic scattering factors.

2. **Fresnel Propagation:** The transmitted wave propagates through free space of distance
   $\Delta z$ to the next slice:
   $$
   \psi_{n+1}(x, y) = \left[\psi_n(x, y) \cdot t_n(x, y)\right] * p(x, y, \Delta z)
   $$ (eq-hrem-prop)
   In reciprocal space, this convolution becomes a simple multiplication by the propagator:
   $$
   P(q_x, q_y, \Delta z) = \exp\left(-i \pi \lambda (q_x^2 + q_y^2) \Delta z\right)
   $$ (eq-hrem-prop-k)

After traversing all slices, the complex wave at the exit surface is the **exit wave**
$\psi_{\mathrm{exit}}(x, y)$.

### Image Synthesis

The objective lens acts on the exit wave in reciprocal space:

$$
\Psi_{\mathrm{image}}(\mathbf{q}) = \mathcal{F}\{\psi_{\mathrm{exit}}(\mathbf{r})\} \cdot T(\mathbf{q}) \cdot e^{-i \chi(\mathbf{q})}
$$ (eq-hrem-image-wave)

The recorded real-space image intensity is the squared modulus of the image wave:

$$
I(x, y) = |\mathcal{F}^{-1}\{\Psi_{\mathrm{image}}(\mathbf{q})\}|^2
$$ (eq-hrem-intensity)

---

## Power Spectra and Thon Rings

For an amorphous specimen (such as carbon foil), the scattering potential is spatially disordered
with an essentially white Fourier spectrum. The 2D fast Fourier transform (FFT) power spectrum of
the image intensity $I(x, y)$ displays concentric rings of zero intensity known as **Thon rings**
(Thon, 1966):

$$
|\mathcal{F}\{I(x, y) - \langle I \rangle\}|^2 \propto |T(\mathbf{q})|^2 = E(q)^2 \sin^2\chi(\mathbf{q})
$$ (eq-hrem-thon)

The radii of the Thon rings correspond directly to the zero crossings of $\sin\chi(q)$, where
$\chi(q) = n\pi$. Ellipticity in the Thon rings directly reveals 2-fold astigmatism $C_{12}$,
while the radial decay of ring contrast measures the coherence envelopes $E_c(q)$ and $E_s(q)$.

---

## PyTex Implementation Architecture

The HREM capability is partitioned into decoupled layers:

1. **`pytex.diffraction.hrem`:**
   - Canonical physical constants and relativistic transformations (`relativistic_wavelength_angstrom`, `relativistic_interaction_parameter_inv_v_angstrom`).
   - `DoubleCorrectionMode`: enumeration of imaging regimes (`UNCORRECTED`, `CS_CORRECTED`, `DOUBLE_CORRECTED`, `NCSI`).
   - `MicroscopeAberrations`: self-describing optics state with higher-order coefficients, Scherzer evaluation, envelopes, and `.describe()`.
   - `CTF1D`: 1D radial transfer evaluation, zero crossings, and information limit.
   - `AtomicSnapshot`: atomic coordinate model supporting crystalline slabs, vacancies, Volterra dislocations, and amorphous packing, with XYZ and ASE conversion.
   - `HREMSimulationResult`: 2D micrograph intensity, exit wave, FFT power spectrum, contrast metrics, line profiles, base64 rendering, and `.describe()`.
   - `pure_python_phase_object_simulation`: fallback phase-object propagation.

2. **`pytex.adapters.abtem`:**
   - Adapter bridging PyTex structures with the open-source `abtem` multislice engine and `ase.Atoms`.
   - Automatic dispatch: uses multislice if `abtem` is present; falls back gracefully to pure-Python phase-object simulation if absent.

3. **`pytex.cli`:**
   - Command-line interfaces: `pytex hrem simulate` and `pytex hrem ctf`.

4. **`pytex.app.services.tem_hrem` and `src/pytex/app/static/js/panels/hrem.js`:**
   - Interactive GUI submodule integrated into the TEM Analysis workspace.

---

## References

### Adapter convention and validation

`to_abtem_ctf` preserves every coefficient represented by `MicroscopeAberrations`:
defocus, twofold astigmatism, coma, trefoil, and third- and fifth-order spherical
aberration. Amplitudes are converted to angstroms and azimuths from degrees to radians.
PyTex's defocus is the coefficient multiplying the positive quadratic term in
`wave_aberration`, so it maps directly to abTEM's `C10`; abTEM's `defocus` alias has
the opposite sign. These conventions follow the
[abTEM CTF documentation](https://abtem.readthedocs.io/en/latest/user_guide/walkthrough/contrast_transfer_function.html).

`tests/unit/test_abtem_aberration_parity.py` compares the complex phase transfer
over a frequency/azimuth grid for coma, trefoil, tiny astigmatism and combined
aberrations. It also checks coefficient and angle conversion directly. These optional
integration tests require abTEM; the analytic core tests and
{doc}`computed HREM examples </examples/generated/hrem-simulation-and-ctf>` do not.
This establishes optical coefficient agreement, not equality of the phase-object and
multislice specimen models or validation against measured micrographs.

For imported structures, `AtomicSnapshot.from_xyz` accepts a multiline single-frame XYZ
string or a file path. A blank comment line is valid. Atom counts must match the frame,
and coordinates must be finite. Supply the simulation cell explicitly when its dimensions
are known: a plain XYZ file carries Cartesian angstrom coordinates, not a periodic lattice.
`tests/unit/test_hrem_xyz.py` checks text/file equivalence and malformed input independently
of ASE and operating-system filename limits.

### Scientific sources

1. Scherzer, O. (1949). The theoretical resolution limit of the electron microscope. *J. Appl. Phys.* **20**, 20–29.
2. Cowley, J. M. & Moodie, A. F. (1957). The scattering of electrons by atoms and crystals. I. A new theoretical approach. *Acta Crystallogr.* **10**, 609–619.
3. Thon, F. (1966). Zur Defokussierungsabhängigkeit des Phasenkontrastes im elektronenmikroskopischen Abbildung. *Z. Naturforsch. A* **21**, 476–478.
4. Frank, J. (1973). An envelope for the transfer function of the electron microscope. *Optik* **38**, 519–536.
5. Haider, M., Uhlemann, S., Schwan, E., Rose, H., Kabius, B. & Urban, K. (1998). Electron microscopy image with a spherical-aberration-corrected objective lens. *Nature* **392**, 768–770.
6. Kabius, B., Hartel, P., Haider, M., Müller, H., Uhlemann, S., Loebau, F., Zach, J. & Rose, H. (2009). First application of an un-monochromated, Cs and Cc corrected transmission electron microscope. *Microsc. Microanal.* **15** (Suppl 2), 1150–1151.
7. Urban, K. W., Jia, C. L., Houben, L., Lentzen, M., Mi, S. B. & Thust, A. (2009). Negative spherical aberration imaging in transmission electron microscopy. *Phil. Trans. R. Soc. A* **367**, 3735–3753.
8. Kirkland, E. J. (2010). *Advanced Computing in Electron Microscopy*, 2nd ed., Springer.
9. Madsen, J., Susi, T. & Sader, K. et al. (2021). abTEM: An open-source framework for simulation of transmission electron microscopy. *ChemPhysChem* **22**, 1–13.
