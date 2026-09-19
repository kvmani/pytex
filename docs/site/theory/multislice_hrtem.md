# Multislice HRTEM Simulation: Theory and Mathematics

An HRTEM image of a crystal thicker than a few ångströms cannot be read as a projection of its
atoms. The fast electron is scattered many times on its way through the specimen, and the wave
that leaves the exit surface is the outcome of that multiple (dynamical) scattering. The
multislice method of Cowley and Moodie (1957) computes that exit wave by cutting the specimen into
thin slices and alternating two exactly solvable steps: a phase shift in each slice and free
propagation between slices. The objective lens then forms the image from the exit wave.

This note derives the method as `pytex.diffraction.multislice` implements it. The implementation
follows abTEM (Madsen & Susi, 2021) step for step, so that the two agree number for number, and
it is validated against abTEM, against Bloch waves and against closed forms (see
[Validation](#validation)). The objective-lens optics themselves — the wave aberration
$\chi$, the Scherzer and NCSI conditions, the coherence envelopes and Thon rings — are developed in
{doc}`hrem_multislice_and_ctf`; this note uses them. The algorithm, its steps, costs and failure
modes are stated in {doc}`../algorithms/multislice_hrtem`, and the numbers are computed live in
{doc}`the multislice worked examples </examples/generated/multislice-hrtem>`.

**Notation.** Symbols follow the
{doc}`terminology and symbol registry <../standards/terminology_and_symbol_registry>`. Here
$\mathbf r = (x, y)$ is a position across the beam, $z$ the depth along it, and $\mathbf g$ a
*lateral spatial frequency* in Å⁻¹ with no factor $2\pi$ — at a Bragg beam, the reciprocal-lattice
vector. (The CTF note writes the same frequency as $q$.) The interaction parameter is written
$\sigma_e$, the transmission function $\tau_n$ and the wave at the entrance of slice $n$ as
$\psi_n$, so that none collides with a registered symbol.

---

## 1. From the Schrödinger equation to the paraxial wave equation

A beam electron of kinetic energy $eV$ has the relativistic wavelength and interaction parameter

$$
\lambda = \frac{h}{\sqrt{2 m_0 e V\left(1 + \dfrac{eV}{2 m_0 c^2}\right)}},
\qquad
\sigma_e = \frac{2\pi m e \lambda}{h^2} = \frac{2\pi}{\lambda V}\,\frac{m_0c^2 + eV}{2m_0c^2 + eV},
$$ (eq-ms-lambda-sigma)

with $m = \gamma m_0$ the relativistic mass. At 200 kV, $\lambda = 0.02508$ Å and
$\sigma_e = 7.288\times10^{-4}$ V⁻¹ Å⁻¹. In the crystal potential $V(\mathbf r, z)$ (a positive
electrostatic potential, in volts) the wave function satisfies the relativistically corrected
Schrödinger equation

$$
\left[\nabla^2 + 4\pi^2 k^2 + \frac{4\pi\sigma_e}{\lambda} V(\mathbf r, z)\right]\Psi = 0,
\qquad k = 1/\lambda .
$$ (eq-ms-schroedinger)

Writing $\Psi = \psi(\mathbf r, z)\,e^{2\pi i z/\lambda}$ separates the fast carrier from a slowly
varying envelope. Because the electron energy (hundreds of keV) dwarfs the potential (tens of eV),
$\psi$ changes little over one wavelength along $z$, and the term $\partial^2\psi/\partial z^2$ can
be dropped against $(4\pi i/\lambda)\,\partial\psi/\partial z$. This **paraxial** (high-energy,
forward-scattering) approximation leaves a first-order equation in $z$:

$$
\frac{\partial\psi}{\partial z}
= \left[\frac{i\lambda}{4\pi}\nabla_{\!\perp}^{2} + i\sigma_e V(\mathbf r, z)\right]\psi .
$$ (eq-ms-paraxial)

It has the form of a time-dependent Schrödinger equation with $z$ in the role of time. Back
scattering is excluded, which is sound for fast electrons but is the first stated limitation of
the method.

## 2. Operator splitting: why the method has two steps

Over a slice from $z_n$ to $z_n + \Delta z$ the formal solution of {eq}`eq-ms-paraxial` is

$$
\psi(z_n + \Delta z) = \exp\!\left[\frac{i\lambda\Delta z}{4\pi}\nabla_{\!\perp}^{2}
  + i\sigma_e v_n(\mathbf r)\right]\psi(z_n) ,
\qquad
v_n(\mathbf r) = \int_{z_n}^{z_n+\Delta z} V(\mathbf r, z)\,dz ,
$$ (eq-ms-formal)

where $v_n$ is the **projected potential** of the slice, in V Å. The exponent is the sum of two
operators that do not commute: the Laplacian is diagonal in Fourier space, the potential in real
space. Splitting the exponential (the Baker–Campbell–Hausdorff expansion) gives

$$
\psi_{n+1} = \mathcal{P}_{\Delta z} \ast \bigl[\tau_n\,\psi_n\bigr] + O(\Delta z^2),
\qquad
\tau_n(\mathbf r) = e^{i\sigma_e v_n(\mathbf r)},
$$ (eq-ms-split)

the product of a **transmission function** $\tau_n$ — a pure phase grating, exact for a slice of
zero thickness — and the **Fresnel propagator** $\mathcal P_{\Delta z}$, exact for vacuum. The
neglected term is proportional to the commutator $[\nabla_\perp^2, v_n]$: it grows with the slice
thickness and with the lateral gradient of the potential. Slices of 1–2 Å are therefore standard, and
a slice thickness that divides the crystal period along the beam keeps every slice identical.

### 2.1 The Fresnel propagator

$\nabla_\perp^2$ becomes $-4\pi^2 g^2$ in Fourier space, so free propagation over $\Delta z$ is a
multiplication of the spectrum:

$$
\mathcal P_{\Delta z}(\mathbf g) = \exp\!\left(-i\pi\lambda g^2\Delta z\right) .
$$ (eq-ms-propagator)

Two propagators compose, $\mathcal P_{\Delta z_1}\mathcal P_{\Delta z_2} = \mathcal P_{\Delta z_1 +
\Delta z_2}$, and $|\mathcal P| = 1$, so propagation conserves $\int|\psi|^2$. The factor
$e^{-i\pi\lambda g^2\Delta z}$ is exactly the phase a beam $\mathbf g$ accumulates from its
excitation error $s_g = -\lambda g^2/2$ at a zone axis: the Ewald-sphere curvature is built in.

### 2.2 Beam tilt

A small tilt $\boldsymbol\theta_{\mathrm{tilt}} = (\theta_x, \theta_y)$ of the incident beam shears
the propagation, and to first order adds a linear phase ramp to the propagator (Kirkland 2010,
sec. 6.9; the form abTEM uses):

$$
\mathcal P_{\Delta z}(\mathbf g) = \exp\!\left[-i\pi\lambda g^2\Delta z
  - 2\pi i\,\Delta z\,(g_x\tan\theta_x + g_y\tan\theta_y)\right].
$$ (eq-ms-tilt)

Tilt breaks the symmetry of the zone-axis pattern — Friedel pairs acquire different intensities —
and, in the image, shifts and blurs columns in proportion to the thickness.

## 3. The projected potential

### 3.1 From scattering factors to potential

In the independent-atom model the crystal potential is a superposition of isolated, spherical,
neutral-atom potentials. The electron scattering factor is the Fourier transform of one of them,

$$
f_e(\mathbf g) = \frac{2\pi m_0 e}{h^2}\int V_{\mathrm{atom}}(\mathbf r)\,
  e^{-2\pi i\,\mathbf g\cdot\mathbf r}\,d^3r ,
$$ (eq-ms-fe)

in Å. Inverting and integrating along $z$ over an atom that lies wholly inside the slice (the
**infinite projection**) sets $g_z = 0$, so the projected potential of the slice, periodic on a
cell of area $A = L_xL_y$, has the Fourier coefficients

$$
\tilde v_n(\mathbf g) = \frac{h^2}{2\pi m_0 e}\,\frac{1}{A}
  \sum_{j \in n} f_{e,j}(g)\,e^{-2\pi^{2}u_{\mathrm{th},j}^{2}g^{2}}\,
  e^{-2\pi i\,\mathbf g\cdot\mathbf r_j} ,
\qquad
\frac{h^2}{2\pi m_0 e} = 2\pi a_0 e = 47.878\ \text{V Å}^2 ,
$$ (eq-ms-potential)

the sum running over the atoms whose centres lie in slice $n$. The optional factor with the
one-axis RMS displacement $u_{\mathrm{th}}$ is the static Debye–Waller damping of a time-averaged
vibrating atom. PyTex evaluates {eq}`eq-ms-potential` on the Fourier grid with **exact** phase
factors: the structure factor of each species factorizes into an $x$ part and a $y$ part and is one
matrix product. abTEM instead places each atom by bilinear interpolation of a delta function and
corrects approximately with a sinc; the two agree to single precision when atoms sit on grid
points and to about a percent otherwise (tests below).

Two consequences are exact and are tested:

- The integral of $v_n$ over the cell is its zero-frequency coefficient times $A$:
  $\int v\,d^2r = 2\pi a_0 e\sum_j f_{e,j}(0)$.
- For a slab filling a box of height $L_z$, the **mean inner potential** is
  $V_0 = 2\pi a_0 e\sum_j f_{e,j}(0)/(A L_z)$ — about 14 V for silicon in this model, above the
  measured value near 12 V because bonding contracts the valence charge. That gap is a property of
  the independent-atom model, not of the propagation.

### 3.2 Parametrizations

`PotentialParametrization` selects the fit to $f_e$:

| Choice | Form of $f_e(g)$ | Source |
| --- | --- | --- |
| `lobato` (default) | $\sum_{i=1}^{5} a_i\,\dfrac{2 + b_i g^2}{(1 + b_i g^2)^2}$ | Lobato & Van Dyck (2014), fitted to Hartree–Fock densities with the correct large-$g$ asymptote |
| `kirkland` | $\sum_{i=1}^{3}\dfrac{a_i}{g^2 + b_i} + c_i\,e^{-d_i g^2}$ | Kirkland (2010), Appendix C |
| `mott_bethe` | $\dfrac{Z - f_x(s)}{8\pi^2 a_0 s^2}$, $s = g/2$ | PyTex's X-ray table ({doc}`dynamical_cbed_and_symmetry_determination`) |

The first two are generated into `electron_potential_parametrizations.json` from abTEM's own
tables, which transcribe the papers; the third is what the Bloch-wave solver uses, and so is the
one to choose when the two methods are compared. For silicon the three agree within 3 % below
$g = 2$ Å⁻¹.

The Kirkland form has a closed-form infinite projection in real space, obtained by transforming
each term in two dimensions — $\int e^{2\pi i\mathbf g\cdot\mathbf r}/(g^2 + b)\,d^2g =
2\pi K_0(2\pi r\sqrt b)$ and $\int e^{-dg^2}e^{2\pi i\mathbf g\cdot\mathbf r}\,d^2g =
(\pi/d)\,e^{-\pi^2r^2/d}$:

$$
v(r) = 4\pi^2 a_0 e\sum_i a_i K_0\!\left(2\pi r\sqrt{b_i}\right)
  + 2\pi^2 a_0 e\sum_i \frac{c_i}{d_i}\,e^{-\pi^2 r^2/d_i} .
$$ (eq-ms-kirkland-real)

The test suite compares the grid synthesis of {eq}`eq-ms-potential` with {eq}`eq-ms-kirkland-real`
for gold. The grid truncates $f_e$ at the Nyquist frequency, which rings about the closed form
(the Gibbs phenomenon of the logarithmic core singularity); at 0.01 Å sampling the ringing is
below 1 % of the potential 0.5 Å from the nucleus.

## 4. Sampling and the band limit

On an $N_x\times N_y$ grid with pixel $\Delta x$ the representable frequencies reach the Nyquist
frequency $g_N = 1/(2\Delta x)$. The transmission step multiplies two functions, and a product of
spectra limited to $g_{\max}$ has spectrum reaching $2g_{\max}$. On a periodic grid everything
beyond $g_N$ wraps back by $2g_N$, landing at $|g| \ge 2g_N - 2g_{\max}$. Those aliases stay out of
the retained band $|g| \le g_{\max}$ when $2g_N - 2g_{\max} \ge g_{\max}$, that is

$$
g_{\max} = \tfrac{2}{3}\,g_N = \frac{1}{3\,\Delta x},
\qquad
\alpha_{\max} = \lambda g_{\max} .
$$ (eq-ms-bandlimit)

PyTex, like abTEM, multiplies both $\tau_n$ and $\mathcal P_{\Delta z}$ by a circular aperture at
$g_{\max}$ with a raised-cosine edge of width $0.01/\Delta x$. Scattering beyond $\alpha_{\max}$ is
removed, so $\int|\psi|^2$ falls below 1 by exactly the intensity scattered out of the band: the
**retained intensity** that every run reports is therefore a direct check on the sampling. HRTEM
needs only the frequencies the objective aperture passes, but the *propagation* needs the
high-angle scattering too, because it feeds back into the low-order beams; a sampling of
0.05–0.1 Å ($\alpha_{\max}\approx 80$–170 mrad at 200 kV) is typical.

## 5. A periodic specimen

The grid is periodic, so the calculation is for an infinite lateral repetition of the box. For a
crystal the box must therefore be an exact lattice period in $x$ and $y$; otherwise every edge is a
seam, the Bragg beams fall between Fourier-grid points, and kinematically forbidden reflections
appear. `zone_axis_cell` constructs such a box for a zone axis $[uvw]$: $\mathbf t_3$ is the
primitive lattice vector along $[uvw]$, and $\mathbf t_1$, $\mathbf t_2$ are the shortest lattice
vectors perpendicular to it and to each other, found by searching integer vectors. For cubic
$[110]$ this is $[001]$, $[1\bar10]$, $[110]$; for hexagonal $[001]$ it is $[100]$, $[120]$,
$[001]$. `periodic_slab` fills $m_1\mathbf t_1\times m_2\mathbf t_2\times m_3\mathbf t_3$ with every
atom whose fractional box coordinates lie in $[0, 1)$ and checks the count against
$|\det[\mathbf t_1\mathbf t_2\mathbf t_3]|\,m_1m_2m_3$ times the atoms per cell. A beam $\mathbf g$
then sits exactly at grid index $(\mathbf g\cdot m_1\mathbf t_1,\ \mathbf g\cdot m_2\mathbf t_2)$.

A low-symmetry lattice at a general zone axis may have no perpendicular lattice pair within reach,
and PyTex then refuses rather than straining the cell, because a strained box simulates a
different crystal.

## 6. From the exit wave to the image

### 6.1 The coherent image

The objective lens multiplies the spectrum of the exit wave by the transfer function
$H(\mathbf g) = A(g)\,e^{-i\chi(\mathbf g)}$, with $A$ the aperture and $\chi$ the wave aberration
of {doc}`hrem_multislice_and_ctf`. The image intensity is

$$
I(\mathbf r) = \Bigl|\mathcal F^{-1}\bigl\{\tilde\psi_{\mathrm{exit}}(\mathbf g)\,H(\mathbf g)\bigr\}\Bigr|^2 .
$$ (eq-ms-image)

The defocus term of $\chi$ is $\pi\lambda\Delta f g^2$, so $e^{-i\chi}$ restricted to it is
exactly $\mathcal P_{\Delta f}$ of {eq}`eq-ms-propagator`: **imaging at defocus $\Delta f$ is
Fresnel propagation of the exit wave by $\Delta f$** (negative $\Delta f$, underfocus, propagates
back towards the specimen). This identity is exact, is a worked example, and is
the reason a focal series costs one multislice run.

### 6.2 Partial coherence

A real source is neither monochromatic nor a point. The **temporal** incoherence — a Gaussian
spread of defocus of standard deviation $\Delta$ from chromatic aberration and energy spread — and
the **spatial** incoherence — a cone of illumination directions of semi-angle $\alpha$ — make the
recorded image an incoherent average of coherent images.

**Quasi-coherent approximation.** Averaging $e^{-i\pi\lambda\delta g^2}$ over
$\delta\sim\mathcal N(0,\Delta^2)$ gives Frank's envelope

$$
E_c(g) = \exp\!\left(-\tfrac12\pi^2\lambda^2\Delta^2 g^4\right),
$$ (eq-ms-ec)

and a first-order expansion in the illumination angle gives
$E_s(g) = \exp[-\pi^2\alpha^2(\Delta f\,g + C_s\lambda^2 g^3)^2]$. Multiplying the wave transfer by
$E_cE_s$ is exact for the **linear** image terms — the interference of each scattered beam with
the unscattered one — and is what abTEM applies. It misstates the **non-linear** terms, the
interference between two scattered beams $\mathbf g$ and $\mathbf g'$, which are damped by
$\chi(\mathbf g) - \chi(\mathbf g')$ rather than by $\chi(\mathbf g)$ alone. For a weak object the
non-linear terms are negligible; for heavy columns they are not.

**Focal integration.** `TemporalCoherence.FOCAL_INTEGRATION` computes the average itself, as a
Gauss–Hermite quadrature over the defocus distribution:

$$
I(\mathbf r) = \sum_k \frac{w_k}{\sqrt\pi}\,
  \Bigl|\mathcal F^{-1}\{\tilde\psi_{\mathrm{exit}}\,H_{\Delta f + \sqrt2\,\Delta\,x_k}\}\Bigr|^2 ,
$$ (eq-ms-focal-integration)

with nodes $x_k$ and weights $w_k$. It is exact for every term, provided the quadrature resolves
the phases it averages. The defocus phase of the interference between beams $\mathbf g$ and
$\mathbf g'$ varies across the spread as $e^{iax}$ with $a = \sqrt2\,\pi\lambda\Delta\,|g^2 - g'^2|$,
and the exact average of that is $e^{-a^2/4}$ — for $\mathbf g' = 0$ precisely $E_c(g)$. An
$n$-node Gauss–Hermite rule integrates $e^{iax}$ to better than $10^{-8}$ only for
$a \le 2\sqrt{n-16}$ (checked numerically to 512 nodes); beyond that it returns spurious transfer
where the true transfer is nil. PyTex therefore chooses $n$ from the largest $a$ the calculation
needs — set by the highest frequency the objective aperture and the grid pass — up to 512 nodes
($a \le 44.5$), and excludes only frequencies beyond what that resolves, where $E_c < e^{-495}$ and
nothing but the interference of two such beams is lost.

The difference between the two treatments is physical, and it is large exactly where it matters.
Two *equivalent* beams, $|\mathbf g| = |\mathbf g'|$, interfere with $a = 0$: the focal spread does
not damp their interference at all, whereas the envelope damps it by $E_c(g)^2$. This is the
non-linear image contribution that lets a crystal image show spacings finer than the linear
information limit. The tests assert agreement within 1 % of the image range for a hydrogen atom
(a weak object) and a difference above 5 % for a gold atom.

### 6.3 Frozen phonons

Thermal vibration does two things: it smears the time-averaged potential (the static Debye–Waller
factor of {eq}`eq-ms-potential`) and it scatters electrons *incoherently* between the Bragg beams —
thermal diffuse scattering (TDS). The frozen-phonon model (Loane, Xu & Silcox 1991) captures both:
each of $N_{\mathrm{fp}}$ configurations displaces every atom by an independent Gaussian of RMS
$u_{\mathrm{th}}$ per axis (an Einstein model, $u_{\mathrm{th}}^2 = B/8\pi^2$), a full multislice is
run per configuration, and **intensities** — not waves — are averaged:

$$
I(\mathbf r) = \frac{1}{N_{\mathrm{fp}}}\sum_{c=1}^{N_{\mathrm{fp}}}
  \Bigl|\mathcal F^{-1}\{\tilde\psi^{(c)}_{\mathrm{exit}}H\}\Bigr|^2 .
$$ (eq-ms-frozen-phonon)

The configuration-averaged *wave* is the elastic (coherent) wave; the rest of the averaged
intensity is the TDS, which appears in the diffraction pattern as a diffuse background between the
Bragg beams (a test places a probe midway between beams and finds it empty for a static lattice and
populated with phonons). A few to a few tens of configurations are typical.

## 7. Series: thickness, focus and the defocus–thickness map

**Thickness series.** The wave at depth $z$ is the wave after the slices above $z$, so the exit
wave at any set of thicknesses comes from one pass (`exit_depths_angstrom`). Along a zone axis the
intensity oscillates between the columns and the background with depth — **channelling**, the
real-space counterpart of the Pendellösung oscillation of the Bragg beams — so HRTEM contrast
reverses with thickness as well as with focus.

**Focal series.** The exit wave does not depend on the lens, and by §6.1 changing the defocus is
free propagation. A through-focus series $I(\mathbf r; \Delta f_j)$ therefore costs one multislice
and one pair of FFTs per image (`MultisliceExitWave.focal_series`). Experimentally, a focal series
is the input of exit-wave reconstruction, which inverts exactly this forward model to recover
$\psi_{\mathrm{exit}}$ free of the lens (Coene et al. 1992; Thust et al. 1996); simulated, it shows
how column contrast reverses at the Scherzer and Lichte conditions and how the weakest contrast
marks Gaussian focus for a weak object.

**Defocus–thickness map.** Neither the thickness nor the defocus of an experimental image is known
a priori. The standard practice (O'Keefe & Kilaas 1988) is to simulate a tableau of images over a
grid of thicknesses (rows) and defoci (columns) and to match the experimental image against it;
`MultisliceExitWave.defocus_thickness_map` produces the tableau from one run.

(multislice-validation)=

## Validation

Every expected value below comes from a closed form or from a different method.

| Check | Expected | Where |
| --- | --- | --- |
| $\int v\,d^2r = 2\pi a_0e\,f_e(0)$ | exact | worked example, `test_multislice.py` |
| Kirkland real-space closed form {eq}`eq-ms-kirkland-real` | < 1 % (Gibbs) | `test_multislice.py` |
| Vacuum leaves the plane wave unchanged; propagators compose | exact | `test_multislice.py` |
| One slice is the band-limited phase object | exact | `test_multislice.py` |
| Thickness series in one pass equals separate thinner runs | exact | `test_multislice.py` |
| Multislice equals Bloch waves on one potential (Si [001] ZOLZ, 54–201 Å, 8 slices per period) | $<1.5\times10^{-3}$ | worked example, `test_multislice.py` |
| Forbidden (200) of diamond: dark in the ZOLZ, weak double diffraction via HOLZ | $<10^{-15}$ / $10^{-9}$–$10^{-3}I_{220}$ | `test_multislice.py` |
| Tilt mirrors the pattern of a mirror-symmetric crystal | exact | `test_multislice.py` |
| Defocus equals Fresnel propagation | exact | worked example |
| Focal integration equals Frank's envelope for a linear fringe | exact | worked example |
| Tables, wavelength, $\sigma_e$ identical to abTEM | exact | `test_multislice_abtem_parity.py` |
| Exit wave vs abTEM, atoms on grid / off grid | $<2\times10^{-3}$ / $<2\times10^{-2}$ | `test_multislice_abtem_parity.py` |
| Image through the lens vs abTEM | $<0.5$ % of the image range | `test_multislice_abtem_parity.py` |

The abTEM tests run only where abTEM is installed; everything else runs in every lane.

## Limitations

- **Independent-atom model.** Neutral spherical atoms; no bonding or ionicity, so the mean inner
  potential is overestimated by 10–20 %.
- **Elastic and paraxial.** No inelastic (plasmon, core-loss) scattering, no absorption potential,
  no back-scattering.
- **Thermal scattering** only through frozen phonons, in an Einstein model with uncorrelated
  displacements.
- **Coherence.** Spatial coherence is always quasi-coherent; temporal coherence is quasi-coherent
  or integrated exactly up to the frequency 512 Gauss–Hermite nodes resolve.
- **Detector.** No modulation transfer function and no shot noise are applied.
- **Box.** Periodic, orthogonal and lateral-lattice-periodic for crystals. A structure imported from
  a file must supply such a box itself.

## References

1. Cowley, J. M. & Moodie, A. F. (1957). The scattering of electrons by atoms and crystals. I. A
   new theoretical approach. *Acta Cryst.* **10**, 609–619. doi:10.1107/S0365110X57002194
2. Ishizuka, K. & Uyeda, N. (1977). A new theoretical and practical approach to the multislice
   method. *Acta Cryst.* **A33**, 740–749. doi:10.1107/S0567739477001879
3. Frank, J. (1973). The envelope of electron microscopic transfer functions for partially coherent
   illumination. *Optik* **38**, 519–536.
4. O'Keefe, M. A. & Kilaas, R. (1988). Advances in high-resolution image simulation. *Scanning
   Microscopy Suppl.* **2**, 225–244.
5. Loane, R. F., Xu, P. & Silcox, J. (1991). Thermal vibrations in convergent-beam electron
   diffraction. *Acta Cryst.* **A47**, 267–278. doi:10.1107/S0108767391000375
6. Coene, W., Janssen, G., Op de Beeck, M. & Van Dyck, D. (1992). Phase retrieval through focus
   variation for ultra-resolution in field-emission transmission electron microscopy. *Phys. Rev.
   Lett.* **69**, 3743–3746. doi:10.1103/PhysRevLett.69.3743
7. Thust, A., Coene, W. M. J., Op de Beeck, M. & Van Dyck, D. (1996). Focal-series reconstruction
   in HRTEM: simulation studies on non-periodic objects. *Ultramicroscopy* **64**, 211–230.
   doi:10.1016/0304-3991(96)00011-2
8. Kirkland, E. J. (2010). *Advanced Computing in Electron Microscopy*, 2nd ed. Springer.
   doi:10.1007/978-1-4419-6533-2
9. Lobato, I. & Van Dyck, D. (2014). An accurate parameterization for scattering factors, electron
   densities and electrostatic potentials for neutral atoms that obey all physical constraints.
   *Acta Cryst.* **A70**, 636–649. doi:10.1107/S205327331401643X
10. Madsen, J. & Susi, T. (2021). The abTEM code: transmission electron microscopy from first
    principles. *Open Research Europe* **1**, 24. doi:10.12688/openreseurope.13015.2
