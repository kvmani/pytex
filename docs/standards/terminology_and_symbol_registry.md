# Terminology And Symbol Registry

This document fixes the repository-wide vocabulary and symbol policy for PyTex. It
is the single source of truth for nomenclature: the same term and the same symbol
keep the same meaning across every documentation surface and in code explanations.

## Purpose

PyTex is large enough that symbol drift and term drift become real scientific
risks. This registry exists so the same term and the same symbol keep the same
meaning across:

- Sphinx concept pages
- workflow guides
- notebook tutorials
- executable worked examples
- LaTeX theory and algorithm notes
- code explanations and docstrings where symbols are discussed
- canonical SVG figures

## Policy

- Stable scientific terms should be defined once here and reused elsewhere.
- Stable mathematical symbols should be introduced here or in the closest
  canonical theory note, then reused consistently.
- If a page needs a local symbol extension, it should state that extension
  explicitly and keep the core registry symbols unchanged.
- Pages that rely on registry terms should link back to the user-facing glossary
  page and, when needed, to this standards document.
- When a canonical SVG figure labels one of these terms or symbols, it should use
  the same wording and symbol form as this registry unless the figure explicitly
  documents a local teaching simplification.
- Executable worked examples must reference symbols from this registry; a symbol
  used in a worked example that is not yet registered must be added here first.
  See [Executable Worked Examples](executable_examples.md).

## Symbol Registration Policy

New notation is introduced through this registry, not ad hoc in a single page:

1. Prefer an existing registered symbol.
2. If a new symbol is genuinely needed, add it to the appropriate table below with
   a one-line fixed meaning before using it in prose, math, figures, notebooks,
   worked examples, or docstrings.
3. Do not reuse a registered symbol for a different meaning; choose a distinct
   symbol or an explicitly scoped local extension.
4. Keep storage-order and normalization conventions (quaternion order, reciprocal
   normalization, Euler labeling) identical to those fixed in
   [Notation and Conventions](notation_and_conventions.md).

## Core Terms

| Term | Fixed meaning |
| --- | --- |
| reference frame | A named, domain-typed coordinate frame such as crystal, specimen, map, detector, laboratory, or reciprocal. |
| canonical Cartesian reference | The right-handed $(X, Y, Z)$ frame in which every `ReferenceFrame`'s axis vectors are expressed. |
| axis vectors | Components of a frame's three labelled axes in the canonical Cartesian reference; dimensionless orientation only, not lattice edge lengths. |
| sample frame | The rolling-geometry specimen frame with axes $\mathrm{RD}$ (rolling direction), $\mathrm{TD}$ (transverse direction), $\mathrm{ND}$ (normal direction). |
| frame transform | A typed rigid map between exactly two named frames, converting source-frame components into target-frame components. |
| frame graph | A registry of frames and declared transforms that resolves any connected pair by composing the shortest declared chain. |
| axis correspondence | A declaration, in axis labels, of which target axis each source axis lands on; equivalent to a signed permutation matrix. |
| frame triad | The drawable three-arrow representation of a reference frame, shared by 3D scenes, corner gizmos, and documentation SVG. |
| orientation | A crystal-to-specimen mapping carried by an explicit `Orientation` object. |
| rotation | A geometric active rotation that does not by itself define crystallographic source and target meaning. |
| misorientation | The orientation mapping between two orientations, before symmetry reduction. |
| disorientation | The symmetry-reduced misorientation of minimal rotation angle in the fundamental zone. |
| symmetry | Point-group-facing operator set used for orientation and direction reduction. |
| space group | Structure-facing crystallographic identity used for phases and CIF-backed construction. |
| pole figure | Distribution of crystal directions or plane normals expressed relative to specimen directions. |
| multiples of a random distribution (m.r.d.) | The scale on which a pole density or ODF value is 1 where the distribution is random. Defined by the solid-angle-weighted mean over the sampled region being exactly 1, so it is a property of an integral and not of a maximum or a sum. |
| pole-figure sampling | Whether a pole figure's intensities are per-pole weights of a cloud of poles (**scattered poles**) or densities already evaluated at the given directions (**sampled density**). The two require different resampling estimators, so the reading is recorded rather than inferred. |
| random-standard defocusing calibration | A reflection-specific radial intensity-loss curve measured from an untextured reference specimen under the specimen scan's instrument conditions. Background is subtracted before azimuthal ring reduction and normalization to the lowest calibrated tilt; it is not a texture model. |
| pole-figure difference | The signed field $\Delta P$ obtained by subtracting one pole figure from another on a shared support. Signed, and therefore not itself a pole figure. |
| residual pole figure | A pole-figure difference between a measurement and the figure an ODF recalculates for it; the spatial form of a goodness-of-fit check. |
| solid-angle weights | Per-point integration weights making a sum over a sampled support approximate an integral over the sphere. Required whenever a mean is taken over a tilt/rotation raster, which over-samples the specimen normal. |
| inverse pole figure | Distribution of specimen directions expressed in crystal coordinates and reduced by symmetry where appropriate. |
| ODF | Orientation distribution function over orientation space. |
| fundamental zone | The symmetry-reduced subset of orientation space used for canonical orientation keys. |
| zone axis | Direct-space crystallographic direction defining an electron-diffraction viewing or incidence condition. |
| zone-axis order | The number of distinct Kikuchi band centre lines crossing a zone axis. The operational measure of how conspicuous the axis is on the screen: a four-band crossing is unmistakable, a two-band one is a guess. Distinct from the axis's indices, and the quantity a Kikuchi map sorts by. |
| Kossel cone | Either of the two cones of semi-angle $90^{\circ} - \theta_{B}$ about a plane normal, on which incoherently scattered electrons satisfy that plane's Bragg condition. Their traces are the two edges of a Kikuchi band. |
| Kikuchi map | The band-and-zone-axis network of a phase drawn on a stereographic projection of the crystal sphere, used to plan a tilt route between zone axes. Distinct from a Kikuchi *pattern*, which is one orientation on one detector in gnomonic coordinates. |
| band-followable route | A tilt path between zone axes in which every leg lies on a single Kikuchi band, so the operator has a line to track. Coincides with the geodesic, because both are defined by the plane the two axes span. |
| holder frame | The specimen-domain frame for TEM work: rigidly attached to the specimen and the holder cartridge, coinciding with the laboratory frame at zero tilt. Not a new frame domain. |
| crystal-to-holder orientation | The rotation $\mathbf{U}$ taking crystal-frame vectors into holder-frame components; absorbs the separately unobservable specimen-mounting rotation. |
| holder tilt | A stage rotation about one of the two holder axes, $\alpha$ about the rod and $\beta$ about the cradle carried inside it. |
| tilt envelope | The set of stage positions a holder can physically reach, as a predicate with a continuous margin rather than a pair of limits. |
| diffraction rotation | The azimuth of the recorded diffraction pattern relative to the laboratory frame; instrument-specific, hysteretic in the lens settings, and not present in file metadata. |
| observation stabilizer | The rotations of the Laue class that map an observed zone plane to itself; the group up to which a single kinematic pattern determines an orientation. |
| ambiguity family | A competing hypothesis about the true orientation, of which only one is correct; distinct from a symmetry equivalent, which is a free choice. |
| powder pattern | Grid-sampled XRD spectrum built from discrete reflections and an optional broadening model. |
| multiplicity | Number of symmetry-equivalent members of a plane or direction family under the phase point group. |
| crystal scene | Reusable geometry bundle for 3D crystal rendering. |
| reciprocal star | The $*$ marking a reciprocal-space **basis vector** or reciprocal-frame axis label ($a^{*}$); never applied to Miller indices, which are already reciprocal-basis components. |
| symmetry family | The symmetry-related orbit of a plane or direction, written $\{hkl\}$ or $\langle uvw \rangle$; a single member is written $(hkl)$ or $[uvw]$. |
| overbar notation | A negative index rendered with a bar over the digit rather than a leading minus, in publication-facing output. |
| zone law | The condition $hu + kv + lw = 0$ for a direction $[uvw]$ lying in a plane $(hkl)$. |
| Bloch wave | An eigenvector of the dynamical structure matrix: a wavefield that propagates through a perfect crystal unchanged in shape, attenuated at its own rate. The many-beam solution is a superposition of them. |
| projection approximation | A dynamical calculation confined to the zeroth Laue zone. It samples the potential *projected* along the beam, whose symmetry is at least as high as the crystal's and is often strictly higher; every CBED symmetry conclusion drawn under it is unsafe. |
| anomalous absorption | The unequal absorption of the Bloch-wave branches, from the off-diagonal imaginary potential. It makes the bright-field rocking curve asymmetric in $s_{g}$ while leaving the dark-field one symmetric. |

## Core Symbols

### Frames, vectors, and rotations

| Symbol | Meaning |
| --- | --- |
| $\mathbf{v}$ | Generic vector in an explicitly named frame. |
| $\hat{\mathbf{v}}$ | Unit vector (normalized) in an explicitly named frame. |
| $\mathbf{R}$ | Rotation matrix acting actively on vectors. |
| $\mathbf{T}$ | Rigid placement of geometry into a shared world frame: $\mathbf{T}(\mathbf{x}) = \mathbf{R}\,\mathbf{x} + \mathbf{t}$ (the `Transform3D` visualization primitive). |
| $q$ | Unit quaternion in `w, x, y, z` storage order. |
| $(\phi_1, \Phi, \phi_2)$ | Bunge Euler angles. |
| $(\mathbf{n}, \omega)$ | Axis-angle pair: rotation axis $\mathbf{n}$ and angle $\omega$. |
| $\boldsymbol{\rho}$ | Rodrigues vector $\hat{\mathbf{n}}\tan(\omega/2)$; the chart in which symmetry fundamental zones are convex polyhedra. |
| $\boldsymbol{\rho}_{F}$ | Rodrigues-Frank homogeneous form $(\hat{\mathbf{n}}, \tan(\omega/2))$, whose magnitude is projective so $\omega = \pi$ stays representable. |
| $\mathbf{h}$ | Homochoric vector $\hat{\mathbf{n}}[\tfrac{3}{4}(\omega - \sin\omega)]^{1/3}$: the equal-volume chart of SO(3), a ball of radius $R_{1}$. |
| $\mathbf{c}$ | Cubochoric coordinate: the equal-volume chart mapped onto a cube of edge $a_{p} = \pi^{2/3}$, so a uniform Cartesian grid is a uniform grid of orientations. |
| $R_{1}$ | Radius of the homochoric ball, $(3\pi/4)^{1/3} \approx 1.3307$. |
| $a_{p}$ | Edge of the cubochoric cube, $\pi^{2/3} \approx 2.1450$; the cube and the ball both enclose the volume $\pi^{2}$ of SO(3). |
| $\mathbf{B}$ | Frame basis matrix whose columns are a frame's axis vectors in the canonical Cartesian reference: $\mathbf{x} = \mathbf{B}\,\mathbf{v}$. |
| $\mathrm{RD}, \mathrm{TD}, \mathrm{ND}$ | Rolling, transverse, and normal directions of the sample frame. |
| $\{hkl\}$ | Symmetry-related family of lattice planes. |
| $\langle uvw \rangle$ | Symmetry-related family of lattice directions. |

### Lattice and reciprocal lattice

| Symbol | Meaning |
| --- | --- |
| $\mathbf{a}, \mathbf{b}, \mathbf{c}$ | Direct-lattice basis vectors. |
| $a, b, c, \alpha, \beta, \gamma$ | Lattice parameters (edge lengths and angles). |
| $\mathbf{a}^{*}, \mathbf{b}^{*}, \mathbf{c}^{*}$ | Reciprocal-lattice basis vectors under the PyTex normalization rule. |
| $\mathbf{g}_{hkl}$ | Reciprocal-lattice vector associated with Miller indices $(hkl)$. |
| $\mathbf{G}$ | Direct-space metric tensor. |
| $d_{hkl}$ | Interplanar spacing for the $(hkl)$ family. |

### Miller indices and crystallographic geometry

| Symbol | Meaning |
| --- | --- |
| $(hkl)$ | Miller plane indices; $(hkil)$ for the four-index hexagonal form. |
| $[uvw]$ | Miller direction indices; $[uvtw]$ for the four-index hexagonal form. |
| $\mathbf{n}$ | Plane normal direction. |
| $\angle(\mathbf{n}_1, \mathbf{n}_2)$ | Angle between two plane normals (interplanar angle). |
| $\angle(\mathbf{d}_1, \mathbf{d}_2)$ | Angle between two lattice directions. |
| $m_{\{hkl\}}$ | Symmetry multiplicity of a plane family under the phase point group. |
| $\mathbf{M}$ | Direction-index correspondence matrix of an orientation relationship: $\mathbf{u}_{c} = \mathbf{M}\,\mathbf{u}_{p}$ with $\mathbf{M} = \mathbf{A}_{c}^{-1}\mathbf{R}\,\mathbf{A}_{p}$ built from the direct structure matrices. |
| $\mathbf{M}^{*}$ | Plane-index correspondence matrix of an orientation relationship: $\mathbf{h}_{c} = \mathbf{M}^{*}\,\mathbf{h}_{p}$, with $\mathbf{M}^{*} = \mathbf{M}^{-\mathsf{T}}$ so the zone law $\mathbf{h} \cdot \mathbf{u}$ is preserved. |
| $G_p$, $G_c$ | Parent and child crystal point groups, as sets of proper rotation operators. |
| $G_c \left(R\,G_p\,R^{\mathsf{T}}\right) G_c$ | Same-parent boundary fingerprint: the set of misorientations two child grains of one parent can exhibit, since $\mathbf{C}_i^{\mathsf{T}}\mathbf{C}_j = \mathbf{V}_i \mathbf{V}_j^{\mathsf{T}}$ and $\mathbf{V}_i = R\,S_{p,i}$. Each child orientation is defined only up to its own crystal symmetry, hence the double coset. |
| $\mathbf{V}_i = \mathbf{C}_i^{\mathsf{T}} \mathbf{P}_i$ | Measured parent-to-child rotation of one orientation pair, in the canonical crystal-to-specimen convention $\mathbf{C} = \mathbf{P}\,\mathbf{V}^{\mathsf{T}}$. The quantity every OR fitting and determination surface averages. |
| $G_c \mathbf{V} G_p$ | Double coset of one measured pair: all symmetry-equivalent descriptions of its parent-to-child rotation. Its maximum-trace (minimum-angle) element is the pair's disorientation description of the relationship, and is the seed for determining an OR with no nominal supplied. Not unique when the relationship's own rotation is symmetric. |

### Orientation and misorientation

| Symbol | Meaning |
| --- | --- |
| $g$ | An orientation (crystal-to-specimen mapping). |
| $\Delta g$ | A misorientation between two orientations. |
| $\omega$ | Disorientation angle: minimal misorientation angle over the symmetry group. |
| $\Sigma$ | Coincidence-site-lattice index of a boundary (for example $\Sigma 3$). |

### Pole figures and texture intensity

| Symbol | Meaning |
| --- | --- |
| $P_{hkl}(\mathbf{y})$ | Pole density of the plane family $\{hkl\}$ at specimen direction $\mathbf{y}$, in multiples of a random distribution. |
| $\Delta P$ | Signed difference of two pole densities on a shared support. |
| $w_i$ | Solid-angle integration weight of sampled direction $i$; weights sum to 1 over the sampled region. |
| $\psi$ | Polar (tilt) angle of a specimen direction from the specimen-frame $+Z$ axis. |
| $d_i$ | Positive random-standard defocusing factor at specimen tilt $i$, normalized to 1 at the lowest calibrated tilt. A specimen intensity is corrected as $(I_i-b)/d_i$, never $I_i/d_i-b$. |
| $\ell$ | Degree of a generalized spherical harmonic term. Diffraction pole figures determine only **even** $\ell$, which is the ghost problem and removes close to half the ODF basis. |
| $P_{\mathrm{rand}}$ | Kernel-density response a random texture produces. Dividing a raw KDE response by it is what puts a pole figure on the m.r.d. scale; it depends on the kernel halfwidth, so un-normalized figures from different kernels are not comparable. |
| $\mathcal{T}$ | Fundamental sector: the spherical region holding one representative of each direction orbit. The standard stereographic triangle, generalized to any point group. |
| $\mathbf{K}$ | Colour-key basis: the $3\times3$ matrix whose **columns** are the fundamental sector's corner directions. |
| $\boldsymbol{\beta}$ | Barycentric weights of a direction in the sector-corner basis, $\mathbf{K}\boldsymbol{\beta} = \hat{\mathbf{d}}_{\mathcal{T}}$. Non-negative exactly inside the sector, so the membership test and the colour are the same inequalities. Distinct from the solid-angle weight $w_i$. |
| $\mathbf{C}_{\mathrm{rgb}}$ | Corner colour matrix whose rows are the RGB triples assigned to the sector corners; the identity gives the standard red/green/blue key. |
| $\gamma_{s}$ | IPF saturation parameter. Channels are raised to $1/\gamma_{s}$, so the default $\gamma_{s} = 0.5$ squares them. A contrast control with no crystallographic content; distinct from the lattice angle $\gamma$ and the relativistic factor. |

### The Kearns parameter

| Symbol | Meaning |
| --- | --- |
| $f$, $f_{\mathrm{RD}}$, $f_{\mathrm{TD}}$, $f_{\mathrm{ND}}$ | Kearns orientation parameter along a specimen direction: the volume-weighted mean of $\cos^{2}$ of the angle between each crystal's basal pole and that direction, and hence the effective fraction of basal poles aligned with it. Bounded to $[0,1]$; exactly $1/3$ for a random texture. Distinct from the ODF density $f(g)$, which is a function on $SO(3)$. |
| $\mathbf{A}$ | Pole orientation tensor $\langle \mathbf{c}\,\mathbf{c}^{\mathsf{T}}angle$ of the basal-pole distribution in the specimen frame, so that $f(\mathbf{d}) = \mathbf{d}^{\mathsf{T}}\mathbf{A}\mathbf{d}$. Unit trace, which is why the Kearns parameters of an orthonormal triad sum identically to 1. The weighted counterpart of the directional-statistics orientation tensor $oldsymbol{\Theta}$. |
| $\phi$ | Tilt of a crystal's basal pole $[0001]$ from the specimen reference direction. Distinct from the first Euler angle $arphi_1$ and from the pole-figure azimuth. |
| $I(\phi)$ | Basal-pole density averaged over the full $360^{\circ}$ of azimuth about the reference direction — the *tilt profile*. Kearns' Eq. (5) integrates it against $\sin\phi\cos^{2}\phi$. Scale-free: only ratios of $I$ enter $f$. |
| $V_{\Delta\phi}$ | Volume fraction of crystals whose basal pole lies in the tilt band $\Delta\phi$, equal to $I(\phi)\sin\phi\,\Delta\phi$ normalized over $[0,\pi/2]$. Vanishes at $\phi = 0$ however intense the pole is there, because the band has no area. |
| $\phi_{hkil}$ | Fixed angle between the $(hkil)$ plane normal and $[0001]$, from the phase's reciprocal metric. What lets a $	heta$-$2	heta$ peak intensity be read as a basal-pole density at a known tilt. |
| $ho$, $eta$ | ODF kernel shrinkage: $ho = \langle\cos^{2}etaangle$ of a pole smeared by the kernel, and $eta = (3ho-1)/2$ the factor by which convolution scales every departure of $\mathbf{A}$ from $\mathbf{I}/3$. Here $eta$ is a scalar shrinkage factor, distinct from the lattice angle $eta$ and the holder tilt $eta$. |

### Directional statistics

| Symbol | Meaning |
| --- | --- |
| $oldsymbol{\Theta}$ | Orientation tensor $	frac{1}{n}\sum_i \hat{\mathbf{v}}_i\hat{\mathbf{v}}_i^{\mathsf{T}}$ of a direction set. Invariant under $\hat{\mathbf{v}} \mapsto -\hat{\mathbf{v}}$, so it is the correct summary for axial data where the resultant cancels. Unit trace. Distinct from the rigid placement $\mathbf{T}$. |
| $\lambda_1 \le \lambda_2 \le \lambda_3$ | Eigenvalues of $oldsymbol{\Theta}$; non-negative and summing to 1. $(	frac13,	frac13,	frac13)$ uniform, $(0,	frac12,	frac12)$ girdle, $(0,0,1)$ cluster. A mean axis is identified only when $\lambda_3 > \lambda_2$. |

### Elastic properties

| Symbol | Meaning |
| --- | --- |
| $C_{ijkl}$, $S_{ijkl}$ | Rank-four elastic stiffness and compliance tensors in crystal-frame Cartesian coordinates. The primary representation; the $6\times6$ Voigt matrices are views of them. |
| $C_{mn}$, $S_{mn}$ | Voigt-compressed stiffness and compliance. Stiffness carries no factors; **compliance carries a factor 2 per shear index**, so $S_{mn} = f_m f_n S_{ijkl}$ with $f_m = 2$ for $m > 3$. Dropping the factors silently corrupts every shear modulus. |
| $E(\hat{\mathbf{n}})$ | Young's modulus along crystal direction $\hat{\mathbf{n}}$, $1/E = n_i n_j n_k n_l S_{ijkl}$. |
| $J$ | Cubic orientation factor $n_1^2n_2^2 + n_2^2n_3^2 + n_3^2n_1^2$; the only way direction enters cubic $E(\hat{\mathbf{n}})$. Runs $0$ along $\langle 100\rangle$ to $1/3$ along $\langle 111\rangle$. |
| $A_{Z}$ | Zener anisotropy ratio $2C_{44}/(C_{11}-C_{12})$; exactly 1 for an elastically isotropic cubic crystal. |
| $K_{V}, K_{R}$ | Voigt and Reuss aggregate bulk moduli. **Identical for cubic symmetry**, so a cubic aggregate's bulk modulus is exact rather than bounded. |
| $\mu_{V}, \mu_{R}, \mu_{H}$ | Voigt, Reuss and Hill aggregate shear moduli. Distinct from the invariant measure $\mathrm{d}\mu$ on $SO(3)$. $\mu_{H}$ is an average of bounds and is not itself a bound. |

### Slip and plasticity

| Symbol | Meaning |
| --- | --- |
| $m$ | Schmid factor $\cos\phi\cos\lambda$ of a slip system under uniaxial stress. Bounded above by exactly $1/2$, attained at $\phi = \lambda = 45^{\circ}$. |
| $\mathbf{N}^{(s)}$ | Symmetric Schmid tensor of slip system $s$, $	frac{1}{2}(\hat{\mathbf{d}}\hat{\mathbf{n}}^{\mathsf{T}} + \hat{\mathbf{n}}\hat{\mathbf{d}}^{\mathsf{T}})$, expressed in the sample frame. |
| $\gamma^{(s)}$ | Slip amount on system $s$; non-negative, with the two shear senses carried as separate systems. Distinct from the IPF saturation $\gamma_{s}$. |
| $\Gamma$ | Total slip $\sum_s \gamma^{(s)}$ minimised by the full-constraint Taylor problem. |
| $M$ | Taylor factor $\Gamma / arepsilon_{\mathrm{eq}}$. Unique even when the slip combination attaining it is not. |
| $arepsilon_{\mathrm{eq}}$ | Von Mises equivalent strain $\sqrt{	frac{2}{3}oldsymbol{arepsilon}:oldsymbol{arepsilon}}$. |

### Diffraction

| Symbol | Meaning |
| --- | --- |
| $\lambda$ | Radiation wavelength. |
| $\theta$ | Bragg half-angle. |
| $2\theta$ | Powder-diffraction scattering angle reported in XRD plots. |
| $F_{hkl}$ | Reflection structure-factor quantity or current PyTex proxy where explicitly stated. |
| $R_{p}$ | Unweighted powder-profile agreement factor, $\sum_i |I_{\mathrm{obs},i}-I_{\mathrm{calc},i}| / \sum_i I_{\mathrm{obs},i}$, following the IUCr pdCIF definition. |
| $R_{wp}$ | Weighted powder-profile agreement factor, $[\sum_i w_i(I_{\mathrm{obs},i}-I_{\mathrm{calc},i})^2 / \sum_i w_i I_{\mathrm{obs},i}^2]^{1/2}$, following the IUCr pdCIF definition. Distinct from an expected R factor or a full Rietveld goodness of fit. |
| $\hat{\mathbf{z}}$ | Unit zone-axis direction in direct space. |
| $u, v$ | Detector-plane plotting coordinates in SAED or detector geometry contexts. |
| $s_{g}$ | Excitation error of reflection $\mathbf{g}$: its deviation from the exact Bragg condition, measured along the zone axis in reciprocal angstrom. |
| $\xi_{g}$ | Two-beam extinction distance of reflection $\mathbf{g}$: the depth period of the intensity exchange between the transmitted and diffracted beams, $\pi V_{c}\cos\theta_{B}/(\lambda|F_{g}|)$. |
| $f_{e}(s)$ | Electron atomic scattering factor in angstrom, from the X-ray form factor by Mott-Bethe: $(Z - f_{x})/(8\pi^{2}a_{0}s^{2})$. |
| $\alpha$ | Convergence semi-angle of the illumination cone in CBED; sets the disc radius $\alpha/\lambda$ and the Kossel-Moellenstedt threshold. |
| $t$ | Foil thickness along the beam. |
| $H$ | Reciprocal-lattice layer spacing along a zone axis, $1/|\mathbf{r}_{uvw}|$; measured by the HOLZ ring radii $G_{n} \simeq \sqrt{2nH/\lambda}$. |
| $\theta_{B}$ | Bragg angle of a specific reflection, $\arcsin(\lambda/2d)$. Written with the subscript wherever it must not be confused with the generic scattering half-angle $\theta$. |
| $2\theta_{B}$ | Angular width of a Kikuchi band, exactly twice the Bragg angle and independent of detector geometry. Since $2\theta_{B} = 2\arcsin(\lambda/2d) \approx \lambda/d$, width is a *decreasing* function of spacing: the widest bands come from high-index planes and the strongest from low-index ones. |
| $\psi$ | Polar angle from the centre of a projection. Azimuthal projection radii: $\tan\psi$ gnomonic, $\tan(\psi/2)$ stereographic, $2\sin(\psi/2)$ equal-area. |
| $\delta\varphi$ | Uncertainty in the diffraction rotation: the azimuthal calibration error between the recorded pattern and the stage axes. A tilt of $\theta$ planned under this uncertainty misses by $2\arcsin(\sin(\delta\varphi/2)\sin\theta) \approx \delta\varphi\sin\theta$. |
| $\nu_{g}$ | Complex Fourier coefficient of the scaled lattice potential, $\lambda F_{g}/(\pi V_{c}\cos\theta_{g})$; the off-diagonal element of the dynamical structure matrix, with $|\nu_{g}| = 1/\xi_{g}$. Its **phase** is what distinguishes a centrosymmetric structure from a non-centrosymmetric one. |
| $\xi'_{g}$ | Absorption distance of reflection $\mathbf{g}$: the imaginary partner of $\xi_{g}$ in the complex optical potential. $\xi'_{0}$ is normal absorption (a scalar $e^{-2\pi t/\xi'_{0}}$); $\xi'_{g}$ for $g \neq 0$ is anomalous absorption. |
| $\mathbf{A}$ | Dynamical structure matrix of the coupled beam equations $\mathrm{d}\psi/\mathrm{d}z = i\pi\mathbf{A}\psi$: diagonal $2s_{g} + i/\xi'_{0}$, off-diagonal $\nu_{g-h} + i/\xi'_{g-h}$. |
| $\gamma_{j}$ | Bloch-wave excitation: the $j$-th eigenvalue of $\mathbf{A}$. Its real part locates the dispersion surface; its imaginary part is that branch's absorption coefficient. |
| $I_{\mathrm{dd}}$ | Double-diffraction intensity assigned to a kinematically forbidden reflection reached as $\mathbf{g}_{1} + \mathbf{g}_{2}$. An observability estimate, never a kinematic intensity: the kinematic intensity of such a reflection is exactly zero. |
| $c$ | Double-diffraction coupling constant scaling $I_{\mathrm{dd}}$. Absorbs what a kinematic treatment cannot supply — beam coupling strength and specimen thickness. Dimensionless, in $(0, 1]$. |

## References

### Normative

- [Notation and Conventions](notation_and_conventions.md)
- [Executable Worked Examples](executable_examples.md)
- [Reference Canon](reference_canon.md)

### Informative

- <a href="../site/concepts/technical_glossary_and_symbols.md">Technical Glossary and Symbols</a>
