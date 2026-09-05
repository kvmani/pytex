# Precise Lattice-Parameter Determination

This note derives the methods PyTex uses to determine a unit cell from a measured powder
diffractogram, and explains why the obvious method fails. It covers peak detection and profile
fitting, the aberrations of a laboratory diffractometer, K$\alpha_2$ treatment, indexing and its
figures of merit, and the two determination methods — Cohen least squares on the reciprocal metric
tensor, and Le Bail whole-pattern decomposition.

The implementation lives in `pytex.diffraction.xrd_peaks`, `pytex.diffraction.xrd_corrections`,
`pytex.diffraction.xrd_indexing` and `pytex.diffraction.xrd_lattice_parameter`.

**This is determination, not refinement.** Nothing here varies an atomic coordinate, a thermal
parameter, or a site occupancy. The structure is held fixed and only the cell and the errors of the
instrument are determined. That restriction is the point: it is what stops texture and an imperfect
structural model from leaking into the answer. For structure refinement see
{doc}`../algorithms/index` and `pytex.diffraction.rietveld`.

## 1. Why averaging over reflections does not work

The intuitive method is to compute a lattice parameter from each reflection and average. It is
worth understanding precisely, because it fails for a reason that no amount of extra data repairs.

Differentiating Bragg's law $\lambda = 2 d \sin\theta$ at fixed $\lambda$,

$$
0 = 2\,\Delta d \, \sin\theta + 2 d \cos\theta \, \Delta\theta
\qquad\Longrightarrow\qquad
\frac{\Delta d}{d} = -\cot\theta \, \Delta\theta .
$$

Two consequences follow, and they pull in opposite directions.

The **favourable** one is that $\cot\theta \to 0$ as $\theta \to 90^{\circ}$. A given angular error
does progressively less damage to the spacing as the reflection moves towards back-reflection, so
high-angle reflections are intrinsically more informative. This is why the weighting in
{ref}`section 6 <lattice-cohen>` is not a convenience but a statement of physics.

The **unfavourable** one is that the error in $d$ therefore *depends on $\theta$*. The errors that
dominate a laboratory scan are systematic, not random — a misaligned detector zero, a specimen
sitting a few tens of micrometres off the diffractometer axis, a beam that penetrates before it
diffracts. Averaging $N$ reflections divides the random scatter by $\sqrt{N}$ and leaves the
systematic part exactly where it was. The result is a lattice parameter good to about one part in
$10^{3}$, roughly two orders of magnitude short of what thermal-expansion and strain work needs.

The library keeps the averaging method (`method="average"`) precisely so that this comparison can
be made on the reader's own data. It is available for cubic cells only, and that restriction is
itself instructive: outside the cubic system a *lattice parameter per reflection* does not exist,
because one reflection cannot determine two cell parameters.

## 2. The aberrations, and their angular signatures

Four systematic effects move peaks on a Bragg–Brentano diffractometer. They can be told apart, and
therefore modelled, only because their angular dependences differ.

| Aberration | $\Delta(2\theta)$ | Typical size |
| --- | --- | --- |
| Detector zero $z$ | $z$ (constant) | $0.01$–$0.05^{\circ}$ |
| Specimen displacement $s$ | $-2 s \cos\theta / R$ | $50\ \mu\mathrm{m}$ gives $0.03^{\circ}$ |
| Specimen transparency | $-\sin 2\theta / (2\mu R)$ | large for low $\mu$ |
| Refraction | $4\delta / \sin 2\theta$ | about $1$ part in $10^{5}$ |

The displacement term is the one that governs precise work on a well-made instrument. Its
$\cos\theta$ dependence means it **vanishes at back-reflection**, which is the entire basis of the
extrapolation methods below. The transparency term peaks at $2\theta = 90^{\circ}$ and vanishes at
both ends; refraction does the opposite, growing towards low angle.

The refraction decrement follows the classical Drude result away from an absorption edge,

$$
\delta = \frac{r_e \lambda^{2} n_e}{2\pi},
\qquad n_e = \rho N_A \frac{Z}{A},
$$

with $Z/A \approx 0.5$ for every element but hydrogen. For nickel at Cu K$\alpha$ this is
$\delta \approx 2.9 \times 10^{-5}$ — the same order as the precision a Cohen determination
reaches, so it is not negligible there even though it is invisible in phase identification.

**Zero and displacement are not refined together from one specimen scan.** They are separable in
principle, being constant against $\cos\theta$, and badly correlated in practice over the angular
range of a single pattern. The convention here is that zero belongs to a calibrated instrument
(`InstrumentBroadening` carries one) and displacement belongs to the specimen.

## 3. Finding the peaks, and how well their positions are known

Everything downstream needs two numbers per reflection: where the peak is, and how well that is
known. A lattice parameter quoted without an uncertainty cannot be compared with anything.

### Detection matched to the instrument

A peak is not "a local maximum above a threshold". It is a feature whose *width* matches what the
diffractometer can produce at that angle. Detection therefore proceeds by:

1. removing the background with SNIP;
2. applying the Anscombe transform $y = 2\sqrt{I + 3/8}$, which makes Poisson noise homoscedastic at
   every count level, so one threshold applies across a pattern whose strongest and weakest peaks
   differ by decades;
3. convolving with zero-mean, unit-$L^{2}$-norm Ricker kernels (the negative second derivative of a
   Gaussian) over a grid of scales, and interpolating the response at each angle to *that angle's*
   expected FWHM from the calibrated Caglioti curve.

Unit $L^{2}$ norm is what makes the threshold meaningful: white noise of unit variance produces a
response of unit variance, so the detection threshold is quotable in noise standard deviations
rather than in counts. A feature narrower than the instrument can produce — a cosmic-ray spike, a
dead channel — or broader than it — background structure, an amorphous halo — gives a weak response
and is rejected with no threshold expressed in angle units at all.

### The K$\alpha_2$ line is a peak

Both lines of the doublet diffract from the same planes, so at fixed $d$,

$$
\sin\theta_{2} = \frac{\lambda_{2}}{\lambda_{1}} \sin\theta_{1},
\qquad
\Delta(2\theta) = 2\,\frac{\Delta\lambda}{\lambda}\tan\theta .
$$

The separation grows as $\tan\theta$: the pair is unresolved at low angle and cleanly split above
roughly $90^{\circ}$, where the filter reports the $\alpha_2$ line as a peak in its own right.
Admitting it to the candidate list is not a cosmetic problem — it is a reflection placed at the
wrong $d$ spacing, and therefore a lattice-parameter bias. Each accepted candidate therefore
suppresses the angle its own $\alpha_2$ partner occupies.

### Fitting, and the position uncertainty

Each candidate is fitted in its own window with a pseudo-Voigt or split pseudo-Voigt profile, a
straight local background, and — when the radiation declares a second line — an $\alpha_2$ partner
whose position is *not free*, being fixed by the relation above. The reported uncertainty is
$\sqrt{[(\mathbf{J}^{\mathsf{T}}\mathbf{W}\mathbf{J})^{-1}]_{11}}$ scaled by the window's reduced
$\chi^{2}$, the goodness-of-fit-scaled convention used throughout crystallography.

## 4. K$\alpha_2$: model it, do not strip it

Rachinger's method removes the $\alpha_2$ contribution by the recursion

$$
I_{1}(2\theta) = I(2\theta) - r\,I_{1}(2\theta'),
\qquad
\sin\theta' = \frac{\lambda_{1}}{\lambda_{2}} \sin\theta ,
$$

sweeping upward in angle, with $r = I_2/I_1 \approx 1/2$. Because $\lambda_1 < \lambda_2$, the
angle $2\theta'$ always lies below $2\theta$, so the recursion only ever needs values it has
already computed.

Its three assumptions are exactly its three weaknesses. The ratio is taken as exact; the two
profiles are taken as identically shaped, when in truth the $\alpha_2$ profile is slightly wider in
$2\theta$ because $\mathrm{d}(2\theta_2)/\mathrm{d}(2\theta_1) \neq 1$; and each subtraction feeds
the next, so the variance propagates as

$$
\operatorname{var} I_{1}(2\theta) = \operatorname{var} I(2\theta)
+ r^{2} \operatorname{var} I_{1}(2\theta') ,
$$

growing with angle. The residual that survives is re-subtracted at the partner of *that* angle and
again beyond it, producing a decaying ringing train above each strong reflection. On a synthetic
nickel pattern the method removes 86 to 94 per cent of each $\alpha_2$ line — enough to make the
pattern much easier to read, and not enough to fit.

PyTex therefore implements stripping for **display and peak picking**, and models the doublet on
every fitting path. Modelling makes the same two physical assumptions but makes them *inside* the
model, where they are visible in the residual, instead of writing their failure into the data.

## 5. Indexing, and its figures of merit

Until a peak carries an $(hkl)$ it contributes nothing, because
$\sin^{2}\theta = (\lambda^{2}/4)\,\mathbf{h}^{\mathsf{T}}\mathbf{G}^{*}\mathbf{h}$ needs the
$\mathbf{h}$ as much as the angle.

The assignment is **global, not greedy**. Walking the peak list and taking the nearest calculated
line can assign two peaks to the same reflection, or pair a peak with a near neighbour and strand
the true partner — and neither failure is visible in the result it produces. The Hungarian
algorithm minimises the total discrepancy over all one-to-one pairings at once.

Two figures of merit are reported. de Wolff's

$$
M_{N} = \frac{Q_{N}}{2\,\langle|\Delta Q|\rangle\,N_{\mathrm{poss}}},
\qquad Q = 1/d^{2},
$$

and Smith and Snyder's

$$
F_{N} = \frac{N}{\langle|\Delta 2\theta|\rangle\,N_{\mathrm{poss}}} .
$$

In both, $N_{\mathrm{poss}}$ — the number of distinct reflections the cell predicts below the
$N$-th observed line — is the factor that does the work. Halving the discrepancies doubles the
figure of merit, but so does a cell that predicts half as many lines while fitting equally well. A
cell that "explains" every peak because it predicts a reflection everywhere scores badly, and that
warning is as useful for a known phase as for an unknown one.

Both are returned with the $N$ they were computed over, because $M_{7}$ and $M_{20}$ are not
comparable and a bare number invites treating them as if they were.

### Re-indexing is not optional

A starting cell wrong by a fraction $e$ misplaces a reflection by
$\Delta(2\theta) = 2 e \tan\theta$, which grows without bound towards back-reflection. A cell wrong
by three parts in a thousand — well within the range of a real alloy against a tabulated
composition — misplaces a $121^{\circ}$ reflection by $0.6^{\circ}$, far outside any sane indexing
tolerance. Those reflections are then silently dropped, and they are precisely the ones that carry
almost all the precision. The pipeline therefore indexes, determines, and **indexes again against
the cell just determined**, without which the tool would require the answer in advance.

(lattice-cohen)=
## 6. Cohen least squares on the reciprocal metric tensor

### The linear form

Write $\mathbf{h} = (h, k, l)^{\mathsf{T}}$. Then $1/d^{2} = \mathbf{h}^{\mathsf{T}}
\mathbf{G}^{*}\mathbf{h}$ with $\mathbf{G}^{*}$ the reciprocal metric tensor, and Bragg's law gives

$$
\sin^{2}\theta = \frac{\lambda^{2}}{4}\,\mathbf{h}^{\mathsf{T}} \mathbf{G}^{*} \mathbf{h} .
$$

This is **linear in the components of $\mathbf{G}^{*}$**. Expanding the quadratic form as a dot
product,

$$
\mathbf{h}^{\mathsf{T}} \mathbf{G}^{*} \mathbf{h} =
\begin{pmatrix} h^{2} & k^{2} & l^{2} & 2hk & 2hl & 2kl \end{pmatrix}
\begin{pmatrix}
G^{*}_{11} & G^{*}_{22} & G^{*}_{33} & G^{*}_{12} & G^{*}_{13} & G^{*}_{23}
\end{pmatrix}^{\mathsf{T}} .
$$

The crystal system supplies a $6 \times n$ constraint matrix reducing those six components to the
genuinely free parameters — one for cubic, two for tetragonal, trigonal and hexagonal, three for
orthorhombic, four for monoclinic, six for triclinic. One implementation therefore serves every
crystal system, with no starting guess, no local minima, and an analytic covariance matrix.

The hexagonal row is the one worth reading. There $a^{*} = b^{*}$ and $\gamma^{*} = 60^{\circ}$,
so $G^{*}_{12} = a^{*}b^{*}\cos\gamma^{*} = G^{*}_{11}/2$, and

$$
\mathbf{h}^{\mathsf{T}} \mathbf{G}^{*} \mathbf{h}
= G^{*}_{11}(h^{2} + k^{2}) + G^{*}_{33} l^{2} + 2hk \cdot \tfrac{1}{2} G^{*}_{11}
= G^{*}_{11}\left(h^{2} + hk + k^{2}\right) + G^{*}_{33} l^{2},
$$

which is the familiar hexagonal expression $1/d^{2} = \tfrac{4}{3}(h^{2}+hk+k^{2})/a^{2} +
l^{2}/c^{2}$ — obtained here as a *consequence* of the symmetry rather than written in as a special
case.

The direct cell is recovered from $\mathbf{G} = (\mathbf{G}^{*})^{-1}$ by
$a = \sqrt{G_{11}}$ and $\cos\alpha = G_{23}/(bc)$, and so on, so one matrix inversion serves every
system here too.

### The systematic-error term

This is the identity that unifies the classical graphical extrapolations with the least-squares
treatment. Suppose an aberration produces a fractional spacing error
$\Delta d / d = -K f(\theta)$ for some function $f$. Since
$\sin^{2}\theta = \lambda^{2}/(4d^{2})$,

$$
\frac{\Delta(\sin^{2}\theta)}{\sin^{2}\theta} = -2\,\frac{\Delta d}{d} = 2 K f(\theta),
$$

and therefore

$$
\boxed{\;
\sin^{2}\theta_{\mathrm{obs}}
= \frac{\lambda^{2}}{4}\,\mathbf{h}^{\mathsf{T}} \mathbf{G}^{*} \mathbf{h}
+ D\,\sin^{2}\theta\,f(\theta), \qquad D = 2K .
\;}
$$

So *whatever function one would plot $a$ against in the classical graphical method is the same
function that appears, multiplied by $\sin^{2}\theta$, as a design column here.* Because every
admissible $f$ vanishes at $\theta = 90^{\circ}$, the fitted cell is the extrapolated one.

Each aberration has its own matching $f$, obtained by writing its $\Delta(2\theta)$ in this form
using $\Delta(\sin^{2}\theta) = \sin\theta\cos\theta\,\Delta(2\theta)$:

| Aberration | Matching $f(\theta)$ | Derivation |
| --- | --- | --- |
| Detector zero | $\cot\theta$ | $\Delta(\sin^{2}\theta) = \sin\theta\cos\theta\,z = \sin^{2}\theta\cdot\cot\theta\cdot z$ |
| Specimen displacement | $\cos^{2}\theta / \sin\theta$ | exact for $-2s\cos\theta/R$ |
| Absorption and displacement | Nelson–Riley | $\tfrac{1}{2}\left(\dfrac{\cos^{2}\theta}{\sin\theta} + \dfrac{\cos^{2}\theta}{\theta}\right)$ |
| Camera absorption | $\cos^{2}\theta$ (Bradley–Jay) | gives Cohen's classical column |

The zero-error case deserves the derivation shown, because a *constant* angular offset looks at
first as though it could not be extrapolated away at all. It can: the constant in $2\theta$ becomes
$\sin^{2}\theta \cot\theta$ in $\sin^{2}\theta$, and $\cot 90^{\circ} = 0$. The fitted $D$ is then
the zero error itself, in radians.

Taking $f = \cos^{2}\theta$ reproduces Cohen's classical $\sin^{2}(2\theta)$ drift column exactly,
since $\sin^{2}\theta\cos^{2}\theta = \sin^{2}(2\theta)/4$.

**One drift term removes one aberration.** A pattern carrying both a zero error and a specimen
displacement cannot be fully corrected by any single $f$, which is exactly why precise work
calibrates the zero against a standard first and refines only the displacement afterwards.

### Weighting and uncertainties

Observations are weighted by $1/\sigma^{2}(\sin^{2}\theta)$ with

$$
\sigma(\sin^{2}\theta) = \left|\frac{\mathrm{d}\sin^{2}\theta}{\mathrm{d}(2\theta)}\right|
\sigma(2\theta) = \tfrac{1}{2}\sin(2\theta)\,\sigma(2\theta)
$$

propagated from each peak's own fitted position uncertainty. High-angle reflections dominate, which
is the quantitative form of $\Delta d/d = -\cot\theta\,\Delta\theta$.

The covariance is $(\mathbf{X}^{\mathsf{T}}\mathbf{W}\mathbf{X})^{-1}$ scaled by the reduced
$\chi^{2}$, and the cell uncertainties are propagated through the reciprocal-to-direct inversion by
$\sigma^{2} = \mathbf{J}\mathbf{C}\mathbf{J}^{\mathsf{T}}$ with a numerically differentiated
$\mathbf{J}$ — numerical rather than analytic because the analytic form differs per crystal system
and the numerical one demonstrably does not.

## 7. Le Bail whole-pattern decomposition

Single-peak fitting runs out of *resolvable* peaks long before a hexagonal or lower-symmetry
pattern runs out of reflections. Whole-pattern decomposition uses every measured point instead.

The calculated profile is a sum of peaks at the positions the current cell puts them, each scaled
by an intensity that is **not** a refined parameter:

1. every reflection starts with equal intensity;
2. **extraction** — the observed intensity at each point is partitioned between the reflections
   overlapping it, in proportion to their current contributions,

$$
I_{k} \leftarrow \sum_{i} y_{\mathrm{obs}}(i)\,
\frac{I_{k} P_{k}(i)}{\sum_{j} I_{j} P_{j}(i)} ;
$$

3. **refinement** — with the intensities held at their extracted values, the cell, the systematic
   term, the Caglioti coefficients and the mixing parameter are refined by bounded
   Levenberg–Marquardt;
4. repeat.

Step 2 is the whole trick. A Pawley fit instead treats the intensities as free least-squares
parameters, whose normal matrix becomes singular exactly when two reflections overlap completely —
which is the case the method exists to handle. Le Bail's partition is stable there, because two
exactly coincident reflections simply split the intensity in their current ratio and neither the
cell nor the fit notices.

Because the intensities are extracted rather than modelled, **neither texture nor a wrong atomic
basis can bias the cell**. The converse is that Le Bail intensities are fine for describing a
profile and unfit for structural work: for two reflections that overlap completely, the partition
between them is whatever ratio the iteration started with.

Three implementation points are load-bearing, and each of them produced a *wrong cell* rather than
merely a poor fit when got wrong:

- **The profile of one reflection is the whole K$\alpha$ multiplet.** Modelling only $\alpha_1$
  against doublet data leaves a residual as large as the $\alpha_2$ peak itself, which the
  refinement absorbs into the cell.
- **Weights come from the measured counts, never from the background-subtracted profile.**
  Subtraction removes signal, not variance: a background point with 150 counts still has
  $\sigma \approx 12$, and weighting it as though $\sigma = 1$ over-weights everything between the
  peaks by an order of magnitude.
- **The extraction returns an integrated intensity, not an amplitude.** The partition sums observed
  counts, so each reflection's profile must be normalised to unit sum; multiplying a unit-height
  profile by the extracted value overstates every peak by the number of points beneath it.

The background's *shape* is removed with SNIP beforehand, so that a flexible background cannot
absorb peak intensity and shift the cell with it. What SNIP leaves is a small level offset, and that
is refined as a straight line in the reduced angular coordinate — two parameters, far too stiff to
follow a peak.

## 8. What this is not: strain is not stress

A determination reports a cell and, against a reference cell, a lattice strain. It does not report
a stress, and it must not be read as though it did.

A symmetric $\theta$–$2\theta$ scan measures the spacing of planes whose normal is parallel to the
scattering vector — that is, normal to the specimen surface. One such measurement gives one strain
component. Converting to stress requires measuring $d(hkl)$ at several specimen tilts $\psi$ and
using the slope of $d$ against $\sin^{2}\psi$ together with the X-ray elastic constants
$\tfrac{1}{2}S_{2}$ and $S_{1}$ of the particular reflection — which are $hkl$-dependent, because
the diffracting grains are a subset selected by orientation and are elastically anisotropic.

This module supplies the precise spacings that such an analysis consumes. The
`describe()` output of every result says so explicitly.

## 9. Precision reached

On synthetic patterns with a known cell and a deliberately injected aberration:

| Method | Relative error in $a$ |
| --- | --- |
| Single reflection | $\sim 10^{-3}$ |
| Average over reflections | $\sim 10^{-3}$ to $4 \times 10^{-4}$ |
| Joint least squares, no drift term | $\sim 10^{-4}$ |
| Cohen with a mismatched $f$ | $\sim 10^{-5}$ |
| Cohen with the matching $f$ | $\sim 10^{-7}$ |
| Le Bail with the matching systematic term | $\sim 10^{-6}$ |

Elastic strains of engineering interest are $10^{-4}$ to $10^{-3}$, so the difference between the
first row and the last is the difference between measuring a strain and not.

The caveat those numbers carry is that a *real* specimen carries more than one aberration at once,
and no single drift term removes two of different angular form. The recommended practice remains
the one Cullity gives: calibrate the zero against a standard, keep it with the instrument, refine
only the specimen displacement, and use high-angle reflections.

## References

Cullity, B. D. & Stock, S. R., *Elements of X-Ray Diffraction*, 3rd ed., Prentice Hall (2001),
Ch. 11 "Precise Parameter Measurements".

Nelson, J. B. & Riley, D. P., *Proc. Phys. Soc.* **57** (1945) 160–177,
[doi:10.1088/0959-5309/57/3/302](https://doi.org/10.1088/0959-5309/57/3/302).

Bradley, A. J. & Jay, A. H., *Proc. Phys. Soc.* **44** (1932) 563–579,
[doi:10.1088/0959-5309/44/5/305](https://doi.org/10.1088/0959-5309/44/5/305).

Cohen, M. U., *Rev. Sci. Instrum.* **6** (1935) 68–74,
[doi:10.1063/1.1751937](https://doi.org/10.1063/1.1751937).

Rachinger, W. A., *J. Sci. Instrum.* **25** (1948) 254–255,
[doi:10.1088/0950-7671/25/7/125](https://doi.org/10.1088/0950-7671/25/7/125).

Wilson, A. J. C., *Mathematical Theory of X-ray Powder Diffractometry*, Philips Technical Library
(1963).

Le Bail, A., Duroy, H. & Fourquet, J. L., *Mater. Res. Bull.* **23** (1988) 447–452,
[doi:10.1016/0025-5408(88)90019-0](https://doi.org/10.1016/0025-5408(88)90019-0).

Pawley, G. S., *J. Appl. Crystallogr.* **14** (1981) 357–361,
[doi:10.1107/S0021889881009618](https://doi.org/10.1107/S0021889881009618).

de Wolff, P. M., *J. Appl. Crystallogr.* **1** (1968) 108–113,
[doi:10.1107/S002188986800508X](https://doi.org/10.1107/S002188986800508X).

Smith, G. S. & Snyder, R. L., *J. Appl. Crystallogr.* **12** (1979) 60–65,
[doi:10.1107/S002188987901178X](https://doi.org/10.1107/S002188987901178X).

Toraya, H., *J. Appl. Crystallogr.* **19** (1986) 440–447,
[doi:10.1107/S0021889886088982](https://doi.org/10.1107/S0021889886088982).

Thompson, P., Cox, D. E. & Hastings, J. B., *J. Appl. Crystallogr.* **20** (1987) 79–83,
[doi:10.1107/S0021889887087090](https://doi.org/10.1107/S0021889887087090).

Bearden, J. A., *Rev. Mod. Phys.* **39** (1967) 78–124,
[doi:10.1103/RevModPhys.39.78](https://doi.org/10.1103/RevModPhys.39.78).

Du, P., Kibbe, W. A. & Lin, S. M., *Bioinformatics* **22** (2006) 2059–2065,
[doi:10.1093/bioinformatics/btl355](https://doi.org/10.1093/bioinformatics/btl355).

Noyan, I. C. & Cohen, J. B., *Residual Stress: Measurement by Diffraction and Interpretation*,
Springer (1987) — the $\sin^{2}\psi$ analysis this note declines to perform.
