# Phase Identification From Powder Patterns

This note derives what it means for a candidate crystal structure to *account for* a measured
powder diffractogram, and why ranking several candidates is a different problem from indexing one
of them well. It covers the four scoring criteria PyTex uses, why measured intensities are the
weakest of the four kinds of evidence, why one uniform cell dilation may be refined without
letting a wrong structure be stretched into fitting, and the two qualifications that separate
"this candidate won" from "this is the phase".

The implementation lives in `pytex.diffraction.xrd_phase_identification`, on top of
`pytex.diffraction.xrd_peaks` for detection and fitting and `pytex.diffraction.xrd_indexing` for
the per-candidate assignment. The stage-by-stage implementation account is in
{doc}`../algorithms/phase_identification`.

**This is ranking, not retrieval, and not quantification.** The candidates must be supplied — as
CIF files or from the built-in catalogue. Nothing here searches a database for structures nobody
proposed, and nothing here apportions a mixture between phases. Both are named in §7 with what to
do instead.

## 1. Why indexing well is not the same as being right

{doc}`precise_lattice_parameter_determination` derives the two classical figures of merit for an
indexed pattern. de Wolff's

$$
M_N = \frac{Q_N}{2 \langle |\Delta Q| \rangle \, N_{\text{poss}}},
\qquad Q = \frac{1}{d^2},
$$

and Smith and Snyder's

$$
F_N = \frac{N}{\langle |\Delta 2\theta| \rangle \, N_{\text{poss}}}
$$

both ask how closely a cell places the lines that *were* observed, discounted by how many lines the
cell predicts. They were designed for autoindexing, where the question is whether a cell derived
from the data is believable, and for that they are exactly right.

They are the wrong instrument for choosing between candidates, and it is worth being precise about
why. Both are computed over the indexed reflections only. Neither term in either expression
changes when

- a strong measured peak is left unindexed, because it contributes to no $\Delta Q$ and to no $N$;
  or
- a reflection the candidate predicts at full intensity is not observed at all, because
  $N_{\text{poss}}$ counts predicted lines whether or not they appeared, and increasing it lowers
  the figure only in proportion.

Those two are the whole content of the fcc/bcc distinction. A face-centred cubic cell of parameter
$a_F$ and a body-centred cubic cell of parameter $a_I$ produce lines at

$$
\sin^2\theta = \frac{\lambda^2}{4 a^2}\,(h^2+k^2+l^2),
$$

with $h,k,l$ all odd or all even for $F$, and $h+k+l$ even for $I$. Both sequences are
$\lambda^2/(4a^2)$ times a set of integers, so for a suitable ratio $a_I / a_F$ several lines of
one coincide with lines of the other. What never coincides is the *pattern of absences*: $F$ has
no $(100)$, $(110)$, $(210)$; $I$ has no $(100)$, $(111)$, $(210)$. Identification is therefore
mostly an argument about lines that are not there, and a figure of merit computed over the lines
that are cannot make it.

## 2. The four criteria

Each candidate is indexed against the fitted peaks by the global one-to-one assignment of
{doc}`precise_lattice_parameter_determination` §5, and then scored on four quantities. Each is
bounded to $[0,1]$, each is reported separately, and each fails for a different physical reason —
which is why the separate values are more diagnostic than their combination.

### 2.1 Explained intensity

Let $\mathcal{P}$ be the fitted peaks, $\mathcal{I} \subseteq \mathcal{P}$ those the candidate
indexed, and $A_p$ the integrated intensity of peak $p$. Then

$$
E = \frac{\sum_{p \in \mathcal{I}} A_p}{\sum_{p \in \mathcal{P}} A_p} .
$$

This is weighted by intensity rather than counted, and the weighting carries the physics. A strong
unindexed peak is decisive evidence that the specimen contains something the candidate is not; a
weak one is routinely a trace impurity, a detector artefact, or a line from the holder, and
penalising it as heavily as a principal reflection would reject the correct phase of any real
specimen. The strongest unindexed peak is also reported on its own, as a fraction of the strongest
measured peak, because a single large value there is the signature of a second phase rather than of
a poor fit.

### 2.2 Completeness

Let $\mathcal{C}$ be the candidate's reflections inside the measured angular range whose calculated
relative intensity exceeds a threshold $\tau$ (default $\tau = 0.05$), and $\mathcal{O}
\subseteq \mathcal{C}$ those matched to a measured peak. Then

$$
C = \frac{|\mathcal{O}|}{|\mathcal{C}|} .
$$

This is the criterion of §1: it asks whether the candidate's own strong lines appeared. It is
evaluated **strictly inside the measured span**, not over a padded window, because a line predicted
beyond the last measured point was never looked for, and counting it as missing would penalise a
candidate for the operator's choice of scan range.

The threshold $\tau$ is a statement about detectability, not about crystallography. Raise it for a
noisy scan in which weak calculated lines genuinely could not have been seen.

### 2.3 Position agreement

With $\varepsilon$ the matching tolerance in degrees,

$$
P = \max\!\left(0,\; 1 - \frac{\langle |\Delta 2\theta| \rangle}{\varepsilon}\right).
$$

The distinction from indexing is that $P$ asks how far *inside* the window the assigned lines
landed rather than whether they landed inside it. A candidate whose every line sits at the edge of
the tolerance has the wrong cell dimensions even though every peak was formally indexed, and a
criterion that is binary in $\varepsilon$ cannot say so.

Because $\varepsilon$ sets the scale, widening the tolerance to admit more matches judges every one
of them against the laxer standard it was admitted under. That is deliberate: the tolerance is not
a free knob that improves a candidate the further it is opened.

### 2.4 Intensity agreement

Let $o_i$ and $c_i$ be the observed integrated and calculated relative intensities of the indexed
reflections, each normalised to unit sum. Then

$$
S = 1 - \tfrac{1}{2} \sum_i \left| \hat{o}_i - \hat{c}_i \right|,
\qquad \hat{o}_i = \frac{o_i}{\sum_j o_j},
\quad \hat{c}_i = \frac{c_i}{\sum_j c_j},
$$

which is one minus the Bray-Curtis dissimilarity. Normalising both to unit sum makes it invariant
to the arbitrary scale of a measured count rate and of a kinematic relative intensity; the factor
$\tfrac{1}{2}$ bounds it to $[0,1]$, since two unit-sum non-negative vectors differ by at most 2 in
$L^1$. It equals 1 exactly when the two distributions coincide.

This is the criterion that separates candidates the first three cannot: two structures on nearly
the same cell with different atomic bases — an ordered and a disordered alloy, the same framework
with a different cation — produce lines at nearly the same angles and differ almost entirely in how
strong those lines are.

With fewer than two indexed reflections $S$ is **undefined**, not zero: a similarity between two
one-element distributions is identically 1 by construction and would be a claim about nothing.

## 3. Why intensities are the weakest evidence

$S$ carries the least weight of the four by default, and the reason is physical rather than
statistical.

The kinematic intensity of a powder reflection is

$$
I_{hkl} \;\propto\; |F_{hkl}|^2 \, m_{hkl} \, L(\theta) \, P(\theta) \, e^{-2M} ,
$$

with $m$ the multiplicity and $L$, $P$ the Lorentz and polarisation factors. Every term is a
property of the structure and the geometry. A *specimen* interposes several effects that are
properties of neither:

- **Preferred orientation.** A rolled sheet, a coating, or a powder of platelets presents some
  plane normals to the diffraction vector far more often than a random powder would. In the March
  model the intensity of the $hkl$ reflection is scaled by a factor that depends on the angle
  between its normal and the texture axis, and for a strong texture this is a factor of several.
  See {doc}`preferred_orientation_in_powder_intensities`.
- **Microabsorption.** In a multi-phase specimen a strongly absorbing phase shields its own
  interiors, so its measured intensity is systematically low by an amount depending on grain size.
- **Extinction**, in well-crystallised grains, attenuates strong reflections preferentially.
- **Coarse or insufficiently randomised powder**, where too few grains satisfy the diffraction
  condition for the intensity to be a statistical average at all.

Every one of these moves an intensity while moving no peak position whatsoever. A ranking that
weighted $S$ heavily would reject the correct phase of any textured specimen, which is most
engineering specimens. Hence the default weighting

$$
w = (E, C, P, S) \cdot (0.40,\; 0.25,\; 0.20,\; 0.15),
$$

and hence the fact that the weights are a declared, overridable parameter rather than a constant:
they encode a judgement about which evidence is trustworthy *for this specimen*, and a laboratory
measuring rolled sheet is right to set $w_S = 0$.

The combined score is the weighted mean over the criteria that are defined,

$$
\text{score} = \frac{\sum_{k \in \mathcal{D}} w_k \, x_k}{\sum_{k \in \mathcal{D}} w_k},
\qquad
\mathcal{D} = \{k : x_k \text{ is defined}\},
$$

so an undefined criterion is renormalised away rather than counted as a failure. "Not measurable
here" and "measured and bad" are different findings and only the second is evidence against a
candidate.

## 4. One uniform cell dilation, and why it is safe

A CIF records the cell of the specimen whoever deposited it measured. A real specimen is a solid
solution, or sits at a different temperature, or carries a residual stress, and differs from the
tabulated cell by a fraction of a per cent. From §1 of
{doc}`precise_lattice_parameter_determination`, a relative cell error $e$ displaces a line by

$$
\Delta(2\theta) = 2 e \tan\theta ,
$$

which diverges towards back-reflection. Three parts in a thousand — ordinary for an alloy against
a pure-element entry — is 0.6° at $2\theta = 121°$, twice any sensible matching tolerance. Without
correction the true phase loses precisely the high-angle lines that would have confirmed it, and
the identification fails on the specimens it most needs to succeed on.

PyTex therefore refines a single scalar $s$ per candidate before indexing, minimising

$$
\Phi(s) = \sum_{p \in \mathcal{P}} \min\!\left(
\min_{j} \left| 2\theta_p - 2\theta_j(s) \right| ,\; \varepsilon \right),
\qquad
\sin\theta_j(s) = \frac{\lambda}{2 s\, d_j},
$$

over a grid on $[1 - \delta, 1 + \delta]$ with $\delta = 0.02$ by default. A grid rather than a
gradient method because $\Phi$ is piecewise linear with a local minimum at every near-coincidence,
so a descent lands in whichever one it started next to.

**This cannot make a wrong structure fit**, and the argument is exact. Under $d \mapsto s d$ every
ratio $d_{hkl} / d_{h'k'l'}$ is unchanged identically, for every pair, whatever the structure. Those
ratios are what indexing tests: from $\sin^2\theta = (\lambda^2/4)\, \mathbf{h}^{\mathsf T}
\mathbf{G}^{*} \mathbf{h}$, a uniform dilation scales $\mathbf{G}^{*}$ by $s^{-2}$ and leaves the
*relative* pattern of $\sin^2\theta$ values untouched. A candidate whose relative line positions
are wrong is wrong at every $s$; one that a scale factor rescues is the right structure with the
wrong cell size, which is the case this exists for. The identity is checked as an executable
worked example.

The refined $s$ is reported per candidate rather than applied silently, and it is informative in
both directions: a few parts in a thousand is a statement about composition or temperature, while
a value pinned at the edge of the search range means the candidate was stretched as far as it was
permitted and still did not fit.

## 5. Two qualifications on the winner

A ranking always has a winner. That is not the same as having an answer, so two further statements
are made, and made separately because they fail for different reasons and have different remedies.

**Conclusive** — the best score reaches an acceptance threshold (0.55 by default). Below it, the
honest reading is that *none* of the candidates offered accounts for this pattern. The remedies are
to widen the candidate list, to consider that the specimen is a mixture, or to check that the
matching tolerance exceeds the instrument's uncorrected zero-point and displacement errors.

**Decisive** — the winner leads the runner-up by at least a margin (0.05 by default). Below it, the
two are not distinguished *by this scan*, which is a statement about the measurement rather than
about the candidates. The remedies are more measurement, not more computation: a longer count at
high angle where the two candidates' calculated lines diverge fastest (by §4, the divergence grows
as $\tan\theta$), a different wavelength, or independent chemistry.

Reporting only a ranked list would leave both of these for the reader to infer, and the second is
easy to miss: a winner at 0.71 and a runner-up at 0.70 looks like a result.

## 6. What the criteria diagnose

The value of reporting the four separately is that the *pattern of failure* names the fault.

| Symptom | Reading |
| --- | --- |
| High $P$, low $C$ | Right cell metric, wrong centring or basis — the classic fcc-offered-as-bcc |
| High $C$ and $P$, low $E$ | The candidate is present, but so is something else: a second phase |
| Low $P$, moderate everything else | Wrong cell dimensions, or a tolerance narrower than the instrument's aberrations |
| High $E$, $C$, $P$, low $S$ | Right framework, wrong basis — or a textured specimen, so check §3 before concluding |
| $s$ pinned at $\pm\delta$ | Stretched as far as permitted and still not fitting: read the candidate with suspicion |
| Everything high for two candidates | Not distinguished by this scan; see §5 |

## 7. What this deliberately is not

**Not a database search.** Hanawalt, Rinn and Frevel's method indexes hundreds of thousands of
reference patterns by their three strongest $d$ spacings so that candidates can be *retrieved* from
an unknown pattern; every modern search-match system descends from it. Retrieval and ranking are
separate problems, and the retrieval half requires a licensed reference database — the Powder
Diffraction File or an equivalent — that this library does not ship. Here the candidates are
already in hand because the user chose them, and only the ranking is performed. A consequence
worth stating to users: a low best score is as likely to mean the right structure was never offered
as that the scan is poor.

**Not quantitative phase analysis.** When several candidates each explain part of a pattern, the
ranked list plus the unindexed peaks is the honest output. Apportioning the specimen between phases
requires their scale factors refined jointly against the whole profile — a multi-phase Rietveld
refinement, {doc}`../algorithms/rietveld_refinement` — not a larger score. The identification is
what tells that refinement which phases to include.

**Not a structure determination.** Nothing here varies an atomic coordinate, a site occupancy or a
thermal parameter. The candidates are taken as given and only compared.

## 8. References

**Normative.**

- Hanawalt, J. D., Rinn, H. W. & Frevel, L. K., "Chemical Analysis by X-Ray Diffraction",
  *Ind. Eng. Chem. Anal. Ed.* **10** (1938) 457-512,
  [doi:10.1021/ac50125a001](https://doi.org/10.1021/ac50125a001). The search-match method, of which
  the criteria here are the already-have-the-candidates case.
- Smith, G. S. & Snyder, R. L., "$F_N$: a criterion for rating powder diffraction patterns",
  *J. Appl. Crystallogr.* **12** (1979) 60-65,
  [doi:10.1107/S002188987901178X](https://doi.org/10.1107/S002188987901178X).
- de Wolff, P. M., "A simplified criterion for the reliability of a powder pattern indexing",
  *J. Appl. Crystallogr.* **1** (1968) 108-113,
  [doi:10.1107/S002188986800508X](https://doi.org/10.1107/S002188986800508X).
- Dollase, W. A., "Correction of intensities for preferred orientation in powder diffractometry",
  *J. Appl. Crystallogr.* **19** (1986) 267-272,
  [doi:10.1107/S0021889886089458](https://doi.org/10.1107/S0021889886089458). Why §3 weights
  intensities least.

**Informative.**

- Cullity, B. D. & Stock, S. R., *Elements of X-Ray Diffraction*, 3rd ed., Prentice Hall (2001),
  Ch. 14 (Chemical Analysis by Diffraction).
- Gates-Rector, S. & Blanton, T., "The Powder Diffraction File: a quality materials
  characterization database", *Powder Diffr.* **34** (2019) 352-360,
  [doi:10.1017/S0885715619000812](https://doi.org/10.1017/S0885715619000812). The reference
  database §7 declines to replace.
- Kuhn, H. W., "The Hungarian method for the assignment problem", *Naval Res. Logist. Quart.* **2**
  (1955) 83-97, [doi:10.1002/nav.3800020109](https://doi.org/10.1002/nav.3800020109). The
  assignment used per candidate.

## See also

- {doc}`../algorithms/phase_identification` — the implementation, stage by stage
- {doc}`precise_lattice_parameter_determination` — detection, fitting, indexing and the figures of
  merit this note builds on
- {doc}`preferred_orientation_in_powder_intensities` — why §3 distrusts measured intensities
- {doc}`powder_xrd_and_saed` — the forward model the candidates are calculated with
- {doc}`crystal_structures_and_cif_import` — how a candidate arrives from a CIF
- {doc}`../workflows/xrd_generation` — the workbench workflow
