# TEM Tilt Navigation Foundation

**Status:** **proposed — formulation for review.** No implementation exists. Program identifier:
**TN**. This document is the scientific formulation the user asked for before any code is written;
if approved it becomes the input to a formal algorithm and software-requirements specification.

**Scope of this document.** Define the problem, the geometry, the algorithm, the ambiguities, the
calibration burden, the reachability and path-planning logic, the visualization mathematics, the
error model, and the validation strategy — for taking a crystalline TEM specimen from a *known,
indexed* zone axis at a *known* holder position to a *requested* target zone axis, and reporting the
required $\alpha$ and $\beta$ holder tilts.

---

## 1. The question, stated precisely

An operator has reached some zone axis by chasing Kikuchi bands. The pattern is indexed. The stage
reads $(\alpha_c, \beta_c)$. The operator now wants zone axis $[uvw]_{\text{target}}$.

> Given the indexed current pattern, the current stage readout, the phase and its symmetry, the
> stage kinematic model and its calibration, and the holder tilt envelope — find every stage
> position $(\alpha, \beta)$ that places a crystallographically acceptable version of
> $[uvw]_{\text{target}}$ along the electron beam; rank them; validate them by forward simulation;
> plan a path to the best one; and state honestly what is *not* determined by the available data.

### 1.1 What this is not

- Not automated spot detection from raw image data. The current pattern arrives already indexed,
  either by `pytex.diffraction.solving.solve_saed_pattern` (picked spots plus a camera constant) or
  by `pytex.diffraction.models.index_saed_pattern` (calibrated detector).
- Not stage control. The output is numbers and a plan for a human or an external control layer. No
  instrument I/O is proposed.
- Not dynamical diffraction. Where an observable is needed (spot visibility, Kikuchi geometry) the
  existing kinematic engine supplies it, and that limitation is stated at the surface.
- Not a solution for a bent, buckled or polycrystalline region. The formulation assumes one rigid
  crystallite rigidly coupled to the holder over the excursion. §14.3 quantifies when that fails.

### 1.2 Why this is worth a foundation document rather than a function

The naive statement — "solve two angles" — has a three-line closed form (§6.1). Everything hard
lies elsewhere: which of the four mathematical branches and up to 48 symmetry images is the one you
can actually reach; which parts of the answer are determined by physics and which by an
uncalibrated instrument constant; and the fact that *the most damaging error mode produces a
perfectly self-consistent answer that tilts the specimen in exactly the wrong direction*. A
document that only derives the closed form would be a trap. §8 is the centre of gravity here.

---

## 2. Coordinate systems

### 2.1 Constraint inherited from the repository

`docs/standards/notation_and_conventions.md` fixes six frame domains — crystal, specimen, map,
detector, laboratory, reciprocal — and states that no subsystem may invent another or silently
collapse two. **TN introduces no new domain.** The mapping is:

| TN frame | Symbol | pytex `FrameDomain` | Definition |
| --- | --- | --- | --- |
| Crystal Cartesian | $C$ | `CRYSTAL` | Right-handed Cartesian tied to the lattice through the structure matrix $\mathbf{A}$ (columns $\mathbf{a},\mathbf{b},\mathbf{c}$). $[uvw] \mapsto \mathbf{A}[u,v,w]^{\mathsf T}$; $(hkl) \mapsto \mathbf{A}^{*}[h,k,l]^{\mathsf T}$ with $\mathbf{A}^{*} = (\mathbf{A}^{-1})^{\mathsf T}$ (Kronecker normalization). |
| Holder | $H$ | `SPECIMEN` | Rigidly attached to the specimen **and** to the holder cartridge. $\hat{\mathbf{x}}_H$ along the holder rod, $\hat{\mathbf{y}}_H$ along the nominal $\beta$ axis, $\hat{\mathbf{z}}_H = \hat{\mathbf{x}}_H \times \hat{\mathbf{y}}_H$. Coincides with the microscope frame at $(\alpha,\beta) = (0,0)$. |
| Microscope / laboratory | $L$ | `LABORATORY` | Fixed to the column. $\hat{\mathbf{z}}_L$ points **up the column, toward the gun**; electrons propagate along $-\hat{\mathbf{z}}_L$. $\hat{\mathbf{x}}_L$ along the nominal $\alpha$ axis (the rod). $\hat{\mathbf{y}}_L = \hat{\mathbf{z}}_L \times \hat{\mathbf{x}}_L$. |
| Pattern / detector | $P$ | `DETECTOR` | Axes of the **stored image array** of the diffraction pattern: $\hat{\mathbf{x}}_P$ along increasing column index, $\hat{\mathbf{y}}_P$ along the direction the display treats as "up", $\hat{\mathbf{z}}_P$ completing a right-handed set. Origin at the transmitted beam. |
| Image | $I$ | `DETECTOR` | Same construction for the real-space image. Related to $P$ by the image–diffraction rotation. |
| Reciprocal | — | `RECIPROCAL` | Unchanged; supplies $\mathbf{g}_{hkl}$. |

The `MAP` domain is unused by TN.

### 2.2 The specimen/holder collapse — a deliberate, load-bearing simplification

A textbook treatment inserts a *specimen* frame between crystal and holder, connected by a
**mounting rotation** $\mathbf{M}_{H \leftarrow S}$ (how the disc sits in the cup), and by the
crystal orientation $\mathbf{M}_{S \leftarrow C}$.

**These two are not separately identifiable from diffraction data, ever.** No observable in this
problem depends on them individually — only on the product. Carrying them separately therefore
manufactures an unobservable degree of freedom, which is the classic route to a sign error that
hides for months.

TN carries the single rotation

$$\mathbf{U} \;\equiv\; \mathbf{M}_{H \leftarrow S}\,\mathbf{M}_{S \leftarrow C} \;:\; C \rightarrow H,$$

the **crystal-to-holder orientation**. Because pytex defines an `Orientation` as a crystal→specimen
mapping, and because the holder frame is declared to be the specimen-domain frame for TEM work,
$\mathbf{U}$ *is* a pytex `Orientation` with no new concept required. Specimen mounting rotation
appears in §10 only as a nuisance parameter affecting the *tilt envelope* (a grid bar shadows at a
mounting-dependent azimuth), never as an orientation unknown.

### 2.3 Beam and alignment condition

The beam direction is $-\hat{\mathbf{z}}_L$. A crystal direction $\mathbf{t}$ is *on axis* when its
Cartesian image in the lab is parallel **or antiparallel** to $\hat{\mathbf{z}}_L$. Antiparallel is
admitted by default because a kinematic SAED pattern cannot distinguish the two senses (§8.1); a
`require_sense` option exists for polarity-sensitive work and is honoured only when the sense has
actually been determined.

---

## 3. Stage kinematics

### 3.1 The double-tilt holder, derived rather than asserted

The most commonly mis-stated point in this problem is whether the two tilt axes are fixed or
moving. Physically:

- $\alpha$ rotates the **holder rod** about its own axis. That axis is fixed in the column.
- $\beta$ rotates a cradle **carried inside the rod**. Its axis is fixed in the *holder*, and
  therefore moves when $\alpha$ changes.

Write $\mathbf{R}_{\hat{\mathbf n}}(\theta)$ for a right-handed rotation by $\theta$ about
$\hat{\mathbf n}$. In the lab frame the $\beta$ axis is $\mathbf{R}_x(\alpha)\hat{\mathbf y}$, so
the physical composition "first $\beta$ about the moving axis, then $\alpha$ about the fixed rod" is

$$
\mathbf{R}_{\text{stage}}
= \underbrace{\Big[\mathbf{R}_x(\alpha)\,\mathbf{R}_y(\beta)\,\mathbf{R}_x(\alpha)^{\mathsf T}\Big]}_{\beta \text{ about its lab-frame axis}} \mathbf{R}_x(\alpha)
= \mathbf{R}_x(\alpha)\,\mathbf{R}_y(\beta).
$$

$$\boxed{\;\mathbf{R}_{\text{stage}}(\alpha,\beta) = \mathbf{R}_x(\alpha)\,\mathbf{R}_y(\beta) \;:\; H \rightarrow L\;}$$

The moving-axis subtlety cancels exactly, which is *why* the naive fixed-axis composition happens to
give the right answer — a coincidence worth recording, because it does **not** survive a
non-orthogonal or mis-set axis pair (§3.3), and code that relies on the coincidence without knowing
it will break silently the first time a real calibration is applied.

### 3.2 The master equation

$$
\boxed{\;\mathbf{v}_L(\alpha,\beta) \;=\; \mathbf{R}_{\text{stage}}(\alpha,\beta)\,\mathbf{U}\,\mathbf{v}_C\;}
\tag{M}
$$

Everything in this document is a consequence of (M): reconstruction inverts it for $\mathbf{U}$
(§5), navigation inverts it for $(\alpha,\beta)$ (§6), reachability asks for its image (§10),
visualization plots it (§12), validation re-evaluates it forward (§14).

### 3.3 The general (calibrated) stage model

The ideal model above is the $\varepsilon \to 0$ limit of

$$
\mathbf{R}_{\text{stage}}(\alpha,\beta) = \mathbf{R}_{\hat{\mathbf a}}\!\big(s_\alpha(\alpha_{\text{read}} - \alpha_0)\big)\;
\mathbf{R}_{\hat{\mathbf b}(\alpha)}\!\big(s_\beta(\beta_{\text{read}} - \beta_0) + k\,\alpha\big),
$$

with

| Parameter | Meaning | Ideal | Source |
| --- | --- | --- | --- |
| $\hat{\mathbf a}$ | $\alpha$-axis direction in $L$ | $\hat{\mathbf x}_L$ | calibration (§9) |
| $\hat{\mathbf b}(\alpha)$ | $\beta$-axis direction in the $\alpha$-rotated holder | $\hat{\mathbf y}_H$ | calibration |
| $\varepsilon$ | non-orthogonality, $\hat{\mathbf a}\cdot\hat{\mathbf b}_0 = \sin\varepsilon$ | $0$ | calibration |
| $s_\alpha, s_\beta \in \{+1,-1\}$ | readout sign conventions | $+1$ | calibration (§9.3) |
| $\alpha_0, \beta_0$ | readout zero offsets | $0$ | calibration |
| $k$ | linear $\alpha \to \beta$ mechanical coupling | $0$ | calibration |

The engine must be written against the general model with the ideal model as its default instance.
The ideal model then supplies the **analytic seed** (§6.1) and the general model the **refined
solution** (§6.3). This two-stage structure is the single most important architectural decision in
the algorithm: it keeps a closed form available for reasoning, testing and visualization while
never letting the closed form be the answer a user acts on.

### 3.4 Other holder types

| Holder | Model |
| --- | --- |
| Double tilt | $\mathbf{R}_x(\alpha)\mathbf{R}_y(\beta)$ — above |
| Tilt–rotate | $\mathbf{R}_x(\alpha)\mathbf{R}_z(\theta)$ — $\theta$ about the holder normal |
| Single tilt | $\mathbf{R}_x(\alpha)$; the reachable set is a great circle, so §10 reduces to a 1-D search over the symmetry orbit |
| High-tilt tomography | $\mathbf{R}_x(\alpha)$ with $|\alpha| \le 70^\circ$ and a strongly $\alpha$-dependent envelope |

All four are instances of "an ordered product of two parameterized rotations about calibrated axes".
The stage model must be a small polymorphic interface, not an `if holder_type ==` chain.

### 3.5 Beam direction in holder coordinates — the working formula

Invert (M) for the beam:

$$
\hat{\mathbf b}_H(\alpha,\beta) \;=\; \mathbf{R}_{\text{stage}}(\alpha,\beta)^{\mathsf T}\hat{\mathbf z}_L
\;=\;\big(-\cos\alpha\,\sin\beta,\;\; \sin\alpha,\;\; \cos\alpha\,\cos\beta\big).
\tag{B}
$$

Equation (B) is remarkably clean and does a great deal of work below. Read it as spherical
coordinates **whose pole is the $\beta$-tilt axis $\hat{\mathbf y}_H$**: $\alpha$ is the latitude
measured from the $x$–$z$ plane, $\beta$ the longitude. Immediate consequences:

- The map $(\alpha,\beta) \mapsto \hat{\mathbf b}_H$ has Jacobian $\cos\alpha$; the solid angle
  swept by a tilt rectangle is $\Omega = (\beta_{\max}-\beta_{\min})(\sin\alpha_{\max} - \sin\alpha_{\min})$ (§10.2).
- Curves of constant $\beta$ are **great circles through $\pm\hat{\mathbf y}_H$**; curves of
  constant $\alpha$ are **small circles about $\hat{\mathbf y}_H$**. The accessible region is a
  geodesic quadrilateral in that coordinate system (§10.1, §12).
- $\lVert \partial\hat{\mathbf b}_H/\partial\alpha \rVert = 1$ but
  $\lVert \partial\hat{\mathbf b}_H/\partial\beta \rVert = \cos\alpha$: **a degree of $\beta$ buys
  only $\cos\alpha$ degrees of crystal rotation.** At $\alpha = 40^\circ$ the $\beta$ knob is 23 %
  less effective; at the $\hat{\mathbf y}_H$ pole it does nothing. This is the double-tilt holder's
  gimbal lock, and it is the origin of the conditioning result in §13.2.

---

## 4. Orientation representation

| Quantity | Representation | Rationale |
| --- | --- | --- |
| $\mathbf{U}$, stage rotations, symmetry operators | unit quaternion $(w,x,y,z)$, matrix on demand | pytex canonical storage; composition and interpolation are cheap and numerically stable |
| Solutions reported to the user | $(\alpha,\beta)$ in degrees | the only numbers an operator can act on |
| Path samples | quaternion **and** $(\alpha,\beta)$ **and** $\hat{\mathbf b}_C$ | each consumer needs a different one; storing all three is trivial and prevents re-derivation drift |
| Uncertainty | $3\times3$ covariance in the rotation-vector (exponential-map) tangent at the estimate | additive, composable, and the only representation in which "$2^\circ$ about the beam axis but $0.3^\circ$ transverse" — the actual anisotropy here — can be stated |
| Diagnostics only | axis–angle, Bunge Euler | human reading |

Euler angles are **never** used internally. The gimbal singularity of Bunge $\Phi \to 0$ is
unrelated to the physical gimbal lock of §3.5, and mixing the two in one debugging session is a
known source of wasted days.

---

## 5. Determining the current orientation $\mathbf{U}$

Three modes, in increasing order of reliability. The engine must support all three and **must
report which was used**, because the ambiguity content differs radically between them (§8).

### 5.1 Mode A — one indexed pattern plus a calibrated pattern rotation

`PatternSolution.orientation` already supplies $\mathbf{R}_{P \leftarrow C}$: the rotation carrying
crystal Cartesian vectors into the pattern frame, with $\hat{\mathbf z}_P$ along the zone axis
toward the viewer. The pattern frame relates to the lab by

$$\mathbf{R}_{L \leftarrow P} \;=\; \mathbf{F}\,\mathbf{R}_z(\varphi_D),$$

where $\varphi_D$ is the **diffraction rotation** (a function of camera length and accelerating
voltage) and $\mathbf{F} \in \{\mathbf{I},\ \mathrm{diag}(1,-1,1)\}$ encodes whether the stored
image is a mirrored rendering of the physical pattern. Then

$$
\boxed{\;\mathbf{U} = \mathbf{R}_{\text{stage}}(\alpha_c,\beta_c)^{\mathsf T}\;\mathbf{F}\,\mathbf{R}_z(\varphi_D)\;\mathbf{R}_{P \leftarrow C}\;}
\tag{A}
$$

Mode A is one line of code and is the mode most users will reach for. It is also the mode that
carries the entire instrumental ambiguity of §8.2–8.3. **If $\mathbf{F}$ is wrong, $\mathbf{U}$
comes out improper and $\det\mathbf{U} = -1$** — a free, exact, self-check that the implementation
must perform and must escalate rather than silently orthogonalize away.

### 5.2 Mode B — two indexed zone axes at two stage positions (recommended)

Suppose zone axes $\hat{\mathbf n}_1, \hat{\mathbf n}_2$ (crystal Cartesian unit vectors) were
indexed at stage positions $(\alpha_1,\beta_1)$ and $(\alpha_2,\beta_2)$. By (B),

$$\mathbf{U}\,\hat{\mathbf n}_i = \hat{\mathbf b}_H(\alpha_i,\beta_i), \qquad i = 1,2 .$$

Two non-parallel vector correspondences determine a rotation uniquely (Wahba's problem; solve by
Kabsch/SVD on the two triads $\{\hat{\mathbf n}_1, \hat{\mathbf n}_2, \hat{\mathbf n}_1\times\hat{\mathbf n}_2\}$
and $\{\hat{\mathbf b}_1,\hat{\mathbf b}_2,\hat{\mathbf b}_1\times\hat{\mathbf b}_2\}$).

**Mode B needs no $\varphi_D$ and no $\mathbf{F}$ whatsoever.** It uses only zone-axis identities and
stage readouts. This is a strong and, in the author's reading of the practical literature,
under-exploited result: the single hardest calibration constant in the problem is simply not
required if the operator has visited two zones — which, having chased Kikuchi bands to get here,
they almost always have.

Mode B also supplies a **calibration-free consistency test**:

$$
\Delta_{\text{consistency}} \;=\; \Big|\;\angle(\hat{\mathbf n}_1,\hat{\mathbf n}_2) \;-\; \angle\big(\hat{\mathbf b}_H(\alpha_1,\beta_1),\,\hat{\mathbf b}_H(\alpha_2,\beta_2)\big)\Big|.
\tag{C}
$$

The interzonal angle is a crystallographic invariant; the beam-direction angle depends only on the
stage model. A non-zero $\Delta_{\text{consistency}}$ therefore indicts, in order of likelihood: a
wrong zone-axis *sign* assignment, a wrong stage sign convention $s_\alpha$/$s_\beta$, a wrong
indexing, or a genuinely non-ideal stage. Testing (C) against all four sign combinations of
$(\pm\hat{\mathbf n}_1, \pm\hat{\mathbf n}_2)$ resolves the relative sense in most cases at zero
experimental cost.

**Residual ambiguity of Mode B.** Flipping *both* zone-axis senses admits a second proper rotation
$\mathbf{U}' = \mathbf{U}\mathbf{Q}$ where $\mathbf{Q}$ is the $180^\circ$ rotation about
$\hat{\mathbf n}_1 \times \hat{\mathbf n}_2$. Proof: $\mathbf{Q}\hat{\mathbf n}_i = -\hat{\mathbf n}_i$
for $i=1,2$ and $\mathbf{Q}$ fixes their cross product, so both triads map consistently and
$\det \mathbf{Q} = +1$.

This ambiguity is **harmless exactly when $\mathbf{Q}$ is an element of the crystal's proper point
group**, which is a one-line check the code must perform and report. For a cubic crystal with
$\hat{\mathbf n}_1 = [001]$, $\hat{\mathbf n}_2 = [110]$, $\mathbf{Q}$ is the two-fold about
$[\bar1 1 0]$, which *is* in $432$ — no ambiguity survives. The check is cheap and the result is
often "no ambiguity", which is a far better user experience than a blanket warning.

### 5.3 Mode C — one pattern plus a deliberate small tilt excursion

Take the current pattern, apply a known $\Delta\alpha$ (and separately $\Delta\beta$), record again,
and measure how the pattern moved. This determines $\varphi_D$ and $\mathbf{F}$ directly (§9.2) and
thereby upgrades Mode A to full determinacy. It costs two extra exposures and about a minute.

Mode C is the recommended bootstrap on an uncalibrated instrument; once run, its output is a
persistent `StageCalibration` record and Mode A becomes sufficient for the rest of the session
(subject to the diffraction-lens hysteresis caveat in §9.5).

### 5.4 Over-determined fitting

With $N \ge 2$ indexed zones the problem becomes least squares over $\mathbf{U}$ and, optionally,
over calibration parameters $\{\varphi_D, \varepsilon, \alpha_0, \beta_0, k\}$:

$$
\min_{\mathbf{U} \in SO(3),\;\boldsymbol\theta} \; \sum_{i=1}^{N} w_i \, \angle^2\!\Big(\mathbf{R}_{\text{stage}}(\alpha_i,\beta_i;\boldsymbol\theta)\,\mathbf{U}\,\hat{\mathbf n}_i,\;\hat{\mathbf z}_L\Big).
$$

Parameterize $\mathbf{U}$ in the exponential map about a current estimate to keep the problem
unconstrained; use Levenberg–Marquardt. Report the parameter covariance. **Do not fit calibration
parameters below $N = 5$** — the problem is under-determined and will produce confidently wrong
axis tilts that absorb indexing error.

---

## 6. Solving for the target tilts

### 6.1 Closed form (ideal stage)

Let $\hat{\mathbf t}_C$ be the target direction in crystal Cartesian coordinates and

$$\hat{\mathbf w} = \mathbf{U}\,\hat{\mathbf t}_C = (w_1, w_2, w_3)^{\mathsf T} \quad\text{(holder frame, known once } \mathbf{U} \text{ is)}.$$

We require $\mathbf{R}_{\text{stage}}(\alpha,\beta)\hat{\mathbf w} = \pm\hat{\mathbf z}_L$, i.e.
$\hat{\mathbf w} = \pm\hat{\mathbf b}_H(\alpha,\beta)$. Matching against (B):

$$
\boxed{\;\beta^{*} = \operatorname{atan2}(-w_1,\, w_3), \qquad
\alpha^{*} = \operatorname{atan2}\!\big(w_2,\; \sqrt{w_1^2 + w_3^2}\,\big)\;}
\tag{S}
$$

aligns $\hat{\mathbf w}$ with $+\hat{\mathbf z}_L$. Verification: with $\rho = \sqrt{w_1^2+w_3^2}$,
$\sin\beta^{*} = -w_1/\rho$ and $\cos\beta^{*} = w_3/\rho$, so
$\mathbf{R}_y(\beta^{*})^{\mathsf T}$ carries $\hat{\mathbf w}$ to $(0, w_2, \rho)$, and
$\mathbf{R}_x(\alpha^{*})^{\mathsf T}$ carries that to $(0,0,\sqrt{w_2^2+\rho^2}) = (0,0,1)$. $\square$

### 6.2 The four branches

The alignment condition has exactly four solutions on $(-180^\circ, 180^\circ]^2$:

| Branch | $(\alpha,\beta)$ | Lands on |
| --- | --- | --- |
| 1 | $\big(\operatorname{atan2}(w_2,\rho),\;\; \operatorname{atan2}(-w_1,w_3)\big)$ | $+\hat{\mathbf z}_L$ |
| 2 | $\big(\operatorname{atan2}(w_2,-\rho),\;\; \operatorname{atan2}(-w_1,w_3) \pm 180^\circ\big)$ | $+\hat{\mathbf z}_L$ |
| 3 | $\big(\operatorname{atan2}(-w_2,-\rho),\;\; \operatorname{atan2}(-w_1,w_3)\big)$ | $-\hat{\mathbf z}_L$ |
| 4 | $\big(\operatorname{atan2}(-w_2,\rho),\;\; \operatorname{atan2}(-w_1,w_3) \pm 180^\circ\big)$ | $-\hat{\mathbf z}_L$ |

Branches 2 and 4 require $|\beta| > 90^\circ$ and are never reachable on a real holder; the engine
enumerates them anyway so that the reachability report can say *why* nothing is available, rather
than returning an empty list. Branch 3 differs from branch 1 by roughly $180^\circ$ in $\alpha$ and
is likewise normally out of range — so in practice, **for each candidate target direction there is
exactly one relevant branch**, and the combinatorial richness of the answer comes entirely from
symmetry (§7) and ambiguity (§8), not from the trigonometry.

**Degeneracy.** If $\rho = 0$ — the target lies along $\pm\hat{\mathbf y}_H$, i.e. along the $\beta$
axis at $\alpha = 0$ — then $\beta$ is arbitrary and $\alpha = \pm 90^\circ$. The engine must detect
$\rho < \rho_{\text{tol}}$ and report a one-parameter family rather than an arbitrary $\operatorname{atan2}(0,0)$.
The configuration is unreachable on any real holder, but the numerics must not produce nonsense on
the way to saying so.

### 6.3 General stage: seed and refine

For the calibrated model of §3.3 there is no closed form. Use (S) as the seed and minimize

$$ f(\alpha,\beta) \;=\; \angle\Big(\mathbf{R}_{\text{stage}}(\alpha,\beta;\boldsymbol\theta)\,\mathbf{U}\,\hat{\mathbf t}_C,\;\; \hat{\mathbf z}_L\Big) $$

by Gauss–Newton with the analytic $2\times3$ Jacobian from (B). Two dimensions, a good seed, and a
smooth objective: convergence in three to five iterations to $10^{-10}$ degrees is expected. Because
the seed is analytic and the residual is a physical angle, a failure to converge is diagnostic (the
stage model is pathological) rather than a numerical accident — the engine should say so.

---

### 6.4 The complete algorithm, in pseudocode

Shown here, early, because every later section fills in one step of it. Forward references are
marked.

```
plan_tilt_to_zone_axis(current, target, stage, tolerance, ranking):

    # -- 1. current orientation ------------------------------------- (§5)
    U, cov_U, mode = reconstruct(current)          # Mode A / B / C
    assert det(U) == +1                            # parity self-check   (§5.1, §8.3)
    if mode is B:
        report interzonal_consistency(U, current)  # eq. (C)             (§5.2)

    # -- 2. candidate target directions ----------------------------- (§7, §8)
    G_proper = proper_rotation_group(phase)
    G_obs    = laue_rotation_group(phase) stabilizing current zone plane
    families = cosets of G_obs modulo (G_obs & G_proper)   # index <= 2  (§8.1)

    candidates = []
    for Q in families:                    # each is a distinct physical hypothesis
        for S in G_proper:                # each is a free, equivalent choice
            for sense in (+1, -1):
                candidates.append(sense * S @ Q @ t_hat_crystal)
    deduplicate on the unit sphere at 1e-8

    # -- 3. solve, per candidate ------------------------------------ (§6)
    solutions = []
    for t in candidates:
        w = U @ t
        for (a, b) in four_branches(w):            # eq. (S), table §6.2
            if stage.is_general:
                a, b = gauss_newton_refine(a, b, w, stage)     # (§6.3)

            # -- 4. reachability -------------------------------------- (§10)
            verdict, margin = stage.envelope.classify(a, b)
            if verdict is UNREACHABLE:
                record best-miss for the NEAREST_APPROACH report; continue

            # -- 5. forward validation -- mandatory, never skipped ---- (§14.1)
            v_lab    = stage.forward(a, b) @ U @ t
            residual = min(angle(v_lab, +z), angle(v_lab, -z))
            if residual > tolerance:
                reject                     # rejected BEFORE ranking, not down-weighted

            # -- 6. uncertainty and path ------------------------------ (§13, §11)
            sigma = propagate(cov_U, stage.calibration_cov, jacobian_of_S(w))
            path  = plan_path(current.stage_position, (a, b), stage, strategy)
            if not path.is_valid:          # limits, margin floor, excursions
                continue

            solutions.append(TiltSolution(a, b, verdict, residual, sigma, path,
                                          family=Q, orbit_member=t,
                                          kikuchi_band=rationalize(n_cur x t)))

    # -- 7. rank and report ----------------------------------------- (§13.4, §8.4)
    if solutions is empty:
        return nearest_approach_report(best_miss)              # (§10.4)
    rank(solutions, ranking)
    return TiltPlanReport(solutions,
                          ambiguity=describe_families(families),   # one entry per
                          discriminating_experiments=...)          # hypothesis (§8.5)
```

Two structural points the pseudocode is meant to make visible. First, **families and symmetry
operators enter the same loop but mean opposite things**: a symmetry operator offers a free
alternative the operator may take at will, a family is a competing hypothesis about reality of which
only one is true. They must be tagged differently in the output even though they are generated
together. Second, **forward validation sits inside the loop, before ranking** — a solution that has
not been re-derived through the calibrated forward model never reaches the user, so a solver bug
cannot present itself as an answer.

---

## 7. Symmetry

### 7.1 The target is an orbit, not a vector

The user asks for $[uvw]$; the crystal offers the whole orbit
$\{\,\pm \mathbf{S}_j\hat{\mathbf t}_C : \mathbf{S}_j \in G_{\text{proper}}\,\}$, every member of
which produces an identical diffraction pattern. For cubic $432$ that is up to 48 directions (24
operators $\times$ the $\pm$ sense). Enumerate via the existing `SymmetrySpec.equivalent_vectors`
with `antipodal=True`; deduplicate on the unit sphere at $10^{-8}$; then run §6 on each.

Two properties follow that are worth stating because they invert the intuition of an operator who
thinks in single poles:

1. **Reachability is a property of the orbit.** A target that looks hopeless can be routine because
   an equivalent sits $12^\circ$ away.
2. **The choice among reachable equivalents is free with respect to the pattern but not with respect
   to the specimen.** All equivalents give the same pattern; they give *different in-plane pattern
   rotations relative to specimen features*. The engine must therefore report, for each solution,
   the predicted azimuth of a nominated reference $\mathbf{g}$ in the pattern frame, so an operator
   comparing a boundary trace or a facet against the pattern can confirm which equivalent they got.

### 7.2 Which symmetry group

Use the **proper (rotation) subgroup** of the phase's point group for orbit generation: improper
operators do not correspond to physically attainable orientations of a chiral specimen. Use the
**Laue class rotation group** for the ambiguity analysis of §8.1 — and note carefully that these two
groups differ for exactly ten of the thirty-two point groups, tabulated in §8.1.
`pytex.core.point_groups` already exposes
`proper_subgroup_symbol_for` and `laue_class_symbol_for`, so both are available without new
crystallography.

### 7.3 Hexagonal and trigonal

Targets given as $[uvtw]$ must be converted through
`docs/standards/hexagonal_and_trigonal_conventions.md` before Cartesian conversion; the redundancy
constraint $u+v+t=0$ is checked, not assumed. Zone axes reported back to the user are rendered in
four-index form for hexagonal phases via `pytex.core.notation`, matching the rest of the repository.

---

## 8. Ambiguity — the heart of the problem

The prompt treats "the classical $180^\circ$ ambiguity" as one thing. It is three things with
different causes, different consequences, and different remedies. Conflating them is why this
failure mode has a reputation for being intractable. The decomposition below is the main analytical
claim of this document.

### 8.1 Layer 1 — crystallographic (Friedel / Laue)

**What one indexed SAED pattern determines uniquely.** The pattern is the intersection of the
weighted reciprocal lattice with a plane through the origin normal to the beam. From spot positions
alone one obtains: the **zone-axis line** (the unordered pair $\pm[uvw]$), the 2-D lattice of the
zone-plane section, and the assignment of $(hkl)$ to spots *up to the symmetry of the observation*.

**What is not determined.** Kinematic intensities obey Friedel's law $|F(\mathbf g)| = |F(-\mathbf g)|$,
so the recorded pattern is centrosymmetric **even when the crystal is not**. Consequently
$\mathbf g$ and $-\mathbf g$ are indistinguishable, and with them the sense of the zone axis.

**The exact statement.** Let $\mathbf{R}_{P\leftarrow C}$ be a reconstruction consistent with the
data. Any other consistent reconstruction is $\mathbf{R}_{P\leftarrow C}\mathbf{Q}$ with

$$
\mathbf{Q} \in G_{\text{obs}} = \Big\{\, \mathbf{Q} \in L \cap SO(3) \;:\; \mathbf{Q}\ \text{maps the zone plane to itself} \,\Big\},
$$

where $L$ is the **Laue class** of the phase. Two consequences follow that are not widely stated in
this form:

- **For a centrosymmetric phase, Friedel adds nothing.** $L \cap SO(3)$ equals the crystal's own
  proper point group, so every element of $G_{\text{obs}}$ is a genuine symmetry operation and every
  "ambiguous" reconstruction is a *symmetry-equivalent description of the same physical
  orientation*. For cubic $m\bar3m$ down $[001]$, $G_{\text{obs}} = 422 \subset 432$: the four-fold
  about $[001]$ and the two-folds about $\langle100\rangle$ and $\langle110\rangle$ in the zone
  plane. **All are crystal symmetries. The much-feared $180^\circ$ ambiguity is, at this layer,
  entirely benign for a centrosymmetric crystal.**
- **For exactly ten of the thirty-two point groups it is real.** $L \cap SO(3)$ strictly exceeds the
  proper point group precisely when the group contains improper operations *other than* inversion —
  e.g. point group $m$ has Laue class $2/m$, whose rotation subgroup $\{E, 2\}$ contains a two-fold
  the crystal does not have. Enumerating with `pytex.core.point_groups` gives the complete list, and
  the enlargement is a factor of **exactly two** in every case:

  | | Point groups | $|L \cap SO(3)| / |G_{\text{proper}}|$ |
  | --- | --- | --- |
  | Centrosymmetric (11) | $\bar1$, $2/m$, $mmm$, $4/m$, $4/mmm$, $\bar3$, $\bar3m$, $6/m$, $6/mmm$, $m\bar3$, $m\bar3m$ | 1 — no ambiguity |
  | Enantiomorphic (11) | $1$, $2$, $222$, $4$, $422$, $3$, $32$, $6$, $622$, $23$, $432$ | 1 — **no ambiguity** |
  | Improper, non-centrosymmetric (10) | $m$, $mm2$, $\bar4$, $4mm$, $\bar42m$, $3m$, $\bar6$, $6mm$, $\bar6m2$, $\bar43m$ | **2** |

  Note the middle row: a chiral crystal such as quartz ($32$) is non-centrosymmetric yet suffers
  **no** layer-1 ambiguity, because its point group is already all rotations. "Non-centrosymmetric"
  is therefore the wrong test, and a design that used it would warn on eleven groups without cause.

  The genuinely distinct solutions are counted by the cosets
  $G_{\text{obs}} / (G_{\text{obs}} \cap G_{\text{proper}})$, whose index is **at most two** and is
  often one even for the affected groups, because $G_{\text{obs}}$ for a *particular* zone may
  contain no offending operator. The engine computes the index per zone and reports it; when it is
  1, it says so, and does not warn.

**Consequence for tilting.** If the reconstruction is off by $\mathbf{Q} \in G_{\text{obs}}$, the
tilts computed for target $\mathbf t$ deliver the operator to $\mathbf{Q}\mathbf t$ instead. When
$\mathbf{Q}$ is a crystal symmetry, that is the same physical pole and nothing is lost. When it is
only a Laue symmetry, the destination is a *Friedel partner*: the SAED pattern there is
indistinguishable from the intended one, but the physical orientation differs, and every
polarity-sensitive measurement — CBED/HOLZ symmetry, polar-axis determination, defect
$\mathbf g\cdot\mathbf b$ analysis, trace analysis — differs with it. This is exactly the class of
error the prompt is right to fear, and it is confined to the ten affected point groups **and** to
polarity-sensitive purposes. Saying that precisely — rather than hoisting a universal warning
banner — is what lets the warning be believed when it does appear.

### 8.2 Layer 2 — instrumental rotation: the ambiguity that actually costs time

Mode A (§5.1) requires $\varphi_D$. Suppose it is wrong by $\delta\varphi$. Then
$\mathbf{U}_{\text{est}} = \mathbf{R}_z(\delta\varphi)\mathbf{U}_{\text{true}}$ (at the reference
stage position), and the operator, following the computed tilts, lands not on $\hat{\mathbf z}_L$
but at an angular error

$$
\boxed{\;\Delta_{\text{residual}} \;=\; 2\arcsin\!\big(\sin(\delta\varphi/2)\,\sin\theta_c\big) \;\approx\; \delta\varphi\,\sin\theta_c\;}
\tag{E}
$$

where $\theta_c$ is the angle between the target zone and the **current** zone. Derivation: with
$\mathbf c = \mathbf{U}_{\text{true}}\hat{\mathbf t}_C$, the tilts are chosen so that
$\mathbf{R}_{\text{stage}}\mathbf{R}_z(\delta\varphi)\mathbf c = \hat{\mathbf z}_L$, while the
specimen actually delivers $\mathbf{R}_{\text{stage}}\mathbf c$; the residual is
$\angle(\mathbf c, \mathbf{R}_z(\delta\varphi)\mathbf c)$, which is the stated expression. $\square$

Equation (E) is the most operationally useful result in this document:

- **The $\delta\varphi = 180^\circ$ case is catastrophic and self-consistent.** $\beta \to -\beta$,
  $\alpha \to -\alpha$: the operator tilts *exactly the wrong way in both axes*, the calculation
  reports a clean zero residual, and nothing appears. This — not Friedel — is the "$180^\circ$
  ambiguity" that wastes a session. It is an instrument-calibration ambiguity, it is **not absorbed
  by crystal symmetry**, and no amount of crystallographic care will detect it.
- **Error scales with the length of the hop.** A $5^\circ$ error in $\varphi_D$ costs $0.44^\circ$
  for a $5^\circ$ hop but $5^\circ$ for a $90^\circ$ hop. Therefore: **for long excursions, plan a
  multi-hop route through intermediate low-index zones and re-solve $\mathbf U$ at each** (§11.5).
  Each hop is short, so each is robust, and re-indexing at every waypoint makes the whole procedure
  self-correcting rather than open-loop. This is what a skilled operator does by instinct; here it
  falls out of the error model.

### 8.3 Layer 3 — parity and sign conventions

Three independent binary unknowns: the image mirroring $\mathbf F$, and the readout signs
$s_\alpha, s_\beta$. Each reflects the computed trajectory rather than rotating it, so unlike layer
2 they cannot be absorbed into any single angle. Detection:

- $\mathbf F$: enforce $\det \mathbf U = +1$ (§5.1); a mirrored recording makes the Mode A product
  improper and is caught for free.
- $s_\alpha, s_\beta$: caught by the Mode B consistency test (C), or directly by the tilt-excursion
  calibration of §9.2, which measures the sign as a by-product.

### 8.4 What the engine must and must not do

**Must not:** silently select one member of an ambiguity class and present it as determined.

**Must:**
1. Classify the ambiguity by layer and report the count of genuinely distinct solution families,
   after quotienting by crystal symmetry.
2. Emit **each** distinct family as a separate ranked solution with its own tilts and path, never a
   single answer with a footnote.
3. State, per family, the specific experimental observation that discriminates it — with the
   predicted outcome for each alternative, so the observation is a decisive test and not a hint.
4. Say "unambiguous" when it is, in as many words. Warning fatigue is itself a failure mode.

### 8.5 Discriminating observations, ranked by cost

| Observation | Resolves | Cost | Decisiveness |
| --- | --- | --- | --- |
| Known small $\alpha$ then $\beta$ excursion; measure Kikuchi motion azimuth (§9.2) | Layer 2 **and** Layer 3, completely | 2 exposures, ~1 min | **Decisive**; also yields $\varphi_D$ numerically |
| A second indexed zone axis (Mode B) | Layers 2 and 3 without any calibration | already available in most sessions | Decisive up to the §5.2 two-fold |
| Interzonal-angle test (C) across sign combinations | Relative zone-axis senses | free, uses existing data | Decisive when the candidate angles differ by more than the indexing error |
| CBED / HOLZ symmetry, or dynamical intensity asymmetry | Layer 1 for non-centrosymmetric phases | convergent probe, thin clean area | Decisive; the only route to absolute polarity |
| Calibrated image–diffraction rotation, applied to specimen morphology | Layer 2 | prior calibration | Decisive if the calibration is current (§9.5) |
| Known specimen feature (facet, growth direction, interface trace) | Layer 1 or 2 depending on the feature | free when present | Case-dependent |
| Ask the user which way the pattern moved | Layers 2 and 3 | one question | Decisive if the user answers accurately; the engine must state the two predicted outcomes so the question is answerable |

The tabulated ranking is deliberate: **the cheapest decisive test is the small-tilt excursion**, and
the design should make it the default first-run experience rather than an advanced option.

---

## 9. Calibration

### 9.1 What comes from metadata and what does not

| Quantity | Source |
| --- | --- |
| $\alpha, \beta$ readouts | **Metadata** (Velox/EMD, DM3/DM4, TIA `.ser`, JEOL) |
| Accelerating voltage, nominal camera length, magnification | **Metadata** |
| Detector pixel size, binning | **Metadata** |
| Diffraction rotation $\varphi_D(L_{\text{cam}}, V)$ | **Calibration.** Rarely stored, never reliably |
| Image rotation $\varphi_I(\text{mag})$, and $\varphi_{ID} = \varphi_I - \varphi_D$ | **Calibration** (classical MoO$_3$ / asbestos procedure) |
| Image parity $\mathbf F$ | **Calibration.** Depends on lens crossover count and on the acquisition software's array conventions |
| STEM scan rotation | Metadata *if* stored, else calibration; it is an independent additive term |
| $\alpha$-axis azimuth in the pattern frame | **Calibration** (§9.2) |
| $s_\alpha, s_\beta$, $\alpha_0, \beta_0$ | **Calibration** |
| Non-orthogonality $\varepsilon$, coupling $k$ | **Calibration** (§9.4); usually small enough to leave at zero with a stated residual budget |
| Backlash amplitude | **Calibration** (§9.6); not correctable, only avoidable |
| Tilt envelope | Holder datasheet + **experimental** verification |
| Eucentric height error | Metadata (stage $Z$) + operator procedure |

The honest summary: **the instrument tells you the two numbers you already knew, and nothing about
the transformation between them and the pattern.** Every quantity that matters must be measured.
This is the strongest argument for making Mode B the recommended path, since it needs none of them.

### 9.2 Primary procedure — the two-excursion tilt calibration

The core measurement. At a stage position where a Kikuchi pattern or a recognizable pattern feature
is visible:

1. Record the reference pattern at $(\alpha_c, \beta_c)$.
2. Apply $+\Delta\alpha$ (5–10°, large enough to exceed measurement noise, small enough to keep the
   same features in view). Record. Measure the in-pattern displacement direction of a tracked
   Kikuchi pole or band intersection; call its azimuth in the stored image $\psi_\alpha$.
3. Return, then apply $+\Delta\beta$. Record. Measure $\psi_\beta$.

**Predictions.** A rigid crystal rotation $\mathbf R$ carries a Kikuchi pole at the pattern centre
to $\mathbf R\hat{\mathbf z}_L$. At $(\alpha,\beta)=(0,0)$: $\mathbf R_x(\Delta\alpha)\hat{\mathbf z}_L
= (0, -\sin\Delta\alpha, \cos\Delta\alpha)$, a displacement along $-\hat{\mathbf y}_L$;
$\mathbf R_y(\Delta\beta)\hat{\mathbf z}_L = (\sin\Delta\beta, 0, \cos\Delta\beta)$, along
$+\hat{\mathbf x}_L$. Hence, with $\mathbf R_{L\leftarrow P} = \mathbf F \mathbf R_z(\varphi_D)$,

$$
\boxed{\;\varphi_D = -\psi_\beta, \qquad \psi_\alpha - \psi_\beta = -90^\circ \;\;(\mathbf F = \mathbf I) \quad\text{or}\quad +90^\circ \;\;(\text{mirrored})\;}
$$

Three results from two exposures: the diffraction rotation numerically, the parity from a sign, and
a redundancy check (the two azimuths must be $90^\circ$ apart — if they are not, the stage axes are
not orthogonal, or the tracked feature was misidentified, and the procedure has diagnosed itself).

The **magnitude** of the displacement is a bonus: a pole must move by exactly $\Delta\alpha$ in
angle, so the same two exposures calibrate the angular scale of the Kikuchi pattern and thereby
check the camera length.

### 9.3 Sign conventions

Run §9.2 with a *negative* excursion as well. $s_\alpha = +1$ if the displacement reverses;
$s_\alpha = -1$ if it does not behave as predicted for the positive case. Vendors differ, software
versions differ, and this must never be assumed from the manual.

### 9.4 Axis geometry from multiple observations

With $N \ge 5$ indexed zones at varied stage positions, fit
$\{\hat{\mathbf a}, \hat{\mathbf b}_0, \alpha_0, \beta_0, k\}$ jointly with $\mathbf U$ per §5.4.
Practical guidance: report the fitted $\varepsilon$ with its uncertainty and **do not adopt it unless
it is more than two standard errors from zero**. A spuriously fitted $1^\circ$ non-orthogonality is
worse than an assumed orthogonal stage, because it looks like knowledge.

### 9.5 Hysteresis in the lenses

$\varphi_D$ depends on the *history* of the projector and diffraction lenses, not only on the
nominal camera length. Mitigation is procedural, not computational: normalize the lenses, and always
approach a camera length from the same direction. The `StageCalibration` record must carry the
camera length and voltage it was measured at, and the engine must **refuse to apply a $\varphi_D$
measured at a different camera length** rather than interpolate. Interpolating a hysteretic
quantity manufactures a plausible wrong number, which is exactly the layer-2 failure of §8.2.

### 9.6 Backlash

Measure by approaching the same nominal stage position from opposite directions and recording the
difference in the resulting beam direction (via a re-indexed pattern). Typical magnitudes are a few
tenths of a degree, which matters when the target tolerance is $0.5^\circ$. Backlash is not
correctable in open loop; it is *avoided* by always approaching the final position from a consistent
direction (§11.4).

### 9.7 Eucentric height — a clarification

The prompt lists eucentric-height error among the calibration parameters. It belongs in a different
category and should be labelled as such: **to first order, being off eucentric height introduces no
orientation error at all.** Tilting about a displaced axis is a rotation composed with a
translation; the rotation — which is all that (M) depends on — is unchanged. The real cost is that
the region of interest translates out of the selected-area aperture during the tilt, so the operator
loses the crystal rather than mis-orients it. This is a serious practical failure mode and belongs
in the path planner (§11.6: bound the excursion between re-centring steps), not in the orientation
model. Stating this correctly saves an implementer from inventing a spurious height term.

### 9.8 STEM scan rotation

An additive, independently settable rotation between the scan raster and the detector. It affects
the *image* frame only; the diffraction pattern is unaffected. It must be recorded separately from
$\varphi_I$ and never folded into it, or a change of scan rotation will silently corrupt a
calibration that was correct.

---

## 10. Reachability

### 10.1 The accessible set, in closed form

From (B), as $(\alpha,\beta)$ ranges over the tilt envelope, the beam direction in holder
coordinates traces a region on the unit sphere bounded by:

- two **great circles** through $\pm\hat{\mathbf y}_H$ (constant $\beta$), and
- two **small circles** about $\hat{\mathbf y}_H$ (constant $\alpha$).

In crystal coordinates this region is $\mathbf U^{\mathsf T}$ applied to it. **Testing reachability
is therefore not a search**: for a candidate direction $\hat{\mathbf w} = \mathbf U\hat{\mathbf t}_C$,

$$
\alpha_{\text{req}} = \arcsin(w_2), \qquad \beta_{\text{req}} = \operatorname{atan2}(-w_1, w_3),
$$

and the target is reachable iff $(\alpha_{\text{req}}, \beta_{\text{req}})$ lies in the envelope —
plus the $-\hat{\mathbf w}$ test for the reversed sense. Constant time per candidate, exact.

### 10.2 How much of orientation space a holder actually reaches

$$\Omega = (\beta_{\max}-\beta_{\min})\,(\sin\alpha_{\max} - \sin\alpha_{\min}) \ \ \text{steradians.}$$

| Envelope | $\Omega$ (sr) | Fraction of $4\pi$ | Counting $\pm$ beam |
| --- | --- | --- | --- |
| $\pm30^\circ / \pm30^\circ$ | 1.047 | 8.3 % | 16.7 % |
| $\pm40^\circ / \pm30^\circ$ | 1.346 | 10.7 % | 21.4 % |
| $\pm70^\circ / 0$ (single tilt, tomography) | 0 (a curve) | 0 | 0 |

A single-tilt holder reaches a **set of measure zero** — a great circle — so an exact zone axis is
reachable only by coincidence, and the correct output there is "nearest approach", never "reachable".
This is a genuine qualitative difference the API must express, not a numerical edge case.

Taking a $\pm30^\circ/\pm30^\circ$ holder and treating orbit members as independently placed (a
crude approximation — real orbits are correlated by symmetry, and the engine computes the exact
answer rather than this estimate):

| Cubic family | Distinct lines | Rough $P(\text{some member reachable})$ |
| --- | --- | --- |
| $\langle 100\rangle$ | 3 | ~42 % |
| $\langle 111\rangle$ | 4 | ~52 % |
| $\langle 110\rangle$ | 6 | ~67 % |
| $\langle 112\rangle$ | 12 | ~89 % |

These figures justify the tool's existence in one line: **for a low-multiplicity target the answer
is often "no", and knowing that in five seconds instead of forty minutes is the deliverable.**

### 10.3 Envelope model

The envelope must be a **predicate**, not a box, because real holders reduce one range as the other
increases and because grid bars shadow at mounting-dependent azimuths:

```
TiltEnvelope.contains(alpha, beta) -> bool
TiltEnvelope.margin(alpha, beta)   -> float   # degrees to the nearest boundary
```

Concrete implementations: `RectangularEnvelope`, `EllipticalEnvelope`
($(\alpha/\alpha_m)^2 + (\beta/\beta_m)^2 \le 1$, a good fit to many double-tilt cartridges),
`PolygonEnvelope` (digitized from a datasheet or measured), and `MaskedEnvelope` (any of the above
minus operator-marked shadowing regions). `margin` is required, not optional — the ranking function
and the path planner both need a continuous distance-to-limit, not a boolean.

### 10.4 Four distinct verdicts

The API must distinguish, and the prompt is right to insist on this:

| Verdict | Meaning |
| --- | --- |
| `EXACT` | Forward-validated residual below `tolerance_deg`, comfortably inside the envelope |
| `WITHIN_TOLERANCE` | Reachable, residual within tolerance but not negligible — states the residual and its dominant cause |
| `NEAREST_APPROACH` | No orbit member is reachable; returns the envelope point minimizing angle to the nearest member, **with that angle**, so the operator can judge whether a partial tilt is useful (often it is — a $6^\circ$ miss still shows the Kikuchi pole) |
| `UNREACHABLE` | No orbit member within a useful distance; states the best miss angle and, where applicable, that re-mounting or a different holder would be required |

---

## 11. Path planning

### 11.1 Why the endpoint is not the answer

A straight line in $(\alpha,\beta)$ space is *not* a straight line on the crystal sphere: by (B) the
$\beta$ contribution is scaled by $\cos\alpha$ and the two rotations do not commute. Different
interpolations produce genuinely different trajectories through orientation space, and they differ
in whether the operator can follow along.

### 11.2 The three candidate paths

| Path | Definition | Property |
| --- | --- | --- |
| **Sequential** ($\beta$ then $\alpha$, or $\alpha$ then $\beta$) | one axis at a time | Traces a great circle through $\pm\hat{\mathbf y}_H$ then a small circle. Easiest to execute by hand; longest travel |
| **Linear in $(\alpha,\beta)$** | $(1-s)(\alpha_c,\beta_c) + s(\alpha_t,\beta_t)$ | Simple; a rhumb-line-like curve on the sphere; no special virtue except simplicity |
| **Geodesic in beam space** | shortest great-circle arc from $\hat{\mathbf b}_C(\text{current})$ to $\hat{\mathbf b}_C(\text{target})$, inverted through (B) at each sample | Minimum crystal travel |

### 11.3 The geodesic path is the Kikuchi band — the practice and the mathematics agree

The geodesic between two zone axes on the crystal sphere lies in the plane spanned by them. That
plane has normal $\hat{\mathbf n}_1 \times \hat{\mathbf n}_2$, which rationalizes to a low-index
$(hkl)$ — and the Kikuchi band of that reflection is precisely the band connecting the two poles.

**So the mathematically optimal path is the one experienced operators already follow by eye.** This
is a strong endorsement of the formulation and it yields the single most useful user-facing output
in the whole system:

> *"Follow the $(1\bar10)$ Kikuchi band from $[001]$ toward $[111]$; total travel $54.7^\circ$."*

`pytex.diffraction.composite.rationalize_zone_axis` already performs the rationalization, applied
here to the cross product rather than to a zone axis. The default path planner should be the
geodesic, and the reported band should be part of the standard output.

### 11.4 Backlash-aware approach

Terminate every path with a consistent approach direction: overshoot the final position by
$\delta_{\text{approach}}$ (default $2^\circ$, from the §9.6 measurement) on each axis in a fixed
sense and return. Adds one leg, removes the dominant repeatability error. When the natural approach
already comes from the preferred side, the planner omits the overshoot rather than adding a
pointless excursion.

### 11.5 Waypoints and multi-hop routing

For $\theta_c$ beyond a threshold (default $30^\circ$, from equation (E) and the calibration
uncertainty budget), route through intermediate **low-index zone axes** that lie near the geodesic
and inside the envelope. At each waypoint the operator re-indexes and the engine re-solves
$\mathbf U$, so accumulated calibration error is reset rather than compounded. Waypoint selection:
enumerate zone axes up to a maximum index, score by (perpendicular distance from the geodesic,
envelope margin, index lowness — low-index zones are easier to recognize and index), take a small
ordered subset.

This converts an open-loop calculation into a closed-loop procedure, and it is the recommended mode
for any excursion where the calibration is not fresh.

### 11.6 Constraints the planner must enforce

- Every sampled point inside the envelope, with margin above a floor (default $2^\circ$) — a path
  that grazes a mechanical limit will hit it once backlash is included.
- No excursion beyond the bounding box of the endpoints unless a waypoint requires it.
- Segment length bounded so the operator can re-centre the eucentric height and the area of
  interest (§9.7) — default a re-centring prompt every $15^\circ$.
- Simultaneous $\alpha$–$\beta$ motion permitted only when the stage supports it; otherwise the
  planner emits a sequence of single-axis legs, because that is what the operator will actually do
  and the plan must match reality.

### 11.7 Path output

At every sample: $(\alpha,\beta)$, the beam direction in crystal coordinates $\hat{\mathbf b}_C$,
its nearest low-index zone axis and the angle to it (so the operator knows what they should be
seeing), the envelope margin, and the cumulative travel. This is the table the operator works from
and the array the visualization animates — the same array, which is what makes §12 a cross-check
rather than an illustration.

---

## 12. Stereographic projection and animation

### 12.1 Formulation

Work in the **crystal frame** so the stereogram is a crystallographic object the user can compare
with a standard projection. The beam direction in crystal coordinates along the path is

$$
\hat{\mathbf b}_C(\alpha,\beta) = \mathbf U^{\mathsf T}\,\hat{\mathbf b}_H(\alpha,\beta)
= \mathbf U^{\mathsf T}\big(-\cos\alpha\sin\beta,\;\sin\alpha,\;\cos\alpha\cos\beta\big)^{\mathsf T},
$$

projected by the equal-angle (stereographic) rule from the south pole,
$(X, Y) = (b_x, b_y)/(1 + b_z)$ for $b_z > 0$, with lower-hemisphere points marked distinctly rather
than folded silently. `pytex.diffraction.stereonets` already provides the projection machinery,
great- and small-circle sampling, and the projection-boundary radius; `pytex.plotting.spherical`
provides the rendering. **No new projection mathematics is required** — which is the point: the
visualization must consume the same functions as the engine.

### 12.2 What the accessible region looks like, and why that is a good thing

From §10.1 the reachable region is bounded by two great circles and two small circles about
$\hat{\mathbf y}_H$. Both curve families are already first-class in `stereonets`
(`sample_great_circle`, `sample_small_circle`), so the reachable region draws exactly, as curves,
with no polygonal approximation and no sampling artefacts. An operator can then *see* whether a pole
falls inside — which is the whole question of §10 rendered as a picture.

### 12.3 Layers

1. Standard stereographic net and the phase's fundamental sector boundary
   (`SymmetrySpec.fundamental_sector`).
2. Low-index poles labelled in the repository's notation, with zone circles for major families.
3. Current beam direction (filled marker) and the requested target (open marker).
4. Every symmetry-equivalent target, distinguished as reachable / unreachable.
5. The accessible region, as the exact curve boundary of §12.2, with a margin contour.
6. The planned trajectory, sampled from the **same array** the planner emitted, with
   $(\alpha,\beta)$ annotated at intervals.
7. Alternative trajectories for each distinct ambiguity family (§8.4), visually distinct, so a
   $180^\circ$ layer-2 error is *seen* as a trajectory heading the opposite way — the single most
   valuable diagnostic in the system.
8. The final residual, drawn to scale as a small circle about the target.

### 12.4 The independence requirement

The prompt asks the animation to be an independent cross-check. Two design rules make this real
rather than aspirational:

- **The renderer consumes `TiltPath` samples and does no kinematics of its own.** It may not
  recompute $\hat{\mathbf b}$ from $(\alpha,\beta)$ by a second code path; it plots what the engine
  produced. A drawing that agrees because it re-implements the same formula proves nothing.
- **Independence is obtained in the test suite, not the renderer**: a validation test recomputes the
  plotted trajectory by an independent route — accumulating small finite rotations by matrix
  exponentials rather than by the closed-form (B) — and asserts agreement to $10^{-9}$. That is a
  real cross-check; a second copy of the same equation in the plotting layer is not.

### 12.5 Animation mechanics — an honest blocker

`pytex.plotting` currently has **no animation infrastructure**: it is Matplotlib figures and
hand-built SVG, all deterministic and headless-safe. Adding animation risks importing
non-determinism and a heavy dependency into a repository whose figure tests are exact.

Recommendation: make the **frame sequence the primitive**. `TiltPath.sample(n)` yields a
deterministic list of states; `render_tilt_stereogram(state)` renders one frame. Then:

- static multi-panel contact sheet (default; testable exactly, embeddable in the docs);
- animated SVG via SMIL — self-contained, no JavaScript, renders in the Sphinx site, and diffable as
  text, which fits the repository's existing SVG practice;
- GIF via `matplotlib.animation` behind an optional extra, never a required dependency.

Test the trajectory and the per-frame geometry; never test pixels.

### 12.6 The teaching mode the prompt asks for

Expose an interactive/parametric mode that sweeps $\alpha$ alone and $\beta$ alone from the current
position and animates $\hat{\mathbf b}_C$. This answers "what does positive $\alpha$ actually do to
my pattern" directly and visually, and — because it uses the calibrated model — it doubles as a
qualitative check on the calibration itself before any real tilt is attempted.

---

## 13. Error propagation

### 13.1 Sources and typical magnitudes

| Source | Typical | Enters as |
| --- | --- | --- |
| Indexing residual $\to$ $\delta\mathbf U$ | $0.3$–$1.0^\circ$ | Rotation covariance on $\mathbf U$, strongly anisotropic: well constrained transverse to the beam, poorly about it |
| $\delta\varphi_D$ | $2$–$5^\circ$ uncalibrated; $<1^\circ$ after §9.2 | Rotation about $\hat{\mathbf z}_L$; residual $\approx \delta\varphi_D\sin\theta_c$ (E) |
| Stage readout | $0.1$–$0.5^\circ$ | Directly on $(\alpha,\beta)$ |
| Backlash | $0.1$–$0.5^\circ$ | Directly, sign-dependent; removed by §11.4 |
| Non-orthogonality $\varepsilon$ | $0$–$1^\circ$ | Second order at small tilt, first order at large |
| Specimen bending | $0.1$–$2^\circ$ | Not modellable; the dominant term in a buckled foil |

Note the anisotropy of the first row: **a single pattern constrains the orientation about the beam
axis much more weakly than transverse to it**, because the in-plane constraint comes from spot
azimuths on a small pattern. A scalar "orientation uncertainty" would misrepresent this by a factor
of several, which is why §4 specifies a full tangent-space covariance.

### 13.2 Conditioning of the inverse

Differentiating (B):

$$
\delta\alpha = \frac{\delta w_2}{\cos\alpha}, \qquad
\delta\beta \approx \frac{\delta w_\perp}{\cos\alpha}.
$$

Both blow up as $|\alpha| \to 90^\circ$ — the gimbal lock of §3.5. Within a real envelope
($|\alpha| \le 40^\circ$) the amplification is at most $1.31$, which is benign; the engine should
nonetheless report the factor $1/\cos\alpha$ per solution and prefer, all else equal, solutions at
smaller $|\alpha|$. Stating the bound is more useful than a warning: it tells a reviewer the
parameterization is safe over the operating range, and exactly where it stops being safe.

### 13.3 Propagation to a reported uncertainty

Two routes, both implemented:

- **Linear:** propagate the tangent-space covariance of $\mathbf U$ and the calibration covariance
  through the analytic Jacobian of (S) to a $2\times2$ covariance on $(\alpha,\beta)$. Cheap,
  adequate away from the singularity, and the default.
- **Monte Carlo:** sample $\mathbf U$ and the calibration parameters, re-solve, and report the
  empirical spread of $(\alpha,\beta)$ and of the forward residual. Used for the reported
  robustness score, and for the §14.6 statistical validation.

Every solution reports $\sigma_\alpha$, $\sigma_\beta$, and $\sigma_{\text{residual}}$. A solution
whose $\sigma_{\text{residual}}$ exceeds the requested tolerance is flagged: *the tilts are right,
but the inputs are not good enough to promise the target.* That distinction — between an imprecise
calculation and imprecise data — is one the operator needs and is rarely given.

### 13.4 Ranking

$$
\text{cost} = w_1 \Theta_{\text{travel}} + w_2 \max(|\Delta\alpha|,|\Delta\beta|) + w_3\,\Phi_{\text{margin}} + w_4\,\sigma_{\text{residual}} + w_5\,(1/\cos\alpha_t - 1) + w_6\,\Phi_{\text{path}}
$$

with $\Theta_{\text{travel}}$ the actual specimen rotation angle
$\angle\big(\mathbf R_{\text{stage}}(\alpha_t,\beta_t)\mathbf R_{\text{stage}}(\alpha_c,\beta_c)^{\mathsf T}\big)$
— not the naive $\sqrt{\Delta\alpha^2+\Delta\beta^2}$, which is not an angle of anything — and
$\Phi_{\text{margin}} = \exp(-m/m_0)$ a soft barrier on the envelope margin. Weights are
user-overridable and their defaults must be documented with reasons, per the repository's tolerance
governance standard. Solutions that fail forward validation (§14.1) are **rejected before ranking**,
never merely down-weighted.

---

## 14. Validation

### 14.1 Mandatory forward validation of every solution

Non-negotiable, and cheap:

1. Take the estimated $\mathbf U$.
2. Build $\mathbf R_{\text{stage}}(\alpha_t, \beta_t)$ from the **calibrated** model — the same
   object the solver used, so this catches solver error, not model error.
3. Compute $\mathbf v_L = \mathbf R_{\text{stage}}\mathbf U \hat{\mathbf t}_C$.
4. Residual $= \min\big(\angle(\mathbf v_L, \hat{\mathbf z}_L), \angle(\mathbf v_L, -\hat{\mathbf z}_L)\big)$
   unless a sense is required and determined.
5. Reject above `tolerance_deg`.

A solution that has not passed step 5 must not be returned in the ranked list at all.

### 14.2 Test matrix

| Class | Cases | Assertion |
| --- | --- | --- |
| Round trip | Random $\mathbf U$, random reachable target, ideal stage | Residual $< 10^{-9}$ deg |
| Forward/reverse | Plan A$\to$B, then B$\to$A from the resulting state | Returns to the original $(\alpha,\beta)$ within $10^{-9}$ |
| Known cubic pairs | $[001]\!\to\![011]$ $45^\circ$; $[001]\!\to\![111]$ $54.736^\circ$; $[001]\!\to\![112]$ $35.264^\circ$; $[011]\!\to\![111]$ $35.264^\circ$; $[111]\!\to\![\bar111]$ $70.529^\circ$ | Computed $\Theta_{\text{travel}}$ matches the interzonal angle when the geodesic path is used and both poles are reachable |
| Hexagonal | $[0001]\!\to\!\langle11\bar20\rangle$ exactly $90^\circ$ by symmetry; $c/a$-dependent pairs against the metric tensor | Matches `CrystalDirection.angle_to` |
| Tetragonal, monoclinic, triclinic | Low-symmetry orbits, small stabilizers | Orbit size matches the point-group order; no spurious equivalents |
| Symmetry orbit | Every orbit member solved independently | All reachable members produce identical *patterns* under forward simulation; in-plane azimuths differ as predicted |
| Near limits | Targets at margin $0.1^\circ$ inside and outside | Verdict flips at the boundary; no solution is returned with negative margin |
| Single-tilt holder | Any off-circle target | Verdict is `NEAREST_APPROACH`, never `EXACT` |
| Gimbal | Target at $\pm\hat{\mathbf y}_H$ ($\rho \to 0$) | One-parameter family reported; no NaN; no arbitrary $\beta$ |
| Ambiguity — centrosymmetric | $m\bar3m$, $6/mmm$ | Engine reports **one** family; no spurious warning |
| Ambiguity — enantiomorphic | $32$ (quartz), $422$, $23$ | **One** family — non-centrosymmetric but ambiguity-free (§8.1); asserts the engine does not warn on chirality |
| Ambiguity — improper non-centrosymmetric | $\bar43m$ (GaAs), $6mm$ (ZnO), $4mm$, $mm2$ | Coset index 2 where the zone admits it; distinct families emitted with discriminating experiments |
| Bad camera rotation | Inject $\delta\varphi_D \in \{5,30,180\}^\circ$ | Measured residual matches (E) to within $10^{-6}$ — this validates the *error model*, not just the code |
| Parity error | Mirrored $\mathbf F$ | Mode A raises on $\det\mathbf U = -1$; Mode B consistency test (C) fires |
| Non-orthogonal axes | $\varepsilon \in \{1,3,5\}^\circ$ | Seed-and-refine converges; ignoring $\varepsilon$ produces the predicted residual, and the prediction is asserted |
| Noisy indexing | Perturb spot positions at realistic amplitude, 1000 trials | Reported $\sigma_{\text{residual}}$ covers the actual error at $68 \pm 3$ % — a calibration of the uncertainty estimate itself, not merely of the answer |

The last row is the one that makes the uncertainty reporting trustworthy, and it is the row most
often omitted from work of this kind.

### 14.3 Independent cross-checks

- **Hand construction.** Reproduce two or three classical Wulff-net tilt constructions for cubic
  from the TEM literature and assert agreement to the reading precision of the net ($\sim1^\circ$).
- **Existing engine.** Forward-simulate the destination pattern with
  `pytex.diffraction.saed.generate_saed_pattern` and confirm the zone axis it reports is the
  intended one — a genuinely independent path through the codebase.
- **Matrix-exponential integration.** Per §12.4, integrate the trajectory by accumulating
  infinitesimal rotations and compare with the closed form.
- **Rigid-crystal assumption.** Using a measured multi-zone dataset, report the Mode B consistency
  residual (C) across all pairs; a systematic non-zero value bounds the specimen bending, which is
  the assumption's failure mode, quantified rather than assumed away.

---

## 15. Inputs and outputs

### 15.1 Inputs

```python
plan_tilt_to_zone_axis(
    current: CurrentState,            # U (or the data to reconstruct it) + (alpha_c, beta_c)
    target: ZoneAxis | CrystalDirection,
    stage: StageModel,                # kinematics + calibration + envelope
    *,
    tolerance_deg: float = 0.5,
    require_sense: bool = False,      # honoured only if the sense is determined
    path: PathStrategy = GEODESIC,
    max_solutions: int = 5,
    ranking: RankingWeights = DEFAULT_RANKING,
) -> TiltPlanReport
```

`CurrentState` is constructible three ways, matching §5: from a `PatternSolution` plus
$(\alpha_c,\beta_c)$ plus a `StageCalibration` (Mode A); from two indexed zones at two stage
positions (Mode B, no calibration needed); or from an externally known `Orientation` (synthetic and
test use). The constructor records **which** mode was used, and the ambiguity analysis reads that
record — the mode is not re-inferred downstream.

### 15.2 Outputs

`TiltPlanReport` carries: the reconstruction mode and its residual; the ambiguity classification by
layer with the count of distinct families and, per family, the discriminating experiment and its
predicted outcomes; the ranked `TiltSolution` list; and, per solution — target tilts, verdict
(§10.4), forward-validated residual, $\sigma$ values, travel, envelope margin, conditioning factor,
which orbit member it lands on, the predicted in-plane azimuth of a reference $\mathbf g$, the
`TiltPath` with per-sample state, and the connecting Kikuchi band $(hkl)$.

Per repository convention every one of these objects implements `describe()` returning prose an
operator can read at the microscope, and `to_json_dict()` against a registered schema
(`pytex.tilt_plan_report/1`) for the manifest system.

---

## 16. Proposed software architecture

New subpackage `src/pytex/tem/` — a new domain of use (instrument operation) rather than a new
crystallography, so it is a peer of `diffraction` and `ebsd`, not a module inside either.

| Module | Contents |
| --- | --- |
| `tem/stage.py` | `StageModel` protocol; `DoubleTiltStage`, `TiltRotateStage`, `SingleTiltStage`; `StageCalibration`; `TiltEnvelope` hierarchy; forward kinematics (M) and (B) |
| `tem/reconstruction.py` | `CurrentState`; Modes A/B/C; over-determined fit; covariance |
| `tem/ambiguity.py` | $G_{\text{obs}}$, the Laue/proper coset count, family enumeration, discriminating-experiment catalogue |
| `tem/navigation.py` | Closed form (S), branch enumeration, symmetry orbit, refinement, forward validation, ranking, `plan_tilt_to_zone_axis` |
| `tem/path.py` | `PathStrategy` implementations, waypoint search, backlash approach, envelope checking, `TiltPath` |
| `tem/calibration.py` | The §9.2 two-excursion procedure as a fitting routine; §9.4 multi-zone fit; YAML contract for a persisted calibration |
| `plotting/tilt_stereogram.py` | Single-frame renderer, contact sheet, SMIL animation; consumes `TiltPath` only |
| `contracts.py` (extended) | `pytex.stage_calibration/1`, `pytex.tilt_plan_report/1` |

Dependency direction: `tem` depends on `core` and `diffraction`; **nothing depends on `tem`**. The
plotting module depends on `tem` data classes but no `tem` module imports plotting — matching the
existing layering.

### 16.1 Suggested phasing

| Phase | Deliverable | Gate |
| --- | --- | --- |
| TN0 | This document, reviewed and approved | Scientific sign-off |
| TN1 | `stage.py` + envelopes + forward kinematics | §14.2 round-trip and gimbal rows pass |
| TN2 | `navigation.py`: closed form, branches, symmetry orbit, forward validation | Known-pair and orbit rows pass |
| TN3 | `reconstruction.py` Modes A/B/C + `ambiguity.py` | Centrosymmetric/non-centrosymmetric ambiguity rows pass |
| TN4 | `path.py` incl. geodesic/Kikuchi-band reporting and waypoints | Forward/reverse and limit rows pass |
| TN5 | `calibration.py` + uncertainty (§13) | Noisy-indexing coverage row passes at $68\pm3$ % |
| TN6 | Visualization + worked notebook + docs | Independent-integration cross-check passes |

TN1–TN2 alone deliver a useful tool for a synthetically known orientation, which makes the first
increment independently valuable — worth preserving through review.

---

## 17. Limitations, blockers, and honest failure modes

1. **Calibration is the binding constraint, not mathematics.** The geometry is exact; $\varphi_D$
   and the parity bits are not, and they are not in the metadata. Mode B exists precisely to route
   around this and should be the documented default.
2. **Diffraction-lens hysteresis makes $\varphi_D$ history-dependent.** Mitigated procedurally
   (§9.5); the engine refuses to extrapolate rather than interpolating a hysteretic quantity.
3. **Specimen bending breaks the rigid-crystal assumption** and is not modellable from a single
   pattern. Quantified by the Mode B consistency residual (§14.3), not assumed away.
4. **No animation infrastructure exists** in `pytex.plotting` (§12.5). Frame-sequence primitive plus
   SMIL SVG is the recommended answer; a heavyweight animation dependency should be resisted.
5. **Absolute polarity for the ten affected point groups is out of reach** from kinematic SAED alone
   (§8.1). CBED/HOLZ is the only route, and it is out of scope; the engine's duty is to say so, not
   to guess.
6. **Off-eucentric operation loses the region, not the orientation** (§9.7). Handled by path
   segmentation; there is no orientation-model term for it, and inventing one would be wrong.
7. **Backlash is avoided, not corrected** (§9.6, §11.4).
8. **No instrument I/O.** The output is a plan; execution is the operator's or an external layer's.
9. **Single crystallite assumed.** A selected area containing more than one grain must be indexed
   per grain upstream; TN does not disentangle overlapping patterns.
10. **The rough reachability probabilities of §10.2 assume independent orbit placement** and are
    illustrative only. The engine computes the exact answer per case; the table must never be quoted
    as a result.

---

## 18. Why the method works, and when it does not

**Why it works.** The whole problem is one rigid-body equation, (M). Every input is a constraint on
one unknown rotation $\mathbf U$ plus a small set of instrument constants. The alignment condition
inverts in closed form because the double-tilt holder's kinematics reduce, by the exact cancellation
of §3.1, to a spherical coordinate system (B) whose pole is the $\beta$ axis. Symmetry, reachability,
path planning and visualization are then all statements about one region on one sphere, computed
with machinery the repository already has.

**When it fails.** In order of practical likelihood:

1. $\varphi_D$ wrong (worst at $180^\circ$): a clean, confident, exactly-backwards answer. Detected
   only by §9.2 calibration or by Mode B. **This is the failure mode the design must be organized
   around**, and it is the reason Mode B is recommended over the more obvious Mode A.
2. Parity/sign convention wrong: reflected trajectory. Detected by $\det\mathbf U$ and by (C).
3. Specimen bent or the area polycrystalline: the rigid assumption fails and no reconstruction is
   valid. Detected by (C) across pairs.
4. Target orbit genuinely outside the envelope: not a failure but an answer, and the tool's most
   common useful output.
5. Target near the $\beta$-axis pole: ill-conditioned by $1/\cos\alpha$, bounded at $1.31$ within a
   real envelope, reported per solution.
6. One of the ten improper non-centrosymmetric point groups combined with a polarity-sensitive
   purpose: irreducible from SAED; reported, not guessed.

---

## 19. Verification status of this document

The derivations here were checked numerically before the document was submitted for review; the
results are recorded so a reviewer knows which statements are asserted and which are demonstrated.

| Claim | §  | Check | Result |
| --- | --- | --- | --- |
| Moving-$\beta$-axis composition equals $\mathbf R_x(\alpha)\mathbf R_y(\beta)$ | 3.1 | Direct matrix comparison | Exact |
| Beam-direction formula (B) | 3.5 | Against $\mathbf R_{\text{stage}}^{\mathsf T}\hat{\mathbf z}$ | Exact |
| Closed form (S) | 6.1 | 20 000 random directions | Max error $4.6\times10^{-16}$ |
| All four branches land where the table says | 6.2 | Direct evaluation | All four confirmed |
| $\lVert\partial\hat{\mathbf b}/\partial\alpha\rVert = 1$, $\lVert\partial\hat{\mathbf b}/\partial\beta\rVert = \cos\alpha$ | 3.5 | Finite differences at $\alpha = 0, 23, 40^\circ$ | Matches to 6 digits |
| Residual law (E), $\Delta = 2\arcsin(\sin(\delta\varphi/2)\sin\theta_c)$ | 8.2 | 3000 random cases at $\delta\varphi = 5, 30, 180^\circ$ | Max deviation $1.3\times10^{-10}$ deg |
| $\delta\varphi = 180^\circ$ negates **both** $\alpha$ and $\beta$ | 8.2 | Worked example | $(-10.32, -52.65) \to (+10.32, +52.65)$ |
| Mode B residual ambiguity is the two-fold about $\hat{\mathbf n}_1\times\hat{\mathbf n}_2$ | 5.2 | Constructed $\mathbf Q$; checked action and $\det$ | Confirmed, $\det = +1$ |
| Solid angles and coverage percentages | 10.2 | Direct integration | $1.0472$ sr / 8.33 %; $1.3463$ sr / 10.71 % |
| Interzonal angles in the test matrix | 14.2 | Metric computation | $45.000$, $54.736$, $35.264$, $35.264$, $70.529$ |
| Laue-vs-proper enlargement over all 32 point groups | 8.1 | `pytex.core.point_groups` enumeration | Factor exactly 2 for the 10 listed groups, 1 for the other 22 |
| $|G_{\text{obs}}| = 8$ for cubic $m\bar3m$ down $[001]$ | 8.1 | Stabilizer enumeration | Confirmed (the group $422$) |

The point-group enumeration corrected a drafting error: quartz ($32$) was initially listed as a
layer-1 ambiguity case. It is enantiomorphic, so its proper group already equals its Laue rotation
group and it has **no** such ambiguity. The corrected criterion — improper operations other than
inversion, not merely absence of a centre — is now stated in §8.1 and tested for in §14.2.

---

## 20. References

### Normative

- Hahn, Th. (ed.), *International Tables for Crystallography, Volume A*, IUCr / Springer, DOI: <https://doi.org/10.1107/97809553602060000100>.
- `docs/standards/notation_and_conventions.md` — frame domains, handedness, quaternion order.
- `docs/architecture/reference_frame_foundation.md` — the frame and transform model TN builds on.
- `docs/standards/benchmark_and_tolerance_governance.md` — the standard every default tolerance in §13.4 must satisfy.

### Informative

- De Graef, M., *Introduction to Conventional Transmission Electron Microscopy*, Cambridge University Press, DOI: <https://doi.org/10.1017/CBO9780511615092> — stage geometry, diffraction rotation, Kikuchi navigation.
- Williams, D. B. and Carter, C. B., *Transmission Electron Microscopy*, Springer, DOI: <https://doi.org/10.1007/978-0-387-76501-3> — indexing, tilt practice, calibration procedures.
- Edington, J. W., *Practical Electron Microscopy in Materials Science*, Macmillan — the classical ratio/angle indexing the existing solver implements.
- Wahba, G., *SIAM Review* 7 (1965) 409, DOI: <https://doi.org/10.1137/1007077> — the two-vector attitude problem underlying Mode B.
- Morawiec, A., *Orientations and Rotations*, Springer, DOI: <https://doi.org/10.1007/978-3-662-09156-2> — rotation representations, symmetry orbits, covariance in the tangent space.
