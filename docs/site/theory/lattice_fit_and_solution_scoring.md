# Fitting The Pattern Lattice, And Scoring The Solutions

Two steps sit either side of indexing a selected-area pattern, and neither is
crystallography. Before: the picked spots lie on a plane lattice, which
over-determines the transmitted beam and exposes any pick that does not belong.
After: several candidate indexings need to be compared, which requires saying
what disagrees, by how much, and how much each kind of disagreement counts.

The indexing between them is documented in
[Ratio/angle indexing of a measured SAED pattern](saed_ratio_angle_indexing.md);
the workflow that strings all three together is
[Indexing a pattern and choosing where to tilt next](../workflows/tem_pattern_indexing.md).

## 1. The plane lattice of a zone-axis pattern

A zone-axis pattern is a plane through the reciprocal lattice, so every spot sits
at

$$ \mathbf{p}_i = \mathbf{c} + m_i\,\mathbf{a} + n_i\,\mathbf{b} $$

with $m_i, n_i$ integers, $\mathbf{a}, \mathbf{b}$ two basis vectors the crystal
fixes, and $\mathbf{c}$ the transmitted beam. Nothing here needs a phase, a
camera constant, or a zone axis. Only the picks.

### Why the centre is worth solving for

With $N$ picked spots the model has $2N$ observations and $6$ parameters — two
each for $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$ — once the integers are known.
Four spots make it over-determined, and the least-squares centre then uses every
spot rather than the single click that placed the beam.

This matters more than the arithmetic suggests. The camera equation
$r = L\lambda / d$ measures $r$ *from the beam*, so an error in $\mathbf{c}$
biases every spacing in the pattern at once, in the same direction. The result
is not a noisy answer but a **self-consistent wrong one**: a plausible lattice
parameter for a material that is not there.

### Estimating it

With the integers fixed the model is linear in every fitted parameter, so one
`lstsq` on the design matrix $[\,\mathbf{1}\;\;m\;\;n\,]$ solves for all three
vectors at once, $x$ and $y$ sharing the design and differing only in the
right-hand side. There is no optimizer and no starting-value sensitivity. The
integers come from rounding $B^{-1}(\mathbf{p}_i - \mathbf{c})$, and the two
steps alternate until the integers stop moving — usually two or three rounds.

### Three failure modes, and what answers them

**Seeding from the centre subdivides the cell.** The trial basis must come from
*differences between spots*, never from offsets to the picked centre. A
difference of two spots is a lattice vector however badly $\mathbf{c}$ was
picked; an offset is one only if the pick was already right. Seed from offsets
and a centre half a spacing out makes some of them spuriously short, the shortest
pair generates a sub-lattice, and the fit explains every spot perfectly by
halving the cell around the wrong centre — machine-precision residuals, an
uncorrected centre, and no diagnostic anywhere.

**One bad pick poisons the seed.** A spot clicked forty pixels off manufactures
short differences that are not lattice vectors either. So the seed is a small
search: the shortest few independent difference pairs are each fitted, and the
best *fit* is kept.

**Counting inliers always prefers the finer lattice.** Halving a cell explains
every spot it explained before, plus the mis-picked one, so inlier count rewards
exactly the models that assert most and predict least. Candidates are therefore
ranked by evidence. A cell of area $A$ scattered with nodes puts a fraction
$\pi t^2 / A$ of the plane within tolerance $t$ of some node, so one inlier
carries

$$ \log\frac{A}{\pi t^{2}} $$

nats of evidence, and a fit carries that times its inlier count. A lattice half
as coarse quadruples its node density and pays $\log 4$ per inlier for the
privilege — which one extra explained spot rarely covers.

The tolerance must be **the same for every candidate**, or the comparison is
empty: a tolerance proportional to each candidate's own basis makes the density
ratio scale-invariant and the criterion blind to the very thing it is meant to
detect. It is therefore a fixed fraction of the shortest observed separation
between two picks, which is a stand-in for picking precision.

### Two limits that cannot be engineered away

**A centre wrong by an exact lattice vector is undetectable.** Every spot is
still an exact node about the displaced origin; the residuals are zero and
nothing is geometrically wrong. What identifies the transmitted beam is that it
is the brightest thing on the plate — a judgement about intensity, not position.

**Beyond half a spacing, refinement becomes relabelling.** A fit free to move the
origin further will settle on a self-consistent assignment around a *different*
node, with residuals small enough to look convincing. The refinement is therefore
leashed at half the shortest spacing and reports when it reaches the leash.

### Reporting the basis

A lattice has infinitely many bases: $(\mathbf{a}, \mathbf{b})$ and
$(\mathbf{a}, \mathbf{b} - \mathbf{a})$ generate the same nodes. Any quantity
read off an arbitrary fitted basis therefore describes the arithmetic rather than
the crystal — a square lattice can come back as two vectors 135° apart, which is
correct and useless. Gauss reduction removes the freedom: repeatedly subtract
from the longer vector the integer multiple of the shorter that shortens it most,
swapping when the order changes. The reduced basis is unique up to sign and
exchange, its included angle lies in $[60°, 120°]$, and 90° then means the
lattice really is rectangular.

## 2. Scoring a candidate indexing

The solver ranks candidates by matched fraction and then residual. That is a sort
key, and it is deliberately not offered as a quality: it orders solutions and
says nothing about whether the best one is any good. Comparing candidates needs
the evidence and the policy kept apart.

### The evidence

For every indexed spot, the measured spacing against the calculated one:

$$ \frac{\Delta d}{d} = \frac{d_{\text{measured}} - d_{\text{calculated}}}{d_{\text{calculated}}} $$

For every pair of indexed spots, the measured angle between the two *measured*
$\mathbf{g}$ vectors against the angle between the corresponding *calculated*
ones. Neither comparison passes through the camera constant.

### The policy

Three terms, each mapped to $[0, 1]$ by

$$ A(x) = \frac{1}{1 + (x/t)^{s}} $$

which scores $\tfrac12$ at the tolerance $t$. A soft curve rather than a
threshold, because a threshold is discontinuous exactly where candidates are
hardest to tell apart; polynomial rather than Gaussian, so two badly wrong
solutions stay comparable instead of both underflowing to zero. The fused score
is the weighted mean of the three terms, normalized by the total weight, so 1
means perfect agreement on everything picked and 0.5 means disagreement at the
stated tolerances.

### Why angles outweigh lengths

A camera constant $L\lambda$ enters every measured length as a common factor and
every measured angle not at all. Set it 5% high and each measured $d$ rises by
exactly 5% while every interspot angle is unchanged. So:

- a **length** disagreement may be evidence about the *instrument*;
- an **angle** disagreement is evidence about the *crystallography*.

Hence the default weights: angle 1.5 against length 1.0. Coverage carries 2.0 —
the largest — because an unindexed spot is unexplained evidence that precision
elsewhere does not answer.

Two consequences worth stating. A solution with one indexed spot has no pair to
measure an angle between; that is *missing evidence*, not disagreement, so the
angle term is held neutral at 0.5 rather than scored zero. And the same relative
deviation appearing on **every** spot is the signature of a calibration error,
while a scatter of deviations is the signature of a wrong indexing — which is why
the per-spot table reports $\Delta d$ and not only an r.m.s.

### What the score cannot do

It cannot resolve what one pattern cannot. A zone axis and its reverse give
identical spot positions, so both score identically and always will. It does not
use intensities: relative intensities in a real pattern are dynamical and vary
with thickness and tilt, so a solution scored on them would be scored on the
specimen rather than on the crystallography.

## Computed checks

Both surfaces are pinned by executable worked examples that recompute on every
test run — see
[simulated SAED plates and the zone-axis atlas](../examples/generated/saed_practice_patterns.md):

- a beam centre displaced by 30 px in each direction is recovered **exactly**
  from eight exact lattice nodes;
- measured vectors stretched by 5% report a relative length deviation of exactly
  $1/1.05 - 1$ on every spot, with the angle term unmoved.

## Sources

- Williams, D. B. and Carter, C. B., *Transmission Electron Microscopy*, 2nd ed.,
  Springer, DOI: 10.1007/978-0-387-76501-3 — chapter 18, the camera equation and
  the role of the beam position.
- Edington, J. W., *Practical Electron Microscopy in Materials Science*,
  Macmillan (1975) — measuring and interpreting zone-axis patterns.
- Nguyen, P. Q. and Vallée, B. (eds), *The LLL Algorithm*, Springer (2010) —
  Gauss reduction as the two-dimensional case of lattice basis reduction.
