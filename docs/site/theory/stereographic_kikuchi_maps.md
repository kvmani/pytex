# Stereographic Kikuchi Maps And Zone-Axis Routing

This note fixes the geometry, the combinatorics and the error analysis behind
`pytex.diffraction.kikuchi_map`:
`StereographicKikuchiMap`, `KikuchiMapBand`,
`KikuchiMapZoneAxis`, `KikuchiRoute`,
`compute_kikuchi_map` and `plan_kikuchi_route`, with the figure
produced by `pytex.plotting.plot_kikuchi_map`.

It is the companion to the note on Kikuchi bands and the gnomonic projection.
That note treats one pattern on one detector at one known orientation; this one
treats the whole band network of the crystal, on the sphere, as the object a TEM
operator navigates by. The user-facing pages are tutorial 30 and the workflow page
on Kikuchi geometry; the executable worked example is
`diffraction-kikuchi-map-zone-axis-tilt-angles`.

## The Band On The Sphere

The origin of a band is unchanged: inelastic scattering illuminates every
direction, and for a lattice plane of spacing $d$ the directions making the Bragg
angle with the plane diffract, where $\sin\theta_B = \lambda/2d$. Those directions
form two cones of semi-angle $90^\circ - \theta_B$ about the plane normal
$\mathbf{g}$ — the Kossel cones. On the unit sphere of directions $S^2$ this is a
band bounded by two small circles,

$$
\hat{\mathbf{r}} \cdot \hat{\mathbf{g}} = \pm \sin\theta_B
$$

whose midline is the great circle $\hat{\mathbf{r}} \cdot \hat{\mathbf{g}} = 0$:
the trace of the plane. The angular width is exactly $2\theta_B$, and since

$$
2\theta_B = 2\arcsin\!\left(\frac{\lambda}{2d}\right) \approx \frac{\lambda}{d}
$$ (eq-band-width)

width is a decreasing function of spacing. At $200$ kV, where
$\lambda = 0.02508$ Å, a $2$ Å spacing gives $2\theta_B = 0.72^\circ$. Two
consequences are worth separating explicitly, because they are routinely
conflated:

- the *widest* bands of a map come from the *smallest* spacings, hence from high-index planes;
- the *strongest* bands come from the largest spacings, hence from low-index planes, because $|F_{\mathbf{g}}|$ falls with $|\mathbf{g}|$.

For nickel at $200$ kV and $|h|,|k|,|l| \le 4$, PyTex reports $\{111\}$ as the
strongest band at $0.706^\circ$ wide, and $\{442\}$ as the widest at
$2.447^\circ$ carrying $5$ percent of the $\{111\}$ intensity.

## Why The Projection Must Be Stereographic

An atlas has to display a hemisphere of directions. Write $\psi$ for the polar
angle from the map centre. The three candidate azimuthal projections give radii

$$
r_{\text{gnomonic}} = \tan\psi, \qquad
  r_{\text{stereographic}} = \tan\frac{\psi}{2}, \qquad
  r_{\text{equal-area}} = 2\sin\frac{\psi}{2}
$$

The gnomonic radius diverges as $\psi \to 90^\circ$: the gnomonic image of a
hemisphere is the entire plane, so a gnomonic atlas does not exist. At
$\psi = 89.99^\circ$ the gnomonic radius is $5729.6$ and the stereographic radius
is $0.99983$. That is why `pytex.diffraction.kikuchi` uses gnomonic
coordinates — where band centre lines are exactly straight, which is what reading
a physical detector wants — and this module cannot.

Between the two bounded projections, stereographic is chosen because it is
conformal. Every number an operator takes off a Kikuchi map is an angle the
goniometer must turn through, and conformality means an angle subtended near the
rim and the same angle near the centre are drawn the same. The equal-area
projection preserves solid angle instead, which is the correct trade for a texture
pole figure whose quantity of interest is a density per unit solid angle, and the
wrong one here. `plot_kikuchi_map` accepts either and defaults to
stereographic.

## Zone Axes As Exact Combinatorics

A direction $[uvw]$ lies on the centre line of $(hkl)$ if and only if

$$
hu + kv + lw = 0
$$ (eq-zone-law)

the Weiss zone law, which is the perpendicularity condition
$\mathbf{g}\cdot\mathbf{r} = 0$ written in the dual bases: the reciprocal-basis
components of the plane contract with the direct-basis components of the direction
without any metric tensor. Equation {eq}`eq-zone-law` is an integer condition,
so the incidence structure of a Kikuchi map is exact combinatorics rather than a
numerical near-coincidence, and PyTex tests it as such.

The zone axis at the crossing of two bands is
$\mathbf{u} \propto \mathbf{g}_1 \times \mathbf{g}_2$, and the *order* of an
axis is the number of reflections in the map satisfying
equation {eq}`eq-zone-law` for it. Order, not index, is the operationally
relevant quantity: it is the number of bands that meet at that point on the
screen, so a four-band intersection is unmistakable and a two-band one is a guess.
The map orders its axes by decreasing order and then by increasing polar angle.

**One band per plane trace.**

Every order of a reflection has the same centre line: $(222)$ traces the same
great circle as $(111)$. Retaining all of them would draw coincident lines and
would multiply the order of every axis they both pass through, corrupting the one
number an operator reads. `compute_kikuchi_map` therefore collapses each
family of collinear triples onto its lowest *allowed* order. The
representative need not be coprime: in a face-centred lattice the $\{100\}$ trace
is drawn by $(200)$, because $(100)$ is extinguished by the centring condition.
`include_higher_orders` restores the unfolded set for studies of
higher-order line positions, where the differing Bragg angles are the point.

**Antipodal identification.**

A zone axis is a *line*, not a direction: the beam traverses the crystal one
way or the other and the pattern is the same, so $[uvw]$ and
$[\bar{u}\bar{v}\bar{w}]$ denote one axis. The map stores one representative per
antipodal pair. Taking the upper hemisphere is not sufficient on its own, because
an axis on the equator has $z = 0$ in both senses; the tie is broken
lexicographically on $x$ then $y$, and a numerically-zero $z$ is snapped to $+0$
first, since a one-hemisphere projection folds on the sign of $z$ and $-0$ carries
a set sign bit.

## Routing: The Band Is The Geodesic

Two zone axes $\mathbf{u}_1, \mathbf{u}_2$ lie on a common band exactly when some
reflection is perpendicular to both, that is when the plane they span is a
rational lattice plane. The shortest arc on $S^2$ between two directions also lies
in the plane they span. Hence

> the great circle an experienced operator follows by eye, tracking a single
> Kikuchi band, is the geodesic between the two zone axes.

This is not a coincidence to be admired but a consequence of both objects being
defined by the same plane. For $[001]$ and $[111]$ in a cubic lattice the spanning
plane is $(1\bar{1}0)$ and the arc length is $\arccos(1/\sqrt{3}) =
54.7356^\circ$.

When no single band joins the endpoints, `plan_kikuchi_route` searches
the zone-axis network. Nodes are zone axes; an edge joins two axes that share at
least one band and are separated by no more than `max_leg_deg`; edge
weight is the angular separation. Dijkstra's algorithm then minimizes total stage
travel among band-followable routes. Splitting an arc at points lying on the same
great circle does not lengthen it, so multi-hop routing along one band is free in
travel.

## Why Several Short Hops Beat One Long One

The reason to prefer several hops is an exact statement about error propagation,
not a preference for tidiness.

An operator at a zone axis knows the beam direction but not the azimuthal
orientation of the crystal about it: the rotation between the diffraction pattern
and the stage axes is calibrated only to a degree or two. Write $\delta\varphi$
for that uncertainty. The intended tilt is by $\theta$ about an axis $\mathbf{n}$
perpendicular to the current zone axis $\mathbf{u}_0$, but the axis actually
driven is $\mathbf{n}' = R(\mathbf{u}_0, \delta\varphi)\,\mathbf{n}$. Conjugation
gives

$$
R(\mathbf{n}', \theta)
  = R(\mathbf{u}_0, \delta\varphi)\,R(\mathbf{n}, \theta)\,
    R(\mathbf{u}_0, \delta\varphi)^{-1}
$$

and since $R(\mathbf{u}_0, \cdot)$ fixes $\mathbf{u}_0$, applying this to
$\mathbf{u}_0$ yields $R(\mathbf{u}_0, \delta\varphi)\,\mathbf{u}_1$. The achieved
direction is therefore the intended one rotated about the *starting* axis by
the calibration error, and the angular miss is exactly

$$
\Delta = 2\arcsin\!\left(\sin\frac{\delta\varphi}{2}\,\sin\theta\right)
         \;\approx\; \delta\varphi\,\sin\theta .
$$ (eq-tilt-miss)

The miss scales with the sine of the *hop length*. Split the excursion into
$n$ hops and re-index the orientation at each one: every hop contributes
$\delta\varphi\sin(\theta/n)$ with an independent error, so the expected total
goes as $\sqrt{n}\,\sin(\theta/n)$. For the $[001] \to [111]$ excursion with
$\delta\varphi = 2^\circ$ this is $1.63^\circ$ for one hop, $1.08^\circ$ for
three, and $0.78^\circ$ for six, at identical total travel. The gain saturates
because the $\sqrt{n}$ accumulation eventually offsets the shorter legs, so an
optimum exists near four hops rather than a monotone improvement.
`DEFAULT_ROUTE_MAX_LEG_DEG` is $30^\circ$ for this reason.

Equation {eq}`eq-tilt-miss` is the quantitative form of the argument in
§11 of the TEM tilt-navigation foundation document, and it is verified against a
direct rotation simulation in tutorial 30.

## Intensities, And What The Map Does Not Model

Band intensities use the Mott–Bethe electron structure factor of
`pytex.diffraction.scattering`, not the atomic-number proxy in
`pytex.diffraction.kinematic`. The proxy replaces $f_e(s)$ by $Z$, which is
independent of $s$; on a monatomic phase without Debye–Waller factors it
therefore assigns every allowed reflection the same intensity, and the band
ordering — which is the whole visual grammar of the map — would carry no
information. With the Mott–Bethe factor, nickel orders as $\{111\}, \{200\},
\{220\}$ and zirconium as $\{10\bar{1}1\}, \{11\bar{2}0\}$, as the reflection
tables have them.

The map remains kinematic and geometric. It does not predict the
excess–deficiency asymmetry across a band, dynamical contrast, or higher-order
Laue zone effects; the note on dynamical CBED covers what a Bloch-wave calculation
adds. It carries no foil thickness, no absorption and no detector: it is a map of
the crystal, and it says which way to turn rather than what the screen will look
like. Route angles are angles between crystal directions, so converting them into
stage $\alpha$ and $\beta$ and checking them against the holder envelope is the
job of `pytex.tem.navigation` and `pytex.tem.path`, which
`KikuchiRoute.describe` states explicitly. Finally, the search minimizes
travel rather than risk: equation {eq}`eq-tilt-miss` would sometimes prefer a
longer path through more intersections, and that trade is not modelled.

## References

- N. K. Levine, W. L. Bell and G. Thomas, *Further applications of Kikuchi diffraction patterns; electron diffraction pattern maps*, J. Appl. Phys.\ **37** (1966) 2141. The original montaged maps.
- J. W. Edington, *Practical Electron Microscopy in Materials Science*, Monograph 2, Philips Technical Library (1975). Kikuchi maps as an operating technique.
- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy*, 2nd ed., Springer (2009), Ch. 19. Kikuchi line geometry and its use for tilting.
- A. J. C. Wilson and E. Prince (eds.), *International Tables for Crystallography*, Vol. C, 2nd ed., Kluwer (1999). Reciprocal-lattice conventions and the zone law.
- J. P. Snyder, *Map Projections — A Working Manual*, USGS Professional Paper 1395 (1987), §§20–21. The azimuthal projections and what each preserves.
