# Orientation Representations And The Equal-Volume Maps

*Relationship to the existing note.*
{doc}`/theory/orientation_representations_and_plane_direction_construction`
already records the classical conversions — matrix, quaternion, axis–angle, Rodrigues and
Rodrigues–Frank — and the construction of an orientation from a plane and a direction.
This note does not repeat them. It covers what that note does not: the two *equal-volume*
charts (homochoric and cubochoric), the measure argument that motivates them, and the
inversion that recovers $(hkl)[uvw]$ indices *from* an orientation.

A rotation is one geometric object; the numbers that denote it are many. This note fixes
the definitions PyTex uses for each representation, derives the two equal-volume charts
(homochoric and cubochoric) from first principles rather than quoting them, and states the
degeneracies a caller must know about. The implementation is
`pytex.core.representations`.

## Definitions

Let $\hat{n}$ be a unit vector and $\omega \in [0,\pi]$ the rotation angle about it, in
the active convention $v' = R v$. PyTex stores rotations as unit quaternions and derives
everything else:

\begin{align}
q &= \left(\cos\tfrac{\omega}{2},\; \hat{n}\sin\tfrac{\omega}{2}\right), &
  \boldsymbol{\rho} &= \hat{n}\tan\tfrac{\omega}{2}, \\
  \boldsymbol{\rho}_{F} &= \left(\hat{n},\; \tan\tfrac{\omega}{2}\right), &
  \mathbf{h} &= \hat{n}\left[\tfrac{3}{4}\left(\omega - \sin\omega\right)\right]^{1/3}.
\end{align}

Here $\boldsymbol{\rho}$ is the Rodrigues vector and $\boldsymbol{\rho}_{F}$ its
homogeneous (Frank) form, in which the magnitude is a projective coordinate: at
$\omega = \pi$ it is the point at infinity, a well-defined and exactly invertible
representation, whereas the 3-vector $\boldsymbol{\rho}$ overflows there and loses its
axis in the product. Finally $\mathbf{h}$ is the homochoric vector. Euler angles follow the axis sequences of
{doc}`/theory/euler_convention_handling`: $ZXZ$ for Bunge
$(\varphi_1, \Phi, \varphi_2)$, $ZYZ$ for Matthies/ABG $(\alpha, \beta, \gamma)$.

Two degeneracies are intrinsic and are reported rather than hidden. First, $q$ and $-q$
denote the same rotation; PyTex returns the representative with non-negative scalar part.
Second, at $\Phi = 0$ or $\Phi = \pi$ only $\varphi_1 \pm \varphi_2$ is determined, and
PyTex resolves it by setting the third angle to zero.

## Why an equal-volume chart is needed

The bi-invariant (Haar) measure on $SO(3)$, in axis-angle coordinates, is

$$
\mathrm{d}\mu = \frac{1}{\pi^{2}}\,(1 - \cos\omega)\,\mathrm{d}\omega\,\mathrm{d}\Omega_{\hat{n}}
$$ (eq-haar)

normalised so that $\int \mathrm{d}\mu = 1$. The factor $(1-\cos\omega)$ is the reason a
cloud of points drawn uniformly in $(\hat{n}, \omega)$ is *not* a uniform sample of
orientations, and the reason uniformly sampled Euler angles are not either. We therefore look
for a radial function $f(\omega)$ such that the Euclidean volume element of
$\mathbf{h} = f(\omega)\hat{n}$ reproduces {eq}`eq-haar`, that is

$$
f^{2}\,\frac{\mathrm{d}f}{\mathrm{d}\omega} \;\propto\; 1 - \cos\omega
  \quad\Longrightarrow\quad
  f(\omega) = \left[\tfrac{3}{4}\left(\omega - \sin\omega\right)\right]^{1/3}
$$

with the constant chosen so that $f$ has unit slope at the origin, $f(\omega) \approx \omega/2$. The whole of $SO(3)$ then fills the ball of radius

$$
R_{1} = f(\pi) = \left(\tfrac{3\pi}{4}\right)^{1/3} \approx 1.330670
$$

whose volume is $\tfrac{4}{3}\pi R_{1}^{3} = \pi^{2}$, the volume of $SO(3)$ under
 {eq}`eq-haar` before normalisation. Antipodal points of the bounding sphere denote the same
rotation — a rotation by $+\pi$ and by $-\pi$ about the same axis are equal — which is
what makes the ball a model of $SO(3)$ rather than of the unit quaternions.

As a check that is used as a test rather than as prose, the mean rotation angle of a uniform
sample follows from {eq}`eq-haar`:

$$
\langle \omega \rangle = \int_{0}^{\pi} \omega\,\frac{1-\cos\omega}{\pi}\,\mathrm{d}\omega
  = \frac{\pi}{2} + \frac{2}{\pi} \approx 126.4756^{\circ}.
$$ (eq-meanangle)

## The cubochoric cube and the map onto the ball

A ball is equal-volume but awkward: a regular grid inside it is not a regular grid of
anything, and it has no product structure. Cubochoric coordinates fix that by mapping the
ball onto a *cube* of the same volume $\pi^{2}$, hence of edge

$$
a_{p} = \pi^{2/3} \approx 2.145029
$$

so that a uniform Cartesian grid in the cube is a uniform grid in orientation space. The map
is that of Roşca, Morawiec and De Graef. PyTex derives it from two conditions, both of
which are asserted directly in the test suite:

**Condition 1: nested surfaces.** The sub-cube of half-edge $z$ maps onto the sphere
enclosing the same volume,

$$
(2z)^{3} = \tfrac{4}{3}\pi r^{3}
  \quad\Longrightarrow\quad
  r(z) = z\left(\tfrac{6}{\pi}\right)^{1/3}
$$ (eq-nested)

which at $z = a_{p}/2$ returns $r = (6\pi)^{1/3}/2 = R_{1}$, as it must.

**Condition 2: each face maps to its own spherical sector.** The six images must tile
the sphere and the map commutes with the octahedral symmetry, so the boundary between the
$+z$ and $+x$ images lies on the mirror plane $x = z$. The $+z$ face therefore maps to
the curvilinear square $\{n_{z} \ge |n_{x}|, |n_{y}|\}$, of solid angle $2\pi/3$, and
never to a spherical cap — six caps of that solid angle would overlap.

Rescale by $(\pi/6)^{1/6}$ so that {eq}`eq-nested` reads $r = \sqrt{6/\pi}\,z$. Writing
the volume element as $\mathrm{d}V = r^{2}\,\mathrm{d}r\,\mathrm{d}\Omega$ and the cube
element in the face coordinates $X = x/z$, $Y = y/z$ as $z^{2}\,\mathrm{d}X\,\mathrm{d}Y \,\mathrm{d}z$, equality of volumes gives $r^{2}r'\,J_{n} = z^{2}$, so the direction map
$\hat{n}(X,Y)$ must have constant Jacobian $J_{n} = \pi/6$: it is an area-preserving map
of the square $[-1,1]^{2}$ onto the spherical sector.

That map factors as a planar wedge followed by a Lambert azimuthal equal-area lift. For
$|y| \le |x|$ within the $+z$ pyramid, put

```{math}
:label: eq-wedge
\begin{aligned}
  \alpha &= \frac{\pi}{12}\,\frac{y}{x}, &
  k &= \frac{2^{1/4}\sqrt{6/\pi}\;x}{\sqrt{\sqrt{2}-\cos\alpha}}, \\
  (T_{1}, T_{2}) &= k\left(\sqrt{2}\cos\alpha - 1,\; \sqrt{2}\sin\alpha\right),
\end{aligned}
```

with $x$ and $y$ exchanging roles on the far side of the face diagonal. The map
 {eq}`eq-wedge` carries the octant onto the plane wedge of half-angle $\pi/4$: at
$y = 0$ the azimuth is zero, and at $y = x$ it is
$\arctan\!\big(\sqrt{2}\sin\tfrac{\pi}{12} / (\sqrt{2}\cos\tfrac{\pi}{12}-1)\big) = \pi/4$.
The prefactor $2^{1/4}\sqrt{6/\pi}$ is not fitted: requiring the face edge at azimuth
$\varphi$ to land on the sector boundary $\theta_{\max}(\varphi) = \arctan(1/\cos\varphi)$,
whose Lambert radius is $2R\sin(\theta_{\max}/2)$, gives at $\varphi = 0$

$$
\text{prefactor}
  = \frac{2\sqrt{6/\pi}\,\sin(\pi/8)}{\sqrt{\sqrt{2}-1}}
  = \sqrt{\tfrac{6}{\pi}}\sqrt{\frac{2-\sqrt{2}}{\sqrt{2}-1}}
  = 2^{1/4}\sqrt{\tfrac{6}{\pi}}
$$

using $(2-\sqrt{2})/(\sqrt{2}-1) = \sqrt{2}$; and the same value then satisfies the boundary
condition at every other azimuth, which is the check that the ansatz is the right one.

The lift onto the sphere of radius $R = \sqrt{6/\pi}\,z$ is Lambert azimuthal equal-area,

$$
\mathbf{b} = \left(T_{1}\sqrt{1 - \tfrac{\rho^{2}}{4R^{2}}},\;
                     T_{2}\sqrt{1 - \tfrac{\rho^{2}}{4R^{2}}},\;
                     R - \tfrac{\rho^{2}}{2R}\right),
  \qquad \rho^{2} = T_{1}^{2}+T_{2}^{2}
$$ (eq-lambert)

which satisfies $|\mathbf{b}| = R$ identically and preserves area. The remaining five
pyramids follow from the octahedral equivariance: a cyclic permutation of the axes is a
rotation of both cube and ball, and the mirror $z \to -z$ negates the third component of the
image and nothing else.

## The inverse

Both steps invert in closed form, so no root finding is needed. From {eq}`eq-lambert`,
$\rho^{2} = 2R(R - b_{3})$ with $R = |\mathbf{b}|$, and $(T_{1},T_{2})$ is
$(b_{1},b_{2})$ rescaled to that radius. The angular part of {eq}`eq-wedge` inverts
through $\sqrt{2}\sin(\alpha - \varphi) = -\sin\varphi$, that is

$$
\alpha = \varphi - \arcsin\!\left(\frac{\sin\varphi}{\sqrt{2}}\right)
$$

after folding the wedge onto its positive-principal half, and the radial part gives
$x = \rho\big/\big(\text{prefactor}\,\sqrt{(3-2\sqrt{2}\cos\alpha)/(\sqrt{2}-\cos\alpha)}\big)$
and $y = 12\alpha x/\pi$. The homochoric-to-axis-angle inverse is the one step with no
closed form: $\omega - \sin\omega = \tfrac{4}{3}\|\mathbf{h}\|^{3}$ is solved by a
vectorized bisection on $[0,\pi]$, where the left-hand side is strictly increasing.

## Naming an orientation: $(hkl)[uvw]$

A texture component is named by the crystal plane lying in the sheet plane and the crystal
direction along the rolling direction. Given the orientation matrix $g$
(crystal-to-specimen) and the specimen axes $\mathbf{n}_{S}$ (ND) and $\mathbf{d}_{S}$
(RD), the crystal-frame Cartesian images are $g^{\mathsf{T}}\mathbf{n}_{S}$ and
$g^{\mathsf{T}}\mathbf{d}_{S}$; their components in the reciprocal and direct bases
respectively are the (generally irrational) indices. PyTex reports the nearest integer triple
within a search bound *together with the residual angle*, because only a measure-zero set
of orientations has an exact ideal label, and a component name quoted without its deviation is
a claim the data does not support.

## Limits and cross-references

- The equal-volume charts are charts, not algebras: rotations are composed as quaternions, never by adding homochoric or cubochoric coordinates.
- Neither chart is symmetry-reduced. Fundamental-zone reduction is the separate concern of {doc}`/theory/fundamental_region_reduction`.
- Rodrigues and homochoric vectors are both bare triples; PyTex rejects an out-of-ball homochoric input rather than clipping it, because that is the only way the mix-up can be caught.
- Euler conventions and their aliases: {doc}`/theory/euler_convention_handling`.
- Orientation-space distance and disorientation: {doc}`/theory/orientation_space_and_disorientation`.

## Normative And Informative Sources

- H.-J. Bunge, *Texture Analysis in Materials Science*, Butterworths, 1982 — the $ZXZ$ Euler convention and the orientation-distribution formalism.
- F. C. Frank, “Orientation mapping”, *Metallurgical Transactions A* **19** (1988) 403–408 — Rodrigues–Frank space.
- A. Morawiec, *Orientations and Rotations*, Springer, 2004 — the invariant measure {eq}`eq-haar` and the representation catalogue.
- D. Roşca, A. Morawiec and M. De Graef, “A new method of constructing a grid in the space of 3D rotations and its applications to texture analysis”, *Modelling and Simulation in Materials Science and Engineering* **22** (2014) 075013, `doi:10.1088/0965-0393/22/7/075013` — the equal-volume cube-to-ball map.
- M. De Graef, *Introduction to Conventional Transmission Electron Microscopy*, Cambridge University Press, 2003 — rotation conventions used throughout the diffraction chapters.
