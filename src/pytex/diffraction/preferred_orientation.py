"""Preferred-orientation corrections for powder diffraction intensities.

A powder pattern's relative intensities assume the crystallites are randomly
oriented. Real specimens rarely are: rolled sheet, pressed powder, and anything
with a plate-like or needle-like habit develop texture, and texture redistributes
intensity between reflections. Uncorrected, that redistribution is easily
misread as a wrong structure, a wrong phase fraction, or a wrong site occupancy.

Two corrections are provided, and they answer different questions.

**March-Dollase** models a single fibre texture with one adjustable parameter.
It is the standard Rietveld correction, it needs no texture measurement, and it
is what to use when preferred orientation is a nuisance to be fitted away.

**ODF weighting** drives the intensities from a measured or modelled orientation
distribution instead. This is the physically direct route: the intensity of a
reflection is proportional to the pole density along the scattering vector, and
an ODF gives exactly that. It requires a texture measurement, and in exchange it
makes no assumption of fibre symmetry, handles arbitrarily complex textures, and
introduces no fitted parameter. Coupling the PyTex texture core into powder
intensities this way is the reason both live in one library.

Both models express the same quantity: a per-reflection multiplicative factor in
multiples of the random-powder intensity. A factor of 1 means "as a random
powder would give"; above 1 means enhanced by texture, below 1 suppressed.

Normalization
-------------
Both models are constructed so that a texture-free specimen gives a factor of
exactly 1 for every reflection: ``march_coefficient = 1`` for March-Dollase, and
a uniform ODF for the pole-density model. Preferred orientation *redistributes*
diffracted intensity between reflections; it does not create or destroy it.

See Also
--------
pytex.texture.ODF : The orientation distribution the pole-density model reads.
pytex.diffraction.xrd.generate_xrd_pattern : Applies a model to a whole pattern.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np
from numpy.typing import ArrayLike

from pytex.core._arrays import FloatArray, freeze_array, normalize_vector
from pytex.core.miller import MillerPlane
from pytex.core.notation import format_plane_family_indices
from pytex.core.provenance import ProvenanceRecord
from pytex.texture.models import KernelSpec, random_pole_density

if TYPE_CHECKING:  # pragma: no cover - import cycle guard for typing only
    from pytex.texture.harmonics import HarmonicODF
    from pytex.texture.models import ODF

#: Angles closer than this to a pole are treated as coincident when the March
#: function is evaluated, which matters only for reporting, not for the value.
_ANGLE_ATOL = 1e-12

#: Quadrature nodes used to evaluate the spherical mean of a smoothing kernel.
#: The integrand is smooth in ``cos(omega)``, so Gauss-Legendre converges fast;
#: 512 nodes is far past convergence for any physically sensible halfwidth.
def _random_pole_density(kernel: KernelSpec) -> float:
    """The pole density a *random* texture produces under this kernel.

    Thin alias for :func:`pytex.texture.models.random_pole_density`, which is
    the single implementation; see there for why the correction is needed.
    """

    return random_pole_density(kernel)


@runtime_checkable
class PreferredOrientationModel(Protocol):
    """The contract a preferred-orientation model satisfies.

    A model turns crystal planes into multiplicative intensity factors in
    multiples of the random-powder intensity, and can explain itself in prose.
    Implementing this protocol is all that is required to plug a new texture
    model into the powder-pattern path.
    """

    def factors(self, planes: Sequence[MillerPlane]) -> FloatArray:
        """Per-plane intensity factors, in multiples of random."""
        ...

    def describe(self) -> str:
        """Convention-explicit prose summary of the model."""
        ...


def march_dollase_factors(
    angles_rad: ArrayLike,
    march_coefficient: float,
) -> FloatArray:
    """The March distribution evaluated at given angles.

    Purpose
    -------
    The single-parameter fibre-texture model of March (1932) as applied to
    diffraction by Dollase (1986):

    .. math::

        P(\\alpha) = \\left( r^{2}\\cos^{2}\\alpha
                     + \\frac{\\sin^{2}\\alpha}{r} \\right)^{-3/2}

    where :math:`\\alpha` is the angle between a plane normal and the
    preferred-orientation axis and :math:`r` is the March coefficient.

    Interpreting ``r``
    ------------------
    ``r = 1`` is a random powder and gives exactly 1 at every angle. ``r < 1``
    describes a **plate-like** habit whose plate normals cluster along the
    preferred-orientation axis, so reflections from those planes are enhanced
    (:math:`P(0) = r^{-3} > 1`). ``r > 1`` describes a **needle-like** habit and
    suppresses them. At the extremes :math:`P(0) = r^{-3}` and
    :math:`P(\\pi/2) = r^{3/2}`.

    Parameters
    ----------
    angles_rad : ArrayLike
        Angles between plane normals and the preferred-orientation axis.
    march_coefficient : float
        The March coefficient ``r``; strictly positive.

    Returns
    -------
    FloatArray
        Factors, read-only. The result is always at least one-dimensional, so a
        scalar angle yields a length-1 array rather than a bare float.

    Notes
    -----
    The distribution is exactly normalized over the sphere: its average over a
    uniform distribution of directions is 1 for every ``r``, since

    .. math::

        \\int_{0}^{1}\\bigl( (r^{2} - r^{-1})u^{2} + r^{-1}\\bigr)^{-3/2}\\,du
          = \\frac{r}{r} = 1 .

    Preferred orientation therefore redistributes intensity rather than creating
    it, and this identity is pinned as a test.
    """

    if not np.isfinite(march_coefficient) or march_coefficient <= 0.0:
        raise ValueError("march_coefficient must be finite and strictly positive.")
    angles = np.asarray(angles_rad, dtype=np.float64)
    if np.any(~np.isfinite(angles)):
        raise ValueError("angles_rad must be finite.")
    if march_coefficient == 1.0:
        # The general expression evaluates cos^2 + sin^2, which is 1 only to
        # within an ulp. The random-powder case is the default assumption of
        # every uncorrected pattern, so it returns exactly 1 and "no preferred
        # orientation" leaves intensities bit-identical rather than nearly so.
        return freeze_array(np.ones_like(angles))
    cosines = np.cos(angles)
    sines = np.sin(angles)
    radicand = march_coefficient**2 * cosines**2 + sines**2 / march_coefficient
    return freeze_array(np.ascontiguousarray(radicand ** (-1.5)))


@dataclass(frozen=True, slots=True)
class MarchDollaseModel:
    """Single-parameter fibre-texture correction for powder intensities.

    Purpose
    -------
    The standard Rietveld preferred-orientation correction. It assumes one
    fibre texture, described by the crystal plane whose normals align with the
    specimen axis and by a single strength parameter, and needs no texture
    measurement.

    Method
    ------
    For each reflection, the March function is averaged over the reflection's
    full symmetry family, since every symmetry-equivalent plane contributes at
    the same Bragg angle and each sits at its own angle to the
    preferred-orientation axis. This family averaging is what makes the
    correction symmetry-consistent rather than dependent on which family
    representative happened to be enumerated.

    Attributes
    ----------
    preferred_orientation : MillerPlane
        The plane whose normals cluster along the specimen axis — ``(001)`` for
        a basal-textured sheet, ``(100)`` for many pressed powders.
    march_coefficient : float
        Strictly positive. ``1`` is a random powder; below 1 is plate-like and
        enhances the preferred reflections; above 1 is needle-like and
        suppresses them. Values outside roughly ``0.3`` to ``3`` describe
        textures far stronger than a nuisance correction should be asked to
        absorb, and usually indicate that a measured ODF is the honest choice.
    provenance : ProvenanceRecord, optional

    Notes
    -----
    Fibre symmetry about the specimen axis is *assumed*, not checked. A sheet
    texture with distinct rolling and transverse behaviour violates it, and
    March-Dollase will absorb the discrepancy into a fitted ``r`` that has no
    physical meaning. Use :class:`ODFPreferredOrientationModel` there.

    See Also
    --------
    march_dollase_factors : The underlying distribution.
    """

    preferred_orientation: MillerPlane
    march_coefficient: float
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.march_coefficient) or self.march_coefficient <= 0.0:
            raise ValueError(
                "MarchDollaseModel.march_coefficient must be finite and strictly positive."
            )

    @property
    def is_random(self) -> bool:
        """Whether the model describes an untextured powder (``r = 1``).

        When true every factor is exactly 1, so the correction can be skipped
        without changing any result.
        """

        return bool(np.isclose(self.march_coefficient, 1.0, rtol=0.0, atol=_ANGLE_ATOL))

    def factors(self, planes: Sequence[MillerPlane]) -> FloatArray:
        """Family-averaged March-Dollase factors for the given reflections.

        Parameters
        ----------
        planes : sequence of MillerPlane
            Reflections to correct. Each must be on the same phase as
            :attr:`preferred_orientation`, so the angle between their normals is
            a meaningful crystallographic quantity; a mismatch raises.

        Returns
        -------
        FloatArray
            ``(n,)`` factors in multiples of random, read-only.
        """

        axis = self.preferred_orientation.normal_cartesian
        values = np.empty(len(planes), dtype=np.float64)
        for index, plane in enumerate(planes):
            if plane.phase != self.preferred_orientation.phase:
                raise ValueError(
                    "MarchDollaseModel.factors requires every plane to share the "
                    "preferred-orientation phase."
                )
            family = plane.symmetry_equivalents()
            normals = family.normals_cartesian()
            # A plane and its opposite normal are the same plane, so the angle is
            # taken to the nearer of the two poles.
            cosines = np.abs(normals @ axis)
            angles = np.arccos(np.clip(cosines, -1.0, 1.0))
            values[index] = float(
                np.mean(march_dollase_factors(angles, self.march_coefficient))
            )
        return freeze_array(np.ascontiguousarray(values))

    def describe(self) -> str:
        """Convention-explicit prose summary of the correction.

        Names the preferred-orientation plane family, the coefficient and the
        habit it implies, the family-averaging convention, and the fibre-symmetry
        assumption that the model does not verify.
        """

        indices = tuple(int(value) for value in self.preferred_orientation.indices)
        family = format_plane_family_indices(indices, style="plain")
        if self.is_random:
            habit = "a random powder, so every factor is exactly 1"
        elif self.march_coefficient < 1.0:
            habit = (
                f"a plate-like habit: {family} normals cluster along the specimen axis, so "
                f"those reflections are enhanced by up to {self.march_coefficient**-3.0:.3f}x"
            )
        else:
            habit = (
                f"a needle-like habit: {family} normals avoid the specimen axis, so those "
                f"reflections are suppressed to as little as {self.march_coefficient**-3.0:.3f}x"
            )
        return (
            f"March-Dollase preferred-orientation correction about {family} with March "
            f"coefficient r = {self.march_coefficient:.4f}, describing {habit}. Factors are "
            "in multiples of the random-powder intensity and are averaged over each "
            "reflection's full symmetry family. The model assumes a single fibre texture "
            "with rotational symmetry about the specimen axis; that assumption is not "
            "verified here, and a sheet texture violating it will be absorbed into a "
            "physically meaningless r. The March distribution integrates to 1 over the "
            "sphere, so the correction redistributes intensity rather than creating it."
        )


@dataclass(frozen=True, slots=True)
class ODFPreferredOrientationModel:
    """Powder intensities driven by a measured or modelled texture.

    Purpose
    -------
    The physically direct preferred-orientation correction, and the reason the
    PyTex texture and diffraction layers belong in one library. The intensity of
    a powder reflection is proportional to the density of ``{hkl}`` poles lying
    along the scattering vector; an orientation distribution function gives that
    pole density directly, in multiples of random, with no fitted parameter and
    no assumption of fibre symmetry.

    Method
    ------
    For each reflection the ODF is evaluated for pole density along the
    scattering direction, over the whole ``{hkl}`` symmetry family — which is
    what a powder reflection actually contains.

    That evaluation returns a kernel-weighted response rather than a value in
    multiples of random: the smoothing kernel peaks at 1 rather than integrating
    to 1, so a uniform texture yields the kernel's spherical mean, not unity.
    The response is therefore divided by that mean, computed in closed form by
    quadrature. This is what makes a uniform texture give a factor of exactly 1
    and the correction interpretable as multiples of random.

    Geometry
    --------
    In symmetric Bragg-Brentano reflection geometry the scattering vector lies
    along the specimen normal at every angle, which is why the default
    ``scattering_direction`` is ND. That is an approximation for asymmetric or
    transmission geometries, where the scattering direction moves with
    ``2*theta``; supply the appropriate direction explicitly in that case.

    Attributes
    ----------
    odf : ODF or HarmonicODF
        The orientation distribution. Its specimen frame defines what the
        scattering direction is expressed in.
    scattering_direction : np.ndarray
        Specimen-frame direction of the scattering vector; normalized on
        construction. Defaults to the specimen normal ``(0, 0, 1)``.
    provenance : ProvenanceRecord, optional

    Notes
    -----
    The correction is only as good as the ODF. An ODF estimated from too few
    pole figures, or over-smoothed by a wide kernel, will flatten real texture
    towards 1 and under-correct; check the reconstruction residuals before
    trusting the factors.
    """

    odf: ODF | HarmonicODF
    scattering_direction: np.ndarray = field(
        default_factory=lambda: np.array([0.0, 0.0, 1.0], dtype=np.float64)
    )
    provenance: ProvenanceRecord | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "scattering_direction", normalize_vector(self.scattering_direction)
        )

    @property
    def specimen_frame_name(self) -> str:
        """Name of the specimen frame the scattering direction is expressed in.

        Read from whichever ODF representation was supplied: the discrete ODF
        carries it on its support orientations, the harmonic ODF exposes it
        directly.
        """

        frame = getattr(self.odf, "specimen_frame", None)
        if frame is None:
            frame = self.odf.orientations.specimen_frame  # type: ignore[union-attr]
        return str(frame.name)

    def factors(self, planes: Sequence[MillerPlane]) -> FloatArray:
        """Pole-density factors for the given reflections.

        Parameters
        ----------
        planes : sequence of MillerPlane
            Reflections to correct, on the phase the ODF describes.

        Returns
        -------
        FloatArray
            ``(n,)`` pole densities along the scattering direction, in multiples
            of random, read-only. Negative values — possible only from a
            truncated harmonic ODF, where they are expansion artefacts rather
            than densities — are clipped to zero, since a negative diffracted
            intensity has no meaning.
        """

        direction = np.asarray(self.scattering_direction, dtype=np.float64)[None, :]
        reference = _random_pole_density(self._pole_kernel())
        values = np.empty(len(planes), dtype=np.float64)
        for index, plane in enumerate(planes):
            density = self.odf.evaluate_pole_density(
                plane.to_crystal_plane(),
                direction,
                include_symmetry_family=True,
            )
            values[index] = float(np.asarray(density).reshape(-1)[0]) / reference
        return freeze_array(np.ascontiguousarray(np.clip(values, 0.0, None)))

    def _pole_kernel(self) -> KernelSpec:
        """The smoothing kernel the ODF uses for pole-density evaluation."""

        kernel = getattr(self.odf, "pole_kernel", None)
        if kernel is None:
            kernel = self.odf.kernel  # type: ignore[union-attr]
        return kernel

    def describe(self) -> str:
        """Convention-explicit prose summary of the correction.

        Names the scattering direction and its frame, states that factors are
        pole densities in multiples of random over the whole ``{hkl}`` family,
        and records the Bragg-Brentano geometry assumption behind the default
        direction.
        """

        direction = np.asarray(self.scattering_direction, dtype=np.float64)
        frame = self.specimen_frame_name
        return (
            "ODF-weighted preferred-orientation correction: each reflection is scaled by the "
            f"pole density of its full {{hkl}} family along the scattering direction "
            f"({direction[0]:.4f}, {direction[1]:.4f}, {direction[2]:.4f}) in the "
            f"'{frame}' specimen frame. Factors are in multiples of a random distribution, so "
            "a uniform texture gives exactly 1 and no correction. The default direction is the "
            "specimen normal, which is the scattering direction at every angle in symmetric "
            "Bragg-Brentano reflection geometry but not in asymmetric or transmission "
            "geometries. No fibre symmetry is assumed; the correction is only as reliable as "
            "the ODF it reads."
        )


def preferred_orientation_factor_table(
    planes: Sequence[MillerPlane],
    model: PreferredOrientationModel,
) -> tuple[tuple[str, float], ...]:
    """A labelled table of the factors a model assigns to given reflections.

    Purpose
    -------
    Make a correction auditable. Reporting only the corrected pattern hides
    which reflections were enhanced and which suppressed, and by how much; that
    breakdown is what distinguishes a plausible texture correction from one
    absorbing a structural error.

    Returns
    -------
    tuple of (str, float)
        ``({hkl} label, factor)`` pairs in the input order, with the family
        brackets that a powder reflection warrants.
    """

    factors = model.factors(planes)
    return tuple(
        (
            format_plane_family_indices(
                tuple(int(value) for value in plane.indices), style="plain"
            ),
            float(factor),
        )
        for plane, factor in zip(planes, factors, strict=True)
    )


__all__ = [
    "MarchDollaseModel",
    "ODFPreferredOrientationModel",
    "PreferredOrientationModel",
    "march_dollase_factors",
    "preferred_orientation_factor_table",
]
