# ruff: noqa: RUF001, RUF002
"""The one declaration of "which X-rays", shared by every diffractometer operation.

Purpose
-------
Every operation that turns an angle into a spacing needs the wavelength, and
before this module each declared its own short list of tube lines. That list
could not describe a synchrotron beamline, a Kα1-monochromated source or a
liquid-metal jet, and the lists had begun to differ from one view to the next.
Here the choice is declared once: the characteristic doublets of the common
anodes, and **a monochromatic beam of any wavelength**, with the polarization
that matters for a synchrotron.

A request names its radiation with three parameters, always under the same
names so a scan analysed in one view means the same thing in the next:

``radiation``
    A tube line, or ``"monochromatic"``.
``wavelength_angstrom``
    Used for ``"monochromatic"`` only. ``λ[Å] = 12.3984 / E[keV]``.
``polarization_fraction``
    Fraction of the beam polarized perpendicular to the scattering plane, for
    ``"monochromatic"`` only: 0.5 unpolarized, about 0.95 for a synchrotron with
    a vertical scattering plane, about 0.05 for a horizontal one.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

from pytex.app.errors import InvalidInputError
from pytex.app.registry import ChoiceParameter, NumberParameter, Parameter
from pytex.diffraction.xrd import HC_KEV_ANGSTROM, RadiationSpec

__all__ = [
    "MONOCHROMATIC",
    "RADIATION_OPTIONS",
    "radiation_from_request",
    "radiation_parameters",
]

#: The option that reads the wavelength from the request.
MONOCHROMATIC = "monochromatic"

_TUBES: Mapping[str, Callable[[], RadiationSpec]] = {
    "cu_ka_doublet": RadiationSpec.cu_ka_doublet,
    "cu_ka": RadiationSpec.cu_ka,
    "mo_ka_doublet": RadiationSpec.mo_ka_doublet,
    "co_ka_doublet": RadiationSpec.co_ka,
    "cr_ka_doublet": RadiationSpec.cr_ka,
    "fe_ka_doublet": RadiationSpec.fe_ka,
}

#: Every radiation a diffractometer operation offers, as (key, label, description).
RADIATION_OPTIONS: tuple[tuple[str, str, str], ...] = (
    ("cu_ka_doublet", "Cu Kα1/Kα2", "Common laboratory copper doublet."),
    ("cu_ka", "Cu Kα (single averaged line)", "One copper line without splitting."),
    ("mo_ka_doublet", "Mo Kα1/Kα2", "Short-wavelength molybdenum doublet."),
    ("co_ka_doublet", "Co Kα1/Kα2", "Reduces Fe fluorescence."),
    ("cr_ka_doublet", "Cr Kα1/Kα2", "Long wavelength; the classical choice for ferrite stress."),
    ("fe_ka_doublet", "Fe Kα1/Kα2", "Avoids Cr and Mn fluorescence."),
    (
        MONOCHROMATIC,
        "Monochromatic / synchrotron",
        "One wavelength of your choice, entered below: a synchrotron beamline or a "
        "Kα1-monochromated source. No Kα2 line is modelled.",
    ),
)


def radiation_parameters(
    *,
    help_text: str,
    default: str = "cu_ka_doublet",
    group: str | None = None,
    advanced: bool = False,
) -> tuple[Parameter, ...]:
    """Return the radiation choice and the monochromatic-beam controls.

    Parameters
    ----------
    help_text
        What the radiation does in *this* operation; the monochromatic option
        and the conversion from energy are explained after it.
    default
        The preselected option.
    group, advanced
        Placement of the choice. The wavelength sits beside it on one row; the
        polarization is always an advanced control, since only a polarized
        (synchrotron) beam needs it changed.
    """

    return (
        ChoiceParameter(
            name="radiation",
            label="Radiation",
            help_text=(
                help_text
                + " Choose Monochromatic / synchrotron for any single wavelength, entered beside "
                "it."
            ),
            options=RADIATION_OPTIONS,
            default=default,
            group=group,
            advanced=advanced,
            row="radiation",
        ),
        NumberParameter(
            name="wavelength_angstrom",
            label="Monochromatic wavelength",
            help_text=(
                "The wavelength of a monochromatic or synchrotron beam; needed, and used, only "
                f"when Radiation is Monochromatic / synchrotron. From the photon energy, λ[Å] = "
                f"{HC_KEV_ANGSTROM:.4f} / E[keV]: 30 keV is 0.4133 Å, 60 keV is 0.2066 Å."
            ),
            units="Å",
            required=False,
            minimum=0.01,
            maximum=5.0,
            exclusive_minimum=False,
            symbol="wavelength",
            group=group,
            advanced=advanced,
            row="radiation",
        ),
        NumberParameter(
            name="polarization_fraction",
            label="Beam polarization perpendicular to the scattering plane",
            help_text=(
                "Fraction f of the monochromatic beam polarized perpendicular to the scattering "
                "plane, which sets the polarization factor f + (1 − f) cos²2θ. 0.5 for an "
                "unpolarized beam; about 0.95 at a synchrotron with a vertical scattering "
                "plane, where the factor is almost 1; about 0.05 with a horizontal one. Used "
                "only when Radiation is Monochromatic / synchrotron."
            ),
            default=0.5,
            minimum=0.0,
            maximum=1.0,
            advanced=True,
            group=group,
        ),
    )


def radiation_from_request(request: Mapping[str, Any]) -> RadiationSpec:
    """Resolve a request's radiation controls into a :class:`RadiationSpec`."""

    key = str(request.get("radiation") or "cu_ka_doublet")
    if key == MONOCHROMATIC:
        wavelength = request.get("wavelength_angstrom")
        if wavelength is None:
            raise InvalidInputError(
                "A monochromatic beam needs its wavelength.",
                field="wavelength_angstrom",
                hint="Enter the wavelength in Å; λ[Å] = 12.3984 / E[keV].",
            )
        raw_fraction = request.get("polarization_fraction")
        fraction = 0.5 if raw_fraction is None else float(raw_fraction)
        try:
            return RadiationSpec.monochromatic(
                float(wavelength), polarization_perpendicular_fraction=fraction
            )
        except ValueError as error:
            raise InvalidInputError(str(error), field="wavelength_angstrom") from error
    try:
        return _TUBES[key]()
    except KeyError:
        raise InvalidInputError(
            f"Unknown radiation {key!r}.",
            field="radiation",
            hint="Choose one of the listed tube lines, or Monochromatic / synchrotron.",
        ) from None
