"""How a phase enters the application, and the catalogue it starts with.

Every panel needs a material. The library's :class:`~pytex.core.lattice.Phase`
is the right object to *compute* with, but it is not the right object to put on a
wire: it carries reference frames, symmetry operator arrays, and provenance that
a browser has no use for and cannot construct. :class:`PhaseSpec` is the wire
form — six cell parameters, a point group, optionally a space group and an atomic
basis — and :meth:`PhaseSpec.to_phase` is the one place it becomes a real phase.

Two properties are deliberate.

**No optional dependency is required for the catalogue or manual path.** A user
can describe any phase of any space group by typing numbers, and the built-in
catalogue below is literal parameters in Python rather than parsed CIF files.
The shared GUI phase control also accepts CIF text through
:func:`phase_from_request` where pymatgen is installed; without it the rest of
the application remains fully usable — which matters because the deployment
target is an intranet host that may never have seen PyPI.

**Cell parameters are checked against the crystal system.** Asking for a cubic
point group with ``a != b`` is a mistake, and catching it at construction with a
sentence naming the offending parameter is worth far more than a silently
non-cubic "cubic" phase propagating into an indexing result.
"""

from __future__ import annotations

import logging
import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.uploads import uploaded_name_and_text
from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import AtomicSite, Lattice, Phase, SpaceGroupSpec, UnitCell
from pytex.core.point_groups import PointGroup, all_point_group_symbols
from pytex.core.provenance import ProvenanceRecord
from pytex.core.symmetry import SymmetrySpec

__all__ = [
    "APP_CRYSTAL_FRAME",
    "BUILTIN_PHASES",
    "PhaseSpec",
    "SiteSpec",
    "builtin_phase",
    "list_builtin_phases",
    "phase_from_request",
]

#: The crystal frame every application-constructed phase is stated in.
#:
#: One frame, shared by every phase the app builds, so two phases entering the
#: same calculation cannot disagree about what ``a``, ``b`` and ``c`` mean. A
#: user who needs a different crystal-frame convention is doing library work,
#: not application work.
APP_CRYSTAL_FRAME = ReferenceFrame(
    name="crystal",
    domain=FrameDomain.CRYSTAL,
    axes=("a", "b", "c"),
    handedness=Handedness.RIGHT,
)

_ANGLE_TOLERANCE_DEG = 1e-6
_LENGTH_TOLERANCE = 1e-9
_LOGGER = logging.getLogger("pytex.app.phases")


@dataclass(frozen=True)
class SiteSpec:
    """One atom of the basis, in fractional coordinates.

    Attributes
    ----------
    species : str
        Element symbol, e.g. ``"Fe"``. Sets the scattering factor, so it must be
        an element rather than a free-text label.
    x, y, z : float
        Fractional coordinates along ``a``, ``b`` and ``c``.
    occupancy : float
        Site occupancy in ``(0, 1]``.
    label : str, optional
        CIF-style site label. Defaults to the species plus an index.
    b_iso : float, optional
        Isotropic displacement parameter in Å².
    """

    species: str
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    label: str | None = None
    b_iso: float | None = None

    def to_json(self) -> dict[str, Any]:
        """Return the wire form of this site."""

        payload: dict[str, Any] = {
            "species": self.species,
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "occupancy": self.occupancy,
        }
        if self.label is not None:
            payload["label"] = self.label
        if self.b_iso is not None:
            payload["b_iso"] = self.b_iso
        return payload

    @classmethod
    def from_json(cls, payload: Mapping[str, Any], *, index: int = 0) -> SiteSpec:
        """Build a site from its wire form, reporting bad input by field."""

        species = payload.get("species")
        if not isinstance(species, str) or not species.strip():
            raise InvalidInputError(
                "Each atomic site needs an element symbol.",
                field=f"sites[{index}].species",
                hint="For example: Fe, Ni, O.",
            )
        coordinates = tuple(
            _require_float(payload.get(axis, 0.0), field=f"sites[{index}].{axis}")
            for axis in ("x", "y", "z")
        )
        occupancy = _require_float(payload.get("occupancy", 1.0), field=f"sites[{index}].occupancy")
        if not 0.0 < occupancy <= 1.0:
            raise InvalidInputError(
                f"Site occupancy must lie in (0, 1]; got {occupancy:g}.",
                field=f"sites[{index}].occupancy",
            )
        b_iso_raw = payload.get("b_iso")
        b_iso = (
            None if b_iso_raw is None else _require_float(b_iso_raw, field=f"sites[{index}].b_iso")
        )
        label = payload.get("label")
        return cls(
            species=species.strip(),
            x=coordinates[0],
            y=coordinates[1],
            z=coordinates[2],
            occupancy=occupancy,
            label=label if isinstance(label, str) and label.strip() else None,
            b_iso=b_iso,
        )

    def to_atomic_site(self, *, index: int) -> AtomicSite:
        """Convert to the library's site object."""

        return AtomicSite(
            label=self.label or f"{self.species}{index + 1}",
            species=self.species,
            fractional_coordinates=np.asarray([self.x, self.y, self.z], dtype=float),
            occupancy=self.occupancy,
            b_iso=self.b_iso,
        )


@dataclass(frozen=True)
class PhaseSpec:
    """A phase as the user describes it: parameters, symmetry, optional basis.

    Purpose
    -------
    The single input every application operation takes when it needs a
    material. JSON-round-trippable in both directions, so a saved analysis
    carries the exact phase it was computed with rather than a name that may
    mean something different next year.

    When to use
    -----------
    Construct one from a request payload with :meth:`from_json`, or take a
    catalogue entry with :func:`builtin_phase`. Call :meth:`to_phase` once, at
    the top of a service, and pass the resulting :class:`Phase` down.

    Attributes
    ----------
    name : str
        Display name, free text.
    a, b, c : float
        Cell edge lengths in Å.
    alpha, beta, gamma : float
        Cell angles in degrees.
    point_group : str
        Hermann-Mauguin symbol, one of the 32 crystallographic point groups.
        This — not the space group — is what symmetry operations come from.
    space_group_symbol, space_group_number : str, int, optional
        Supplied together or not at all. The *centring letter* of the symbol is
        what determines systematic absences, so a diffraction result computed
        without a space group will list reflections that a real pattern does not
        show. The application warns rather than guesses.
    sites : tuple of SiteSpec
        Atomic basis. Required for structure factors and for a crystal viewer
        that shows atoms; without it only lattice geometry is available.
    source : str, optional
        Where the numbers came from, carried into provenance and shown in the
        UI. Catalogue entries cite a reference here.

    Returns
    -------
    PhaseSpec
        A validated specification; construction raises
        :class:`~pytex.app.errors.InvalidInputError` on any inconsistency.

    Examples
    --------
    >>> spec = PhaseSpec(name="Ni", a=3.52387, b=3.52387, c=3.52387,
    ...                  alpha=90.0, beta=90.0, gamma=90.0, point_group="m-3m")
    >>> phase = spec.to_phase()
    >>> round(float(phase.lattice.a), 5)
    3.52387
    """

    name: str
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    point_group: str
    space_group_symbol: str | None = None
    space_group_number: int | None = None
    sites: tuple[SiteSpec, ...] = ()
    source: str | None = None

    def __post_init__(self) -> None:
        symbol = self.point_group.strip()
        if symbol not in all_point_group_symbols():
            raise InvalidInputError(
                f"{symbol!r} is not one of the 32 crystallographic point groups.",
                field="point_group",
                hint="Choose from: " + ", ".join(all_point_group_symbols()) + ".",
            )
        object.__setattr__(self, "point_group", symbol)
        for field_name in ("a", "b", "c"):
            value = float(getattr(self, field_name))
            if not math.isfinite(value) or value <= 0.0:
                raise InvalidInputError(
                    f"Cell edge {field_name} must be a positive length; got {value:g} Å.",
                    field=field_name,
                )
            object.__setattr__(self, field_name, value)
        for field_name in ("alpha", "beta", "gamma"):
            value = float(getattr(self, field_name))
            if not math.isfinite(value) or not 0.0 < value < 180.0:
                raise InvalidInputError(
                    f"Cell angle {field_name} must lie strictly between 0 and 180 degrees; "
                    f"got {value:g}.",
                    field=field_name,
                )
            object.__setattr__(self, field_name, value)
        if (self.space_group_symbol is None) != (self.space_group_number is None):
            raise InvalidInputError(
                "Give the space-group symbol and number together, or neither.",
                field="space_group_symbol",
                hint="For example: Fm-3m with number 225.",
            )
        if self.space_group_number is not None and not 1 <= self.space_group_number <= 230:
            raise InvalidInputError(
                f"Space-group number must lie in [1, 230]; got {self.space_group_number}.",
                field="space_group_number",
            )
        object.__setattr__(self, "sites", tuple(self.sites))
        self._check_crystal_system()
        self._check_space_group_system()

    @property
    def crystal_system(self) -> str:
        """Crystal system implied by the point group."""

        return str(PointGroup.from_symbol(self.point_group).crystal_system)

    def _check_crystal_system(self) -> None:
        """Reject cell parameters that contradict the point group.

        The metric constraints of each crystal system are not conventions to be
        relaxed: a tetragonal point group applied to a cell with ``a != b``
        produces symmetry operations that do not map the lattice to itself, and
        every downstream orbit, multiplicity, and angle is then wrong. Better to
        refuse at the point of entry, naming the parameter to fix.
        """

        system = self.crystal_system
        equal_edges: tuple[tuple[str, str], ...]
        fixed_angles: tuple[tuple[str, float], ...]
        if system == "cubic":
            equal_edges = (("a", "b"), ("a", "c"))
            fixed_angles = (("alpha", 90.0), ("beta", 90.0), ("gamma", 90.0))
        elif system == "tetragonal":
            equal_edges = (("a", "b"),)
            fixed_angles = (("alpha", 90.0), ("beta", 90.0), ("gamma", 90.0))
        elif system in {"hexagonal", "trigonal"}:
            # Trigonal phases are accepted in the hexagonal setting, which is
            # what CIFs, International Tables, and every downstream
            # Miller-Bravais conversion in PyTex assume. The rhombohedral
            # setting is a different cell, not a looser check.
            equal_edges = (("a", "b"),)
            fixed_angles = (("alpha", 90.0), ("beta", 90.0), ("gamma", 120.0))
        elif system == "orthorhombic":
            equal_edges = ()
            fixed_angles = (("alpha", 90.0), ("beta", 90.0), ("gamma", 90.0))
        elif system == "monoclinic":
            equal_edges = ()
            fixed_angles = (("alpha", 90.0), ("gamma", 90.0))
        else:  # triclinic imposes nothing
            equal_edges = ()
            fixed_angles = ()
        for first, second in equal_edges:
            if abs(getattr(self, first) - getattr(self, second)) > _LENGTH_TOLERANCE:
                raise InvalidInputError(
                    f"A {system} phase requires {first} = {second}; got "
                    f"{getattr(self, first):g} Å and {getattr(self, second):g} Å.",
                    field=second,
                    hint=f"Either correct {second} or choose a lower-symmetry point group.",
                )
        for angle_name, expected in fixed_angles:
            if abs(getattr(self, angle_name) - expected) > _ANGLE_TOLERANCE_DEG:
                raise InvalidInputError(
                    f"A {system} phase requires {angle_name} = {expected:g}°; got "
                    f"{getattr(self, angle_name):g}°.",
                    field=angle_name,
                    hint=f"Either correct {angle_name} or choose a lower-symmetry point group.",
                )

    def _check_space_group_system(self) -> None:
        """Reject a space group that belongs to a different crystal system.

        The space group is not decoration: its centring letter decides which
        reflections are systematically absent, so a cubic F number carried onto
        a tetragonal phase silently deletes most of its pattern. The pairing
        arises the moment a user edits a catalogue phase — the cell and the
        point group change, and the space group inherited from the entry they
        started from does not.

        Only the system is checked, because that is what the number alone can
        settle. The ranges are the standard sequence of International Tables
        Volume A, where space groups are numbered in order of increasing
        symmetry.
        """

        if self.space_group_number is None:
            return
        ranges = {
            "triclinic": (1, 2),
            "monoclinic": (3, 15),
            "orthorhombic": (16, 74),
            "tetragonal": (75, 142),
            "trigonal": (143, 167),
            "hexagonal": (168, 194),
            "cubic": (195, 230),
        }
        system = self.crystal_system
        low, high = ranges[system]
        if not low <= self.space_group_number <= high:
            actual = next(
                (name for name, (a, b) in ranges.items() if a <= self.space_group_number <= b),
                "another system",
            )
            symbol = self.space_group_symbol or f"number {self.space_group_number}"
            raise InvalidInputError(
                f"Space group {symbol} (number {self.space_group_number}) is {actual}, but "
                f"point group {self.point_group} is {system}.",
                field="space_group_number",
                hint=(
                    f"A {system} phase takes a space group numbered {low}-{high}. Clear the "
                    "space group to compute without systematic absences, or correct it to one "
                    "of this system."
                ),
            )

    @property
    def uses_miller_bravais(self) -> bool:
        """True when four-index notation is the natural one for this phase."""

        return self.crystal_system in {"hexagonal", "trigonal"}

    @property
    def has_structure(self) -> bool:
        """True when an atomic basis is present, so intensities are available."""

        return bool(self.sites)

    def to_json(self) -> dict[str, Any]:
        """Return the wire form, suitable for saving beside a result."""

        payload: dict[str, Any] = {
            "name": self.name,
            "a": self.a,
            "b": self.b,
            "c": self.c,
            "alpha": self.alpha,
            "beta": self.beta,
            "gamma": self.gamma,
            "point_group": self.point_group,
            "crystal_system": self.crystal_system,
        }
        if self.space_group_symbol is not None:
            payload["space_group_symbol"] = self.space_group_symbol
            payload["space_group_number"] = self.space_group_number
        if self.sites:
            payload["sites"] = [site.to_json() for site in self.sites]
        if self.source is not None:
            payload["source"] = self.source
        return payload

    @classmethod
    def from_phase(cls, phase: Phase, *, source: str | None = None) -> PhaseSpec:
        """Build the application wire form from a canonical phase.

        Purpose
        -------
        Keep imported structures on the same service contract as catalogue and
        manually entered phases. This is the application boundary from a fully
        semantic :class:`~pytex.core.lattice.Phase` to its JSON-ready form; it
        does not expose a parser-specific structure object.

        Parameters
        ----------
        phase : Phase
            A canonical phase carrying lattice, symmetry, and optionally an
            atomic unit cell and space group.
        source : str, optional
            User-facing provenance text to retain in results and exports.

        Returns
        -------
        PhaseSpec
            A validated phase specification preserving the imported structure.
        """

        lattice = phase.lattice
        sites = (
            tuple(
                SiteSpec(
                    species=site.species,
                    x=float(site.fractional_coordinates[0]),
                    y=float(site.fractional_coordinates[1]),
                    z=float(site.fractional_coordinates[2]),
                    occupancy=float(site.occupancy),
                    label=site.label,
                    b_iso=site.b_iso,
                )
                for site in phase.unit_cell.sites
            )
            if phase.unit_cell is not None
            else ()
        )
        return cls(
            name=phase.name,
            a=lattice.a,
            b=lattice.b,
            c=lattice.c,
            alpha=lattice.alpha_deg,
            beta=lattice.beta_deg,
            gamma=lattice.gamma_deg,
            point_group=phase.symmetry.point_group,
            space_group_symbol=phase.space_group_symbol,
            space_group_number=phase.space_group_number,
            sites=sites,
            source=source,
        )

    @classmethod
    def from_json(cls, payload: Mapping[str, Any]) -> PhaseSpec:
        """Build a specification from a request payload.

        Accepts either a full description or ``{"builtin": "ni_fcc"}``, so the
        frontend can send a catalogue reference without duplicating its numbers.
        Both forms may be combined: a ``builtin`` key with further keys applies
        those as overrides, which is how "start from nickel and change ``a``"
        works without the user retyping the structure.
        """

        if not isinstance(payload, Mapping):
            raise InvalidInputError("The phase must be an object.", field="phase")
        overrides = dict(payload)
        builtin_id = overrides.pop("builtin", None)
        base: PhaseSpec | None = None
        if builtin_id is not None:
            if not isinstance(builtin_id, str):
                raise InvalidInputError(
                    "The built-in phase identifier must be a string.", field="phase.builtin"
                )
            base = builtin_phase(builtin_id)
            if not overrides:
                return base
        overrides.pop("crystal_system", None)
        sites_payload = overrides.pop("sites", None)
        sites: tuple[SiteSpec, ...] | None = None
        if sites_payload is not None:
            if not isinstance(sites_payload, Sequence) or isinstance(sites_payload, str):
                raise InvalidInputError("Atomic sites must be a list.", field="phase.sites")
            sites = tuple(
                SiteSpec.from_json(entry, index=index)
                for index, entry in enumerate(sites_payload)
                if isinstance(entry, Mapping)
            )
        if base is not None:
            known = {
                "name",
                "a",
                "b",
                "c",
                "alpha",
                "beta",
                "gamma",
                "point_group",
                "space_group_symbol",
                "space_group_number",
                "source",
            }
            unexpected = sorted(set(overrides) - known)
            if unexpected:
                raise InvalidInputError(
                    f"Unknown phase field(s): {', '.join(unexpected)}.",
                    field="phase",
                    hint=f"Accepted fields: {', '.join(sorted(known | {'sites'}))}.",
                )
            changes: dict[str, Any] = {key: overrides[key] for key in overrides}
            for key in ("a", "b", "c", "alpha", "beta", "gamma"):
                if key in changes:
                    changes[key] = _require_float(changes[key], field=f"phase.{key}")
            if sites is not None:
                changes["sites"] = sites
            return replace(base, **changes)
        missing = [
            key
            for key in ("a", "b", "c", "alpha", "beta", "gamma", "point_group")
            if key not in overrides
        ]
        if missing:
            raise InvalidInputError(
                f"The phase is missing: {', '.join(missing)}.",
                field="phase",
                hint="Give the six cell parameters and a point group, or name a built-in phase.",
            )
        point_group = overrides["point_group"]
        if not isinstance(point_group, str):
            raise InvalidInputError(
                "The point group must be a Hermann-Mauguin symbol.", field="phase.point_group"
            )
        number = overrides.get("space_group_number")
        return cls(
            name=str(overrides.get("name") or "phase"),
            a=_require_float(overrides["a"], field="phase.a"),
            b=_require_float(overrides["b"], field="phase.b"),
            c=_require_float(overrides["c"], field="phase.c"),
            alpha=_require_float(overrides["alpha"], field="phase.alpha"),
            beta=_require_float(overrides["beta"], field="phase.beta"),
            gamma=_require_float(overrides["gamma"], field="phase.gamma"),
            point_group=point_group,
            space_group_symbol=(
                str(overrides["space_group_symbol"])
                if overrides.get("space_group_symbol") is not None
                else None
            ),
            space_group_number=int(number) if number is not None else None,
            sites=sites or (),
            source=str(overrides["source"]) if overrides.get("source") is not None else None,
        )

    def to_phase(self) -> Phase:
        """Build the library phase this specification describes.

        Returns
        -------
        Phase
            Stated in :data:`APP_CRYSTAL_FRAME`, carrying a unit cell when an
            atomic basis was given and a space group when one was named.
        """

        lattice = Lattice(
            self.a,
            self.b,
            self.c,
            self.alpha,
            self.beta,
            self.gamma,
            crystal_frame=APP_CRYSTAL_FRAME,
        )
        symmetry = SymmetrySpec.from_point_group(
            self.point_group, reference_frame=APP_CRYSTAL_FRAME
        )
        unit_cell = (
            UnitCell(
                lattice=lattice,
                sites=tuple(
                    site.to_atomic_site(index=index) for index, site in enumerate(self.sites)
                ),
            )
            if self.sites
            else None
        )
        space_group = (
            SpaceGroupSpec(
                symbol=self.space_group_symbol,
                number=int(self.space_group_number),
                reference_frame=APP_CRYSTAL_FRAME,
                crystal_system=self.crystal_system,
            )
            if self.space_group_symbol is not None and self.space_group_number is not None
            else None
        )
        return Phase(
            self.name,
            lattice=lattice,
            symmetry=symmetry,
            crystal_frame=APP_CRYSTAL_FRAME,
            unit_cell=unit_cell,
            space_group=space_group,
            space_group_symbol=self.space_group_symbol,
            space_group_number=self.space_group_number,
        )


def _require_float(value: Any, *, field: str) -> float:
    """Coerce a JSON value to a finite float, or raise with the field named."""

    if isinstance(value, bool) or not isinstance(value, int | float | str):
        raise InvalidInputError(f"{field} must be a number.", field=field)
    try:
        number = float(value)
    except ValueError as exc:
        raise InvalidInputError(f"{field} must be a number; got {value!r}.", field=field) from exc
    if not math.isfinite(number):
        raise InvalidInputError(f"{field} must be finite.", field=field)
    return number


def _fcc_sites(species: str) -> tuple[SiteSpec, ...]:
    """The four face-centred positions of a monatomic cubic-F cell."""

    return (
        SiteSpec(species, 0.0, 0.0, 0.0),
        SiteSpec(species, 0.0, 0.5, 0.5),
        SiteSpec(species, 0.5, 0.0, 0.5),
        SiteSpec(species, 0.5, 0.5, 0.0),
    )


def _bcc_sites(species: str) -> tuple[SiteSpec, ...]:
    """The two body-centred positions of a monatomic cubic-I cell."""

    return (SiteSpec(species, 0.0, 0.0, 0.0), SiteSpec(species, 0.5, 0.5, 0.5))


def _hcp_sites(species: str) -> tuple[SiteSpec, ...]:
    """The two positions of the hexagonal close-packed basis."""

    return (
        SiteSpec(species, 1.0 / 3.0, 2.0 / 3.0, 0.25),
        SiteSpec(species, 2.0 / 3.0, 1.0 / 3.0, 0.75),
    )


#: The phases the application offers out of the box.
#:
#: Literal parameters with a cited source, not parsed files: the catalogue must
#: work on an installation with no optional dependencies and no data directory.
#: Room-temperature values unless stated otherwise.
BUILTIN_PHASES: dict[str, PhaseSpec] = {
    "ni_fcc": PhaseSpec(
        name="Nickel (fcc)",
        a=3.52387,
        b=3.52387,
        c=3.52387,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Fm-3m",
        space_group_number=225,
        sites=_fcc_sites("Ni"),
        source="Wyckoff, Crystal Structures Vol. 1 (1963); COD 9008476.",
    ),
    "al_fcc": PhaseSpec(
        name="Aluminium (fcc)",
        a=4.0495,
        b=4.0495,
        c=4.0495,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Fm-3m",
        space_group_number=225,
        sites=_fcc_sites("Al"),
        source="Wyckoff, Crystal Structures Vol. 1 (1963).",
    ),
    "cu_fcc": PhaseSpec(
        name="Copper (fcc)",
        a=3.6149,
        b=3.6149,
        c=3.6149,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Fm-3m",
        space_group_number=225,
        sites=_fcc_sites("Cu"),
        source="Wyckoff, Crystal Structures Vol. 1 (1963).",
    ),
    "austenite_fcc": PhaseSpec(
        name="Austenite (fcc Fe)",
        a=3.60,
        b=3.60,
        c=3.60,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Fm-3m",
        space_group_number=225,
        sites=_fcc_sites("Fe"),
        source="Conventional austenite cell parameter used in OR analysis.",
    ),
    "fe_bcc": PhaseSpec(
        name="Ferrite (bcc Fe)",
        a=2.8665,
        b=2.8665,
        c=2.8665,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Im-3m",
        space_group_number=229,
        sites=_bcc_sites("Fe"),
        source="Wyckoff, Crystal Structures Vol. 1 (1963); COD 9008536.",
    ),
    "zr_bcc_beta": PhaseSpec(
        name="Zirconium (bcc, beta at 863 °C)",
        a=3.6090,
        b=3.6090,
        c=3.6090,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Im-3m",
        space_group_number=229,
        sites=_bcc_sites("Zr"),
        source=(
            "Zuzek et al., Bull. Alloy Phase Diagrams 11 (1990) 385; "
            "Maimaitiyili et al., J. Synchrotron Rad. 22 (2015) 995, "
            "doi:10.1107/S1600577515009054 (a = 3.6090 Å at 863 °C)."
        ),
    ),
    "w_bcc": PhaseSpec(
        name="Tungsten (bcc)",
        a=3.1652,
        b=3.1652,
        c=3.1652,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Im-3m",
        space_group_number=229,
        sites=_bcc_sites("W"),
        source="Wyckoff, Crystal Structures Vol. 1 (1963).",
    ),
    "si_diamond": PhaseSpec(
        name="Silicon (diamond cubic)",
        a=5.43102,
        b=5.43102,
        c=5.43102,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Fd-3m",
        space_group_number=227,
        sites=(
            SiteSpec("Si", 0.0, 0.0, 0.0),
            SiteSpec("Si", 0.0, 0.5, 0.5),
            SiteSpec("Si", 0.5, 0.0, 0.5),
            SiteSpec("Si", 0.5, 0.5, 0.0),
            SiteSpec("Si", 0.25, 0.25, 0.25),
            SiteSpec("Si", 0.25, 0.75, 0.75),
            SiteSpec("Si", 0.75, 0.25, 0.75),
            SiteSpec("Si", 0.75, 0.75, 0.25),
        ),
        source="CODATA lattice parameter of silicon; COD 9008565.",
    ),
    "diamond": PhaseSpec(
        name="Diamond (C)",
        a=3.56679,
        b=3.56679,
        c=3.56679,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Fd-3m",
        space_group_number=227,
        sites=(
            SiteSpec("C", 0.0, 0.0, 0.0),
            SiteSpec("C", 0.0, 0.5, 0.5),
            SiteSpec("C", 0.5, 0.0, 0.5),
            SiteSpec("C", 0.5, 0.5, 0.0),
            SiteSpec("C", 0.25, 0.25, 0.25),
            SiteSpec("C", 0.25, 0.75, 0.75),
            SiteSpec("C", 0.75, 0.25, 0.75),
            SiteSpec("C", 0.75, 0.75, 0.25),
        ),
        source="Straumanis & Aka, J. Am. Chem. Soc. 73 (1951) 5643; COD 9012304.",
    ),
    "nacl": PhaseSpec(
        name="Halite (NaCl)",
        a=5.6402,
        b=5.6402,
        c=5.6402,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="m-3m",
        space_group_symbol="Fm-3m",
        space_group_number=225,
        sites=(
            SiteSpec("Na", 0.0, 0.0, 0.0),
            SiteSpec("Na", 0.0, 0.5, 0.5),
            SiteSpec("Na", 0.5, 0.0, 0.5),
            SiteSpec("Na", 0.5, 0.5, 0.0),
            SiteSpec("Cl", 0.5, 0.5, 0.5),
            SiteSpec("Cl", 0.5, 0.0, 0.0),
            SiteSpec("Cl", 0.0, 0.5, 0.0),
            SiteSpec("Cl", 0.0, 0.0, 0.5),
        ),
        source="Wyckoff, Crystal Structures Vol. 1 (1963).",
    ),
    "ti_hcp": PhaseSpec(
        name="Titanium (hcp, alpha)",
        a=2.9508,
        b=2.9508,
        c=4.6855,
        alpha=90.0,
        beta=90.0,
        gamma=120.0,
        point_group="6/mmm",
        space_group_symbol="P6_3/mmc",
        space_group_number=194,
        sites=_hcp_sites("Ti"),
        source="Wyckoff, Crystal Structures Vol. 1 (1963).",
    ),
    "zr_hcp": PhaseSpec(
        name="Zirconium (hcp, alpha)",
        a=3.232,
        b=3.232,
        c=5.147,
        alpha=90.0,
        beta=90.0,
        gamma=120.0,
        point_group="6/mmm",
        space_group_symbol="P6_3/mmc",
        space_group_number=194,
        sites=_hcp_sites("Zr"),
        source="Wyckoff, Crystal Structures Vol. 1 (1963); COD 9008523.",
    ),
    "mg_hcp": PhaseSpec(
        name="Magnesium (hcp)",
        a=3.2094,
        b=3.2094,
        c=5.2108,
        alpha=90.0,
        beta=90.0,
        gamma=120.0,
        point_group="6/mmm",
        space_group_symbol="P6_3/mmc",
        space_group_number=194,
        sites=_hcp_sites("Mg"),
        source="Wyckoff, Crystal Structures Vol. 1 (1963).",
    ),
    "alpha_u": PhaseSpec(
        name="Alpha uranium (orthorhombic)",
        a=2.8537,
        b=5.8695,
        c=4.9548,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        point_group="mmm",
        space_group_symbol="Cmcm",
        space_group_number=63,
        sites=(
            SiteSpec("U", 0.0, 0.1025, 0.25),
            SiteSpec("U", 0.0, 0.8975, 0.75),
            SiteSpec("U", 0.5, 0.6025, 0.25),
            SiteSpec("U", 0.5, 0.3975, 0.75),
        ),
        source="Jacob & Warren, J. Am. Chem. Soc. 59 (1937) 2588.",
    ),
    "quartz_alpha": PhaseSpec(
        name="Alpha quartz (SiO2)",
        a=4.9134,
        b=4.9134,
        c=5.4052,
        alpha=90.0,
        beta=90.0,
        gamma=120.0,
        point_group="32",
        space_group_symbol="P3_221",
        space_group_number=154,
        sites=(
            SiteSpec("Si", 0.4697, 0.0, 0.0),
            SiteSpec("Si", 0.0, 0.4697, 2.0 / 3.0),
            SiteSpec("Si", 0.5303, 0.5303, 1.0 / 3.0),
            SiteSpec("O", 0.4135, 0.2669, 0.1191),
            SiteSpec("O", 0.2669, 0.4135, 0.5476),
            SiteSpec("O", 0.7331, 0.1466, 0.7857),
            SiteSpec("O", 0.5865, 0.8534, 0.2143),
            SiteSpec("O", 0.8534, 0.5865, 0.4524),
            SiteSpec("O", 0.1466, 0.7331, 0.8809),
        ),
        source="Levien, Prewitt & Weidner, Am. Mineral. 65 (1980) 920.",
    ),
}


def list_builtin_phases() -> tuple[dict[str, Any], ...]:
    """Return catalogue entries for the phase picker.

    Each entry carries the identifier, display name, crystal system, space
    group, and source, which is everything the picker shows before a phase is
    chosen. The full specification arrives with the selection.
    """

    entries: list[dict[str, Any]] = []
    for identifier, spec in BUILTIN_PHASES.items():
        entry = spec.to_json()
        entry["id"] = identifier
        entries.append(entry)
    return tuple(entries)


def builtin_phase(identifier: str) -> PhaseSpec:
    """Return one catalogue entry by identifier."""

    try:
        return BUILTIN_PHASES[identifier]
    except KeyError:
        raise InvalidInputError(
            f"No built-in phase named {identifier!r}.",
            field="phase.builtin",
            hint="Available: " + ", ".join(sorted(BUILTIN_PHASES)) + ".",
        ) from None


def phase_from_request(payload: Mapping[str, Any] | None) -> tuple[PhaseSpec, Phase]:
    """Resolve a request's phase payload into both forms at once.

    Services need the specification (to echo into the result and the export
    manifest) and the phase (to compute with). Returning both from one call
    keeps them from being derived separately and drifting.
    """

    if payload is None:
        raise InvalidInputError(
            "This operation needs a phase.",
            field="phase",
            hint=(
                "Choose a built-in phase, load a .cif file, or enter cell parameters and a "
                "point group."
            ),
        )
    cif_payload = payload.get("cif")
    if cif_payload is not None:
        if set(payload) != {"cif"}:
            unexpected = ", ".join(sorted(str(key) for key in payload if key != "cif"))
            raise InvalidInputError(
                f"A CIF phase cannot also override phase field(s): {unexpected}.",
                field="phase",
                hint="Load the CIF by itself, then edit the imported cell as a custom phase.",
            )
        name, text = uploaded_name_and_text(cif_payload, field="phase", suffixes=(".cif",))
        provenance = ProvenanceRecord(
            source_system="cif",
            source_identifier=name,
            metadata={"reader": "pymatgen.Structure.from_str", "format": "cif"},
        )
        try:
            with warnings.catch_warnings(record=True) as parser_notes:
                warnings.simplefilter("always")
                phase = Phase.from_cif_string(
                    text,
                    crystal_frame=APP_CRYSTAL_FRAME,
                    provenance=provenance,
                )
        except (IndexError, KeyError, RuntimeError, TypeError, ValueError) as error:
            raise InvalidInputError(
                f"PyTex could not read {name} as a crystallographic structure: {error}",
                field="phase",
                hint=(
                    "Open a valid text CIF containing cell parameters and atomic sites. "
                    "If another crystallography program also rejects it, repair or export it again."
                ),
            ) from error
        for note in dict.fromkeys(str(item.message) for item in parser_notes):
            _LOGGER.warning("CIF parser note for %s: %s", name, note)
        spec = PhaseSpec.from_phase(
            phase,
            source=f"CIF file {name}; parsed by pymatgen into canonical PyTex phase semantics.",
        )
        return spec, phase
    spec = PhaseSpec.from_json(payload)
    return spec, spec.to_phase()
