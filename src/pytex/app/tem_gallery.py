"""The practice-pattern gallery: three SAED plates with known answers.

Why a gallery
-------------
The indexing half of the TEM panel used to require a micrograph. That is exactly
backwards for learning it and for testing it: a researcher wanting to check that
their calibration workflow is right, or a student meeting reciprocal space for
the first time, has to find a pattern before they can do anything at all — and
whatever they find comes with no answer key, so a mistake in the workflow looks
identical to a mistake in the crystallography.

These three entries remove that barrier. Each is a real crystallographic
calculation — the same one :mod:`pytex.tem.synthetic` performs for any phase and
axis — pinned to a specific material, a specific zone axis, and a specific,
stated instrument setting, so the pattern on screen is what that microscope would
record from that specimen. Because the pattern was built from a zone axis, the
answer is known, and the panel can tell a user whether they indexed it correctly.

The three were chosen to be the three cases a microscopist meets first, and to be
different from one another in the way that matters:

``fcc_al_001``
    Face-centred cubic down [001]: a square, four-fold pattern where every
    reflection has unmixed indices. The reference case for the ratio-and-angle
    method.
``bcc_fe_110``
    Body-centred cubic down [110]: a centred rectangle, two-fold, with
    ``h + k + l`` even. The pattern whose rectangle is routinely mistaken for a
    square when the camera constant is wrong.
``hcp_zr_2-1-10``
    Hexagonal close-packed down [2̄110]: a rectangle whose *aspect ratio is the
    axial ratio*, so it measures ``c/a`` directly and is the standard way to tell
    one hcp metal from another.

Instrument setting
------------------
The accelerating voltage and camera length are **not** properties of an entry.
They belong to the microscope, one set of them applies to whatever specimen is
in the column, and the panel exposes them as controls shared by all three
plates — defaulting to 200 kV and 400 mm. The camera constant is then computed
as ``L·λ`` with ``λ`` from :func:`electron_wavelength_angstrom`, rather than
typed in, which is the relation the calibration field exists to teach: shorten
the camera and more of reciprocal space fits on the plate.

What *is* per-entry is everything the specimen and the exposure decide: the
phase, the zone axis, where the beam falls on the detector, how the pattern is
rolled about it, and the centroiding scatter.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pytex.app.errors import InvalidInputError
from pytex.app.phases import builtin_phase

__all__ = [
    "GALLERY",
    "GalleryEntry",
    "SuggestedTarget",
    "gallery_entry",
    "gallery_options",
]


@dataclass(frozen=True, slots=True)
class SuggestedTarget:
    """A zone axis worth tilting to from a gallery pattern, and why.

    Attributes
    ----------
    indices : tuple of int
        The target, in the phase's three-index basis.
    reason : str
        What the trip buys — the reason a microscopist would make it.
    """

    indices: tuple[int, ...]
    reason: str

    def to_json(self) -> dict[str, Any]:
        """The suggestion as JSON-ready data."""

        return {"indices": list(self.indices), "reason": self.reason}


@dataclass(frozen=True, slots=True)
class GalleryEntry:
    """One practice pattern: a phase, an axis, and the microscope that recorded it.

    Attributes
    ----------
    identifier : str
        Stable id used by the manifest and the browser.
    title, summary : str
        What the entry is, for the gallery card.
    teaches : str
        What a user should take away from indexing it. The gallery is a teaching
        surface, and an entry without a lesson is a picture.
    phase_key : str
        Catalogue identifier of the phase, from :data:`pytex.app.phases.BUILTIN_PHASES`.
    zone_axis : tuple of int
        Beam direction in the phase's three-index basis.
    width_px, height_px, pixel_size_mm : int, int, float
        The detector raster.
    centre_px : tuple of float
        Where the transmitted beam falls. Deliberately off the geometric centre
        on some entries, because on a real plate it is.
    in_plane_rotation_deg : float
        Roll of the pattern about the beam. Non-zero everywhere: a pattern
        aligned with the detector axes would let a wrong workflow pass.
    jitter_px : float
        Centroiding error added to each spot, so indexing residuals are
        realistic.
    rng_seed : int
        Fixes the jitter, so the entry looks the same every time it is opened.
    targets : tuple of SuggestedTarget
        Where to go next, and why.
    """

    identifier: str
    title: str
    summary: str
    teaches: str
    phase_key: str
    zone_axis: tuple[int, ...]
    width_px: int = 1024
    height_px: int = 1024
    pixel_size_mm: float = 0.024
    centre_px: tuple[float, float] = (512.0, 512.0)
    in_plane_rotation_deg: float = 0.0
    jitter_px: float = 0.9
    rng_seed: int = 20260814
    targets: tuple[SuggestedTarget, ...] = ()

    def phase_spec(self) -> Any:
        """The catalogue phase this entry is built from."""

        return builtin_phase(self.phase_key)


GALLERY: tuple[GalleryEntry, ...] = (
    GalleryEntry(
        identifier="fcc_al_001",
        title="Aluminium (fcc) down [001]",
        summary="A square four-fold pattern; the reference case for ratio-and-angle indexing.",
        teaches=(
            "Every reflection here has unmixed indices — 200, 020, 220, 400 — because face "
            "centring extinguishes the mixed ones, and that alone identifies the lattice before "
            "any measurement. The pattern is four-fold, so the two shortest independent spots are "
            "at 90° and their lengths are equal; the next ring out is at 45° and longer by √2. "
            "Index "
            "the innermost four spots first: the ratio 1 : 1 : √2 with a 45° angle is the "
            "signature of a cubic ⟨001⟩ zone and does not depend on the camera constant at all. "
            "The camera constant only enters when you want the lattice parameter, which is why a "
            "wrong one gives a self-consistent pattern of the wrong material."
        ),
        phase_key="al_fcc",
        zone_axis=(0, 0, 1),
        in_plane_rotation_deg=17.0,
        centre_px=(512.0, 512.0),
        rng_seed=1001,
        targets=(
            SuggestedTarget(
                (0, 1, 1),
                "45° away and six-fold degenerate — the standard next stop, and close enough "
                "that a wide holder reaches it in one move.",
            ),
            SuggestedTarget(
                (1, 1, 1),
                "54.74° away: the six-fold pattern, and the axis where a stacking fault or a twin "
                "on {111} shows its streaking most clearly.",
            ),
            SuggestedTarget(
                (1, 1, 2),
                "35.26° away, and the cheapest useful move on a narrow holder — a second axis "
                "this close is usually enough to fix the rotation about the beam.",
            ),
        ),
    ),
    GalleryEntry(
        identifier="bcc_fe_110",
        title="Ferrite (bcc Fe) down [110]",
        summary="A centred rectangle, two-fold; the pattern most often mistaken for a square.",
        teaches=(
            "Body centring extinguishes every reflection with h + k + l odd, and the zone law "
            "h + k = 0 admits only half of what survives that, so the two shortest vectors here "
            "are 11̄0 — a member of {110} — and 002. They are perpendicular but *not* equal in "
            "length: |g| is √2/a and 2/a respectively, so their ratio is exactly √2 for any bcc "
            "metal, whatever a is. The rectangle's aspect "
            "ratio is a lattice-independent check on your calibration: measure it, and if it is "
            "not √2 the beam centre is misplaced or the pattern is not what you think. This is "
            "the classic trap — a bcc ⟨110⟩ read as a cubic ⟨001⟩ square, which then indexes to a "
            "plausible and entirely wrong lattice parameter."
        ),
        phase_key="fe_bcc",
        zone_axis=(1, 1, 0),
        in_plane_rotation_deg=-26.0,
        centre_px=(486.0, 534.0),
        rng_seed=1002,
        targets=(
            SuggestedTarget(
                (1, 1, 2),
                "30° to the nearest member of ⟨112⟩ — the cheapest useful hop from here, and the "
                "standard axis for imaging ½⟨111⟩ screw dislocations edge-on.",
            ),
            SuggestedTarget(
                (1, 1, 1),
                "35.26° away: the six-fold pattern, and the axis that separates the bcc variants "
                "of a Kurdjumov-Sachs relationship.",
            ),
            SuggestedTarget(
                (0, 0, 1),
                "45° to the nearest ⟨100⟩ member — the four-fold pattern, and ⟨110⟩ plus ⟨100⟩ "
                "together fix the orientation completely.",
            ),
        ),
    ),
    GalleryEntry(
        identifier="hcp_zr_2-1-10",
        title="Zirconium (hcp) down [2̄110]",
        summary="A rectangle whose aspect ratio is the axial ratio c/a — measured, not assumed.",
        teaches=(
            "This is the prism zone, written [2̄110] in Miller-Bravais and [100] in three indices. "
            "Two reflections dominate it: 0002 along the c* direction and 011̄0 perpendicular to "
            "it. Their lengths are 2/c and 2/(√3·a), so the rectangle's aspect ratio is √3·a/c — "
            "1.088 for zirconium, 1.091 for titanium, 1.067 for magnesium — fixed by the axial "
            "ratio alone. This one pattern therefore measures c/a and separates the hcp metals "
            "from each other with no camera-constant calibration at all, because a ratio of two "
            "lengths on the same plate does not care what the calibration is. Note also that "
            "0001 is absent "
            "while 0002 is present: the hcp basis extinguishes odd l on the 000l row, and a "
            "0001 spot appearing in a real plate is double diffraction, not a lattice reflection."
        ),
        phase_key="zr_hcp",
        # [2-1-10] in Miller-Bravais is [100] in the three-index basis PyTex
        # computes in: u = 2U + V = 3, v = 2V + U = 0, w = W = 0, reduced.
        zone_axis=(1, 0, 0),
        in_plane_rotation_deg=8.0,
        centre_px=(528.0, 498.0),
        rng_seed=1003,
        targets=(
            SuggestedTarget(
                (1, -1, 0),
                "[11̄00] in Miller-Bravais: 30° away, the second prism zone, and the pair "
                "separates prismatic from pyramidal slip traces.",
            ),
            SuggestedTarget(
                (1, 0, 1),
                "[21̄1̄3], 57.87° away — the pyramidal zone where ⟨c+a⟩ dislocations are imaged, "
                "and beyond a single move on a conventional holder.",
            ),
            SuggestedTarget(
                (0, 0, 1),
                "The basal axis [0001], exactly 90° away whatever c/a is — unreachable in one "
                "move on any conventional holder, and the standard demonstration of why an "
                "intermediate axis is needed.",
            ),
        ),
    ),
)


def gallery_options() -> tuple[tuple[str, str, str], ...]:
    """The gallery as ``(value, label, help)`` triples for a choice control."""

    return tuple((entry.identifier, entry.title, entry.summary) for entry in GALLERY)


def gallery_entry(identifier: str) -> GalleryEntry:
    """One gallery entry by identifier.

    Raises
    ------
    InvalidInputError
        For an unknown identifier, naming the ones that exist — the frontend
        sends this value, so a stale bookmark must produce a usable message
        rather than a stack trace.
    """

    for entry in GALLERY:
        if entry.identifier == identifier:
            return entry
    raise InvalidInputError(
        f"There is no gallery pattern called {identifier!r}.",
        field="pattern",
        hint="Available: " + ", ".join(item.identifier for item in GALLERY) + ".",
    )
