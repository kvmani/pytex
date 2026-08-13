"""Diffraction patterns, including the composite pattern of a two-phase crystal.

What it does
    Simulates the selected-area pattern down a chosen zone axis, and — the
    reason this panel exists — the *composite* pattern of a parent phase
    together with the product variants an orientation relationship generates
    from it, all on one detector.

When to use it
    Before a session, to know what a pattern down a given zone axis should look
    like. After one, to decide which of the twenty-four martensite variants a
    set of extra spots belongs to. The composite pattern is the thing a
    two-phase micrograph actually shows, and reading it by hand means keeping
    twenty-four rotated reciprocal lattices in your head at once.

What every spot carries
    Each spot travels with its full row — indices, d, |g|, relative intensity,
    which phase, which variant, whether it is a double-diffraction spot — taken
    from :func:`pytex.diffraction.export.composite_reflection_table`. That is
    the same table the CSV export writes, so what a hover says and what the file
    contains cannot disagree.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pytex.app.errors import InvalidInputError
from pytex.app.phases import phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
    NumberParameter,
    TextParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import (
    _RELATIONSHIP_CONSTRUCTORS,
    _RELATIONSHIPS,
    direction_label,
    phase_parameter,
    relationship_name,
)

__all__: tuple[str, ...] = ()

_CITATION_WILLIAMS = (
    "Williams & Carter, Transmission Electron Microscopy, 2nd ed., Part 2 (diffraction)."
)
_CITATION_MORITO = "Morito et al., Acta Mater. 51 (2003) 1789 (variant numbering)."

#: Columns of the spot table, shared by the on-screen hover card and every
#: export. Declared once so the two cannot describe a spot differently.
_SPOT_COLUMNS: tuple[Column, ...] = (
    Column("origin", "Source", help_text="Parent, or the product variant the spot belongs to."),
    Column("phase", "Phase"),
    Column("hkl_label", "Reflection"),
    Column(
        "zone_axis_indexed",
        "Zone axis",
        help_text=(
            "For a variant, the nearest low-index zone axis of the product crystal, with how far "
            "the exact axis lies from it. The exact axis is irrational: a rational parent axis "
            "maps to an irrational child axis under a real orientation relationship, and rounding "
            "it silently would be the more misleading choice."
        ),
    ),
    Column("d_angstrom", "d", units="Å", numeric=True, digits=5),
    Column("g_inv_angstrom", "|g|", units="Å⁻¹", numeric=True, digits=5),
    Column("detector_x_mm", "x", units="mm", numeric=True, digits=4),
    Column("detector_y_mm", "y", units="mm", numeric=True, digits=4),
    Column(
        "relative_intensity",
        "Intensity",
        numeric=True,
        digits=4,
        help_text="Kinematic, relative to the strongest spot in the pattern.",
    ),
    Column(
        "excitation_error_inv_angstrom",
        "s",
        units="Å⁻¹",
        numeric=True,
        digits=5,
        help_text="Excitation error: how far the reflection sits from the Ewald sphere.",
    ),
    Column("double_diffraction", "Double diffraction"),
)


def _relationship(name: str, parent: Any, child: Any) -> Any:
    from pytex.core.transformation import OrientationRelationship

    constructor = getattr(OrientationRelationship, _RELATIONSHIP_CONSTRUCTORS[name])
    try:
        return constructor(parent_phase=parent, child_phase=child)
    except (ValueError, TypeError) as error:
        raise InvalidInputError(
            f"The {relationship_name(name)} relationship does not apply to these phases: {error}",
            field="relationship",
            hint=(
                "The fcc-to-bcc relationships need a cubic parent and a cubic child; Burgers "
                "needs a cubic parent and a hexagonal child."
            ),
        ) from error


def _variant_selection(text: str | None, available: int) -> tuple[int, ...] | None:
    """Parse a variant selection such as ``"1 2 5"`` or ``"1-6"``.

    Returning ``None`` for an empty selection means "all variants", which is the
    honest default: the point of a composite pattern is that a real
    micrograph contains whichever variants happen to be present, and hiding some
    by default would teach the wrong lesson.
    """

    if not text or not str(text).strip():
        return None
    chosen: list[int] = []
    for token in str(text).replace(",", " ").split():
        if "-" in token[1:]:
            start_text, _, end_text = token.partition("-")
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise InvalidInputError(
                    f"{token!r} is not a variant range.",
                    field="variants",
                    hint="Use numbers and ranges, for example '1 2 5' or '1-6'.",
                ) from error
            chosen.extend(range(start, end + 1))
        else:
            try:
                chosen.append(int(token))
            except ValueError as error:
                raise InvalidInputError(
                    f"{token!r} is not a variant number.",
                    field="variants",
                    hint="Use numbers and ranges, for example '1 2 5' or '1-6'.",
                ) from error
    out_of_range = sorted({value for value in chosen if not 1 <= value <= available})
    if out_of_range:
        raise InvalidInputError(
            f"This relationship has {available} variants; {out_of_range} are outside that range.",
            field="variants",
            hint=f"Choose from 1 to {available}, or leave the field empty for all of them.",
        )
    return tuple(sorted(set(chosen)))


@REGISTRY.operation(
    "diffraction.composite_saed",
    title="Composite SAED pattern",
    summary="Parent and product-variant spots on one detector, as a real pattern shows them.",
    help_text=(
        "Simulates the selected-area pattern of a parent crystal down a chosen zone axis, "
        "together with the patterns of the product variants an orientation relationship "
        "generates from it — all on the shared detector, in the correct relative orientation.\n\n"
        "This is what a two-phase micrograph actually contains, and it is the pattern that is "
        "hard to read by hand: twenty-four Kurdjumov-Sachs variants put twenty-four rotated "
        "reciprocal lattices on the same plate, and deciding which variant a given extra spot "
        "belongs to means testing all of them. Hovering a spot answers that directly.\n\n"
        "Start with all variants to see the full complexity, then restrict the selection to one "
        "or two and watch the pattern resolve — that comparison is the fastest way to learn to "
        "read a martensite pattern.\n\n"
        "**What the intensities do and do not say.** They are kinematic and relative to the "
        "strongest spot, and the scattering model uses the atomic number with no angular "
        "dependence. Two consequences are worth knowing before reading anything into them. For a "
        "**monatomic** phase — iron, nickel, zirconium — every reflection the centring allows "
        "comes out at the same intensity, so the pattern shows which spots exist rather than how "
        "bright they are. Intensity differences appear only where the structure factor genuinely "
        "differs between reflections, as in rock salt, where the 111 is weak and the 200 strong. "
        "And in every case dynamical scattering makes measured intensities differ from kinematic "
        "ones, sometimes by a lot."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text="The phase the zone axis is stated in — austenite, or beta for Burgers.",
            builtin="austenite_fcc",
        ),
        ChoiceParameter(
            name="relationship",
            label="Orientation relationship",
            help_text="Which relationship generates the product variants.",
            options=_RELATIONSHIPS,
            default="kurdjumov_sachs",
        ),
        IndicesParameter(
            name="zone_axis",
            label="Parent zone axis [uvw]",
            help_text=(
                "The beam direction, in the parent crystal. Low-index axes give the patterns "
                "worth working from: [001], [011], [111]."
            ),
            default=(0, 0, 1),
        ),
        TextParameter(
            name="variants",
            label="Variants to include",
            help_text=(
                "Numbers or ranges, for example '1 2 5' or '1-6'. Leave empty for every variant, "
                "which is what a real two-phase area may contain."
            ),
            required=False,
            placeholder="all",
        ),
        BooleanParameter(
            name="include_parent",
            label="Include the parent pattern",
            help_text="Show the parent's own reflections alongside the product spots.",
            default=True,
        ),
        NumberParameter(
            name="max_g_inv_angstrom",
            label="Largest |g| to include",
            help_text=(
                "Reflections beyond this reciprocal-space radius are dropped. Raising it fills "
                "the pattern with high-order spots that a real aperture would exclude anyway."
            ),
            units="Å⁻¹",
            default=1.6,
            minimum=0.2,
            maximum=6.0,
            advanced=True,
        ),
        NumberParameter(
            name="min_intensity",
            label="Weakest spot to draw",
            help_text=(
                "Relative to the strongest spot. Zero draws everything, including reflections "
                "far too weak to see on a plate. For a monatomic phase this control does nothing, "
                "because the kinematic model gives every allowed reflection the same intensity."
            ),
            default=0.005,
            minimum=0.0,
            maximum=1.0,
            advanced=True,
        ),
        NumberParameter(
            name="beam_energy_kev",
            label="Accelerating voltage",
            help_text=(
                "Sets the electron wavelength and so the curvature of the Ewald sphere. 200 kV "
                "is the common instrument; at 300 kV the sphere is flatter and more reflections "
                "are excited at once."
            ),
            units="kV",
            default=200.0,
            minimum=20.0,
            maximum=1000.0,
            group="Instrument",
        ),
        NumberParameter(
            name="camera_constant_mm_angstrom",
            label="Camera constant",
            help_text=(
                "The instrument constant relating a spot's radius on the plate to 1/d: "
                "r = (camera constant) / d. Set it to your microscope's calibrated value to get "
                "millimetres you can measure against a real plate."
            ),
            units="mm·Å",
            default=180.0,
            minimum=1.0,
            maximum=5000.0,
            group="Instrument",
        ),
        IntegerParameter(
            name="max_index",
            label="Index limit",
            help_text="Largest |h|, |k| or |l| the simulation enumerates.",
            default=4,
            minimum=1,
            maximum=8,
            advanced=True,
        ),
        phase_parameter(
            name="child_phase",
            label="Product phase",
            help_text="The product phase — ferrite or martensite, or alpha for Burgers.",
            builtin="fe_bcc",
        ),
    ),
    returns="One row per spot; detector geometry and per-source grouping under `data`.",
    panel="diffraction",
    citations=(_CITATION_WILLIAMS, _CITATION_MORITO),
    tags=(
        "SAED",
        "diffraction",
        "composite",
        "variant",
        "martensite",
        "zone axis",
        "electron",
    ),
)
def _composite_saed(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.lattice import CrystalDirection
    from pytex.diffraction.composite import simulate_composite_saed
    from pytex.diffraction.export import composite_reflection_table
    from pytex.diffraction.kinematic import KinematicSimulationConfig

    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    relationship = _relationship(str(request["relationship"]), parent_phase, child_phase)
    available = len(relationship.generate_variants())
    variants = _variant_selection(request.get("variants"), available)
    zone_indices = tuple(request["zone_axis"])

    config = KinematicSimulationConfig(
        max_index=int(request["max_index"]),
        g_max_inv_angstrom=float(request["max_g_inv_angstrom"]),
        beam_energy_kev=float(request["beam_energy_kev"]),
        camera_constant_mm_angstrom=float(request["camera_constant_mm_angstrom"]),
    )
    pattern = simulate_composite_saed(
        relationship,
        CrystalDirection(coordinates=np.asarray(zone_indices, dtype=float), phase=parent_phase),
        variant_indices=list(variants) if variants else None,
        include_parent=bool(request["include_parent"]),
        config=config,
    )
    table = composite_reflection_table(pattern, intensity_threshold=float(request["min_intensity"]))
    rows = tuple(dict(record) for record in table.to_records())
    # The library's row carries the exact (irrational) child zone axis. A hover
    # card wants the nearest indexable axis and the deviation from it, which is
    # what a microscopist would write down, so both are attached here.
    nearest: dict[int, str] = {}
    for variant_pattern in pattern.variant_patterns:
        rationalized = variant_pattern.nearest_zone_axis
        indices = " ".join(str(int(value)) for value in rationalized.indices)
        nearest[int(variant_pattern.variant.variant_index)] = (
            f"[{indices}] ({rationalized.deviation_deg:.2f}° off)"
        )
    for row in rows:
        if row["source"] == "parent":
            row["origin"] = "Parent"
            row["zone_axis_indexed"] = str(row["zone_axis"])
            continue
        variant_index = int(row["variant"])
        row["origin"] = f"Variant {variant_index}"
        row["zone_axis_indexed"] = nearest.get(variant_index, str(row["zone_axis"]))
    if not rows:
        raise InvalidInputError(
            "No reflection survives these limits, so there is nothing to draw.",
            field="min_intensity",
            hint="Lower the weakest-spot threshold, or raise the largest |g|.",
        )

    radius = max(abs(float(row["detector_r_mm"])) for row in rows) or 1.0
    # Group by (source, variant), not by source: the whole point of the panel is
    # to distinguish variant 3 from variant 17, and a single "child" bucket makes
    # exactly the distinction the user came for impossible.
    sources: dict[tuple[str, Any], dict[str, Any]] = {}
    for row in rows:
        key = (str(row["source"]), row["variant"])
        sources.setdefault(
            key,
            {
                "source": key[0],
                "variant": key[1],
                "label": "Parent" if key[0] == "parent" else f"Variant {key[1]}",
                "spots": [],
            },
        )
        sources[key]["spots"].append(row)

    notes = [
        "Detector coordinates use the camera constant reported under `data`. A spot's radius is "
        "set by the *in-plane* component of g, so a reflection with a non-zero excitation error "
        "sits marginally inside r = (camera constant) / d.",
    ]
    # Say it here, where the numbers are, rather than only in the help: a
    # monatomic phase produces a pattern in which every allowed spot has
    # intensity 1.0, and a user who reads that as a prediction will be misled.
    monatomic = [
        spec.name
        for spec in (parent_spec, child_spec)
        if len({site.species for site in spec.sites}) <= 1
    ]
    if monatomic:
        notes.append(
            f"{' and '.join(monatomic)} "
            + ("is" if len(monatomic) == 1 else "are")
            + " monatomic, so the kinematic model gives every allowed reflection the same "
            "intensity. Read the pattern for which spots are present, not for how bright they are."
        )

    zone_text = direction_label(zone_indices, spec=parent_spec)
    variant_text = (
        f"variants {', '.join(str(value) for value in variants)}"
        if variants
        else f"all {available} variants"
    )
    result = AppResult(
        title=f"Composite SAED down {zone_text} of {parent_spec.name}",
        summary=(
            f"{len(rows)} spots down {zone_text}: {parent_spec.name} with {variant_text} of "
            f"{child_spec.name} under the "
            f"{relationship_name(str(request['relationship']))} relationship, on one detector. "
            "Positions are exact; intensities are kinematic and should be read with the note "
            "below."
        ),
        table=ResultTable(
            columns=_SPOT_COLUMNS,
            rows=rows,
            caption=f"Reflections in the composite pattern down {zone_text}.",
        ),
        data={
            "spots": rows,
            "sources": list(sources.values()),
            "detector_radius_mm": radius,
            "camera_constant_mm_angstrom": table.camera_constant_mm_angstrom,
            "wavelength_angstrom": table.wavelength_angstrom,
            "beam_energy_kev": table.beam_energy_kev,
            "zone_axis_label": zone_text,
            "variant_count": available,
            "columns": [column.to_json() for column in _SPOT_COLUMNS],
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": request["relationship"],
            "zone_axis": list(zone_indices),
            "variants": request.get("variants"),
            "include_parent": bool(request["include_parent"]),
            "max_g_inv_angstrom": float(request["max_g_inv_angstrom"]),
            "min_intensity": float(request["min_intensity"]),
            "max_index": int(request["max_index"]),
            "beam_energy_kev": float(request["beam_energy_kev"]),
            "camera_constant_mm_angstrom": float(request["camera_constant_mm_angstrom"]),
        },
        notes=tuple(notes),
        citations=(_CITATION_WILLIAMS, _CITATION_MORITO),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="diffraction.example.ks_001",
            title="Why a martensite pattern looks crowded",
            panel="diffraction",
            summary="Austenite [001] with all 24 Kurdjumov-Sachs variants.",
            teaches=(
                "Every extra spot belongs to one of 24 rotated reciprocal lattices, and hovering "
                "says which. This is the pattern that makes variant identification by hand "
                "impractical — and the reason the next example is worth running immediately "
                "after."
            ),
            operation="diffraction.composite_saed",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "zone_axis": [0, 0, 1],
            },
        ),
        ExampleScenario(
            id="diffraction.example.ks_single",
            title="The same pattern with one variant",
            panel="diffraction",
            summary="Austenite [001] with a single Kurdjumov-Sachs variant.",
            teaches=(
                "With one variant the pattern is simple enough to index by eye, and the "
                "relationship between the parent square array and the product spots becomes "
                "visible. Switch back to all 24 and the same structure is still there, 24 times "
                "over."
            ),
            operation="diffraction.composite_saed",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "zone_axis": [0, 0, 1],
                "variants": "1",
            },
        ),
        ExampleScenario(
            id="diffraction.example.nw",
            title="Nishiyama-Wassermann, for comparison",
            panel="diffraction",
            summary="The same zone axis under the 12-variant relationship.",
            teaches=(
                "N-W puts 12 variants where K-S puts 24, and the missing ones are the "
                "5.26°-rotated partners. Comparing the two patterns down the same axis is how "
                "the two relationships are told apart experimentally."
            ),
            operation="diffraction.composite_saed",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "nishiyama_wassermann",
                "zone_axis": [0, 0, 1],
            },
        ),
        ExampleScenario(
            id="diffraction.example.burgers",
            title="Burgers: bcc beta to hexagonal alpha",
            panel="diffraction",
            summary="A bcc [111] zone axis with the 12 Burgers variants of hcp zirconium.",
            teaches=(
                "The hexagonal product puts six-fold arrays inside the parent's three-fold one. "
                "This is the pattern seen in quenched titanium and zirconium, and the same "
                "machinery that handled steel handles it with no change of method."
            ),
            operation="diffraction.composite_saed",
            request={
                "phase": {"builtin": "fe_bcc"},
                "child_phase": {"builtin": "zr_hcp"},
                "relationship": "burgers",
                "zone_axis": [1, 1, 1],
            },
        ),
    )
)
