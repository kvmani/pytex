"""Transformation variants: where they point, and how they differ from each other.

What it does
    Takes an orientation relationship and shows the *set* of child orientations
    one parent grain produces. Two views, because the two questions a variant
    set is asked are different questions.

    The **pole figure** answers "where do the child poles land?". Every variant
    contributes its whole symmetry family of a chosen child plane, projected
    into the parent frame, so the figure is directly comparable with a measured
    pole figure from the same parent grain. Variants can be coloured by packet,
    which is what makes the structure visible: under Kurdjumov-Sachs the 24
    variants fall into 4 packets of 6, one per member of the parent {111}
    family, and the packet is what a lath martensite micrograph shows as a
    block.

    The **intervariant misorientation spectrum** answers "what boundaries do
    these variants make with each other?". Two child grains descended from one
    parent can only meet at certain misorientations, and that discrete spectrum
    is what an EBSD misorientation histogram is compared against when deciding
    whether a region shares a parent.

When to use it
    Before interpreting a variant-selection argument, to know what no selection
    looks like. After indexing an EBSD map, to compare the measured boundary
    spectrum against the admissible one. When teaching why one parent grain
    gives 24 child orientations and only 4 apparent plate directions.

What every pole carries
    Its variant, its packet, its indices in the child crystal, its position on
    the projection and its polar and azimuthal angles — the same row the CSV
    export writes, so a hover and a file cannot disagree.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from pytex.app.errors import DependencyMissingError, InvalidInputError
from pytex.app.phases import PhaseSpec, phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesParameter,
    IntegerParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import (
    _RELATIONSHIP_CONSTRUCTORS,
    _RELATIONSHIPS,
    direction_label,
    phase_parameter,
    plane_label,
    relationship_name,
)

__all__: tuple[str, ...] = ()

_CITATION_MORITO = (
    "Morito, Tanaka, Konishi, Furuhara & Maki, Acta Mater. 51 (2003) 1789 "
    "(packet and block structure of lath martensite)."
)
_CITATION_BUNGE = "Bunge, Texture Analysis in Materials Science (1982), chapter 2."
_CITATION_RANDLE = (
    "Randle & Engler, Introduction to Texture Analysis, 2nd ed., chapter 2 (projections)."
)

#: Columns of the pole table, shared by the hover card and every export.
_POLE_COLUMNS: tuple[Column, ...] = (
    Column("variant", "Variant", numeric=True, help_text="Index in `generate_variants()` order."),
    Column(
        "packet",
        "Packet",
        numeric=True,
        help_text=(
            "Variants sharing the parent plane they carry into exact parallelism. Under "
            "Kurdjumov-Sachs and the {111} family this is the packet of lath martensite."
        ),
    ),
    Column("pole", "Child plane", help_text="The family member of the plotted child plane."),
    Column(
        "x",
        "x",
        numeric=True,
        digits=5,
        help_text="Projected coordinate in the unit disc, in the parent frame.",
    ),
    Column("y", "y", numeric=True, digits=5),
    Column(
        "polar_deg",
        "Polar angle",
        units="°",
        numeric=True,
        digits=3,
        help_text="Angle from the projection pole (the parent z axis). 0° is the centre.",
    ),
    Column(
        "azimuth_deg",
        "Azimuth",
        units="°",
        numeric=True,
        digits=3,
        help_text="Angle around the projection, measured from the parent x axis.",
    ),
)

#: Columns of the intervariant misorientation table.
_PAIR_COLUMNS: tuple[Column, ...] = (
    Column("variant_a", "Variant A", numeric=True),
    Column("variant_b", "Variant B", numeric=True),
    Column(
        "angle_deg",
        "Disorientation",
        units="°",
        numeric=True,
        digits=4,
        help_text="Minimal angle over the child point group, which is what EBSD reports.",
    ),
    Column(
        "axis",
        "Axis",
        help_text=(
            "The nearest low-index direction of the child crystal to the rotation axis of the "
            "minimal representative. These axes are not all rational, so the deviation column "
            "beside it says how far the label is from the exact axis."
        ),
    ),
    Column(
        "axis_deviation_deg",
        "Axis deviation",
        units="°",
        numeric=True,
        digits=3,
        help_text="How far the exact axis lies from the low-index label. 0 means it is exact.",
    ),
    Column(
        "axis_x",
        "Axis x",
        numeric=True,
        digits=5,
        help_text="The exact axis, in Cartesian child-crystal components, for re-plotting.",
    ),
    Column("axis_y", "Axis y", numeric=True, digits=5),
    Column("axis_z", "Axis z", numeric=True, digits=5),
    Column(
        "same_packet",
        "Same packet",
        help_text="Whether the two variants share a parent habit plane.",
    ),
)


#: Packet colours, as hex, matching `PACKET_COLORS` in `js/panels/variants.js`.
#:
#: The two lists are the same colours in two languages, which is a duplication
#: worth naming: the browser cannot read a Python constant and matplotlib cannot
#: read a CSS one. `test_app_variants.py` compares them, so a change to either
#: without the other is a test failure rather than a figure that quietly stops
#: matching the screen it was published from.
_PACKET_COLORS: tuple[str, ...] = (
    "#2c7fdd",
    "#e87917",
    "#259d71",
    "#ce3b98",
    "#8c56d2",
    "#b8a41e",
)

#: The parent's own poles: achromatic, because they are the reference.
_PARENT_COLOR = "#333333"


def _unit_from_angles(polar_deg: float, azimuth_deg: float) -> list[float]:
    """Rebuild a unit vector from the polar and azimuthal angles in a pole row.

    The figure is drawn from the table rather than from the vectors that
    produced it, so that the published figure and the exported numbers are the
    same numbers. The round trip is exact to floating point: the table stores
    the angles of the folded vector, and the projection depends on nothing else.
    """

    polar = math.radians(polar_deg)
    azimuth = math.radians(azimuth_deg)
    return [
        math.sin(polar) * math.cos(azimuth),
        math.sin(polar) * math.sin(azimuth),
        math.cos(polar),
    ]


#: Largest |u|, |v| or |w| considered when naming a disorientation axis.
#:
#: Eight rather than six because the martensite intervariant axes are quoted in
#: the literature up to ⟨5 5 7⟩ (Morito et al., Table 2), and a limit of six
#: cannot reach them: it labelled the 14.88° pair as ⟨2 5 0⟩ nearly four degrees
#: away, which is a label worse than none. Beyond eight an axis is irrational
#: for practical purposes, and the deviation column says so rather than the
#: label growing longer.
_AXIS_MAX_INDEX = 8


def _canonical_sign(triple: tuple[int, int, int]) -> tuple[int, int, int]:
    """Return the representative of an axis with its first nonzero index positive.

    A rotation axis is a line, not a ray: [1 -1 1] and [-1 1 -1] name the same
    axis, and the sign that comes out of a symmetry reduction is arbitrary.
    Without this the same physical axis appears under two labels in one table,
    and the reader is left to work out that they are the same thing.
    """

    for value in triple:
        if value != 0:
            return triple if value > 0 else (-triple[0], -triple[1], -triple[2])
    return triple


def _nearest_axis(axis_cartesian: np.ndarray, phase: Any, spec: PhaseSpec) -> tuple[str, float]:
    """Name the nearest low-index lattice direction to a Cartesian axis.

    A disorientation axis comes out of the symmetry reduction as three Cartesian
    components, and three decimals is not how anyone quotes one: the interesting
    fact about the 60° Kurdjumov-Sachs pairs is that they are about ⟨111⟩, and
    ``+0.577 +0.577 +0.577`` is that fact in a form that has to be decoded.

    The residual is returned with the label rather than hidden, because these
    axes are *not* all rational — an orientation relationship carrying a
    rational parent axis produces irrational child axes in general — and a
    label with no error attached would silently claim otherwise.

    Returns the label and the angular residual in degrees. An axis further than
    a few degrees from every candidate is reported as irrational by the caller.
    """

    from itertools import product
    from math import gcd

    # Columns of `matrix` are the Cartesian images of a, b, c, so `matrix @ uvw`
    # is the Cartesian vector of the direction [uvw].
    basis = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    target = np.asarray(axis_cartesian, dtype=float)
    norm = float(np.linalg.norm(target))
    if norm < 1e-12:
        return "—", 180.0
    target = target / norm

    # Candidate triples, primitive and sign-canonical: an axis and its reverse
    # describe the same rotation axis once the angle's sign is fixed, and a
    # non-primitive triple names the same direction as its reduced form.
    candidates: list[tuple[int, int, int]] = []
    span = range(-_AXIS_MAX_INDEX, _AXIS_MAX_INDEX + 1)
    for triple in product(span, span, span):
        if triple == (0, 0, 0):
            continue
        if gcd(gcd(abs(triple[0]), abs(triple[1])), abs(triple[2])) != 1:
            continue
        candidates.append(triple)

    images = np.asarray(candidates, dtype=float) @ basis.T
    units = images / np.linalg.norm(images, axis=1)[:, None]
    cosines = units @ target
    best = int(np.argmax(np.abs(cosines)))
    residual = float(np.degrees(np.arccos(min(abs(float(cosines[best])), 1.0))))
    return direction_label(_canonical_sign(candidates[best]), spec=spec), residual


def _relationship(name: str, parent: Any, child: Any) -> Any:
    """Build a named relationship, or explain why these phases cannot carry it."""

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


def _packet_labels(relationship: Any, parent_spec: PhaseSpec, indices: tuple[int, ...]) -> Any:
    """Group the variants by the parent plane each carries into parallelism.

    Raises an :class:`InvalidInputError` naming the packet-plane field rather
    than letting a library ``ValueError`` reach the user as an internal message:
    choosing a plane that is not the relationship's defining family is an easy
    and reasonable thing for a user to try.
    """

    from pytex.core.lattice import CrystalPlane, MillerIndex
    from pytex.core.transformation import variant_close_packed_groups

    parent_phase = relationship.parent_phase
    plane = CrystalPlane(
        MillerIndex(np.asarray(indices, dtype=int), phase=parent_phase),
        phase=parent_phase,
    )
    try:
        return variant_close_packed_groups(relationship, plane)
    except (ValueError, TypeError) as error:
        raise InvalidInputError(
            f"Variants cannot be grouped by {plane_label(indices, spec=parent_spec)} of "
            f"{parent_spec.name}: {error}",
            field="packet_plane",
            hint=(
                "Use the family the relationship is defined on: {111} for the fcc-to-bcc "
                "relationships, {110} for Burgers."
            ),
        ) from error


def _project(normals: np.ndarray, method: str) -> np.ndarray:
    """Project unit normals onto the unit disc, whichever projection is asked for.

    :func:`pytex.texture.projections.project_directions` returns each projection
    in its own natural radius: the stereographic net reaches 1 at the equator,
    the equal-area net reaches √2. Both are correct and the difference is
    invisible until something draws them in the same circle — at which point the
    equal-area figure spills over the rim by 41%.

    Dividing by the equatorial radius makes the rim mean "90° from the
    projection axis" in both, which is what a reader assumes and what makes the
    exported x and y comparable between the two. It changes no angle and no
    relative spacing: it is a choice of unit, fixed here so that the panel, the
    CSV and any figure built from the same rows cannot disagree about it.
    """

    from pytex.texture.projections import project_directions

    projected = np.asarray(project_directions(normals, method=method), dtype=float)
    equatorial_radius = math.sqrt(2.0) if method == "equal_area" else 1.0
    return projected / equatorial_radius


def _child_family(child_phase: Any, indices: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    """Return the symmetry family of a child plane, as indices and unit normals.

    Both are needed and they must correspond row for row: the normals are what
    is projected, the indices are what the hover card and the CSV name.
    """

    from pytex.core.miller import MillerPlaneSet

    single = MillerPlaneSet.from_hkl(np.asarray([indices], dtype=int), phase=child_phase)
    members, mask = single.symmetry_equivalent_indices()
    valid = np.asarray(members[0])[np.asarray(mask[0], dtype=bool)]
    if valid.size == 0:
        raise InvalidInputError(
            "That plane has no symmetry family; check the indices.", field="pole"
        )
    family = MillerPlaneSet.from_hkl(np.asarray(valid, dtype=int), phase=child_phase)
    normals = np.asarray(family.normals_cartesian(), dtype=float)
    return np.asarray(valid, dtype=int), normals


@REGISTRY.operation(
    "variants.pole_figure",
    title="Variant pole figure",
    summary="Where the child poles of every transformation variant land, in the parent frame.",
    help_text=(
        "Plots the chosen child plane family of every transformation variant, projected into "
        "the parent crystal frame. This is the figure that makes a variant set legible: one "
        "parent grain under Kurdjumov-Sachs produces 24 child orientations, and asking where "
        "each one puts its {100} is the question a measured pole figure of that grain answers "
        "too, so the two are directly comparable.\n\n"
        "**Colour by packet to see the structure.** The 24 Kurdjumov-Sachs variants fall into 4 "
        "groups of 6, one per member of the parent {111} family, because each variant carries "
        "exactly one {111} into exact parallelism with a child {110}. That group is the *packet* "
        "of lath martensite, and it is what a micrograph shows as a block of parallel laths. "
        "Burgers gives 6 groups of 2 on the parent {110} family instead.\n\n"
        "**Reading the projection.** The centre of the disc is the parent z axis and the rim is "
        "the equator; poles in the lower hemisphere are folded up, which is standard for a plane "
        "normal because a plane and its opposite normal are the same plane. Equal-area is the "
        "projection to use when density matters, because it does not concentrate area toward the "
        "rim; stereographic is the one that preserves angles, and is what a Wulff net measures.\n\n"
        "**What this does not claim.** Every variant is drawn with equal weight. A real "
        "microstructure shows variant *selection* — some variants far more often than others — "
        "and that is a property of the transformation conditions, not of the crystallography. "
        "This figure is the no-selection baseline a selection argument is measured against."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text="The phase that transforms, and the frame the figure is drawn in.",
            builtin="austenite_fcc",
        ),
        phase_parameter(
            name="child_phase",
            label="Child phase",
            help_text="The product phase — ferrite or martensite, or alpha for Burgers.",
            builtin="fe_bcc",
        ),
        ChoiceParameter(
            name="relationship",
            label="Orientation relationship",
            help_text="Which relationship generates the variants.",
            options=_RELATIONSHIPS,
            default="kurdjumov_sachs",
        ),
        IndicesParameter(
            name="pole",
            label="Child plane to plot",
            help_text=(
                "The child plane whose symmetry family is projected. (100) is the usual choice "
                "for a cubic product: three poles per variant, so the figure stays readable."
            ),
            default=(1, 0, 0),
        ),
        IndicesParameter(
            name="packet_plane",
            label="Parent plane defining packets",
            help_text=(
                "The parent family whose members the variants are grouped by. Use the family "
                "the relationship is built on: (111) for the fcc-to-bcc relationships, (110) "
                "for Burgers."
            ),
            default=(1, 1, 1),
        ),
        ChoiceParameter(
            name="projection",
            label="Projection",
            help_text="How the sphere is flattened onto the disc.",
            options=(
                (
                    "stereographic",
                    "Stereographic",
                    "Preserves angles; what a Wulff net measures.",
                ),
                (
                    "equal_area",
                    "Equal area (Schmidt)",
                    "Preserves area; use when pole density is the point.",
                ),
            ),
            default="stereographic",
        ),
        BooleanParameter(
            name="include_parent",
            label="Include the parent's own poles",
            help_text=(
                "Draws the same indices in the parent crystal, unrotated. It gives the figure a "
                "reference frame: the parent poles are where a pole would sit with no "
                "transformation at all."
            ),
            default=True,
        ),
    ),
    returns="One row per plotted pole; per-variant and per-packet summaries under `data`.",
    panel="variants",
    citations=(_CITATION_MORITO, _CITATION_RANDLE, _CITATION_BUNGE),
    tags=(
        "variant",
        "pole figure",
        "packet",
        "martensite",
        "orientation relationship",
        "OR",
        "stereographic",
        "projection",
    ),
)
def _variant_pole_figure(request: dict[str, Any]) -> dict[str, Any]:
    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    name = str(request["relationship"])
    relationship = _relationship(name, parent_phase, child_phase)
    variants = relationship.generate_variants()
    packets = _packet_labels(relationship, parent_spec, tuple(request["packet_plane"]))
    method = str(request["projection"])

    pole_indices, pole_normals = _child_family(child_phase, tuple(request["pole"]))

    # Every variant contributes the whole family at once. `parent_to_child`
    # carries a parent vector into the child frame, so its transpose is what
    # brings a child pole back into the parent frame the figure is drawn in.
    matrices = np.stack(
        [variant.parent_to_child_rotation.as_matrix() for variant in variants], axis=0
    )
    in_parent = np.einsum("vji,pj->vpi", matrices, pole_normals, optimize=True)

    flat = in_parent.reshape(-1, 3)
    projected = _project(flat, method)
    # `project_directions` folds to the upper hemisphere before projecting, and
    # the polar angle must be read from the folded vector, not the original, or
    # a pole at 175° would be reported as such while being drawn at 5°.
    folded = np.where(flat[:, 2:3] < 0.0, -flat, flat)
    polar_deg = np.degrees(np.arccos(np.clip(folded[:, 2], -1.0, 1.0)))
    azimuth_deg = np.degrees(np.arctan2(folded[:, 1], folded[:, 0])) % 360.0

    rows: list[dict[str, Any]] = []
    position = 0
    for variant_position, variant in enumerate(variants):
        for member in pole_indices:
            rows.append(
                {
                    "variant": int(variant.variant_index),
                    "packet": int(packets[variant_position]) + 1,
                    "pole": plane_label(tuple(int(value) for value in member), spec=child_spec),
                    "x": float(projected[position, 0]),
                    "y": float(projected[position, 1]),
                    "polar_deg": float(polar_deg[position]),
                    "azimuth_deg": float(azimuth_deg[position]),
                }
            )
            position += 1

    parent_rows: list[dict[str, Any]] = []
    if bool(request["include_parent"]):
        parent_indices, parent_normals = _child_family(parent_phase, tuple(request["pole"]))
        parent_projected = _project(parent_normals, method)
        parent_folded = np.where(parent_normals[:, 2:3] < 0.0, -parent_normals, parent_normals)
        for member, point, vector in zip(
            parent_indices, parent_projected, parent_folded, strict=True
        ):
            parent_rows.append(
                {
                    "variant": 0,
                    "packet": 0,
                    "pole": plane_label(tuple(int(value) for value in member), spec=parent_spec),
                    "x": float(point[0]),
                    "y": float(point[1]),
                    "polar_deg": float(np.degrees(np.arccos(np.clip(vector[2], -1.0, 1.0)))),
                    "azimuth_deg": float(np.degrees(np.arctan2(vector[1], vector[0])) % 360.0),
                }
            )

    packet_sizes: dict[int, int] = {}
    for label in packets:
        packet_sizes[int(label) + 1] = packet_sizes.get(int(label) + 1, 0) + 1
    # Uneven packets mean the chosen parent family is not the one this
    # relationship is built on — grouping Pitsch, which is defined on {100}, by
    # {111} gives 4/3/4/1 rather than four equal groups. The grouping is still
    # computed, because it is a real property of a real question; but a figure
    # whose colours mean nothing is worse than one that says so.
    uneven = len(set(packet_sizes.values())) > 1
    variant_rows = [
        {
            "variant": int(variant.variant_index),
            "packet": int(packets[position]) + 1,
            "angle_deg": float(variant.parent_to_child_rotation.angle_deg),
        }
        for position, variant in enumerate(variants)
    ]

    display = relationship_name(name)
    pole_text = plane_label(tuple(int(v) for v in request["pole"]), spec=child_spec)
    packet_text = plane_label(tuple(int(v) for v in request["packet_plane"]), spec=parent_spec)
    packet_count = len(packet_sizes)
    sizes = sorted(set(packet_sizes.values()))
    size_text = f"{sizes[0]}" if len(sizes) == 1 else "/".join(str(value) for value in sizes)
    result = AppResult(
        title=f"{display} variant pole figure: {pole_text} of {child_spec.name}",
        summary=(
            f"{len(variants)} variants of {child_spec.name} in {parent_spec.name}, each "
            f"contributing its {len(pole_indices)}-member {pole_text} family, so "
            f"{len(rows)} poles on a "
            f"{'stereographic' if method == 'stereographic' else 'equal-area'} projection. "
            f"Grouping by {packet_text} of the parent gives {packet_count} packets of "
            f"{size_text} variants."
        ),
        table=ResultTable(
            columns=_POLE_COLUMNS,
            rows=tuple(rows + parent_rows),
            caption=(
                f"{pole_text} poles of every {display} variant, projected into the "
                f"{parent_spec.name} frame."
            ),
        ),
        data={
            "poles": rows,
            "parent_poles": parent_rows,
            "variants": variant_rows,
            "packet_sizes": packet_sizes,
            "packet_count": packet_count,
            "variant_count": len(variants),
            "family_size": len(pole_indices),
            "projection": method,
            "pole_label": pole_text,
            "packet_plane_label": packet_text,
            "columns": [column.to_json() for column in _POLE_COLUMNS],
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": name,
            "pole": [int(value) for value in request["pole"]],
            "packet_plane": [int(value) for value in request["packet_plane"]],
            "projection": method,
            "include_parent": bool(request["include_parent"]),
        },
        notes=(
            (
                f"The packets are uneven ({size_text} variants), which means "
                f"{packet_text} is not the family the {display} relationship is built on. The "
                "grouping is still the nearest-parallel parent plane of each variant, but it "
                "does not carry the packet meaning it does for the defining family."
            )
            if uneven
            else (
                f"Every packet holds {size_text} variants, as it must when the grouping plane "
                "is the family the relationship is defined on."
            ),
            "Every variant is drawn with equal weight. A real microstructure shows variant "
            "selection, which is a property of the transformation conditions rather than of the "
            "crystallography; this figure is the no-selection baseline.",
            "Poles in the lower hemisphere are folded onto the upper one, because a plane and "
            "its opposite normal describe the same plane.",
        ),
        citations=(_CITATION_MORITO, _CITATION_RANDLE),
    )
    return result.to_json()


@REGISTRY.operation(
    "variants.intervariant_misorientations",
    title="Intervariant misorientation spectrum",
    summary="Every boundary two child grains of one parent can make with each other.",
    help_text=(
        "Lists the disorientation of every pair of transformation variants. Two child grains "
        "that grew from the same parent cannot meet at an arbitrary misorientation: the set of "
        "possible boundaries is fixed by the relationship and by the two point groups, and it "
        "is discrete.\n\n"
        "That is what makes this table useful rather than merely descriptive. A measured "
        "misorientation histogram from an EBSD map of prior-austenite grains should show peaks "
        "at these angles and nowhere else; peaks away from them are boundaries between "
        "*different* parent grains. It is the same reasoning parent-grain reconstruction rests "
        "on.\n\n"
        "**Same-packet pairs are the low-angle ones.** Variants sharing a parent habit plane "
        "differ by small rotations about that plane normal, so they make the sub-block "
        "boundaries a micrograph shows as faint lines within a packet. Cross-packet pairs give "
        "the high-angle boundaries that delimit blocks.\n\n"
        "**Disorientation, not misorientation.** The angle reported is the minimum over the "
        "child point group, which is the convention EBSD software reports and the only one under "
        "which a histogram is comparable between tools."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text="The phase that transforms.",
            builtin="austenite_fcc",
        ),
        phase_parameter(
            name="child_phase",
            label="Child phase",
            help_text="The product phase, whose symmetry reduces each misorientation.",
            builtin="fe_bcc",
        ),
        ChoiceParameter(
            name="relationship",
            label="Orientation relationship",
            help_text="Which relationship generates the variants.",
            options=_RELATIONSHIPS,
            default="kurdjumov_sachs",
        ),
        IndicesParameter(
            name="packet_plane",
            label="Parent plane defining packets",
            help_text=(
                "Used only to label each pair as within or across a packet. (111) for the "
                "fcc-to-bcc relationships, (110) for Burgers."
            ),
            default=(1, 1, 1),
        ),
        BooleanParameter(
            name="merge_equal_angles",
            label="Report the spectrum, not every pair",
            help_text=(
                "On, pairs at the same angle are collapsed to one row carrying the count, which "
                "is the discrete spectrum to compare a histogram against. Off lists all pairs, "
                "which is what is wanted when a particular pair is the question."
            ),
            default=False,
            advanced=True,
        ),
    ),
    returns="One row per variant pair, or per distinct angle; the spectrum under `data`.",
    panel="variants",
    citations=(_CITATION_MORITO, _CITATION_BUNGE),
    tags=(
        "variant",
        "misorientation",
        "disorientation",
        "boundary",
        "packet",
        "block",
        "EBSD",
        "parent reconstruction",
    ),
)
def _intervariant_misorientations(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.core.transformation import intervariant_misorientations

    parent_spec, parent_phase = phase_from_request(request["phase"])
    child_spec, child_phase = phase_from_request(request["child_phase"])
    name = str(request["relationship"])
    relationship = _relationship(name, parent_phase, child_phase)
    variants = relationship.generate_variants()
    if len(variants) < 2:
        raise InvalidInputError(
            f"The {relationship_name(name)} relationship gives only {len(variants)} variant, so "
            "there are no variant pairs to report.",
            field="relationship",
            hint="Bain gives 3 variants; Kurdjumov-Sachs gives 24.",
        )
    packets = _packet_labels(relationship, parent_spec, tuple(request["packet_plane"]))
    by_variant = {
        variant.variant_index: int(packets[position]) + 1
        for position, variant in enumerate(variants)
    }

    pairs = intervariant_misorientations(relationship, variants=variants)
    rows: list[dict[str, Any]] = []
    for pair in pairs:
        axis = np.asarray(pair.axis_child_frame, dtype=float)
        same = by_variant[pair.variant_a] == by_variant[pair.variant_b]
        label, deviation = _nearest_axis(axis, child_phase, child_spec)
        rows.append(
            {
                "variant_a": int(pair.variant_a),
                "variant_b": int(pair.variant_b),
                "angle_deg": float(pair.angle_deg),
                "axis": label,
                "axis_deviation_deg": deviation,
                "axis_x": float(axis[0]),
                "axis_y": float(axis[1]),
                "axis_z": float(axis[2]),
                "same_packet": "yes" if same else "no",
            }
        )
    rows.sort(key=lambda row: float(row["angle_deg"]))

    # The spectrum: distinct angles with their multiplicity. Rounding to three
    # decimals merges representatives that differ only at the floating-point
    # floor of chained matrix products, which is the whole point of asking for
    # a spectrum rather than for a list of pairs.
    spectrum: dict[float, dict[str, Any]] = {}
    for row in rows:
        key = round(float(row["angle_deg"]), 3)
        entry = spectrum.setdefault(
            key, {"angle_deg": key, "pairs": 0, "same_packet": 0, "cross_packet": 0}
        )
        entry["pairs"] += 1
        if row["same_packet"] == "yes":
            entry["same_packet"] += 1
        else:
            entry["cross_packet"] += 1
    spectrum_rows = sorted(spectrum.values(), key=lambda entry: float(entry["angle_deg"]))

    if bool(request["merge_equal_angles"]):
        columns: tuple[Column, ...] = (
            Column("angle_deg", "Disorientation", units="°", numeric=True, digits=3),
            Column("pairs", "Pairs", numeric=True),
            Column("same_packet", "Within a packet", numeric=True),
            Column("cross_packet", "Across packets", numeric=True),
        )
        table_rows: tuple[dict[str, Any], ...] = tuple(spectrum_rows)
        caption = f"Distinct disorientations among the {len(variants)} variants."
    else:
        columns = _PAIR_COLUMNS
        table_rows = tuple(rows)
        caption = f"Disorientation of every pair of the {len(variants)} variants."

    angles = np.asarray([row["angle_deg"] for row in rows], dtype=float)
    low_angle = int(np.count_nonzero(angles < 15.0))
    low_angle_within = sum(
        1 for row in rows if float(row["angle_deg"]) < 15.0 and row["same_packet"] == "yes"
    )
    # The within-packet spectrum is the interesting subset and is smaller than
    # the whole: for Kurdjumov-Sachs it is three angles out of ten. Saying "the
    # low-angle pairs are the within-packet ones" would be the tidy sentence and
    # the wrong one — 24 of the 48 pairs below 15 degrees cross a packet.
    within = sorted(
        {entry["angle_deg"] for entry in spectrum_rows if entry["same_packet"]},
    )
    within_pairs = sum(int(entry["same_packet"]) for entry in spectrum_rows)
    within_text = (
        "no pair shares a packet under this grouping"
        if not within
        else (
            f"the {within_pairs} pairs that share a packet take only "
            + ("one disorientation, " if len(within) == 1 else f"{len(within)} of them, ")
            + ", ".join(f"{value:g}°" for value in within)
        )
    )
    display = relationship_name(name)
    result = AppResult(
        title=f"{display}: misorientations among {len(variants)} variants",
        summary=(
            f"{len(rows)} variant pairs of {child_spec.name} from one {parent_spec.name} grain, "
            f"falling at {len(spectrum_rows)} distinct disorientations between "
            f"{angles.min():.2f}° and {angles.max():.2f}° — a discrete spectrum, not a spread. "
            f"Of those, {within_text}. "
            f"{low_angle} pairs sit below 15°, of which {low_angle_within} share a packet — "
            "sharing a packet and being low-angle are related but not the same thing."
        ),
        table=ResultTable(columns=columns, rows=table_rows, caption=caption),
        data={
            "pairs": rows,
            "spectrum": spectrum_rows,
            "variant_count": len(variants),
            "pair_count": len(rows),
            "distinct_angles": len(spectrum_rows),
            "min_angle_deg": float(angles.min()),
            "max_angle_deg": float(angles.max()),
            "low_angle_pairs": low_angle,
            "packets": by_variant,
            "columns": [column.to_json() for column in columns],
        },
        inputs={
            "phase": parent_spec.to_json(),
            "child_phase": child_spec.to_json(),
            "relationship": name,
            "packet_plane": [int(value) for value in request["packet_plane"]],
            "merge_equal_angles": bool(request["merge_equal_angles"]),
        },
        notes=(
            "Angles are disorientations: the minimum over the child point group, which is the "
            "convention EBSD software reports.",
            "This is the set of boundaries variants of one parent *can* make. Which of them a "
            "microstructure actually shows depends on variant selection and on which variants "
            "happen to be neighbours.",
        ),
        citations=(_CITATION_MORITO, _CITATION_BUNGE),
    )
    return result.to_json()


@REGISTRY.operation(
    "variants.render",
    title="Publication figure of the pole figure",
    summary="Render the variant pole figure through the publication renderer.",
    help_text=(
        "Draws the same poles as the interactive figure through "
        "`pytex.plotting.spherical` in the journal style, with the Wulff net, the packet "
        "colours, and the parent poles marked.\n\n"
        "The screen figure and this one are the same numbers seen through two renderers: the "
        "browser draws from the pole table, and so does this. Nothing is recomputed, so the "
        "figure cannot disagree with the CSV.\n\n"
        "**Choosing a format.** A pole figure is line art and a few dozen markers, so SVG is "
        "small, stays editable, and is the better artefact — the opposite of the crystal "
        "viewer, where a sphere mesh makes SVG the expensive choice. PNG at 600 dpi is offered "
        "for submission systems that will not take vector art."
    ),
    parameters=(
        phase_parameter(
            label="Parent phase",
            help_text="The phase that transforms, and the frame the figure is drawn in.",
            builtin="austenite_fcc",
        ),
        phase_parameter(
            name="child_phase",
            label="Child phase",
            help_text="The product phase.",
            builtin="fe_bcc",
        ),
        ChoiceParameter(
            name="relationship",
            label="Orientation relationship",
            help_text="Which relationship generates the variants.",
            options=_RELATIONSHIPS,
            default="kurdjumov_sachs",
        ),
        IndicesParameter(
            name="pole",
            label="Child plane to plot",
            help_text="The child plane whose symmetry family is projected.",
            default=(1, 0, 0),
        ),
        IndicesParameter(
            name="packet_plane",
            label="Parent plane defining packets",
            help_text="The parent family the variants are grouped and coloured by.",
            default=(1, 1, 1),
        ),
        ChoiceParameter(
            name="projection",
            label="Projection",
            help_text="How the sphere is flattened onto the disc.",
            options=(
                ("stereographic", "Stereographic", "Preserves angles."),
                ("equal_area", "Equal area (Schmidt)", "Preserves area."),
            ),
            default="stereographic",
        ),
        BooleanParameter(
            name="include_parent",
            label="Include the parent's own poles",
            help_text="Draws the same indices in the parent crystal, unrotated.",
            default=True,
        ),
        ChoiceParameter(
            name="format",
            label="Format",
            help_text="SVG for an editable vector figure; PNG for a raster one.",
            options=(
                ("svg", "SVG", "Vector, editable, and small for line art like this."),
                ("png", "PNG", "Raster, for submission systems that refuse vector art."),
            ),
            default="svg",
        ),
        IntegerParameter(
            name="dpi",
            label="Resolution",
            help_text="Pixels per inch, for PNG only.",
            default=600,
            minimum=72,
            maximum=1200,
            advanced=True,
        ),
    ),
    returns="The encoded image under `data.image`, with its format and encoding.",
    panel="variants",
    citations=(_CITATION_MORITO, _CITATION_RANDLE),
    tags=("variant", "pole figure", "figure", "export", "publication"),
)
def _variant_render(request: dict[str, Any]) -> dict[str, Any]:
    try:
        import matplotlib
    except ImportError as error:
        raise DependencyMissingError(
            "matplotlib", purpose="Rendering a publication figure", extra="plotting"
        ) from error

    matplotlib.use("Agg", force=False)
    import base64
    import io

    import matplotlib.pyplot as plt

    from pytex.plotting.spherical import plot_stereographic_vectors

    # The figure is drawn from the same computation the table reports, called
    # here rather than reimplemented, so the published figure and the exported
    # numbers cannot describe different poles.
    computed = _variant_pole_figure(
        {
            "phase": request["phase"],
            "child_phase": request["child_phase"],
            "relationship": request["relationship"],
            "pole": request["pole"],
            "packet_plane": request["packet_plane"],
            "projection": request["projection"],
            "include_parent": request["include_parent"],
        }
    )
    data = computed["data"]
    image_format = str(request["format"])
    dpi = int(request["dpi"])

    vectors: list[list[float]] = []
    colors: list[str] = []
    for pole in data["poles"]:
        vectors.append(_unit_from_angles(pole["polar_deg"], pole["azimuth_deg"]))
        colors.append(_PACKET_COLORS[(int(pole["packet"]) - 1) % len(_PACKET_COLORS)])
    for pole in data["parent_poles"]:
        vectors.append(_unit_from_angles(pole["polar_deg"], pole["azimuth_deg"]))
        colors.append(_PARENT_COLOR)

    axes = plot_stereographic_vectors(
        np.asarray(vectors, dtype=float),
        colors=colors,
        method=str(request["projection"]),
        render="pole",
        include_wulff_net=True,
        title=computed["title"],
    )
    figure = axes.get_figure()
    buffer = io.BytesIO()
    try:
        figure.savefig(buffer, format=image_format, dpi=dpi, bbox_inches="tight")
    finally:
        # A leaked figure is a defect here and a memory leak in a long-running
        # server, which is the same rule the crystal renderer follows.
        plt.close(figure)
    payload = buffer.getvalue()
    image = (
        payload.decode("utf-8")
        if image_format == "svg"
        else base64.b64encode(payload).decode("ascii")
    )

    result = AppResult(
        title=f"{computed['title']}: publication figure",
        summary=(
            f"{image_format.upper()} figure of {len(vectors)} poles "
            f"({len(data['poles'])} from {data['variant_count']} variants in "
            f"{data['packet_count']} packets, {len(data['parent_poles'])} from the parent), "
            f"rendered through the journal style"
            + (f" at {dpi} dpi" if image_format == "png" else "")
            + f" ({len(payload) / 1024:.0f} kB). The poles are the same rows the table exports."
        ),
        data={
            "image": image,
            "format": image_format,
            "encoding": "text" if image_format == "svg" else "base64",
            "bytes": len(payload),
            "dpi": dpi,
            "pole_count": len(vectors),
        },
        inputs=dict(computed["inputs"], format=image_format, dpi=dpi),
        citations=(_CITATION_MORITO, _CITATION_RANDLE),
    )
    return result.to_json()


REGISTRY.add_examples(
    (
        ExampleScenario(
            id="variants.example.ks_packets",
            title="The 24 variants, and the 4 packets they fall into",
            panel="variants",
            summary="Kurdjumov-Sachs {100} poles of austenite-to-ferrite, coloured by packet.",
            teaches=(
                "Twenty-four variants put seventy-two poles on the figure, and they are not "
                "scattered: colouring by packet shows four groups of six, one per member of the "
                "parent {111} family. That grouping is what a lath martensite micrograph shows "
                "as a block, and it is why one parent grain gives twenty-four orientations but "
                "only four apparent plate directions."
            ),
            operation="variants.pole_figure",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "pole": [1, 0, 0],
                "packet_plane": [1, 1, 1],
            },
        ),
        ExampleScenario(
            id="variants.example.nw_twelve",
            title="Nishiyama-Wassermann: half as many variants",
            panel="variants",
            summary="The same figure under the 12-variant relationship, for comparison.",
            teaches=(
                "Run this straight after the Kurdjumov-Sachs figure. The same parent, the same "
                "poles plotted, half the variants — and the packets still number four, because "
                "packets are counted by the parent {111} family, which has four members "
                "whatever relationship sits on it. Nishiyama-Wassermann simply puts three "
                "variants in each rather than six."
            ),
            operation="variants.pole_figure",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "nishiyama_wassermann",
                "pole": [1, 0, 0],
                "packet_plane": [1, 1, 1],
            },
        ),
        ExampleScenario(
            id="variants.example.bain_three",
            title="Bain: three variants and nothing to choose between",
            panel="variants",
            summary="The pure-strain path, whose three variants are the cube axes themselves.",
            teaches=(
                "The Bain path has no shear and only three variants, and its {100} poles land "
                "on the parent's own cube poles. That is why Bain is the reference the "
                "shear-carrying relationships are read against: it shows what the lattice "
                "correspondence alone requires, before any rotation is added to make a plane "
                "and a direction match."
            ),
            operation="variants.pole_figure",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "bain",
                "pole": [1, 0, 0],
                "packet_plane": [1, 1, 1],
            },
        ),
        ExampleScenario(
            id="variants.example.burgers_poles",
            title="Burgers: six packets of two in zirconium",
            panel="variants",
            summary="The bcc-to-hcp path, with the basal pole of every variant plotted.",
            teaches=(
                "Burgers is grouped on the parent {110} family, which has six members, so the "
                "twelve variants fall into six packets of two rather than four of six. Plotting "
                "the child (0001) shows it directly: each packet contributes one basal pole "
                "position, because both of its variants share the parent plane the basal plane "
                "lies on."
            ),
            operation="variants.pole_figure",
            request={
                "phase": {"builtin": "zr_bcc_beta"},
                "child_phase": {"builtin": "zr_hcp"},
                "relationship": "burgers",
                "pole": [0, 0, 1],
                "packet_plane": [1, 1, 0],
            },
        ),
        ExampleScenario(
            id="variants.example.ks_spectrum",
            title="The boundaries 24 variants can make with each other",
            panel="variants",
            summary="The Kurdjumov-Sachs intervariant disorientation spectrum.",
            teaches=(
                "The 276 variant pairs do not spread over the angle range: they collapse onto a "
                "handful of discrete disorientations. That discreteness is the test — a "
                "measured misorientation histogram from prior-austenite grains should show "
                "peaks here and nowhere else, and a peak somewhere else is a boundary between "
                "two different parent grains."
            ),
            operation="variants.intervariant_misorientations",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "packet_plane": [1, 1, 1],
                "merge_equal_angles": True,
            },
        ),
        ExampleScenario(
            id="variants.example.ks_pairs",
            title="Which pairs are the low-angle ones",
            panel="variants",
            summary="Every Kurdjumov-Sachs variant pair, sorted by disorientation.",
            teaches=(
                "Sorted by angle, the top of the table is almost entirely same-packet pairs: "
                "variants sharing a parent {111} differ by small rotations about it, and those "
                "are the sub-block boundaries within a packet. Scroll down and the same-packet "
                "column turns to 'no' — those are the block boundaries a micrograph shows as "
                "sharp lines."
            ),
            operation="variants.intervariant_misorientations",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "child_phase": {"builtin": "fe_bcc"},
                "relationship": "kurdjumov_sachs",
                "packet_plane": [1, 1, 1],
                "merge_equal_angles": False,
            },
        ),
    )
)
