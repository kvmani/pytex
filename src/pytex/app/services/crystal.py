"""The crystal viewer: a structure you can turn in your hands.

What it does
    Builds the three-dimensional scene of a structure — atoms, bonds, the unit
    cell, and any number of superimposed planes and directions — and sends it to
    the browser as geometry in Cartesian angstrom. The browser holds a camera and
    redraws on drag; nothing round-trips while the mouse is down.

When to use it
    To see what a plane or a direction actually *is* in a structure: which atoms
    a slip plane cuts, how a Burgers direction lies in it, how a family of planes
    is arranged around the cell. The calculator gives the numbers; this gives the
    picture the numbers describe.

The division of labour
    Every crystallographic decision is made here, in Python: which atoms exist,
    where the plane polygon is clipped, what a direction arrow's endpoints are,
    what colour a species takes. The browser receives finished vertices and
    applies a rotation and an orthographic divide — under fifty lines of viewing
    arithmetic, no crystallography. When a figure is exported, the camera comes
    back and :mod:`pytex.plotting.crystal3d` renders the same scene at
    publication quality, so the picture on screen and the figure in the paper are
    the same geometry seen the same way.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from pytex.app.errors import DependencyMissingError, InvalidInputError
from pytex.app.phases import PhaseSpec, phase_from_request
from pytex.app.registry import (
    REGISTRY,
    BooleanParameter,
    ChoiceParameter,
    ExampleScenario,
    IndicesListParameter,
    IntegerParameter,
    NumberParameter,
    TextParameter,
)
from pytex.app.results import AppResult, Column, ResultTable
from pytex.app.services.calculator import direction_label, phase_parameter, plane_label
from pytex.core.lattice import CrystalDirection, CrystalPlane, MillerIndex, Phase

__all__ = ["camera_angles_from_matrix", "scene_payload"]

_CITATION_VESTA = (
    "Momma & Izumi, VESTA 3, J. Appl. Crystallogr. 44 (2011) 1272 (visual conventions)."
)

_RENDER_STYLES = (
    (
        "ball_and_stick",
        "Ball and stick",
        "Spheres at covalent-radius scale joined by bonds. The default, and the clearest for "
        "seeing how planes cut a structure.",
    ),
    (
        "space_filling",
        "Space filling",
        "Atomic-radius spheres, bonds suppressed. Shows packing and free volume.",
    ),
    ("stick", "Stick", "Uniform thin bonds with small atoms. Best for open framework structures."),
    ("wireframe", "Wireframe", "Line bonds and marker atoms. Lightest to read behind overlays."),
    (
        "polyhedral",
        "Polyhedral",
        "Coordination polyhedra instead of individual bonds. The mineralogist's view.",
    ),
)


def _plane_overlay(indices: tuple[int, ...], phase: Phase) -> CrystalPlane:
    return CrystalPlane(
        miller=MillerIndex(indices=np.asarray(indices, dtype=int), phase=phase), phase=phase
    )


def _direction_overlay(indices: tuple[int, ...], phase: Phase) -> CrystalDirection:
    return CrystalDirection(coordinates=np.asarray(indices, dtype=float), phase=phase)


def _overlay_label(default: str, replacements: Sequence[str] | None, index: int) -> str:
    """The caller's label for overlay ``index``, or the scene's own."""

    if replacements is None or index >= len(replacements):
        return str(default)
    return str(replacements[index])


def camera_angles_from_matrix(matrix: Any) -> tuple[float, float]:
    """Convert a browser camera matrix into the renderer's elevation and azimuth.

    Purpose
    -------
    The interactive view holds its camera as a 3x3 rotation, because that is
    what a drag composes cleanly; the publication renderer is driven by an
    elevation and an azimuth. This is the one conversion between them, and it
    lives here rather than in JavaScript so that the exported figure and the
    on-screen view cannot disagree through two slightly different derivations.

    The third row of the camera matrix is the direction the viewer looks along,
    expressed in crystal coordinates. Elevation is its angle above the xy-plane
    and azimuth is its bearing within that plane, which is precisely how
    matplotlib's 3-D axes define ``elev`` and ``azim``.

    Parameters
    ----------
    matrix : array-like
        Nine numbers in row-major order, or a ``(3, 3)`` array.

    Returns
    -------
    tuple of (float, float)
        Elevation and azimuth in degrees.

    Raises
    ------
    InvalidInputError
        If the matrix is not nine finite numbers, or its view row is degenerate.
    """

    values = np.asarray(matrix, dtype=float).reshape(-1)
    if values.size != 9 or not np.all(np.isfinite(values)):
        raise InvalidInputError(
            "The camera matrix must be nine finite numbers in row-major order.",
            field="camera_matrix",
            hint="Reset the view and try again.",
        )
    forward = values.reshape(3, 3)[2]
    norm = float(np.linalg.norm(forward))
    if norm < 1e-9:
        raise InvalidInputError(
            "The camera matrix has a degenerate viewing direction.",
            field="camera_matrix",
            hint="Reset the view and try again.",
        )
    forward = forward / norm
    elevation = float(np.degrees(np.arcsin(np.clip(forward[2], -1.0, 1.0))))
    azimuth = float(np.degrees(np.arctan2(forward[1], forward[0])))
    return elevation, azimuth


def scene_payload(
    scene: Any,
    *,
    spec: PhaseSpec,
    plane_labels: Sequence[str] | None = None,
    direction_labels: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Convert a :class:`CrystalScene` into the JSON the browser draws.

    Only what a renderer needs crosses the wire: positions, radii, colours,
    labels, and the polygons already clipped to the cell. Everything is in
    Cartesian angstrom in the crystal frame, so the browser's camera is a pure
    rotation with no crystallography hidden inside it.

    Each atom and each bond carries enough of its own description to answer a
    hover — species, occupancy, bond length — because a scene the user cannot
    interrogate is a picture rather than an instrument.

    Overlay labels are replaced by ``plane_labels`` and ``direction_labels``
    when given, in scene order. The scene builder labels for matplotlib, whose
    mathtext ``$(1\\bar{1}0)$`` is markup rather than text; the browser draws
    the string literally, so it must be handed the plain form the rest of the
    application already uses in tables and prose.
    """

    bounds = np.asarray(scene.bounds(), dtype=float)
    centre = (bounds[0] + bounds[1]) / 2.0
    radius = float(np.linalg.norm(bounds[1] - bounds[0]) / 2.0) or 1.0

    atoms = [
        {
            "position": [float(value) for value in atom.position_angstrom],
            "radius": float(atom.radius_angstrom),
            "color": str(atom.color),
            "species": str(atom.species),
            "occupancy": float(atom.occupancy),
            "label": atom.label,
        }
        for atom in scene.atoms
    ]
    bonds = [
        {
            "start": [float(value) for value in bond.start_angstrom],
            "end": [float(value) for value in bond.end_angstrom],
            "color": str(bond.color),
            "start_color": bond.start_color,
            "end_color": bond.end_color,
            "length": float(bond.length_angstrom),
            "species": f"{bond.start_species or '?'}-{bond.end_species or '?'}",
        }
        for bond in scene.bonds
    ]
    cell_edges = [
        [[float(value) for value in point] for point in edge] for edge in scene.lattice_edges
    ]
    for cell in scene.cells:
        cell_edges.extend(
            [[float(value) for value in point] for point in edge] for edge in cell.edges_angstrom
        )
    planes = [
        {
            "vertices": [[float(value) for value in point] for point in plane.vertices_angstrom],
            "normal": [float(value) for value in plane.normal_angstrom],
            "color": str(plane.color),
            "alpha": float(plane.alpha),
            "label": _overlay_label(plane.label, plane_labels, index),
        }
        for index, plane in enumerate(scene.planes)
    ]
    directions = [
        {
            "start": [float(value) for value in direction.start_angstrom],
            "end": [float(value) for value in direction.end_angstrom],
            "color": str(direction.color),
            "label": _overlay_label(direction.label, direction_labels, index),
        }
        for index, direction in enumerate(scene.directions)
    ]
    return {
        "atoms": atoms,
        "bonds": bonds,
        "cell_edges": cell_edges,
        "planes": planes,
        "directions": directions,
        "centre": [float(value) for value in centre],
        "radius": radius,
        "bounds": [[float(value) for value in row] for row in bounds],
        "phase": spec.to_json(),
        "axes": _axis_arrows(spec),
    }


def _axis_arrows(spec: PhaseSpec) -> list[dict[str, Any]]:
    """The a, b, c axes as arrows, so the view always says which way is which.

    Without them a rotated structure is unreadable: every cubic cell looks the
    same from every direction, and the whole point of turning it is to know
    where you now are.
    """

    phase = spec.to_phase()
    basis = np.asarray(phase.lattice.direct_basis().matrix, dtype=float)
    return [
        {
            "label": name,
            "vector": [float(value) for value in basis[:, index]],
        }
        for index, name in enumerate(("a", "b", "c"))
    ]


@REGISTRY.operation(
    "crystal.scene",
    title="Crystal structure",
    summary="The 3D structure, with any number of planes and directions superimposed.",
    help_text=(
        "Builds the structure as a scene you can rotate: drag to turn it, scroll to zoom, and "
        "hover any atom or bond for its identity and its bond length.\n\n"
        "Enter as many planes and directions as you like — each plane is clipped to the cell and "
        "drawn as a translucent polygon with its indices, and each direction is drawn as a "
        "labelled arrow from the cell origin. Superimposing a slip plane and its Burgers "
        "direction is the canonical use: the arrow must lie *in* the polygon, and seeing that it "
        "does is worth more than checking that the zone-law integer is zero.\n\n"
        "The render style changes the whole visual system at once. Ball-and-stick reads best "
        "behind plane overlays; space-filling shows packing; polyhedral is the mineralogist's "
        "view of a framework structure."
    ),
    parameters=(
        phase_parameter(),
        ChoiceParameter(
            name="render_style",
            label="Style",
            help_text="The visual system for atoms and bonds.",
            options=_RENDER_STYLES,
            default="ball_and_stick",
        ),
        IntegerParameter(
            name="repeat_a",
            label="Repeats along a",
            help_text="How many cells to build along a. More cells show longer-range structure.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
        ),
        IntegerParameter(
            name="repeat_b",
            label="Repeats along b",
            help_text="How many cells to build along b.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
        ),
        IntegerParameter(
            name="repeat_c",
            label="Repeats along c",
            help_text="How many cells to build along c.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
        ),
        IndicesListParameter(
            name="planes",
            label="Planes to superimpose (hkl)",
            help_text=(
                "One plane per row. Each is drawn as a translucent polygon clipped to the cell "
                "and labelled with its indices. Leave empty for none."
            ),
            required=False,
            group="Overlays",
        ),
        IndicesListParameter(
            name="directions",
            label="Directions to superimpose [uvw]",
            help_text="One direction per row, drawn as a labelled arrow from the cell origin.",
            required=False,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_bonds",
            label="Draw bonds",
            help_text=(
                "Bonds are inferred from covalent radii plus a tolerance, not read from a file, "
                "so they are an aid to reading the picture rather than a statement about "
                "chemistry."
            ),
            default=True,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_unit_cells",
            label="Outline every cell",
            help_text="Draw the edges of each repeated cell, not only the outer box.",
            default=False,
            group="Overlays",
        ),
        ChoiceParameter(
            name="atom_labels",
            label="Atom labels",
            help_text="Label every atom with its element or its site name.",
            options=(
                ("none", "None", "No labels; cleanest for dense structures."),
                ("species", "Element", "The element symbol on each atom."),
                ("site", "Site", "The crystallographic site label, as in the CIF."),
            ),
            default="none",
            advanced=True,
            group="Overlays",
        ),
        NumberParameter(
            name="bond_tolerance_angstrom",
            label="Bond tolerance",
            help_text=(
                "Added to the sum of covalent radii when deciding whether two atoms are bonded. "
                "Raise it to catch long bonds; raise it too far and everything bonds to "
                "everything."
            ),
            units="Å",
            default=0.45,
            minimum=0.0,
            maximum=2.0,
            advanced=True,
            group="Overlays",
        ),
    ),
    returns="The scene under `data.scene`; the atom list as the table.",
    panel="crystal",
    citations=(_CITATION_VESTA,),
    tags=("crystal", "structure", "3D", "viewer", "plane", "direction", "unit cell", "bonds"),
)
def _crystal_scene(request: dict[str, Any]) -> dict[str, Any]:
    from pytex.plotting.crystal3d import build_crystal_scene

    spec, phase = phase_from_request(request["phase"])
    if not spec.has_structure:
        raise InvalidInputError(
            f"{spec.name} carries no atomic basis, so there is nothing to draw.",
            field="phase",
            hint=(
                "Choose a built-in phase, or add atomic sites to the phase description. Lattice "
                "geometry alone supports the calculator but not the viewer."
            ),
        )
    repeats = (int(request["repeat_a"]), int(request["repeat_b"]), int(request["repeat_c"]))
    plane_rows = tuple(request.get("planes") or ())
    direction_rows = tuple(request.get("directions") or ())

    scene = build_crystal_scene(
        phase,
        repeats=repeats,
        render_style=str(request["render_style"]),
        show_bonds=bool(request["show_bonds"]),
        bond_tolerance_angstrom=float(request["bond_tolerance_angstrom"]),
        show_unit_cells=bool(request["show_unit_cells"]),
        atom_label_mode=str(request["atom_labels"]),
        plane_overlays=tuple(_plane_overlay(row, phase) for row in plane_rows),
        direction_overlays=tuple(_direction_overlay(row, phase) for row in direction_rows),
    )
    plane_texts = tuple(plane_label(row, spec=spec) for row in plane_rows)
    direction_texts = tuple(direction_label(row, spec=spec) for row in direction_rows)
    payload = scene_payload(
        scene, spec=spec, plane_labels=plane_texts, direction_labels=direction_texts
    )

    rows = tuple(
        {
            "index": index + 1,
            "species": atom["species"],
            "x": atom["position"][0],
            "y": atom["position"][1],
            "z": atom["position"][2],
            "occupancy": atom["occupancy"],
        }
        for index, atom in enumerate(payload["atoms"])
    )
    bond_summary = {
        f"{pair[0]}-{pair[1]}": stats for pair, stats in scene.bond_length_summary().items()
    }
    overlay_text = ", ".join([*plane_texts, *direction_texts])
    result = AppResult(
        title=f"{spec.name}: crystal structure",
        summary=(
            f"{len(rows)} atoms in a {repeats[0]}x{repeats[1]}x{repeats[2]} block of "
            f"{spec.name} ({spec.crystal_system}, {spec.point_group})"
            + (f", with {overlay_text} superimposed" if overlay_text else "")
            + ". Positions are Cartesian angstrom in the crystal frame. Bonds are inferred from "
            "covalent radii plus the stated tolerance, so they aid reading rather than assert "
            "chemistry."
        ),
        table=ResultTable(
            columns=(
                Column("index", "#", numeric=True),
                Column("species", "Element"),
                Column("x", "x", units="Å", numeric=True, digits=5),
                Column("y", "y", units="Å", numeric=True, digits=5),
                Column("z", "z", units="Å", numeric=True, digits=5),
                Column("occupancy", "Occupancy", numeric=True, digits=3),
            ),
            rows=rows,
            caption=f"Atom positions in the displayed block of {spec.name}.",
        ),
        data={
            "scene": payload,
            "bond_summary": bond_summary,
            "repeats": list(repeats),
        },
        inputs={
            "phase": spec.to_json(),
            "render_style": request["render_style"],
            "repeat_a": repeats[0],
            "repeat_b": repeats[1],
            "repeat_c": repeats[2],
            "planes": [list(row) for row in plane_rows],
            "directions": [list(row) for row in direction_rows],
            "show_bonds": bool(request["show_bonds"]),
            "show_unit_cells": bool(request["show_unit_cells"]),
            "atom_labels": request["atom_labels"],
            "bond_tolerance_angstrom": float(request["bond_tolerance_angstrom"]),
        },
        citations=(_CITATION_VESTA,),
    )
    return result.to_json()


@REGISTRY.operation(
    "crystal.render",
    title="Publication figure of the structure",
    summary="Render the current view through the publication renderer.",
    help_text=(
        "Takes the camera you set by dragging and renders the same scene through "
        "`pytex.plotting.crystal3d` in the journal style, with lit spheres and correct "
        "depth ordering.\n\n"
        "The browser view and this figure are the same geometry seen from the same direction, "
        "because the camera angles travel with the request. What you framed on screen is what "
        "the figure shows.\n\n"
        "**Choosing a format.** A lit sphere is a mesh, and a vector format must write every "
        "facet of it: an SVG of a two-cell rock-salt block runs to megabytes and will defeat "
        "most journal submission systems. PNG at 600 dpi is therefore the default, and is the "
        "better publication artefact for ball-and-stick and space-filling views. Choose SVG "
        "when the figure must stay editable, or for the stick and wireframe styles, which are "
        "line art and stay small; the sphere mesh is coarsened automatically for SVG so the "
        "file remains openable."
    ),
    parameters=(
        phase_parameter(),
        ChoiceParameter(
            name="render_style",
            label="Style",
            help_text="The visual system for atoms and bonds.",
            options=_RENDER_STYLES,
            default="ball_and_stick",
        ),
        IntegerParameter(
            name="repeat_a",
            label="Repeats along a",
            help_text="Cells along a.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
        ),
        IntegerParameter(
            name="repeat_b",
            label="Repeats along b",
            help_text="Cells along b.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
        ),
        IntegerParameter(
            name="repeat_c",
            label="Repeats along c",
            help_text="Cells along c.",
            default=1,
            minimum=1,
            maximum=6,
            group="Extent",
        ),
        IndicesListParameter(
            name="planes",
            label="Planes (hkl)",
            help_text="Planes to superimpose, one per row.",
            required=False,
            group="Overlays",
        ),
        IndicesListParameter(
            name="directions",
            label="Directions [uvw]",
            help_text="Directions to superimpose, one per row.",
            required=False,
            group="Overlays",
        ),
        NumberParameter(
            name="elevation_deg",
            label="Camera elevation",
            help_text="Angle above the horizon, matching the on-screen view.",
            units="°",
            default=22.0,
            minimum=-90.0,
            maximum=90.0,
            group="Camera",
        ),
        NumberParameter(
            name="azimuth_deg",
            label="Camera azimuth",
            help_text="Rotation about the vertical, matching the on-screen view.",
            units="°",
            default=34.0,
            group="Camera",
        ),
        BooleanParameter(
            name="show_bonds",
            label="Draw bonds",
            help_text="Include inferred bonds.",
            default=True,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_unit_cells",
            label="Outline every cell",
            help_text="Draw the edges of each repeated cell, not only the outer box.",
            default=False,
            group="Overlays",
        ),
        ChoiceParameter(
            name="atom_labels",
            label="Atom labels",
            help_text="Label every atom with its element or its site name.",
            options=(
                ("none", "None", "No labels; cleanest for dense structures."),
                ("species", "Element", "The element symbol on each atom."),
                ("site", "Site", "The crystallographic site label, as in the CIF."),
            ),
            default="none",
            advanced=True,
            group="Overlays",
        ),
        NumberParameter(
            name="bond_tolerance_angstrom",
            label="Bond tolerance",
            help_text=(
                "Added to the sum of covalent radii when deciding whether two atoms are bonded. "
                "Must match the viewer's value, or the figure will show a different set of "
                "bonds from the one on screen."
            ),
            units="Å",
            default=0.45,
            minimum=0.0,
            maximum=2.0,
            advanced=True,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_legend",
            label="Species legend",
            help_text="Add a colour key for the elements present.",
            default=True,
            group="Overlays",
        ),
        BooleanParameter(
            name="show_frame_indicator",
            label="Axis gizmo",
            help_text=(
                "Add a small indicator showing where a, b and c point in the rendered view, so "
                "the figure states its own orientation."
            ),
            default=True,
            group="Overlays",
        ),
        TextParameter(
            name="camera_matrix",
            label="Camera matrix",
            help_text=(
                "Nine numbers, row-major, describing the on-screen camera. Sent automatically by "
                "the viewer so the figure shows exactly the view you framed; it overrides the "
                "elevation and azimuth above when present."
            ),
            required=False,
            advanced=True,
            group="Camera",
        ),
        ChoiceParameter(
            name="format",
            label="Format",
            help_text="PNG for lit spheres; SVG when the figure must stay editable.",
            options=(
                ("png", "PNG", "Raster at the chosen resolution. Compact, and exact for spheres."),
                ("svg", "SVG", "Vector and editable. Large for sphere-based styles."),
            ),
            default="png",
            group="Output",
        ),
        IntegerParameter(
            name="dpi",
            label="Resolution",
            help_text=(
                "Dots per inch for PNG. 300 is the usual journal minimum; 600 is safe for "
                "anything. Ignored for SVG."
            ),
            units="dpi",
            default=600,
            minimum=72,
            maximum=1200,
            group="Output",
        ),
    ),
    returns="The image under `data.image`: markup for SVG, base64 for PNG.",
    panel="crystal",
    citations=(_CITATION_VESTA,),
    tags=("export", "figure", "SVG", "publication", "crystal"),
)
def _crystal_render(request: dict[str, Any]) -> dict[str, Any]:
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

    from pytex.plotting.crystal3d import plot_crystal_structure_3d

    spec, phase = phase_from_request(request["phase"])
    repeats = (int(request["repeat_a"]), int(request["repeat_b"]), int(request["repeat_c"]))
    plane_rows = tuple(request.get("planes") or ())
    direction_rows = tuple(request.get("directions") or ())
    image_format = str(request["format"])
    dpi = int(request["dpi"])
    elevation = float(request["elevation_deg"])
    azimuth = float(request["azimuth_deg"])
    if request.get("camera_matrix"):
        elevation, azimuth = camera_angles_from_matrix(
            [float(value) for value in str(request["camera_matrix"]).replace(",", " ").split()]
        )

    # A sphere is drawn as a quad mesh, and a vector format writes every facet:
    # at the renderer's default resolution a two-cell block becomes tens of
    # thousands of polygons and an SVG of several megabytes. Coarsening the mesh
    # for SVG keeps the file openable, at a cost invisible at figure size. PNG
    # keeps the full mesh, because there the cost is a fixed pixel count either
    # way.
    style_overrides = (
        {"crystal": {"atom_surface_resolution": 18, "bond_surface_resolution": 16}}
        if image_format == "svg"
        else None
    )

    axes = plot_crystal_structure_3d(
        phase,
        repeats=repeats,
        render_style=str(request["render_style"]),
        show_bonds=bool(request["show_bonds"]),
        bond_tolerance_angstrom=float(request["bond_tolerance_angstrom"]),
        show_unit_cells=bool(request["show_unit_cells"]),
        atom_label_mode=str(request["atom_labels"]),
        plane_overlays=tuple(_plane_overlay(row, phase) for row in plane_rows),
        direction_overlays=tuple(_direction_overlay(row, phase) for row in direction_rows),
        elev_deg=elevation,
        azim_deg=azimuth,
        show_legend=bool(request["show_legend"]),
        show_frame_indicator=bool(request["show_frame_indicator"]),
        style_overrides=style_overrides,
    )
    figure = axes.get_figure()
    buffer = io.BytesIO()
    try:
        figure.savefig(buffer, format=image_format, dpi=dpi, bbox_inches="tight")
    finally:
        # Tests treat a leaked figure as a defect, and a server that leaked one
        # per export would exhaust memory in an afternoon.
        plt.close(figure)
    payload = buffer.getvalue()
    image = (
        payload.decode("utf-8")
        if image_format == "svg"
        else base64.b64encode(payload).decode("ascii")
    )

    result = AppResult(
        title=f"{spec.name}: publication figure",
        summary=(
            f"{image_format.upper()} figure of {spec.name} at elevation {elevation:.1f}° and "
            f"azimuth {azimuth:.1f}°, rendered through the journal style"
            + (f" at {dpi} dpi" if image_format == "png" else "")
            + f" ({len(payload) / 1024:.0f} kB). The geometry is identical to the interactive "
            "view; only the renderer differs."
        ),
        notes=(
            (
                "The sphere mesh is coarsened for SVG so the file stays openable. For a "
                "sphere-based style, PNG at 600 dpi is both smaller and closer to what the "
                "renderer intends.",
            )
            if image_format == "svg"
            else ()
        ),
        data={
            "image": image,
            "format": image_format,
            "encoding": "text" if image_format == "svg" else "base64",
            "bytes": len(payload),
            "dpi": dpi,
            "elevation_deg": elevation,
            "azimuth_deg": azimuth,
        },
        inputs={
            "phase": spec.to_json(),
            "render_style": request["render_style"],
            "repeat_a": repeats[0],
            "repeat_b": repeats[1],
            "repeat_c": repeats[2],
            "planes": [list(row) for row in plane_rows],
            "directions": [list(row) for row in direction_rows],
            "elevation_deg": elevation,
            "azimuth_deg": azimuth,
            "format": image_format,
            "dpi": dpi,
        },
        citations=(_CITATION_VESTA,),
    )
    return result.to_json()


# --------------------------------------------------------------------------
# Canonical examples
# --------------------------------------------------------------------------

REGISTRY.add_examples(
    (
        ExampleScenario(
            id="crystal.example.nacl",
            title="Rock salt: two interpenetrating fcc lattices",
            panel="crystal",
            summary="NaCl with the (100) cube face and the (111) close-packed plane drawn in.",
            teaches=(
                "Rotate until you look down [100]: sodium and chlorine alternate along every "
                "cube edge. The (111) plane, by contrast, cuts through atoms of only one species "
                "at a time — which is why the 111 reflection of NaCl is weak while 200 is strong."
            ),
            operation="crystal.scene",
            request={
                "phase": {"builtin": "nacl"},
                "planes": [[1, 0, 0], [1, 1, 1]],
                "atom_labels": "species",
            },
        ),
        ExampleScenario(
            id="crystal.example.fcc_slip",
            title="The fcc slip system, drawn",
            panel="crystal",
            summary="Austenite with (111) and the [1-10] slip direction superimposed.",
            teaches=(
                "The arrow lies *in* the polygon — that is what h·u + k·v + l·w = 0 looks like. "
                "Turn the structure until the (111) plane is edge-on and the atoms in it line up "
                "into the close-packed rows the dislocation actually moves along."
            ),
            operation="crystal.scene",
            request={
                "phase": {"builtin": "austenite_fcc"},
                "planes": [[1, 1, 1]],
                "directions": [[1, -1, 0]],
                "repeat_a": 2,
                "repeat_b": 2,
                "repeat_c": 2,
            },
        ),
        ExampleScenario(
            id="crystal.example.bcc",
            title="Ferrite and its {110} planes",
            panel="crystal",
            summary="bcc iron with two members of the {110} family drawn in.",
            teaches=(
                "The body-centring atom sits in every (110) plane, not between them — which is "
                "why {110} is the most densely packed plane of a bcc metal and carries most of "
                "its slip, despite bcc having no truly close-packed plane at all."
            ),
            operation="crystal.scene",
            request={
                "phase": {"builtin": "fe_bcc"},
                "planes": [[1, 1, 0], [1, -1, 0]],
                "repeat_a": 2,
                "repeat_b": 2,
                "repeat_c": 2,
            },
        ),
        ExampleScenario(
            id="crystal.example.zr_basal",
            title="Zirconium: basal and prism planes together",
            panel="crystal",
            summary="hcp zirconium with (0001) and a prism plane, plus the a-direction.",
            teaches=(
                "The basal plane is close-packed and the prism plane is not, and the picture "
                "shows why: turn it so the basal plane is edge-on and the ABAB stacking is "
                "visible as two offset layers. That stacking, not symmetry alone, is what makes "
                "prism slip compete with basal slip in zirconium."
            ),
            operation="crystal.scene",
            request={
                "phase": {"builtin": "zr_hcp"},
                "planes": [[0, 0, 1], [1, 0, 0]],
                "directions": [[1, 0, 0]],
                "repeat_a": 2,
                "repeat_b": 2,
                "repeat_c": 2,
            },
        ),
    )
)
