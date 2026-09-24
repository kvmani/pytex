# ruff: noqa: RUF001
"""Generate the canonical figure of the sin^2(psi) measurement geometry.

The figure is a *generated asset*: the scattering vector, its projection onto
the surface and both angle arcs are computed from
:func:`pytex.diffraction.xrd_residual_stress.measurement_direction`, the same
function the analysis uses, so the drawn convention cannot drift from the
implemented one. The SVG is written by hand rather than by matplotlib, so the
bytes depend on nothing but this script and are compared exactly by
``tests/unit/test_residual_stress_figures.py``.

Usage::

    python scripts/generate_residual_stress_figures.py

Output (tracked as a canonical documentation asset):

- ``docs/figures/sin2psi_geometry.svg`` -- the specimen frame S1 S2 S3, the
  scattering vector m at azimuth phi and tilt psi, the diffracting planes
  normal to it, and the two angles.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from pytex.diffraction.xrd_residual_stress import measurement_direction

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = REPO_ROOT / "docs" / "figures"
GEOMETRY_FIGURE = FIGURES_DIR / "sin2psi_geometry.svg"

#: The drawn example direction. Any values would do; these keep the labels apart.
PHI_DEG = 60.0
PSI_DEG = 50.0

_WIDTH = 640
_HEIGHT = 440
_ORIGIN = np.array([300.0, 300.0])
_SCALE = 190.0

_INK = "#111827"
_AXIS = "#374151"
_VECTOR = "#b45309"
_PROJECTION = "#1d4ed8"
_PLANE = "#0f766e"
_SURFACE_FILL = "#e5e7eb"


def _project(point: np.ndarray) -> np.ndarray:
    """An oblique view: S1 towards the viewer and left, S2 right, S3 up."""

    x, y, z = (float(value) for value in point)
    return _ORIGIN + _SCALE * np.array([0.92 * y - 0.55 * x, -0.95 * z + 0.30 * x + 0.12 * y])


def _fmt(value: float) -> str:
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text


def _points(points: list[np.ndarray]) -> str:
    return " ".join(f"{_fmt(p[0])},{_fmt(p[1])}" for p in points)


def _line(
    start: np.ndarray,
    end: np.ndarray,
    color: str,
    *,
    width: float = 1.6,
    dash: str | None = None,
    arrow: str | None = None,
) -> str:
    a, b = _project(start), _project(end)
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    head = f' marker-end="url(#{arrow})"' if arrow else ""
    return (
        f'<line x1="{_fmt(a[0])}" y1="{_fmt(a[1])}" x2="{_fmt(b[0])}" y2="{_fmt(b[1])}" '
        f'stroke="{color}" stroke-width="{_fmt(width)}"{extra}{head}/>'
    )


def _text(
    point: np.ndarray,
    text: str,
    color: str,
    *,
    size: int = 15,
    anchor: str = "middle",
    italic: bool = False,
    offset: tuple[float, float] = (0.0, 0.0),
) -> str:
    p = _project(point) + np.array(offset)
    style = ' font-style="italic"' if italic else ""
    return (
        f'<text x="{_fmt(p[0])}" y="{_fmt(p[1])}" fill="{color}" font-size="{size}" '
        f'text-anchor="{anchor}"{style}>{text}</text>'
    )


def _arrow_marker(identifier: str, color: str) -> str:
    return (
        f'<marker id="{identifier}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="11" '
        f'markerHeight="11" markerUnits="userSpaceOnUse" orient="auto">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>'
    )


def geometry_svg() -> str:
    """Return the geometry figure as SVG text."""

    m = measurement_direction(PHI_DEG, PSI_DEG)
    in_plane = np.array([m[0], m[1], 0.0])
    in_plane_unit = in_plane / np.linalg.norm(in_plane)
    e1, e2, e3 = np.eye(3)

    # The specimen surface, a square patch through the origin.
    corners = [
        np.array(c)
        for c in ((-0.75, -0.75, 0), (0.95, -0.75, 0), (0.95, 0.95, 0), (-0.75, 0.95, 0))
    ]
    surface = (
        f'<polygon points="{_points([_project(c) for c in corners])}" fill="{_SURFACE_FILL}" '
        f'fill-opacity="0.7" stroke="{_AXIS}" stroke-width="0.8"/>'
    )

    # phi arc in the surface, from S1 to the projection of m.
    phi_arc = [
        _project(0.32 * np.array([math.cos(t), math.sin(t), 0.0]))
        for t in np.linspace(0.0, math.radians(PHI_DEG), 24)
    ]
    # psi arc in the tilt plane, from S3 to m.
    psi_arc = [
        _project(0.42 * (math.cos(t) * e3 + math.sin(t) * in_plane_unit))
        for t in np.linspace(0.0, math.radians(PSI_DEG), 24)
    ]

    # The diffracting planes: short traces normal to m, drawn through points on m.
    normal_a = np.cross(m, e3)
    normal_a /= np.linalg.norm(normal_a)
    normal_b = np.cross(m, normal_a)
    planes = []
    for fraction in (0.62, 0.74, 0.86):
        centre = fraction * m
        quad = [
            centre + 0.13 * (sa * normal_a + sb * normal_b)
            for sa, sb in ((-1, -1), (1, -1), (1, 1), (-1, 1))
        ]
        planes.append(
            f'<polygon points="{_points([_project(q) for q in quad])}" fill="{_PLANE}" '
            f'fill-opacity="0.16" stroke="{_PLANE}" stroke-width="0.9"/>'
        )

    body = [
        f'<rect width="{_WIDTH}" height="{_HEIGHT}" fill="#ffffff"/>',
        surface,
        _line(-0.0 * e1, 1.05 * e1, _AXIS, arrow="axis"),
        _line(-0.0 * e2, 1.05 * e2, _AXIS, arrow="axis"),
        _line(-0.0 * e3, 1.05 * e3, _AXIS, arrow="axis"),
        _text(1.13 * e1, "S1", _INK, offset=(-6.0, 10.0)),
        _text(1.12 * e2, "S2", _INK, offset=(10.0, 5.0)),
        _text(1.10 * e3, "S3 (surface normal)", _INK, anchor="start", offset=(8.0, 0.0)),
        _line(np.zeros(3), in_plane, _PROJECTION, width=1.3, dash="6 4"),
        _line(in_plane, m, _PROJECTION, width=0.9, dash="3 3"),
        *planes,
        _line(np.zeros(3), m, _VECTOR, width=2.4, arrow="vector"),
        _text(1.06 * m, "m (scattering vector)", _VECTOR, anchor="start", offset=(8.0, -2.0)),
        f'<polyline points="{_points(phi_arc)}" fill="none" stroke="{_PROJECTION}" '
        f'stroke-width="1.5"/>',
        _text(
            0.40
            * np.array(
                [math.cos(math.radians(PHI_DEG / 2)), math.sin(math.radians(PHI_DEG / 2)), 0.0]
            ),
            "φ",
            _PROJECTION,
            size=18,
            italic=True,
            offset=(0.0, 14.0),
        ),
        f'<polyline points="{_points(psi_arc)}" fill="none" stroke="{_VECTOR}" '
        f'stroke-width="1.5"/>',
        _text(
            0.50
            * (
                math.cos(math.radians(PSI_DEG / 2)) * e3
                + math.sin(math.radians(PSI_DEG / 2)) * in_plane_unit
            ),
            "ψ",
            _VECTOR,
            size=18,
            italic=True,
            offset=(8.0, 0.0),
        ),
        _text(
            0.78 * m + 0.16 * normal_a,
            "diffracting (hkl) planes",
            _PLANE,
            size=12,
            anchor="start",
            offset=(34.0, 40.0),
        ),
        f'<text x="20" y="30" fill="{_INK}" font-size="15" font-weight="bold">'
        "The sin²ψ measurement geometry</text>",
        f'<text x="20" y="52" fill="{_AXIS}" font-size="12">'
        "m = (cos φ sin ψ, sin φ sin ψ, cos ψ);  ε(φ,ψ) = (d − d₀)/d₀ = m·ε·m</text>",
        f'<text x="20" y="{_HEIGHT - 18}" fill="{_AXIS}" font-size="12">'
        f"Drawn at φ = {PHI_DEG:g}°, ψ = {PSI_DEG:g}°. φ is measured in the surface from S1 "
        "towards S2; ψ from S3.</text>",
    ]
    markers = _arrow_marker("axis", _AXIS) + _arrow_marker("vector", _VECTOR)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {_WIDTH} {_HEIGHT}" '
        f'width="{_WIDTH}" height="{_HEIGHT}" role="img" '
        'aria-labelledby="sin2psi-title" font-family="DejaVu Sans, Arial, sans-serif">\n'
        '<title id="sin2psi-title">The sin²ψ measurement geometry: the specimen frame, the '
        "scattering vector at azimuth φ and tilt ψ, and the planes that diffract</title>\n"
        f"<defs>{markers}</defs>\n" + "\n".join(body) + "\n</svg>\n"
    )


def main() -> int:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    GEOMETRY_FIGURE.write_text(geometry_svg(), encoding="utf-8", newline="\n")
    print(f"wrote {GEOMETRY_FIGURE.relative_to(REPO_ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
