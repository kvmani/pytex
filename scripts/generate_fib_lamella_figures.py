# ruff: noqa: RUF001
"""Generate the canonical frame-convention figure of FIB lamella planning.

The figure is a *generated asset*: the drawn lamella --- its normal ``n_L``, its
long axis ``t_L`` and the target direction ``d*`` at ``eps*`` out of the surface ---
is computed by :func:`pytex.fib.lamella_geometry` for a constructed orientation,
so the convention on the page cannot drift from the one the code implements. The
SVG is written by hand rather than by matplotlib, so its bytes depend on nothing
but this script and are compared exactly by
``tests/unit/test_fib_lamella_figures.py``.

Usage::

    python scripts/generate_fib_lamella_figures.py

Output (tracked as a canonical documentation asset):

- ``docs/figures/fib_lamella_frames.svg`` -- the sample frame, the vertical
  lamella with the in-plane constraint on its normal, the mill direction, and the
  SEM image, FIB ion-view and TEM holder frames the azimuth is carried into.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from pytex.core import frame_catalog
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import lamella_geometry, target_orbit

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = REPO_ROOT / "docs" / "figures"
FRAMES_FIGURE = FIGURES_DIR / "fib_lamella_frames.svg"

#: The drawn example: the target rises EPS_DEG out of the surface, its in-plane
#: part at AZIMUTH_DEG from X_s. Chosen so every label has room.
EPS_DEG = 20.0
AZIMUTH_DEG = 140.0

_WIDTH = 1000
_HEIGHT = 690
_ORIGIN = np.array([280.0, 330.0])
_SCALE = 190.0

_INK = "#07122f"
_MUTED = "#40506f"
_PAPER = "#fbfdff"
_BLUE = "#2563eb"
_TEAL = "#0f9f9f"
_VIOLET = "#7c3aed"
_AMBER = "#f59e0b"
_ROSE = "#e11d48"
_SURFACE = "#e2e8f0"
_SLAB = "#93c5fd"


def example_orientation() -> np.ndarray:
    """A crystal-to-specimen rotation putting ``[001]`` at the drawn direction."""

    eps, azimuth = math.radians(EPS_DEG), math.radians(AZIMUTH_DEG)
    target = np.array(
        [math.cos(eps) * math.cos(azimuth), math.cos(eps) * math.sin(azimuth), math.sin(eps)]
    )
    z = np.array([0.0, 0.0, 1.0])
    axis = np.cross(z, target)
    angle = math.acos(float(z @ target))
    return np.asarray(Rotation.from_axis_angle(axis / np.linalg.norm(axis), angle).as_matrix())


def example_geometry():  # type: ignore[no-untyped-def]
    """The lamella the figure draws, solved by the library."""

    frame = frame_catalog.crystal_frame()
    phase = Phase(
        "tetragonal",
        lattice=Lattice(4.0, 4.0, 5.0, 90.0, 90.0, 90.0, crystal_frame=frame),
        symmetry=SymmetrySpec.from_point_group("4/mmm", reference_frame=frame),
        crystal_frame=frame,
    )
    return lamella_geometry(example_orientation(), target_orbit(phase, (0, 0, 1)))


def _project(point: np.ndarray) -> np.ndarray:
    """An oblique view: X_s towards the viewer and left, Y_s right, Z_s up."""

    x, y, z = (float(value) for value in point)
    return _ORIGIN + _SCALE * np.array([0.92 * y - 0.55 * x, -0.95 * z + 0.30 * x + 0.12 * y])


def _fmt(value: float) -> str:
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text


def _points(points: list[np.ndarray]) -> str:
    return " ".join(f"{_fmt(p[0])},{_fmt(p[1])}" for p in points)


def _polygon(corners: list[np.ndarray], fill: str, *, opacity: float, stroke: str) -> str:
    return (
        f'<polygon points="{_points([_project(c) for c in corners])}" fill="{fill}" '
        f'fill-opacity="{_fmt(opacity)}" stroke="{stroke}" stroke-width="1"/>'
    )


def _line(
    start: np.ndarray,
    end: np.ndarray,
    color: str,
    *,
    width: float = 1.8,
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


def _flat_line(
    a: tuple[float, float],
    b: tuple[float, float],
    color: str,
    *,
    arrow: str | None = None,
    width: float = 1.8,
    dash: str | None = None,
) -> str:
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
    return _flat_text(
        (float(p[0]), float(p[1])), text, color, size=size, anchor=anchor, italic=italic
    )


def _flat_text(
    at: tuple[float, float],
    text: str,
    color: str,
    *,
    size: int = 15,
    anchor: str = "start",
    italic: bool = False,
    bold: bool = False,
) -> str:
    style = ' font-style="italic"' if italic else ""
    weight = ' font-weight="bold"' if bold else ""
    return (
        f'<text x="{_fmt(at[0])}" y="{_fmt(at[1])}" fill="{color}" font-size="{size}" '
        f'font-family="Arial" text-anchor="{anchor}"{style}{weight}>{text}</text>'
    )


def _arrow_marker(identifier: str, color: str) -> str:
    return (
        f'<marker id="{identifier}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="11" '
        f'markerHeight="11" markerUnits="userSpaceOnUse" orient="auto">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>'
    )


def _card(x: float, y: float, width: float, height: float, title: str) -> list[str]:
    return [
        f'<rect x="{_fmt(x)}" y="{_fmt(y)}" width="{_fmt(width)}" height="{_fmt(height)}" rx="8" '
        f'fill="#ffffff" stroke="#cbd5e1" stroke-width="1.2"/>',
        _flat_text((x + 14, y + 26), title, _INK, size=17, bold=True),
    ]


def frames_svg() -> str:
    """Return the frame-convention figure as SVG text."""

    geometry = example_geometry()
    n = np.asarray(geometry.normal_sample)
    t = np.asarray(geometry.long_axis_sample)
    d = np.asarray(geometry.direction_sample)
    if float(d @ n) < 0.0:
        d = -d
    z = np.array([0.0, 0.0, 1.0])
    ex, ey = np.eye(3)[0], np.eye(3)[1]

    # The surface patch, and the lamella plane (spanned by t and Z) standing in it.
    surface = [
        np.array(c) for c in ((-0.9, -0.9, 0), (0.9, -0.9, 0), (0.9, 0.9, 0), (-0.9, 0.9, 0))
    ]
    half_l, half_w, depth = 0.72, 0.06, 0.7
    lamella_plane = [half_l * t, -half_l * t, -half_l * t - depth * z, half_l * t - depth * z]
    footprint = [
        s1 * half_l * t + s2 * half_w * n for s1, s2 in ((-1, -1), (1, -1), (1, 1), (-1, 1))
    ]
    theta = math.radians(geometry.theta_sample_deg)
    theta_arc = [
        _project(0.3 * np.array([math.cos(a), math.sin(a), 0.0]))
        for a in np.linspace(0.0, theta, 24)
    ]
    eps = math.radians(geometry.eps_deg)
    eps_arc = [
        _project(0.62 * (math.cos(a) * n + math.sin(a) * z)) for a in np.linspace(0.0, eps, 16)
    ]
    mid = theta / 2.0
    left = [
        _polygon(surface, _SURFACE, opacity=0.75, stroke=_MUTED),
        _polygon(lamella_plane, _SLAB, opacity=0.45, stroke=_BLUE),
        _polygon(footprint, _BLUE, opacity=0.35, stroke=_BLUE),
        _line(np.zeros(3), 1.0 * ex, _MUTED, arrow="axis"),
        _line(np.zeros(3), 1.0 * ey, _MUTED, arrow="axis"),
        _line(np.zeros(3), 0.92 * z, _MUTED, arrow="axis"),
        _text(1.08 * ex, "X", _INK, italic=True, offset=(-4.0, 14.0)),
        _text(1.06 * ey, "Y", _INK, italic=True, offset=(12.0, 6.0)),
        _text(0.95 * z, "Z (outward normal)", _INK, anchor="end", italic=True, offset=(-10.0, 0.0)),
        _line(np.zeros(3), 0.95 * n, _ROSE, width=2.6, arrow="normal"),
        _text(0.95 * n, "n (TEM beam)", _ROSE, anchor="start", offset=(10.0, 16.0)),
        _line(np.zeros(3), 0.78 * t, _TEAL, width=2.0, arrow="long"),
        _text(0.8 * t, "t = Z × n", _TEAL, anchor="end", offset=(-8.0, -6.0)),
        _line(np.zeros(3), 0.9 * d, _VIOLET, width=2.2, dash="7 4", arrow="target"),
        _text(0.92 * d, "d* (target zone axis)", _VIOLET, anchor="start", offset=(8.0, -8.0)),
        f'<polyline points="{_points(theta_arc)}" fill="none" stroke="{_ROSE}" '
        'stroke-width="1.5"/>',
        _text(
            0.38 * np.array([math.cos(mid), math.sin(mid), 0.0]),
            "θ",
            _ROSE,
            size=18,
            italic=True,
            offset=(0.0, 20.0),
        ),
        f'<polyline points="{_points(eps_arc)}" fill="none" stroke="{_VIOLET}" '
        'stroke-width="1.5"/>',
        _text(
            0.66 * (math.cos(eps / 2) * n + math.sin(eps / 2) * z),
            "ε*",
            _VIOLET,
            size=18,
            italic=True,
            offset=(14.0, 6.0),
        ),
        _line(0.3 * n + 1.25 * z, 0.3 * n + 0.08 * z, _AMBER, width=3.0, arrow="ion"),
        _text(
            0.3 * n + 1.25 * z, "ion beam mills along −Z", _INK, anchor="start", offset=(10.0, 4.0)
        ),
        _text(
            -half_l * t - 0.5 * depth * z,
            "lamella plane",
            _BLUE,
            size=14,
            anchor="end",
            offset=(-6.0, 0.0),
        ),
        _text(
            -half_l * t - depth * z, "depth D", _MUTED, size=14, anchor="end", offset=(-6.0, 4.0)
        ),
    ]
    left += [
        _flat_text(
            (24, 36),
            "FIB lamella planning: the frames and the in-plane constraint",
            _INK,
            size=20,
            bold=True,
        ),
        _flat_text(
            (24, 60),
            "A strictly vertical mill leaves n in the surface plane, so the "
            "TEM beam can only lie along in-plane directions.",
            _MUTED,
            size=15,
        ),
        _flat_text(
            (24, _HEIGHT - 52),
            f"Drawn for a target rising ε* = {geometry.eps_deg:.0f}° "
            f"out of the surface with the lamella normal at θ = "
            f"{geometry.theta_sample_deg:.0f}°, as pytex.fib.lamella_geometry solves it.",
            _MUTED,
            size=15,
        ),
        _flat_text(
            (24, _HEIGHT - 28),
            "ε* = arcsin |d*·Z| is the residual tilt the TEM "
            "holder must supply; (n, t, Z) is the right-handed lamella frame.",
            _MUTED,
            size=15,
        ),
    ]

    # Right column: the three frames the azimuth is carried into.
    x0, width = 650.0, 330.0
    right: list[str] = []
    right += _card(x0, 90, width, 140, "SEM image frame I")
    right += [
        _flat_line((x0 + 30, 138), (x0 + 96, 138), _INK, arrow="axis2"),
        _flat_line((x0 + 30, 138), (x0 + 30, 196), _INK, arrow="axis2"),
        _flat_text((x0 + 102, 143), "x", _INK, italic=True),
        _flat_text((x0 + 36, 214), "y (rows)", _INK, italic=True),
        _flat_text((x0 + 130, 146), "q = A p + t", _INK, size=15),
        _flat_text((x0 + 130, 168), "identity by default;", _MUTED, size=14),
        _flat_text((x0 + 130, 188), "an affine fit reports", _MUTED, size=14),
        _flat_text((x0 + 130, 208), "its residual", _MUTED, size=14),
    ]
    right += _card(x0, 246, width, 170, "FIB ion view P (stage at T)")
    right += [
        _flat_text((x0 + 16, 288), "θion = sR (θ + Rstage) + R0", _INK, size=15),
        _flat_text((x0 + 16, 314), "sR and R0 are calibrated once with", _ROSE, size=14),
        _flat_text((x0 + 16, 334), "fiducial trenches, never assumed.", _ROSE, size=14),
        _flat_text((x0 + 16, 358), "T = 52° or 54°: at that tilt the ion view", _MUTED, size=14),
        _flat_text((x0 + 16, 378), "is the plan view, and the SEM view is", _MUTED, size=14),
        _flat_text((x0 + 16, 398), "foreshortened by cos T.", _MUTED, size=14),
    ]
    right += _card(x0, 432, width, 150, "TEM holder frame H")
    cx, cy = x0 + 62, 512.0
    arc = [
        (cx + 32 * math.cos(a), cy + 32 * math.sin(a))
        for a in np.linspace(0.3, 2 * math.pi - 0.3, 30)
    ]
    right += [
        f'<circle cx="{_fmt(cx)}" cy="{_fmt(cy)}" r="7" fill="{_ROSE}"/>',
        f'<polyline points="{" ".join(f"{_fmt(a)},{_fmt(b)}" for a, b in arc)}" fill="none" '
        f'stroke="{_VIOLET}" stroke-width="1.8" marker-end="url(#phi)"/>',
        _flat_text((cx, cy + 56), "φ unknown", _VIOLET, size=14, anchor="middle"),
        _flat_text((x0 + 120, 478), "n lands on the beam; the", _INK, size=15),
        _flat_text((x0 + 120, 500), "in-plane rotation φ and the", _MUTED, size=14),
        _flat_text((x0 + 120, 520), "front/back flip are unknown", _MUTED, size=14),
        _flat_text((x0 + 120, 540), "until it is mounted, so", _MUTED, size=14),
        _flat_text((x0 + 120, 560), "feasibility is over all φ.", _MUTED, size=14),
    ]

    markers = "".join(
        _arrow_marker(identifier, color)
        for identifier, color in (
            ("axis", _MUTED),
            ("axis2", _INK),
            ("normal", _ROSE),
            ("long", _TEAL),
            ("target", _VIOLET),
            ("ion", _AMBER),
            ("phi", _VIOLET),
        )
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {_WIDTH} {_HEIGHT}" '
        f'width="{_WIDTH}" height="{_HEIGHT}" role="img" aria-labelledby="fib-frames-title '
        'fib-frames-desc" font-family="Arial">\n'
        '<title id="fib-frames-title">FIB lamella planning: sample, lamella, SEM image, FIB ion '
        "view and TEM holder frames</title>\n"
        '<desc id="fib-frames-desc">The sample frame X, Y in the surface and Z the outward normal; '
        "a vertically milled lamella whose normal n lies in the surface at azimuth theta, with "
        "long axis t = Z x n; the ion beam milling along -Z; the target zone axis d* at eps* "
        "out of the surface, which is the residual TEM tilt; and the SEM image, calibrated "
        "FIB ion-view and TEM holder frames the azimuth is carried into, the last with the "
        "unknown mounting rotation phi.</desc>\n"
        f"<defs>{markers}</defs>\n"
        f'<rect width="{_WIDTH}" height="{_HEIGHT}" fill="{_PAPER}"/>\n'
        + "\n".join(left + right)
        + "\n</svg>\n"
    )


def main() -> int:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    FRAMES_FIGURE.write_text(frames_svg(), encoding="utf-8", newline="\n")
    print(f"wrote {FRAMES_FIGURE.relative_to(REPO_ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
