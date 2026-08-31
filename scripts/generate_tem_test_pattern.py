"""Render the tracked TEM test pattern: a simulated hcp plate, as a PNG on disk.

Why a raster fixture exists at all
----------------------------------
The workbench's practice plates are sent to the browser as *coordinates* and
drawn as vectors, which is the right way to deliver a simulation. It means,
however, that nothing in the repository exercises the other path a user takes:
opening an image file from disk, calibrating it, and picking spots off pixels.
That path has its own failure modes — the file never appears, the camera is
inherited from the previous pattern, the pixel read that snaps a click to a
centroid fails — and they are invisible to every test that starts from the
gallery.

So this script writes one small, deterministic micrograph of a pattern whose
answer is known, and a JSON sidecar carrying that answer. The browser test opens
the PNG exactly as a user would, indexes it, and compares the result with the
sidecar; the unit test regenerates the PNG and compares it byte for byte with
the tracked one, so the fixture cannot drift away from the simulation that
produced it.

Why zirconium down the basal axis
---------------------------------
It is the pattern an hcp session starts from and the one a microscopist
recognises on sight: six-fold, with the six prism {10-10} reflections innermost
and the six {11-20} at sqrt(3) times their radius — a ratio fixed by the
hexagonal net rather than by the material. A plate that indexes to anything else
is wrong in a way a reader of the test can see.

The pattern is rolled twelve degrees about the beam on purpose. A real plate is
never aligned with the detector axes, and a workflow tested only at zero roll has
not been tested.

Run manually with::

    python scripts/generate_tem_test_pattern.py

It writes ``fixtures/tem/zr_hcp_basal_saed.png`` and its ``.json`` sidecar, and
prints what it wrote.
"""

from __future__ import annotations

import json
import math
import struct
import sys
import zlib
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

#: Where the fixture lands. Named from what it is, so a second pattern of a
#: different phase or zone sits beside it without renaming this one.
OUTPUT_DIR = REPO_ROOT / "fixtures" / "tem"
STEM = "zr_hcp_basal_saed"

#: The exposure. A 400 mm camera length at 200 kV on a 24 micron detector pitch,
#: which is an ordinary setting rather than a convenient one.
PHASE_ID = "zr_hcp"
ZONE_AXIS = (0, 0, 1)
CAMERA_LENGTH_MM = 400.0
BEAM_ENERGY_KEV = 200.0
DETECTOR_PX = 512
PIXEL_SIZE_MM = 0.048
IN_PLANE_ROTATION_DEG = 12.0

#: Rendering. The background is not black and the spots are not white discs:
#: a plate has a fog level and the spots have wings, and a fixture that is
#: cleaner than any real exposure would let a centroiding bug pass.
BACKGROUND_LEVEL = 14
FOG_AMPLITUDE = 6.0
BEAM_RADIUS_PX = 7.0
SPOT_RADIUS_PX = 3.4


def _write_png(path: Path, image: np.ndarray) -> None:
    """Write an 8-bit greyscale PNG, with no dependency beyond the standard library.

    Pillow is not a dependency of this project and adding one to write a
    twelve-kilobyte fixture would be a poor trade. The format needed here is the
    simplest PNG there is: colour type 0, bit depth 8, one filter byte of zero
    per row.
    """

    if image.dtype != np.uint8 or image.ndim != 2:
        raise ValueError("_write_png takes a 2-D uint8 array.")
    height, width = image.shape

    def chunk(kind: bytes, payload: bytes) -> bytes:
        return (
            struct.pack(">I", len(payload))
            + kind
            + payload
            + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF)
        )

    header = struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0)
    # One zero filter byte in front of every scanline: filter type 0, "None".
    raw = b"".join(b"\x00" + image[row].tobytes() for row in range(height))
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", header)
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )


def _render(image_model: Any) -> np.ndarray:
    """Paint the simulated spots onto a raster, as an exposure rather than a plot.

    Each reflection is a Gaussian of a width set by the illumination rather than
    by its own strength, with the intensity carried by amplitude — which is what
    a plate records. The transmitted beam is wider and saturated, because it is
    orders of magnitude brighter than any reflection and is what a beam stop
    exists for.
    """

    width = int(image_model.raster.width_px)
    height = int(image_model.raster.height_px)
    columns, rows = np.meshgrid(np.arange(width, dtype=float), np.arange(height, dtype=float))

    # A slow, smooth background gradient: a real plate is not uniformly dark, and
    # a centroid taken over a window must survive that.
    canvas = np.full((height, width), float(BACKGROUND_LEVEL), dtype=float)
    canvas += FOG_AMPLITUDE * np.exp(
        -(((columns - width * 0.35) ** 2 + (rows - height * 0.6) ** 2) / (2.0 * (width * 0.7) ** 2))
    )

    centre_x, centre_y = image_model.centre_px
    canvas += 255.0 * np.exp(
        -(((columns - centre_x) ** 2 + (rows - centre_y) ** 2) / (2.0 * BEAM_RADIUS_PX**2))
    )

    for spot in image_model.spots:
        x, y = (float(value) for value in spot.position_px)
        amplitude = 60.0 + 175.0 * float(spot.relative_intensity)
        canvas += amplitude * np.exp(
            -(((columns - x) ** 2 + (rows - y) ** 2) / (2.0 * SPOT_RADIUS_PX**2))
        )

    return np.clip(np.rint(canvas), 0, 255).astype(np.uint8)


#: Significant digits every float in the sidecar is written to.
#:
#: The sidecar is a tracked baseline compared against a fresh generation, so it
#: has to be reproducible on any machine, and a float64 straight from the
#: simulation is not: the last one or two digits of a spot intensity differ
#: between BLAS builds, which failed the comparison on Linux while passing on
#: Windows. Twelve digits is far more precision than a diffraction intensity
#: carries and enough slack to absorb that.
_SIDECAR_SIGNIFICANT_DIGITS = 12


def _stable(value: float) -> float:
    """A float rounded to `_SIDECAR_SIGNIFICANT_DIGITS` significant digits."""

    number = float(value)
    if number == 0.0 or not math.isfinite(number):
        return number
    exponent = math.floor(math.log10(abs(number)))
    return round(number, _SIDECAR_SIGNIFICANT_DIGITS - 1 - exponent)


def build() -> tuple[np.ndarray, dict[str, Any]]:
    """Simulate the pattern and render it. Returns the raster and its answer.

    Every float in the returned sidecar is rounded to
    `_SIDECAR_SIGNIFICANT_DIGITS`; the raster is rendered from the unrounded
    simulation.
    """

    from pytex.app.phases import builtin_phase
    from pytex.core.lattice import ZoneAxis
    from pytex.diffraction.kinematic import electron_wavelength_angstrom
    from pytex.tem.synthetic import DetectorRaster, synthesize_saed_image

    spec = builtin_phase(PHASE_ID)
    phase = spec.to_phase()
    camera_constant = float(CAMERA_LENGTH_MM * electron_wavelength_angstrom(BEAM_ENERGY_KEV))
    image_model = synthesize_saed_image(
        phase,
        ZoneAxis(indices=np.asarray(ZONE_AXIS, dtype=int), phase=phase),
        camera_constant_mm_angstrom=camera_constant,
        raster=DetectorRaster(
            width_px=DETECTOR_PX,
            height_px=DETECTOR_PX,
            pixel_size_mm=PIXEL_SIZE_MM,
        ),
        in_plane_rotation_deg=IN_PLANE_ROTATION_DEG,
    )
    if not image_model.spots:
        raise RuntimeError("the chosen exposure puts every reflection off the plate")

    truth: dict[str, Any] = {
        "schema": "pytex.tem_test_pattern/1",
        "description": (
            "A simulated selected-area pattern of alpha zirconium with the beam along [0001], "
            "written as an 8-bit greyscale PNG so that opening a pattern from disk can be "
            "tested the way a user does it. Generated by scripts/generate_tem_test_pattern.py; "
            "regenerate rather than edit."
        ),
        "image": f"{STEM}.png",
        "phase": {"builtin": PHASE_ID},
        "phase_name": spec.name,
        "zone_axis": list(ZONE_AXIS),
        "beam_energy_kev": BEAM_ENERGY_KEV,
        "camera_length_mm": CAMERA_LENGTH_MM,
        "camera_constant_mm_angstrom": _stable(camera_constant),
        "pixel_size_mm": PIXEL_SIZE_MM,
        "in_plane_rotation_deg": IN_PLANE_ROTATION_DEG,
        "width_px": DETECTOR_PX,
        "height_px": DETECTOR_PX,
        "centre_px": [_stable(value) for value in image_model.centre_px],
        "crystal_to_pattern": [
            _stable(value) for value in image_model.crystal_to_pattern().reshape(-1)
        ],
        "spots": [
            {
                "hkl": [int(value) for value in spot.miller_indices],
                "x": _stable(spot.position_px[0]),
                "y": _stable(spot.position_px[1]),
                "d_angstrom": _stable(spot.d_spacing_angstrom),
                "g_inv_angstrom": _stable(spot.g_inv_angstrom),
                "intensity": _stable(spot.relative_intensity),
            }
            for spot in image_model.spots
        ],
        # Three strong, mutually non-collinear reflections: what a test (or a
        # user) should click to seed the indexing.
        "seed_spots": [
            {"x": _stable(spot.position_px[0]), "y": _stable(spot.position_px[1])}
            for spot in image_model.independent_seed_spots(3)
        ],
        "describe": image_model.describe(),
    }
    return _render(image_model), truth


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    raster, truth = build()
    png_path = OUTPUT_DIR / f"{STEM}.png"
    json_path = OUTPUT_DIR / f"{STEM}.json"
    _write_png(png_path, raster)
    json_path.write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {png_path.relative_to(REPO_ROOT)} ({png_path.stat().st_size} bytes)")
    print(f"wrote {json_path.relative_to(REPO_ROOT)} with {len(truth['spots'])} reflections")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
