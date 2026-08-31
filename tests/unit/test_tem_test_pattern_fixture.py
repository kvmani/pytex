"""The tracked TEM test pattern: still the pattern its generator produces.

A raster fixture is a pinned baseline, and a pinned baseline that nothing checks
is a file that quietly stops meaning what its name says. Three things are
asserted here, in increasing order of how badly they would mislead if wrong:

1. **The PNG is what the generator writes today.** Byte for byte, from a
   regeneration into a temporary directory. The simulation is deterministic —
   no jitter, no random seed in play — so a difference is a real change in the
   pattern rather than in the weather.
2. **The sidecar describes that PNG.** Size, phase, zone axis and calibration
   agree with the file, so a test that reads the answer from the JSON is reading
   the answer to the picture it opened.
3. **The rendered image really has spots where the answer says.** The brightest
   pixel in a window about each stated reflection is at that reflection, which is
   what makes clicking on the image and getting the indexed answer possible at
   all.

The pattern itself is checked against crystallography rather than against stored
numbers: down [0001] the hexagonal net puts the {11-20} ring at exactly sqrt(3)
times the radius of the {10-10} ring, whatever the axial ratio is.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from scripts.generate_tem_test_pattern import (
    OUTPUT_DIR,
    STEM,
    _write_png,
    build,
)

PNG_PATH = OUTPUT_DIR / f"{STEM}.png"
JSON_PATH = OUTPUT_DIR / f"{STEM}.json"


@pytest.fixture(scope="module")
def truth() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def test_the_tracked_files_exist_and_are_small() -> None:
    """A fixture in git forever should be kilobytes, not megabytes."""

    assert PNG_PATH.is_file(), "run scripts/generate_tem_test_pattern.py"
    assert JSON_PATH.is_file()
    assert PNG_PATH.stat().st_size < 200_000


def test_the_png_is_byte_identical_to_a_fresh_generation(tmp_path: Path) -> None:
    raster, _ = build()
    regenerated = tmp_path / "regenerated.png"
    _write_png(regenerated, raster)
    assert regenerated.read_bytes() == PNG_PATH.read_bytes(), (
        "the tracked pattern is no longer what scripts/generate_tem_test_pattern.py "
        "produces; regenerate it in the same commit as the change that moved it"
    )


def test_the_sidecar_answers_for_the_image_beside_it(truth: dict) -> None:
    _, regenerated = build()
    assert truth == regenerated
    assert truth["image"] == PNG_PATH.name
    assert truth["phase"] == {"builtin": "zr_hcp"}
    assert truth["zone_axis"] == [0, 0, 1]
    assert truth["camera_constant_mm_angstrom"] > 0.0
    assert len(truth["seed_spots"]) == 3


def test_the_basal_net_places_the_two_inner_rings_a_root_three_apart(truth: dict) -> None:
    radii = sorted({round(spot["g_inv_angstrom"], 6) for spot in truth["spots"]})
    assert len(radii) >= 2
    # A property of the hexagonal net, not of zirconium: |g_11-20| / |g_10-10| =
    # sqrt(3) for any c/a. The tolerance is a part in a hundred thousand because
    # the pinned cell is a measured one, whose a and b agree to the digits the
    # CIF states.
    assert radii[1] / radii[0] == pytest.approx(math.sqrt(3.0), rel=1e-5)


def test_every_stated_reflection_is_actually_bright_in_the_picture(truth: dict) -> None:
    """The answer and the pixels must be about the same pattern."""

    raster, _ = build()
    assert raster.shape == (truth["height_px"], truth["width_px"])

    centre = np.asarray(truth["centre_px"], dtype=float)
    background = float(np.median(raster))
    half = 3
    checked = 0
    for spot in truth["spots"]:
        column = round(spot["x"])
        row = round(spot["y"])
        if not (half <= column < raster.shape[1] - half and half <= row < raster.shape[0] - half):
            continue
        # Skip the few reflections close enough to the beam that its own glow
        # dominates the window: there the brightest pixel is honestly the beam.
        if math.hypot(spot["x"] - centre[0], spot["y"] - centre[1]) < 20.0:
            continue
        window = raster[row - half : row + half + 1, column - half : column + half + 1]
        peak = np.unravel_index(int(np.argmax(window)), window.shape)
        assert abs(int(peak[0]) - half) <= 1 and abs(int(peak[1]) - half) <= 1
        assert float(window.max()) > background + 30.0
        checked += 1
    assert checked > 20, "the fixture should carry a plateful of checkable reflections"
