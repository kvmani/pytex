"""Analysing a user's own EBSD scan, not only the practice datasets.

The claim under test is that an imported scan is not a second-class dataset: it
reaches the same segmentation, the same KAM and GROD, the same boundary network
and the same colourings as a constructed map, because it arrives as the same
`CrystalMap`. The scans here are small enough to have answers worked out by hand,
which is what lets the assertions be about numbers rather than about shapes.
"""

from __future__ import annotations

import base64

import numpy as np
import pytest

from pytex.app import REGISTRY
from pytex.app.errors import InvalidInputError
from pytex.app.services.ebsd import _EMPTY_CELL

# Four points, two rows, two columns. The first three share an orientation and
# the fourth is 60 degrees about [001] away from them, so a 5 degree grain
# threshold must find exactly two grains.
SQUARE_ANG = """\
# Phase 1
# MaterialName  	Nickel
# Formula     	Ni
# Symmetry              43
# LatticeConstants      3.520 3.520 3.520  90.000  90.000  90.000
#
# GRID: SqrGrid
# XSTEP: 1.000000
# YSTEP: 1.000000
# NCOLS_ODD: 2
# NCOLS_EVEN: 2
# NROWS: 2
#
   0.00000   0.00000   0.00000      0.00000      0.00000  60.0  0.950  0  1  0.500
   0.00000   0.00000   0.00000      1.00000      0.00000  55.0  0.900  0  1  0.600
   0.00000   0.00000   0.00000      0.00000      1.00000  50.0  0.850  0  1  0.700
   1.04720   0.00000   0.00000      1.00000      1.00000  45.0  0.800  0  1  0.800
"""

SQUARE_CTF = (
    "Channel Text File\n"
    "Prj\timport-test\n"
    "JobMode\tGrid\n"
    "XCells\t2\n"
    "YCells\t2\n"
    "XStep\t0.5000\n"
    "YStep\t0.5000\n"
    "AcqE1\t0\n"
    "Euler angles refer to Sample Coordinate system (CS0)!\tMag\t500\n"
    "Phases\t1\n"
    "3.524;3.524;3.524\t90;90;90\tNickel\t11\t225\t\tcomment\n"
    "Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\n"
    "1\t0.0\t0.0\t10\t0\t0.0\t0.0\t0.0\t0.40\t180\t200\n"
    "1\t0.5\t0.0\t10\t0\t0.0\t0.0\t0.0\t0.40\t180\t200\n"
    "1\t0.0\t0.5\t10\t0\t0.0\t0.0\t0.0\t0.40\t180\t200\n"
    "1\t0.5\t0.5\t11\t0\t60.0\t0.0\t0.0\t0.30\t190\t210\n"
)

HEX_ANG = "fixtures/ebsd/synthetic_hex_grid.ang"  # a staggered scan, read from the repository


def analyse(name: str, text: str, **overrides: object) -> dict:
    request: dict[str, object] = {"scan_file": {"name": name, "text": text}}
    request.update(overrides)
    return REGISTRY.call("ebsd.map", request)


def image_pixels(result: dict) -> np.ndarray:
    image = result["data"]["image"]
    raw = base64.b64decode(image["data"])
    return np.frombuffer(raw, dtype=np.uint8).reshape(image["height"], image["width"], 3)


def test_a_square_ang_scan_is_analysed_like_any_other_map() -> None:
    result = analyse("nickel.ang", SQUARE_ANG)

    assert result["data"]["dataset"]["id"] == "file:nickel.ang"
    assert result["title"].endswith("nickel.ang")
    # The construction: three points together, one 60 degrees away.
    assert result["data"]["grain_count"] == 2
    assert result["data"]["grid_shape"] == [2, 2]
    assert image_pixels(result).shape == (2, 2, 3)
    sizes = sorted(row["size"] for row in result["data"]["grains"])
    assert sizes == [1, 3]


def test_a_ctf_scan_reaches_the_same_place() -> None:
    result = analyse("nickel.ctf", SQUARE_CTF)
    assert result["data"]["grain_count"] == 2
    assert result["data"]["grid_shape"] == [2, 2]
    assert sorted(row["size"] for row in result["data"]["grains"]) == [1, 3]
    # The .ctf step is 0.5 um, and every reported length follows from it.
    assert result["data"]["step_um"] == pytest.approx(0.5)
    assert result["data"]["extent_um"] == pytest.approx([0.0, 0.0, 0.5, 0.5])


def test_the_boundary_between_the_two_grains_carries_its_own_misorientation() -> None:
    """60 degrees about [001] in a cubic phase is 90 degrees reduced.

    The disorientation of a 60 degree rotation about a four-fold axis is 30
    degrees, because the symmetry offers a shorter equivalent. Asserting the
    reduced angle rather than the one written in the file is the point: the
    imported scan goes through the same symmetry-aware machinery as everything
    else.
    """

    result = analyse("nickel.ang", SQUARE_ANG)
    summary = {row["character"]: row for row in result["data"]["boundary_summary"]}
    total = sum(row["count"] for row in summary.values())
    assert total > 0
    angles = [row["mean_misorientation_deg"] for row in summary.values() if row["count"]]
    # 1e-3, not machine epsilon: the file stores the Euler angle as 1.04720
    # radians, which is 60.00008 degrees, so the answer cannot be exact and
    # pretending otherwise would be testing the rounding of the fixture.
    assert angles
    assert all(angle == pytest.approx(30.0, abs=1e-3) for angle in angles)


def test_a_staggered_scan_keeps_its_offset_rows_and_invents_nothing() -> None:
    """A hexagonal scan has no rectangular shape, and is drawn as what it is.

    EDAX writes hexagonal scans by default, so refusing them would refuse most
    ``.ang`` files. They are placed on a half-step raster, which is the pitch the
    stagger lives on: every measurement gets its own cell, the cells between them
    stay empty, and nothing is interpolated into the gaps.
    """

    from pathlib import Path

    text = Path(HEX_ANG).read_text(encoding="utf-8")
    result = analyse("hex.ang", text)

    pixels = image_pixels(result)
    # Three rows of a 3/2/3 staggered scan on a half-step raster: 5 columns.
    assert pixels.shape == (3, 5, 3)
    filled = np.any(pixels != _EMPTY_CELL, axis=-1)
    assert int(filled.sum()) == 8, "every measurement must land on exactly one cell"
    # The middle row is the short one, and its two points sit between the others.
    assert filled[0].tolist() == [True, False, True, False, True]
    assert filled[1].tolist() == [False, True, False, True, False]
    assert "staggered scan" in result["summary"]


def test_an_imported_map_refuses_to_claim_a_known_answer() -> None:
    """The practice datasets are checkable because they were constructed.

    A measurement is not, and presenting one as though it carried the same
    guarantee would be the most misleading thing this panel could do.
    """

    dataset = analyse("nickel.ang", SQUARE_ANG)["data"]["dataset"]
    assert "None" in dataset["known_answer"]
    assert "measurement, not a construction" in dataset["known_answer"]


def test_the_quality_channels_come_from_the_file() -> None:
    result = analyse("nickel.ang", SQUARE_ANG, colouring="confidence_index")
    scale = result["data"]["colour_scale"]
    assert scale["label"] == "Confidence index"
    # The four CI values in the file are 0.95, 0.90, 0.85 and 0.80.
    assert scale["minimum"] == pytest.approx(0.80)
    assert scale["maximum"] == pytest.approx(0.95)


def test_an_open_file_wins_over_the_practice_dataset() -> None:
    """Someone who has opened their data is looking at their data."""

    result = analyse("nickel.ang", SQUARE_ANG, dataset="equiaxed_polycrystal")
    assert result["data"]["dataset"]["id"] == "file:nickel.ang"
    assert result["data"]["grain_count"] == 2


@pytest.mark.parametrize(
    ("payload", "fragment"),
    [
        ({"name": "map.txt", "text": "1 2 3"}, "not a file this reads"),
        ({"name": "map.ang", "text": "   "}, "empty"),
        ({"name": "map.ang", "text": "not an ang file at all\n"}, "could not be read"),
    ],
)
def test_an_unreadable_file_is_refused_beside_its_own_control(payload: dict, fragment: str) -> None:
    with pytest.raises(InvalidInputError) as raised:
        REGISTRY.call("ebsd.map", {"scan_file": payload})
    assert fragment in str(raised.value)
    assert raised.value.details["field"] == "scan_file"


def test_nothing_is_left_on_disk_afterwards() -> None:
    """The temporary file exists only while the reader has it open."""

    import tempfile
    from pathlib import Path

    before = set(Path(tempfile.gettempdir()).glob("*.ang"))
    analyse("nickel.ang", SQUARE_ANG)
    after = set(Path(tempfile.gettempdir()).glob("*.ang"))
    assert after <= before
