"""Computed images as files, at the resolution and precision they were computed at.

A simulated micrograph is an array of numbers. Drawing it into a figure and
rasterizing the figure resamples it and quantizes it to eight bits of a colour
map; neither is what a user who wants to measure the image, or to put it next
to an experimental one, should get. The workbench therefore offers every
computed image in three forms, each without resampling:

- **PNG**, one pixel per computed pixel, eight bits through the colour map the
  panel shows: the picture as seen, for slides and papers.
- **TIFF**, one pixel per computed pixel, 32-bit floating point: the numbers
  themselves, which ImageJ/Fiji, DigitalMicrograph, Python and MATLAB all read.
- **ZIP**, for a series: every image in both forms, a native-resolution montage,
  and a table of what each image is.

Orientation is fixed once, here. PyTex image arrays put row 0 at ``y = 0``, the
bottom of the specimen, with ``y`` increasing upwards; image files put their
first row at the top. Every writer below therefore flips the rows, so a file
opens the same way up as the panel draws it.
"""

from __future__ import annotations

import base64
import io
import zipfile
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

#: Media types of the raster forms, for data URLs and downloads.
PNG_TYPE = "image/png"
TIFF_TYPE = "image/tiff"
ZIP_TYPE = "application/zip"


def _top_down(array: np.ndarray) -> np.ndarray:
    values = np.asarray(array)
    if values.ndim not in (2, 3):
        raise ValueError(f"An image needs two dimensions (or three for RGB), got {values.shape}.")
    return np.ascontiguousarray(values[::-1])


def data_url(payload: bytes, media_type: str) -> str:
    """``payload`` as a base64 data URL of the given media type."""
    return f"data:{media_type};base64," + base64.b64encode(payload).decode("ascii")


def float_tiff_bytes(array: np.ndarray) -> bytes:
    """A 2-D array as a single-channel 32-bit floating-point TIFF, top row first.

    Values are written as computed (cast to float32); no scaling or colour map
    is applied, so the file is the data.
    """

    from PIL import Image

    values = _top_down(np.asarray(array, dtype=np.float32))
    if values.ndim != 2:
        raise ValueError("A float TIFF holds one channel.")
    buffer = io.BytesIO()
    Image.fromarray(values).save(buffer, format="TIFF")
    return buffer.getvalue()


def gray_png_bytes(array: np.ndarray, *, colormap: str = "gray") -> bytes:
    """A 2-D array as an 8-bit PNG through ``colormap``, scaled to its own range, top row first."""

    import matplotlib

    values = np.asarray(array, dtype=np.float64)
    if values.ndim != 2:
        raise ValueError("A colour-mapped PNG is drawn from one channel.")
    low, high = float(np.min(values)), float(np.max(values))
    scaled = (values - low) / (high - low) if high > low else np.zeros_like(values)
    rgba = matplotlib.colormaps[colormap](_top_down(scaled), bytes=True)
    return rgb_png_bytes(rgba[..., :3], top_down=True)


def rgb_png_bytes(pixels: np.ndarray, *, top_down: bool = False) -> bytes:
    """An ``(ny, nx, 3)`` uint8 array as a PNG; flipped to top-down unless it already is."""

    from PIL import Image

    values = np.asarray(pixels, dtype=np.uint8)
    if values.ndim != 3 or values.shape[2] != 3:
        raise ValueError(f"RGB pixels must have shape (ny, nx, 3), got {values.shape}.")
    buffer = io.BytesIO()
    Image.fromarray(values if top_down else _top_down(values)).save(buffer, format="PNG")
    return buffer.getvalue()


def montage(
    tiles: np.ndarray, *, gap_px: int = 4, colormap: str = "gray"
) -> np.ndarray:
    """Tiles ``(rows, columns, ny, nx)`` as one RGB image at their native pixels.

    Each tile is scaled to its own range, as the workbench draws a tableau, and
    tiles are separated by ``gap_px`` dark pixels. Row 0 of the tableau is at the
    top of the result; within a tile ``y`` still increases upwards, so the
    result is returned top-down and must be written with ``top_down=True``.
    """

    import matplotlib

    stack = np.asarray(tiles, dtype=np.float64)
    rows, columns, ny, nx = stack.shape
    height = rows * ny + (rows - 1) * gap_px
    width = columns * nx + (columns - 1) * gap_px
    canvas = np.full((height, width, 3), 24, dtype=np.uint8)
    cmap = matplotlib.colormaps[colormap]
    for i in range(rows):
        for j in range(columns):
            tile = stack[i, j]
            low, high = float(tile.min()), float(tile.max())
            scaled = (tile - low) / (high - low) if high > low else np.zeros_like(tile)
            rgb = cmap(_top_down(scaled), bytes=True)[..., :3]
            y0, x0 = i * (ny + gap_px), j * (nx + gap_px)
            canvas[y0 : y0 + ny, x0 : x0 + nx] = rgb
    return canvas


def zip_bytes(files: Mapping[str, bytes]) -> bytes:
    """Files as a deflated ZIP archive, in the order given."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, payload in files.items():
            archive.writestr(name, payload)
    return buffer.getvalue()


def raster(
    label: str, payload: bytes, media_type: str, filename: str, shape: Sequence[int]
) -> dict[str, Any]:
    """One downloadable raster in the form the workbench's Save menu lists.

    ``shape`` is ``(height, width)`` in pixels, for the menu's size note.
    """

    return {
        "label": label,
        "data": data_url(payload, media_type),
        "filename": filename,
        "width": int(shape[1]),
        "height": int(shape[0]),
    }


def tiff_raster(array: np.ndarray, *, name: str, stem: str, quantity: str) -> dict[str, Any]:
    """The 32-bit float TIFF download of one computed image.

    ``name`` says what the image is (``"micrograph"``), ``stem`` is the file-name
    stem and ``quantity`` what the numbers are (``"intensity"``), for the label.
    The panel offers the PNG itself from the image it already draws.
    """

    values = np.asarray(array)
    ny, nx = values.shape
    return raster(
        f"Download TIFF ({name}, 32-bit {quantity}, {nx} × {ny} px)",
        float_tiff_bytes(values),
        TIFF_TYPE,
        f"{stem}-{nx}x{ny}px-float32.tif",
        values.shape,
    )
