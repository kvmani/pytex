"""XYZ input boundaries independent of platform filename limits and optional ASE."""

from pathlib import Path

import numpy as np
import pytest

from pytex.diffraction.hrem import AtomicSnapshot


def test_multiline_xyz_never_probes_the_filesystem(monkeypatch: pytest.MonkeyPatch) -> None:
    def unexpected_probe(*args: object, **kwargs: object) -> None:
        raise AssertionError("XYZ content must not be interpreted as a path")

    monkeypatch.setattr(Path, "read_text", unexpected_probe)
    monkeypatch.setattr(Path, "is_file", unexpected_probe)
    result = AtomicSnapshot.from_xyz("2\n\nC 0 1 2\nO 3 4 5\n", cell=np.eye(3) * 10)
    assert result.species == ("C", "O")
    np.testing.assert_array_equal(result.positions, [[0, 1, 2], [3, 4, 5]])
    assert result.label == "Imported XYZ snapshot"


@pytest.mark.parametrize("as_string", [False, True])
def test_xyz_file_input(tmp_path: Path, as_string: bool) -> None:
    path = tmp_path / "atoms.xyz"
    path.write_text("1\nSingle carbon\nC 0 0 0\n", encoding="utf-8")
    result = AtomicSnapshot.from_xyz(str(path) if as_string else path)
    assert result.natoms == 1
    assert result.label == "Single carbon"


@pytest.mark.parametrize("content", [
    "2\ntruncated\nC 0 0 0\n",
    "1\nextra frame\nC 0 0 0\n1\nsecond\nC 1 1 1\n",
    "0\ninvalid count\nC 0 0 0\n",
    "1\nnonfinite\nC nan 0 0\n",
    "1\nshort coordinates\nC 0 0\n",
])
def test_invalid_xyz_is_rejected(content: str) -> None:
    with pytest.raises(ValueError):
        AtomicSnapshot.from_xyz(content)


def test_missing_xyz_file_is_reported(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        AtomicSnapshot.from_xyz(tmp_path / "missing.xyz")
