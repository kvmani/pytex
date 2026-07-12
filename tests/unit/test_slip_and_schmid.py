from __future__ import annotations

import math

import numpy as np
import pytest

from pytex import (
    CrystalMap,
    FrameDomain,
    Handedness,
    Lattice,
    Orientation,
    OrientationSet,
    Phase,
    ReferenceFrame,
    Rotation,
    SlipSystem,
    SymmetrySpec,
    bcc_110_slip,
    fcc_octahedral_slip,
)
from pytex.core.miller import MillerDirection, MillerPlane


def _cubic_phase(name: str = "fcc") -> tuple[Phase, ReferenceFrame, ReferenceFrame]:
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=crystal)
    lattice = Lattice(3.6, 3.6, 3.6, 90.0, 90.0, 90.0, crystal_frame=crystal)
    phase = Phase(name=name, lattice=lattice, symmetry=symmetry, crystal_frame=crystal)
    return phase, crystal, specimen


def test_fcc_and_bcc_families_have_twelve_systems() -> None:
    phase, _, _ = _cubic_phase()
    assert fcc_octahedral_slip(phase).count == 12
    assert bcc_110_slip(phase).count == 12


def test_slip_system_requires_direction_in_plane() -> None:
    phase, _, _ = _cubic_phase()
    with pytest.raises(ValueError, match="lie in the slip plane"):
        SlipSystem(
            plane=MillerPlane.from_hkl((1, 1, 1), phase=phase),
            direction=MillerDirection.from_uvw((1, 0, 0), phase=phase),
        )


def test_max_schmid_factor_matches_analytic_001_value() -> None:
    phase, crystal, specimen = _cubic_phase()
    symmetry = phase.symmetry
    orientation = Orientation(
        Rotation.identity(),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        phase=phase,
    )
    family = fcc_octahedral_slip(phase)
    # The classic maximum Schmid factor for <110>{111} under [001] is sqrt(6)/6.
    expected = math.sqrt(6.0) / 6.0
    assert family.max_schmid_factor(orientation, [0.0, 0.0, 1.0]) == pytest.approx(expected)
    factors = family.schmid_factors(orientation, [0.0, 0.0, 1.0])
    assert factors.shape == (12,)
    assert np.all(factors <= 0.5 + 1e-9)


def test_schmid_factors_are_orientation_batched() -> None:
    phase, crystal, specimen = _cubic_phase()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            ),
            Orientation(
                Rotation.from_bunge_euler(30.0, 20.0, 10.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            ),
        ]
    )
    family = fcc_octahedral_slip(phase)
    factors = family.schmid_factors(orientations, [0.0, 0.0, 1.0])
    assert factors.shape == (2, 12)
    maxima = family.max_schmid_factor(orientations, [0.0, 0.0, 1.0])
    assert maxima.shape == (2,)
    assert maxima[0] == pytest.approx(math.sqrt(6.0) / 6.0)


def test_crystal_map_schmid_factor_map_reshapes_to_grid() -> None:
    phase, crystal, specimen = _cubic_phase()
    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=phase.symmetry,
                phase=phase,
            )
            for _ in range(4)
        ]
    )
    crystal_map = CrystalMap(
        coordinates=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        orientations=orientations,
        map_frame=specimen,
        grid_shape=(2, 2),
        step_sizes=(1.0, 1.0),
    )
    family = fcc_octahedral_slip(phase)
    schmid = crystal_map.schmid_factor_map(family, "z")
    assert schmid.shape == (2, 2)
    assert np.allclose(schmid, math.sqrt(6.0) / 6.0)
