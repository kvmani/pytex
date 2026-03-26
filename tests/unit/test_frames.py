from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose

from pytex.core import FrameDomain, FrameTransform, Handedness, ReferenceFrame


def make_frame(name: str, domain: FrameDomain) -> ReferenceFrame:
    return ReferenceFrame(
        name=name,
        domain=domain,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
    )


def test_identity_transform_round_trip() -> None:
    specimen = make_frame("specimen", FrameDomain.SPECIMEN)
    transform = FrameTransform.identity(specimen)
    vectors = np.array([[1.0, 2.0, 3.0], [-4.0, 5.0, 6.0]])
    assert_allclose(transform.apply_to_vectors(vectors), vectors)
    assert_allclose(transform.inverse().apply_to_vectors(vectors), vectors)


def test_transform_composition_matches_stepwise_application() -> None:
    crystal = make_frame("crystal", FrameDomain.CRYSTAL)
    specimen = make_frame("specimen", FrameDomain.SPECIMEN)
    lab = make_frame("lab", FrameDomain.LABORATORY)

    first = FrameTransform(
        source=crystal,
        target=specimen,
        rotation_matrix=np.eye(3),
        translation_vector=np.array([1.0, 0.0, 0.0]),
    )
    second = FrameTransform(
        source=specimen,
        target=lab,
        rotation_matrix=np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        translation_vector=np.array([0.0, 2.0, 0.0]),
    )
    chained = second.compose(first)
    vector = np.array([[1.0, 2.0, 3.0]])
    stepwise = second.apply_to_vectors(first.apply_to_vectors(vector))
    assert_allclose(chained.apply_to_vectors(vector), stepwise)
