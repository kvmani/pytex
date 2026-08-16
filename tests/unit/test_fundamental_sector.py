from __future__ import annotations

import matplotlib
import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.core import FrameDomain, Handedness, ReferenceFrame, SymmetrySpec
from pytex.plotting.ipf import IPFColorKey, plot_ipf_key

matplotlib.use("Agg")

PROPER_GROUPS = ("1", "2", "222", "4", "422", "3", "32", "6", "622", "23", "432")


def make_crystal_frame() -> ReferenceFrame:
    return ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
    )


def random_unit_vectors(count: int, *, seed: int = 7) -> np.ndarray:
    rng = np.random.default_rng(seed)
    vectors = rng.normal(size=(count, 3))
    return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)


@pytest.mark.parametrize("proper_group", PROPER_GROUPS)
def test_sector_area_fraction_matches_group_order(proper_group: str) -> None:
    spec = SymmetrySpec.from_point_group(proper_group)
    sector = spec.fundamental_sector(antipodal=True)
    samples = random_unit_vectors(60000)
    fraction = float(np.count_nonzero(sector.contains(samples))) / samples.shape[0]
    expected = 1.0 / (2.0 * spec.order)
    tolerance = 5.0 * float(np.sqrt(expected * (1.0 - expected) / samples.shape[0]))
    assert abs(fraction - expected) < tolerance


@pytest.mark.parametrize("proper_group", PROPER_GROUPS)
def test_reduction_lands_inside_the_exact_sector(proper_group: str) -> None:
    spec = SymmetrySpec.from_point_group(proper_group)
    sector = spec.fundamental_sector(antipodal=True)
    samples = random_unit_vectors(200, seed=11)
    reduced = spec.reduce_vectors_to_fundamental_sector(samples, antipodal=True)
    assert np.all(sector.contains(reduced, atol=1e-6))
    # Reduction must be idempotent on the exact sector.
    twice = spec.reduce_vectors_to_fundamental_sector(reduced, antipodal=True)
    assert_allclose(np.asarray(twice), np.asarray(reduced), atol=1e-10)


@pytest.mark.parametrize("proper_group", PROPER_GROUPS)
def test_boundary_trace_lies_on_sphere_and_sector(proper_group: str) -> None:
    sector = SymmetrySpec.from_point_group(proper_group).fundamental_sector(antipodal=True)
    trace = sector.boundary_trace(samples_per_edge=32)
    assert trace.shape[0] >= 32
    assert_allclose(np.linalg.norm(trace, axis=1), 1.0, atol=1e-10)
    assert np.all(sector.contains(trace, atol=1e-8))
    center = sector.center()
    assert sector.contains(center, atol=1e-8)
    assert_allclose(np.linalg.norm(center), 1.0, atol=1e-12)


def test_cyclic_groups_get_larger_sectors_than_dihedral_groups() -> None:
    samples = random_unit_vectors(30000, seed=3)

    def fraction(proper_group: str) -> float:
        sector = SymmetrySpec.from_point_group(proper_group).fundamental_sector(antipodal=True)
        return float(np.count_nonzero(sector.contains(samples))) / samples.shape[0]

    assert fraction("2") == pytest.approx(2.0 * fraction("222"), rel=0.1)
    assert fraction("3") == pytest.approx(2.0 * fraction("32"), rel=0.1)
    assert fraction("6") == pytest.approx(2.0 * fraction("622"), rel=0.1)
    assert fraction("23") == pytest.approx(2.0 * fraction("432"), rel=0.1)


def test_tetrahedral_sector_differs_from_octahedral_sector() -> None:
    tetrahedral = SymmetrySpec.from_point_group("23").fundamental_sector(antipodal=True)
    octahedral = SymmetrySpec.from_point_group("432").fundamental_sector(antipodal=True)
    assert not np.allclose(tetrahedral.vertices, octahedral.vertices)
    # The 100 axis is a corner of the tetrahedral sector but the octahedral
    # sector corner along that edge is 101.
    assert_allclose(tetrahedral.vertices[1], [1.0, 0.0, 0.0], atol=1e-12)
    assert_allclose(octahedral.vertices[1], np.array([1.0, 0.0, 1.0]) / np.sqrt(2.0), atol=1e-12)


def test_hemisphere_sector_boundary_is_the_equator() -> None:
    sector = SymmetrySpec.from_point_group("-1").fundamental_sector(antipodal=True)
    assert sector.vertices.shape[0] == 0
    trace = sector.boundary_trace(samples_per_edge=16)
    assert_allclose(trace[:, 2], 0.0, atol=1e-12)


def test_monoclinic_wedge_covers_half_hemisphere() -> None:
    sector = SymmetrySpec.from_point_group("2/m").fundamental_sector(antipodal=True)
    assert sector.contains([0.0, 1.0, 0.0])
    assert sector.contains([-0.5, 0.5, 0.5])
    assert not sector.contains([0.0, -1.0, 0.0])
    assert not sector.contains([0.0, 0.5, -0.5])


BOUNDARY_DIRECTIONS = np.array(
    [
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [1.0, 1.0, 2.0],
        [-1.0, 2.0, 3.0],
        [0.0, 0.0, -1.0],
        [-1.0, -1.0, -1.0],
    ]
)


def reduce_by_orbit_search(
    spec: SymmetrySpec,
    vector: np.ndarray,
    *,
    antipodal: bool,
) -> np.ndarray:
    """A deliberately naive reduction, used as the reference answer.

    Enumerates the orbit one operator at a time and picks the in-sector member
    with the largest ``(z, x, y)``, falling back to the canonical representative
    when the sector admits none. This is the definition the vectorized
    :meth:`SymmetrySpec.reduce_vectors_to_fundamental_sector` must reproduce, so
    it is written independently of that implementation.
    """

    sector = spec.fundamental_sector(antipodal=antipodal)
    candidates = spec.equivalent_vectors(vector, antipodal=antipodal)
    inside = [
        candidate
        for candidate in candidates
        if bool(np.all(sector.edge_normals @ candidate >= -1e-8))
    ]

    def order_key(item: np.ndarray) -> tuple[float, float, float]:
        rounded = np.round(item, decimals=12)
        return (float(rounded[2]), float(rounded[0]), float(rounded[1]))

    if inside:
        return np.asarray(max(inside, key=order_key), dtype=float)
    return np.asarray(spec.canonicalize_vector(vector, antipodal=antipodal), dtype=float)


@pytest.mark.parametrize("proper_group", PROPER_GROUPS)
@pytest.mark.parametrize("antipodal", [True, False])
def test_batch_reduction_matches_an_independent_orbit_search(
    proper_group: str,
    antipodal: bool,
) -> None:
    spec = SymmetrySpec.from_point_group(proper_group)
    samples = np.concatenate(
        [
            random_unit_vectors(120, seed=23),
            BOUNDARY_DIRECTIONS / np.linalg.norm(BOUNDARY_DIRECTIONS, axis=1, keepdims=True),
        ]
    )
    reduced = np.asarray(spec.reduce_vectors_to_fundamental_sector(samples, antipodal=antipodal))
    expected = np.stack(
        [reduce_by_orbit_search(spec, vector, antipodal=antipodal) for vector in samples]
    )
    assert_allclose(reduced, expected, atol=1e-12)
    # The single-vector entry point must agree with its own batch form.
    for index in (0, 5, len(samples) - 1):
        single = spec.reduce_vector_to_fundamental_sector(samples[index], antipodal=antipodal)
        assert_allclose(single, reduced[index], atol=1e-12)


def test_batch_reduction_is_independent_of_the_chunk_boundary() -> None:
    from pytex.core import symmetry as symmetry_module

    spec = SymmetrySpec.from_point_group("432")
    samples = random_unit_vectors(symmetry_module._SECTOR_REDUCTION_CHUNK + 37, seed=5)
    whole = np.asarray(spec.reduce_vectors_to_fundamental_sector(samples))
    split = np.concatenate(
        [
            np.asarray(spec.reduce_vectors_to_fundamental_sector(samples[:1000])),
            np.asarray(spec.reduce_vectors_to_fundamental_sector(samples[1000:])),
        ]
    )
    assert whole.shape == samples.shape
    assert_allclose(whole, split, atol=1e-14)


def test_batch_reduction_accepts_an_empty_set() -> None:
    spec = SymmetrySpec.from_point_group("m-3m")
    reduced = np.asarray(spec.reduce_vectors_to_fundamental_sector(np.zeros((0, 3))))
    assert reduced.shape == (0, 3)


def test_cubic_reduction_sends_symmetry_axes_to_the_sector_corners() -> None:
    spec = SymmetrySpec.from_point_group("m-3m")
    axes = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 1.0], [-1.0, 1.0, -1.0]])
    reduced = np.asarray(spec.reduce_vectors_to_fundamental_sector(axes))
    expected = np.array(
        [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0] / np.sqrt(2.0),
            [1.0, 1.0, 1.0] / np.sqrt(3.0),
        ]
    )
    assert_allclose(reduced, expected, atol=1e-12)


def make_cubic_key() -> IPFColorKey:
    symmetry = SymmetrySpec.from_point_group("m-3m", reference_frame=make_crystal_frame())
    return IPFColorKey(crystal_symmetry=symmetry, specimen_direction="ND")


def test_legend_mesh_produces_valid_colors_and_points() -> None:
    key = make_cubic_key()
    points, colors = key.legend_mesh(resolution_deg=3.0)
    assert points.shape[0] == colors.shape[0]
    assert points.shape[1] == 2
    assert np.all(np.isfinite(points))
    assert np.all((colors >= 0.0) & (colors <= 1.0))
    boundary = key.boundary_points_2d()
    assert boundary.shape[1] == 2
    assert np.all(np.isfinite(boundary))
    with pytest.raises(ValueError):
        key.legend_mesh(resolution_deg=0.0)


LAUE_CLASSES = ("-1", "2/m", "mmm", "4/m", "4/mmm", "-3", "-3m", "6/m", "6/mmm", "m-3", "m-3m")


def test_legend_mesh_works_for_every_laue_class() -> None:
    frame = make_crystal_frame()
    for point_group in LAUE_CLASSES:
        symmetry = SymmetrySpec.from_point_group(point_group, reference_frame=frame)
        key = IPFColorKey(crystal_symmetry=symmetry, specimen_direction="ND")
        points, colors = key.legend_mesh(resolution_deg=6.0)
        assert points.shape[0] > 0
        assert np.all((colors >= 0.0) & (colors <= 1.0))


def test_plot_ipf_key_renders_scatter_and_boundary() -> None:
    import matplotlib.pyplot as plt

    key = make_cubic_key()
    fig, axis = plot_ipf_key(key, resolution_deg=4.0)
    try:
        assert len(axis.collections) == 1
        assert len(axis.lines) == 1
        assert axis.get_title() == "IPF key m-3m"
    finally:
        plt.close(fig)
