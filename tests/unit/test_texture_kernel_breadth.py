from __future__ import annotations

import numpy as np
import pytest

from pytex import AbelPoissonKernel, GaussianSO3Kernel
from pytex.texture import KernelSpec

KERNELS = [GaussianSO3Kernel, AbelPoissonKernel]


@pytest.mark.parametrize("kernel_class", KERNELS)
def test_kernels_are_normalized_and_satisfy_the_halfwidth_property(kernel_class) -> None:
    kernel = kernel_class(10.0)
    coefficients = kernel.chebyshev_coefficients(6)
    # A_0 = 1 exactly: the kernel integrates to one over SO(3).
    assert coefficients[0] == pytest.approx(1.0, abs=1e-12)
    peak = float(kernel.evaluate(np.zeros(1))[0])
    at_halfwidth = float(kernel.evaluate(np.array([np.deg2rad(10.0)]))[0])
    assert at_halfwidth / peak == pytest.approx(0.5, abs=1e-6)
    assert peak > 0.0


@pytest.mark.parametrize("kernel_class", KERNELS)
def test_kernel_coefficients_round_trip_through_quadrature(kernel_class) -> None:
    """Quadrature of the evaluated kernel recovers the closed-form spectrum."""

    kernel = kernel_class(15.0)
    omega = np.linspace(0.0, np.pi, 20001)
    psi = kernel.evaluate(omega)
    half = omega / 2.0
    spacing = float(omega[1] - omega[0])
    expected = kernel.chebyshev_coefficients(6)
    for order in range(7):
        integrand = psi * (2.0 / np.pi) * np.sin(half) * np.sin((2 * order + 1) * half)
        recovered = spacing * (
            0.5 * integrand[0] + integrand[1:-1].sum() + 0.5 * integrand[-1]
        )
        assert recovered == pytest.approx(expected[order], abs=2e-4)


@pytest.mark.parametrize("kernel_class", KERNELS)
def test_kernel_bandwidth_grows_as_halfwidth_shrinks(kernel_class) -> None:
    broad = kernel_class(20.0).bandwidth()
    sharp = kernel_class(5.0).bandwidth()
    assert sharp > broad > 0


def test_gaussian_decays_faster_than_abel_poisson_at_equal_halfwidth() -> None:
    orders = 24
    gaussian = GaussianSO3Kernel(10.0).chebyshev_coefficients(orders)
    abel = AbelPoissonKernel(10.0).chebyshev_coefficients(orders)
    # The heat-kernel spectrum decays super-geometrically; Abel-Poisson only
    # geometrically, so its tail dominates.
    assert gaussian[orders] < abel[orders]


def test_kernel_spec_routes_new_kernel_names() -> None:
    gaussian = KernelSpec(name="gaussian", halfwidth_deg=12.0).as_so3_kernel()
    abel = KernelSpec(name="abel_poisson", halfwidth_deg=12.0).as_so3_kernel()
    assert isinstance(gaussian, GaussianSO3Kernel)
    assert isinstance(abel, AbelPoissonKernel)
    values = KernelSpec(name="gaussian", halfwidth_deg=12.0).evaluate(
        np.array([0.0, np.deg2rad(12.0)])
    )
    assert values[0] == pytest.approx(1.0, abs=1e-9)
    assert values[1] == pytest.approx(0.5, abs=1e-6)
    with pytest.raises(ValueError, match="Kernel name must be one of"):
        KernelSpec(name="bogus", halfwidth_deg=10.0)
