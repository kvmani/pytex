"""Worked examples: averaging axes, where the vector mean does not work.

A crystal direction and its negative are the same physical axis, so axial
data carries an arbitrary sign and the arithmetic mean partially cancels.
The orientation tensor T = (1/n) sum v v^T is invariant under v -> -v, which
is why it is the right summary, and its eigenvalues have exact closed forms
at the three limiting distributions:

    uniform (1/3, 1/3, 1/3)    girdle (0, 1/2, 1/2)    cluster (0, 0, 1)

These examples check the closed forms and demonstrate the failure the tensor
exists to avoid: axes tightly clustered about z but randomly signed, where
the normalized resultant comes out wrong in sign and several degrees off
axis while the tensor recovers z.

See :doc:`../../theory/directional_statistics_and_mean_axes`.
"""

from __future__ import annotations

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

SPHERE_SETUP = """
import numpy as np
from pytex.core import crystal_frame
from pytex.core.sphere import SphericalVectorSet

frame = crystal_frame()


def axes(values, antipodal=True):
    return SphericalVectorSet.from_vectors(
        values, reference_frame=frame, antipodal=antipodal
    )
"""

_THETA = SymbolUse(
    r"\boldsymbol{\Theta}",
    "Orientation tensor (1/n) sum v v^T; the second moment of a direction set.",
)
_LAMBDA = SymbolUse(
    r"\lambda_{1} \le \lambda_{2} \le \lambda_{3}",
    "Eigenvalues of the orientation tensor; non-negative and summing to one.",
)

_THEORY = SeeAlso(
    "Directional statistics and mean axes",
    "../../theory/directional_statistics_and_mean_axes",
)


GIRDLE_AND_CLUSTER_EIGENVALUES = WorkedExample(
    id="directional-orientation-tensor-limiting-eigenvalues",
    title="A girdle gives eigenvalues (0, 1/2, 1/2) and a cluster (0, 0, 1)",
    domain="core",
    scenario=(
        "Build two limiting direction sets - one spread uniformly around a "
        "great circle in the xy-plane, one with every axis parallel to z - and "
        "read the eigenvalues of their orientation tensors. The closed forms "
        "are exact: a girdle has <cos^2> = <sin^2> = 1/2 around the circle so "
        "its eigenvalues are (0, 1/2, 1/2), and a perfect cluster gives "
        "(0, 0, 1). Together with the uniform case (1/3, 1/3, 1/3) these are "
        "the three corners of the eigenvalue triangle, and they let three "
        "numbers classify a fabric without contouring anything. The example "
        "returns both eigenvalue triples."
    ),
    setup=SPHERE_SETUP,
    code=(
        "rng = np.random.default_rng(4)\n"
        "angle = rng.random(100000) * 2.0 * np.pi\n"
        "girdle = np.stack(\n"
        "    [np.cos(angle), np.sin(angle), np.zeros_like(angle)], axis=-1\n"
        ")\n"
        "cluster = np.tile([0.0, 0.0, 1.0], (1000, 1))\n"
        "result = np.array(\n"
        "    [\n"
        "        np.linalg.eigvalsh(np.asarray(axes(girdle).orientation_tensor())),\n"
        "        np.linalg.eigvalsh(np.asarray(axes(cluster).orientation_tensor())),\n"
        "    ]\n"
        ")"
    ),
    expected=[[0.0, 0.5, 0.5], [0.0, 0.0, 1.0]],
    unit="",
    tolerance=5e-3,
    reference=(
        "Analytic: a uniform girdle in a plane has second moments "
        "<cos^2> = <sin^2> = 1/2 and zero out of plane, giving (0, 1/2, 1/2); "
        "identical parallel axes give (0, 0, 1)."
    ),
    citation=(
        "Woodcock, Specification of fabric shapes using an eigenvalue method, "
        "GSA Bulletin 88 (1977) 1231-1236."
    ),
    symbols=(_THETA, _LAMBDA),
    see_also=(_THEORY,),
    result_format="{:.4f}",
)


ORIENTATION_TENSOR_HAS_UNIT_TRACE = WorkedExample(
    id="directional-orientation-tensor-unit-trace",
    title="The orientation tensor of unit vectors has trace exactly one",
    domain="core",
    scenario=(
        "Compute the trace of the orientation tensor for a set of random unit "
        "directions. Because each vector is normalized, the diagonal sums to "
        "the mean of |v|^2, which is exactly one - so the eigenvalues always "
        "sum to one and live on a triangle. It is a free check on any "
        "implementation of the tensor, and it is what makes the three "
        "eigenvalues directly comparable between datasets of different size."
    ),
    setup=SPHERE_SETUP,
    code=(
        "rng = np.random.default_rng(4)\n"
        "values = rng.normal(size=(50000, 3))\n"
        "values = values / np.linalg.norm(values, axis=1, keepdims=True)\n"
        "tensor = np.asarray(axes(values).orientation_tensor())\n"
        "result = float(np.trace(tensor))"
    ),
    expected=1.0,
    unit="",
    tolerance=1e-12,
    reference=(
        "Analytic identity: tr((1/n) sum v v^T) = (1/n) sum |v|^2 = 1 for unit "
        "vectors, independent of the distribution."
    ),
    citation=(
        "Fisher, Lewis and Embleton, Statistical Analysis of Spherical Data "
        "(CUP 1987)."
    ),
    symbols=(_THETA,),
    see_also=(_THEORY,),
    result_format="{:.12f}",
)


RANDOMLY_SIGNED_AXES_DEFEAT_THE_VECTOR_MEAN = WorkedExample(
    id="directional-mean-axis-of-randomly-signed-axes",
    title="The tensor recovers z from randomly signed axes; the resultant does not",
    domain="core",
    scenario=(
        "Take 3000 axes tightly clustered about z - an unambiguous fibre - and "
        "give each an independent random sign, which is what a real axial "
        "measurement delivers. The normalized resultant comes out wrong in "
        "sign and several degrees off axis, and would change if the signs were "
        "redrawn. The orientation tensor is blind to sign, because "
        "(-v)(-v)^T = v v^T, so its principal eigenvector recovers z. The "
        "example returns the absolute z-component of the tensor mean, which "
        "must be 1."
    ),
    setup=SPHERE_SETUP,
    code=(
        "rng = np.random.default_rng(4)\n"
        "tight = np.stack(\n"
        "    [\n"
        "        0.2 * rng.normal(size=3000),\n"
        "        0.2 * rng.normal(size=3000),\n"
        "        np.ones(3000),\n"
        "    ],\n"
        "    axis=-1,\n"
        ")\n"
        "tight = tight / np.linalg.norm(tight, axis=1, keepdims=True)\n"
        "signed = tight * rng.choice([-1.0, 1.0], size=(3000, 1))\n"
        "result = float(abs(axes(signed).mean_direction()[2]))"
    ),
    expected=1.0,
    unit="",
    tolerance=5e-3,
    reference=(
        "Analytic: the axes are drawn about z, so the principal eigenvector of "
        "the orientation tensor is z and its z-component is 1. Random signs "
        "cannot affect it, since the tensor is invariant under v -> -v."
    ),
    citation=(
        "Mardia and Jupp, Directional Statistics (Wiley 2000) - axial data and "
        "the failure of the resultant."
    ),
    symbols=(_THETA,),
    see_also=(_THEORY,),
    result_format="{:.6f}",
)


GROUP = ExampleGroup(
    slug="directional-statistics",
    title="Directional statistics and mean axes",
    summary=(
        "Averaging axes rather than vectors: the orientation tensor has unit "
        "trace, its eigenvalues take exact values at the girdle and cluster "
        "limits, and it recovers a fibre axis from randomly signed data where "
        "the vector resultant fails outright."
    ),
    examples=(
        GIRDLE_AND_CLUSTER_EIGENVALUES,
        ORIENTATION_TENSOR_HAS_UNIT_TRACE,
        RANDOMLY_SIGNED_AXES_DEFEAT_THE_VECTOR_MEAN,
    ),
)
