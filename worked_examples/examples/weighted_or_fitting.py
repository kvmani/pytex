"""Measured-pair weighting checked against analytic identities."""

from ..framework import ExampleGroup, SeeAlso, SymbolUse, WorkedExample

SETUP = """
import numpy as np
from pytex import OrientationRelationship, OrientationSet, Rotation, specimen_frame
from pytex.core import fit_orientation_relationship
from pytex.app.phases import builtin_phase

parent = builtin_phase('austenite_fcc').to_phase()
child = builtin_phase('fe_bcc').to_phase()
nominal = OrientationRelationship.from_kurdjumov_sachs_correspondence(
    parent_phase=parent, child_phase=child,
)
angles = np.radians([0.0, 2.0])
measured = np.stack([
    Rotation.from_axis_angle([1, 0, 0], angle).as_matrix()
    @ nominal.parent_to_child_rotation.as_matrix() for angle in angles
])
parents = OrientationSet.from_matrices(
    np.repeat(np.eye(3)[None], 2, axis=0), specimen_frame=specimen_frame(), phase=parent,
)
children = OrientationSet.from_matrices(
    measured.transpose(0, 2, 1), specimen_frame=specimen_frame(), phase=child,
)
"""

_CITATION = "Markley et al. (2007), Averaging Quaternions, doi:10.2514/1.28949."
_LINKS = (SeeAlso("OR determination", "../../theory/orientation_relationship_determination"),)
_SYMBOLS = (SymbolUse(r"w_i", "Normalized scalar evidence weight of measured pair i."),)

GROUP = ExampleGroup(
    slug="weighted-or-fitting",
    title="Weighted orientation-relationship fitting",
    summary=("Retain every measured pair while controlling its influence with declared evidence "
           "weights. These checks use a common-axis circular identity and an exact exclusion."),
    examples=(
        WorkedExample(
            id="or-weighted-fit-circular-identity",
            title="A weighted quaternion fit agrees with the circular mean",
            domain="transformation",
            scenario="Two pairs differ by two degrees about one axis. Give the first three times the evidence weight and compare the fitted displacement against the analytic circular mean.",
            setup=SETUP,
            code="""
report = fit_orientation_relationship(parents, children, nominal, pair_weights=[3, 1])
analytic_angle = np.degrees(np.arctan2(np.sin(angles[1]), 3 + np.cos(angles[1])))
result = abs(report.residuals_deg[0] - analytic_angle)
""",
            expected=0.0, unit="deg", tolerance=1e-10,
            reference="For common-axis rotations the weighted chordal optimum has angle arg(sum(w*exp(i*angle))); the two independent formulas must agree.",
            citation=_CITATION, symbols=_SYMBOLS, see_also=_LINKS,
        ),
        WorkedExample(
            id="or-weighted-fit-excluded-residual",
            title="Excluding a pair leaves its residual visible",
            domain="transformation",
            scenario="Exclude the perturbed second pair while retaining it for review. The first exact pair fixes the rotation; the second must still report its imposed two-degree residual.",
            setup=SETUP,
            code="""
report = fit_orientation_relationship(parents, children, nominal, pair_weights=[1, 0])
assert report.weighted_mean_residual_deg < 1e-10
assert report.effective_pair_count == 1.0
result = report.residuals_deg[1]
""",
            expected=2.0, unit="deg", tolerance=1e-10,
            reference="The excluded pair was constructed by a two-degree rotation relative to the only positive-weight pair. It contributes no term to the scatter matrix.",
            citation=_CITATION, symbols=_SYMBOLS, see_also=_LINKS,
        ),
    ),
)
