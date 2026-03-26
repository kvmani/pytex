from __future__ import annotations

from pytex import (
    ODF,
    CrystalPlane,
    FrameDomain,
    Handedness,
    InversePoleFigure,
    KernelSpec,
    Lattice,
    MillerIndex,
    Orientation,
    OrientationSet,
    Phase,
    PoleFigure,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)


def main() -> None:
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
    phase = Phase(name="fcc-demo", lattice=lattice, symmetry=symmetry, crystal_frame=crystal)

    orientations = OrientationSet.from_orientations(
        [
            Orientation(
                Rotation.identity(),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
                phase=phase,
            ),
            Orientation(
                Rotation.from_bunge_euler(45.0, 35.0, 15.0),
                crystal_frame=crystal,
                specimen_frame=specimen,
                symmetry=symmetry,
                phase=phase,
            ),
        ]
    )
    pole = CrystalPlane(miller=MillerIndex([1, 0, 0], phase=phase), phase=phase)

    pole_figure = PoleFigure.from_orientations(orientations, pole)
    inverse_pole_figure = InversePoleFigure.from_orientations(orientations, [1.0, 0.0, 0.0])
    odf = ODF.from_orientations(
        orientations,
        kernel=KernelSpec(name="de_la_vallee_poussin", halfwidth_deg=12.0),
    )

    print("Pole figure projected coordinates shape:", pole_figure.project().shape)
    print("Inverse pole figure projected coordinates shape:", inverse_pole_figure.project().shape)
    print("ODF self-density:", odf.evaluate(orientations))


if __name__ == "__main__":
    main()
