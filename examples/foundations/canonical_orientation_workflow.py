from __future__ import annotations

from pytex import (
    FrameDomain,
    Handedness,
    Orientation,
    OrientationSet,
    ProvenanceRecord,
    ReferenceFrame,
    Rotation,
    SymmetrySpec,
)


def main() -> None:
    provenance = ProvenanceRecord.minimal("example", note="Foundational orientation workflow.")
    crystal = ReferenceFrame(
        name="crystal",
        domain=FrameDomain.CRYSTAL,
        axes=("a", "b", "c"),
        handedness=Handedness.RIGHT,
        provenance=provenance,
    )
    specimen = ReferenceFrame(
        name="specimen",
        domain=FrameDomain.SPECIMEN,
        axes=("x", "y", "z"),
        handedness=Handedness.RIGHT,
        provenance=provenance,
    )
    symmetry = SymmetrySpec.identity(reference_frame=crystal)
    o1 = Orientation(
        rotation=Rotation.from_bunge_euler(0.0, 0.0, 0.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        provenance=provenance,
    )
    o2 = Orientation(
        rotation=Rotation.from_bunge_euler(45.0, 35.0, 15.0),
        crystal_frame=crystal,
        specimen_frame=specimen,
        symmetry=symmetry,
        provenance=provenance,
    )

    orientation_set = OrientationSet.from_orientations([o1, o2])
    misorientation = o1.misorientation_to(o2)

    print("Orientation set length:", len(orientation_set))
    print("Misorientation matrix:")
    print(misorientation.as_matrix())


if __name__ == "__main__":
    main()
