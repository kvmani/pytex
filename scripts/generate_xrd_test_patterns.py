"""Generate realistic synthetic experimental XRD pattern fixtures for testing.

Produces:
- fixtures/diffraction/experimental_ni_fcc_pattern.xy
- fixtures/diffraction/experimental_ni_fcc_pattern.xrdml

Features applied:
- Ni FCC standard phase (a = 3.52387 Å, Fm-3m)
- Laboratory Cu Ka radiation doublet (λ1 = 1.540598 Å, λ2 = 1.544426 Å, ratio 0.5)
- Instrument broadening via Caglioti relation: U = 0.005, V = -0.002, W = 0.008 (deg²)
- Pseudo-Voigt peak profile (η = 0.6 Lorentzian fraction)
- Designed curved background: low-angle air-scattering exponential + quadratic baseline
- Realistic Poisson counting noise with fixed seed for determinism.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pytex.app.phases import builtin_phase
from pytex.diffraction.xrd import RadiationSpec, generate_xrd_pattern


def main() -> None:
    phase = builtin_phase("ni_fcc").to_phase()
    caglioti = (0.005, -0.002, 0.008)

    pattern = generate_xrd_pattern(
        phase,
        radiation=RadiationSpec.cu_ka_doublet(),
        two_theta_range_deg=(20.0, 100.0),
        resolution_deg=0.02,
        caglioti_uvw=caglioti,
        profile="pseudo_voigt",
        pseudo_voigt_eta=0.6,
    )

    grid = pattern.two_theta_grid_deg
    # Background: low-angle air scattering + sample holder curvature
    bg = 220.0 * np.exp(-grid / 15.0) + 80.0 + 0.8 * (grid - 50.0) + 0.01 * (grid - 50.0) ** 2
    # Scale peak profile so max peak is ~12,000 counts
    max_p = float(np.max(pattern.intensity_grid))
    peak_counts = (pattern.intensity_grid / max_p) * 12000.0
    total = peak_counts + bg

    rng = np.random.default_rng(20260906)
    noisy = rng.poisson(total).astype(int)

    out_dir = Path("fixtures/diffraction")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Write .xy file
    xy_path = out_dir / "experimental_ni_fcc_pattern.xy"
    with xy_path.open("w", encoding="utf-8") as f:
        f.write("# pytex-powder-pattern: 1\n")
        f.write("# name: experimental-ni-fcc-standard\n")
        f.write("# specimen: Polycrystalline Nickel (fcc, a=3.52387 A)\n")
        f.write("# radiation: Cu Ka doublet (1.540598 A / 1.544426 A)\n")
        f.write(
            "# instrument_broadening: Caglioti U=0.005 V=-0.002 W=0.008 (pseudo-Voigt eta=0.6)\n"
        )
        f.write(
            "# background: Designed exponential air-scattering + quadratic sample curvature\n"
        )
        f.write("# intensity_unit: counts\n")
        f.write("# columns: two_theta_deg intensity\n")
        for angle, count in zip(grid, noisy, strict=True):
            f.write(f"{angle:.2f}  {count:d}\n")

    print(f"Generated {xy_path} ({len(grid)} points, max={int(np.max(noisy))})")

    # 2. Write .xrdml file
    xrdml_path = out_dir / "experimental_ni_fcc_pattern.xrdml"
    counts_str = " ".join(str(c) for c in noisy)
    xrdml_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<xrdMeasurements xmlns="http://www.panalytical.com/XRDML/1.0" version="1.0">
  <comment>
    <entry>Synthesized experimental powder XRD pattern of Ni FCC with Caglioti broadening.
    </entry>
  </comment>
  <sample>
    <id>Ni-FCC-Powder</id>
    <name>Nickel FCC standard</name>
  </sample>
  <xrdMeasurement measurementType="Scan">
    <usedWavelength>
      <kAlpha1 unit="Angstrom">1.540598</kAlpha1>
      <kAlpha2 unit="Angstrom">1.544426</kAlpha2>
      <ratioKAlpha2KAlpha1>0.5</ratioKAlpha2KAlpha1>
    </usedWavelength>
    <scan scanAxis="2Theta" status="Completed">
      <header>
        <startTimeStamp>2026-09-06T10:00:00</startTimeStamp>
        <author>
          <name>PyTex Synthetic Measurement Generator</name>
        </author>
      </header>
      <dataPoints>
        <positions axis="2Theta" unit="deg">
          <startPosition>{grid[0]:.4f}</startPosition>
          <endPosition>{grid[-1]:.4f}</endPosition>
        </positions>
        <commonCountingTime unit="seconds">1.00</commonCountingTime>
        <counts unit="counts">{counts_str}</counts>
      </dataPoints>
    </scan>
  </xrdMeasurement>
</xrdMeasurements>
"""
    xrdml_path.write_text(xrdml_content, encoding="utf-8")
    print(f"Generated {xrdml_path}")


if __name__ == "__main__":
    main()
