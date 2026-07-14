"""Regenerate the pinned X-ray form-factor coefficient table from pymatgen.

Writes ``src/pytex/diffraction/_data/xray_scattering_factors.json`` with the
De Graef & McHenry ``(a_i, b_i)`` coefficients bundled in pymatgen's
``ATOMIC_SCATTERING_PARAMS`` (the table pymatgen's XRDCalculator uses), so
PyTex ``intensity_model="xray_tabulated"`` structure factors share provenance
with the pinned pymatgen external baselines. Requires the ``adapters`` extra
(pymatgen). Run whenever pymatgen updates its table:

    python scripts/generate_xray_scattering_factor_table.py
"""

from __future__ import annotations

import json
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def main() -> int:
    from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS

    coefficients = {
        symbol: [[float(a), float(b)] for a, b in rows]
        for symbol, rows in sorted(ATOMIC_SCATTERING_PARAMS.items())
    }
    payload = {
        "description": (
            "X-ray atomic form-factor coefficients (a_i, b_i) for the "
            "De Graef & McHenry parametrization "
            "f(s) = Z - 41.78214 s^2 sum_i a_i exp(-b_i s^2)."
        ),
        "source": "pymatgen ATOMIC_SCATTERING_PARAMS (pymatgen.analysis.diffraction.xrd)",
        "generator": "scripts/generate_xray_scattering_factor_table.py",
        "coefficients": coefficients,
    }
    target = _repo_root() / "src/pytex/diffraction/_data/xray_scattering_factors.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(coefficients)} element entries to {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
