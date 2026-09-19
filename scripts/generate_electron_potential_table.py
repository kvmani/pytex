"""Regenerate the pinned electron-potential parametrization table used by multislice.

Writes ``src/pytex/diffraction/_data/electron_potential_parametrizations.json``
holding the two independent-atom electron scattering-factor fits that the
multislice engine in `pytex.diffraction.multislice` offers:

- **Lobato & Van Dyck (2014)**, Acta Cryst. A70, 636-649: five hydrogen-like
  terms ``f_e(g) = sum_i a_i (2 + b_i g^2) / (1 + b_i g^2)^2``, stored as the
  ``[a_1..a_5], [b_1..b_5]`` rows of the published table.
- **Kirkland (2010)**, Advanced Computing in Electron Microscopy, 2nd ed.,
  Appendix C: three Lorentzians and three Gaussians
  ``f_e(g) = sum_i a_i / (g^2 + b_i) + c_i exp(-d_i g^2)``, stored as the
  ``[a], [b], [c], [d]`` rows of the published table.

``g`` is the magnitude of the scattering vector in 1/angstrom (``g = 2 sin(theta) /
lambda``) and ``f_e`` is in angstrom.

The coefficients are copied from the JSON files abTEM distributes (abTEM is
GPL-3.0-or-later, like PyTex), which transcribe the published tables. Taking
them from abTEM rather than re-typing them from the papers is deliberate: the
multislice engine is validated against abTEM, and the parity test
``tests/unit/test_multislice_abtem_parity.py`` asserts the tables are identical
whenever abTEM is installed. Requires abTEM. Run when abTEM updates its tables:

    python scripts/generate_electron_potential_table.py
"""

from __future__ import annotations

import json
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def main() -> int:
    import abtem

    source = Path(abtem.__file__).resolve().parent / "parametrizations" / "data"
    tables: dict[str, dict[str, list[list[float]]]] = {}
    for name in ("lobato", "kirkland"):
        raw = json.loads((source / f"{name}.json").read_text(encoding="utf-8"))
        tables[name] = {
            symbol: [[float(value) for value in row] for row in rows]
            for symbol, rows in raw.items()
        }
    payload = {
        "description": (
            "Independent-atom electron scattering-factor parametrizations f_e(g) in angstrom, "
            "g = |scattering vector| in 1/angstrom."
        ),
        "parametrizations": {
            "lobato": {
                "formula": "f_e(g) = sum_i a_i (2 + b_i g^2) / (1 + b_i g^2)^2, i = 1..5",
                "rows": ["a", "b"],
                "reference": (
                    "Lobato, I. & Van Dyck, D. (2014). Acta Cryst. A70, 636-649. "
                    "doi:10.1107/S205327331401643X"
                ),
            },
            "kirkland": {
                "formula": "f_e(g) = sum_i a_i / (g^2 + b_i) + c_i exp(-d_i g^2), i = 1..3",
                "rows": ["a", "b", "c", "d"],
                "reference": (
                    "Kirkland, E. J. (2010). Advanced Computing in Electron Microscopy, "
                    "2nd ed., Springer, Appendix C. doi:10.1007/978-1-4419-6533-2"
                ),
            },
        },
        "source": f"abTEM {abtem.__version__} parametrizations/data/{{lobato,kirkland}}.json",
        "generator": "scripts/generate_electron_potential_table.py",
        "coefficients": tables,
    }
    target = (
        _repo_root() / "src/pytex/diffraction/_data/electron_potential_parametrizations.json"
    )
    target.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(f"Wrote {', '.join(f'{k}: {len(v)}' for k, v in tables.items())} entries to {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
