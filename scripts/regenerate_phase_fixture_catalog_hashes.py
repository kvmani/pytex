from __future__ import annotations

import hashlib
import json
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> int:
    repo_root = _repo_root()
    catalog_path = repo_root / "fixtures/phases/catalog.json"
    payload = json.loads(catalog_path.read_text(encoding="utf-8"))
    for entry in payload["fixtures"]:
        artifact_path = repo_root / entry["artifact_path"]
        metadata_path = repo_root / entry["metadata_path"]
        entry["artifact_sha256"] = _sha256(artifact_path)
        entry["metadata_sha256"] = _sha256(metadata_path)
    catalog_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Updated fixture hashes in {catalog_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
