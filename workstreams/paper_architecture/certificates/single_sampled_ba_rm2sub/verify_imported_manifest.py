#!/usr/bin/env python3
"""Verify the relocated one-sampled BA/RM2Sub certificate snapshot."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[3]
MANIFEST = HERE / "SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(entry: dict[str, str], path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    actual = digest(path)
    if actual != entry["sha256"]:
        raise ValueError(
            f"hash mismatch for {path}: expected {entry['sha256']}, got {actual}"
        )


def main() -> None:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))

    for entry in payload["local_files"]:
        verify(entry, HERE / Path(entry["path"]).name)

    for entry in payload["companion_linear_time_files"]:
        verify(entry, HERE / "linear_time_audit" / Path(entry["path"]).name)

    for entry in payload["frozen_dependencies"]:
        verify(entry, REPO_ROOT / entry["path"])

    count = sum(
        len(payload[key])
        for key in (
            "local_files",
            "companion_linear_time_files",
            "frozen_dependencies",
        )
    )
    print(
        json.dumps(
            {
                "claim": payload["claim"],
                "manifest_sha256": digest(MANIFEST),
                "status": "verified",
                "verified_entries": count,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
