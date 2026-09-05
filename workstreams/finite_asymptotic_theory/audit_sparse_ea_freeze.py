#!/usr/bin/env python3
"""Verify the frozen sparse-EA artifact inventory and certificate chain."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "SPARSE_EA_FREEZE_MANIFEST.json"


def main() -> int:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if payload.get("status") != "FROZEN":
        raise AssertionError("freeze manifest is not frozen")
    for row in payload["artifacts"]:
        path = ROOT / row["file"]
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        if digest != row["sha256"]:
            raise AssertionError(f"hash mismatch: {row['file']}")
        if len(data) != int(row["bytes"]):
            raise AssertionError(f"size mismatch: {row['file']}")
    subprocess.run(
        [sys.executable, str(ROOT / "audit_one_stage_sparse_ea_certificate.py")],
        cwd=ROOT,
        check=True,
    )
    print(f"freeze artifacts: PASS ({len(payload['artifacts'])} files)")
    print("sparse-EA semantic certificate: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
