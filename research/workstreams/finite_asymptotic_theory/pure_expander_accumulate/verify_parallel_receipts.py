#!/usr/bin/env python3
"""Compare serial and parallel dominant-character receipts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("serial", type=Path)
    parser.add_argument("parallel", type=Path)
    args = parser.parse_args()
    serial = json.loads(args.serial.read_text(encoding="utf-8"))
    parallel = json.loads(args.parallel.read_text(encoding="utf-8"))
    if serial["shells"] != parallel["shells"]:
        raise AssertionError("serial and parallel shell rows differ")
    if serial["claim"] != parallel["claim"]:
        raise AssertionError("serial and parallel claims differ")
    result = {
        "shells_equal": True,
        "claim_equal": True,
        "serial_elapsed_seconds": serial["elapsed_seconds"],
        "parallel_elapsed_seconds": parallel["elapsed_seconds"],
        "speedup": serial["elapsed_seconds"] / parallel["elapsed_seconds"],
        "parallel_worker_elapsed_seconds": parallel["worker_elapsed_seconds"],
        "serial_sha256": sha256(args.serial),
        "parallel_sha256": sha256(args.parallel),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
