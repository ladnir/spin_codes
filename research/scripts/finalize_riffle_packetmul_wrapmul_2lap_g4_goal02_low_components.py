#!/usr/bin/env python3
"""Authenticate the native Goal 02 low-component spectrum receipt."""

from __future__ import annotations

import hashlib
import json
from math import prod
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
SOURCE = (
    ROOT
    / "scripts"
    / "certify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components.cpp"
)
BINARY = SOURCE.with_suffix(".exe")
RAW = CANDIDATE / "receipts" / "goal02_low_components_primary_raw.json"
OUTPUT = CANDIDATE / "receipts" / "goal02_low_components_primary.json"
DEGREES = (1, 2, 4, 9, 10, 18, 20)
DIMENSION_CAP = 22


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    payload = json.loads(RAW.read_text())
    expected_masks = [
        mask
        for mask in range(1, 1 << len(DEGREES))
        if sum(DEGREES[index] for index in range(len(DEGREES)) if (mask >> index) & 1)
        <= DIMENSION_CAP
    ]
    rows = payload["rows"]
    if [int(row["component_support_mask_hex"], 16) for row in rows] != expected_masks:
        raise RuntimeError("low components finalize: support schedule mismatch")
    for row in rows:
        expected_count = prod((1 << degree) - 1 for degree in row["component_degrees"])
        if row["exact_nonzero_state_count"] != expected_count:
            raise RuntimeError("low components finalize: exact state count mismatch")
        if row["visited_state_count"] != expected_count:
            raise RuntimeError("low components finalize: incomplete orbit coverage")
    total_states = sum(row["exact_nonzero_state_count"] for row in rows)
    total_low = sum(row["low_response_state_count"] for row in rows)
    global_minimum = min(row["minimum_response_weight"] for row in rows)
    if (
        payload["support_count"] != len(expected_masks)
        or payload["total_exact_nonzero_states"] != total_states
        or payload["total_low_response_states"] != total_low
        or payload["global_minimum_response_weight"] != global_minimum
    ):
        raise RuntimeError("low components finalize: aggregate mismatch")
    payload.update(
        {
            "source_sha256": digest(SOURCE),
            "executable_sha256": digest(BINARY),
            "raw_receipt_sha256": digest(RAW),
            "finalizer_source_sha256": digest(Path(__file__).resolve()),
            "compiler": "MSVC 14.50.35717",
            "compiler_flags": [
                "/std:c++20",
                "/O2",
                "/EHsc",
                "/W4",
                "/permissive-",
            ],
            "coverage_statement": (
                "Every nonzero initial state whose exact irreducible-component "
                "support has total dimension at most 22 is enumerated once. "
                "Orbit rolling sums evaluate all 32772-node response windows."
            ),
        }
    )
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"support_count={len(rows)}")
    print(f"total_states={total_states}")
    print(f"global_minimum={global_minimum}")
    print(f"status={payload['result']}")


if __name__ == "__main__":
    main()
