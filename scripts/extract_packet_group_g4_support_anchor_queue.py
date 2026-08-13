#!/usr/bin/env python3
"""Extract distinct feasible support-mesh anchors as a tuning queue."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hull-only", action="store_true")
    args = parser.parse_args()
    ledger = json.loads(args.ledger.read_text(encoding="utf-8"))
    rows = []
    seen = set()
    for stratum in ledger["support_strata"]:
        if not stratum.get("feasible"):
            continue
        mask = int(stratum["support_mask"])
        source = (
            stratum.get("exact_hull_vertices", [])
            if args.hull_only
            else stratum["anchor_profiles"]
        )
        for index, raw in enumerate(source):
            profile = tuple(int(value) for value in raw)
            if profile in seen:
                continue
            seen.add(profile)
            rows.append({"name": f"s{mask:02x}_anchor_{index:04d}", "profile": profile})
    output = {
        "schema": "packet-group-g4-support-anchor-tuning-queue-v1",
        "group_bits": 4,
        "hull_only": args.hull_only,
        "uncovered_integer_residuals": rows,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"support_anchor_profiles={len(rows)}")


if __name__ == "__main__":
    main()
