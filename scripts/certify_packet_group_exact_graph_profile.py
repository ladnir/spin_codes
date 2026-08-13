#!/usr/bin/env python3
"""Outward-certify one g=4 witness with an exact graph-averaged outer branch."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import certify_packet_group_triangle_ledger as certifier


def rows(path: Path) -> list[dict]:
    value = json.loads(path.read_text(encoding="utf-8"))
    return value if isinstance(value, list) else value.get("rows", [])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument(
        "--outer-type",
        choices=("exact_graph_linear_bl", "exact_graph_total_spectrum"),
        required=True,
    )
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--target-log2", type=float, default=-111.41506501642425)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profile = tuple(int(value) for value in args.profile.split(","))
    matching = [row for row in rows(args.artifact) if tuple(row.get("profile", [])) == profile]
    if not matching:
        raise ValueError("exact graph profile certifier: profile is absent")
    row = dict(min(matching, key=lambda value: float(value["combined_log2"])))
    row["outer_type"] = args.outer_type
    certifier.GROUP_BITS = 4
    certifier.BLOCK_ATOMS = 16
    certifier._OUTWARD_NORMALIZATION_CACHE.clear()
    hardened = certifier.harden_witness(
        row.get("name", "selected"),
        row,
        args.artifact,
        certifier.split_cap_table(),
        args.iterations,
    )
    value = certifier.evaluate_vertex(profile, hardened)
    report = {
        "status": "OUTWARD_CERTIFIED_G4_EXACT_GRAPH_PROFILE_WITNESS",
        "profile": profile,
        "source_artifact": str(args.artifact),
        "source_artifact_sha256": certifier._file_digest(args.artifact.resolve()),
        "source_witness_name": row.get("name"),
        "outer_type": args.outer_type,
        "combined_log2_interval": [str(value.lo), str(value.hi)],
        "target_log2": args.target_log2,
        "target_margin_lower_bits": str(value.lo.__class__(str(args.target_log2)) - value.hi),
        "closes_target": bool(value.hi <= value.lo.__class__(str(args.target_log2))),
        "hardening_report": hardened["report"],
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
