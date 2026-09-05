#!/usr/bin/env python3
"""Audit domination by the fully split profile after PacketMul.

For each recorded occupation, compare every region matrix with the profile
having one active outer block per packet group.  Also check every recorded
elementary packing edge in the reversed direction: the less packed endpoint
should dominate the more packed endpoint entrywise.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from audit_riffle_ldpcsplitstate_region_packing import packed_neighbors


ROOT = Path("constructions/riffle_ldpcsplitstate_packetmul_g4_t256_s64/receipts")
DEFAULT_OUTPUT = ROOT / "unpacking_domination_recorded_profiles.json"


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    rows = []
    total_edges = 0
    total_edge_violations = 0
    for path in args.receipts:
        payload = json.loads(path.read_text(encoding="utf-8"))
        for occupation_row in payload["occupation_rows"]:
            occupation = int(occupation_row["occupation"])
            unique = {
                tuple(row["profile_n1_n2_n3_n4"]): row
                for row in occupation_row["profiles"]
            }
            distinct_profile = (occupation, 0, 0, 0)
            distinct = np.asarray(
                unique[distinct_profile]["region_matrix"], dtype=np.float64
            )
            profile_violations = 0
            maximum_profile_ratio = 0.0
            maximum_profile_witness = None
            for profile, row in unique.items():
                matrix = np.asarray(row["region_matrix"], dtype=np.float64)
                ratio = float(np.max(matrix / distinct))
                if profile != distinct_profile and ratio > maximum_profile_ratio:
                    maximum_profile_ratio = ratio
                    maximum_profile_witness = list(profile)
                if bool(
                    np.any(
                        matrix
                        > args.absolute_tolerance
                        + (1.0 + args.relative_tolerance) * distinct
                    )
                ):
                    profile_violations += 1

            edges = 0
            edge_violations = 0
            maximum_edge_ratio = 0.0
            maximum_edge_witness = None
            for less_packed, less_packed_row in unique.items():
                less_packed_matrix = np.asarray(
                    less_packed_row["region_matrix"], dtype=np.float64
                )
                for more_packed in packed_neighbors(less_packed):
                    if more_packed not in unique:
                        continue
                    edges += 1
                    more_packed_matrix = np.asarray(
                        unique[more_packed]["region_matrix"], dtype=np.float64
                    )
                    ratio = float(np.max(more_packed_matrix / less_packed_matrix))
                    if ratio > maximum_edge_ratio:
                        maximum_edge_ratio = ratio
                        maximum_edge_witness = {
                            "less_packed": list(less_packed),
                            "more_packed": list(more_packed),
                        }
                    if bool(
                        np.any(
                            more_packed_matrix
                            > args.absolute_tolerance
                            + (1.0 + args.relative_tolerance) * less_packed_matrix
                        )
                    ):
                        edge_violations += 1
            total_edges += edges
            total_edge_violations += edge_violations
            rows.append(
                {
                    "occupation": occupation,
                    "recorded_unique_profiles": len(unique),
                    "profile_domination_violations": profile_violations,
                    "maximum_profile_to_distinct_ratio": maximum_profile_ratio,
                    "maximum_profile_ratio_witness": maximum_profile_witness,
                    "recorded_elementary_edges": edges,
                    "reversed_edge_violations": edge_violations,
                    "maximum_more_to_less_packed_ratio": maximum_edge_ratio,
                    "maximum_edge_ratio_witness": maximum_edge_witness,
                    "source": str(path),
                }
            )
            print(
                f"occupation,{occupation},profiles,{len(unique)},"
                f"profile_violations,{profile_violations},edges,{edges},"
                f"edge_violations,{edge_violations},"
                f"max_edge_ratio,{maximum_edge_ratio:.12g}",
                flush=True,
            )
    return {
        "schema": "riffle-ldpcsplitstate-packetmul-unpacking-audit-v1",
        "candidate": "Riffle LDPCSplitState PacketMul g=4 t=256 s=64",
        "rows": rows,
        "summary": {
            "recorded_elementary_edges": total_edges,
            "reversed_edge_violations": total_edge_violations,
        },
        "scope": (
            "Float64 comparison of exact recorded profile matrices. The "
            "result is diagnostic, not outward-rounded, and covers only the "
            "profiles in the supplied receipts."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--receipts",
        type=Path,
        nargs="+",
        default=[ROOT / "shared_group_low_occupation_profiles.json"],
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-300)
    parser.add_argument("--relative-tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
