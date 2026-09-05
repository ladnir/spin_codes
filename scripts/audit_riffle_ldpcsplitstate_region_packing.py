#!/usr/bin/env python3
"""Audit region-matrix domination by the maximally packed profile.

The audit reads exact profile receipts that share one Chernoff point.  For
each occupation, it compares every recorded region matrix with the matrix of
the maximally packed profile.  The audit is finite evidence for the proposed
region-level packing lemma; it does not cover unrecorded profiles.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path("constructions/riffle_ldpcsplitstate_g4_t256_s64/receipts")
DEFAULT_OUTPUT = ROOT / "region_packing_recorded_profiles_audit.json"


def packed_profile(occupation: int) -> list[int]:
    quads, remainder = divmod(occupation, 4)
    result = [0, 0, 0, quads]
    if remainder:
        result[remainder - 1] = 1
    return result


def packed_neighbors(profile: tuple[int, int, int, int]):
    """Yield profiles after one elementary move toward parts of width four."""
    n1, n2, n3, n4 = profile
    if n1 >= 2:
        yield (n1 - 2, n2 + 1, n3, n4)
    if n1 >= 1 and n2 >= 1:
        yield (n1 - 1, n2 - 1, n3 + 1, n4)
    if n1 >= 1 and n3 >= 1:
        yield (n1 - 1, n2, n3 - 1, n4 + 1)
    if n2 >= 2:
        yield (n1, n2 - 2, n3, n4 + 1)
    if n2 >= 1 and n3 >= 1:
        yield (n1 + 1, n2 - 1, n3 - 1, n4 + 1)
    if n3 >= 2:
        yield (n1, n2 + 1, n3 - 2, n4 + 1)


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    by_occupation: dict[int, list[dict[str, object]]] = {}
    sources: dict[int, list[str]] = {}
    for path in args.receipts:
        payload = json.loads(path.read_text(encoding="utf-8"))
        occupation_rows = payload.get("occupation_rows")
        if occupation_rows is None:
            occupation_rows = [
                {"occupation": payload["parameters"]["occupation"], "profiles": payload["rows"]}
            ]
        for occupation_row in occupation_rows:
            occupation = int(occupation_row["occupation"])
            by_occupation.setdefault(occupation, []).extend(
                occupation_row["profiles"]
            )
            sources.setdefault(occupation, []).append(str(path))

    rows = []
    for occupation, profiles in sorted(by_occupation.items()):
        target_profile = packed_profile(occupation)
        packed_rows = [
            row
            for row in profiles
            if row["profile_n1_n2_n3_n4"] == target_profile
        ]
        if not packed_rows:
            raise ValueError(f"packed profile missing at occupation {occupation}")
        packed = np.asarray(packed_rows[0]["region_matrix"], dtype=np.float64)
        unique: dict[tuple[int, ...], dict[str, object]] = {}
        for row in profiles:
            unique.setdefault(tuple(row["profile_n1_n2_n3_n4"]), row)

        maximum_excess = (-float("inf"), None)
        maximum_ratio = (0.0, None)
        violations = 0
        for profile, row in unique.items():
            matrix = np.asarray(row["region_matrix"], dtype=np.float64)
            excess = float(np.max(matrix - packed))
            ratio = float(np.max(matrix / packed))
            if excess > maximum_excess[0]:
                maximum_excess = (excess, profile)
            if list(profile) != target_profile and ratio > maximum_ratio[0]:
                maximum_ratio = (ratio, profile)
            violation_mask = matrix > (
                args.absolute_tolerance
                + (1.0 + args.relative_tolerance) * packed
            )
            if bool(np.any(violation_mask)):
                violations += 1

        edge_count = 0
        edge_violations = 0
        maximum_edge_ratio = (0.0, None)
        for source_profile, source_row in unique.items():
            source_matrix = np.asarray(source_row["region_matrix"], dtype=np.float64)
            for target_profile_tuple in packed_neighbors(source_profile):
                target_row = unique.get(target_profile_tuple)
                if target_row is None:
                    continue
                edge_count += 1
                target_matrix = np.asarray(
                    target_row["region_matrix"], dtype=np.float64
                )
                ratio = float(np.max(source_matrix / target_matrix))
                excess = float(np.max(source_matrix - target_matrix))
                violation_mask = source_matrix > (
                    args.absolute_tolerance
                    + (1.0 + args.relative_tolerance) * target_matrix
                )
                if bool(np.any(violation_mask)):
                    edge_violations += 1
                if ratio > maximum_edge_ratio[0]:
                    maximum_edge_ratio = (
                        ratio,
                        {
                            "source": list(source_profile),
                            "target": list(target_profile_tuple),
                            "maximum_entrywise_ratio": ratio,
                            "maximum_absolute_excess": excess,
                        },
                    )
        rows.append(
            {
                "occupation": occupation,
                "recorded_unique_profiles": len(unique),
                "packed_profile_n1_n2_n3_n4": target_profile,
                "entrywise_violations": violations,
                "maximum_absolute_excess": maximum_excess[0],
                "maximum_absolute_excess_profile": list(maximum_excess[1]),
                "maximum_entrywise_ratio": maximum_ratio[0],
                "maximum_ratio_profile": list(maximum_ratio[1]),
                "recorded_elementary_packing_edges": edge_count,
                "elementary_edge_violations": edge_violations,
                "maximum_elementary_edge_ratio": maximum_edge_ratio[0],
                "maximum_elementary_edge_witness": maximum_edge_ratio[1],
                "sources": sources[occupation],
            }
        )
        print(
            f"occupation,{occupation},profiles,{len(unique)},"
            f"violations,{violations},max_ratio,{maximum_ratio[0]:.12f}",
            f",edges,{edge_count},edge_violations,{edge_violations},"
            f"max_edge_ratio,{maximum_edge_ratio[0]:.12f}",
            flush=True,
        )
    return {
        "schema": "riffle-ldpcsplitstate-region-packing-recorded-audit-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "rows": rows,
        "scope": (
            "Exact floating-point comparison of region matrices already "
            "recorded at a common Chernoff point for each occupation. The "
            "result is not an outward-rounded proof and does not cover "
            "profiles absent from the source receipts."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-300)
    parser.add_argument("--relative-tolerance", type=float, default=1e-12)
    parser.add_argument(
        "--receipts",
        type=Path,
        nargs="+",
        default=[
            ROOT / "shared_group_low_occupation_profiles.json",
            ROOT / "shared_group_occupation16_profiles.json",
            ROOT / "shared_group_occupation128_near_packed.json",
            ROOT / "termination_fugacity_occupation128.json",
        ],
    )
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
