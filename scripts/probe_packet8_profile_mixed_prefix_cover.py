#!/usr/bin/env python3
"""Hierarchical mixed-witness cover of all ordered profile chambers.

A prefix ``(p0,...,pr)`` represents the convex union of all ordered chambers
beginning with that prefix.  Its candidate vertices are the proper uniform
prefixes and every uniform superset of the full prefix.  For prefixes starting
at class zero, the minimum-weight clipping points are included as well.

At each node an LP searches for one convex mixture of fixed logged witnesses
that is safe at every candidate vertex.  Convexity then covers the entire
prefix region.  Failed nodes split on the next packet class.  This is a
binary64 discovery/diagnostic; the selected mixtures still require outward
certification.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from probe_packet8_profile_ordered_chambers import (
    add_full_bijection,
    clipped_from_zero,
    uniform_subset,
)
from probe_packet8_profile_simplex_landscape import PROFILE_COUNT_LOG2, normalization_log2
from probe_packet8_shared_witness_io import load_shared_witness_arrays


ALL_CLASSES = (1 << 9) - 1


def prefix_point_keys(prefix: tuple[int, ...]) -> list[tuple[str, int]]:
    keys: set[tuple[str, int]] = set()
    mask = 0
    for packet_class in prefix:
        mask |= 1 << packet_class
        keys.add(("uniform", mask))

    remaining = ALL_CLASSES ^ mask
    addition = remaining
    while True:
        keys.add(("uniform", mask | addition))
        if addition == 0:
            break
        addition = (addition - 1) & remaining

    if prefix[0] == 0:
        keys.discard(("uniform", 1))
        keys.update(
            ("clipped", subset)
            for kind, subset in tuple(keys)
            if kind == "uniform" and subset != 1
        )
    return sorted(keys)


def solve_mixture(values: np.ndarray):
    witness_count = values.shape[1]
    objective = np.zeros(witness_count + 1)
    objective[-1] = 1.0
    result = linprog(
        objective,
        A_ub=np.hstack((values, -np.ones((len(values), 1)))),
        b_ub=np.zeros(len(values)),
        A_eq=np.asarray([[1.0] * witness_count + [0.0]]),
        b_eq=np.asarray([1.0]),
        bounds=[(0.0, None)] * witness_count + [(None, None)],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    weights = result.x[:-1]
    return float(np.max(values @ weights)), weights


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upgraded-cache", type=Path, required=True)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    parser.add_argument("--minimum-outer-weight", type=int, default=21)
    parser.add_argument("--max-nodes", type=int, default=10000)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--show-unresolved", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.max_nodes <= 0 or args.show_unresolved <= 0:
        raise SystemExit("mixed prefix cover: invalid node/row limit")

    names, constants, charges = load_shared_witness_arrays(
        args.upgraded_cache, args.extra_shared_report
    )
    full = add_full_bijection([])[0]
    names.append(full.name)
    constants = np.append(constants, full.constant_log2)
    charges = np.vstack((charges, full.linear_charge))
    target = -40.0 - PROFILE_COUNT_LOG2

    point_values: dict[tuple[str, int], np.ndarray] = {}
    for subset in range(1, 1 << 9):
        point = uniform_subset(subset)
        point_values[("uniform", subset)] = (
            constants - charges @ point - normalization_log2(point[None, :])[0]
        )
        if subset & 1 and subset != 1:
            point = clipped_from_zero(point, args.minimum_outer_weight)
            point_values[("clipped", subset)] = (
                constants - charges @ point - normalization_log2(point[None, :])[0]
            )

    stack = [(packet_class,) for packet_class in reversed(range(9))]
    processed = 0
    covered_chambers = 0
    covered_nodes = []
    unresolved = []
    print("packet-8 hierarchical mixed-witness prefix cover", flush=True)
    print(
        f"witnesses={len(constants)} target={target:.12f} max_nodes={args.max_nodes}",
        flush=True,
    )
    while stack and processed < args.max_nodes:
        prefix = stack.pop()
        keys = prefix_point_keys(prefix)
        values = np.vstack([point_values[key] for key in keys])
        score, weights = solve_mixture(values)
        processed += 1
        if score <= target:
            chambers = math.factorial(9 - len(prefix))
            covered_chambers += chambers
            active = np.flatnonzero(weights > 1e-9)
            covered_nodes.append(
                {
                    "prefix": list(prefix),
                    "chambers": chambers,
                    "score_log2": score,
                    "margin_bits": target - score,
                    "points": len(keys),
                    "mixture": [
                        {"name": names[int(index)], "weight": float(weights[index])}
                        for index in active
                    ],
                }
            )
        elif len(prefix) == 9:
            unresolved.append({"prefix": list(prefix), "score_log2": score})
        else:
            used = set(prefix)
            for packet_class in reversed(range(9)):
                if packet_class not in used:
                    stack.append(prefix + (packet_class,))

        if args.progress_every and processed % args.progress_every == 0:
            print(
                f"processed_nodes={processed} open_nodes={len(stack)} "
                f"covered_chambers={covered_chambers} "
                f"unresolved_leaf_chambers={len(unresolved)}",
                flush=True,
            )

    total_chambers = math.factorial(9)
    report = {
        "status": "DIAGNOSTIC_BINARY64_HIERARCHICAL_CONVEX_MIXTURES",
        "target_log2": target,
        "minimum_outer_weight": args.minimum_outer_weight,
        "processed_nodes": processed,
        "open_nodes": [list(prefix) for prefix in stack],
        "covered_chambers": covered_chambers,
        "total_chambers": total_chambers,
        "covered_nodes": covered_nodes,
        "unresolved_leaf_chambers": unresolved,
        "complete_cover": not stack and not unresolved and covered_chambers == total_chambers,
    }
    print(f"processed_nodes={processed}")
    print(f"covered_nodes={len(covered_nodes)}")
    print(f"open_nodes={len(stack)}")
    print(f"covered_chambers={covered_chambers}")
    print(f"total_chambers={total_chambers}")
    print(f"unresolved_leaf_chambers={len(unresolved)}")
    for row in sorted(unresolved, key=lambda item: item["score_log2"], reverse=True)[
        : args.show_unresolved
    ]:
        print(
            f"unresolved_score={row['score_log2']:.12f} prefix="
            + ",".join(map(str, row["prefix"]))
        )
    print(
        "complete_mixed_prefix_cover="
        + ("PASS" if report["complete_cover"] else "NO")
    )
    print("status=DIAGNOSTIC_BINARY64_HIERARCHICAL_CONVEX_MIXTURES")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
