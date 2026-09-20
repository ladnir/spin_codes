#!/usr/bin/env python3
"""Test sparse-type isolation against learned packet-type clustering."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_riffle_small_exact_vs_type_bound import (  # noqa: E402
    histogram_sequence_count,
    log_fraction,
    optimize_coefficient_bound,
    outer_histogram_distribution,
)
from probe_riffle_small_clustered_tilts import optimize_rows  # noqa: E402
from probe_riffle_small_typewise_bound import DEFAULT_RECEIPTS  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=4)
    parser.add_argument("--parity-symbols", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--distance", type=int, required=True)
    parser.add_argument(
        "--regions", type=int, nargs="+", default=[1, 2, 4, 8, 16, 32]
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    clustered_path = DEFAULT_RECEIPTS / (
        f"goal20_clustered_b{args.data_blocks}_p{args.parity_symbols}"
        f"_d{args.distance}.json"
    )
    clustered = json.loads(clustered_path.read_text(encoding="utf-8"))
    outer, _signatures = outer_histogram_distribution(
        args.data_blocks, args.parity_symbols
    )
    rows = tuple(
        (
            histogram,
            log_fraction(outer[histogram]),
            math.log(histogram_sequence_count(histogram)),
        )
        for histogram in sorted(outer)
    )
    cached = clustered["typewise_parameters"]
    if [list(row[0]) for row in rows] != [
        row["histogram_h0_through_h4"] for row in cached
    ]:
        raise RuntimeError("the cached type order does not match the outer distribution")

    type_parameters = np.asarray([row["parameters"] for row in cached], dtype=float)
    type_logs = np.asarray([row["log_bound_contribution"] for row in cached], dtype=float)
    packet_positions = int(clustered["packet_positions"])
    common_log, common_parameters, common_success = optimize_coefficient_bound(
        args.distance, packet_positions, rows, None
    )
    typewise_log = float(logsumexp(type_logs))
    order = np.argsort(-type_logs)

    region_results = []
    for regions in sorted(set(args.regions)):
        if not 1 <= regions <= len(rows):
            raise ValueError("invalid region count")
        isolated = order[: max(0, regions - 1)]
        remainder = order[max(0, regions - 1) :]
        contribution_logs = [float(type_logs[index]) for index in isolated]
        remainder_success = True
        if len(remainder):
            remainder_rows = tuple(rows[int(index)] for index in remainder)
            raw_weights = np.exp(type_logs[remainder] - float(np.max(type_logs[remainder])))
            mean_parameters = np.average(
                type_parameters[remainder], axis=0, weights=raw_weights
            )
            remainder_log, _parameters, remainder_success = optimize_rows(
                args.distance,
                packet_positions,
                remainder_rows,
                [mean_parameters, common_parameters],
            )
            contribution_logs.append(remainder_log)
        total_log = float(logsumexp(contribution_logs))
        recovered = common_log - total_log
        available = common_log - typewise_log
        region_results.append(
            {
                "regions": regions,
                "isolated_packet_types": len(isolated),
                "bound_log2": total_log / math.log(2.0),
                "improvement_over_common_bits": recovered / math.log(2.0),
                "fraction_of_typewise_improvement_recovered": (
                    recovered / available if available > 0 else 1.0
                ),
                "isolated_typewise_bound_mass": float(
                    np.sum(np.exp(type_logs[isolated] - typewise_log))
                ),
                "remainder_optimizer_success": remainder_success,
            }
        )

    binary_weights = np.asarray(
        [sum(index * count for index, count in enumerate(row[0])) for row in rows]
    )
    minimum_binary_weight = int(np.min(binary_weights))
    minimum_indices = np.flatnonzero(binary_weights == minimum_binary_weight)
    remainder = np.flatnonzero(binary_weights != minimum_binary_weight)
    minimum_logs = [float(type_logs[index]) for index in minimum_indices]
    if len(remainder):
        remainder_rows = tuple(rows[int(index)] for index in remainder)
        raw_weights = np.exp(type_logs[remainder] - float(np.max(type_logs[remainder])))
        mean_parameters = np.average(
            type_parameters[remainder], axis=0, weights=raw_weights
        )
        remainder_log, _parameters, remainder_success = optimize_rows(
            args.distance,
            packet_positions,
            remainder_rows,
            [mean_parameters, common_parameters],
        )
        minimum_logs.append(remainder_log)
    else:
        remainder_success = True
    minimum_shell_log = float(logsumexp(minimum_logs))

    payload = {
        "schema": "riffle-small-sparse-isolation-probe-v1",
        "evidence_label": "NUMERICAL_SPARSE_TYPE_ISOLATION_BOUNDS",
        "data_blocks": args.data_blocks,
        "parity_symbols": args.parity_symbols,
        "binary_output_length": 4 * packet_positions,
        "distance": args.distance,
        "packet_types": len(rows),
        "common_tilt_bound_log2": common_log / math.log(2.0),
        "typewise_tilt_bound_log2": typewise_log / math.log(2.0),
        "common_optimizer_success": common_success,
        "top_type_isolation": region_results,
        "minimum_binary_weight_isolation": {
            "minimum_outer_binary_weight": minimum_binary_weight,
            "isolated_packet_types": len(minimum_indices),
            "regions_including_bulk": len(minimum_indices) + int(bool(len(remainder))),
            "bound_log2": minimum_shell_log / math.log(2.0),
            "improvement_over_common_bits": (
                common_log - minimum_shell_log
            ) / math.log(2.0),
            "fraction_of_typewise_improvement_recovered": (
                (common_log - minimum_shell_log) / (common_log - typewise_log)
            ),
            "typewise_bound_mass_in_minimum_shell": float(
                np.sum(np.exp(type_logs[minimum_indices] - typewise_log))
            ),
            "bulk_optimizer_success": remainder_success,
        },
        "scope": (
            "The isolated types use their cached numerical typewise tilts. "
            "The remaining types share one re-optimized tilt. The selection by "
            "optimized contribution is diagnostic. Isolation by minimum binary "
            "weight is an explicit construction-level rule."
        ),
    }
    output = args.output or (
        DEFAULT_RECEIPTS
        / (
            f"goal20_sparse_isolation_b{args.data_blocks}_p{args.parity_symbols}"
            f"_d{args.distance}.json"
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output": str(output), **payload}, indent=2))


if __name__ == "__main__":
    main()
