#!/usr/bin/env python3
"""Combine an exact minimum outer shell with clustered bulk tilt bounds."""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_riffle_small_exact_vs_type_bound import (  # noqa: E402
    histogram_sequence_count,
    log2_fraction,
    log_fraction,
    optimize_coefficient_bound,
    outer_histogram_distribution,
    target_accumulator_enumerators,
)
from probe_riffle_small_clustered_tilts import (  # noqa: E402
    cluster_bound,
    weighted_kmeans,
)
from probe_riffle_small_typewise_bound import DEFAULT_RECEIPTS  # noqa: E402


def logaddexp2(left: float, right: float) -> float:
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    maximum = max(left, right)
    return maximum + math.log2(2 ** (left - maximum) + 2 ** (right - maximum))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=4)
    parser.add_argument("--parity-symbols", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--distance", type=int, required=True)
    parser.add_argument(
        "--bulk-clusters", type=int, nargs="+", default=[1, 2, 4, 8, 16]
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
    ordered_histograms = sorted(outer)
    rows = tuple(
        (
            histogram,
            log_fraction(outer[histogram]),
            math.log(histogram_sequence_count(histogram)),
        )
        for histogram in ordered_histograms
    )
    cached = clustered["typewise_parameters"]
    if [list(histogram) for histogram in ordered_histograms] != [
        row["histogram_h0_through_h4"] for row in cached
    ]:
        raise RuntimeError("the cached type order does not match the outer distribution")

    packet_positions = 2 * (args.data_blocks + args.parity_symbols)
    accumulator = target_accumulator_enumerators(
        packet_positions, set(outer), args.distance
    )
    exact_by_input_weight: dict[int, Fraction] = {}
    for histogram in ordered_histograms:
        input_weight = sum(index * count for index, count in enumerate(histogram))
        inner_bad = sum(accumulator[histogram].values())
        contribution = (
            outer[histogram]
            * inner_bad
            / histogram_sequence_count(histogram)
        )
        exact_by_input_weight[input_weight] = (
            exact_by_input_weight.get(input_weight, Fraction(0)) + contribution
        )
    exact_total = sum(exact_by_input_weight.values(), Fraction(0))
    expected_exact_log2 = float(clustered["exact_cumulative_log2"])
    if not math.isclose(
        log2_fraction(exact_total), expected_exact_log2, rel_tol=1e-11, abs_tol=1e-11
    ):
        raise RuntimeError("the exact shell decomposition disagrees with the receipt")

    minimum_input_weight = min(
        sum(k * histogram[k] for k in range(5))
        for histogram in ordered_histograms
    )
    input_weights = np.asarray(
        [
            sum(k * histogram[k] for k in range(5))
            for histogram in ordered_histograms
        ],
        dtype=np.int64,
    )
    sparse_indices = np.array(
        [
            index
            for index, histogram in enumerate(ordered_histograms)
            if input_weights[index] == minimum_input_weight
        ],
        dtype=np.int64,
    )
    sparse_index_set = set(map(int, sparse_indices))
    bulk_indices = np.array(
        [
            index
            for index in range(len(rows))
            if index not in sparse_index_set
            and input_weights[index] <= 2 * args.distance
        ],
        dtype=np.int64,
    )
    impossible_indices = np.flatnonzero(input_weights > 2 * args.distance)
    exact_sparse = exact_by_input_weight[minimum_input_weight]
    exact_sparse_log2 = log2_fraction(exact_sparse)

    type_parameters = np.asarray([row["parameters"] for row in cached], dtype=float)
    type_logs = np.asarray([row["log_bound_contribution"] for row in cached], dtype=float)
    bulk_rows = tuple(rows[int(index)] for index in bulk_indices)
    bulk_common_log, bulk_common_parameters, bulk_common_success = (
        optimize_coefficient_bound(args.distance, packet_positions, bulk_rows, None)
    )
    raw_weights = np.exp(
        type_logs[bulk_indices] - float(np.max(type_logs[bulk_indices]))
    )
    type_weights = raw_weights / float(np.sum(raw_weights))
    bulk_parameters = type_parameters[bulk_indices]
    feature_mean = np.average(bulk_parameters, axis=0, weights=type_weights)
    feature_variance = np.average(
        (bulk_parameters - feature_mean) ** 2, axis=0, weights=type_weights
    )
    feature_scale = np.sqrt(np.maximum(feature_variance, 1e-12))
    features = (bulk_parameters - feature_mean) / feature_scale

    results = []
    best_hybrid = math.inf
    requested_cluster_counts = sorted(
        set(min(cluster_count, len(bulk_rows)) for cluster_count in args.bulk_clusters)
    )
    if not requested_cluster_counts or requested_cluster_counts[0] < 1:
        raise ValueError("invalid bulk cluster count")
    for cluster_count in requested_cluster_counts:
        if cluster_count == 1:
            bulk_log = bulk_common_log
            summaries = [
                {
                    "cluster": 0,
                    "packet_types": len(bulk_rows),
                    "normalized_typewise_bound_mass": 1.0,
                    "log_bound": bulk_log,
                    "parameters": [float(value) for value in bulk_common_parameters],
                    "optimizer_success": bulk_common_success,
                }
            ]
        else:
            assignments = weighted_kmeans(
                features, type_weights, cluster_count
            )
            bulk_log, summaries = cluster_bound(
                args.distance,
                packet_positions,
                bulk_rows,
                assignments,
                bulk_parameters,
                type_weights,
                bulk_common_parameters,
            )
        hybrid_log2 = logaddexp2(exact_sparse_log2, bulk_log / math.log(2.0))
        best_hybrid = min(best_hybrid, hybrid_log2)
        results.append(
            {
                "bulk_clusters": cluster_count,
                "bulk_bound_log2": bulk_log / math.log(2.0),
                "hybrid_bound_log2": hybrid_log2,
                "best_hybrid_with_at_most_this_many_bulk_clusters_log2": best_hybrid,
                "loss_over_complete_exact_tail_bits": hybrid_log2
                - expected_exact_log2,
                "cluster_summaries": summaries,
            }
        )

    payload = {
        "schema": "riffle-small-exact-shell-clustered-bulk-v1",
        "evidence_label": "EXACT_RATIONAL_SHELL_PLUS_NUMERICAL_BULK_BOUNDS",
        "data_blocks": args.data_blocks,
        "parity_symbols": args.parity_symbols,
        "binary_output_length": 4 * packet_positions,
        "distance": args.distance,
        "exact_complete_tail_log2": expected_exact_log2,
        "minimum_outer_binary_weight": minimum_input_weight,
        "minimum_shell_packet_types": len(sparse_indices),
        "eligible_nonsparse_packet_types": len(bulk_indices),
        "discarded_impossible_packet_types": len(impossible_indices),
        "deterministic_cutoff": "outer binary weight H <= 2D",
        "exact_minimum_shell_log2": exact_sparse_log2,
        "exact_minimum_shell_fraction_of_complete_tail": float(
            exact_sparse / exact_total
        ),
        "exact_tail_by_outer_binary_weight": {
            str(weight): {
                "numerator": str(contribution.numerator),
                "denominator": str(contribution.denominator),
                "log2": log2_fraction(contribution),
            }
            for weight, contribution in sorted(exact_by_input_weight.items())
            if contribution
        },
        "bulk_results": results,
        "scope": (
            "The minimum-outer-weight contribution is exact rational. Types "
            "with input weight above 2D are impossible because H<=2W. The "
            "remaining bulk uses numerical coefficient upper bounds on a "
            "learned finite partition. The hybrid values are valid for the "
            "evaluated parameters but are not yet certified against "
            "floating-point error."
        ),
    }
    output = args.output or (
        DEFAULT_RECEIPTS
        / (
            f"goal20_exact_shell_hybrid_b{args.data_blocks}_p{args.parity_symbols}"
            f"_d{args.distance}.json"
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "exact_complete_tail_log2": expected_exact_log2,
                "exact_minimum_shell_log2": exact_sparse_log2,
                "exact_minimum_shell_fraction": float(exact_sparse / exact_total),
                "bulk_results": [
                    {
                        key: row[key]
                        for key in (
                            "bulk_clusters",
                            "hybrid_bound_log2",
                            "loss_over_complete_exact_tail_bits",
                        )
                    }
                    for row in results
                ],
                "status": "PASS",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
