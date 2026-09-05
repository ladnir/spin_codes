#!/usr/bin/env python3
"""Cluster packet types and optimize one coefficient tilt per cluster.

The calculation is a numerical diagnostic.  Every reported cluster value is
still a valid coefficient upper bound because its tilt is evaluated on the
complete set of packet types assigned to that cluster.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from analyze_riffle_small_exact_vs_type_bound import (
    ROOT,
    coefficient_bound_objective,
    histogram_sequence_count,
    log_fraction,
    optimize_coefficient_bound,
    outer_histogram_distribution,
)
from probe_riffle_small_typewise_bound import (
    DEFAULT_RECEIPTS,
    exact_log_from_receipt,
    single_type_objective,
)


PARAMETER_BOUNDS = [(-12.0, 3.0)] + [(-14.0, 14.0)] * 4
PARAMETER_LOWER = np.array([bound[0] for bound in PARAMETER_BOUNDS])
PARAMETER_UPPER = np.array([bound[1] for bound in PARAMETER_BOUNDS])


def optimize_rows(
    distance: int,
    packet_positions: int,
    rows: tuple[tuple[tuple[int, ...], float, float], ...],
    starts: list[np.ndarray],
) -> tuple[float, np.ndarray, bool]:
    results = []
    for start in starts:
        start = np.clip(start, PARAMETER_LOWER, PARAMETER_UPPER)
        result = minimize(
            coefficient_bound_objective,
            start,
            args=(distance, packet_positions, rows),
            method="Nelder-Mead",
            bounds=PARAMETER_BOUNDS,
            options={"maxiter": 3000, "xatol": 1e-9, "fatol": 1e-9},
        )
        if not result.success:
            fallback = minimize(
                coefficient_bound_objective,
                result.x,
                args=(distance, packet_positions, rows),
                method="Powell",
                bounds=PARAMETER_BOUNDS,
                options={"maxiter": 2000, "xtol": 1e-9, "ftol": 1e-9},
            )
            if float(fallback.fun) <= float(result.fun):
                result = fallback
        results.append(result)
    best = min(results, key=lambda result: float(result.fun))
    parameters = np.clip(best.x, PARAMETER_LOWER, PARAMETER_UPPER)
    value = coefficient_bound_objective(
        parameters, distance, packet_positions, rows
    )
    return float(value), parameters, bool(best.success)


def optimize_individual_types(
    distance: int,
    packet_positions: int,
    rows: tuple[tuple[tuple[int, ...], float, float], ...],
    common_parameters: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[bool]]:
    parameters = np.empty((len(rows), 5), dtype=np.float64)
    log_contributions = np.empty(len(rows), dtype=np.float64)
    successes = []
    for index, (histogram, log_outer, log_total) in enumerate(rows):
        result = minimize(
            single_type_objective,
            common_parameters,
            args=(distance, packet_positions, histogram, log_total),
            method="Nelder-Mead",
            bounds=PARAMETER_BOUNDS,
            options={"maxiter": 2500, "xatol": 1e-9, "fatol": 1e-9},
        )
        if not result.success:
            fallback = minimize(
                single_type_objective,
                result.x,
                args=(distance, packet_positions, histogram, log_total),
                method="Powell",
                bounds=PARAMETER_BOUNDS,
                options={"maxiter": 2000, "xtol": 1e-9, "ftol": 1e-9},
            )
            if float(fallback.fun) <= float(result.fun):
                result = fallback
        parameters[index] = np.clip(result.x, PARAMETER_LOWER, PARAMETER_UPPER)
        log_contributions[index] = log_outer + single_type_objective(
            parameters[index], distance, packet_positions, histogram, log_total
        )
        successes.append(bool(result.success))
    return parameters, log_contributions, successes


def weighted_kmeans(
    features: np.ndarray,
    weights: np.ndarray,
    cluster_count: int,
    maximum_iterations: int = 200,
) -> np.ndarray:
    """Deterministic weighted Lloyd iteration with farthest-point seeding."""
    if not 1 <= cluster_count <= len(features):
        raise ValueError("invalid cluster count")
    weighted_mean = np.average(features, axis=0, weights=weights)
    distances = np.sum((features - weighted_mean) ** 2, axis=1)
    first = int(np.argmax(weights * distances))
    centers = [features[first].copy()]
    while len(centers) < cluster_count:
        nearest = np.min(
            np.stack(
                [np.sum((features - center) ** 2, axis=1) for center in centers]
            ),
            axis=0,
        )
        candidate = int(np.argmax(weights * nearest))
        if any(np.array_equal(features[candidate], center) for center in centers):
            candidate = int(np.argmax(nearest))
        centers.append(features[candidate].copy())
    centers_array = np.stack(centers)

    assignments = np.full(len(features), -1, dtype=np.int64)
    for _ in range(maximum_iterations):
        squared = np.stack(
            [np.sum((features - center) ** 2, axis=1) for center in centers_array],
            axis=1,
        )
        following = np.argmin(squared, axis=1)
        if np.array_equal(following, assignments):
            break
        assignments = following
        for cluster in range(cluster_count):
            mask = assignments == cluster
            if np.any(mask):
                centers_array[cluster] = np.average(
                    features[mask], axis=0, weights=weights[mask]
                )
    if len(set(map(int, assignments))) != cluster_count:
        raise RuntimeError("weighted k-means produced an empty cluster")
    return assignments


def cluster_bound(
    distance: int,
    packet_positions: int,
    rows: tuple[tuple[tuple[int, ...], float, float], ...],
    assignments: np.ndarray,
    type_parameters: np.ndarray,
    type_weights: np.ndarray,
    common_parameters: np.ndarray,
) -> tuple[float, list[dict[str, object]]]:
    cluster_logs = []
    summaries = []
    for cluster in sorted(set(map(int, assignments))):
        indices = np.flatnonzero(assignments == cluster)
        cluster_rows = tuple(rows[int(index)] for index in indices)
        mean_parameters = np.average(
            type_parameters[indices], axis=0, weights=type_weights[indices]
        )
        log_bound, parameters, success = optimize_rows(
            distance,
            packet_positions,
            cluster_rows,
            [mean_parameters, common_parameters],
        )
        cluster_logs.append(log_bound)
        summaries.append(
            {
                "cluster": cluster,
                "packet_types": len(indices),
                "normalized_typewise_bound_mass": float(np.sum(type_weights[indices])),
                "log_bound": log_bound,
                "parameters": [float(value) for value in parameters],
                "optimizer_success": success,
            }
        )
    return float(logsumexp(cluster_logs)), summaries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=4)
    parser.add_argument("--parity-symbols", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--distance", type=int, required=True)
    parser.add_argument(
        "--clusters", type=int, nargs="+", default=[1, 2, 4, 8, 16, 32]
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

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
    packet_positions = 2 * (args.data_blocks + args.parity_symbols)
    common_log, common_parameters, common_success = optimize_coefficient_bound(
        args.distance, packet_positions, rows, None
    )
    type_parameters, type_logs, type_successes = optimize_individual_types(
        args.distance, packet_positions, rows, common_parameters
    )
    typewise_log = float(logsumexp(type_logs))

    raw_weights = np.exp(type_logs - float(np.max(type_logs)))
    type_weights = raw_weights / float(np.sum(raw_weights))
    feature_mean = np.average(type_parameters, axis=0, weights=type_weights)
    feature_variance = np.average(
        (type_parameters - feature_mean) ** 2, axis=0, weights=type_weights
    )
    feature_scale = np.sqrt(np.maximum(feature_variance, 1e-12))
    features = (type_parameters - feature_mean) / feature_scale

    requested = sorted(set(args.clusters))
    if requested[0] != 1:
        requested.insert(0, 1)
    if requested[-1] > len(rows):
        raise ValueError("cluster count exceeds the number of packet types")

    cluster_rows = []
    best_at_most = math.inf
    for cluster_count in requested:
        if cluster_count == 1:
            log_bound = common_log
            summaries = [
                {
                    "cluster": 0,
                    "packet_types": len(rows),
                    "normalized_typewise_bound_mass": 1.0,
                    "log_bound": common_log,
                    "parameters": [float(value) for value in common_parameters],
                    "optimizer_success": common_success,
                }
            ]
        else:
            assignments = weighted_kmeans(
                features, type_weights, cluster_count
            )
            log_bound, summaries = cluster_bound(
                args.distance,
                packet_positions,
                rows,
                assignments,
                type_parameters,
                type_weights,
                common_parameters,
            )
        best_at_most = min(best_at_most, log_bound)
        recovered = common_log - best_at_most
        available = common_log - typewise_log
        cluster_rows.append(
            {
                "clusters": cluster_count,
                "raw_cluster_bound_log2": log_bound / math.log(2.0),
                "best_bound_with_at_most_this_many_clusters_log2": (
                    best_at_most / math.log(2.0)
                ),
                "improvement_over_common_bits": recovered / math.log(2.0),
                "fraction_of_typewise_improvement_recovered": (
                    recovered / available if available > 0 else 1.0
                ),
                "cluster_summaries": summaries,
            }
        )

    exact_log2 = exact_log_from_receipt(
        args.data_blocks, args.parity_symbols, args.distance
    )
    if typewise_log / math.log(2.0) + 1e-8 < exact_log2:
        raise RuntimeError("the typewise bound fell below the exact tail")
    if any(
        row["best_bound_with_at_most_this_many_clusters_log2"] + 1e-8 < exact_log2
        for row in cluster_rows
    ):
        raise RuntimeError("a clustered bound fell below the exact tail")

    payload = {
        "schema": "riffle-small-clustered-tilt-probe-v1",
        "evidence_label": "NUMERICAL_CLUSTERED_COEFFICIENT_UPPER_BOUNDS",
        "data_blocks": args.data_blocks,
        "parity_symbols": args.parity_symbols,
        "packet_positions": packet_positions,
        "binary_output_length": 4 * packet_positions,
        "distance": args.distance,
        "packet_types": len(rows),
        "exact_cumulative_log2": exact_log2,
        "common_tilt_bound_log2": common_log / math.log(2.0),
        "typewise_tilt_bound_log2": typewise_log / math.log(2.0),
        "available_typewise_improvement_bits": (
            common_log - typewise_log
        ) / math.log(2.0),
        "typewise_optimizer_failures": sum(not value for value in type_successes),
        "feature": (
            "Weighted k-means on standardized per-type optimal log eta and "
            "four log packet fugacities. Weights are proportional to each "
            "type's optimized bound contribution."
        ),
        "cluster_results": cluster_rows,
        "typewise_parameters": [
            {
                "histogram_h0_through_h4": list(rows[index][0]),
                "parameters": [float(value) for value in type_parameters[index]],
                "log_bound_contribution": float(type_logs[index]),
                "optimizer_success": type_successes[index],
            }
            for index in range(len(rows))
        ],
        "scope": (
            "The cluster partition is learned from per-type numerical optima. "
            "Every evaluated cluster tilt gives a valid upper bound for its "
            "assigned finite type set. The learned partition is diagnostic, "
            "not yet a scalable or certified partition rule."
        ),
    }
    output = args.output or (
        DEFAULT_RECEIPTS
        / (
            f"goal20_clustered_b{args.data_blocks}_p{args.parity_symbols}"
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
                "exact_log2": exact_log2,
                "common_log2": common_log / math.log(2.0),
                "typewise_log2": typewise_log / math.log(2.0),
                "cluster_results": [
                    {
                        key: row[key]
                        for key in (
                            "clusters",
                            "best_bound_with_at_most_this_many_clusters_log2",
                            "improvement_over_common_bits",
                            "fraction_of_typewise_improvement_recovered",
                        )
                    }
                    for row in cluster_rows
                ],
                "status": "PASS",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
