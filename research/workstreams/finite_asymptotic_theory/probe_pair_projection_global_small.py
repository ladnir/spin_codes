#!/usr/bin/env python3
"""Globally optimize the iterated pair-projection envelope at length eight.

The earlier projection probe greedily minimizes each backward factorization.
Here every intermediate factor is one variable in a single convex log-domain
program.  SLSQP supplies a nearest-binary64 diagnostic for whether greediness,
rather than the projection inequality itself, caused the large loss.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from analyze_accumulator_pair_chain_small import compositions4, rank_two
from probe_pair_projection_holder_bound_small import (
    factorize,
    log_inverse_geometric_overlap,
)
from probe_triangle_holder_pair_bound_small import (
    direct_sum_spectrum,
    one_word_transition,
    triple_weights,
)


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "pair_projection_global_young_B8_probe.json"


def expand(values: np.ndarray, support: list[int], length: int) -> tuple[np.ndarray, np.ndarray]:
    size = len(support)
    log_f = np.full(length, -math.inf)
    log_g = np.full(length, -math.inf)
    for index, weight in enumerate(support):
        log_f[weight - 1] = values[index]
        log_g[weight - 1] = values[size + index]
    return log_f, log_g


def transform(log_values: np.ndarray, log_kernel: np.ndarray) -> np.ndarray:
    return logsumexp(log_kernel + log_values[None, :], axis=1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stages", type=int, default=3)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if not 1 <= args.stages <= 3:
        parser.error("this diagnostic supports one through three stages")

    length = 8
    messages = 1 << 4
    spectrum = direct_sum_spectrum(1)
    base = np.asarray(spectrum, dtype=float) / messages
    binary = one_word_transition(length)
    log_kernel = np.full_like(binary, -math.inf)
    positive = binary > 0.0
    log_kernel[positive] = (2.0 / 3.0) * np.log(binary[positive])
    distributions = [base]
    for _ in range(args.stages):
        distributions.append(distributions[-1] @ binary)
    supports = [
        [weight for weight in range(1, length + 1) if distributions[stage][weight] > 0]
        for stage in range(args.stages)
    ]
    types = [kind for kind in compositions4(length) if rank_two(kind)]
    penalties = np.asarray(
        [log_inverse_geometric_overlap(kind, length) for kind in types]
    )

    offsets = [0]
    for support in supports:
        offsets.append(offsets[-1] + 2 * len(support))

    rows = []
    for target_weight in range(1, length + 1):
        # Build the greedy feasible point from the existing factorizer.
        next_f = np.full(length, -math.inf)
        next_f[target_weight - 1] = 0.0
        next_g = np.zeros(length)
        greedy_by_stage: list[tuple[np.ndarray, np.ndarray] | None] = [None] * args.stages
        for stage in range(args.stages - 1, -1, -1):
            transformed_f = transform(next_f, log_kernel[1:, 1:])
            transformed_g = transform(next_g, log_kernel[1:, 1:])
            lower = np.asarray(
                [
                    transformed_f[triple_weights(kind)[0] - 1]
                    + transformed_f[triple_weights(kind)[1] - 1]
                    + transformed_g[triple_weights(kind)[2] - 1]
                    + penalties[index]
                    for index, kind in enumerate(types)
                ]
            )
            factor = factorize(lower, types, distributions[stage], "young")
            if factor is None:
                raise RuntimeError("greedy factorization has empty support")
            next_f = np.asarray(factor["f_log"])
            next_g = np.asarray(factor["g_log"])
            greedy_by_stage[stage] = (next_f.copy(), next_g.copy())

        initial_parts = []
        for stage, support in enumerate(supports):
            pair = greedy_by_stage[stage]
            assert pair is not None
            initial_parts.append(
                np.asarray(
                    [pair[0][weight - 1] for weight in support]
                    + [pair[1][weight - 1] for weight in support]
                )
            )
        initial = np.concatenate(initial_parts)

        def unpack(point: np.ndarray, stage: int) -> tuple[np.ndarray, np.ndarray]:
            return expand(point[offsets[stage] : offsets[stage + 1]], supports[stage], length)

        def objective(point: np.ndarray) -> float:
            log_f, log_g = unpack(point, 0)
            support = supports[0]
            log_prob = np.log(np.asarray([base[weight] for weight in support]))
            f_values = np.asarray([log_f[weight - 1] for weight in support])
            g_values = np.asarray([log_g[weight - 1] for weight in support])
            return float(
                logsumexp(log_prob + 2.0 * f_values)
                + logsumexp(log_prob + g_values)
            )

        constraint_stage_rows: list[list[int]] = []
        for stage, support in enumerate(supports):
            support_set = set(support)
            constraint_stage_rows.append(
                [
                    index
                    for index, kind in enumerate(types)
                    if all(weight in support_set for weight in triple_weights(kind))
                ]
            )

        def inequalities(point: np.ndarray) -> np.ndarray:
            residuals = []
            for stage in range(args.stages):
                current_f, current_g = unpack(point, stage)
                if stage + 1 == args.stages:
                    following_f = np.full(length, -math.inf)
                    following_f[target_weight - 1] = 0.0
                    following_g = np.zeros(length)
                else:
                    following_f, following_g = unpack(point, stage + 1)
                transformed_f = transform(following_f, log_kernel[1:, 1:])
                transformed_g = transform(following_g, log_kernel[1:, 1:])
                for index in constraint_stage_rows[stage]:
                    first, second, difference = triple_weights(types[index])
                    left = (
                        current_f[first - 1]
                        + current_f[second - 1]
                        + current_g[difference - 1]
                    )
                    right = (
                        transformed_f[first - 1]
                        + transformed_f[second - 1]
                        + transformed_g[difference - 1]
                        + penalties[index]
                    )
                    if not math.isfinite(float(right)):
                        continue
                    residuals.append(left - right)
            return np.asarray(residuals)

        bounds = []
        for support in supports:
            for index in range(2 * len(support)):
                bounds.append((0.0, 0.0) if index == 0 else (None, None))
        greedy_objective = objective(initial)
        initial_slack = inequalities(initial)
        if not np.isfinite(initial).all() or not np.isfinite(initial_slack).all():
            raise RuntimeError(
                f"nonfinite greedy point at weight {target_weight}: "
                f"point={np.isfinite(initial).sum()}/{initial.size}, "
                f"slack={np.isfinite(initial_slack).sum()}/{initial_slack.size}"
            )
        if float(np.min(initial_slack)) < -1e-7:
            raise RuntimeError(
                f"greedy point is infeasible at weight {target_weight}: "
                f"minimum slack {float(np.min(initial_slack))}"
            )
        result = minimize(
            objective,
            initial,
            method="SLSQP",
            bounds=bounds,
            constraints={"type": "ineq", "fun": inequalities},
            options={"ftol": 1e-11, "maxiter": 3000, "disp": False},
        )
        if not result.success:
            raise RuntimeError(f"global optimization failed at weight {target_weight}: {result.message}")
        minimum_slack = float(np.min(inequalities(result.x)))
        marginal = float(distributions[args.stages][target_weight])
        mean = messages * marginal
        off_diagonal = messages * messages * math.exp(float(result.fun))
        variance_upper = max(0.0, mean + off_diagonal - mean * mean)
        row = {
            "target_weight": target_weight,
            "mean": mean,
            "greedy_objective_log2": greedy_objective / math.log(2.0),
            "global_objective_log2": float(result.fun) / math.log(2.0),
            "global_variance_over_mean_upper": variance_upper / mean if mean else None,
            "minimum_constraint_slack": minimum_slack,
            "iterations": int(result.nit),
        }
        rows.append(row)
        print(
            f"weight,{target_weight},greedy,{row['greedy_objective_log2']:.9f},"
            f"global,{row['global_objective_log2']:.9f},"
            f"variance_ratio,{row['global_variance_over_mean_upper']:.9f}",
            flush=True,
        )

    worst = max(rows, key=lambda row: float(row["global_variance_over_mean_upper"] or 0.0))
    payload = {
        "schema": "pair-projection-global-young-small-probe-v1",
        "status": "EXACT_KERNEL_MAJORANT_WITH_BINARY64_GLOBAL_CONVEX_OPTIMIZATION",
        "parameters": {
            "base_code": "RM(1,3) [8,4,4]",
            "length": length,
            "accumulator_stages": args.stages,
            "averaging_bound": "Fourier convolution E[f(U)f(V)g(U+V)] <= E[f^2]E[g]",
        },
        "worst": worst,
        "rows": rows,
        "limitations": [
            "The pair-projection kernel majorant is exact.",
            "SLSQP and the reported objectives use nearest binary64 arithmetic.",
            "Convexity does not make an unaudited SLSQP result an outward certificate.",
            "The length-eight result does not prove the length-512 target.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"worst": worst}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
