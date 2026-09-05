#!/usr/bin/env python3
"""Combine the scalable pair projection with genus-two MacWilliams constraints.

The character-projection envelope is a pointwise upper bound on the common-
interleaver pair kernel.  Averaging that envelope from ordinary marginals was
far too loose.  This probe instead maximizes it over pair enumerators whose
dual MacWilliams transform is nonnegative and has the same ordinary spectrum.

All combinatorial kernels are exact before conversion to binary64.  The
factorizations and LP optima are diagnostics, not outward certificates.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_accumulator_pair_chain_small import compositions4, direct_sum_words, rank_two
from probe_genus2_macwilliams_lp_small import (
    macwilliams_matrix,
    marginal_rows,
    pair_counts,
    pair_transition,
    rank_one_bounds,
    solve_objective,
)
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
OUTPUT = WORKSTREAM / "projection_macwilliams_lp_B16_probe.json"


def propagate(log_values: np.ndarray, log_kernel: np.ndarray) -> np.ndarray:
    return logsumexp(log_kernel + log_values[None, :], axis=1)


def projection_envelope(
    *,
    length: int,
    stages: int,
    target_weight: int,
    rank_two_types: list[tuple[int, int, int, int]],
    overlap_penalty: np.ndarray,
    distributions: list[np.ndarray],
    log_binary_power: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, object]]]:
    log_f = np.full(length, -math.inf)
    log_f[target_weight - 1] = 0.0
    log_g = np.zeros(length)
    steps = []
    for backward_step in range(1, stages + 1):
        transformed_f = propagate(log_f, log_binary_power[1:, 1:])
        transformed_g = propagate(log_g, log_binary_power[1:, 1:])
        lower = np.asarray(
            [
                transformed_f[triple_weights(kind)[0] - 1]
                + transformed_f[triple_weights(kind)[1] - 1]
                + transformed_g[triple_weights(kind)[2] - 1]
                + overlap_penalty[index]
                for index, kind in enumerate(rank_two_types)
            ]
        )
        input_stage = stages - backward_step
        result = factorize(
            lower,
            rank_two_types,
            distributions[input_stage],
            "young",
        )
        if result is None:
            raise RuntimeError("projection factorization had empty support")
        log_f = np.asarray(result["f_log"])
        log_g = np.asarray(result["g_log"])
        steps.append(
            {
                "input_stage": input_stage,
                "objective_log2": result["objective_log2"],
                "minimum_constraint_slack": result["minimum_constraint_slack"],
                "constraints": result["constraints"],
            }
        )
    return log_f, log_g, steps


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blocks", type=int, choices=(1, 2), default=2)
    parser.add_argument("--maximum-stages", type=int, default=3)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    length = 8 * args.blocks
    dimension = 4 * args.blocks
    messages = 1 << dimension
    spectrum = np.asarray(direct_sum_spectrum(args.blocks), dtype=float)
    types = compositions4(length)
    rank_two_types = [kind for kind in types if rank_two(kind)]
    rank_two_index = {kind: index for index, kind in enumerate(rank_two_types)}
    words = direct_sum_words(args.blocks)
    actual_counts = pair_counts(words, length, types)

    print(f"building_macwilliams length={length} types={len(types)}", flush=True)
    transform = macwilliams_matrix(types)
    marginal_matrix, marginal_rhs = marginal_rows(
        types, [int(value) for value in spectrum], messages
    )
    bounds = rank_one_bounds(types, [int(value) for value in spectrum])
    transition = pair_transition(types)
    binary = one_word_transition(length)
    log_binary_power = np.full_like(binary, -math.inf)
    positive = binary > 0
    log_binary_power[positive] = (2.0 / 3.0) * np.log(binary[positive])
    overlap_penalty = np.asarray(
        [log_inverse_geometric_overlap(kind, length) for kind in rank_two_types]
    )

    distributions = [spectrum / messages]
    for _ in range(args.maximum_stages):
        distributions.append(distributions[-1] @ binary)
    pair_power = np.eye(len(types))
    rows = []
    for stages in range(1, args.maximum_stages + 1):
        pair_power = pair_power @ transition
        for target_weight in range(1, length + 1):
            mean = messages * float(distributions[stages][target_weight])
            if mean <= 1e-14:
                continue
            terminal = np.asarray(
                [
                    float(
                        kind[2] + kind[3] == target_weight
                        and kind[1] + kind[3] == target_weight
                    )
                    for kind in types
                ]
            )
            exact_objective = pair_power @ terminal
            log_f, log_g, steps = projection_envelope(
                length=length,
                stages=stages,
                target_weight=target_weight,
                rank_two_types=rank_two_types,
                overlap_penalty=overlap_penalty,
                distributions=distributions,
                log_binary_power=log_binary_power,
            )
            envelope = exact_objective.copy()
            for index, kind in enumerate(types):
                if kind not in rank_two_index:
                    continue
                first, second, difference = triple_weights(kind)
                value = log_f[first - 1] + log_f[second - 1] + log_g[difference - 1]
                envelope[index] = min(1.0, math.exp(value)) if math.isfinite(value) else 0.0
            admissible = np.asarray(
                [
                    all(spectrum[weight] > 0 for weight in triple_weights(kind))
                    for kind in types
                ]
            )
            pointwise_defect = float(
                np.max((exact_objective - envelope)[admissible])
            )
            actual_second = float(actual_counts @ exact_objective)
            actual_ratio = max(0.0, actual_second - mean * mean) / mean if mean else None
            direct_second = float(actual_counts @ envelope)
            direct_ratio = max(0.0, direct_second - mean * mean) / mean if mean else None
            formal_second, formal_status = solve_objective(
                envelope,
                transform,
                messages,
                marginal_matrix,
                marginal_rhs,
                bounds,
                "formal-dual",
            )
            formal_ratio = (
                max(0.0, formal_second - mean * mean) / mean
                if mean and math.isfinite(formal_second)
                else None
            )
            row = {
                "accumulator_stages": stages,
                "target_weight": target_weight,
                "mean": mean,
                "admissible_support_maximum_exact_minus_envelope": pointwise_defect,
                "actual_exact_variance_over_mean": actual_ratio,
                "actual_projection_variance_over_mean_upper": direct_ratio,
                "formal_dual_projection_variance_over_mean_upper": formal_ratio,
                "formal_dual_status": formal_status,
                "factorization_steps": steps,
            }
            rows.append(row)
            print(
                f"stage,{stages},weight,{target_weight},actual,{actual_ratio},"
                f"actual_projection,{direct_ratio},formal_projection,{formal_ratio},"
                f"defect,{pointwise_defect}",
                flush=True,
            )

    fields = (
        "actual_exact_variance_over_mean",
        "actual_projection_variance_over_mean_upper",
        "formal_dual_projection_variance_over_mean_upper",
    )
    worst = {
        field: max(
            (row for row in rows if row[field] is not None),
            key=lambda row: float(row[field]),
        )
        for field in fields
    }
    payload = {
        "schema": "projection-macwilliams-lp-small-probe-v1",
        "status": "EXACT_KERNEL_MAJORANTS_WITH_BINARY64_FACTORIZATION_AND_LP",
        "parameters": {
            "base_code": f"direct sum of {args.blocks} RM(1,3) [8,4,4] blocks",
            "length": length,
            "dimension": dimension,
            "pair_types": len(types),
            "rank_two_pair_types": len(rank_two_types),
        },
        "worst": worst,
        "rows": rows,
        "limitations": [
            "The character-projection kernel inequality is exact.",
            "The factor witnesses and LP optima use nearest binary64 arithmetic.",
            "The LP constrains an unknown genus-two enumerator but does not enumerate message pairs.",
            "The small-length result does not prove the EBCH128 or length-512 target.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"worst": worst}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
