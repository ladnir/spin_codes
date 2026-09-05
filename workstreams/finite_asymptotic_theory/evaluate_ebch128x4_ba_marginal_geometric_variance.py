#!/usr/bin/env python3
"""Bound BA shell variances through three one-word marginal events.

If both transformed words have weight w, their transformed difference has
weight at most 2w.  For a fixed rank-two input pair with initial weights
(h1,h2,h3), the joint probability is therefore at most the geometric mean
of the three corresponding one-word probabilities.  A Triangle--Holder
inequality then averages this product using only the ordinary base spectrum.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from evaluate_ebch128x4_ba_expected_spectrum import (
    B,
    K,
    SOURCE,
    direct_sum_spectrum,
    read_wd,
)
from probe_triangle_holder_pair_bound_small import one_word_transition


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch128x4_ba_marginal_geometric_variance_probe.json"
LN2 = math.log(2.0)


def logdiffexp(log_large: float, log_small: float) -> float:
    if log_small == -math.inf:
        return log_large
    if log_small > log_large + 1e-10:
        raise ArithmeticError("negative variance after geometric bound")
    if log_small >= log_large:
        return -math.inf
    return log_large + math.log1p(-math.exp(log_small - log_large))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-stages", type=int, default=16)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    base_log2_counts = direct_sum_spectrum(read_wd(SOURCE))
    base_log_probabilities = base_log2_counts * LN2 - K * LN2
    base_log_probabilities[0] = -math.inf
    transition = one_word_transition(B)
    power = np.eye(B + 1)
    rows = []
    stage_summaries = []
    for stages in range(1, args.maximum_stages + 1):
        power = power @ transition
        marginal_distribution = np.exp(base_log_probabilities) @ power
        marginal_distribution[0] = 0.0
        stage_rows = []
        for weight in range(1, B + 1):
            one_shell = power[:, weight].copy()
            difference_limit = 2 * min(weight, B - weight)
            difference_tail = power[:, 1 : difference_limit + 1].sum(axis=1)
            one_shell[0] = 0.0
            difference_tail[0] = 0.0
            positive_shell = one_shell > 0.0
            positive_tail = difference_tail > 0.0
            if not positive_shell.any() or marginal_distribution[weight] == 0.0:
                continue

            log_one_shell = np.full(B + 1, -math.inf)
            log_one_shell[positive_shell] = np.log(one_shell[positive_shell])
            log_difference_tail = np.full(B + 1, -math.inf)
            log_difference_tail[positive_tail] = np.log(difference_tail[positive_tail])

            def factor_logs(one_word_exponent: float) -> tuple[float, float]:
                difference_exponent = 1.0 - 2.0 * one_word_exponent
                if one_word_exponent == 0.0:
                    shell_terms = base_log_probabilities
                else:
                    shell_terms = (
                        base_log_probabilities
                        + 2.0 * one_word_exponent * log_one_shell
                    )
                if difference_exponent == 0.0:
                    tail_terms = base_log_probabilities
                else:
                    tail_terms = (
                        base_log_probabilities
                        + 2.0 * difference_exponent * log_difference_tail
                    )
                return float(logsumexp(shell_terms)), float(logsumexp(tail_terms))

            def factor_objective(one_word_exponent: float) -> float:
                shell_value, tail_value = factor_logs(one_word_exponent)
                return shell_value + 0.5 * tail_value

            candidates = [(0.0, factor_objective(0.0)), (0.5, factor_objective(0.5))]
            optimized = minimize_scalar(
                factor_objective,
                bounds=(0.0, 0.5),
                method="bounded",
                options={"xatol": 1e-12},
            )
            candidates.append((float(optimized.x), float(optimized.fun)))
            one_word_exponent, _ = min(candidates, key=lambda item: item[1])
            difference_exponent = 1.0 - 2.0 * one_word_exponent
            log_shell_factor, log_tail_factor = factor_logs(one_word_exponent)
            log_off_diagonal = (
                2.0 * K * LN2 + log_shell_factor + 0.5 * log_tail_factor
            )
            log_mean = K * LN2 + math.log(float(marginal_distribution[weight]))
            log_second_upper = float(np.logaddexp(log_mean, log_off_diagonal))
            log_mean_square = 2.0 * log_mean
            log_variance_upper = logdiffexp(log_second_upper, log_mean_square)
            log_ratio = log_variance_upper - log_mean
            row = {
                "accumulator_stages": stages,
                "target_weight": weight,
                "mean_log2": log_mean / LN2,
                "variance_over_mean_log2_upper": log_ratio / LN2,
                "variance_over_mean_upper": math.exp(log_ratio) if log_ratio < 700 else math.inf,
                "shell_factor_log2": log_shell_factor / LN2,
                "difference_tail_factor_log2": log_tail_factor / LN2,
                "difference_weight_upper": difference_limit,
                "one_word_intersection_exponent": one_word_exponent,
                "difference_intersection_exponent": difference_exponent,
            }
            rows.append(row)
            stage_rows.append(row)

        worst = max(stage_rows, key=lambda row: float(row["variance_over_mean_log2_upper"]))
        below_512 = [
            row for row in stage_rows if float(row["variance_over_mean_log2_upper"]) <= 9.0
        ]
        summary = {
            "accumulator_stages": stages,
            "worst_weight": worst["target_weight"],
            "worst_variance_over_mean_log2_upper": worst["variance_over_mean_log2_upper"],
            "shells_at_most_512": len(below_512),
            "active_shells": len(stage_rows),
        }
        stage_summaries.append(summary)
        print(
            f"stage,{stages},worst_weight,{worst['target_weight']},"
            f"worst_log2_ratio,{worst['variance_over_mean_log2_upper']:.9f},"
            f"shells_le_512,{len(below_512)},{len(stage_rows)}",
            flush=True,
        )

    result = {
        "schema": "ebch128x4-ba-marginal-geometric-variance-probe-v1",
        "status": "EXACT_REDUCTION_WITH_BINARY64_ONE_WORD_CHAIN",
        "construction": {
            "base": "direct sum of four fixed extended BCH [128,64,22] codes",
            "maximum_accumulator_stages": args.maximum_stages,
            "sampler": "one independent uniform 512-coordinate permutation before each accumulator",
        },
        "bound": {
            "pointwise": "Phi_t,w(h1,h2,h3) <= q_t,w(h1)^a q_t,w(h2)^a q_t,<=2 min(w,B-w)(h3)^(1-2a), optimized over 0<=a<=1/2",
            "averaging": "E[f(H(U))f(H(V))g(H(U+V))] <= E[f(H)^2] sqrt(E[g(H)^2])",
            "degenerate_pairs": "f(0)=g(0)=0; the diagonal contribution is added as the shell mean",
        },
        "stage_summaries": stage_summaries,
        "rows": rows,
        "limitations": [
            "The geometric intersection inequality and Triangle--Holder reduction are exact.",
            "The one-word transition powers and reported logarithms use nearest binary64 arithmetic.",
            "A passing diagnostic would still require outward arithmetic and the complete conditional transfer.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
