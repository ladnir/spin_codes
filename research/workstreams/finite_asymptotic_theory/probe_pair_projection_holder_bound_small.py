#!/usr/bin/env python3
"""Probe a scalable pair-kernel bound using its three binary projections.

For a rank-two pair type n with three nonzero-character weights h_i, let
r_ij(n) be the hypergeometric probability of that pair type under independent
binary interleavers of components i and j. The exact pair kernel P and the
binary accumulator weight kernel Q satisfy

    P(n,m) <= Q(h_1,w_1) Q(h_2,w_2) / r_12(n),

and the two analogous inequalities. Their geometric mean gives a product of
three one-dimensional kernels. This program propagates a factored envelope
through that bound and tests its loss on direct sums of RM(1,3).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, minimize
from scipy.special import gammaln, logsumexp

from analyze_accumulator_pair_chain_small import compositions4, rank_two
from probe_triangle_holder_pair_bound_small import (
    direct_sum_spectrum,
    one_word_transition,
    triple_weights,
)


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "pair_projection_holder_B8_probe.json"


def log_choose(total: int, selected: int) -> float:
    return float(
        gammaln(total + 1)
        - gammaln(selected + 1)
        - gammaln(total - selected + 1)
    )


def log_inverse_geometric_overlap(
    kind: tuple[int, int, int, int], block_bits: int
) -> float:
    h1, h2, h3 = triple_weights(kind)
    log_multinomial = float(gammaln(block_bits + 1)) - sum(
        float(gammaln(value + 1)) for value in kind
    )
    return (
        2.0
        / 3.0
        * (
            log_choose(block_bits, h1)
            + log_choose(block_bits, h2)
            + log_choose(block_bits, h3)
        )
        - log_multinomial
    )


def factorize(
    lower: np.ndarray,
    types: list[tuple[int, int, int, int]],
    probabilities: np.ndarray,
    averaging_bound: str,
) -> dict[str, object] | None:
    block_bits = len(probabilities) - 1
    support = [weight for weight in range(1, block_bits + 1) if probabilities[weight] > 0]
    support_index = {weight: index for index, weight in enumerate(support)}
    coefficient_rows = []
    retained_lower = []
    for bound, kind in zip(lower, types):
        first, second, difference = triple_weights(kind)
        if not math.isfinite(float(bound)):
            continue
        if any(weight not in support_index for weight in (first, second, difference)):
            continue
        row = np.zeros(2 * len(support))
        row[support_index[first]] += 1.0
        row[support_index[second]] += 1.0
        row[len(support) + support_index[difference]] += 1.0
        coefficient_rows.append(row)
        retained_lower.append(bound)
    coefficients = np.asarray(coefficient_rows)
    lower = np.asarray(retained_lower)
    if not retained_lower:
        return None

    log_probabilities = np.log(np.asarray([probabilities[w] for w in support]))

    def objective(values: np.ndarray) -> float:
        first = logsumexp(log_probabilities + 2.0 * values[: len(support)])
        if averaging_bound == "holder":
            second = logsumexp(
                log_probabilities + 2.0 * values[len(support) :]
            )
            return float(first + 0.5 * second)
        second = logsumexp(log_probabilities + values[len(support) :])
        return float(first + second)

    gauge = np.zeros((1, 2 * len(support)))
    gauge[0, 0] = 1.0
    constraints = [
        LinearConstraint(coefficients, lower, np.full_like(lower, np.inf)),
        LinearConstraint(gauge, np.zeros(1), np.zeros(1)),
    ]
    result = minimize(
        objective,
        np.zeros(2 * len(support)),
        method="SLSQP",
        constraints=constraints,
        options={"ftol": 1e-11, "maxiter": 3000},
    )
    if not result.success:
        raise RuntimeError(f"factor optimization failed: {result.message}")
    slack = coefficients @ result.x - lower
    full_f = np.full(block_bits, -math.inf)
    full_g = np.full(block_bits, -math.inf)
    for weight, index in support_index.items():
        full_f[weight - 1] = result.x[index]
        full_g[weight - 1] = result.x[len(support) + index]
    return {
        "objective_log2": objective(result.x) / math.log(2.0),
        "minimum_constraint_slack": float(slack.min()),
        "f_log": full_f,
        "g_log": full_g,
        "factor_support": support,
        "constraints": len(lower),
    }


def propagate_log(
    log_f: np.ndarray, log_kernel: np.ndarray
) -> np.ndarray:
    return logsumexp(log_kernel + log_f[None, :], axis=1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blocks", type=int, choices=(1, 2, 3), default=1)
    parser.add_argument("--maximum-stages", type=int, default=3)
    parser.add_argument(
        "--averaging-bound", choices=("holder", "young"), default="holder"
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    block_bits = 8 * args.blocks
    message_count = 1 << (4 * args.blocks)
    spectrum = direct_sum_spectrum(args.blocks)
    base = np.asarray(spectrum, dtype=float) / message_count
    binary = one_word_transition(block_bits)
    log_binary_power = np.full_like(binary, -math.inf)
    positive = binary > 0
    log_binary_power[positive] = (2.0 / 3.0) * np.log(binary[positive])

    types = [kind for kind in compositions4(block_bits) if rank_two(kind)]
    overlap_penalty = np.asarray(
        [log_inverse_geometric_overlap(kind, block_bits) for kind in types]
    )
    distributions = [base]
    for _ in range(args.maximum_stages):
        distributions.append(distributions[-1] @ binary)

    rows = []
    for stages in range(1, args.maximum_stages + 1):
        for target_weight in range(1, block_bits + 1):
            log_f = np.full(block_bits, -math.inf)
            log_f[target_weight - 1] = 0.0
            log_g = np.zeros(block_bits)
            steps = []
            valid = True
            for backward_step in range(1, stages + 1):
                transformed_f = propagate_log(log_f, log_binary_power[1:, 1:])
                transformed_g = propagate_log(log_g, log_binary_power[1:, 1:])
                lower = np.empty(len(types))
                for index, kind in enumerate(types):
                    h1, h2, h3 = triple_weights(kind)
                    lower[index] = (
                        transformed_f[h1 - 1]
                        + transformed_f[h2 - 1]
                        + transformed_g[h3 - 1]
                        + overlap_penalty[index]
                    )
                input_stage = stages - backward_step
                factor = factorize(
                    lower, types, distributions[input_stage], args.averaging_bound
                )
                if factor is None:
                    valid = False
                    break
                log_f = np.asarray(factor.pop("f_log"))
                log_g = np.asarray(factor.pop("g_log"))
                steps.append(factor)

            if not valid:
                continue
            log_prob = np.full(block_bits, -math.inf)
            positive_base = base[1:] > 0
            log_prob[positive_base] = np.log(base[1:][positive_base])
            if args.averaging_bound == "holder":
                objective = float(
                    logsumexp(log_prob + 2.0 * log_f)
                    + 0.5 * logsumexp(log_prob + 2.0 * log_g)
                )
            else:
                objective = float(
                    logsumexp(log_prob + 2.0 * log_f)
                    + logsumexp(log_prob + log_g)
                )
            mean = message_count * float(distributions[stages][target_weight])
            off_diagonal = message_count * message_count * math.exp(objective)
            variance_upper = max(0.0, mean + off_diagonal - mean * mean)
            rows.append(
                {
                    "accumulator_stages": stages,
                    "target_weight": target_weight,
                    "mean": mean,
                    "variance_over_mean_upper": variance_upper / mean if mean else None,
                    "factor_objective_log2": objective / math.log(2.0),
                    "steps": steps,
                }
            )
        active = [
            row
            for row in rows
            if row["accumulator_stages"] == stages
            and row["variance_over_mean_upper"] is not None
        ]
        worst = max(active, key=lambda row: row["variance_over_mean_upper"])
        print(
            f"stage,{stages},worst_weight,{worst['target_weight']},"
            f"variance_ratio,{worst['variance_over_mean_upper']:.9g}",
            flush=True,
        )

    result = {
        "schema": "pair-projection-holder-small-probe-v1",
        "status": "BINARY64_CONVEX_DIAGNOSTIC",
        "parameters": {
            "base_code": f"direct sum of {args.blocks} RM(1,3) [8,4,4] blocks",
            "block_bits": block_bits,
            "dimension": 4 * args.blocks,
            "maximum_stages": args.maximum_stages,
            "averaging_bound": args.averaging_bound,
        },
        "kernel_majorant": "P(n,m) <= product_i Q(h_i,w_i)^(2/3) / (r_12(n) r_13(n) r_23(n))^(1/3)",
        "averaging_inequality": (
            "E[f(U)f(V)g(U+V)] <= E[f(U)^2] E[g(U)]"
            if args.averaging_bound == "young"
            else "E[f(U)f(V)g(U+V)] <= E[f(U)^2] sqrt(E[g(U)^2])"
        ),
        "rows": rows,
        "limitations": [
            "The geometric-mean kernel inequality is exact, but the factor witnesses use nearest binary64 optimization.",
            "The factorization is optimized greedily at each stage rather than globally.",
            "A small-block result does not prove the length-512 variance target.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
