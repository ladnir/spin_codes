#!/usr/bin/env python3
"""Test a triangle-Holder second-moment bound on the length-eight model.

For independent messages U and V, the codewords U, V, and U+V are pairwise
independent and have the same ordinary weight distribution.  The program
majorizes the exact pair-kernel shell probability by

    f(wt(U)) f(wt(V)) g(wt(U+V))

and applies a triangle Holder inequality.  A convex optimization finds the
best such factorization on the ordinary-spectrum support.  This is a small
proof-of-concept for avoiding a genus-two weight enumerator.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import argparse
from collections import Counter

import numpy as np
from scipy.optimize import LinearConstraint, minimize
from scipy.special import logsumexp

from analyze_accumulator_pair_chain_small import compositions4, rm13_words
from verify_accumulator_pair_type_kernel import kernel_counts, multinomial, pair_type


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "triangle_holder_pair_bound_B8_probe.json"
LOCAL_BITS = 8


def triple_weights(kind: tuple[int, int, int, int]) -> tuple[int, int, int]:
    # Symbol order is 00, 01, 10, 11.
    return kind[2] + kind[3], kind[1] + kind[3], kind[1] + kind[2]


def complete_chain(types: list[tuple[int, int, int, int]], block_bits: int) -> np.ndarray:
    index = {kind: position for position, kind in enumerate(types)}
    result = np.zeros((len(types), len(types)))
    for row, input_type in enumerate(types):
        denominator = multinomial(input_type)
        for output_type, count in kernel_counts(input_type).items():
            result[row, index[output_type]] = count / denominator
    if not np.allclose(result.sum(axis=1), 1.0, atol=2e-15):
        raise AssertionError("pair transition is not stochastic")
    return result


def direct_sum_pair_distribution(
    types: list[tuple[int, int, int, int]], blocks: int
) -> np.ndarray:
    index = {kind: position for position, kind in enumerate(types)}
    local_words = rm13_words()
    local_counts: Counter[tuple[int, int, int, int]] = Counter(
        pair_type(first, second, LOCAL_BITS)
        for first in local_words
        for second in local_words
    )
    counts: Counter[tuple[int, int, int, int]] = Counter({(0, 0, 0, 0): 1})
    for _ in range(blocks):
        following: Counter[tuple[int, int, int, int]] = Counter()
        for left, left_count in counts.items():
            for right, right_count in local_counts.items():
                combined = tuple(a + b for a, b in zip(left, right))
                following[combined] += left_count * right_count
        counts = following
    result = np.zeros(len(types))
    for kind, count in counts.items():
        result[index[kind]] = count
    message_count = 1 << (4 * blocks)
    if int(result.sum()) != message_count * message_count:
        raise AssertionError("direct-sum pair distribution has the wrong mass")
    return result / (message_count * message_count)


def direct_sum_spectrum(blocks: int) -> np.ndarray:
    local_words = rm13_words()
    local = np.bincount(
        [word.bit_count() for word in local_words], minlength=LOCAL_BITS + 1
    ).astype(object)
    result = np.asarray([1], dtype=object)
    for _ in range(blocks):
        result = np.convolve(result, local)
    return result


def one_word_transition(block_bits: int) -> np.ndarray:
    result = np.zeros((block_bits + 1, block_bits + 1))
    result[0, 0] = 1.0
    for input_weight in range(1, block_bits + 1):
        denominator = math.comb(block_bits, input_weight)
        half_down = input_weight // 2
        half_up = (input_weight + 1) // 2
        for output_weight in range(1, block_bits + 1):
            if half_down > block_bits - output_weight or half_up > output_weight:
                continue
            result[input_weight, output_weight] = (
                math.comb(block_bits - output_weight, half_down)
                * math.comb(output_weight - 1, half_up - 1)
                / denominator
            )
    if not np.allclose(result.sum(axis=1), 1.0, atol=2e-15):
        raise AssertionError("one-word transition is not stochastic")
    return result


def optimize_factorization(
    bad_probability: np.ndarray,
    types: list[tuple[int, int, int, int]],
    support: list[int],
    probabilities: np.ndarray,
    averaging_bound: str,
) -> dict[str, object] | None:
    support_index = {weight: index for index, weight in enumerate(support)}
    coefficient_rows = []
    lower_bounds = []
    constraint_types = []
    for probability, kind in zip(bad_probability, types):
        if probability <= 0:
            continue
        first, second, difference = triple_weights(kind)
        if first not in support_index or second not in support_index or difference not in support_index:
            continue
        row = np.zeros(2 * len(support))
        row[support_index[first]] += 1.0
        row[support_index[second]] += 1.0
        row[len(support) + support_index[difference]] += 1.0
        coefficient_rows.append(row)
        lower_bounds.append(math.log(probability))
        constraint_types.append((first, second, difference))
    if not coefficient_rows:
        return None
    coefficients = np.asarray(coefficient_rows)
    lower = np.asarray(lower_bounds)

    def objective(values: np.ndarray) -> float:
        first = logsumexp(np.log(probabilities) + 2.0 * values[: len(support)])
        if averaging_bound == "holder":
            second = logsumexp(
                np.log(probabilities) + 2.0 * values[len(support) :]
            )
            return float(first + 0.5 * second)
        second = logsumexp(np.log(probabilities) + values[len(support) :])
        return float(first + second)

    # Fix one gauge coordinate. The objective and constraints are invariant
    # under f <- exp(c) f and g <- exp(-2c) g.
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
        options={"ftol": 1e-12, "maxiter": 2000},
    )
    if not result.success:
        raise RuntimeError(f"factor optimization failed: {result.message}")
    slack = coefficients @ result.x - lower
    active = [
        {
            "triple_weights": list(kind),
            "slack": float(value),
        }
        for kind, value in zip(constraint_types, slack)
        if value <= 1e-7
    ]
    return {
        "objective_log2": objective(result.x) / math.log(2.0),
        "minimum_constraint_slack": float(slack.min()),
        "f_log2": [float(value / math.log(2.0)) for value in result.x[: len(support)]],
        "g_log2": [float(value / math.log(2.0)) for value in result.x[len(support) :]],
        "constraints": len(lower),
        "active_constraints": active,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blocks", type=int, choices=(1, 2, 3), default=1)
    parser.add_argument(
        "--averaging-bound", choices=("holder", "young"), default="holder"
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    block_bits = LOCAL_BITS * args.blocks
    message_count = 1 << (4 * args.blocks)
    types = compositions4(block_bits)
    transition = complete_chain(types, block_bits)
    initial = direct_sum_pair_distribution(types, args.blocks)
    spectrum = direct_sum_spectrum(args.blocks)
    # The zero message and diagonal pairs are handled exactly. Setting both
    # factor functions to zero at weight zero makes the product vanish when
    # U=0, V=0, or U=V. Only positive ordinary-spectrum weights remain.
    support = [weight for weight, count in enumerate(spectrum) if weight and count]
    probabilities = np.asarray([spectrum[weight] / message_count for weight in support])

    rows = []
    terminal = np.zeros((len(types), block_bits + 1))
    for index, kind in enumerate(types):
        first, second, _difference = triple_weights(kind)
        if first == second:
            terminal[index, first] = 1.0
    backward = terminal
    word_distribution = np.asarray(spectrum, dtype=float) / message_count
    word_transition = one_word_transition(block_bits)
    for stages in range(1, 4):
        backward = transition @ backward
        word_distribution = word_distribution @ word_transition
        for target_weight in range(1, block_bits + 1):
            bad_probability = backward[:, target_weight]
            if not np.any(bad_probability):
                continue
            factor = optimize_factorization(
                bad_probability,
                types,
                support,
                probabilities,
                args.averaging_bound,
            )
            if factor is None:
                continue
            marginal = float(word_distribution[target_weight])
            joint = float(initial @ bad_probability)
            mean = message_count * marginal
            exact_second = message_count * message_count * joint
            exact_variance = max(0.0, exact_second - mean * mean)
            off_diagonal_upper = message_count * message_count * 2.0 ** float(
                factor["objective_log2"]
            )
            second_upper = mean + off_diagonal_upper
            variance_upper = max(0.0, second_upper - mean * mean)
            rows.append(
                {
                    "accumulator_stages": stages,
                    "target_weight": target_weight,
                    "mean": mean,
                    "exact_variance": exact_variance,
                    "exact_variance_over_mean": exact_variance / mean if mean else None,
                    "holder_second_moment_upper": second_upper,
                    "holder_off_diagonal_second_moment_upper": off_diagonal_upper,
                    "holder_variance_upper": variance_upper,
                    "holder_variance_over_mean_upper": variance_upper / mean if mean else None,
                    "factorization": factor,
                }
            )
        finite = [
            row
            for row in rows
            if row["accumulator_stages"] == stages
            and row["holder_variance_over_mean_upper"] is not None
        ]
        worst = max(finite, key=lambda row: float(row["holder_variance_over_mean_upper"]))
        print(
            f"stage,{stages},worst_weight,{worst['target_weight']},"
            f"holder_var_ratio,{worst['holder_variance_over_mean_upper']:.9g}",
            flush=True,
        )

    result = {
        "schema": "triangle-pair-bound-small-probe-v3",
        "status": "EXACT_PAIR_KERNEL_WITH_BINARY64_CONVEX_FACTORIZATION",
        "parameters": {
            "base_code": f"direct sum of {args.blocks} RM(1,3) [8,4,4] blocks",
            "block_bits": block_bits,
            "dimension": 4 * args.blocks,
            "ordinary_spectrum": {
                str(weight): int(count)
                for weight, count in enumerate(spectrum)
                if count
            },
            "factor_support": support,
            "degenerate_pair_handling": "U=0, V=0, and U=V are excluded by f(0)=g(0)=0; the diagonal second-moment term is added exactly as the mean",
            "accumulator_stages": [1, 2, 3],
            "averaging_bound": args.averaging_bound,
        },
        "inequality": (
            "E[f(U) f(V) g(U+V)] <= E[f(U)^2] E[g(U)] by Fourier convolution"
            if args.averaging_bound == "young"
            else "E[f(U) f(V) g(U+V)] <= E[f(U)^2] sqrt(E[g(U)^2]) by Cauchy--Schwarz"
        ),
        "rows": rows,
        "limitations": [
            "The factorization optimization uses nearest binary64 arithmetic.",
            "The small-block result does not prove a length-512 factorization.",
            "A length-512 implementation still needs a non-materialized way to bound the three-stage pair kernel.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
