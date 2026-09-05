#!/usr/bin/env python3
"""Diagnose the joint Golay--BA-3 and RM2Sub exponent at distance 0.11.

This program keeps the exact finite ensemble-expected BA spectrum in log
space.  For one active-block weight class, it conditions an iid Bernoulli
row on that weight and uses the three-state RM2Sub envelope.  The resulting
optimization is a diagnostic for the transfer-weighted outer theorem.  It
does not yet cover mixtures of several active-block weight classes and does
not use outward arithmetic.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import gammaln, logsumexp

from analyze_rm2sub_dense_occupation import (
    DEFAULT_SELECTION,
    DenseEnvelope,
    association_matrix,
    state_code_weights,
    verify_transpose_columns,
)


GOLAY_LENGTH = 24
GOLAY_DIMENSION = 12
GOLAY_SPECTRUM = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}


def log_binomial(length: int, weight: np.ndarray | int) -> np.ndarray:
    values = np.asarray(weight)
    answer = np.full(values.shape, -math.inf, dtype=np.float64)
    valid = (0 <= values) & (values <= length)
    answer[valid] = (
        gammaln(length + 1)
        - gammaln(values[valid] + 1)
        - gammaln(length - values[valid] + 1)
    )
    return answer


def log_binomial_vector(lengths: np.ndarray, weight: int) -> np.ndarray:
    answer = np.full(lengths.shape, -math.inf, dtype=np.float64)
    valid = (lengths >= 0) & (weight >= 0) & (weight <= lengths)
    answer[valid] = (
        gammaln(lengths[valid] + 1)
        - gammaln(weight + 1)
        - gammaln(lengths[valid] - weight + 1)
    )
    return answer


def golay_direct_sum_log_probability(length: int) -> np.ndarray:
    if length % GOLAY_LENGTH:
        raise ValueError("B must be divisible by 24")
    constituent = np.zeros(GOLAY_LENGTH + 1, dtype=np.float64)
    for weight, count in GOLAY_SPECTRUM.items():
        constituent[weight] = count / (1 << GOLAY_DIMENSION)
    distribution = np.asarray([1.0], dtype=np.float64)
    for _ in range(length // GOLAY_LENGTH):
        distribution = np.convolve(distribution, constituent)
    result = np.full(length + 1, -math.inf, dtype=np.float64)
    positive = distribution > 0.0
    result[positive] = np.log(distribution[positive])
    return result


def accumulator_log_transition(length: int) -> np.ndarray:
    transition = np.full((length + 1, length + 1), -math.inf)
    transition[0, 0] = 0.0
    output_weights = np.arange(length + 1)
    for input_weight in range(1, length + 1):
        runs = (input_weight + 1) // 2
        counts = log_binomial_vector(output_weights - 1, runs - 1)
        counts += log_binomial_vector(length - output_weights, input_weight - runs)
        transition[input_weight] = counts - float(
            log_binomial(length, input_weight)
        )
    return transition


def expected_ba_log_spectrum(length: int) -> np.ndarray:
    log_probability = golay_direct_sum_log_probability(length)
    transition = accumulator_log_transition(length)
    after_one = logsumexp(log_probability[:, None] + transition, axis=0)
    after_two = logsumexp(after_one[:, None] + transition, axis=0)
    return after_two + (length // 2) * math.log(2.0)


def binary_entropy_natural(value: float) -> float:
    if value <= 0.0 or value >= 1.0:
        return 0.0
    return -value * math.log(value) - (1.0 - value) * math.log1p(-value)


def golay_exponent(alpha: float) -> float:
    if alpha <= 0.0:
        return 0.0
    if alpha >= 1.0:
        return 0.0

    def objective(log_r: float) -> float:
        terms = [
            math.log(count) + weight * log_r
            for weight, count in GOLAY_SPECTRUM.items()
        ]
        scale = max(terms)
        log_polynomial = scale + math.log(
            sum(math.exp(term - scale) for term in terms)
        )
        return log_polynomial / GOLAY_LENGTH - alpha * log_r

    result = minimize_scalar(
        objective,
        bounds=(-80.0, 80.0),
        method="bounded",
        options={"xatol": 1.0e-13},
    )
    return float(result.fun)


def accumulator_exponent(input_density: float, output_density: float) -> float:
    if input_density < 0.0 or input_density > 1.0:
        return -math.inf
    if output_density < input_density / 2.0:
        return -math.inf
    if output_density > 1.0 - input_density / 2.0:
        return -math.inf
    if input_density == 0.0:
        return 0.0 if output_density in (0.0, 1.0) else -math.inf
    left = input_density / (2.0 * output_density) if output_density else 0.0
    right = (
        input_density / (2.0 * (1.0 - output_density))
        if output_density < 1.0
        else 0.0
    )
    return (
        output_density * binary_entropy_natural(left)
        + (1.0 - output_density) * binary_entropy_natural(right)
        - binary_entropy_natural(input_density)
    )


def ba_asymptotic_exponent(output_density: float) -> float:
    beta_max = min(1.0, 2.0 * output_density, 2.0 * (1.0 - output_density))
    if beta_max <= 0.0:
        return -math.inf

    def decode(point: np.ndarray) -> tuple[float, float]:
        beta = beta_max * logistic(float(point[0]))
        alpha_max = 2.0 * min(beta, 1.0 - beta)
        alpha = alpha_max * logistic(float(point[1]))
        return alpha, beta

    def negative(point: np.ndarray) -> float:
        alpha, beta = decode(point)
        return -(
            golay_exponent(alpha)
            + accumulator_exponent(alpha, beta)
            + accumulator_exponent(beta, output_density)
        )

    starts = []
    for beta_fraction in (0.1, 0.25, 0.5, 0.75, 0.9):
        for alpha_fraction in (0.1, 0.35, 0.65, 0.9):
            starts.append(
                np.asarray([logit(beta_fraction), logit(alpha_fraction)])
            )
    best = None
    for start in starts:
        result = minimize(
            negative,
            start,
            method="Nelder-Mead",
            options={"maxiter": 800, "xatol": 1.0e-10, "fatol": 1.0e-12},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    return -float(best.fun)


def expurgated_random_hull(output_density: float, cutoff: float = 0.104) -> float:
    """Return the ideal expurgated-random concave hull used diagnostically."""
    entropy = binary_entropy_natural

    def tangent_equation(point: float) -> float:
        derivative = math.log((1.0 - point) / point)
        secant = (entropy(point) - math.log(2.0) / 2.0) / (point - cutoff)
        return derivative - secant

    tangent = brentq(tangent_equation, cutoff + 1.0e-10, 0.499999)
    slope = math.log((1.0 - tangent) / tangent)
    if output_density < tangent:
        return (output_density - cutoff) * slope
    if output_density > 1.0 - tangent:
        return (1.0 - cutoff - output_density) * slope
    return entropy(output_density) - math.log(2.0) / 2.0


def logit(probability: float) -> float:
    return math.log(probability / (1.0 - probability))


def logistic(value: float) -> float:
    epsilon = 1.0e-15
    if value >= 0.0:
        inverse = math.exp(-value)
        return min(1.0 - epsilon, 1.0 / (1.0 + inverse))
    direct = math.exp(value)
    return max(epsilon, direct / (1.0 + direct))


def kl_natural(alpha: float, probability: float) -> float:
    if alpha == 0.0:
        return -math.log1p(-probability)
    if alpha == 1.0:
        return -math.log(probability)
    return (
        alpha * math.log(alpha / probability)
        + (1.0 - alpha) * math.log((1.0 - alpha) / (1.0 - probability))
    )


def joint_objective(
    envelope: DenseEnvelope,
    *,
    alpha: float,
    row_density: float,
    outer_exponent_natural: float,
    candidate_probability: float,
    value_probability: float,
    surprisal: float,
    transfer_name: str,
) -> tuple[float, dict[str, float]]:
    actual_bit_probability = candidate_probability * value_probability
    if not 0.0 < actual_bit_probability < 1.0:
        return math.inf, {}
    transfer_function = (
        envelope.transfer
        if transfer_name == "three-state"
        else envelope.transfer_small_density
    )
    transfer, details = transfer_function(
        candidate_probability=2.0 * actual_bit_probability,
        surprisal=surprisal,
    )
    radius = float(np.max(np.abs(np.linalg.eigvals(transfer))))
    value = (
        alpha * outer_exponent_natural
        + kl_natural(alpha, candidate_probability)
        + alpha * kl_natural(row_density, value_probability)
        + math.log(radius) / envelope.step_bits
        + envelope.delta * surprisal
    )
    return value, {
        "actual_bit_probability": actual_bit_probability,
        "value_probability": value_probability,
        "perron_radius": radius,
        "z": math.exp(-surprisal),
    }


def optimize_reference(
    envelope: DenseEnvelope,
    *,
    alpha: float,
    row_density: float,
    outer_exponent_natural: float,
    transfer_name: str,
) -> dict[str, float | bool]:
    epsilon = 1.0e-9
    starts = [
        (alpha, row_density, max(0.02, 1.6 * alpha * row_density)),
        (alpha, 0.5, max(0.02, 0.8 * alpha)),
        (min(0.999, 1.5 * alpha), row_density, 0.2),
        (min(0.999, 2.0 * alpha), 0.5, 0.8),
        (0.25, row_density, 0.8),
        (0.5, 0.5, 1.2),
        (0.8, row_density, 2.0),
    ]
    best = None
    for probability, value_probability, surprisal in starts:
        probability = min(1.0 - epsilon, max(epsilon, probability))
        value_probability = min(
            1.0 - epsilon, max(epsilon, value_probability)
        )
        result = minimize(
            lambda point: joint_objective(
                envelope,
                alpha=alpha,
                row_density=row_density,
                outer_exponent_natural=outer_exponent_natural,
                candidate_probability=logistic(float(point[0])),
                value_probability=logistic(float(point[1])),
                surprisal=math.exp(float(point[2])),
                transfer_name=transfer_name,
            )[0],
            np.asarray(
                [
                    logit(probability),
                    logit(value_probability),
                    math.log(surprisal),
                ]
            ),
            method="Nelder-Mead",
            options={"maxiter": 500, "xatol": 1.0e-9, "fatol": 1.0e-12},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    probability = logistic(float(best.x[0]))
    value_probability = logistic(float(best.x[1]))
    surprisal = math.exp(float(best.x[2]))
    value, details = joint_objective(
        envelope,
        alpha=alpha,
        row_density=row_density,
        outer_exponent_natural=outer_exponent_natural,
        candidate_probability=probability,
        value_probability=value_probability,
        surprisal=surprisal,
        transfer_name=transfer_name,
    )
    return {
        "objective_natural_per_output_bit": value,
        "candidate_probability": probability,
        "value_probability": value_probability,
        "surprisal": surprisal,
        "optimizer_success": bool(best.success),
        **details,
    }


def load_envelope(selection: Path, delta: float) -> DenseEnvelope:
    payload = json.loads(selection.read_text(encoding="utf-8"))
    selected = payload["selected"]
    generator_words = [
        int(value, 16) for value in selected["A_generator_words_hex"]
    ]
    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    verify_transpose_columns(
        generator_words=generator_words,
        columns=columns,
        output_bits=len(columns),
    )
    weights = state_code_weights(generator_words, len(columns))
    classes, counts, association = association_matrix(weights)
    return DenseEnvelope(
        classes=classes,
        counts=counts,
        association=association,
        state_bits=len(generator_words),
        step_bits=len(columns),
        delta=delta,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--length", type=int, default=480)
    parser.add_argument(
        "--spectrum-mode",
        choices=["finite", "asymptotic", "expurgated-random-hull"],
        default="finite",
    )
    parser.add_argument("--delta", type=float, default=0.11)
    parser.add_argument(
        "--alphas",
        nargs="+",
        type=float,
        default=[0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0],
    )
    parser.add_argument("--weight-points", type=int, default=65)
    parser.add_argument("--selection-polynomial", type=float, default=2.0)
    parser.add_argument(
        "--transfer", choices=["three-state", "four-state"], default="three-state"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    length = args.length
    envelope = load_envelope(args.selection, args.delta)
    if args.spectrum_mode == "finite":
        log_spectrum = expected_ba_log_spectrum(length)
        lower = math.ceil(0.104 * length)
        upper = math.floor(0.896 * length)
        candidates = np.unique(
            np.linspace(lower, upper, args.weight_points).round().astype(int)
        )
        candidate_rows = [
            (
                int(weight),
                weight / length,
                (log_spectrum[weight] + args.selection_polynomial * math.log(length))
                / length,
            )
            for weight in candidates
            if log_spectrum[weight] + args.selection_polynomial * math.log(length)
            >= 0.0
        ]
    elif args.spectrum_mode == "asymptotic":
        densities = np.linspace(0.104, 0.896, args.weight_points)
        candidate_rows = [
            (None, float(density), ba_asymptotic_exponent(float(density)))
            for density in densities
        ]
    else:
        densities = np.linspace(0.104, 0.896, args.weight_points)
        candidate_rows = [
            (None, float(density), expurgated_random_hull(float(density)))
            for density in densities
        ]
    rows = []
    for alpha in args.alphas:
        weight_rows = []
        for weight, row_density, outer_exponent in candidate_rows:
            if outer_exponent < 0.0:
                continue
            row = optimize_reference(
                envelope,
                alpha=alpha,
                row_density=row_density,
                outer_exponent_natural=outer_exponent,
                transfer_name=args.transfer,
            )
            weight_rows.append(
                {
                    "weight": weight,
                    "row_density": row_density,
                    "outer_exponent_natural": outer_exponent,
                    **row,
                }
            )
        worst = max(
            weight_rows,
            key=lambda item: float(item["objective_natural_per_output_bit"]),
        )
        rows.append({"alpha": alpha, "worst": worst, "weights": weight_rows})

    result = {
        "schema": "golay-ba3-rm2sub-joint-diagnostic-v1",
        "status": "BINARY64_SINGLE_WEIGHT_TYPE_DIAGNOSTIC",
        "length": length,
        "spectrum_mode": args.spectrum_mode,
        "delta": args.delta,
        "selection_polynomial": args.selection_polynomial,
        "transfer": args.transfer,
        "weight_points": len(candidate_rows),
        "rows": rows,
        "all_sampled_closed": all(
            float(row["worst"]["objective_natural_per_output_bit"]) < 0.0
            for row in rows
        ),
        "limitations": [
            "Binary64 optimization is not an outward certificate.",
            "The diagnostic restricts all active BA rows to one weight class.",
            "Mixtures of BA row-weight classes require a concave-envelope proof.",
            "The finite selected-spectrum factor is B^selection_polynomial.",
            "The three-state RM2Sub envelope is used; vanishing density needs the four-state proof.",
        ],
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
