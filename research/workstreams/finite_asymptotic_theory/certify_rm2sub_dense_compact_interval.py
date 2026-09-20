#!/usr/bin/env python3
"""Outward interval certificate for RM2Sub density 10^-4 <= alpha <= 1."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from mpmath import iv

from analyze_rm2sub_dense_occupation import (
    DEFAULT_SELECTION,
    DenseEnvelope,
    association_matrix,
    state_code_weights,
    verify_transpose_columns,
)


iv.dps = 100


def as_interval(value: Fraction):
    return iv.mpf(value.numerator) / value.denominator


def upper_float(value) -> float:
    return math.nextafter(float(value.b), math.inf)


def point_fraction(value: float, digits: int = 17) -> Fraction:
    return Fraction(format(value, f".{digits}g"))


def interval_three_state_transfer(
    *,
    candidate_probability: Fraction,
    z: Fraction,
    classes: list[int],
    counts: list[int],
    association: list[list[int]],
    state_space: int,
    step_bits: int,
) -> list[list[object]]:
    p = as_interval(candidate_probability)
    z_iv = as_interval(z)
    beta = p / 2
    u = 1 - beta + beta * z_iv
    v = beta + (1 - beta) * z_iv
    signed = 1 - beta - beta * z_iv
    live = [u ** (step_bits - w) * v**w for w in classes]
    syndrome = [u ** (step_bits - w) * signed**w for w in classes]
    zero_to_zero = sum(
        count * moment
        for count, moment in zip(counts, syndrome, strict=True)
    ) / state_space
    zero_total = u**step_bits
    zero_to_deterministic = zero_total - zero_to_zero
    paired_total = iv.mpf(0)
    for row, syndrome_moment in enumerate(syndrome):
        paired_total += syndrome_moment * sum(
            association[row][column] * live[column]
            for column in range(len(classes))
        )
    paired_total /= state_space
    paired_nonzero = paired_total - zero_to_zero * live[0]
    deterministic = paired_nonzero / zero_to_deterministic
    live_states = state_space - 1
    uniform_live = (
        sum(
            count * moment
            for count, moment in zip(counts, live, strict=True)
        )
        - live[0]
    ) / live_states
    punctured = as_interval(Fraction(live_states, live_states - 1)) * uniform_live
    return [
        [zero_to_zero, zero_to_deterministic, iv.mpf(0)],
        [deterministic / live_states, iv.mpf(0), deterministic],
        [punctured / live_states, iv.mpf(0), punctured],
    ]


def collatz_upper(
    transfer: list[list[object]], vector: list[Fraction]
) -> Fraction:
    ratios = []
    for row, denominator in zip(transfer, vector, strict=True):
        numerator = sum(
            entry * as_interval(coordinate)
            for entry, coordinate in zip(row, vector, strict=True)
        )
        ratios.append(numerator / as_interval(denominator))
    outward = math.nextafter(max(upper_float(value) for value in ratios), math.inf)
    return point_fraction(outward)


def exponent_upper(
    *,
    alpha: Fraction,
    candidate_probability: Fraction,
    z: Fraction,
    radius_upper: Fraction,
) -> float:
    alpha_iv = as_interval(alpha)
    p = as_interval(candidate_probability)
    result = alpha_iv * iv.log(iv.mpf(2)) / 2
    if alpha > 0:
        result += alpha_iv * iv.log(alpha_iv / p)
    if alpha < 1:
        complement = 1 - alpha_iv
        result += complement * iv.log(complement / (1 - p))
    result += iv.log(as_interval(radius_upper)) / 128
    result -= as_interval(Fraction(11, 100)) * iv.log(as_interval(z))
    return upper_float(result)


@dataclass(frozen=True)
class AlphaBox:
    lower: Fraction
    upper: Fraction

    def split(self) -> tuple["AlphaBox", "AlphaBox"]:
        midpoint = (self.lower + self.upper) / 2
        return (
            AlphaBox(self.lower, midpoint),
            AlphaBox(midpoint, self.upper),
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--maximum-depth", type=int, default=30)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "rm2sub_dense_compact_interval_d11.json"
        ),
    )
    args = parser.parse_args()

    payload = json.loads(args.selection.read_text(encoding="utf-8"))
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
    classes_array, counts_array, association_array = association_matrix(weights)
    classes = [int(value) for value in classes_array]
    counts = [int(value) for value in counts_array]
    association = [[int(value) for value in row] for row in association_array]
    state_space = 1 << len(generator_words)
    diagnostic = DenseEnvelope(
        classes=classes_array,
        counts=counts_array,
        association=association_array,
        state_bits=len(generator_words),
        step_bits=len(columns),
        delta=0.11,
    )

    stack = [(AlphaBox(Fraction(1, 10_000), Fraction(1)), 0)]
    accepted: list[dict[str, object]] = []
    processed = 0
    deepest = 0
    unresolved = None
    while stack:
        box, depth = stack.pop()
        processed += 1
        deepest = max(deepest, depth)
        midpoint = (box.lower + box.upper) / 2
        optimized = diagnostic.optimize_variant(float(midpoint), "three_state")
        p = point_fraction(float(optimized["candidate_probability"]))
        z = point_fraction(math.exp(-float(optimized["surprisal"])))
        numeric_endpoints = [
            diagnostic.objective(
                alpha=float(endpoint),
                candidate_probability=float(p),
                surprisal=-math.log(float(z)),
                envelope="three_state",
            )[0]
            for endpoint in (box.lower, box.upper)
        ]
        if max(numeric_endpoints) >= -1e-10:
            if depth >= args.maximum_depth:
                unresolved = {
                    "lower": str(box.lower),
                    "upper": str(box.upper),
                    "depth": depth,
                    "numeric_endpoint_values": numeric_endpoints,
                }
                break
            left, right = box.split()
            stack.append((right, depth + 1))
            stack.append((left, depth + 1))
            continue

        float_transfer, _ = diagnostic.transfer(
            candidate_probability=float(p),
            surprisal=-math.log(float(z)),
        )
        eigenvalues, eigenvectors = np.linalg.eig(float_transfer)
        index = int(np.argmax(np.abs(eigenvalues)))
        float_vector = np.abs(np.real(eigenvectors[:, index]))
        float_vector /= np.max(float_vector)
        vector = [
            max(Fraction(1, 10**30), point_fraction(float(value)))
            for value in float_vector
        ]
        transfer = interval_three_state_transfer(
            candidate_probability=p,
            z=z,
            classes=classes,
            counts=counts,
            association=association,
            state_space=state_space,
            step_bits=len(columns),
        )
        radius_upper = collatz_upper(transfer, vector)
        endpoint_uppers = [
            exponent_upper(
                alpha=endpoint,
                candidate_probability=p,
                z=z,
                radius_upper=radius_upper,
            )
            for endpoint in (box.lower, box.upper)
        ]
        if max(endpoint_uppers) < 0:
            accepted.append(
                {
                    "lower": str(box.lower),
                    "upper": str(box.upper),
                    "candidate_probability": str(p),
                    "z": str(z),
                    "collatz_vector": [str(value) for value in vector],
                    "radius_upper": str(radius_upper),
                    "endpoint_exponent_uppers_natural": endpoint_uppers,
                }
            )
            continue
        if depth >= args.maximum_depth:
            unresolved = {
                "lower": str(box.lower),
                "upper": str(box.upper),
                "depth": depth,
                "interval_endpoint_uppers": endpoint_uppers,
            }
            break
        left, right = box.split()
        stack.append((right, depth + 1))
        stack.append((left, depth + 1))

    accepted.sort(key=lambda cell: Fraction(str(cell["lower"])))
    result = {
        "schema": "rm2sub-dense-compact-interval-v1",
        "status": "proved" if unresolved is None else "unresolved",
        "claim": (
            "For delta=0.11, the three-state random-outer RM2Sub exponent "
            "is strictly negative for every 10^-4 <= alpha <= 1."
        ),
        "interval_dps": iv.dps,
        "processed_boxes": processed,
        "accepted_boxes": len(accepted),
        "deepest_box": deepest,
        "largest_accepted_upper_natural": max(
            value
            for cell in accepted
            for value in cell["endpoint_exponent_uppers_natural"]
        ) if accepted else None,
        "convexity_rule": (
            "For fixed p,z and Collatz radius bound, the exponent is convex "
            "in alpha; endpoint negativity certifies the whole box."
        ),
        "boxes": accepted,
        "unresolved": unresolved,
        "limitations": [
            "Binary64 optimization proposes witnesses but proves no box.",
            "All accepted inequalities use 100-digit outward intervals.",
            "The interval begins at alpha=10^-4.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "boxes"},
            indent=2,
            sort_keys=True,
        )
    )
    if unresolved is not None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
