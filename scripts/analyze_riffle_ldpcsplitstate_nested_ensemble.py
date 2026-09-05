#!/usr/bin/env python3
"""Exact expected spectrum of the layered nested-code ensemble.

The ensemble samples each column of P independently and uniformly from the
weight-6 vectors in F_2^128.  It samples the S and R columns independently
and uniformly from the weight-3 vectors in F_2^64.  A uniform 128-position
interleaver separates the two accumulator layers.

The result is exact for this with-replacement ensemble.  It is a constituent
calculation, not an end-to-end distance certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from explore_riffle_ldpcsplitstate_local import accumulator_weight_transition


DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/nested_depth2_ensemble_spectrum.json"
)


def xor_fixed_weight_transition(bits: int, column_weight: int) -> np.ndarray:
    transition = np.zeros((bits + 1, bits + 1), dtype=np.float64)
    denominator = math.comb(bits, column_weight)
    for current_weight in range(bits + 1):
        minimum_intersection = max(0, column_weight - (bits - current_weight))
        maximum_intersection = min(column_weight, current_weight)
        for intersection in range(minimum_intersection, maximum_intersection + 1):
            following_weight = (
                current_weight + column_weight - 2 * intersection
            )
            transition[current_weight, following_weight] += (
                math.comb(current_weight, intersection)
                * math.comb(bits - current_weight, column_weight - intersection)
                / denominator
            )
    return transition


def iterated_distributions(
    transition: np.ndarray, maximum_steps: int
) -> np.ndarray:
    result = np.zeros((maximum_steps + 1, transition.shape[0]), dtype=np.float64)
    result[0, 0] = 1.0
    for step in range(1, maximum_steps + 1):
        result[step] = result[step - 1] @ transition
    return result


def expected_spectrum(
    *, p_column_weight: int, parity_column_weight: int, accumulator_depth: int
) -> np.ndarray:
    injection_transition = xor_fixed_weight_transition(128, p_column_weight)
    injection_by_input_weight = iterated_distributions(injection_transition, 64)

    accumulator_transition = accumulator_weight_transition(128)
    auxiliary_by_input_weight = injection_by_input_weight
    for _ in range(accumulator_depth):
        auxiliary_by_input_weight = (
            auxiliary_by_input_weight @ accumulator_transition
        )

    parity_transition = xor_fixed_weight_transition(64, parity_column_weight)
    parity_by_selected_columns = iterated_distributions(parity_transition, 192)

    spectrum = np.zeros(257, dtype=np.float64)
    for information_weight in range(65):
        information_count = math.comb(64, information_weight)
        for auxiliary_weight, auxiliary_probability in enumerate(
            auxiliary_by_input_weight[information_weight]
        ):
            if auxiliary_probability == 0.0:
                continue
            parity_distribution = parity_by_selected_columns[
                information_weight + auxiliary_weight
            ]
            base_weight = information_weight + auxiliary_weight
            spectrum[
                base_weight : base_weight + len(parity_distribution)
            ] += information_count * auxiliary_probability * parity_distribution
    return spectrum


def chernoff_scan(
    spectrum: np.ndarray, *, epochs: int, target_distance: int
) -> dict[str, object]:
    nonzero_spectrum = spectrum.copy()
    nonzero_spectrum[0] -= 1.0
    denominator = (1 << 64) - 1
    grid = np.linspace(-2.0, 1.25, 131)
    best: dict[str, float] | None = None
    rows = []
    weights = np.arange(len(spectrum))
    for log_surprisal in grid:
        surprisal = math.exp(float(log_surprisal))
        z = math.exp(-surprisal)
        moment = float(nonzero_spectrum @ (z**weights)) / denominator
        bound = (
            epochs * math.log2(moment)
            + target_distance * surprisal / math.log(2.0)
        )
        row = {
            "log_surprisal": float(log_surprisal),
            "z": z,
            "live_moment_log2": math.log2(moment),
            "repeated_target_log2_bound": bound,
        }
        rows.append(row)
        if best is None or bound < best["repeated_target_log2_bound"]:
            best = row
    assert best is not None
    return {
        "best": best,
        "repeated_live_suppression_bits": -best["repeated_target_log2_bound"],
        "rows": rows,
    }


def expurgated_chernoff_scans(
    spectrum: np.ndarray,
    *,
    epochs: int,
    target_distance: int,
    structural_event_probability: float = 1.0,
    structural_event_name: str = "none",
) -> list[dict[str, object]]:
    rows = []
    for cutoff in (23, 31, 39, 47, 55, 63):
        expected_bad_words = float(np.sum(spectrum[1 : cutoff + 1]))
        joint_probability = structural_event_probability - expected_bad_words
        if joint_probability <= 0.0:
            rows.append(
                {
                    "cutoff": cutoff,
                    "expected_bad_words": expected_bad_words,
                    "structural_event": structural_event_name,
                    "structural_event_probability": structural_event_probability,
                    "joint_good_probability_lower_bound": 0.0,
                    "available": False,
                }
            )
            continue
        conditional_upper_spectrum = spectrum.copy()
        conditional_upper_spectrum[1 : cutoff + 1] = 0.0
        conditional_upper_spectrum[cutoff + 1 :] /= joint_probability
        scan = chernoff_scan(
            conditional_upper_spectrum,
            epochs=epochs,
            target_distance=target_distance,
        )
        rows.append(
            {
                "cutoff": cutoff,
                "expected_bad_words": expected_bad_words,
                "structural_event": structural_event_name,
                "structural_event_probability": structural_event_probability,
                "joint_good_probability_lower_bound": joint_probability,
                "available": True,
                "conditional_expected_moment_upper_bound": scan,
                "claim": (
                    "There exists a sampled code with no nonzero word at "
                    f"weight at most {cutoff} and live-state moment no larger "
                    "than this conditional-expectation bound."
                ),
            }
        )
    return rows


def distinct_sampling_probability(population: int, draws: int) -> float:
    probability = 1.0
    for used in range(draws):
        probability *= (population - used) / population
    return probability


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    spectrum = expected_spectrum(
        p_column_weight=args.p_column_weight,
        parity_column_weight=args.parity_column_weight,
        accumulator_depth=args.accumulator_depth,
    )
    target_distance = math.floor(
        args.relative_distance * 2 * args.message_bits
    )
    epochs = (2 * args.message_bits) // 256
    b_distinct_probability = distinct_sampling_probability(
        math.comb(64, args.parity_column_weight), 192
    )
    expected_low_counts = {
        str(threshold): float(np.sum(spectrum[1 : threshold + 1]))
        for threshold in (23, 31, 39, 47, 55, 63)
    }
    positive = np.flatnonzero(spectrum[1:] > 0.0) + 1
    return {
        "schema": "riffle-ldpcsplitstate-nested-ensemble-spectrum-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "ensemble": {
            "P_columns": (
                f"64 independent uniform weight-{args.p_column_weight} "
                "vectors in F2^128"
            ),
            "accumulator_depth": args.accumulator_depth,
            "interleavers": (
                "independent uniform 128-position permutation between "
                "adjacent accumulator layers"
            ),
            "S_R_columns": (
                f"192 independent uniform weight-{args.parity_column_weight} "
                "vectors in F2^64"
            ),
            "sampling": "with replacement",
        },
        "parameters": {
            "message_bits": args.message_bits,
            "relative_distance": args.relative_distance,
            "epochs": epochs,
            "target_distance": target_distance,
        },
        "checks": {
            "expected_codewords": float(np.sum(spectrum)),
            "target_codewords": float(1 << 64),
            "expected_zero_weight_words": float(spectrum[0]),
            "first_positive_expected_weight": int(positive[0]),
        },
        "expected_low_weight_codewords": expected_low_counts,
        "expected_weight_spectrum": {
            str(weight): float(count)
            for weight, count in enumerate(spectrum)
            if count > 0.0
        },
        "live_state_chernoff": chernoff_scan(
            spectrum, epochs=epochs, target_distance=target_distance
        ),
        "expurgated_live_state_chernoff": expurgated_chernoff_scans(
            spectrum, epochs=epochs, target_distance=target_distance
        ),
        "distinct_B_expurgated_live_state_chernoff": expurgated_chernoff_scans(
            spectrum,
            epochs=epochs,
            target_distance=target_distance,
            structural_event_probability=b_distinct_probability,
            structural_event_name=(
                "all 192 degree-3 S/R columns are distinct; this gives B "
                "no kernel words of weights two or three"
            ),
        ),
        "scope": (
            "Exact expected constituent spectrum for the stated sampled "
            "ensemble. It does not assert concentration for one fixed pair "
            "and does not include zero-state activation, termination, the "
            "outer spectrum, placement, or the complete recurrence."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p-column-weight", type=int, default=6)
    parser.add_argument("--parity-column-weight", type=int, default=3)
    parser.add_argument("--accumulator-depth", type=int, default=2)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checks = payload["checks"]
    print(
        f"expected_codewords,{checks['expected_codewords']:.6f},"
        f"zero_weight,{checks['expected_zero_weight_words']:.6e}"
    )
    for threshold, count in payload["expected_low_weight_codewords"].items():
        print(f"expected_words_le_{threshold},{count:.9e}")
    chernoff = payload["live_state_chernoff"]
    print(
        "live_suppression_bits,"
        f"{chernoff['repeated_live_suppression_bits']:.6f},"
        f"best_log_surprisal,{chernoff['best']['log_surprisal']:.6f}"
    )
    for row in payload["expurgated_live_state_chernoff"]:
        if row["available"]:
            suppression = row["conditional_expected_moment_upper_bound"][
                "repeated_live_suppression_bits"
            ]
            print(
                f"expurgated_through_{row['cutoff']},"
                f"good_probability,{row['joint_good_probability_lower_bound']:.9f},"
                f"live_suppression_bits,{suppression:.6f}"
            )
    for row in payload["distinct_B_expurgated_live_state_chernoff"]:
        if row["available"]:
            suppression = row["conditional_expected_moment_upper_bound"][
                "repeated_live_suppression_bits"
            ]
            print(
                f"distinct_B_expurgated_through_{row['cutoff']},"
                f"joint_probability,{row['joint_good_probability_lower_bound']:.9f},"
                f"live_suppression_bits,{suppression:.6f}"
            )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
