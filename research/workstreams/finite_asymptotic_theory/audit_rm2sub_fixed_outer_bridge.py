#!/usr/bin/env python3
"""Audit the fixed-outer bridge from Toeplitz convolution to RM2Sub-S19.

This is a structural audit and a nearest-binary64 diagnostic.  It verifies
the selected RM2Sub A and B maps, computes the unavoidable silent-subcode
rank budget, and measures an optimistic live-only Chernoff contraction.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[2]
WORKSTREAM = Path(__file__).resolve().parent
SELECTION = ROOT / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts/min_state/s19_rm2sub_selection.json"
)
B_KERNEL = ROOT / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts/min_state/"
    "s19_rm2sub_b_kernel_spectrum.json"
)
A_SPECTRUM = ROOT / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts/min_state/"
    "s19_rm2sub_a_spectrum_audit.json"
)
TOEPLITZ_RECEIPT = WORKSTREAM / "ba3_B240_repeated_toeplitz_prefix_outward.json"
OUTPUT = WORKSTREAM / "rm2sub_fixed_outer_bridge_audit.json"

STEP_BITS = 128
STATE_BITS = 19
STATE_SIZE = 1 << STATE_BITS


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gf2_rank(words: list[int], width: int) -> int:
    basis = [0] * width
    result = 0
    for source in words:
        value = source
        while value:
            pivot = value.bit_length() - 1
            if basis[pivot]:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                result += 1
                break
    return result


def syndrome(word: int, columns: list[int]) -> int:
    result = 0
    while word:
        bit = (word & -word).bit_length() - 1
        result ^= columns[bit]
        word &= word - 1
    return result


def enumerate_a_spectrum(generators: list[int]) -> dict[int, int]:
    spectrum: dict[int, int] = {0: 1}
    word = 0
    previous_gray = 0
    for counter in range(1, 1 << len(generators)):
        gray = counter ^ (counter >> 1)
        changed = gray ^ previous_gray
        word ^= generators[changed.bit_length() - 1]
        weight = word.bit_count()
        spectrum[weight] = spectrum.get(weight, 0) + 1
        previous_gray = gray
    return spectrum


def log_weight_enumerator(theta: float, spectrum: list[tuple[int, int]]) -> float:
    terms = [math.log(count) - theta * weight for weight, count in spectrum]
    maximum = max(terms)
    return maximum + math.log(sum(math.exp(term - maximum) for term in terms))


def optimistic_live_only_bound(
    live_epochs: int,
    distance: int,
    spectrum: list[tuple[int, int]],
) -> dict[str, float | int]:
    """Ignore terminations and deterministic output outside live epochs."""

    def objective(theta: float) -> float:
        log_moment = log_weight_enumerator(theta, spectrum) - math.log(
            STATE_SIZE - 1
        )
        return (distance * theta + live_epochs * log_moment) / math.log(2.0)

    result = minimize_scalar(
        objective,
        bounds=(0.0, 5.0),
        method="bounded",
        options={"xatol": 1.0e-13},
    )
    return {
        "live_epochs": live_epochs,
        "z": math.exp(-float(result.x)),
        "log2_chernoff_bound": float(result.fun),
        "margin_bits": -float(result.fun),
    }


def log2_hamming_ball(length: int, radius: int) -> float:
    """Nearest-binary64 log2 of sum_{w=0}^radius C(length,w)."""
    largest = (
        math.lgamma(length + 1)
        - math.lgamma(radius + 1)
        - math.lgamma(length - radius + 1)
    )
    total = 1.0
    ratio_product = 1.0
    for weight in range(radius, 0, -1):
        ratio_product *= weight / (length - weight + 1)
        total += ratio_product
        if ratio_product < total * 1.0e-18:
            break
    return (largest + math.log(total)) / math.log(2.0)


def dimension_only_weight_bound(dimension: int, distance: int) -> dict[str, float]:
    """Optimize z^-D(1+z)^dimension, a universal systematic bound."""

    def objective(theta: float) -> float:
        return (
            dimension * math.log1p(math.exp(-theta)) + distance * theta
        ) / math.log(2.0)

    result = minimize_scalar(objective, bounds=(0.0, 10.0), method="bounded")
    return {
        "log2_upper": float(result.fun),
        "z": math.exp(-float(result.x)),
    }


def main() -> None:
    selection_payload = json.loads(SELECTION.read_text(encoding="utf-8"))
    selected = selection_payload["selected"]
    a_generators = [int(value, 16) for value in selected["A_generator_words_hex"]]
    b_columns = [int(value, 16) for value in selected["B_columns_hex"]]
    if len(a_generators) != STATE_BITS or len(b_columns) != STEP_BITS:
        raise ArithmeticError("selected RM2Sub dimensions do not match S19/T128")
    if gf2_rank(a_generators, STEP_BITS) != STATE_BITS:
        raise ArithmeticError("A is not injective")
    if gf2_rank(b_columns, STATE_BITS) != STATE_BITS:
        raise ArithmeticError("B is not surjective")
    if any(syndrome(word, b_columns) for word in a_generators):
        raise ArithmeticError("BA is not zero")

    observed_spectrum = enumerate_a_spectrum(a_generators)
    spectrum_payload = json.loads(A_SPECTRUM.read_text(encoding="utf-8"))
    expected_spectrum = {
        int(row["weight"]): int(row["count"])
        for row in spectrum_payload["spectrum"]
    }
    if observed_spectrum != expected_spectrum:
        raise ArithmeticError("enumerated A spectrum does not match its receipt")

    kernel_payload = json.loads(B_KERNEL.read_text(encoding="utf-8"))
    if int(kernel_payload["parameters"]["kernel_dimension"]) != STEP_BITS - STATE_BITS:
        raise ArithmeticError("B-kernel dimension mismatch")
    if int(kernel_payload["checks"]["minimum_kernel_distance"]) != 6:
        raise ArithmeticError("unexpected B-kernel minimum distance")

    toeplitz = json.loads(TOEPLITZ_RECEIPT.read_text(encoding="utf-8"))
    parameters = toeplitz["parameters"]
    output_bits = int(parameters["output_bits"])
    parent_dimension = int(parameters["parent_dimension"])
    target_dimension = int(parameters["target_dimension_after_zero_shortening"])
    distance = int(parameters["distance_cutoff"])
    if output_bits % STEP_BITS:
        raise ArithmeticError("RM2Sub step length does not divide N")
    epochs = output_bits // STEP_BITS
    maximum_syndrome_rank = STATE_BITS * epochs
    parent_silent_lower = parent_dimension - maximum_syndrome_rank
    target_silent_lower = target_dimension - maximum_syndrome_rank
    if target_silent_lower <= 0:
        raise ArithmeticError("unexpected nonpositive silent rank lower bound")

    spectrum = sorted(observed_spectrum.items())
    optimistic_rows = [
        optimistic_live_only_bound(span, distance, spectrum)
        for span in (epochs - 1, 10_000, 5_000, 4_000, 3_700, 3_684, 3_683)
    ]
    ball_log2 = log2_hamming_ball(output_bits, distance)
    random_silent_log2_expected = target_silent_lower - output_bits + ball_log2
    outer_bits = int(parameters["outer_bits"])
    outer_rows = int(parameters["outer_rows"])
    ba_b_squared_cost = 2.0 * outer_rows * math.log2(outer_bits)

    lower_span = 0
    upper_span = epochs - 1
    while lower_span + 1 < upper_span:
        middle = (lower_span + upper_span) // 2
        middle_log2 = optimistic_live_only_bound(
            middle, distance, spectrum
        )["log2_chernoff_bound"]
        if middle_log2 <= -40.0:
            upper_span = middle
        else:
            lower_span = middle
    minimum_span_40 = upper_span

    nonactivation = []
    for row in kernel_payload["by_total_weight"]:
        weight = int(row["total_weight"])
        kernel_words = int(row["kernel_words"])
        shell_size = int(row["shell_size"])
        if weight in (0, 1, 2, 3, 4, 5, 6, 8, 64, 120, 128):
            probability = kernel_words / shell_size
            nonactivation.append(
                {
                    "input_weight": weight,
                    "kernel_words": kernel_words,
                    "shell_size": shell_size,
                    "probability": probability,
                    "log2_probability": (
                        math.log2(probability) if probability else None
                    ),
                }
            )

    payload = {
        "schema": "rm2sub-s19-fixed-outer-bridge-audit-v1",
        "status": "EXACT_STRUCTURE_AND_BINARY64_DIAGNOSTIC",
        "conclusion": {
            "prefix_rank_port_is_not_direct": True,
            "reason": (
                "RM2Sub observes at most 19 syndrome bits per 128-bit epoch. "
                "A large subcode has zero syndrome in every epoch and bypasses "
                "all random multipliers. Raw prefix dimensions do not bound "
                "the weight spectrum of this silent subcode."
            ),
            "next_exact_object": (
                "the weight enumerators of suffix-silent subcodes, together "
                "with a zero/live/punctured RM2Sub transfer"
            ),
        },
        "parameters": {
            "output_bits": output_bits,
            "parent_dimension": parent_dimension,
            "target_dimension": target_dimension,
            "distance_cutoff": distance,
            "step_bits": STEP_BITS,
            "state_bits": STATE_BITS,
            "epochs": epochs,
        },
        "exact_checks": {
            "rank_A": gf2_rank(a_generators, STEP_BITS),
            "rank_B": gf2_rank(b_columns, STATE_BITS),
            "BA_zero": True,
            "A_spectrum_matches_receipt": True,
            "A_minimum_distance": min(weight for weight in observed_spectrum if weight),
            "B_kernel_dimension": STEP_BITS - STATE_BITS,
            "B_kernel_minimum_distance": int(
                kernel_payload["checks"]["minimum_kernel_distance"]
            ),
        },
        "silent_subcode_rank_budget": {
            "maximum_total_syndrome_rank": maximum_syndrome_rank,
            "parent_silent_dimension_lower": parent_silent_lower,
            "parent_silent_rate_lower": parent_silent_lower / output_bits,
            "shortened_silent_dimension_lower": target_silent_lower,
            "shortened_silent_rate_lower": target_silent_lower / output_bits,
            "identity_on_silent_subcode": True,
        },
        "silent_subcode_spectrum_diagnostics": {
            "uniform_random_code_log2_expected_bad_words": random_silent_log2_expected,
            "uniform_random_code_margin_bits": -random_silent_log2_expected,
            "B_squared_selection_cost_all_rows_bits": ba_b_squared_cost,
            "random_benchmark_margin_after_B_squared_cost_bits": (
                -random_silent_log2_expected - ba_b_squared_cost
            ),
            "dimension_only_systematic_enumerator": dimension_only_weight_bound(
                target_silent_lower, distance
            ),
            "interpretation": (
                "A random code at the forced silent dimension has enormous "
                "margin, but dimension alone gives no useful upper bound. "
                "The missing information is the silent subcode's spectrum."
            ),
        },
        "selected_B_nonactivation_probabilities": nonactivation,
        "optimistic_live_only_diagnostic": {
            "scope": (
                "Uses the proved maximum-coset bound W_A(z)/(2^19-1), "
                "but ignores state termination and every deterministic-output "
                "segment. It is not an RM2Sub code bound."
            ),
            "minimum_live_epochs_for_single_word_40_bit_bound": minimum_span_40,
            "rows": optimistic_rows,
        },
        "dependencies": [
            {"path": path.relative_to(ROOT).as_posix(), "sha256": sha256(path)}
            for path in (SELECTION, B_KERNEL, A_SPECTRUM, TOEPLITZ_RECEIPT, Path(__file__))
        ],
    }
    with OUTPUT.open("w", encoding="utf-8", newline="\n") as output_file:
        output_file.write(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"output": str(OUTPUT), **payload["silent_subcode_rank_budget"]}, indent=2))


if __name__ == "__main__":
    main()
