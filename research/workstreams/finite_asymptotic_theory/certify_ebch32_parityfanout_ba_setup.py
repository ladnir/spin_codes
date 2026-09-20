#!/usr/bin/env python3
"""Outward setup/spectrum certificate for EBCH32--PF31x33--BA-3/B=256.

The proof path uses nonnegative arithmetic.  Each binary64 arithmetic result
used as an upper bound is advanced one ULP toward +infinity.  Exact positive
integer denominators are rounded toward -infinity before division.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import platform
import sys


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json"

B = 256
SOURCE_SIZE = 31
TARGET_SIZE = 33
LOWER_WEIGHT = 24
UPPER_WEIGHT = 232
LOCAL_SPECTRUM = {
    0: 1,
    8: 620,
    12: 13_888,
    16: 36_518,
    20: 13_888,
    24: 620,
    32: 1,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def up(value: float) -> float:
    if math.isnan(value):
        raise ArithmeticError("NaN in outward calculation")
    return value if value == math.inf else math.nextafter(value, math.inf)


def down(value: float) -> float:
    if math.isnan(value):
        raise ArithmeticError("NaN in outward calculation")
    return math.nextafter(value, -math.inf)


def add_up(left: float, right: float) -> float:
    if left < 0.0 or right < 0.0:
        raise ArithmeticError("negative addend in positive proof path")
    return up(left + right)


def mul_up(left: float, right: float) -> float:
    if left < 0.0 or right < 0.0:
        raise ArithmeticError("negative factor in positive proof path")
    if left == 0.0 or right == 0.0:
        return 0.0
    return up(left * right)


def div_up(numerator: float, denominator: float) -> float:
    if numerator < 0.0 or denominator <= 0.0:
        raise ArithmeticError("invalid positive outward division")
    return 0.0 if numerator == 0.0 else up(numerator / denominator)


def integer_upper(value: int) -> float:
    if value < 0:
        raise ArithmeticError("negative integer in positive proof path")
    return 0.0 if value == 0 else up(float(value))


def integer_lower(value: int) -> float:
    if value <= 0:
        raise ArithmeticError("nonpositive denominator")
    return down(float(value))


def ratio_upper(numerator: int, denominator: int) -> float:
    return div_up(integer_upper(numerator), integer_lower(denominator))


def direct_sum_counts() -> list[int]:
    counts = [1]
    for _ in range(8):
        updated = [0] * (len(counts) + 32)
        for left_weight, left_count in enumerate(counts):
            for right_weight, right_count in LOCAL_SPECTRUM.items():
                updated[left_weight + right_weight] += left_count * right_count
        counts = updated
    if len(counts) != B + 1 or sum(counts) != 1 << 128:
        raise ArithmeticError("EBCH32 direct-sum mass check failed")
    if counts[B] != 1:
        raise ArithmeticError("expected a unique direct-sum all-ones word")
    return counts


def valid_overlap(total: int, selected: int) -> bool:
    return 0 <= selected <= total


def fanout_spectrum_upper(counts: list[int]) -> list[float]:
    result = [0.0] * (B + 1)
    source_denominator = math.comb(B, SOURCE_SIZE)
    target_denominator = math.comb(B - SOURCE_SIZE, TARGET_SIZE)
    for weight, count in enumerate(counts):
        if count == 0:
            continue
        source_numerators: list[tuple[int, int]] = []
        for source_overlap in range(SOURCE_SIZE + 1):
            left = source_overlap
            right = SOURCE_SIZE - source_overlap
            if not (valid_overlap(weight, left) and valid_overlap(B - weight, right)):
                continue
            numerator = math.comb(weight, left) * math.comb(B - weight, right)
            source_numerators.append((source_overlap, numerator))
        if sum(numerator for _, numerator in source_numerators) != source_denominator:
            raise ArithmeticError("source-overlap Vandermonde check failed")

        multiplicity = integer_upper(count)
        for source_overlap, source_numerator in source_numerators:
            source_mass = ratio_upper(source_numerator, source_denominator)
            if source_overlap % 2 == 0:
                result[weight] = add_up(
                    result[weight], mul_up(multiplicity, source_mass)
                )
                continue

            remaining_weight = weight - source_overlap
            target_numerators: list[tuple[int, int]] = []
            for target_overlap in range(TARGET_SIZE + 1):
                left = target_overlap
                right = TARGET_SIZE - target_overlap
                if not (
                    valid_overlap(remaining_weight, left)
                    and valid_overlap(B - SOURCE_SIZE - remaining_weight, right)
                ):
                    continue
                numerator = math.comb(remaining_weight, left) * math.comb(
                    B - SOURCE_SIZE - remaining_weight, right
                )
                target_numerators.append((target_overlap, numerator))
            if sum(numerator for _, numerator in target_numerators) != target_denominator:
                raise ArithmeticError("target-overlap Vandermonde check failed")
            for target_overlap, target_numerator in target_numerators:
                output_weight = weight + TARGET_SIZE - 2 * target_overlap
                probability = mul_up(
                    source_mass,
                    ratio_upper(target_numerator, target_denominator),
                )
                result[output_weight] = add_up(
                    result[output_weight], mul_up(multiplicity, probability)
                )

    # For the unique all-ones word, a=31 and b=33 deterministically.
    if counts[B] != 1:
        raise ArithmeticError("all-ones input audit failed")
    all_ones_output = B + TARGET_SIZE - 2 * TARGET_SIZE
    if all_ones_output != 223:
        raise ArithmeticError("all-ones fanout-weight audit failed")
    return result


def accumulator_transition_upper(input_weight: int, output_weight: int) -> float:
    if input_weight == 0:
        return 1.0 if output_weight == 0 else 0.0
    runs = (input_weight + 1) // 2
    second_runs = input_weight - runs
    if output_weight < runs or B - output_weight < second_runs:
        return 0.0
    numerator = math.comb(output_weight - 1, runs - 1) * math.comb(
        B - output_weight, second_runs
    )
    return ratio_upper(numerator, math.comb(B, input_weight))


def apply_accumulator_upper(spectrum: list[float]) -> list[float]:
    result = [0.0] * (B + 1)
    for input_weight, multiplicity in enumerate(spectrum):
        if multiplicity == 0.0:
            continue
        for output_weight in range(B + 1):
            transition = accumulator_transition_upper(input_weight, output_weight)
            if transition:
                result[output_weight] = add_up(
                    result[output_weight], mul_up(multiplicity, transition)
                )
    return result


def expected_spectrum_upper() -> list[float]:
    """Return an outward upper bound on every final expected multiplicity."""
    counts = direct_sum_counts()
    fanout = fanout_spectrum_upper(counts)
    return apply_accumulator_upper(apply_accumulator_upper(fanout))


def conditioned_spectrum_upper() -> tuple[list[float], float, float]:
    """Apply the Markov conditioning envelope for weights 24 through 232."""
    after_two = expected_spectrum_upper()
    tail = 0.0
    for weight in list(range(1, LOWER_WEIGHT)) + list(
        range(UPPER_WEIGHT + 1, B + 1)
    ):
        tail = add_up(tail, after_two[weight])
    if not tail < 1.0:
        raise ArithmeticError("Markov lower bound is nonpositive")
    good_lower = down(1.0 - tail)
    conditioned = [0.0] * (B + 1)
    for weight in range(LOWER_WEIGHT, UPPER_WEIGHT + 1):
        conditioned[weight] = div_up(after_two[weight], good_lower)
    return conditioned, tail, good_lower


def certify() -> dict[str, object]:
    if sys.float_info.radix != 2 or sys.float_info.mant_dig != 53:
        raise RuntimeError("the verifier requires IEEE-754 binary64")
    counts = direct_sum_counts()
    after_two = expected_spectrum_upper()
    _, tail, good_lower = conditioned_spectrum_upper()
    tail_rows = []
    for weight in list(range(1, LOWER_WEIGHT)) + list(
        range(UPPER_WEIGHT + 1, B + 1)
    ):
        if after_two[weight] != 0.0:
            tail_rows.append(
                {
                    "weight": weight,
                    "expected_multiplicity_upper_hex": after_two[weight].hex(),
                }
            )
    expected_trials_upper = div_up(1.0, good_lower)

    return {
        "schema": "ebch32-parityfanout31x33-ba3-b256-setup-outward-v1",
        "status": "OUTWARD_CERTIFICATE",
        "claim": {
            "permitted_final_nonzero_weights": [LOWER_WEIGHT, UPPER_WEIGHT],
            "expected_rejected_word_count_upper_hex": tail.hex(),
            "good_event_probability_lower_hex": good_lower.hex(),
            "expected_rejection_trials_per_row_upper_hex": expected_trials_upper.hex(),
            "display_expected_rejected_word_count_upper": tail,
            "display_good_event_probability_lower": good_lower,
            "display_expected_trials_upper": expected_trials_upper,
        },
        "parameters": {
            "row_length": B,
            "row_dimension": 128,
            "constituent_count": 8,
            "constituent": {
                "parameters": [32, 16, 8],
                "weight_spectrum": LOCAL_SPECTRUM,
            },
            "parity_fanout": {
                "source_size": SOURCE_SIZE,
                "target_size": TARGET_SIZE,
                "sampling": "S uniform; T uniform in the complement of S",
            },
            "accumulator_count": 2,
            "independent_uniform_accumulator_interleavers": 2,
        },
        "exact_identity_audits": {
            "direct_sum_total_words": str(sum(counts)),
            "direct_sum_total_words_expected": str(1 << 128),
            "direct_sum_all_ones_multiplicity": counts[B],
            "all_ones_source_overlap": SOURCE_SIZE,
            "all_ones_target_overlap": TARGET_SIZE,
            "all_ones_fanout_output_weight": 223,
            "source_and_target_vandermonde_checks": "passed for every nonzero input-weight coefficient",
        },
        "arithmetic": {
            "format": "IEEE-754 binary64",
            "proof_rule": (
                "Every nonnegative proof-path addition, multiplication, and "
                "division is rounded one ULP toward +infinity; exact positive "
                "integer denominators are rounded toward -infinity."
            ),
            "transcendental_operations_in_proof_path": 0,
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "tail_weight_rows": tail_rows,
        "dependencies": [
            {
                "path": "workstreams/finite_asymptotic_theory/certify_ebch32_parityfanout_ba_setup.py",
                "sha256": sha256(Path(__file__).resolve()),
            },
            {
                "path": "workstreams/finite_asymptotic_theory/ebch32_16_delta8_spectrum.csv",
                "sha256": sha256(WORKSTREAM / "ebch32_16_delta8_spectrum.csv"),
            },
        ],
        "limitations": [
            "This receipt certifies the expected row spectrum tail and the Markov setup bound only.",
            "It does not certify any RM2Sub distance contribution.",
            "An efficient exact or one-sided-safe decision procedure for the accepted-row event remains open.",
        ],
    }


def main() -> None:
    payload = certify()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"receipt={OUTPUT}")


if __name__ == "__main__":
    main()
