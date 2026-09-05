#!/usr/bin/env python3
"""Outward certificate for the finite B=240 occupation-one contribution.

The verifier uses only nonnegative arithmetic in the proof path.  Every
binary64 addition, multiplication, and division is followed by one directed
step toward +infinity.  Exact integer inputs are bracketed before conversion.
The diagnostic receipt selects Chernoff witnesses; its numerical claims are
not trusted.  Any selected binary64 value z in (0,1) gives a valid bound.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
DIAGNOSTIC = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_k20_d11.json"
ACTIVATION = REPO_ROOT / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts/min_state/"
    "s19_rm2sub_b_kernel_spectrum.json"
)
LIVE_SPECTRUM = REPO_ROOT / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
    "preaddmul_rm2sub_t128_s20/receipts/min_state/"
    "s19_rm2sub_a_spectrum_audit.json"
)
OUTPUT = WORKSTREAM / os.environ.get(
    "SPIN_Q1_OUTPUT",
    "golay_ba3_rm2sub_finite_B240_q1_outward_k20_d11.json",
)

B = 240
L = 8832
DISTANCE = 233_164
EPOCHS_PER_REGION = L // 128
STATE_DENOMINATOR = (1 << 19) - 1
GOLAY_SPECTRUM = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
LOWER_WEIGHT = int(os.environ.get("SPIN_LOWER_WEIGHT", "25"))
UPPER_WEIGHT = int(os.environ.get("SPIN_UPPER_WEIGHT", "215"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def up(value: float) -> float:
    if math.isnan(value):
        raise ArithmeticError("NaN in outward calculation")
    if value == math.inf:
        return value
    return math.nextafter(value, math.inf)


def down(value: float) -> float:
    if math.isnan(value):
        raise ArithmeticError("NaN in outward calculation")
    return math.nextafter(value, -math.inf)


def add_up(left: float, right: float) -> float:
    if left < 0.0 or right < 0.0:
        raise ArithmeticError("positive arithmetic received a negative addend")
    return up(left + right)


def mul_up(left: float, right: float) -> float:
    if left < 0.0 or right < 0.0:
        raise ArithmeticError("positive arithmetic received a negative factor")
    if left == 0.0 or right == 0.0:
        return 0.0
    return up(left * right)


def div_up(numerator: float, denominator: float) -> float:
    if numerator < 0.0 or denominator <= 0.0:
        raise ArithmeticError("invalid positive outward division")
    if numerator == 0.0:
        return 0.0
    return up(numerator / denominator)


def integer_upper(value: int) -> float:
    if value < 0:
        raise ArithmeticError("expected a nonnegative integer")
    if value == 0:
        return 0.0
    return up(float(value))


def integer_lower(value: int) -> float:
    if value <= 0:
        raise ArithmeticError("expected a positive integer denominator")
    return down(float(value))


def ratio_upper(numerator: int, denominator: int) -> float:
    if numerator == 0:
        return 0.0
    return div_up(integer_upper(numerator), integer_lower(denominator))


def power_up(base: float, exponent: int) -> float:
    if not 0.0 <= base < math.inf or exponent < 0:
        raise ArithmeticError("invalid outward integer power")
    result = 1.0
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = mul_up(result, factor)
        power >>= 1
        if power:
            factor = mul_up(factor, factor)
    return result


Matrix = tuple[float, float, float, float]
ZERO_MATRIX: Matrix = (0.0, 0.0, 0.0, 0.0)
IDENTITY: Matrix = (1.0, 0.0, 0.0, 1.0)


def matrix_add_up(left: Matrix, right: Matrix) -> Matrix:
    return tuple(add_up(a, b) for a, b in zip(left, right))  # type: ignore[return-value]


def matrix_scale_up(matrix: Matrix, scalar: float) -> Matrix:
    return tuple(mul_up(value, scalar) for value in matrix)  # type: ignore[return-value]


def matrix_divide_up(matrix: Matrix, denominator: int) -> Matrix:
    lower = integer_lower(denominator)
    return tuple(div_up(value, lower) for value in matrix)  # type: ignore[return-value]


def matrix_multiply_up(left: Matrix, right: Matrix) -> Matrix:
    a, b, c, d = left
    e, f, g, h = right
    return (
        add_up(mul_up(a, e), mul_up(b, g)),
        add_up(mul_up(a, f), mul_up(b, h)),
        add_up(mul_up(c, e), mul_up(d, g)),
        add_up(mul_up(c, f), mul_up(d, h)),
    )


def golay_direct_sum_counts() -> list[int]:
    counts = [1]
    for _ in range(B // 24):
        updated = [0] * (len(counts) + 24)
        for left_weight, left_count in enumerate(counts):
            for right_weight, right_count in GOLAY_SPECTRUM.items():
                updated[left_weight + right_weight] += left_count * right_count
        counts = updated
    if sum(counts) != 1 << (B // 2):
        raise ArithmeticError("Golay direct-sum mass check failed")
    return counts


def accumulator_transition_upper(input_weight: int, output_weight: int) -> float:
    if input_weight == 0:
        return 1.0 if output_weight == 0 else 0.0
    runs = (input_weight + 1) // 2
    if output_weight < runs or B - output_weight < input_weight - runs:
        return 0.0
    numerator = math.comb(output_weight - 1, runs - 1) * math.comb(
        B - output_weight, input_weight - runs
    )
    return ratio_upper(numerator, math.comb(B, input_weight))


def ba_expected_spectrum_upper() -> tuple[list[float], float, float]:
    transition = [
        [accumulator_transition_upper(u, h) for h in range(B + 1)]
        for u in range(B + 1)
    ]
    initial = [integer_upper(value) for value in golay_direct_sum_counts()]
    after_one = [0.0] * (B + 1)
    for u, multiplicity in enumerate(initial):
        if multiplicity == 0.0:
            continue
        for h, probability in enumerate(transition[u]):
            if probability:
                after_one[h] = add_up(
                    after_one[h], mul_up(multiplicity, probability)
                )
    after_two = [0.0] * (B + 1)
    for u, multiplicity in enumerate(after_one):
        if multiplicity == 0.0:
            continue
        for h, probability in enumerate(transition[u]):
            if probability:
                after_two[h] = add_up(
                    after_two[h], mul_up(multiplicity, probability)
                )

    tail = 0.0
    for weight in list(range(1, LOWER_WEIGHT)) + list(
        range(UPPER_WEIGHT + 1, B + 1)
    ):
        tail = add_up(tail, after_two[weight])
    if not tail < 1.0:
        raise ArithmeticError("Markov good-event lower bound vanished")
    good_probability_lower = down(1.0 - tail)
    conditioned = [0.0] * (B + 1)
    for weight in range(LOWER_WEIGHT, UPPER_WEIGHT + 1):
        conditioned[weight] = div_up(after_two[weight], good_probability_lower)
    return conditioned, tail, good_probability_lower


def load_constituent_data() -> tuple[int, int, list[tuple[int, int]]]:
    activation = json.loads(ACTIVATION.read_text(encoding="utf-8"))
    by_weight = {
        int(row["total_weight"]): row for row in activation["by_total_weight"]
    }
    for weight in (0, 1):
        row = by_weight[weight]
        if int(row["shell_size"]) != math.comb(128, weight):
            raise ArithmeticError("activation shell-size check failed")
    live = json.loads(LIVE_SPECTRUM.read_text(encoding="utf-8"))
    if not live["checks"]["minimum_distance"] == 48:
        raise ArithmeticError("live constituent distance check failed")
    spectrum = [
        (int(row["weight"]), int(row["count"]))
        for row in live["spectrum"]
        if int(row["weight"]) != 0
    ]
    if sum(count for _, count in spectrum) != STATE_DENOMINATOR:
        raise ArithmeticError("live constituent mass check failed")
    return (
        int(by_weight[1]["kernel_words"]),
        int(by_weight[1]["shell_size"]),
        spectrum,
    )


def support_average_upper(
    z: float, input_weight: int, spectrum: list[tuple[int, int]]
) -> float:
    total = 0.0
    shell_denominator = math.comb(128, input_weight)
    for codeword_weight, count in spectrum:
        kernel = 0.0
        minimum = max(0, input_weight - (128 - codeword_weight))
        maximum = min(input_weight, codeword_weight)
        for intersection in range(minimum, maximum + 1):
            placements = math.comb(codeword_weight, intersection) * math.comb(
                128 - codeword_weight, input_weight - intersection
            )
            probability = ratio_upper(placements, shell_denominator)
            exponent = input_weight + codeword_weight - 2 * intersection
            kernel = add_up(kernel, mul_up(probability, power_up(z, exponent)))
        live_probability = ratio_upper(count, STATE_DENOMINATOR)
        total = add_up(total, mul_up(live_probability, kernel))
    return total


def region_matrices(
    z: float, spectrum: list[tuple[int, int]], nonactivation_one: float
) -> tuple[Matrix, Matrix]:
    punctured = ratio_upper(STATE_DENOMINATOR, STATE_DENOMINATOR - 1)
    live_zero = support_average_upper(z, 0, spectrum)
    live_one = support_average_upper(z, 1, spectrum)
    zero_epoch: Matrix = (
        1.0,
        0.0,
        0.0,
        mul_up(punctured, live_zero),
    )
    averaged_one = mul_up(punctured, live_one)
    impulse_epoch: Matrix = (
        mul_up(nonactivation_one, z),
        z,
        div_up(averaged_one, integer_lower(STATE_DENOMINATOR)),
        averaged_one,
    )

    powers = [IDENTITY]
    for _ in range(EPOCHS_PER_REGION):
        powers.append(matrix_multiply_up(powers[-1], zero_epoch))
    inactive = powers[EPOCHS_PER_REGION]
    active = ZERO_MATRIX
    for epoch in range(EPOCHS_PER_REGION):
        term = matrix_multiply_up(
            matrix_multiply_up(powers[epoch], impulse_epoch),
            powers[EPOCHS_PER_REGION - 1 - epoch],
        )
        active = matrix_add_up(active, term)
    return inactive, matrix_divide_up(active, EPOCHS_PER_REGION)


def weight_conditioned_moments_upper(
    inactive: Matrix, active: Matrix
) -> list[float]:
    coefficients = [ZERO_MATRIX] * (B + 1)
    coefficients[0] = IDENTITY
    for completed in range(B):
        updated = [ZERO_MATRIX] * (B + 1)
        for weight in range(completed + 1):
            source = coefficients[weight]
            updated[weight] = matrix_add_up(
                updated[weight], matrix_multiply_up(source, inactive)
            )
            updated[weight + 1] = matrix_add_up(
                updated[weight + 1], matrix_multiply_up(source, active)
            )
        coefficients = updated
    moments = [0.0] * (B + 1)
    for weight in range(B + 1):
        row_sum = add_up(coefficients[weight][0], coefficients[weight][1])
        moments[weight] = div_up(row_sum, integer_lower(math.comb(B, weight)))
    return moments


def tilt_witnesses() -> dict[int, float]:
    payload = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    rows = payload["one_active"]["weight_rows"]
    result: dict[int, float] = {}
    for row in rows:
        weight = int(row["outer_weight"])
        if not LOWER_WEIGHT <= weight <= UPPER_WEIGHT:
            continue
        # Round the diagnostic log-surprisal to a 0.1 grid.  The calls to exp
        # merely choose a witness.  The proof treats the resulting binary64 z
        # as an exact rational and checks only 0 < z < 1.
        log_surprisal = round(float(row["best_log_surprisal"]) * 10.0) / 10.0
        z = math.exp(-math.exp(log_surprisal))
        if not 0.0 < z < 1.0:
            raise ArithmeticError("invalid Chernoff witness")
        result[weight] = z
    if set(result) != set(range(LOWER_WEIGHT, UPPER_WEIGHT + 1)):
        raise ArithmeticError("diagnostic receipt omitted a permitted weight")
    return result


def certify() -> dict[str, object]:
    if sys.float_info.radix != 2 or sys.float_info.mant_dig != 53:
        raise RuntimeError("the verifier requires IEEE-754 binary64 floats")
    if L % 128 or EPOCHS_PER_REGION != 69:
        raise ArithmeticError("finite parameter consistency check failed")

    conditioned_spectrum, tail_upper, good_lower = ba_expected_spectrum_upper()
    kernel_one, shell_one, live_spectrum = load_constituent_data()
    nonactivation_one = ratio_upper(kernel_one, shell_one)
    witnesses = tilt_witnesses()

    cached_moments: dict[str, list[float]] = {}
    for z in sorted(set(witnesses.values())):
        inactive, active = region_matrices(z, live_spectrum, nonactivation_one)
        cached_moments[z.hex()] = weight_conditioned_moments_upper(
            inactive, active
        )

    terms: list[dict[str, object]] = []
    aggregate = 0.0
    for weight in range(LOWER_WEIGHT, UPPER_WEIGHT + 1):
        z = witnesses[weight]
        moment = cached_moments[z.hex()][weight]
        inverse = div_up(1.0, z)
        chernoff = mul_up(moment, power_up(inverse, DISTANCE))
        probability = min(1.0, chernoff)
        term = mul_up(
            integer_upper(L),
            mul_up(conditioned_spectrum[weight], probability),
        )
        aggregate = add_up(aggregate, term)
        terms.append(
            {
                "outer_weight": weight,
                "z_binary64_hex_exact": z.hex(),
                "conditioned_expected_multiplicity_upper_hex": (
                    conditioned_spectrum[weight].hex()
                ),
                "inner_probability_upper_hex": probability.hex(),
                "contribution_upper_hex": term.hex(),
            }
        )

    certified_integer_margin = 0
    while aggregate <= math.ldexp(1.0, -(certified_integer_margin + 1)):
        certified_integer_margin += 1
    dominant = sorted(
        terms,
        key=lambda row: float.fromhex(str(row["contribution_upper_hex"])),
        reverse=True,
    )[:20]
    return {
        "schema": "golay-ba3-rm2sub-finite-b240-q1-outward-v1",
        "status": "OUTWARD_CERTIFICATE",
        "claim": {
            "probability_space": (
                "One Golay--BA-3 draw per outer row, independently conditioned "
                "on G_240; independent row-to-region permutations and RM2Sub "
                "multipliers; full parent message space before shortening."
            ),
            "occupation": 1,
            "bad_weight_at_most": DISTANCE,
            "expected_bad_word_count_upper_hex": aggregate.hex(),
            "certified_integer_margin_bits": certified_integer_margin,
            "comparison_to_2^-40": aggregate <= math.ldexp(1.0, -40),
            "comparison_to_2^-47": aggregate <= math.ldexp(1.0, -47),
            "display_margin_bits_not_used_by_verifier": -math.log2(aggregate),
        },
        "parameters": {
            "message_bits_after_shortening": 1 << 20,
            "parent_dimension": B * L // 2,
            "shortened_input_coordinates": B * L // 2 - (1 << 20),
            "outer_bits": B,
            "outer_rows": L,
            "parent_output_bits": B * L,
            "inner_epochs": B * L // 128,
            "epochs_per_region": EPOCHS_PER_REGION,
            "permitted_outer_weights": [LOWER_WEIGHT, UPPER_WEIGHT],
        },
        "outer_conditioning": {
            "expected_tail_word_count_upper_hex": tail_upper.hex(),
            "good_event_probability_lower_hex": good_lower.hex(),
            "expected_trials_upper_display": up(1.0 / good_lower),
        },
        "arithmetic": {
            "format": "IEEE-754 binary64",
            "proof_rule": (
                "All nonnegative proof-path operations are rounded one ULP "
                "toward +infinity; exact integer denominators are rounded "
                "toward -infinity before division."
            ),
            "underflow_rule": (
                "A positive operation that rounds to zero is advanced to the "
                "least positive subnormal by the same upward rule."
            ),
            "transcendental_scope": (
                "exp selects legal binary64 Chernoff witnesses only; no "
                "transcendental result enters as a claimed bound."
            ),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "unique_tilts": len(cached_moments),
        },
        "dependencies": [
            {"path": str(DIAGNOSTIC.relative_to(REPO_ROOT)), "sha256": sha256(DIAGNOSTIC)},
            {"path": str(ACTIVATION.relative_to(REPO_ROOT)), "sha256": sha256(ACTIVATION)},
            {"path": str(LIVE_SPECTRUM.relative_to(REPO_ROOT)), "sha256": sha256(LIVE_SPECTRUM)},
        ],
        "dominant_contributions": dominant,
        "weight_witnesses": terms,
        "limitations": [
            "This receipt covers occupation Q=1 only.",
            "The conditioned setup law is specified mathematically; an efficient G_240 decision procedure remains open.",
            "The activation and live-spectrum source claims are hash-bound but require the separate RM2Sub interface audit.",
        ],
    }


def main() -> None:
    payload = certify()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"receipt={OUTPUT}")


if __name__ == "__main__":
    main()
