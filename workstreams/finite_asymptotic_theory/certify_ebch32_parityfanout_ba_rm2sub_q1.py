#!/usr/bin/env python3
"""Outward q=1 certificate for EBCH32--PF31x33--BA/RM2Sub-S19."""

from __future__ import annotations

import json
import math
from pathlib import Path
import platform
import sys

from certify_ebch32_parityfanout_ba_setup import (
    B,
    LOWER_WEIGHT,
    UPPER_WEIGHT,
    WORKSTREAM,
    add_up,
    conditioned_spectrum_upper,
    div_up,
    integer_lower,
    integer_upper,
    mul_up,
    ratio_upper,
    sha256,
    up,
)


REPO_ROOT = WORKSTREAM.parents[1]
DIAGNOSTIC = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_q1_d11_diagnostic.json"
SETUP_VERIFIER = WORKSTREAM / "certify_ebch32_parityfanout_ba_setup.py"
SETUP_RECEIPT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json"
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
FROZEN_MANIFEST = REPO_ROOT / (
    "constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_"
    "splitstate_preaddmul_rm2sub_t128_s19/frozen_source/SOURCE_MANIFEST.json"
)
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_rm2sub_B256_q1_outward_d11.json"

L = 8192
DISTANCE = 230_686
STEP_BITS = 128
EPOCHS_PER_REGION = L // STEP_BITS
STATE_DENOMINATOR = (1 << 19) - 1

Matrix = tuple[float, float, float, float]
ZERO_MATRIX: Matrix = (0.0, 0.0, 0.0, 0.0)
IDENTITY: Matrix = (1.0, 0.0, 0.0, 1.0)


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


def matrix_add_up(left: Matrix, right: Matrix) -> Matrix:
    return tuple(add_up(a, b) for a, b in zip(left, right))  # type: ignore[return-value]


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


def load_constituent_data() -> tuple[int, int, list[tuple[int, int]]]:
    activation = json.loads(ACTIVATION.read_text(encoding="utf-8"))
    by_weight = {
        int(row["total_weight"]): row for row in activation["by_total_weight"]
    }
    for weight in (0, 1):
        if int(by_weight[weight]["shell_size"]) != math.comb(STEP_BITS, weight):
            raise ArithmeticError("activation shell-size check failed")
    live = json.loads(LIVE_SPECTRUM.read_text(encoding="utf-8"))
    if live["checks"]["minimum_distance"] != 48:
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
    shell_denominator = math.comb(STEP_BITS, input_weight)
    z_powers = [1.0]
    for _ in range(STEP_BITS):
        z_powers.append(mul_up(z_powers[-1], z))
    for codeword_weight, count in spectrum:
        kernel = 0.0
        minimum = max(0, input_weight - (STEP_BITS - codeword_weight))
        maximum = min(input_weight, codeword_weight)
        for intersection in range(minimum, maximum + 1):
            placements = math.comb(codeword_weight, intersection) * math.comb(
                STEP_BITS - codeword_weight, input_weight - intersection
            )
            exponent = input_weight + codeword_weight - 2 * intersection
            kernel = add_up(
                kernel,
                mul_up(
                    ratio_upper(placements, shell_denominator),
                    z_powers[exponent],
                ),
            )
        total = add_up(
            total,
            mul_up(ratio_upper(count, STATE_DENOMINATOR), kernel),
        )
    return total


def region_matrices(
    z: float, spectrum: list[tuple[int, int]], nonactivation_one: float
) -> tuple[Matrix, Matrix]:
    punctured = ratio_upper(STATE_DENOMINATOR, STATE_DENOMINATOR - 1)
    live_zero = support_average_upper(z, 0, spectrum)
    live_one = support_average_upper(z, 1, spectrum)
    zero_epoch: Matrix = (1.0, 0.0, 0.0, mul_up(punctured, live_zero))
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
    rows = payload["conditioned_q1"]["weight_rows"]
    result: dict[int, float] = {}
    for row in rows:
        weight = int(row["outer_weight"])
        log_surprisal = round(float(row["best_log_surprisal"]) * 10.0) / 10.0
        z = math.exp(-math.exp(log_surprisal))
        if not 0.0 < z < 1.0:
            raise ArithmeticError("invalid Chernoff witness")
        result[weight] = z
    if set(result) != set(range(LOWER_WEIGHT, UPPER_WEIGHT + 1)):
        raise ArithmeticError("diagnostic omitted a permitted weight")
    return result


def certify() -> dict[str, object]:
    if sys.float_info.radix != 2 or sys.float_info.mant_dig != 53:
        raise RuntimeError("the verifier requires IEEE-754 binary64")
    if EPOCHS_PER_REGION != 64 or B * L != 1 << 21:
        raise ArithmeticError("finite geometry check failed")

    conditioned, tail_upper, good_lower = conditioned_spectrum_upper()
    kernel_one, shell_one, live = load_constituent_data()
    nonactivation_one = ratio_upper(kernel_one, shell_one)
    witnesses = tilt_witnesses()

    cache: dict[str, list[float]] = {}
    for z in sorted(set(witnesses.values())):
        inactive, active = region_matrices(z, live, nonactivation_one)
        cache[z.hex()] = weight_conditioned_moments_upper(inactive, active)

    aggregate = 0.0
    terms = []
    for weight in range(LOWER_WEIGHT, UPPER_WEIGHT + 1):
        z = witnesses[weight]
        moment = cache[z.hex()][weight]
        chernoff = mul_up(moment, power_up(div_up(1.0, z), DISTANCE))
        probability = min(1.0, chernoff)
        term = mul_up(integer_upper(L), mul_up(conditioned[weight], probability))
        aggregate = add_up(aggregate, term)
        terms.append(
            {
                "outer_weight": weight,
                "z_binary64_hex_exact": z.hex(),
                "conditioned_expected_multiplicity_upper_hex": conditioned[weight].hex(),
                "inner_probability_upper_hex": probability.hex(),
                "contribution_upper_hex": term.hex(),
            }
        )

    integer_margin = 0
    while aggregate <= math.ldexp(1.0, -(integer_margin + 1)):
        integer_margin += 1
    dominant = sorted(
        terms,
        key=lambda row: float.fromhex(str(row["contribution_upper_hex"])),
        reverse=True,
    )[:20]
    return {
        "schema": "ebch32-parityfanout31x33-ba3-rm2sub-b256-q1-outward-v1",
        "status": "OUTWARD_CERTIFICATE",
        "claim": {
            "probability_space": (
                "Independent EBCH32--ParityFanout31x33--BA-3 row draws, each "
                "conditioned on G_256; independent region permutations and "
                "RM2Sub-S19 multipliers; full parent message space."
            ),
            "occupation": 1,
            "bad_weight_at_most": DISTANCE,
            "expected_bad_word_count_upper_hex": aggregate.hex(),
            "certified_integer_margin_bits": integer_margin,
            "comparison_to_2^-40": aggregate <= math.ldexp(1.0, -40),
            "display_margin_bits_not_used_by_verifier": -math.log2(aggregate),
        },
        "parameters": {
            "message_bits": 1 << 20,
            "parent_dimension": B * L // 2,
            "shortened_input_coordinates": 0,
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": B * L,
            "distance": DISTANCE,
            "inner_epochs": B * L // STEP_BITS,
            "epochs_per_region": EPOCHS_PER_REGION,
            "permitted_outer_weights": [LOWER_WEIGHT, UPPER_WEIGHT],
        },
        "outer_conditioning": {
            "expected_tail_word_count_upper_hex": tail_upper.hex(),
            "good_event_probability_lower_hex": good_lower.hex(),
        },
        "arithmetic": {
            "format": "IEEE-754 binary64",
            "proof_rule": (
                "All nonnegative proof-path operations are advanced one ULP "
                "toward +infinity; exact integer denominators are rounded "
                "toward -infinity before division."
            ),
            "transcendental_scope": (
                "exp selects legal binary64 Chernoff witnesses only; its "
                "value is treated as an exact rational in the proof path."
            ),
            "unique_tilts": len(cache),
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "dependencies": [
            {"path": str(DIAGNOSTIC.relative_to(REPO_ROOT)), "sha256": sha256(DIAGNOSTIC)},
            {"path": str(SETUP_VERIFIER.relative_to(REPO_ROOT)), "sha256": sha256(SETUP_VERIFIER)},
            {"path": str(SETUP_RECEIPT.relative_to(REPO_ROOT)), "sha256": sha256(SETUP_RECEIPT)},
            {"path": str(ACTIVATION.relative_to(REPO_ROOT)), "sha256": sha256(ACTIVATION)},
            {"path": str(LIVE_SPECTRUM.relative_to(REPO_ROOT)), "sha256": sha256(LIVE_SPECTRUM)},
            {"path": str(FROZEN_MANIFEST.relative_to(REPO_ROOT)), "sha256": sha256(FROZEN_MANIFEST)},
        ],
        "dominant_contributions": dominant,
        "weight_witnesses": terms,
        "limitations": [
            "This receipt covers occupation Q=1 only.",
            "An efficient G_256 setup test remains open.",
            "The frozen implementation does not yet implement this EBCH32--BA outer.",
            "The RM2Sub source receipts are hash-bound but the full interface-equivalence audit remains open.",
        ],
    }


def main() -> None:
    payload = certify()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"receipt={OUTPUT}")


if __name__ == "__main__":
    main()
