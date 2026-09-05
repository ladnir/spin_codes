#!/usr/bin/env python3
"""Outward certificate for finite B=240 occupations 2 through 64."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys

from certify_golay_ba_rm2sub_finite_one_active import (
    ACTIVATION,
    B,
    DISTANCE,
    IDENTITY,
    LIVE_SPECTRUM,
    L,
    LOWER_WEIGHT,
    Matrix,
    REPO_ROOT,
    STATE_DENOMINATOR,
    UPPER_WEIGHT,
    WORKSTREAM,
    ZERO_MATRIX,
    add_up,
    ba_expected_spectrum_upper,
    div_up,
    down,
    integer_lower,
    integer_upper,
    matrix_add_up,
    matrix_multiply_up,
    matrix_scale_up,
    mul_up,
    power_up,
    ratio_upper,
    sha256,
    up,
)


DIAGNOSTIC = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_bulk_q2_64_k20_d11.json"
Q1_VERIFIER = WORKSTREAM / "certify_golay_ba_rm2sub_finite_one_active.py"
OUTPUT = WORKSTREAM / os.environ.get(
    "SPIN_Q2_64_OUTPUT",
    "golay_ba3_rm2sub_finite_B240_q2_64_outward_k20_d11.json",
)
STEP_BITS = 128
EPOCHS_PER_REGION = L // STEP_BITS


# A scaled nonnegative upper bound represents mantissa * 2**exponent.  A
# nonzero mantissa is normalized to [1/2,1).  This prevents the 240-region
# product from underflowing before it is combined with the outer envelope.
Scaled = tuple[float, int]
SCALED_ZERO: Scaled = (0.0, 0)
SCALED_ONE: Scaled = (0.5, 1)
ScaledMatrix = tuple[Scaled, Scaled, Scaled, Scaled]
SCALED_IDENTITY: ScaledMatrix = (
    SCALED_ONE,
    SCALED_ZERO,
    SCALED_ZERO,
    SCALED_ONE,
)


def scaled_normalize(mantissa: float, exponent: int) -> Scaled:
    if mantissa == 0.0:
        return SCALED_ZERO
    if mantissa < 0.0 or not math.isfinite(mantissa):
        raise ArithmeticError("invalid scaled upper bound")
    normalized, shift = math.frexp(mantissa)
    return normalized, exponent + shift


def scaled_from_float(value: float) -> Scaled:
    if value == 0.0:
        return SCALED_ZERO
    if value < 0.0 or not math.isfinite(value):
        raise ArithmeticError("invalid float-to-scaled conversion")
    return math.frexp(value)


def scaled_add_up(left: Scaled, right: Scaled) -> Scaled:
    if left[0] == 0.0:
        return right
    if right[0] == 0.0:
        return left
    exponent = max(left[1], right[1])
    left_mantissa = up(math.ldexp(left[0], left[1] - exponent))
    right_mantissa = up(math.ldexp(right[0], right[1] - exponent))
    return scaled_normalize(add_up(left_mantissa, right_mantissa), exponent)


def scaled_multiply_up(left: Scaled, right: Scaled) -> Scaled:
    if left[0] == 0.0 or right[0] == 0.0:
        return SCALED_ZERO
    return scaled_normalize(
        mul_up(left[0], right[0]), left[1] + right[1]
    )


def scaled_power_up(base: Scaled, exponent: int) -> Scaled:
    if exponent < 0:
        raise ArithmeticError("negative scaled exponent")
    result = SCALED_ONE
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = scaled_multiply_up(result, factor)
        power >>= 1
        if power:
            factor = scaled_multiply_up(factor, factor)
    return result


def scaled_matrix_from_float(matrix: Matrix) -> ScaledMatrix:
    return tuple(scaled_from_float(value) for value in matrix)  # type: ignore[return-value]


def scaled_matrix_multiply_up(
    left: ScaledMatrix, right: ScaledMatrix
) -> ScaledMatrix:
    a, b, c, d = left
    e, f, g, h = right
    return (
        scaled_add_up(scaled_multiply_up(a, e), scaled_multiply_up(b, g)),
        scaled_add_up(scaled_multiply_up(a, f), scaled_multiply_up(b, h)),
        scaled_add_up(scaled_multiply_up(c, e), scaled_multiply_up(d, g)),
        scaled_add_up(scaled_multiply_up(c, f), scaled_multiply_up(d, h)),
    )


def scaled_at_most_power_of_two(value: Scaled, negative_exponent: int) -> bool:
    threshold = (0.5, 1 - negative_exponent)
    if value[0] == 0.0:
        return True
    if value[1] != threshold[1]:
        return value[1] < threshold[1]
    return value[0] <= threshold[0]


def scaled_hex(value: Scaled) -> dict[str, object]:
    return {
        "mantissa_binary64_hex": value[0].hex(),
        "binary_exponent": value[1],
    }


def mul_down(left: float, right: float) -> float:
    if left < 0.0 or right < 0.0:
        raise ArithmeticError("positive arithmetic received a negative factor")
    if left == 0.0 or right == 0.0:
        return 0.0
    return down(left * right)


def power_down(base: float, exponent: int) -> float:
    if not 0.0 < base < math.inf or exponent < 0:
        raise ArithmeticError("invalid downward integer power")
    result = 1.0
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = mul_down(result, factor)
        power >>= 1
        if power:
            factor = mul_down(factor, factor)
    return result


def load_activation_and_live() -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    activation_payload = json.loads(ACTIVATION.read_text(encoding="utf-8"))
    activation = [(0, 0)] * (STEP_BITS + 1)
    for row in activation_payload["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight <= STEP_BITS:
            shell = int(row["shell_size"])
            kernel = int(row["kernel_words"])
            if shell != math.comb(STEP_BITS, weight):
                raise ArithmeticError("activation shell-size check failed")
            activation[weight] = (kernel, shell)
    live_payload = json.loads(LIVE_SPECTRUM.read_text(encoding="utf-8"))
    if live_payload["checks"]["minimum_distance"] != 48:
        raise ArithmeticError("live constituent distance check failed")
    live = [
        (int(row["weight"]), int(row["count"]))
        for row in live_payload["spectrum"]
        if int(row["weight"]) != 0
    ]
    if sum(count for _, count in live) != STATE_DENOMINATOR:
        raise ArithmeticError("live spectrum mass check failed")
    return activation, live


def impulse_matrices_upper(
    z: float,
    activation: list[tuple[int, int]],
    live: list[tuple[int, int]],
) -> list[Matrix]:
    z_powers = [1.0]
    for _ in range(STEP_BITS):
        z_powers.append(mul_up(z_powers[-1], z))
    punctured = ratio_upper(STATE_DENOMINATOR, STATE_DENOMINATOR - 1)
    matrices: list[Matrix] = []
    for input_weight in range(STEP_BITS + 1):
        support_average = 0.0
        shell_denominator = math.comb(STEP_BITS, input_weight)
        for codeword_weight, count in live:
            kernel = 0.0
            minimum = max(0, input_weight - (STEP_BITS - codeword_weight))
            maximum = min(input_weight, codeword_weight)
            for intersection in range(minimum, maximum + 1):
                placements = math.comb(codeword_weight, intersection) * math.comb(
                    STEP_BITS - codeword_weight,
                    input_weight - intersection,
                )
                exponent = input_weight + codeword_weight - 2 * intersection
                kernel = add_up(
                    kernel,
                    mul_up(
                        ratio_upper(placements, shell_denominator),
                        z_powers[exponent],
                    ),
                )
            support_average = add_up(
                support_average,
                mul_up(ratio_upper(count, STATE_DENOMINATOR), kernel),
            )
        if input_weight == 0:
            matrices.append(
                (1.0, 0.0, 0.0, mul_up(punctured, support_average))
            )
            continue
        nonactivation = ratio_upper(*activation[input_weight])
        output = z_powers[input_weight]
        averaged = mul_up(punctured, support_average)
        matrices.append(
            (
                mul_up(nonactivation, output),
                output,
                div_up(averaged, integer_lower(STATE_DENOMINATOR)),
                averaged,
            )
        )
    return matrices


def candidate_epoch_matrices(impulses: list[Matrix], p: float) -> list[Matrix]:
    if not 0.0 < p < 1.0:
        raise ArithmeticError("invalid Bernoulli witness")
    one_minus_p_upper = up(1.0 - p)
    p_powers = [1.0]
    q_powers = [1.0]
    for _ in range(STEP_BITS):
        p_powers.append(mul_up(p_powers[-1], p))
        q_powers.append(mul_up(q_powers[-1], one_minus_p_upper))
    result: list[Matrix] = []
    for candidates in range(STEP_BITS + 1):
        matrix = ZERO_MATRIX
        for active in range(candidates + 1):
            probability = mul_up(
                integer_upper(math.comb(candidates, active)),
                mul_up(
                    p_powers[active], q_powers[candidates - active]
                ),
            )
            matrix = matrix_add_up(
                matrix, matrix_scale_up(impulses[active], probability)
            )
        result.append(matrix)
    return result


def region_matrices_upper(candidates: list[Matrix], maximum_q: int) -> list[Matrix]:
    current = candidates[: maximum_q + 1]
    current_maximum = min(STEP_BITS, maximum_q)
    for completed_epochs in range(1, EPOCHS_PER_REGION):
        next_maximum = min((completed_epochs + 1) * STEP_BITS, maximum_q)
        updated = [ZERO_MATRIX] * (next_maximum + 1)
        for next_count in range(min(STEP_BITS, next_maximum) + 1):
            maximum_source = min(current_maximum, next_maximum - next_count)
            for source_count in range(maximum_source + 1):
                destination = source_count + next_count
                numerator = math.comb(STEP_BITS, next_count) * math.comb(
                    completed_epochs * STEP_BITS, source_count
                )
                denominator = math.comb(
                    (completed_epochs + 1) * STEP_BITS, destination
                )
                weight = ratio_upper(numerator, denominator)
                product = matrix_multiply_up(
                    current[source_count], candidates[next_count]
                )
                updated[destination] = matrix_add_up(
                    updated[destination], matrix_scale_up(product, weight)
                )
        current = updated
        current_maximum = next_maximum
    return current


def total_moments_upper(regions: list[Matrix]) -> list[Scaled]:
    moments = []
    for region in regions:
        product = SCALED_IDENTITY
        scaled_region = scaled_matrix_from_float(region)
        for _ in range(B):
            product = scaled_matrix_multiply_up(product, scaled_region)
        moments.append(scaled_add_up(product[0], product[1]))
    return moments


def bernoulli_envelope_upper(
    conditioned_spectrum: list[float], p: float
) -> tuple[float, int]:
    q_lower = down(1.0 - p)
    best = 0.0
    best_weight = -1
    for weight in range(LOWER_WEIGHT, UPPER_WEIGHT + 1):
        denominator = mul_down(
            integer_lower(math.comb(B, weight)),
            mul_down(
                power_down(p, weight), power_down(q_lower, B - weight)
            ),
        )
        candidate = div_up(conditioned_spectrum[weight], denominator)
        if candidate > best:
            best = candidate
            best_weight = weight
    return best, best_weight


def diagnostic_witnesses() -> dict[int, tuple[float, float]]:
    payload = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    result = {}
    for row in payload["occupation_rows"]:
        occupation = int(row["active_regular_outer_blocks"])
        if not 2 <= occupation <= 64:
            continue
        p = float(row["best_candidate_probability"])
        z = math.exp(-math.exp(float(row["best_log_surprisal"])))
        if not 0.0 < p < 1.0 or not 0.0 < z < 1.0:
            raise ArithmeticError("invalid diagnostic witness")
        result[occupation] = (p, z)
    if set(result) != set(range(2, 65)):
        raise ArithmeticError("diagnostic receipt omitted an occupation")
    return result


def certify() -> dict[str, object]:
    if sys.float_info.radix != 2 or sys.float_info.mant_dig != 53:
        raise RuntimeError("the verifier requires IEEE-754 binary64 floats")
    conditioned_spectrum, tail_upper, good_lower = ba_expected_spectrum_upper()
    activation, live = load_activation_and_live()
    witnesses = diagnostic_witnesses()

    envelopes: dict[str, tuple[float, int]] = {}
    caches: dict[str, list[Scaled]] = {}
    for p, z in sorted(set(witnesses.values())):
        p_key = p.hex()
        if p_key not in envelopes:
            envelopes[p_key] = bernoulli_envelope_upper(
                conditioned_spectrum, p
            )
        impulses = impulse_matrices_upper(z, activation, live)
        candidates = candidate_epoch_matrices(impulses, p)
        regions = region_matrices_upper(candidates, 64)
        caches[f"{p.hex()}|{z.hex()}"] = total_moments_upper(regions)

    rows = []
    aggregate = SCALED_ZERO
    for occupation in range(2, 65):
        p, z = witnesses[occupation]
        envelope, maximizing_weight = envelopes[p.hex()]
        moment = caches[f"{p.hex()}|{z.hex()}"][occupation]
        chernoff = scaled_multiply_up(
            moment,
            scaled_power_up(
                scaled_from_float(div_up(1.0, z)), DISTANCE
            ),
        )
        inner = chernoff if scaled_at_most_power_of_two(chernoff, 0) else SCALED_ONE
        contribution = scaled_multiply_up(
            scaled_from_float(integer_upper(math.comb(L, occupation))),
            scaled_multiply_up(
                scaled_power_up(scaled_from_float(envelope), occupation),
                inner,
            ),
        )
        aggregate = scaled_add_up(aggregate, contribution)
        rows.append(
            {
                "active_outer_rows": occupation,
                "candidate_probability_binary64_hex_exact": p.hex(),
                "z_binary64_hex_exact": z.hex(),
                "bernoulli_envelope_upper_hex": envelope.hex(),
                "envelope_maximizing_weight": maximizing_weight,
                "inner_probability_upper_scaled": scaled_hex(inner),
                "contribution_upper_scaled": scaled_hex(contribution),
            }
        )

    integer_margin = 0
    while scaled_at_most_power_of_two(aggregate, integer_margin + 1):
        integer_margin += 1
    dominant = sorted(
        rows,
        key=lambda row: (
            int(row["contribution_upper_scaled"]["binary_exponent"]),
            float.fromhex(
                str(row["contribution_upper_scaled"]["mantissa_binary64_hex"])
            ),
        ),
        reverse=True,
    )[:20]
    return {
        "schema": "golay-ba3-rm2sub-finite-b240-q2-64-outward-v1",
        "status": "OUTWARD_CERTIFICATE",
        "claim": {
            "probability_space": (
                "Independent Golay--BA-3 row draws, each conditioned on "
                "G_240; independent route and region permutations and "
                "RM2Sub multipliers; full parent space before shortening."
            ),
            "occupations": [2, 64],
            "bad_weight_at_most": DISTANCE,
            "expected_bad_word_count_upper_scaled": scaled_hex(aggregate),
            "certified_integer_margin_bits": integer_margin,
            "comparison_to_2^-40": scaled_at_most_power_of_two(aggregate, 40),
            "comparison_to_2^-80": scaled_at_most_power_of_two(aggregate, 80),
            "display_margin_bits_not_used_by_verifier": (
                -math.log2(aggregate[0]) - aggregate[1]
            ),
        },
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "parent_output_bits": B * L,
            "distance": DISTANCE,
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
                "All positive upper-bound operations advance one ULP toward "
                "+infinity. Denominators and Bernoulli masses used below a "
                "division are rounded toward -infinity."
            ),
            "transcendental_scope": (
                "exp selects legal binary64 witnesses only and is not used "
                "to assert a numerical bound."
            ),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "unique_probability_tilt_pairs": len(caches),
        },
        "dependencies": [
            {"path": str(DIAGNOSTIC.relative_to(REPO_ROOT)), "sha256": sha256(DIAGNOSTIC)},
            {"path": str(Q1_VERIFIER.relative_to(REPO_ROOT)), "sha256": sha256(Q1_VERIFIER)},
            {"path": str(ACTIVATION.relative_to(REPO_ROOT)), "sha256": sha256(ACTIVATION)},
            {"path": str(LIVE_SPECTRUM.relative_to(REPO_ROOT)), "sha256": sha256(LIVE_SPECTRUM)},
        ],
        "dominant_contributions": dominant,
        "occupation_rows": rows,
        "limitations": [
            "This receipt covers only occupations 2 through 64.",
            "Efficient G_240 conditioning remains open.",
            "The separate RM2Sub implementation-equivalence audit remains open.",
        ],
    }


def main() -> None:
    payload = certify()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"receipt={OUTPUT}")


if __name__ == "__main__":
    main()
