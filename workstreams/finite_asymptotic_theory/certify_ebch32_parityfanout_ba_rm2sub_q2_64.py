#!/usr/bin/env python3
"""Outward q=2..64 certificate for EBCH32--PF31x33--BA/RM2Sub-S19."""

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
    conditioned_spectrum_upper,
    down,
    integer_lower,
    integer_upper,
    sha256,
    up,
)
from certify_ebch32_parityfanout_ba_rm2sub_q1 import (
    ACTIVATION,
    DISTANCE,
    FROZEN_MANIFEST,
    LIVE_SPECTRUM,
    L,
    REPO_ROOT,
    STEP_BITS,
)
from certify_golay_ba_rm2sub_finite_q2_64 import (
    SCALED_IDENTITY,
    SCALED_ONE,
    SCALED_ZERO,
    Matrix,
    Scaled,
    candidate_epoch_matrices,
    impulse_matrices_upper,
    load_activation_and_live,
    mul_down,
    power_down,
    ratio_upper,
    scaled_add_up,
    scaled_at_most_power_of_two,
    scaled_from_float,
    scaled_hex,
    scaled_matrix_from_float,
    scaled_matrix_multiply_up,
    scaled_multiply_up,
    scaled_power_up,
)
from certify_golay_ba_rm2sub_finite_one_active import (
    ZERO_MATRIX,
    matrix_add_up,
    matrix_multiply_up,
    matrix_scale_up,
)


DIAGNOSTIC = WORKSTREAM / "ebch32_parityfanout31x33_ba3_rm2sub_B256_q2_64_d11_diagnostic.json"
SETUP_VERIFIER = WORKSTREAM / "certify_ebch32_parityfanout_ba_setup.py"
SETUP_RECEIPT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json"
Q1_VERIFIER = WORKSTREAM / "certify_ebch32_parityfanout_ba_rm2sub_q1.py"
GENERIC_SCALED_VERIFIER = WORKSTREAM / "certify_golay_ba_rm2sub_finite_q2_64.py"
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_rm2sub_B256_q2_64_outward_d11.json"
EPOCHS_PER_REGION = L // STEP_BITS


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
            mul_down(power_down(p, weight), power_down(q_lower, B - weight)),
        )
        candidate = ratio_from_float_upper(conditioned_spectrum[weight], denominator)
        if candidate > best:
            best = candidate
            best_weight = weight
    return best, best_weight


def ratio_from_float_upper(numerator: float, denominator: float) -> float:
    if numerator < 0.0 or denominator <= 0.0:
        raise ArithmeticError("invalid positive ratio")
    return 0.0 if numerator == 0.0 else up(numerator / denominator)


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
        raise ArithmeticError("diagnostic omitted an occupation")
    return result


def certify() -> dict[str, object]:
    if sys.float_info.radix != 2 or sys.float_info.mant_dig != 53:
        raise RuntimeError("the verifier requires IEEE-754 binary64")
    if EPOCHS_PER_REGION != 64 or B * L != 1 << 21:
        raise ArithmeticError("finite geometry check failed")

    conditioned, tail_upper, good_lower = conditioned_spectrum_upper()
    activation, live = load_activation_and_live()
    witnesses = diagnostic_witnesses()

    envelopes: dict[str, tuple[float, int]] = {}
    caches: dict[str, list[Scaled]] = {}
    for p, z in sorted(set(witnesses.values())):
        if p.hex() not in envelopes:
            envelopes[p.hex()] = bernoulli_envelope_upper(conditioned, p)
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
            scaled_power_up(scaled_from_float(ratio_from_float_upper(1.0, z)), DISTANCE),
        )
        inner = chernoff if scaled_at_most_power_of_two(chernoff, 0) else SCALED_ONE
        contribution = scaled_multiply_up(
            scaled_from_float(integer_upper(math.comb(L, occupation))),
            scaled_multiply_up(
                scaled_power_up(scaled_from_float(envelope), occupation), inner
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
            float.fromhex(str(row["contribution_upper_scaled"]["mantissa_binary64_hex"])),
        ),
        reverse=True,
    )[:20]
    return {
        "schema": "ebch32-parityfanout31x33-ba3-rm2sub-b256-q2-64-outward-v1",
        "status": "OUTWARD_CERTIFICATE",
        "claim": {
            "probability_space": (
                "Independent EBCH32--ParityFanout31x33--BA-3 row draws, each "
                "conditioned on G_256; independent region permutations and "
                "RM2Sub-S19 multipliers; full parent message space."
            ),
            "occupations": [2, 64],
            "bad_weight_at_most": DISTANCE,
            "expected_bad_word_count_upper_scaled": scaled_hex(aggregate),
            "certified_integer_margin_bits": integer_margin,
            "comparison_to_2^-40": scaled_at_most_power_of_two(aggregate, 40),
            "display_margin_bits_not_used_by_verifier": -math.log2(aggregate[0]) - aggregate[1],
        },
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": B * L,
            "distance": DISTANCE,
            "epochs_per_region": EPOCHS_PER_REGION,
            "permitted_outer_weights": [LOWER_WEIGHT, UPPER_WEIGHT],
        },
        "outer_conditioning": {
            "expected_tail_word_count_upper_hex": tail_upper.hex(),
            "good_event_probability_lower_hex": good_lower.hex(),
        },
        "arithmetic": {
            "format": "IEEE-754 binary64 with scaled nonnegative products",
            "proof_rule": (
                "All positive upper-bound operations advance one ULP toward "
                "+infinity; denominators and Bernoulli masses used below a "
                "division are rounded toward -infinity."
            ),
            "transcendental_scope": "exp selects exact-binary64 witnesses only.",
            "unique_probability_tilt_pairs": len(caches),
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "dependencies": [
            {"path": str(DIAGNOSTIC.relative_to(REPO_ROOT)), "sha256": sha256(DIAGNOSTIC)},
            {"path": str(SETUP_VERIFIER.relative_to(REPO_ROOT)), "sha256": sha256(SETUP_VERIFIER)},
            {"path": str(SETUP_RECEIPT.relative_to(REPO_ROOT)), "sha256": sha256(SETUP_RECEIPT)},
            {"path": str(Q1_VERIFIER.relative_to(REPO_ROOT)), "sha256": sha256(Q1_VERIFIER)},
            {"path": str(GENERIC_SCALED_VERIFIER.relative_to(REPO_ROOT)), "sha256": sha256(GENERIC_SCALED_VERIFIER)},
            {"path": str(ACTIVATION.relative_to(REPO_ROOT)), "sha256": sha256(ACTIVATION)},
            {"path": str(LIVE_SPECTRUM.relative_to(REPO_ROOT)), "sha256": sha256(LIVE_SPECTRUM)},
            {"path": str(FROZEN_MANIFEST.relative_to(REPO_ROOT)), "sha256": sha256(FROZEN_MANIFEST)},
        ],
        "dominant_contributions": dominant,
        "occupation_rows": rows,
        "limitations": [
            "This receipt covers occupations Q=2 through 64 only.",
            "An efficient G_256 setup test remains open.",
            "The frozen implementation does not yet implement this EBCH32--BA outer.",
            "The full RM2Sub interface-equivalence audit remains open.",
        ],
    }


def main() -> None:
    payload = certify()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"receipt={OUTPUT}")


if __name__ == "__main__":
    main()
