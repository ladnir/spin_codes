#!/usr/bin/env python3
"""Outward mean spectrum and conditional cap event for Block Expand-5.

The regional expander spectrum is computed from exact integer occupancy
counts.  Five accumulator stages use exact binomial transition ratios.  Arb
encloses only divisions, additions, and logarithms.  The cap-event conclusion
remains conditional on Var(A_w) <= 2 E[A_w].
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from flint import arb, ctx, fmpz_poly

from certify_single_random_constituent_dense_outward import (
    exact_caps,
    fraction_of_float,
    upper_float,
)
from evaluate_block_expand_cap_budget import caps_for_variance_factor_and_delta


WORKSTREAM = Path(__file__).resolve().parent
DIAGNOSTIC = WORKSTREAM / "block_expand_accumulate_512_256_d14_t0_5_spectrum.json"
TRANSFER = WORKSTREAM / "block_expand5_F2_conditional_transfer_outward.json"
OUTPUT = WORKSTREAM / "block_expand5_mean_cap_outward.json"
K = 256
B = 512
REGIONS = ((37, 8), (36, 6))
STAGES = 5
VARIANCE_FACTOR = 2
CONDITIONAL_RANDOM_LOG2 = mp.mpf("-51.6439589890")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def occupancy_rows(length: int) -> list[list[int]]:
    rows = [[1] + [0] * length]
    for _balls in range(1, K + 1):
        previous = rows[-1]
        following = [0] * (length + 1)
        for odd, count in enumerate(previous):
            if not count:
                continue
            if odd < length:
                following[odd + 1] += count * (length - odd)
            if odd:
                following[odd - 1] += count * odd
        rows.append(following)
    return rows


def expander_spectrum() -> list[arb]:
    occupancy = {length: occupancy_rows(length) for length, _ in REGIONS}
    spectrum = [arb(0) for _ in range(B + 1)]
    for support in range(K + 1):
        polynomial = fmpz_poly([1])
        denominator = 1
        for length, copies in REGIONS:
            regional = fmpz_poly(occupancy[length][support])
            polynomial *= regional**copies
            denominator *= length ** (copies * support)
        message_count = math.comb(K, support)
        scale = arb(message_count) / denominator
        for weight in range(polynomial.degree() + 1):
            coefficient = int(polynomial[weight])
            if coefficient:
                spectrum[weight] += scale * coefficient
        if support % 16 == 0 or support == K:
            print(f"expander_support,{support},{K}", flush=True)
    return spectrum


def apply_accumulator(spectrum: list[arb]) -> list[arb]:
    result = [arb(0) for _ in range(B + 1)]
    result[0] = spectrum[0]
    for support in range(1, B + 1):
        if spectrum[support] == 0:
            continue
        lower = (support + 1) // 2
        upper = B - support // 2
        denominator = math.comb(B, support)
        down = support // 2
        up = (support + 1) // 2
        for weight in range(lower, upper + 1):
            numerator = math.comb(B - weight, down) * math.comb(
                weight - 1, up - 1
            )
            if numerator:
                result[weight] += spectrum[support] * numerator / denominator
    return result


def log2_bounds(value: arb) -> tuple[float, float]:
    quotient = value.log() / arb(2).log()
    return (
        math.nextafter(float(quotient.lower()), -math.inf),
        math.nextafter(float(quotient.upper()), math.inf),
    )


def main() -> None:
    ctx.prec = 256
    mp.mp.dps = 120
    diagnostic = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    diagnostic_stage = next(
        row for row in diagnostic["stages"] if row["accumulator_stages"] == STAGES
    )
    diagnostic_logs = diagnostic_stage["expected_spectrum_log2_by_weight"]
    frozen_caps, _ = exact_caps()
    diagnostic_kernel = mp.power(
        2, diagnostic_stage["expected_nonzero_kernel_words_log2"]
    )
    diagnostic_tail = sum(
        (
            mp.power(2, diagnostic_logs[weight])
            for weight in range(1, B + 1)
            if not frozen_caps[weight] and diagnostic_logs[weight] is not None
        ),
        mp.mpf(0),
    )
    positive_shells = sum(bool(frozen_caps[weight]) for weight in range(1, B + 1))
    setup_target = mp.power(2, -40) - mp.power(2, CONDITIONAL_RANDOM_LOG2)
    delta = (
        mp.mpf("0.5")
        * (setup_target - diagnostic_kernel - diagnostic_tail)
        / positive_shells
    )
    caps = caps_for_variance_factor_and_delta(
        diagnostic_logs, frozen_caps, VARIANCE_FACTOR, delta
    )

    spectrum = expander_spectrum()
    total_messages = arb(1 << K)
    if not spectrum[0] > 0 or not sum(spectrum, arb(0)).contains(total_messages):
        raise AssertionError("expander spectrum has invalid mass")
    for stage in range(1, STAGES + 1):
        spectrum = apply_accumulator(spectrum)
        if not sum(spectrum, arb(0)).contains(total_messages):
            raise AssertionError(f"stage {stage} spectrum has invalid mass")
        print(f"accumulator_stage,{stage},{STAGES}", flush=True)

    kernel = spectrum[0] - 1
    tail = sum(
        (spectrum[weight] for weight in range(1, B + 1) if not caps[weight]),
        arb(0),
    )
    central = arb(0)
    for weight in range(1, B + 1):
        cap = caps[weight]
        if not cap:
            continue
        mean = spectrum[weight]
        deviation = arb(cap + 1) - mean
        if not deviation > 0:
            raise AssertionError(f"cap does not exceed shell mean at weight {weight}")
        variance_upper = VARIANCE_FACTOR * mean
        central += variance_upper / (variance_upper + deviation * deviation)
    setup_failure = kernel + tail + central
    setup_log2 = log2_bounds(setup_failure)

    transfer = json.loads(TRANSFER.read_text(encoding="utf-8"))
    transfer_log2 = float(transfer["claim"]["conditional_distance_failure_log2_upper"])
    transfer_mass = arb(2) ** (
        arb(fraction_of_float(transfer_log2).numerator)
        / fraction_of_float(transfer_log2).denominator
    )
    combined = setup_failure + transfer_mass
    combined_log2 = log2_bounds(combined)
    if not combined_log2[1] < -40:
        raise AssertionError("outward conditional cap accounting misses 40 bits")

    shell_rows = []
    for weight, mean in enumerate(spectrum):
        lower, upper = log2_bounds(mean) if mean > 0 else (-math.inf, -math.inf)
        diagnostic_value = diagnostic_logs[weight]
        shell_rows.append(
            {
                "weight": weight,
                "mean_log2_lower": lower,
                "mean_log2_upper": upper,
                "binary64_diagnostic_log2": diagnostic_value,
                "cap": str(caps[weight]),
            }
        )

    payload = {
        "schema": "block-expand5-mean-cap-outward-v1",
        "status": "OUTWARD_MEAN_AND_CAP_EVENT_CONDITIONAL_ON_VARIANCE",
        "claim": {
            "mean_spectrum_outward": True,
            "variance_hypothesis": "Var(A_w) <= 2 E[A_w] for every positive shell",
            "outer_event_failure_log2_lower": setup_log2[0],
            "outer_event_failure_log2_upper": setup_log2[1],
            "combined_failure_log2_upper": combined_log2[1],
            "combined_margin_bits_lower": -combined_log2[1],
            "comparison_to_2^-40": True,
        },
        "parameters": {
            "message_bits": K,
            "output_bits": B,
            "regions": [[length, copies] for length, copies in REGIONS],
            "left_degree": sum(copies for _length, copies in REGIONS),
            "accumulator_stages": STAGES,
            "variance_factor": VARIANCE_FACTOR,
            "per_positive_shell_failure_log2": float(mp.log(delta, 2)),
        },
        "components": {
            "kernel_log2_bounds": log2_bounds(kernel),
            "zero_cap_tail_log2_bounds": log2_bounds(tail),
            "central_cantelli_log2_bounds": log2_bounds(central),
            "transfer_log2_upper": transfer_log2,
        },
        "method": (
            "exact integer parity-occupancy polynomials for eight length-37 and six length-36 regions; "
            "five exact binomial accumulator IOWE transitions; 256-bit Arb enclosure"
        ),
        "shells": shell_rows,
        "sources": [
            {"file": DIAGNOSTIC.name, "sha256": sha256(DIAGNOSTIC)},
            {"file": TRANSFER.name, "sha256": sha256(TRANSFER)},
        ],
        "limitations": [
            "The cap event and combined claim remain conditional on the stated length-512 variance inequality.",
            "The RandomStepConv-M22 inner is a mathematical comparator, not the RM2Sub implementation.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(json.dumps(payload["components"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
