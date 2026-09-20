#!/usr/bin/env python3
"""Verify the SPIN transfer for provisional Block Expand-5 shell caps.

The integer-cap transfer uses outward arithmetic.  The cap event itself is
conditional on outward versions of the recorded means and on
Var(A_w) <= 2 E[A_w].  Therefore this receipt is not a distance certificate.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import arb, ctx

import certify_single_random_constituent_sparse_outward as sparse
from certify_single_random_constituent_dense_outward import (
    P_CENTRAL,
    P_DEFECT,
    a,
    exact_band_majorant,
    exact_caps,
    fraction_of_float,
    upper_float,
)
from evaluate_block_expand_cap_budget import (
    caps_for_variance_factor_and_delta,
)


WORKSTREAM = Path(__file__).resolve().parent
SPECTRUM = WORKSTREAM / "block_expand_accumulate_512_256_d14_t0_5_spectrum.json"
RANDOM_SPARSE = WORKSTREAM / "single_random_constituent_B512_sparse_outward_s22.json"
RANDOM_DENSE = WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json"
OUTPUT = WORKSTREAM / "block_expand5_F2_conditional_transfer_outward.json"
B = 512
L = 4096
VARIANCE_FACTOR = 2
CONDITIONAL_RANDOM_LOG2 = mp.mpf("-51.6439589890")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q1_bound_filtered(caps: list[int]) -> tuple[float, dict[str, object]]:
    diagnostic = json.loads(sparse.Q1_DIAGNOSTIC.read_text(encoding="utf-8"))
    groups: dict[float, list[int]] = {}
    for row in diagnostic["weight_rows"]:
        weight = int(row["weight"])
        if caps[weight]:
            groups.setdefault(float(row["log_surprisal"]), []).append(weight)
    pointwise: dict[int, float] = {}
    for group_index, (u, weights) in enumerate(sorted(groups.items())):
        transitions, surprisal = sparse.transition_matrices(u, (Fraction(1),))
        regions, region_exponents = sparse.coefficients_one_dimension(
            transitions[0], 0, transitions[1], 0, L, 1
        )
        coordinates, coordinate_exponents = sparse.coefficients_one_dimension(
            regions[0],
            int(region_exponents[0]),
            regions[1],
            int(region_exponents[1]),
            B,
            B,
        )
        moments = sparse.row_sum_log2_upper(coordinates, coordinate_exponents)
        correction = sparse.correction_log2_upper(surprisal)
        for weight in weights:
            inner = min(0.0, sparse.up(moments[weight] + correction))
            outer = sparse.log2_integer_upper(L * caps[weight])
            pointwise[weight] = sparse.up(outer + inner)
        print(
            f"q1_witness_group,{group_index + 1},{len(groups)}", flush=True
        )
    maximum = max(pointwise.values())
    aggregate = sparse.up(maximum + sparse.log2_integer_upper(len(pointwise)))
    return aggregate, {
        "positive_cap_shells": len(pointwise),
        "maximum_pointwise_log2_upper": maximum,
    }


def exp2_arb(log2_value: float) -> arb:
    return arb(2) ** a(fraction_of_float(log2_value))


def main() -> None:
    ctx.prec = 192
    mp.mp.dps = 120
    spectrum_payload = json.loads(SPECTRUM.read_text(encoding="utf-8"))
    stage = next(
        row
        for row in spectrum_payload["stages"]
        if row["accumulator_stages"] == 5
    )
    log_means = stage["expected_spectrum_log2_by_weight"]
    kernel = mp.power(2, stage["expected_nonzero_kernel_words_log2"])
    frozen_caps, _ = exact_caps()
    tail = mp.mpf(0)
    for weight in range(1, B + 1):
        if not frozen_caps[weight] and log_means[weight] is not None:
            tail += mp.power(2, log_means[weight])
    positive_shells = sum(bool(frozen_caps[weight]) for weight in range(1, B + 1))
    setup_target = mp.power(2, -40) - mp.power(2, CONDITIONAL_RANDOM_LOG2)
    delta = mp.mpf("0.5") * (setup_target - kernel - tail) / positive_shells
    caps = caps_for_variance_factor_and_delta(
        log_means, frozen_caps, VARIANCE_FACTOR, delta
    )

    old_low = exact_band_majorant(frozen_caps, 42, 79, P_DEFECT)
    old_high = exact_band_majorant(frozen_caps, 433, 470, 1 - P_DEFECT)
    old_central = exact_band_majorant(frozen_caps, 80, 432, P_CENTRAL)
    new_low = exact_band_majorant(caps, 42, 79, P_DEFECT)
    new_high = exact_band_majorant(caps, 433, 470, 1 - P_DEFECT)
    new_central = exact_band_majorant(caps, 80, 432, P_CENTRAL)
    domination = {
        "low": new_low <= old_low,
        "high": new_high <= old_high,
        "central": new_central <= old_central,
    }
    if not all(domination.values()):
        raise AssertionError("Block Expand caps do not dominate by the frozen bands")

    q1, q1_details = q1_bound_filtered(caps)
    q2, q2_details = sparse.q2_bound(caps)
    old_sparse = json.loads(RANDOM_SPARSE.read_text(encoding="utf-8"))
    old_dense = json.loads(RANDOM_DENSE.read_text(encoding="utf-8"))
    q3_159 = float(old_sparse["claim"]["q3_159_log2_upper"])
    dense = float(old_dense["claim"]["aggregate_log2_upper"])
    conditional_terms = (q1, q2, q3_159, dense)
    conditional_mass = sum((exp2_arb(value) for value in conditional_terms), arb(0))
    conditional_log2 = upper_float(conditional_mass.log() / arb(2).log())

    central_failure = mp.mpf(0)
    for weight in range(1, B + 1):
        if not caps[weight]:
            continue
        mean = mp.power(2, log_means[weight])
        deviation = mp.mpf(caps[weight] + 1) - mean
        variance = VARIANCE_FACTOR * mean
        central_failure += variance / (variance + deviation * deviation)
    setup_failure = kernel + tail + central_failure
    setup_log2 = mp.log(setup_failure, 2)
    combined = arb(str(setup_failure)) + conditional_mass
    combined_log2 = upper_float(combined.log() / arb(2).log())
    if not combined_log2 < -40:
        raise AssertionError("conditional Block Expand accounting misses 40 bits")

    payload = {
        "schema": "block-expand5-F2-conditional-transfer-outward-v1",
        "status": "OUTWARD_TRANSFER_CONDITIONAL_ON_UNPROVED_CAP_EVENT",
        "claim": {
            "conditional_distance_failure_log2_upper": conditional_log2,
            "diagnostic_outer_event_failure_log2": float(setup_log2),
            "diagnostic_combined_failure_log2_upper": combined_log2,
            "diagnostic_combined_margin_bits": -combined_log2,
            "comparison_to_2^-40": True,
        },
        "hypothesis": (
            "outward shell means equal the recorded values within a certified enclosure and "
            "Var(A_w) <= 2 E[A_w] for every positive shell"
        ),
        "parameters": {
            "message_bits_per_constituent": 256,
            "output_bits_per_constituent": 512,
            "expander_left_degree": 14,
            "accumulator_stages": 5,
            "outer_rows": 4096,
            "inner": "RandomStepConv-M22",
            "positive_cap_shells_before_recentring": positive_shells,
            "positive_cap_shells_after_recentring": sum(bool(value) for value in caps),
            "per_shell_failure_log2": float(mp.log(delta, 2)),
        },
        "transfer": {
            "q1_log2_upper": q1,
            "q1_details": q1_details,
            "q2_log2_upper": q2,
            "q2_details": q2_details,
            "q3_159_log2_upper_reused": q3_159,
            "q160_4096_log2_upper_reused": dense,
            "band_majorant_domination": domination,
            "rule": (
                "Q=1 and Q=2 are reevaluated with the exact integer caps. "
                "For Q>=3, each new band majorant is no larger than the one "
                "used by the cited outward receipts."
            ),
        },
        "setup_accounting_diagnostic": {
            "kernel_log2": float(mp.log(kernel, 2)),
            "zero_cap_tail_log2": float(mp.log(tail, 2)),
            "central_cantelli_union_log2": float(mp.log(central_failure, 2)),
        },
        "sources": [
            {"file": SPECTRUM.name, "sha256": sha256(SPECTRUM)},
            {"file": RANDOM_SPARSE.name, "sha256": sha256(RANDOM_SPARSE)},
            {"file": RANDOM_DENSE.name, "sha256": sha256(RANDOM_DENSE)},
        ],
        "limitations": [
            "The shell means are binary64 diagnostics and are not outward enclosures.",
            "The variance-factor-two hypothesis is open at length 512.",
            "The reported combined margin is conditional accounting, not an unconditional distance theorem.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(json.dumps(payload["transfer"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
