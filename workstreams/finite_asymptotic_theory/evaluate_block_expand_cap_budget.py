#!/usr/bin/env python3
"""Measure the variance budget for Block Expand against the frozen SPIN caps.

This is a high-precision diagnostic.  It imports the exact integer cap vector
from the one-shot random-constituent certificate and asks how large a uniform
bound Var(A_w) <= F E[A_w] may be while retaining a 40-bit combined target.
The Block Expand means are binary64 inputs, so the result is not a proof.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp

from certify_single_random_constituent_dense_outward import (
    P_CENTRAL,
    P_DEFECT,
    exact_band_majorant,
    exact_caps,
)


WORKSTREAM = Path(__file__).resolve().parent
INPUT = WORKSTREAM / "block_expand_accumulate_512_256_d14_t0_5_spectrum.json"
OUTPUT = WORKSTREAM / "block_expand_cap_budget_diagnostic.json"
B = 512
K = 256
CONDITIONAL_DISTANCE_LOG2 = mp.mpf("-51.6439589890")


def exp2(value: float | mp.mpf) -> mp.mpf:
    return mp.power(2, value)


def event_failure(
    log_means: list[float | None], kernel_log2: float, caps: list[int], factor: mp.mpf
) -> tuple[mp.mpf, mp.mpf, mp.mpf, list[int]]:
    kernel = exp2(kernel_log2)
    zero_cap = mp.mpf(0)
    central = mp.mpf(0)
    exceeded = []
    for weight in range(1, B + 1):
        encoded = log_means[weight]
        mean = mp.mpf(0) if encoded is None else exp2(mp.mpf(encoded))
        cap = caps[weight]
        if cap == 0:
            zero_cap += mean
            continue
        deviation = mp.mpf(cap + 1) - mean
        if deviation <= 0:
            exceeded.append(weight)
            central = mp.inf
            continue
        variance = factor * mean
        central += variance / (variance + deviation * deviation)
    return kernel, zero_cap, central, exceeded


def caps_for_variance_factor_and_delta(
    log_means: list[float | None],
    support_caps: list[int],
    factor: int,
    per_shell_failure: mp.mpf,
) -> list[int]:
    """Choose per-shell Cantelli caps with the requested failure target."""
    result = [0] * (B + 1)
    multiplier = (1 - per_shell_failure) / per_shell_failure
    for weight in range(1, B + 1):
        if not support_caps[weight]:
            continue
        mean = exp2(mp.mpf(log_means[weight]))
        deviation = mp.sqrt(mp.mpf(factor) * mean * multiplier)
        result[weight] = max(0, int(mp.ceil(mean + deviation - 1)))
    return result


def caps_for_variance_factor(
    log_means: list[float | None], support_caps: list[int], factor: int
) -> list[int]:
    return caps_for_variance_factor_and_delta(
        log_means, support_caps, factor, mp.mpf(1) / (1 << 51)
    )


def consecutive_ranges(values: list[int]) -> list[list[int]]:
    if not values:
        return []
    result = []
    lower = previous = values[0]
    for value in values[1:]:
        if value != previous + 1:
            result.append([lower, previous])
            lower = value
        previous = value
    result.append([lower, previous])
    return result


def random_mean_log2(weight: int) -> mp.mpf:
    messages = mp.mpf(2) ** K - 1
    return mp.log(messages * math.comb(B, weight) / (mp.mpf(2) ** B), 2)


def main() -> None:
    mp.mp.dps = 100
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    caps, _ = exact_caps()
    bands = (
        ("low", 42, 79, P_DEFECT),
        ("central", 80, 432, P_CENTRAL),
        ("high", 433, 470, 1 - P_DEFECT),
    )
    frozen_band_logs = {
        name: mp.log(exact_band_majorant(caps, lower, upper, probability), 2)
        for name, lower, upper, probability in bands
    }
    rows = []
    for stage in payload["stages"]:
        count = int(stage["accumulator_stages"])
        if count not in (2, 3, 4, 5):
            continue
        log_means = stage["expected_spectrum_log2_by_weight"]
        kernel_log2 = float(stage["expected_nonzero_kernel_words_log2"])
        frozen_terms = event_failure(log_means, kernel_log2, caps, mp.mpf(0))
        ratios = []
        for weight in range(1, B + 1):
            if caps[weight] and log_means[weight] is not None:
                ratios.append((mp.mpf(log_means[weight]) - random_mean_log2(weight), weight))
        maximum_ratio, ratio_weight = max(ratios)
        variance_rows = []
        for factor in (1, 2, 4, 16, 512):
            candidate_caps = caps_for_variance_factor(log_means, caps, factor)
            terms = event_failure(
                log_means, kernel_log2, candidate_caps, mp.mpf(factor)
            )
            cap_inflations = [
                (mp.mpf(candidate_caps[weight]) / caps[weight], weight)
                for weight in range(1, B + 1)
                if caps[weight]
            ]
            maximum_cap_inflation, cap_inflation_weight = max(cap_inflations)
            total = terms[0] + terms[1] + terms[2]
            band_inflations = {}
            for name, lower, upper, probability in bands:
                candidate = exact_band_majorant(
                    candidate_caps, lower, upper, probability
                )
                band_inflations[name] = float(
                    mp.log(candidate, 2) - frozen_band_logs[name]
                )
            variance_rows.append(
                {
                    "assumed_variance_factor": factor,
                    "outer_event_failure_log2": float(mp.log(total, 2)),
                    "central_cantelli_union_log2": float(mp.log(terms[2], 2)),
                    "maximum_cap_inflation_over_frozen_log2": float(
                        mp.log(maximum_cap_inflation, 2)
                    ),
                    "maximum_cap_inflation_weight": cap_inflation_weight,
                    "two_band_majorant_inflation_log2": band_inflations,
                }
            )
        positive_shells = sum(bool(caps[weight]) for weight in range(1, B + 1))
        setup_target = exp2(-40) - exp2(CONDITIONAL_DISTANCE_LOG2)
        fixed_failure = frozen_terms[0] + frozen_terms[1]
        available = setup_target - fixed_failure
        budget_row = None
        if available > 0:
            # Retain half of the available setup mass as diagnostic
            # reserve for outward rounding and bookkeeping.
            delta = mp.mpf("0.5") * available / positive_shells
            budget_caps = caps_for_variance_factor_and_delta(
                log_means, caps, 2, delta
            )
            budget_terms = event_failure(
                log_means, kernel_log2, budget_caps, mp.mpf(2)
            )
            budget_inflations = {}
            for name, lower, upper, probability in bands:
                candidate = exact_band_majorant(
                    budget_caps, lower, upper, probability
                )
                budget_inflations[name] = float(
                    mp.log(candidate, 2) - frozen_band_logs[name]
                )
            budget_row = {
                "assumed_variance_factor": 2,
                "per_positive_shell_failure_log2": float(mp.log(delta, 2)),
                "positive_shells": positive_shells,
                "outer_event_failure_log2": float(
                    mp.log(sum(budget_terms[:3]), 2)
                ),
                "two_band_majorant_inflation_log2": budget_inflations,
            }
        rows.append(
            {
                "accumulator_stages": count,
                "frozen_cap_mean_excess_ranges": consecutive_ranges(
                    frozen_terms[3]
                ),
                "kernel_failure_markov_log2": float(mp.log(frozen_terms[0], 2)),
                "zero_cap_tail_markov_log2": float(mp.log(frozen_terms[1], 2)),
                "maximum_mean_ratio_to_random_log2_on_positive_caps": float(maximum_ratio),
                "maximum_mean_ratio_weight": ratio_weight,
                "variance_scenarios": variance_rows,
                "factor_two_spend_to_combined_40_bit_budget": budget_row,
            }
        )

    result = {
        "schema": "block-expand-cap-budget-diagnostic-v1",
        "status": "BINARY64_MEANS_HIGH_PRECISION_ACCOUNTING_DIAGNOSTIC",
        "construction": "degree-14 regional expander followed by t independently interleaved accumulators",
        "cap_interface": "frozen integer caps from the one-shot random-[512,256] RandomStepConv-M22 certificate",
        "target": "outer-event failure plus the certified conditional distance failure is below 2^-40",
        "rows": rows,
        "limitations": [
            "Block Expand shell means are binary64 diagnostics, not outward enclosures.",
            "The variance scenarios are conditional targets; no variance inequality is proved here.",
            "Each scenario recenters its positive caps by Cantelli at per-shell failure 2^-51 and preserves the frozen zero/nonzero support.",
            "Only the Block Expand-5 factor-two budget-spending caps have been propagated through the complete conditional SPIN transfer.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
