#!/usr/bin/env python3
"""Audit late starts and terminated episodes in one-lap RandomStepConv."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import gammaln, logsumexp
from scipy.stats import binom

from analyze_riffle_rm512_randomstepconv_g4_sigma20 import (
    DATA_GROUPS,
    DEFAULT_SPECTRUM,
    PACKETS_PER_LOCAL,
    ExactRMMoments,
)
from analyze_riffle_randomstepconv_g4_sigma20_goal01 import (
    TARGET_G,
    TARGET_PACKETS,
    TARGET_SIGMA,
    episode_terms,
    log_binomial,
    optimize_tilts,
)
from analyze_riffle_randomstepconv_g4_sigma20_goal02 import optimize_point


LOG2 = math.log(2.0)
DISTANCE = 188766
BULK_SUPPORT = 32767
LOW_SUPPORT = 95
MAX_PARITY_PACKETS = 64
DEFAULT_RM_RECEIPT = Path(
    "constructions/riffle_rm512_randomstepconv_g4_sigma20/receipts/"
    "rm512_randomstepconv_distance09.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_rm512_randomstepconv_g4_sigma20/receipts/"
    "inner_episode_structure_audit.json"
)


def small_trajectory_count_check(maximum_n: int = 7) -> dict[str, object]:
    """Check the surviving-final-episode count by exhaustive state paths."""
    comparisons = 0
    for n in range(1, maximum_n + 1):
        observed: Counter[tuple[int, int, int]] = Counter()
        for support_mask in range(1, 1 << n):
            support = support_mask.bit_count()
            stack = [(0, False, 0, 0)]
            while stack:
                position, state_live, terminations, live_steps = stack.pop()
                if position == n:
                    if state_live:
                        observed[(support, terminations, live_steps)] += 1
                    continue
                active = bool((support_mask >> position) & 1)
                live = state_live or active
                if not live:
                    stack.append(
                        (position + 1, False, terminations, live_steps)
                    )
                    continue
                stack.append(
                    (position + 1, True, terminations, live_steps + 1)
                )
                stack.append(
                    (position + 1, False, terminations + 1, live_steps + 1)
                )
        for support in range(1, n + 1):
            for terminations in range(support):
                for live_steps in range(support, n + 1):
                    if live_steps - terminations - 1 < support - terminations - 1:
                        expected = 0
                    else:
                        expected = (
                            math.comb(live_steps - 1, terminations)
                            * math.comb(
                                n - live_steps + terminations, terminations
                            )
                            * math.comb(
                                live_steps - terminations - 1,
                                support - terminations - 1,
                            )
                        )
                    comparisons += 1
                    if observed[(support, terminations, live_steps)] != expected:
                        raise AssertionError(
                            "surviving-final-episode count identity failed"
                        )
    return {
        "maximum_packet_positions": maximum_n,
        "comparisons": comparisons,
        "status": "EXACT_INTEGER_MATCH",
    }


def rm_outer_coefficient_log2(
    model: ExactRMMoments, support: int
) -> float:
    def centered(log_t: float) -> float:
        _log_normalizer, probabilities = model.tilted_block_distribution(log_t)
        return (
            DATA_GROUPS * float(probabilities @ model.supports) - support
        )

    log_t = brentq(centered, -50.0, 5.0)
    log_normalizer, block = model.tilted_block_distribution(log_t)
    maximum_degree = DATA_GROUPS * PACKETS_PER_LOCAL
    transform_size = 1 << maximum_degree.bit_length()
    distribution = np.fft.irfft(
        np.fft.rfft(np.pad(block, (0, transform_size - len(block))))
        ** DATA_GROUPS,
        transform_size,
    )[: maximum_degree + 1]
    probability = float(distribution[support])
    if probability <= 0.0:
        raise ArithmeticError("FFT failed at requested RM support")
    return (
        DATA_GROUPS * log_normalizer
        + math.log(probability)
        - support * log_t
    ) / LOG2


def occupation_decomposition(
    model: ExactRMMoments, support: int
) -> list[dict[str, object]]:
    active_logs = logsumexp(
        model.log_multiplicities[model.nonzero, None]
        + model.log_support_probabilities[model.nonzero],
        axis=0,
    )

    def active_mean(log_t: float) -> float:
        tilted = active_logs + model.supports * log_t
        normalizer = float(logsumexp(tilted))
        return float(np.exp(tilted - normalizer) @ model.supports)

    log_t = brentq(lambda value: active_mean(value) - support, -50.0, 5.0)
    tilted = active_logs + model.supports * log_t
    log_active_normalizer = float(logsumexp(tilted))
    active_support = np.exp(tilted - log_active_normalizer)
    minimum_support = int(np.flatnonzero(np.isfinite(active_logs))[0])
    maximum_occupation = support // minimum_support
    convolution = np.asarray([1.0])
    rows: list[dict[str, object]] = []
    for occupation in range(1, maximum_occupation + 1):
        convolution = np.convolve(convolution, active_support)[: support + 1]
        probability = (
            float(convolution[support])
            if support < len(convolution)
            else 0.0
        )
        if probability <= 0.0:
            continue
        log_coefficient = (
            gammaln(DATA_GROUPS + 1)
            - gammaln(occupation + 1)
            - gammaln(DATA_GROUPS - occupation + 1)
            + occupation * log_active_normalizer
            - support * log_t
            + math.log(probability)
        ) / LOG2
        rows.append(
            {
                "active_group_count": occupation,
                "log2_coefficient": float(log_coefficient),
            }
        )
    total = float(
        logsumexp([row["log2_coefficient"] * LOG2 for row in rows])
    ) / LOG2
    for row in rows:
        row["log2_fraction_of_support_coefficient"] = (
            float(row["log2_coefficient"]) - total
        )
    return rows


def termination_upper_profile(support: int) -> dict[str, object]:
    optimized = optimize_tilts(
        TARGET_PACKETS, TARGET_G, TARGET_SIGMA, support, DISTANCE
    )
    terms, radii = episode_terms(
        TARGET_PACKETS,
        TARGET_G,
        TARGET_SIGMA,
        support,
        float(optimized["weight_tilt"]),
    )
    normalizer = float(logsumexp(terms))
    probabilities = np.exp(terms - normalizer)
    counts = np.arange(support + 1, dtype=np.float64)
    cumulative = np.cumsum(probabilities)
    common_point = optimize_point(support, DISTANCE)
    base = (
        -log_binomial(TARGET_PACKETS, support)
        - DISTANCE * math.log(float(optimized["weight_tilt"]))
    )
    zero_termination_upper = float(terms[0] + base) / LOG2
    return {
        "termination_separated_log2_upper": optimized[
            "diagnostic_log2_upper"
        ],
        "common_radius_pointwise_log2_upper": common_point[
            "diagnostic_log2_upper"
        ],
        "dominant_termination_count": int(np.argmax(terms)),
        "mean_termination_count_under_upper_tilt": float(
            probabilities @ counts
        ),
        "termination_count_quantiles": {
            "q10": int(np.searchsorted(cumulative, 0.10)),
            "q50": int(np.searchsorted(cumulative, 0.50)),
            "q90": int(np.searchsorted(cumulative, 0.90)),
        },
        "zero_termination_component_log2_upper": zero_termination_upper,
        "dominant_gap_tilt": float(
            radii[int(np.argmax(terms))]
        ),
        "weight_tilt": optimized["weight_tilt"],
    }


def surviving_final_episode_rows(
    support: int,
    maximum_terminations: int,
    precomputed_tail: tuple[np.ndarray, np.ndarray] | None = None,
) -> list[dict[str, object]]:
    if precomputed_tail is None:
        live = np.arange(support, TARGET_PACKETS + 1, dtype=np.int64)
        output_tail = binom.logcdf(DISTANCE, TARGET_G * live, 0.5)
    else:
        all_live, all_tail = precomputed_tail
        offset = support - int(all_live[0])
        live = all_live[offset:]
        output_tail = all_tail[offset:]
    live_float = live.astype(np.float64)
    denominator = (
        gammaln(TARGET_PACKETS + 1)
        - gammaln(support + 1)
        - gammaln(TARGET_PACKETS - support + 1)
    )
    q = math.ldexp(1.0, -TARGET_SIGMA)
    rows: list[dict[str, object]] = []
    for terminations in range(maximum_terminations + 1):
        episodes = terminations + 1
        if episodes > support:
            break
        log_episode_lengths = (
            gammaln(live_float)
            - gammaln(terminations + 1)
            - gammaln(live_float - terminations)
        )
        off_positions = TARGET_PACKETS - live + terminations
        log_off_gaps = (
            gammaln(off_positions + 1)
            - gammaln(terminations + 1)
            - gammaln(off_positions - terminations + 1)
        )
        non_start_live = live - episodes
        non_start_active = support - episodes
        log_active_placements = (
            gammaln(non_start_live + 1)
            - gammaln(non_start_active + 1)
            - gammaln(non_start_live - non_start_active + 1)
        )
        logs = (
            log_episode_lengths
            + log_off_gaps
            + log_active_placements
            - denominator
            + terminations * math.log(q)
            + (live - terminations) * math.log1p(-q)
            + output_tail
        )
        log_total = float(logsumexp(logs))
        probabilities = np.exp(logs - log_total)
        cumulative = np.cumsum(probabilities)
        rows.append(
            {
                "termination_count": terminations,
                "log2_probability": log_total / LOG2,
                "live_steps_mode": int(live[int(np.argmax(logs))]),
                "live_steps_median": int(
                    live[int(np.searchsorted(cumulative, 0.50))]
                ),
            }
        )
    return rows


def family_summary(
    support: int,
    outer_log2: float,
    maximum_terminations: int,
    precomputed_tail: tuple[np.ndarray, np.ndarray] | None = None,
) -> dict[str, object]:
    rows = surviving_final_episode_rows(
        support, maximum_terminations, precomputed_tail
    )
    family_log2 = float(
        logsumexp(
            [float(row["log2_probability"]) * LOG2 for row in rows]
        )
    ) / LOG2
    return {
        "full_packet_support": support,
        "termination_counts_included": [0, maximum_terminations],
        "log2_probability": family_log2,
        "outer_weighted_expected_count_log2": outer_log2 + family_log2,
        "dominant_rows": sorted(
            rows,
            key=lambda row: float(row["log2_probability"]),
            reverse=True,
        )[:8],
        "event": (
            "e terminated live episodes followed by one live episode that "
            "survives through the final position; total live length is summed; "
            "output weight is at most D"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--rm-receipt", type=Path, default=DEFAULT_RM_RECEIPT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    rm_receipt = json.loads(args.rm_receipt.read_text(encoding="utf-8"))
    model = ExactRMMoments(args.spectrum)
    low_outer = rm_outer_coefficient_log2(model, LOW_SUPPORT)
    bulk_outer = float(
        rm_receipt["result"]["representative_support_audit"][
            "relaxed_data_outer_coefficient_log2"
        ]
    )
    bulk_interval_inner = float(
        rm_receipt["result"]["representative_support_audit"][
            "inner_probability_upper_log2"
        ]
    )
    low_upper = termination_upper_profile(LOW_SUPPORT)
    bulk_upper = termination_upper_profile(BULK_SUPPORT)

    live = np.arange(LOW_SUPPORT, TARGET_PACKETS + 1, dtype=np.int64)
    precomputed_tail = (live, binom.logcdf(DISTANCE, TARGET_G * live, 0.5))
    low_no_parity = family_summary(
        LOW_SUPPORT, low_outer, 15, precomputed_tail
    )
    low_max_parity = family_summary(
        LOW_SUPPORT + MAX_PARITY_PACKETS,
        low_outer,
        15,
        precomputed_tail,
    )
    parity_crossing = []
    for parity_packets in range(13):
        summary = family_summary(
            LOW_SUPPORT + parity_packets,
            low_outer,
            15,
            precomputed_tail,
        )
        parity_crossing.append(
            {
                "parity_packets": parity_packets,
                "outer_weighted_expected_count_log2": summary[
                    "outer_weighted_expected_count_log2"
                ],
            }
        )

    payload = {
        "schema": "riffle-randomstepconv-episode-structure-audit-v1",
        "model": "Riffle RM512-RandomStepConv g=4 sigma=20",
        "parameters": {
            "packet_positions": TARGET_PACKETS,
            "packet_bits": TARGET_G,
            "state_bits": TARGET_SIGMA,
            "distance_threshold": DISTANCE,
            "bulk_data_support": BULK_SUPPORT,
            "low_data_support": LOW_SUPPORT,
        },
        "exact_checks": {
            "surviving_final_episode_count": small_trajectory_count_check()
        },
        "bulk_support_audit": {
            "data_support": BULK_SUPPORT,
            "outer_coefficient_log2": bulk_outer,
            "wide_interval_inner_log2_upper": bulk_interval_inner,
            **bulk_upper,
            "outer_plus_pointwise_inner_log2_upper": (
                bulk_outer
                + float(bulk_upper["termination_separated_log2_upper"])
            ),
            "wide_interval_parameter_loss_bits": (
                bulk_interval_inner
                - float(bulk_upper["common_radius_pointwise_log2_upper"])
            ),
            "interpretation": (
                "The wide interval reused the support-16384 parameters at the "
                "support-32767 endpoint. Pointwise optimization removes almost "
                "all of the apparent bulk obstruction."
            ),
        },
        "low_support_audit": {
            "data_support": LOW_SUPPORT,
            "outer_coefficient_log2": low_outer,
            "outer_occupation_decomposition": occupation_decomposition(
                model, LOW_SUPPORT
            ),
            **low_upper,
            "outer_plus_inner_log2_upper": (
                low_outer
                + float(low_upper["termination_separated_log2_upper"])
            ),
            "zero_parity_explicit_family": low_no_parity,
            "max_parity_explicit_family": low_max_parity,
            "parity_support_crossing": parity_crossing,
            "first_negative_explicit_family_parity_support": next(
                row["parity_packets"]
                for row in parity_crossing
                if row["outer_weighted_expected_count_log2"] < 0.0
            ),
        },
        "conclusion": (
            "Late start remains the live-length mechanism. The bulk positive "
            "bound was an interval-parameter artifact. The residual explicit "
            "failure family uses about six live episodes at data support 95 "
            "and survives only when the parity block has at most three active "
            "packets. The proof/refutation question therefore moves to the "
            "joint spectrum of one RM data group and its two field parities."
        ),
        "scope": (
            "The pointwise upper bounds are floating diagnostics from the "
            "previously verified fixed-support identity. The explicit episode "
            "family is a probability lower bound; its trajectory count has an "
            "independent exhaustive small-instance check. This audit does not "
            "bound the parity-kernel RM subcode."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"bulk_interval_loss_bits,{payload['bulk_support_audit']['wide_interval_parameter_loss_bits']:.6f}")
    print(f"bulk_pointwise_combined,{payload['bulk_support_audit']['outer_plus_pointwise_inner_log2_upper']:.6f}")
    print(f"low_support,{LOW_SUPPORT}")
    print(f"low_pointwise_combined,{payload['low_support_audit']['outer_plus_inner_log2_upper']:.6f}")
    print(f"zero_parity_explicit,{low_no_parity['outer_weighted_expected_count_log2']:.6f}")
    print(f"max_parity_explicit,{low_max_parity['outer_weighted_expected_count_log2']:.6f}")
    print(
        "first_negative_parity_support,"
        f"{payload['low_support_audit']['first_negative_explicit_family_parity_support']}"
    )
    print(f"output,{args.output}")


if __name__ == "__main__":
    main()
