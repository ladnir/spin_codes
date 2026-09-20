#!/usr/bin/env python3
"""Audit multiplicity accounting in the fixed-rate Outer256 model.

The audit separates the modeled constituent spectrum, the exact data-block
generating function under that model, and the remaining worst-case
relaxations. It is a floating diagnostic, not a code-spectrum certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom

from analyze_riffle_outer256model_randomstepconv_g4_sigma20 import (
    DATA_PAIRS,
    LOG2,
    PACKETS_PER_BLOCK,
    TARGET_PACKETS,
    ModeledLocalMoments,
    canonical_late_suffix_event,
    log_binomial,
)


DEFAULT_SOURCE = Path(
    "constructions/riffle_outer256model_randomstepconv_g4_sigma20/receipts/"
    "fixed_rate_outer256_model.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_outer256model_randomstepconv_g4_sigma20/receipts/"
    "outer_multiplicity_audit.json"
)


def quantiles(values: np.ndarray, probabilities: np.ndarray) -> dict[str, int]:
    levels = (0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)
    cumulative = np.cumsum(probabilities)
    return {
        f"q{round(100 * level):02d}": int(
            values[min(int(np.searchsorted(cumulative, level)), len(values) - 1)]
        )
        for level in levels
    }


def local_profile(model: ModeledLocalMoments, log_t: float) -> dict[str, object]:
    joint_logs = (
        model.log_multiplicities[:, None]
        + model.log_support_probabilities
        + model.supports[None, :] * log_t
    )[model.nonzero]
    joint = np.exp(joint_logs - float(logsumexp(joint_logs)))
    weights = model.weights[model.nonzero]
    weight_probabilities = joint.sum(axis=1)
    support_probabilities = joint.sum(axis=0)
    return {
        "mean_binary_weight": float(weight_probabilities @ weights),
        "binary_weight_quantiles": quantiles(weights, weight_probabilities),
        "mean_packet_support": float(support_probabilities @ model.supports),
        "packet_support_quantiles": quantiles(
            model.supports.astype(np.int64), support_probabilities
        ),
        "weight_mass": {
            "at_most_70": float(weight_probabilities[weights <= 70].sum()),
            "from_90_through_112": float(
                weight_probabilities[(weights >= 90) & (weights <= 112)].sum()
            ),
            "at_least_114": float(weight_probabilities[weights >= 114].sum()),
        },
    }


def tilted_block_distribution(
    model: ModeledLocalMoments, log_t: float
) -> tuple[float, np.ndarray]:
    active_by_support = logsumexp(
        model.log_multiplicities[model.nonzero, None]
        + model.log_support_probabilities[model.nonzero],
        axis=0,
    )
    logs = active_by_support + model.supports * log_t
    logs[0] = np.logaddexp(logs[0], 0.0)
    log_normalizer = float(logsumexp(logs))
    return log_normalizer, np.exp(logs - log_normalizer)


def convolved_support_distribution(block: np.ndarray) -> np.ndarray:
    maximum_degree = DATA_PAIRS * PACKETS_PER_BLOCK
    transform_size = 1 << maximum_degree.bit_length()
    padded = np.pad(block, (0, transform_size - len(block)))
    distribution = np.fft.irfft(
        np.fft.rfft(padded) ** DATA_PAIRS, transform_size
    )[: maximum_degree + 1]
    minimum_before_clipping = float(distribution.min())
    distribution = np.maximum(distribution, 0.0)
    distribution /= distribution.sum()
    if minimum_before_clipping < -1e-12:
        raise ArithmeticError("FFT support distribution has material negative mass")
    return distribution


def select_baseline(payload: dict[str, object]) -> dict[str, object]:
    matches = [
        row
        for row in payload["results"]
        if int(row["modeled_minimum_distance"]) == 38
        and float(row["low_weight_inflation_bits"]) == 0.0
    ]
    if len(matches) != 1:
        raise ValueError("expected one distance-38, zero-inflation baseline")
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    source = json.loads(args.source.read_text(encoding="utf-8"))
    baseline = select_baseline(source)
    bad_event = baseline["canonical_bad_event"]
    representative_support = int(
        bad_event["outer_saddle_profile"]["representative_total_packet_support"]
    )
    dominant = max(
        baseline["segments"], key=lambda row: float(row["log2_contribution_upper"])
    )
    combined_log_t = (
        float(dominant["log_support_base"])
        + float(dominant["reciprocal_binomial_secant_slope_log2"]) * LOG2
        + float(dominant["support_log_tilt"])
    )
    support_tail_tilt = math.exp(float(dominant["support_log_tilt"]))

    model = ModeledLocalMoments(38, 0.0)
    log_block_normalizer, block_distribution = tilted_block_distribution(
        model, combined_log_t
    )
    total_distribution = convolved_support_distribution(block_distribution)
    saddle_probability = float(total_distribution[representative_support])
    relaxed_coefficient_log2 = (
        DATA_PAIRS * log_block_normalizer
        + math.log(saddle_probability)
        - representative_support * combined_log_t
    ) / LOG2

    log_m1, _log_m2 = model.log_moments(combined_log_t)
    active_probability = 1.0 / (1.0 + math.exp(-log_m1))
    active_mean = DATA_PAIRS * active_probability
    active_stddev = math.sqrt(
        DATA_PAIRS * active_probability * (1.0 - active_probability)
    )

    left = int(dominant["support_left"])
    right = int(dominant["support_right"])
    support_indices = np.arange(left, right + 1)
    chernoff_ratio = float(
        np.sum(
            total_distribution[left : right + 1]
            * np.exp((right - support_indices) * math.log(support_tail_tilt))
        )
    )

    secant_intercept = (
        float(dominant["reciprocal_binomial_secant_intercept_log2"]) * LOG2
    )
    secant_slope = (
        float(dominant["reciprocal_binomial_secant_slope_log2"]) * LOG2
    )
    exact_reciprocal_log = -log_binomial(TARGET_PACKETS, representative_support)
    secant_gap = (
        secant_intercept
        + secant_slope * representative_support
        - exact_reciprocal_log
    ) / LOG2

    inner_upper_log2 = (
        float(dominant["log_constant"])
        + representative_support * float(dominant["log_support_base"])
        - log_binomial(TARGET_PACKETS, representative_support)
    ) / LOG2
    pointwise_upper_log2 = relaxed_coefficient_log2 + inner_upper_log2
    segment_upper_log2 = float(dominant["log2_contribution_upper"])

    explicit_event = canonical_late_suffix_event(
        representative_support, int(baseline["distance"])
    )
    modeled_explicit_event_log2 = (
        relaxed_coefficient_log2 + float(explicit_event["log2_probability_lower"])
    )

    log_q = model.log_q(combined_log_t)
    ordinary_nonzero_parity_log_moment = (
        float(
            logsumexp(
                model.log_multiplicities[model.nonzero] + log_q[model.nonzero]
            )
        )
        - math.log((1 << 128) - 1)
    ) / LOG2

    payload = {
        "schema": "riffle-outer256model-multiplicity-audit-v1",
        "source": str(args.source),
        "question": (
            "Does the dominant outer calculation sum configurations with their "
            "spectrum multiplicities, or replace every active block by its "
            "minimum distance?"
        ),
        "answer": (
            "It sums the full modeled local weight and packet-support spectrum. "
            "Minimum distance only truncates the guessed spectrum. The dominant "
            "bad-event tilt changes the typical active-block weight from 128 to "
            "101.1; it does not change it to the minimum weight 38."
        ),
        "configuration_sum": {
            "identity": (
                "F(x)^8192, where F(x)=1+sum_{w>0,s} A_w P[S=s|w] x^s"
            ),
            "interpretation": (
                "Expanding F(x)^8192 sums every data-block weight/support "
                "configuration with its multinomial placement count and product "
                "of local multiplicities."
            ),
            "status_relative_to_modeled_spectrum": "EXACT",
        },
        "unconditioned_active_block": local_profile(model, 0.0),
        "dominant_conditioned_active_block": local_profile(model, combined_log_t),
        "dominant_occupation": {
            "active_probability_under_tilt": active_probability,
            "mean_active_blocks": active_mean,
            "standard_deviation_active_blocks": active_stddev,
            "active_block_quantiles": {
                "q01": int(binom.ppf(0.01, DATA_PAIRS, active_probability)),
                "q10": int(binom.ppf(0.10, DATA_PAIRS, active_probability)),
                "q50": int(binom.ppf(0.50, DATA_PAIRS, active_probability)),
                "q90": int(binom.ppf(0.90, DATA_PAIRS, active_probability)),
                "q99": int(binom.ppf(0.99, DATA_PAIRS, active_probability)),
            },
        },
        "representative_support_audit": {
            "packet_support": representative_support,
            "combined_packet_support_tilt": math.exp(combined_log_t),
            "saddle_probability_of_exact_support": saddle_probability,
            "relaxed_data_outer_coefficient_log2": relaxed_coefficient_log2,
            "inner_probability_upper_log2": inner_upper_log2,
            "pointwise_first_moment_upper_log2": pointwise_upper_log2,
            "dominant_segment_upper_log2": segment_upper_log2,
            "segment_minus_pointwise_upper_bits": (
                segment_upper_log2 - pointwise_upper_log2
            ),
            "explicit_late_suffix_probability_lower_log2": explicit_event[
                "log2_probability_lower"
            ],
            "modeled_explicit_event_expected_count_log2": (
                modeled_explicit_event_log2
            ),
        },
        "outer_relaxations": {
            "reciprocal_binomial_secant_gap_at_representative_support_bits": (
                secant_gap
            ),
            "support_interval_chernoff_loss_bits": -math.log2(chernoff_ratio),
            "many_pair_parity_treatment": (
                "The parity-pair block is discarded when the support tilt is "
                "below one. This is the only block-level worst-case treatment "
                "in the many-active-block branch."
            ),
            "ordinary_nonzero_block_parity_moment_log2_heuristic": (
                ordinary_nonzero_parity_log_moment
            ),
            "parity_heuristic_scope": (
                "The heuristic treats the determined parity-pair value like a "
                "uniform nonzero modeled input. Correlations prevent using this "
                "number as a proof."
            ),
        },
        "verdict": (
            "The outer bulk multiplicity is fully exploited relative to the "
            "guessed spectrum. The segment aggregation is within 27.3 bits of "
            "the exact relaxed point calculation at h=32748. The positive "
            "first-moment upper bound is not caused by assigning minimum weight "
            "to each block. The remaining large uncertainty is the inner upper "
            "bound at high packet support, plus the unproved local spectrum and "
            "the worst-case parity-pair treatment."
        ),
        "scope": (
            "The local multiplicities are real-valued random-like guesses. The "
            "FFT coefficient is exact up to floating error for that model after "
            "discarding the determined parity-pair support. It is not evidence "
            "that an explicit [256,128] code has this spectrum."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"output,{args.output}")
    print(f"raw_mean_weight,{payload['unconditioned_active_block']['mean_binary_weight']:.6f}")
    print(
        "conditioned_mean_weight,"
        f"{payload['dominant_conditioned_active_block']['mean_binary_weight']:.6f}"
    )
    print(f"relaxed_outer_coefficient_log2,{relaxed_coefficient_log2:.6f}")
    print(f"pointwise_upper_log2,{pointwise_upper_log2:.6f}")
    print(f"segment_minus_pointwise_bits,{segment_upper_log2 - pointwise_upper_log2:.6f}")
    print(f"modeled_explicit_event_log2,{modeled_explicit_event_log2:.6f}")


if __name__ == "__main__":
    main()
