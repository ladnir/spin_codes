#!/usr/bin/env python3
"""Heuristic full-size occupation-three enumerator for Riffle g=4."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import log_ndtr, logsumexp

from analyze_riffle_parity_ladder_zero import DEFAULT_SPECTRUM, load_spectrum


OUTER_SYMBOLS = 16_386
FIELD_NONZERO = (1 << 64) - 1
PACKET_POSITIONS = 524_352
BINARY_LENGTH = 4 * PACKET_POSITIONS
PACKETS_PER_OCCUPATION_THREE_WORD = 96
STATE_WEIGHT_LAW = np.asarray((1, 4, 6, 4, 1), dtype=np.float64) / 16.0
DEFAULT_RELATIVE_DISTANCES = tuple(value / 1000 for value in range(50, 121, 5))
GROWTH_DATA_BLOCKS = (256, 1_024, 4_096, 16_384)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/"
    "target_occ3_bulk_heuristic.json"
)


def one_block_support_distribution() -> np.ndarray:
    spectrum = load_spectrum(DEFAULT_SPECTRUM)
    nonzero_packet = (0, 4, 6, 4, 1)
    powers = [[1]]
    for _ in range(32):
        old = powers[-1]
        following = [0] * (len(old) + 4)
        for left, left_value in enumerate(old):
            for right, right_value in enumerate(nonzero_packet):
                following[left + right] += left_value * right_value
        powers.append(following)

    result = np.zeros(33, dtype=np.float64)
    for support in range(33):
        polynomial = powers[support]
        placements = math.comb(32, support)
        for weight, codewords in spectrum.items():
            if not weight or weight >= len(polynomial):
                continue
            packetizations = placements * polynomial[weight]
            if packetizations:
                result[support] += (
                    codewords
                    * packetizations
                    / math.comb(128, weight)
                    / FIELD_NONZERO
                )
    if abs(float(np.sum(result)) - 1.0) > 1e-12:
        raise RuntimeError("one-block support distribution lost mass")
    return result


def spacing_tail_log(relative_weight: float, active_packets: int) -> float:
    threshold = 4.0 * relative_weight
    coefficients = np.arange(5, dtype=np.float64) - threshold

    def cumulants(tilt: float) -> tuple[float, float, float]:
        denominators = 1.0 - tilt * coefficients
        moment = float(np.sum(STATE_WEIGHT_LAW / denominators))
        first = float(
            np.sum(STATE_WEIGHT_LAW * coefficients / denominators**2)
        )
        second = float(
            np.sum(2.0 * STATE_WEIGHT_LAW * coefficients**2 / denominators**3)
        )
        cumulant = -math.log1p(tilt * threshold) + active_packets * math.log(moment)
        derivative = (
            -threshold / (1.0 + tilt * threshold)
            + active_packets * first / moment
        )
        curvature = (
            threshold**2 / (1.0 + tilt * threshold) ** 2
            + active_packets * (second / moment - (first / moment) ** 2)
        )
        return cumulant, derivative, curvature

    tilt = brentq(
        lambda value: cumulants(value)[1],
        -1.0 / threshold + 1e-12,
        -1e-15,
        xtol=1e-14,
        rtol=1e-14,
    )
    cumulant, _derivative, curvature = cumulants(tilt)
    w = -math.sqrt(-2.0 * cumulant)
    u = tilt * math.sqrt(curvature)
    log_normal = float(log_ndtr(w))
    mills = math.exp(
        -0.5 * w * w - 0.5 * math.log(2.0 * math.pi) - log_normal
    )
    correction = 1.0 + mills * (1.0 / w - 1.0 / u)
    return log_normal + (math.log(correction) if correction > 0.0 else 0.0)


def main() -> None:
    one_block = one_block_support_distribution()
    support_law = np.convolve(np.convolve(one_block, one_block), one_block)
    log_support_law = np.full(support_law.size, -math.inf, dtype=np.float64)
    positive = support_law > 0.0
    log_support_law[positive] = np.log(support_law[positive])
    outer_words = math.comb(OUTER_SYMBOLS, 3) * FIELD_NONZERO
    outer_log = math.log(outer_words)

    def mixture_log(relative_weight: float) -> float:
        terms = [
            log_support_law[support] + spacing_tail_log(relative_weight, support)
            for support in range(1, support_law.size)
            if support_law[support] > 1e-18
        ]
        return float(logsumexp(terms))

    rows = []
    for relative_weight in DEFAULT_RELATIVE_DISTANCES:
        per_word_log = mixture_log(relative_weight)
        rows.append(
            {
                "relative_binary_weight": relative_weight,
                "distance": math.floor(relative_weight * BINARY_LENGTH),
                "heuristic_per_word_log2_probability": per_word_log / math.log(2.0),
                "heuristic_expected_count_log2": (
                    per_word_log + outer_log
                )
                / math.log(2.0),
            }
        )

    crossing = brentq(
        lambda value: mixture_log(value) + outer_log,
        0.05,
        0.15,
        xtol=1e-12,
    )
    growth_rows = []
    for data_blocks in GROWTH_DATA_BLOCKS:
        output_length = 128 * (data_blocks + 2)
        occupation_three_words = (
            math.comb(data_blocks + 2, 3) * FIELD_NONZERO
        )
        occupation_three_log = math.log(occupation_three_words)
        size_crossing = brentq(
            lambda value: mixture_log(value) + occupation_three_log,
            0.03,
            0.15,
            xtol=1e-12,
        )
        growth_rows.append(
            {
                "data_blocks": data_blocks,
                "message_bits": 64 * data_blocks,
                "binary_output_length": output_length,
                "occupation_three_words_log2": (
                    occupation_three_log / math.log(2.0)
                ),
                "heuristic_crossing_relative_weight": size_crossing,
                "heuristic_crossing_distance": math.floor(
                    size_crossing * output_length
                ),
                "heuristic_expected_count_log2_at_9_percent": (
                    mixture_log(0.09) + occupation_three_log
                )
                / math.log(2.0),
            }
        )
    payload = {
        "schema": "riffle-target-occ3-bulk-heuristic-v1",
        "evidence_label": "NONRIGOROUS_BULK_OCCUPATION_THREE_MODEL",
        "packet_positions": PACKET_POSITIONS,
        "binary_length": BINARY_LENGTH,
        "exact_outer_symbol_weight_three_words": str(outer_words),
        "exact_outer_symbol_weight_three_words_log2": outer_log / math.log(2.0),
        "one_block_packet_support_mean": float(
            np.dot(np.arange(one_block.size), one_block)
        ),
        "three_block_packet_support_mean": float(
            np.dot(np.arange(support_law.size), support_law)
        ),
        "three_block_packet_support_standard_deviation": float(
            math.sqrt(
                np.dot(
                    (
                        np.arange(support_law.size)
                        - np.dot(np.arange(support_law.size), support_law)
                    )
                    ** 2,
                    support_law,
                )
            )
        ),
        "heuristic_crossing_relative_weight": crossing,
        "heuristic_crossing_distance": math.floor(crossing * BINARY_LENGTH),
        "growth_rows": growth_rows,
        "rows": rows,
        "assumptions": [
            "Occupation-three coefficient ratios behave like uniform field ratios.",
            "The three nonzero BCH symbols use the product of their exact marginal packet-support laws.",
            "Accumulator state weights use the stationary Binomial(4,1/2) law with negligible short-range correlation.",
            "Uniform finite gaps use the continuous Dirichlet-spacing saddle approximation.",
        ],
        "scope": (
            "The MDS occupation-three word count and the marginal packet-support "
            "law are exact. The joint outer ratio law, accumulator-state law, "
            "and spacing tail are heuristic. This is an engineering estimate, "
            "not a distance certificate."
        ),
    }
    DEFAULT_OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(DEFAULT_OUTPUT),
                "crossing_relative_weight": crossing,
                "crossing_distance": payload["heuristic_crossing_distance"],
                "at_9_percent_log2": next(
                    row["heuristic_expected_count_log2"]
                    for row in rows
                    if row["relative_binary_weight"] == 0.09
                ),
                "growth_rows": growth_rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
