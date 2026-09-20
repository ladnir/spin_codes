#!/usr/bin/env python3
"""Measure BCH-weight correlations induced by alpha_i = gamma^i.

The bulk experiment samples nonzero x in GF(2^64) and coefficient indices i,
then records

    (wt(B(x)), wt(B(gamma^i x))).

The exact-tail experiment reads authenticated low-packet-support BCH messages.
It advances each selected source through all 16,384 coefficient lags and uses
a hash lookup to count targets that also occur in the authenticated list.

Bulk results are statistical evidence.  Tail intersections are exact only
inside the supplied exhaustive packet-support range.  The script prints JSON
and does not modify the workspace.
"""

from __future__ import annotations

import argparse
import json
import math
import struct
import sys
from collections import Counter
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as dp  # noqa: E402


COEFFICIENT_COUNT = dp.USER_BLOCKS
CODE_BITS = 128
PACKET_BITS = 4
PACKET_COUNT = CODE_BITS // PACKET_BITS
DEFAULT_FIXED_LAGS = (0, 1, 2, 3, 7, 15, 31, 63, 127, 255, 1023, 4095, 8191, 16383)
DEFAULT_THRESHOLDS = (32, 36, 40, 44, 48, 52)


def parse_int_list(value: str) -> tuple[int, ...]:
    result = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    if not result:
        raise argparse.ArgumentTypeError("list must not be empty")
    return result


def build_encode_tables() -> tuple[np.ndarray, np.ndarray]:
    """Build eight byte-position tables for the fixed linear BCH encoder."""
    low = np.zeros((8, 256), dtype=np.uint64)
    high = np.zeros((8, 256), dtype=np.uint64)
    for byte_position in range(8):
        for byte_value in range(256):
            codeword = 0
            for bit in range(8):
                if (byte_value >> bit) & 1:
                    codeword ^= dp.GENERATOR_ROWS[8 * byte_position + bit]
            low[byte_position, byte_value] = np.uint64(codeword & dp.FIELD_ORDER)
            high[byte_position, byte_value] = np.uint64(codeword >> 64)
    return low, high


ENCODE_LOW, ENCODE_HIGH = build_encode_tables()


def encode_weights(messages: np.ndarray) -> np.ndarray:
    """Return BCH binary weights for a uint64 message batch."""
    low = np.zeros(messages.size, dtype=np.uint64)
    high = np.zeros(messages.size, dtype=np.uint64)
    for byte_position in range(8):
        values = ((messages >> np.uint64(8 * byte_position)) & np.uint64(0xFF)).astype(
            np.uint8
        )
        low ^= ENCODE_LOW[byte_position, values]
        high ^= ENCODE_HIGH[byte_position, values]
    return np.bitwise_count(low).astype(np.uint16) + np.bitwise_count(high).astype(
        np.uint16
    )


def field_multiply_batch(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Carryless GF(2^64) multiplication on equal-size uint64 batches."""
    result = np.zeros(left.size, dtype=np.uint64)
    left = left.copy()
    right = right.copy()
    reduction = np.uint64(dp.FIELD_REDUCTION)
    one = np.uint64(1)
    for _ in range(64):
        result ^= left * (right & one)
        top = left >> np.uint64(63)
        left = (left << one) ^ (top * reduction)
        right >>= one
    return result


def coefficient_powers() -> np.ndarray:
    powers = np.empty(COEFFICIENT_COUNT, dtype=np.uint64)
    powers[0] = np.uint64(1)
    for index in range(1, COEFFICIENT_COUNT):
        powers[index] = np.uint64(dp.field_multiply_x(int(powers[index - 1])))
    if np.unique(powers).size != COEFFICIENT_COUNT:
        raise RuntimeError("coefficient schedule contains a repeated field element")
    return powers


def random_nonzero_messages(rng: np.random.Generator, count: int) -> np.ndarray:
    return rng.integers(1, 1 << 64, size=count, dtype=np.uint64)


def summarize_pairs(
    source_weights: np.ndarray,
    target_weights: np.ndarray,
    thresholds: tuple[int, ...],
) -> dict[str, object]:
    source_float = source_weights.astype(np.float64)
    target_float = target_weights.astype(np.float64)
    count = source_weights.size
    correlation = float(np.corrcoef(source_float, target_float)[0, 1])
    threshold_rows = []
    for threshold in thresholds:
        source_low = source_weights <= threshold
        target_low = target_weights <= threshold
        source_probability = float(np.mean(source_low))
        target_probability = float(np.mean(target_low))
        joint_probability = float(np.mean(source_low & target_low))
        independent_probability = source_probability * target_probability
        threshold_rows.append(
            {
                "threshold_inclusive": threshold,
                "source_count": int(np.sum(source_low)),
                "target_count": int(np.sum(target_low)),
                "joint_count": int(np.sum(source_low & target_low)),
                "source_probability": source_probability,
                "target_probability": target_probability,
                "joint_probability": joint_probability,
                "independent_probability_from_sample_marginals": independent_probability,
                "joint_to_independent_ratio": (
                    joint_probability / independent_probability
                    if independent_probability > 0.0
                    else None
                ),
            }
        )
    return {
        "samples": count,
        "source_mean": float(np.mean(source_float)),
        "target_mean": float(np.mean(target_float)),
        "source_standard_deviation": float(np.std(source_float)),
        "target_standard_deviation": float(np.std(target_float)),
        "pearson_correlation": correlation,
        "mean_absolute_weight_difference": float(
            np.mean(np.abs(source_float - target_float))
        ),
        "thresholds": threshold_rows,
    }


def bulk_random_lags(
    sample_count: int,
    batch_size: int,
    seed: int,
    powers: np.ndarray,
    thresholds: tuple[int, ...],
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    joint_histogram = np.zeros((CODE_BITS + 1, CODE_BITS + 1), dtype=np.int64)
    source_sum = 0.0
    target_sum = 0.0
    source_square_sum = 0.0
    target_square_sum = 0.0
    cross_sum = 0.0
    absolute_difference_sum = 0.0
    threshold_counts = {
        threshold: [0, 0, 0] for threshold in thresholds
    }
    completed = 0
    while completed < sample_count:
        count = min(batch_size, sample_count - completed)
        messages = random_nonzero_messages(rng, count)
        indices = rng.integers(0, COEFFICIENT_COUNT, size=count, dtype=np.int64)
        targets = field_multiply_batch(powers[indices], messages)
        source_weights = encode_weights(messages)
        target_weights = encode_weights(targets)
        np.add.at(joint_histogram, (source_weights, target_weights), 1)

        source_float = source_weights.astype(np.float64)
        target_float = target_weights.astype(np.float64)
        source_sum += float(np.sum(source_float))
        target_sum += float(np.sum(target_float))
        source_square_sum += float(np.dot(source_float, source_float))
        target_square_sum += float(np.dot(target_float, target_float))
        cross_sum += float(np.dot(source_float, target_float))
        absolute_difference_sum += float(np.sum(np.abs(source_float - target_float)))
        for threshold in thresholds:
            source_low = source_weights <= threshold
            target_low = target_weights <= threshold
            counts = threshold_counts[threshold]
            counts[0] += int(np.sum(source_low))
            counts[1] += int(np.sum(target_low))
            counts[2] += int(np.sum(source_low & target_low))
        completed += count

    source_mean = source_sum / sample_count
    target_mean = target_sum / sample_count
    source_variance = source_square_sum / sample_count - source_mean * source_mean
    target_variance = target_square_sum / sample_count - target_mean * target_mean
    covariance = cross_sum / sample_count - source_mean * target_mean
    correlation = covariance / math.sqrt(source_variance * target_variance)
    threshold_rows = []
    for threshold in thresholds:
        source_count, target_count, joint_count = threshold_counts[threshold]
        source_probability = source_count / sample_count
        target_probability = target_count / sample_count
        joint_probability = joint_count / sample_count
        independent_probability = source_probability * target_probability
        threshold_rows.append(
            {
                "threshold_inclusive": threshold,
                "source_count": source_count,
                "target_count": target_count,
                "joint_count": joint_count,
                "source_probability": source_probability,
                "target_probability": target_probability,
                "joint_probability": joint_probability,
                "independent_probability_from_sample_marginals": independent_probability,
                "joint_to_independent_ratio": (
                    joint_probability / independent_probability
                    if independent_probability > 0.0
                    else None
                ),
            }
        )

    nonzero_bins = np.argwhere(joint_histogram > 0)
    top_bins = sorted(
        (
            {
                "source_weight": int(source_weight),
                "target_weight": int(target_weight),
                "count": int(joint_histogram[source_weight, target_weight]),
            }
            for source_weight, target_weight in nonzero_bins
        ),
        key=lambda item: item["count"],
        reverse=True,
    )[:20]
    return {
        "law": "uniform nonzero x and uniform i in [0,16384)",
        "samples": sample_count,
        "seed": seed,
        "source_mean": source_mean,
        "target_mean": target_mean,
        "source_standard_deviation": math.sqrt(source_variance),
        "target_standard_deviation": math.sqrt(target_variance),
        "pearson_correlation": correlation,
        "mean_absolute_weight_difference": absolute_difference_sum / sample_count,
        "thresholds": threshold_rows,
        "top_joint_bins": top_bins,
    }


def fixed_lag_profiles(
    sample_count: int,
    seed: int,
    powers: np.ndarray,
    lags: tuple[int, ...],
    thresholds: tuple[int, ...],
) -> list[dict[str, object]]:
    rng = np.random.default_rng(seed)
    messages = random_nonzero_messages(rng, sample_count)
    source_weights = encode_weights(messages)
    profiles = []
    for lag in lags:
        if not 0 <= lag < COEFFICIENT_COUNT:
            raise ValueError(f"fixed lag {lag} lies outside the coefficient schedule")
        if lag == 0:
            targets = messages
        else:
            scalar = np.full(sample_count, powers[lag], dtype=np.uint64)
            targets = field_multiply_batch(scalar, messages)
        target_weights = encode_weights(targets)
        profile = summarize_pairs(source_weights, target_weights, thresholds)
        profile["lag"] = lag
        profiles.append(profile)
    return profiles


def read_authenticated_low_messages(
    path: Path, exhaustive_support: int
) -> dict[int, tuple[int, int]]:
    data = path.read_bytes()
    if len(data) % 9:
        raise ValueError("authenticated low-message file has a partial record")
    result: dict[int, tuple[int, int]] = {}
    for offset in range(0, len(data), 9):
        support, message = struct.unpack_from("<BQ", data, offset)
        if support > exhaustive_support:
            raise ValueError("record exceeds the declared exhaustive support")
        codeword = dp.encode_message(message)
        actual_support = sum(
            bool((codeword >> (PACKET_BITS * packet)) & 0xF)
            for packet in range(PACKET_COUNT)
        )
        if actual_support != support:
            raise ValueError("authenticated low-message record failed support check")
        if message in result:
            raise ValueError("authenticated low-message file repeats a message")
        result[message] = (support, codeword.bit_count())
    return result


def exact_low_tail_scan(
    records: dict[int, tuple[int, int]],
    source_support_maximum: int,
    source_binary_weight: int | None,
    lag_start: int,
    lag_count: int,
) -> dict[str, object]:
    sources = [
        (message, support, binary_weight)
        for message, (support, binary_weight) in records.items()
        if support <= source_support_maximum
        and (source_binary_weight is None or binary_weight == source_binary_weight)
    ]
    hits_by_lag = [0] * lag_count
    weight22_hits_by_lag = [0] * lag_count
    support11_hits_by_lag = [0] * lag_count
    pair_histogram: Counter[tuple[int, int, int, int]] = Counter()
    lookup = records.get
    multiply_x = dp.field_multiply_x
    starting_power = dp.field_power(2, lag_start)
    total_hits = 0
    weight22_hits = 0
    support11_hits = 0

    for source, source_support, source_weight in sources:
        target_message = dp.field_multiply(starting_power, source)
        for lag_offset in range(lag_count):
            lag = lag_start + lag_offset
            target = lookup(target_message)
            if target is not None:
                target_support, target_weight = target
                hits_by_lag[lag_offset] += 1
                total_hits += 1
                pair_histogram[
                    (source_support, target_support, source_weight, target_weight)
                ] += 1
                if source_weight == 22 and target_weight == 22:
                    weight22_hits_by_lag[lag_offset] += 1
                    weight22_hits += 1
                if source_support == 11 and target_support == 11:
                    support11_hits_by_lag[lag_offset] += 1
                    support11_hits += 1
            target_message = multiply_x(target_message)

    top_lags = sorted(
        (
            {
                "lag": lag_start + lag_offset,
                "authenticated_target_hits": hits,
                "weight22_to_weight22_hits": weight22_hits_by_lag[lag_offset],
                "support11_to_support11_hits": support11_hits_by_lag[lag_offset],
            }
            for lag_offset, hits in enumerate(hits_by_lag)
            if hits
        ),
        key=lambda item: (
            item["authenticated_target_hits"],
            item["weight22_to_weight22_hits"],
        ),
        reverse=True,
    )[:30]
    top_pairs = [
        {
            "source_packet_support": key[0],
            "target_packet_support": key[1],
            "source_binary_weight": key[2],
            "target_binary_weight": key[3],
            "count_across_sources_and_lags": count,
        }
        for key, count in pair_histogram.most_common(30)
    ]
    source_support_counts = Counter(support for _, support, _ in sources)
    source_weight_counts = Counter(weight for _, _, weight in sources)
    return {
        "source_support_maximum": source_support_maximum,
        "source_binary_weight_filter": source_binary_weight,
        "source_count": len(sources),
        "target_authenticated_count": len(records),
        "coefficient_lag_start": lag_start,
        "coefficient_lags_scanned_per_source": lag_count,
        "source_support_counts": sorted(source_support_counts.items()),
        "source_binary_weight_counts": sorted(source_weight_counts.items()),
        "authenticated_target_hits": total_hits,
        "weight22_to_weight22_hits": weight22_hits,
        "support11_to_support11_hits": support11_hits,
        "lags_with_authenticated_hits": sum(hit > 0 for hit in hits_by_lag),
        "lags_with_weight22_hits": sum(hit > 0 for hit in weight22_hits_by_lag),
        "authenticated_hit_lags": [
            lag_start + lag_offset
            for lag_offset, hit in enumerate(hits_by_lag)
            if hit > 0
        ],
        "weight22_hit_lags": [
            lag_start + lag_offset
            for lag_offset, hit in enumerate(weight22_hits_by_lag)
            if hit > 0
        ],
        "maximum_authenticated_hit_lag": max(
            (
                lag_start + lag_offset
                for lag_offset, hit in enumerate(hits_by_lag)
                if hit > 0
            ),
            default=None,
        ),
        "maximum_weight22_hit_lag": max(
            (
                lag_start + lag_offset
                for lag_offset, hit in enumerate(weight22_hits_by_lag)
                if hit > 0
            ),
            default=None,
        ),
        "top_lags": top_lags,
        "top_pair_profiles": top_pairs,
        "scope": (
            "Exact for sources at or below the selected packet support and "
            "targets inside the supplied exhaustive packet-support list. It does "
            "not exclude target words above that support."
        ),
    }


def validation_gate(powers: np.ndarray) -> str:
    messages = np.asarray(
        [1, 2, 3, 0x123456789ABCDEF0, 0xFFFFFFFFFFFFFFFF], dtype=np.uint64
    )
    weights = encode_weights(messages)
    for message, weight in zip(messages.tolist(), weights.tolist()):
        if weight != dp.encode_message(message).bit_count():
            raise RuntimeError("batched BCH encoder failed scalar validation")
    indices = np.asarray([0, 1, 31, 1023, 16383], dtype=np.int64)
    products = field_multiply_batch(powers[indices], messages)
    for index, message, product in zip(
        indices.tolist(), messages.tolist(), products.tolist()
    ):
        expected = dp.field_multiply(int(powers[index]), message)
        if product != expected:
            raise RuntimeError("batched field multiplication failed scalar validation")
    return "PASS"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-samples", type=int, default=1_000_000)
    parser.add_argument("--batch-size", type=int, default=100_000)
    parser.add_argument("--fixed-samples", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20_260_823)
    parser.add_argument(
        "--fixed-lags", type=parse_int_list, default=DEFAULT_FIXED_LAGS
    )
    parser.add_argument(
        "--thresholds", type=parse_int_list, default=DEFAULT_THRESHOLDS
    )
    parser.add_argument("--low-message-file", type=Path)
    parser.add_argument("--low-exhaustive-support", type=int, default=0)
    parser.add_argument("--low-source-support", type=int, default=12)
    parser.add_argument("--low-source-binary-weight", type=int)
    parser.add_argument("--exact-lag-start", type=int, default=0)
    parser.add_argument("--exact-lag-count", type=int, default=COEFFICIENT_COUNT)
    args = parser.parse_args()

    if args.random_samples < 0 or args.fixed_samples < 0:
        raise ValueError("sample counts must be nonnegative")
    if args.batch_size < 1:
        raise ValueError("batch size must be positive")
    if bool(args.low_message_file) != bool(args.low_exhaustive_support):
        raise ValueError("low-message file and exhaustive support are required together")
    if args.low_source_support < 1:
        raise ValueError("low source support must be positive")
    if args.low_source_binary_weight is not None and not 1 <= args.low_source_binary_weight <= CODE_BITS:
        raise ValueError("low source binary weight lies outside the BCH block")
    if args.exact_lag_start < 0 or args.exact_lag_count < 1:
        raise ValueError("exact lag range must be nonnegative and nonempty")

    powers = coefficient_powers()
    validation = validation_gate(powers)
    random_lags = (
        bulk_random_lags(
            args.random_samples,
            args.batch_size,
            args.seed,
            powers,
            args.thresholds,
        )
        if args.random_samples
        else None
    )
    fixed_lags = (
        fixed_lag_profiles(
            args.fixed_samples,
            args.seed ^ 0xA17A,
            powers,
            args.fixed_lags,
            args.thresholds,
        )
        if args.fixed_samples
        else None
    )
    low_tail = None
    if args.low_message_file:
        records = read_authenticated_low_messages(
            args.low_message_file, args.low_exhaustive_support
        )
        low_tail = exact_low_tail_scan(
            records,
            args.low_source_support,
            args.low_source_binary_weight,
            args.exact_lag_start,
            args.exact_lag_count,
        )

    payload = {
        "schema": "riffle-bch-alpha-joint-spectrum-v1",
        "evidence_label": "BULK_MONTE_CARLO_AND_EXACT_AUTHENTICATED_LOW_TAIL",
        "construction": {
            "field": "GF(2^64) with reduction polynomial X^64+X^4+X^3+X+1",
            "coefficient_schedule": "alpha_i = gamma^i for 0 <= i < 16384",
            "bch": "extended binary BCH [128,64,22]",
        },
        "validation": validation,
        "random_lag_sample": random_lags,
        "fixed_lag_samples": fixed_lags,
        "authenticated_low_tail": low_tail,
        "scope": (
            "Monte Carlo rows are numerical evidence only. Exact low-tail rows "
            "are complete only within the declared packet-support ranges."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
