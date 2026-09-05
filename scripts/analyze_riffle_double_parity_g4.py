#!/usr/bin/env python3
"""First proof diagnostics for the sequential Riffle DP g=4 candidate.

The candidate has 16,384 data symbols over GF(2^64), two field-valued
parity symbols, and one extended BCH [128,64,22] word per outer symbol.
The outer BCH words remain sequential.  A uniform permutation acts on the
resulting four-bit packets before the recursive inner transform.

The exact calculations cover the terminal-placement obstruction and a
small, authenticated part of the local BCH tail.  Optional Monte Carlo
sampling estimates the four block-weight-three support categories.  A
sampled estimate is not an upper bound or an end-to-end certificate.
"""

from __future__ import annotations

import argparse
import math
import struct
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as dp  # noqa: E402
import analyze_riffle_g8_late_suffix as base  # noqa: E402


PACKET_BITS = 4
PACKETS_PER_BCH_WORD = base.LOCAL_LENGTH // PACKET_BITS
OUTER_BLOCKS = dp.USER_BLOCKS + 2
CODE_LENGTH = OUTER_BLOCKS * base.LOCAL_LENGTH
DISTANCE_THRESHOLD = 9 * CODE_LENGTH // 100
INNER_NODES = CODE_LENGTH // 64
TERMINAL_NODES = DISTANCE_THRESHOLD // 64
PACKET_COUNT = CODE_LENGTH // PACKET_BITS
TERMINAL_PACKETS = TERMINAL_NODES * (64 // PACKET_BITS)


def terminal_probability(packet_support: int) -> Fraction:
    """Return the probability that all active packets enter terminal nodes."""
    if packet_support < 0 or packet_support > TERMINAL_PACKETS:
        return Fraction()
    return Fraction(
        math.comb(TERMINAL_PACKETS, packet_support),
        math.comb(PACKET_COUNT, packet_support),
    )


def mds_terminal_lower_bound(block_weight: int) -> Fraction:
    """Lower-bound the terminal contribution from one complete MDS shell."""
    maximum_support = PACKETS_PER_BCH_WORD * block_weight
    return dp.mds_block_weight_count(block_weight) * terminal_probability(
        maximum_support
    )


def nibble_support(codeword: int) -> int:
    return sum(
        bool((codeword >> (PACKET_BITS * packet)) & 0xF)
        for packet in range(PACKETS_PER_BCH_WORD)
    )


def authenticated_local_tail() -> tuple[Counter[int], int]:
    """Map every BCH word on at most eight bytes to its nibble support.

    Any BCH word on at most eight nibbles also lies in this enumerated set.
    The absence of such a word therefore proves nibble support at least nine.
    Counts above eight are only a partial nibble-support enumerator.
    """
    messages = dp.exact_low_support_messages(dp.MAX_EXACT_LOCAL_PACKET_SUPPORT)
    expected = sum(
        base.packet_support_enumerator()[
            1 : dp.MAX_EXACT_LOCAL_PACKET_SUPPORT + 1
        ]
    )
    if len(messages) != expected:
        raise RuntimeError("authenticated byte-support enumeration has wrong mass")
    counts = Counter(
        nibble_support(dp.encode_message(message)) for message in messages
    )
    return counts, len(messages)


def read_exact_low_messages(path: Path) -> dict[int, set[int]]:
    """Read support-byte/message records emitted by the C++ enumerator."""
    data = path.read_bytes()
    if len(data) % 9:
        raise ValueError("low-message file has a partial record")
    messages: dict[int, set[int]] = {}
    for offset in range(0, len(data), 9):
        support, message = struct.unpack_from("<BQ", data, offset)
        actual_support = nibble_support(dp.encode_message(message))
        if support != actual_support:
            raise ValueError("low-message record failed BCH support verification")
        messages.setdefault(support, set()).add(message)
    return messages


def multiply_x_inverse(value: int) -> int:
    if value & 1:
        return ((value ^ dp.FIELD_REDUCTION) >> 1) | (1 << 63)
    return value >> 1


def exact_support_33_tail(
    messages_by_support: dict[int, set[int]],
) -> tuple[dict[str, int], Fraction]:
    """Count every block-weight-three DP word on at most 33 nibbles.

    The input must establish local nibble distance 11 and contain every word
    at support 11.  Three nonzero local words then have total support at least
    33.  Equality requires all three values to belong to the supplied set.
    """
    if any(messages_by_support.get(support) for support in range(1, 11)):
        raise ValueError("unexpected BCH word below nibble support 11")
    low = messages_by_support.get(11, set())
    if not low:
        raise ValueError("support-11 BCH set is empty")

    powers = [1]
    for _ in range(1, dp.USER_BLOCKS):
        powers.append(dp.field_multiply_x(powers[-1]))
    power_log = {value: exponent for exponent, value in enumerate(powers)}

    one = 0
    for parameter in low:
        inverse = dp.field_inverse(parameter)
        for parity1 in low:
            if dp.field_multiply(parity1, inverse) in power_log:
                one += 1

    inverse_power_gap: dict[int, int] = {}
    inverse_power = 1
    for gap in range(1, dp.USER_BLOCKS):
        inverse_power = multiply_x_inverse(inverse_power)
        inverse_power_gap[inverse_power] = gap

    pair_p0 = 0
    for value0 in low:
        inverse0 = dp.field_inverse(value0)
        for value1 in low:
            if (value0 ^ value1) not in low:
                continue
            gap = inverse_power_gap.get(dp.field_multiply(value1, inverse0))
            if gap is not None:
                pair_p0 += dp.USER_BLOCKS - gap

    ratio_counts: Counter[int] = Counter()
    for parameter in low:
        inverse = dp.field_inverse(parameter)
        for parity1 in low:
            ratio_counts[dp.field_multiply(parity1, inverse)] += 1
    one_plus_power_gap = {
        powers[gap] ^ 1: gap for gap in range(1, dp.USER_BLOCKS)
    }
    pair_p1 = 0
    for ratio, multiplicity in ratio_counts.items():
        shifted = ratio
        for index0 in range(dp.USER_BLOCKS - 1):
            gap = one_plus_power_gap.get(shifted)
            if gap is not None and index0 + gap < dp.USER_BLOCKS:
                pair_p1 += multiplicity
            shifted = multiply_x_inverse(shifted)

    scaled_power_logs: dict[int, dict[int, int]] = {}
    for value in low:
        scaled = value
        logs: dict[int, int] = {}
        for exponent in range(dp.USER_BLOCKS):
            logs[scaled] = exponent
            scaled = dp.field_multiply_x(scaled)
        scaled_power_logs[value] = logs

    triple = 0
    for value0 in low:
        for value1 in low:
            value2 = value0 ^ value1
            if value2 not in low:
                continue
            shifted1 = value1
            logs2 = scaled_power_logs[value2]
            for gap1 in range(1, dp.USER_BLOCKS - 1):
                shifted1 = dp.field_multiply_x(shifted1)
                gap2 = logs2.get(value0 ^ shifted1)
                if gap2 is not None and gap2 > gap1:
                    triple += dp.USER_BLOCKS - gap2

    counts = {
        "one": one,
        "pair_p0": pair_p0,
        "pair_p1": pair_p1,
        "triple": triple,
    }
    contribution = sum(counts.values()) * terminal_probability(33)
    return counts, contribution


def weight_three_layer_upper_bound(
    messages_by_support: dict[int, set[int]], total_support: int
) -> tuple[int, Fraction]:
    """Upper-bound one block-weight-three layer from exact local counts."""
    ordered_value_triples = 0
    for support0 in range(11, total_support - 21):
        for support1 in range(11, total_support - support0 - 10):
            support2 = total_support - support0 - support1
            if support2 < 11:
                continue
            ordered_value_triples += (
                len(messages_by_support.get(support0, set()))
                * len(messages_by_support.get(support1, set()))
                * len(messages_by_support.get(support2, set()))
            )

    # For a fixed ordered value triple, at most C(B,2) three-data supports,
    # B supports of each two-data category, and one one-data support realize it.
    maximum_outer_placements = math.comb(dp.USER_BLOCKS, 2) + 2 * dp.USER_BLOCKS + 1
    codeword_bound = ordered_value_triples * maximum_outer_placements
    return codeword_bound, codeword_bound * terminal_probability(total_support)


def sample_weight_three_category(
    category: str,
    sample_count: int,
    rng: np.random.Generator,
    batch_size: int = 50_000,
) -> tuple[Counter[int], float]:
    """Estimate one MDS block-weight-three terminal contribution."""
    generator_low = np.asarray(
        [row & dp.FIELD_ORDER for row in dp.GENERATOR_ROWS], dtype=np.uint64
    )
    generator_high = np.asarray(
        [row >> 64 for row in dp.GENERATOR_ROWS], dtype=np.uint64
    )
    maximum = np.uint64(dp.FIELD_ORDER)

    powers = np.empty(dp.USER_BLOCKS, dtype=np.uint64)
    powers[0] = 1
    for index in range(1, dp.USER_BLOCKS):
        powers[index] = np.uint64(dp.field_multiply_x(int(powers[index - 1])))

    def multiply(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        result = np.zeros(left.size, dtype=np.uint64)
        left = left.copy()
        right = right.copy()
        for _ in range(64):
            result ^= np.where(right & np.uint64(1), left, np.uint64(0))
            left = (left << np.uint64(1)) ^ np.where(
                left >> np.uint64(63),
                np.uint64(dp.FIELD_REDUCTION),
                np.uint64(0),
            )
            right >>= np.uint64(1)
        return result

    def support(messages: np.ndarray) -> np.ndarray:
        low = np.zeros(messages.size, dtype=np.uint64)
        high = np.zeros(messages.size, dtype=np.uint64)
        for bit in range(64):
            active = (messages >> np.uint64(bit)) & np.uint64(1)
            mask = np.where(active != 0, maximum, np.uint64(0))
            low ^= mask & generator_low[bit]
            high ^= mask & generator_high[bit]
        low_bytes = low.view(np.uint8).reshape(-1, 8)
        high_bytes = high.view(np.uint8).reshape(-1, 8)
        return (
            np.count_nonzero(low_bytes & 0xF, axis=1)
            + np.count_nonzero(low_bytes >> 4, axis=1)
            + np.count_nonzero(high_bytes & 0xF, axis=1)
            + np.count_nonzero(high_bytes >> 4, axis=1)
        )

    category_size = {
        "one": 1,
        "pair_p0": 2,
        "pair_p1": 2,
        "triple": 3,
    }[category]
    probabilities = np.asarray(
        [float(terminal_probability(packet_support)) for packet_support in range(97)]
    )
    histogram: Counter[int] = Counter()
    probability_sum = 0.0
    completed = 0
    while completed < sample_count:
        count = min(batch_size, sample_count - completed)
        parameter = rng.integers(1, 1 << 64, size=count, dtype=np.uint64)
        indices = rng.integers(
            0, dp.USER_BLOCKS, size=(count, category_size), dtype=np.int64
        )
        if category_size > 1:
            collision = np.any(
                np.diff(np.sort(indices, axis=1), axis=1) == 0, axis=1
            )
            while np.any(collision):
                indices[collision] = rng.integers(
                    0,
                    dp.USER_BLOCKS,
                    size=(int(np.sum(collision)), category_size),
                )
                collision = np.any(
                    np.diff(np.sort(indices, axis=1), axis=1) == 0, axis=1
                )
            indices.sort(axis=1)
        alpha = [powers[indices[:, position]] for position in range(category_size)]

        if category == "one":
            messages = (parameter, parameter, multiply(alpha[0], parameter))
        elif category == "pair_p0":
            messages = (
                multiply(alpha[1], parameter),
                multiply(alpha[0], parameter),
                multiply(alpha[0] ^ alpha[1], parameter),
            )
        elif category == "pair_p1":
            messages = (
                parameter,
                parameter,
                multiply(alpha[0] ^ alpha[1], parameter),
            )
        else:
            messages = (
                multiply(alpha[1] ^ alpha[2], parameter),
                multiply(alpha[2] ^ alpha[0], parameter),
                multiply(alpha[0] ^ alpha[1], parameter),
            )

        total_support = sum(support(message) for message in messages)
        values, counts = np.unique(total_support, return_counts=True)
        for packet_support, frequency in zip(values.tolist(), counts.tolist()):
            histogram[int(packet_support)] += int(frequency)
            probability_sum += frequency * probabilities[int(packet_support)]
        completed += count

    return histogram, probability_sum / sample_count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samples",
        type=int,
        default=0,
        help="sample this many words from each block-weight-three category",
    )
    parser.add_argument(
        "--skip-local-tail",
        action="store_true",
        help="skip the exact enumeration of BCH words on at most eight bytes",
    )
    parser.add_argument(
        "--exact-low-message-file",
        type=Path,
        help=(
            "binary support/message records from enumerate_bch_g4_low_support; "
            "compute the exact minimum block-weight-three tail"
        ),
    )
    parser.add_argument(
        "--exact-low-max-support",
        type=int,
        default=0,
        help="largest support exhaustively scanned in the exact low-message file",
    )
    args = parser.parse_args()
    if args.samples < 0:
        raise SystemExit("sample count must be nonnegative")
    if bool(args.exact_low_message_file) != bool(args.exact_low_max_support):
        raise SystemExit(
            "exact low-message file and its maximum support must be supplied together"
        )

    print("Riffle DP g=4 first proof diagnostic")
    print(
        f"N={CODE_LENGTH} d={DISTANCE_THRESHOLD} inner_nodes={INNER_NODES} "
        f"packets={PACKET_COUNT} terminal_packets={TERMINAL_PACKETS}"
    )
    print("double_mds_block_weight,max_packet_terminal_lower_log2")
    for block_weight in range(3, 13):
        print(
            f"{block_weight},"
            f"{base.log2_fraction(mds_terminal_lower_bound(block_weight)):.12f}"
        )

    if not args.skip_local_tail:
        counts, message_count = authenticated_local_tail()
        print(f"exact_bch_messages_byte_support_le_8={message_count}")
        print(
            "partial_nibble_support_counts="
            + ",".join(f"{support}:{count}" for support, count in sorted(counts.items()))
        )
        if counts and min(counts) <= 8:
            raise RuntimeError("found a BCH word on at most eight nibbles")
        print("proved_nonzero_bch_nibble_support_at_least=9")

    if args.exact_low_message_file:
        messages_by_support = read_exact_low_messages(args.exact_low_message_file)
        counts, contribution = exact_support_33_tail(messages_by_support)
        print(
            "exact_low_message_counts="
            + ",".join(
                f"{support}:{len(messages)}"
                for support, messages in sorted(messages_by_support.items())
            )
        )
        print(
            "exact_weight_three_support_33_counts="
            + ",".join(f"{category}:{count}" for category, count in counts.items())
        )
        print(
            "exact_weight_three_support_33_total="
            f"{sum(counts.values())}"
        )
        print(
            "exact_weight_three_support_33_terminal_log2="
            f"{base.log2_fraction(contribution):.12f}"
        )
        for total_support in range(34, args.exact_low_max_support + 23):
            codeword_bound, contribution_bound = weight_three_layer_upper_bound(
                messages_by_support, total_support
            )
            print(
                f"weight_three_support_{total_support}_codeword_upper="
                f"{codeword_bound}"
            )
            print(
                f"weight_three_support_{total_support}_terminal_upper_log2="
                f"{base.log2_fraction(contribution_bound):.12f}"
            )

    if args.samples:
        rng = np.random.default_rng(0xD04)
        category_counts = {
            "one": dp.USER_BLOCKS,
            "pair_p0": math.comb(dp.USER_BLOCKS, 2),
            "pair_p1": math.comb(dp.USER_BLOCKS, 2),
            "triple": math.comb(dp.USER_BLOCKS, 3),
        }
        print(f"samples_per_weight_three_category={args.samples}")
        for category, outer_support_count in category_counts.items():
            histogram, sampled_mean = sample_weight_three_category(
                category, args.samples, rng
            )
            estimate = outer_support_count * dp.FIELD_ORDER * sampled_mean
            average_support = sum(
                support * count for support, count in histogram.items()
            ) / args.samples
            print(
                f"sampled_{category}_minimum_support={min(histogram)} "
                f"mean_support={average_support:.6f} "
                f"estimated_terminal_log2={math.log2(estimate):.12f}"
            )

    print("status=DIAGNOSTIC_NOT_END_TO_END_CERTIFICATE")


if __name__ == "__main__":
    main()
