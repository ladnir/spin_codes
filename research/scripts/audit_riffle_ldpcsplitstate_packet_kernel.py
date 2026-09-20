#!/usr/bin/env python3
"""Audit packet-support cancellations in a fixed LDPCSplitState compressor."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from explore_riffle_ldpcsplitstate_local import gf2_rank


DEFAULT_INPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/fixed_nested_pair_depth2_search.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/fixed_nested_pair_depth2_packet_kernel.json"
)


def packet_columns(source: dict[str, object]) -> list[list[int]]:
    s_columns = [int(value, 16) for value in source["matrices"]["S_columns_hex_64"]]
    r_columns = [int(value, 16) for value in source["matrices"]["R_columns_hex_64"]]
    b_columns = [int(value, 16) for value in source["matrices"]["B_columns_hex_64"]]
    identity_columns = b_columns[192:]
    return [
        [
            s_columns[slot],
            r_columns[slot],
            r_columns[64 + slot],
            identity_columns[slot],
        ]
        for slot in range(64)
    ]


def symbol_syndromes(columns: list[list[int]]) -> list[list[tuple[int, int]]]:
    result = []
    for packet in columns:
        symbols = []
        for pattern in range(1, 16):
            syndrome = 0
            for lane in range(4):
                if (pattern >> lane) & 1:
                    syndrome ^= packet[lane]
            symbols.append((pattern, syndrome))
        result.append(symbols)
    return result


def audit(source: dict[str, object]) -> dict[str, object]:
    columns = packet_columns(source)
    symbols = symbol_syndromes(columns)
    ranks = [gf2_rank(packet, 64) for packet in columns]
    zero_symbols = [
        {"slot": slot, "pattern": pattern}
        for slot, packet_symbols in enumerate(symbols)
        for pattern, syndrome in packet_symbols
        if syndrome == 0
    ]

    by_syndrome: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for slot, packet_symbols in enumerate(symbols):
        for pattern, syndrome in packet_symbols:
            by_syndrome[syndrome].append((slot, pattern))
    pair_witnesses = []
    for syndrome, occurrences in by_syndrome.items():
        for left in range(len(occurrences)):
            for right in range(left + 1, len(occurrences)):
                if occurrences[left][0] != occurrences[right][0]:
                    pair_witnesses.append(
                        {
                            "syndrome_hex": f"{syndrome:016x}",
                            "left_slot": occurrences[left][0],
                            "left_pattern": occurrences[left][1],
                            "right_slot": occurrences[right][0],
                            "right_pattern": occurrences[right][1],
                        }
                    )

    triple_witnesses: set[tuple[int, int, int, int, int, int]] = set()
    for left_slot in range(64):
        for right_slot in range(left_slot + 1, 64):
            for left_pattern, left_syndrome in symbols[left_slot]:
                for right_pattern, right_syndrome in symbols[right_slot]:
                    target = left_syndrome ^ right_syndrome
                    for third_slot, third_pattern in by_syndrome.get(target, ()):
                        if third_slot > right_slot:
                            triple_witnesses.add(
                                (
                                    left_slot,
                                    left_pattern,
                                    right_slot,
                                    right_pattern,
                                    third_slot,
                                    third_pattern,
                                )
                            )

    flat_columns = [column for packet in columns for column in packet]
    pairs_by_syndrome: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for left in range(256):
        for right in range(left + 1, 256):
            pairs_by_syndrome[flat_columns[left] ^ flat_columns[right]].append(
                (left, right)
            )
    weight_four_supports: set[tuple[int, int, int, int]] = set()
    for pairs in pairs_by_syndrome.values():
        for first in range(len(pairs)):
            for second in range(first + 1, len(pairs)):
                support = (*pairs[first], *pairs[second])
                if len(set(support)) == 4:
                    weight_four_supports.add(tuple(sorted(support)))
    packet_support_histogram: dict[int, int] = defaultdict(int)
    lane_composition_histogram: dict[tuple[int, int, int, int], int] = defaultdict(int)
    four_packet_witnesses = []
    for support in sorted(weight_four_supports):
        packet_count = len({coordinate // 4 for coordinate in support})
        packet_support_histogram[packet_count] += 1
        if packet_count == 4 and len(four_packet_witnesses) < 32:
            four_packet_witnesses.append(
                [
                    {
                        "slot": coordinate // 4,
                        "lane": coordinate % 4,
                    }
                    for coordinate in support
                ]
            )
        if packet_count == 4:
            lane_counts = [0, 0, 0, 0]
            for coordinate in support:
                lane_counts[coordinate % 4] += 1
            lane_composition_histogram[tuple(lane_counts)] += 1

    placement_rows = []
    placement_denominator = 2048 * 2047 * 2046 * 2045
    for composition, kernel_words in sorted(lane_composition_histogram.items()):
        compatible_labelings = math.prod(math.factorial(count) for count in composition)
        probability = (
            32 * kernel_words * compatible_labelings / placement_denominator
        )
        placement_rows.append(
            {
                "lane_composition": list(composition),
                "kernel_words": kernel_words,
                "compatible_labelings_per_kernel_word": compatible_labelings,
                "failure_probability": probability,
                "failure_bits": -math.log2(probability),
            }
        )
    worst_singleton_placement = max(
        placement_rows, key=lambda row: row["failure_probability"]
    )

    minimum_packet_support_lower_bound = 1
    if not zero_symbols:
        minimum_packet_support_lower_bound = 2
    if minimum_packet_support_lower_bound == 2 and not pair_witnesses:
        minimum_packet_support_lower_bound = 3
    if minimum_packet_support_lower_bound == 3 and not triple_witnesses:
        minimum_packet_support_lower_bound = 4
    exact_minimum_packet_support = (
        4
        if minimum_packet_support_lower_bound == 4
        and packet_support_histogram.get(4, 0) > 0
        else None
    )

    return {
        "schema": "riffle-ldpcsplitstate-packet-kernel-audit-v1",
        "candidate": source["candidate"],
        "source_chosen_seed": source["parameters"]["chosen_seed"],
        "packet_layout": (
            "slot k uses columns S[k], R[k], R[64+k], and identity[k] "
            "for lanes zero through three"
        ),
        "packet_subspace_ranks": {
            "minimum": min(ranks),
            "maximum": max(ranks),
            "histogram": {
                str(rank): ranks.count(rank) for rank in sorted(set(ranks))
            },
        },
        "zero_nonzero_packet_symbols": {
            "count": len(zero_symbols),
            "witnesses": zero_symbols[:32],
        },
        "two_packet_cancellations": {
            "count": len(pair_witnesses),
            "witnesses": pair_witnesses[:32],
        },
        "three_packet_cancellations": {
            "count": len(triple_witnesses),
            "witnesses": [
                {
                    "first_slot": witness[0],
                    "first_pattern": witness[1],
                    "second_slot": witness[2],
                    "second_pattern": witness[3],
                    "third_slot": witness[4],
                    "third_pattern": witness[5],
                }
                for witness in sorted(triple_witnesses)[:32]
            ],
        },
        "coordinate_weight_four_kernel": {
            "count": len(weight_four_supports),
            "packet_support_histogram": {
                str(key): packet_support_histogram[key]
                for key in sorted(packet_support_histogram)
            },
            "lane_composition_histogram": {
                ",".join(map(str, key)): lane_composition_histogram[key]
                for key in sorted(lane_composition_histogram)
            },
            "four_packet_witnesses": four_packet_witnesses,
        },
        "four_singleton_packet_placement": {
            "packet_slots": 2048,
            "epochs": 32,
            "slots_per_epoch": 64,
            "by_lane_composition": placement_rows,
            "worst_lane_composition": worst_singleton_placement,
            "probability_space": (
                "four labeled singleton packets with fixed lanes under one "
                "uniform permutation of 2048 packet groups"
            ),
        },
        "proved_minimum_nonactivation_packet_support_lower_bound": (
            minimum_packet_support_lower_bound
        ),
        "exact_minimum_nonactivation_packet_support": exact_minimum_packet_support,
        "scope": (
            "Exact exhaustive audit through three distinct packet slots for "
            "the recorded fixed B. A lower bound of four does not assert that "
            "a four-packet cancellation exists."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source = json.loads(args.input.read_text(encoding="utf-8"))
    payload = audit(source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "packet_rank_min,"
        f"{payload['packet_subspace_ranks']['minimum']},"
        "zero_symbols,"
        f"{payload['zero_nonzero_packet_symbols']['count']}"
    )
    print(
        "two_packet_cancellations,"
        f"{payload['two_packet_cancellations']['count']}"
    )
    print(
        "three_packet_cancellations,"
        f"{payload['three_packet_cancellations']['count']}"
    )
    print(
        "packet_support_lower_bound,"
        f"{payload['proved_minimum_nonactivation_packet_support_lower_bound']}"
    )
    print(
        "exact_packet_support_minimum,"
        f"{payload['exact_minimum_nonactivation_packet_support']}"
    )
    worst = payload["four_singleton_packet_placement"]["worst_lane_composition"]
    print(
        "four_singleton_worst_failure_bits,"
        f"{worst['failure_bits']:.6f},lane_composition,"
        + ",".join(map(str, worst["lane_composition"]))
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
