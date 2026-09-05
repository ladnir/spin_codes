#!/usr/bin/env python3
"""Count the even-complement low-shell relation core exactly.

Let ``Lw`` contain the messages whose committed extended BCH encoding has
weight ``w``.  This script counts the ordered relations

    x + y + z = 0

for ``L24^3`` and ``L22 x L24^2``.  The affine coordinate group preserves
the code, Hamming weight, and addition.  A partner count is therefore
constant on each affine orbit of the fixed first word.

The implementation classifies both committed shells into affine orbits.  It
then scans one representative per orbit against a pre-encoded target shell.
The hot scan uses chunked NumPy bitwise popcounts and does not re-encode the
target shell for each representative.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as bch  # noqa: E402
from analyze_ebch128_adjacent_triples import (  # noqa: E402
    affine_message_actions,
    transformed_messages_fast,
)
from analyze_riffle_bch_alpha_joint_spectrum import (  # noqa: E402
    ENCODE_HIGH,
    ENCODE_LOW,
)
from enumerate_ebch128_weight22_affine import (  # noqa: E402
    AFFINE_GROUP_SIZE,
    coordinate_maps,
)


ROOT = SCRIPT_DIRECTORY.parent
CONSTRUCTION = (
    ROOT / "constructions" / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
)
DEFAULT_WEIGHT22 = CONSTRUCTION / "ebch128_weight22_messages.bin"
DEFAULT_WEIGHT24 = CONSTRUCTION / "ebch128_weight24_messages.bin"
EXPECTED_SHELL_SIZES = {22: 243_840, 24: 6_855_968}
EXPECTED_ORBIT_HISTOGRAMS = {
    22: {AFFINE_GROUP_SIZE: 15},
    24: {4_064: 147, AFFINE_GROUP_SIZE: 385},
}


def read_messages(path: Path, *, expected: int) -> np.ndarray:
    messages = np.fromfile(path, dtype="<u8")
    if messages.size != expected:
        raise RuntimeError(
            f"{path} contains {messages.size} messages; expected {expected}"
        )
    if np.any(messages[1:] <= messages[:-1]):
        raise RuntimeError(f"{path} is not sorted and unique")
    return messages


def encode_halves(messages: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Encode a uint64 message array into two uint64 codeword halves."""
    low = np.zeros(messages.size, dtype=np.uint64)
    high = np.zeros(messages.size, dtype=np.uint64)
    for byte_position in range(8):
        values = (
            (messages >> np.uint64(8 * byte_position)) & np.uint64(0xFF)
        ).astype(np.uint8)
        low ^= ENCODE_LOW[byte_position, values]
        high ^= ENCODE_HIGH[byte_position, values]
    return low, high


def classify_affine_orbits(
    messages: np.ndarray,
    low: np.ndarray,
    high: np.ndarray,
    *,
    weight: int,
    actions: np.ndarray,
) -> list[dict[str, int]]:
    """Return one message and the orbit size for every affine orbit."""
    covered = np.zeros(messages.size, dtype=np.bool_)
    representatives: list[dict[str, int]] = []
    cursor = 0
    while cursor < messages.size:
        while cursor < messages.size and covered[cursor]:
            cursor += 1
        if cursor == messages.size:
            break
        message = int(messages[cursor])
        codeword = int(low[cursor]) | (int(high[cursor]) << 64)
        if codeword.bit_count() != weight:
            raise RuntimeError(f"L{weight} contains a word of the wrong weight")
        orbit = np.unique(transformed_messages_fast(codeword, actions))
        indices = np.searchsorted(messages, orbit)
        if np.any(indices == messages.size) or np.any(messages[indices] != orbit):
            raise RuntimeError(f"an affine image left L{weight}")
        if np.any(covered[indices]):
            raise RuntimeError(f"affine orbit overlap detected in L{weight}")
        covered[indices] = True
        representatives.append(
            {
                "message": message,
                "codeword_low": int(low[cursor]),
                "codeword_high": int(high[cursor]),
                "orbit_size": int(orbit.size),
            }
        )
    if not np.all(covered):
        raise RuntimeError(f"affine orbits do not cover L{weight}")
    histogram = dict(sorted(Counter(row["orbit_size"] for row in representatives).items()))
    if histogram != EXPECTED_ORBIT_HISTOGRAMS[weight]:
        raise RuntimeError(
            f"L{weight} orbit histogram changed: {histogram}"
        )
    return representatives


def count_partners(
    representatives: list[dict[str, int]],
    target_low: np.ndarray,
    target_high: np.ndarray,
    *,
    xor_weight: int,
    chunk_size: int,
    label: str,
) -> list[int]:
    """Count target words at the requested XOR weight for each representative."""
    counts: list[int] = []
    started = time.perf_counter()
    for representative_index, representative in enumerate(representatives):
        count = 0
        rep_low = np.uint64(representative["codeword_low"])
        rep_high = np.uint64(representative["codeword_high"])
        for begin in range(0, target_low.size, chunk_size):
            end = min(begin + chunk_size, target_low.size)
            weights = np.bitwise_count(target_low[begin:end] ^ rep_low)
            weights += np.bitwise_count(target_high[begin:end] ^ rep_high)
            count += int(np.count_nonzero(weights == xor_weight))
        counts.append(count)
        if (representative_index + 1) % 16 == 0 or (
            representative_index + 1 == len(representatives)
        ):
            elapsed = time.perf_counter() - started
            print(
                f"{label}: {representative_index + 1}/{len(representatives)} "
                f"representatives in {elapsed:.1f}s",
                file=sys.stderr,
                flush=True,
            )
    return counts


def weighted_relation_count(
    representatives: list[dict[str, int]], partner_counts: list[int]
) -> int:
    if len(representatives) != len(partner_counts):
        raise ValueError("representative and partner counts have different lengths")
    return sum(
        row["orbit_size"] * count
        for row, count in zip(representatives, partner_counts)
    )


def orbit_rows(
    representatives: list[dict[str, int]], partner_counts: list[int]
) -> list[dict[str, object]]:
    return [
        {
            "representative_message_hex": hex(row["message"]),
            "orbit_size": row["orbit_size"],
            "partner_count": count,
        }
        for row, count in zip(representatives, partner_counts)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight22", type=Path, default=DEFAULT_WEIGHT22)
    parser.add_argument("--weight24", type=Path, default=DEFAULT_WEIGHT24)
    parser.add_argument("--chunk-size", type=int, default=1 << 20)
    args = parser.parse_args()
    if args.chunk_size <= 0:
        raise ValueError("chunk size must be positive")

    started = time.perf_counter()
    messages22 = read_messages(
        args.weight22, expected=EXPECTED_SHELL_SIZES[22]
    )
    messages24 = read_messages(
        args.weight24, expected=EXPECTED_SHELL_SIZES[24]
    )
    low22, high22 = encode_halves(messages22)
    low24, high24 = encode_halves(messages24)

    coordinate_to_point, point_to_coordinate, multiplication = coordinate_maps()
    actions = affine_message_actions(
        coordinate_to_point, point_to_coordinate, multiplication
    )
    representatives22 = classify_affine_orbits(
        messages22, low22, high22, weight=22, actions=actions
    )
    representatives24 = classify_affine_orbits(
        messages24, low24, high24, weight=24, actions=actions
    )
    classified_seconds = time.perf_counter() - started

    l22_cube_partners = count_partners(
        representatives22,
        low22,
        high22,
        xor_weight=22,
        chunk_size=args.chunk_size,
        label="L22^3 audit",
    )
    l22_cube = weighted_relation_count(representatives22, l22_cube_partners)
    if l22_cube != 1_365_504:
        raise RuntimeError(f"L22^3 audit changed: {l22_cube}")

    mixed_from22_partners = count_partners(
        representatives22,
        low24,
        high24,
        xor_weight=24,
        chunk_size=args.chunk_size,
        label="L22 x L24^2 from L22",
    )
    mixed_from22 = weighted_relation_count(
        representatives22, mixed_from22_partners
    )

    mixed_from24_partners = count_partners(
        representatives24,
        low22,
        high22,
        xor_weight=24,
        chunk_size=args.chunk_size,
        label="L22 x L24^2 from L24",
    )
    mixed_from24 = weighted_relation_count(
        representatives24, mixed_from24_partners
    )
    if mixed_from22 != mixed_from24:
        raise RuntimeError(
            "the two projections of L22 x L24^2 disagree: "
            f"{mixed_from22} != {mixed_from24}"
        )

    l24_cube_partners = count_partners(
        representatives24,
        low24,
        high24,
        xor_weight=24,
        chunk_size=args.chunk_size,
        label="L24^3",
    )
    if any(count & 1 for count in l24_cube_partners):
        raise RuntimeError("an L24 representative has an odd partner count")
    l24_cube = weighted_relation_count(representatives24, l24_cube_partners)
    if l24_cube % 6:
        raise RuntimeError("the ordered L24^3 count is not divisible by six")

    elapsed = time.perf_counter() - started
    payload = {
        "schema": "ebch128-weight24-relation-core-v1",
        "code": "extended primitive binary BCH [128,64,22]",
        "equation": "x+y+z=0",
        "shells": {
            "L22": {
                "size": int(messages22.size),
                "affine_orbits": len(representatives22),
                "orbit_size_histogram": {
                    str(size): count
                    for size, count in EXPECTED_ORBIT_HISTOGRAMS[22].items()
                },
            },
            "L24": {
                "size": int(messages24.size),
                "affine_orbits": len(representatives24),
                "orbit_size_histogram": {
                    str(size): count
                    for size, count in EXPECTED_ORBIT_HISTOGRAMS[24].items()
                },
            },
        },
        "relations": {
            "L22^3_audit": {
                "ordered_relations": l22_cube,
                "log2_ordered_relations": math.log2(l22_cube),
                "orbit_rows": orbit_rows(
                    representatives22, l22_cube_partners
                ),
            },
            "L22_x_L24^2": {
                "ordered_relations": mixed_from22,
                "log2_ordered_relations": math.log2(mixed_from22),
                "from_L22_orbit_rows": orbit_rows(
                    representatives22, mixed_from22_partners
                ),
                "from_L24_orbit_rows": orbit_rows(
                    representatives24, mixed_from24_partners
                ),
            },
            "L24^3": {
                "ordered_relations": l24_cube,
                "unordered_relations": l24_cube // 6,
                "log2_ordered_relations": math.log2(l24_cube),
                "orbit_rows": orbit_rows(
                    representatives24, l24_cube_partners
                ),
            },
        },
        "timing": {
            "classification_seconds": classified_seconds,
            "total_seconds": elapsed,
        },
        "validation": {
            "committed_shells_sorted_and_unique": True,
            "affine_orbits_partition_both_shells": True,
            "orbit_histograms_match_committed_enumerations": True,
            "L22_cube_reproduces_goal11": True,
            "mixed_relation_count_matches_both_shell_projections": True,
            "L24_partner_counts_are_even": True,
            "L24_cube_count_is_divisible_by_six": True,
        },
        "scope": (
            "Exact additive relation counts. This receipt does not yet count "
            "which field ratios match the shifted outer schedule."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
