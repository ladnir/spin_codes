#!/usr/bin/env python3
"""Enumerate the adjacent EBCH occupation-three profiles exactly.

The two profiles are ``(22,22,24)`` and ``(22,106,106)``.  Let ``L22`` be
the complete minimum shell and let ``e`` encode to the all-one word.

For the first profile, enumerate ordered pairs ``x,y`` in ``L22`` for which
``B(x+y)`` has weight 24.  For the second profile, reuse the relations
``x+y+z=0`` in ``L22`` and complement the last two words:

    B(x), B(e+y), B(e+z).

The affine coordinate group preserves every BCH weight and fixes the all-one
word.  The minimum shell is a union of 15 free affine orbits.  It therefore
suffices to test one first-word representative from each orbit.  Optional
binary outputs store exact field-ratio histograms for later outer-support
matching.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as bch  # noqa: E402
from analyze_ebch128_weight22_triples import (  # noqa: E402
    batch_inverse,
    read_messages,
    transformed_messages,
    write_ratio_histogram,
)
from analyze_riffle_bch_alpha_joint_spectrum import encode_weights  # noqa: E402
from analyze_riffle_bch_alpha_joint_spectrum import field_multiply_batch  # noqa: E402
from enumerate_ebch128_weight22_affine import (  # noqa: E402
    AFFINE_GROUP_SIZE,
    EXPECTED_COUNT,
    affine_orbit,
    carryless_multiply_low,
    coordinate_maps,
    generator_inverse,
)


DEFAULT_MESSAGES = (
    SCRIPT_DIRECTORY.parent
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "ebch128_weight22_messages.bin"
)
ALL_ONE_MESSAGE = 0x8E63CF44EFD4FA21


def classify_minimum_shell(
    messages: list[int],
) -> tuple[
    list[int],
    dict[int, int],
    list[int],
    list[int],
    list[list[int]],
]:
    """Return orbit representatives and the authenticated coordinate maps."""
    message_by_codeword: dict[int, int] = {}
    for message in messages:
        codeword = bch.encode_message(message)
        if codeword.bit_count() != 22:
            raise RuntimeError("minimum-shell file contains a non-weight-22 word")
        message_by_codeword[codeword] = message
    if len(message_by_codeword) != EXPECTED_COUNT:
        raise RuntimeError("BCH encoding is not injective on the minimum shell")

    coordinate_to_point, point_to_coordinate, multiplication = coordinate_maps()
    complete = set(message_by_codeword)
    uncovered = set(complete)
    representatives: list[int] = []
    while uncovered:
        representative = min(uncovered)
        orbit = affine_orbit(
            representative,
            coordinate_to_point,
            point_to_coordinate,
            multiplication,
        )
        if len(orbit) != AFFINE_GROUP_SIZE or not orbit.issubset(complete):
            raise RuntimeError("minimum-shell affine orbit validation failed")
        representatives.append(representative)
        uncovered.difference_update(orbit)
    if len(representatives) != 15:
        raise RuntimeError("minimum shell does not have 15 affine orbits")
    return (
        representatives,
        message_by_codeword,
        coordinate_to_point,
        point_to_coordinate,
        multiplication,
    )


def partner_lists(
    representatives: list[int],
    message_by_codeword: dict[int, int],
    messages_array: np.ndarray,
) -> tuple[list[list[int]], list[list[int]], list[dict[str, object]]]:
    """Find representative partners for both adjacent profiles."""
    minimum_codewords = set(message_by_codeword)
    additive_partners: list[list[int]] = []
    weight24_partners: list[list[int]] = []
    rows: list[dict[str, object]] = []
    for representative in representatives:
        representative_message = message_by_codeword[representative]
        minimum = [
            message_by_codeword[other]
            for other in minimum_codewords
            if representative ^ other in minimum_codewords
        ]
        target_weights = encode_weights(
            messages_array ^ np.uint64(representative_message)
        )
        weight24 = [
            int(value) for value in messages_array[target_weights == 24]
        ]
        additive_partners.append(minimum)
        weight24_partners.append(weight24)
        rows.append(
            {
                "representative_message_hex": hex(representative_message),
                "weight22_xor_weight22_partners": len(minimum),
                "weight24_xor_partners": len(weight24),
            }
        )
    return additive_partners, weight24_partners, rows


def ratio_histogram(
    *,
    representatives: list[int],
    partners_by_representative: list[list[int]],
    message_by_codeword: dict[int, int],
    coordinate_to_point: list[int],
    point_to_coordinate: list[int],
    multiplication: list[list[int]],
    complement_second: bool,
) -> dict[int, int]:
    """Count second/first field ratios over every affine image of each pair."""
    histogram: dict[int, int] = {}
    for representative, partners in zip(
        representatives, partners_by_representative
    ):
        first = transformed_messages(
            representative,
            coordinate_to_point,
            point_to_coordinate,
            multiplication,
            message_by_codeword,
        )
        inverses = batch_inverse(first)
        for partner_message in partners:
            partner_codeword = bch.encode_message(partner_message)
            second = transformed_messages(
                partner_codeword,
                coordinate_to_point,
                point_to_coordinate,
                multiplication,
                message_by_codeword,
            )
            if complement_second:
                second = [value ^ ALL_ONE_MESSAGE for value in second]
            for inverse, value in zip(inverses, second):
                ratio = bch.field_multiply(value, inverse)
                histogram[ratio] = histogram.get(ratio, 0) + 1
    return histogram


def histogram_summary(
    histogram: dict[int, int], output: Path | None
) -> dict[str, object]:
    top = sorted(histogram.items(), key=lambda item: item[1], reverse=True)[:12]
    return {
        "total_mass": sum(histogram.values()),
        "distinct_ratios": len(histogram),
        "maximum_multiplicity": max(histogram.values()) if histogram else 0,
        "top_ratios": [
            {"ratio_hex": hex(ratio), "ordered_pairs": count}
            for ratio, count in top
        ],
        "binary_output": str(output) if output else None,
        "binary_record": "little-endian uint64 ratio, uint64 count",
    }


def affine_message_actions(
    coordinate_to_point: list[int],
    point_to_coordinate: list[int],
    multiplication: list[list[int]],
) -> np.ndarray:
    """Build the linear message action of every ordered affine map."""
    inverse = generator_inverse()
    decode_basis = np.zeros(128, dtype=np.uint64)
    for coordinate in range(64):
        decode_basis[coordinate] = np.uint64(
            carryless_multiply_low(1 << coordinate, inverse)
        )
    actions = np.empty((AFFINE_GROUP_SIZE, 128), dtype=np.uint64)
    source_points = coordinate_to_point + [0]
    row = 0
    for products in multiplication:
        for translation in range(128):
            destinations = [
                point_to_coordinate[products[point] ^ translation]
                for point in source_points
            ]
            actions[row] = decode_basis[np.asarray(destinations, dtype=np.uint8)]
            row += 1
    if row != AFFINE_GROUP_SIZE:
        raise RuntimeError("affine message-action table has the wrong height")
    return actions


def transformed_messages_fast(codeword: int, actions: np.ndarray) -> np.ndarray:
    coordinates = [
        coordinate for coordinate in range(128) if (codeword >> coordinate) & 1
    ]
    return np.bitwise_xor.reduce(actions[:, coordinates], axis=1)


def write_raw_ratio_stream_222224(
    path: Path,
    *,
    representatives: list[int],
    partners_by_representative: list[list[int]],
    message_by_codeword: dict[int, int],
    coordinate_to_point: list[int],
    point_to_coordinate: list[int],
    multiplication: list[list[int]],
) -> int:
    """Write one uint64 y/x record for every ordered (22,22,24) relation."""
    actions = affine_message_actions(
        coordinate_to_point, point_to_coordinate, multiplication
    )
    probe = transformed_messages_fast(representatives[0], actions)
    expected_probe = transformed_messages(
        representatives[0],
        coordinate_to_point,
        point_to_coordinate,
        multiplication,
        message_by_codeword,
    )
    if probe.tolist() != expected_probe:
        raise RuntimeError("vectorized affine message action failed exact audit")

    written = 0
    with path.open("wb") as destination:
        for representative, partners in zip(
            representatives, partners_by_representative
        ):
            first = transformed_messages_fast(representative, actions)
            inverses = np.asarray(batch_inverse(first.tolist()), dtype=np.uint64)
            for partner_message in partners:
                second = transformed_messages_fast(
                    bch.encode_message(partner_message), actions
                )
                ratios = field_multiply_batch(second, inverses)
                ratios.tofile(destination)
                written += int(ratios.size)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--messages", type=Path, default=DEFAULT_MESSAGES)
    parser.add_argument("--ratio-22106106", type=Path)
    parser.add_argument("--ratio-222224", type=Path)
    parser.add_argument(
        "--raw-ratio-222224",
        type=Path,
        help="write one little-endian uint64 ratio per ordered relation",
    )
    parser.add_argument(
        "--counts-only",
        action="store_true",
        help="stop after exact global relation counts; do not expand ratios",
    )
    args = parser.parse_args()

    messages = read_messages(args.messages)
    messages_array = np.asarray(messages, dtype=np.uint64)
    (
        representatives,
        message_by_codeword,
        coordinate_to_point,
        point_to_coordinate,
        multiplication,
    ) = classify_minimum_shell(messages)
    additive, weight24, orbit_rows = partner_lists(
        representatives, message_by_codeword, messages_array
    )

    ordered_22106106 = AFFINE_GROUP_SIZE * sum(map(len, additive))
    ordered_222224 = AFFINE_GROUP_SIZE * sum(map(len, weight24))
    if ordered_22106106 != 1_365_504:
        raise RuntimeError("minimum additive-triple count changed")
    if ordered_222224 % 2:
        raise RuntimeError("(22,22,24) ordered pair count is not even")

    ratios_22106106: dict[int, int] = {}
    ratios_222224: dict[int, int] = {}
    raw_ratio_count = None
    if args.raw_ratio_222224:
        raw_ratio_count = write_raw_ratio_stream_222224(
            args.raw_ratio_222224,
            representatives=representatives,
            partners_by_representative=weight24,
            message_by_codeword=message_by_codeword,
            coordinate_to_point=coordinate_to_point,
            point_to_coordinate=point_to_coordinate,
            multiplication=multiplication,
        )
        if raw_ratio_count != ordered_222224:
            raise RuntimeError("raw (22,22,24) ratio stream has the wrong size")
    elif not args.counts_only:
        ratios_22106106 = ratio_histogram(
            representatives=representatives,
            partners_by_representative=additive,
            message_by_codeword=message_by_codeword,
            coordinate_to_point=coordinate_to_point,
            point_to_coordinate=point_to_coordinate,
            multiplication=multiplication,
            complement_second=True,
        )
        ratios_222224 = ratio_histogram(
            representatives=representatives,
            partners_by_representative=weight24,
            message_by_codeword=message_by_codeword,
            coordinate_to_point=coordinate_to_point,
            point_to_coordinate=point_to_coordinate,
            multiplication=multiplication,
            complement_second=False,
        )
        if sum(ratios_22106106.values()) != ordered_22106106:
            raise RuntimeError("(22,106,106) ratio mass mismatch")
        if sum(ratios_222224.values()) != ordered_222224:
            raise RuntimeError("(22,22,24) ratio mass mismatch")
        if args.ratio_22106106:
            write_ratio_histogram(args.ratio_22106106, ratios_22106106)
        if args.ratio_222224:
            write_ratio_histogram(args.ratio_222224, ratios_222224)

    if bch.encode_message(ALL_ONE_MESSAGE) != (1 << 128) - 1:
        raise RuntimeError("committed all-one message changed")
    payload = {
        "schema": "ebch128-adjacent-additive-triples-v1",
        "code": "extended primitive binary BCH [128,64,22]",
        "minimum_shell_size": len(messages),
        "affine_orbits": len(representatives),
        "affine_orbit_size": AFFINE_GROUP_SIZE,
        "all_one_message_hex": hex(ALL_ONE_MESSAGE),
        "representative_rows": orbit_rows,
        "profiles": {
            "22_106_106": {
                "identity": "(x,e+y,e+z), where x+y+z=0 in L22",
                "ordered_relations": ordered_22106106,
                "log2_ordered_relations": math.log2(ordered_22106106),
                "ratio": "(e+y)/x, with the light coordinate first",
                "ratio_histogram": (
                    histogram_summary(ratios_22106106, args.ratio_22106106)
                    if not args.counts_only and not args.raw_ratio_222224
                    else None
                ),
            },
            "22_22_24": {
                "identity": "x,y in L22 and wt(B(x+y))=24",
                "ordered_relations": ordered_222224,
                "unordered_relations": ordered_222224 // 2,
                "log2_ordered_relations": (
                    math.log2(ordered_222224) if ordered_222224 else -math.inf
                ),
                "ratio": "y/x between the two light coordinates",
                "ratio_histogram": (
                    histogram_summary(ratios_222224, args.ratio_222224)
                    if not args.counts_only and not args.raw_ratio_222224
                    else None
                ),
                "raw_ratio_stream": (
                    {
                        "output": str(args.raw_ratio_222224),
                        "records": raw_ratio_count,
                        "binary_record": "little-endian uint64 ratio",
                    }
                    if args.raw_ratio_222224
                    else None
                ),
            },
        },
        "validation": {
            "all_messages_reencoded_at_weight_22": True,
            "affine_orbits_partition_complete_shell": True,
            "all_one_message_reencoded": True,
            "minimum_triple_count_matches_goal11": True,
            "weight24_targets_reencoded_vectorially": True,
            "ratio_histogram_masses_match_relation_counts": (
                True
                if not args.counts_only and not args.raw_ratio_222224
                else "NOT_RUN"
            ),
            "raw_ratio_stream_matches_relation_count": (
                True if args.raw_ratio_222224 else "NOT_RUN"
            ),
        },
        "scope": (
            "Exact global additive relation counts and field-ratio histograms. "
            "No outer support schedule has yet been applied."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
