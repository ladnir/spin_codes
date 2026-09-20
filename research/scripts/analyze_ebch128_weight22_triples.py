#!/usr/bin/env python3
"""Exactly count additive triples in the EBCH weight-22 shell.

Let L22 contain the 64-bit messages whose extended BCH encodings have weight
22.  The script counts ordered pairs (x,y) in L22^2 for which x+y also lies
in L22.  The affine group on 128 coordinates preserves L22.  The committed
shell has 15 free affine orbits, so it suffices to test one representative
from each orbit against the complete shell.
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as bch  # noqa: E402
from enumerate_ebch128_weight22_affine import (  # noqa: E402
    AFFINE_GROUP_SIZE,
    EXPECTED_COUNT,
    EXPECTED_WEIGHT,
    affine_orbit,
    codeword_points,
    coordinate_maps,
)


DEFAULT_MESSAGES = (
    SCRIPT_DIRECTORY.parent
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "ebch128_weight22_messages.bin"
)


def read_messages(path: Path) -> list[int]:
    data = path.read_bytes()
    if len(data) != 8 * EXPECTED_COUNT:
        raise RuntimeError("weight-22 message file has the wrong size")
    messages = [
        record[0] for record in struct.iter_unpack("<Q", data)
    ]
    if any(left >= right for left, right in zip(messages, messages[1:])):
        raise RuntimeError("weight-22 messages are not sorted and unique")
    return messages


def transform_codeword(
    points: list[int],
    products: list[int],
    translation: int,
    point_to_coordinate: list[int],
) -> int:
    image = 0
    for point in points:
        image |= 1 << point_to_coordinate[products[point] ^ translation]
    return image


def transformed_messages(
    codeword: int,
    coordinate_to_point: list[int],
    point_to_coordinate: list[int],
    multiplication: list[list[int]],
    message_by_codeword: dict[int, int],
) -> list[int]:
    points = codeword_points(codeword, coordinate_to_point)
    result = []
    for products in multiplication:
        for translation in range(128):
            image = transform_codeword(
                points, products, translation, point_to_coordinate
            )
            result.append(message_by_codeword[image])
    if len(result) != AFFINE_GROUP_SIZE or len(set(result)) != AFFINE_GROUP_SIZE:
        raise RuntimeError("ordered affine images do not form one free orbit")
    return result


def batch_inverse(values: list[int]) -> list[int]:
    prefixes = []
    product = 1
    for value in values:
        prefixes.append(product)
        product = bch.field_multiply(product, value)
    inverse_product = bch.field_inverse(product)
    result = [0] * len(values)
    for index in range(len(values) - 1, -1, -1):
        result[index] = bch.field_multiply(inverse_product, prefixes[index])
        inverse_product = bch.field_multiply(inverse_product, values[index])
    if any(
        bch.field_multiply(value, inverse) != 1
        for value, inverse in zip(values[::1024], result[::1024])
    ):
        raise RuntimeError("batch field inversion failed probe validation")
    return result


def ratio_histogram(
    representatives: list[int],
    partners_by_representative: list[list[int]],
    coordinate_to_point: list[int],
    point_to_coordinate: list[int],
    multiplication: list[list[int]],
    message_by_codeword: dict[int, int],
) -> dict[int, int]:
    histogram: dict[int, int] = {}
    for representative, partners in zip(
        representatives, partners_by_representative
    ):
        first_messages = transformed_messages(
            representative,
            coordinate_to_point,
            point_to_coordinate,
            multiplication,
            message_by_codeword,
        )
        first_inverses = batch_inverse(first_messages)
        for partner in partners:
            second_messages = transformed_messages(
                partner,
                coordinate_to_point,
                point_to_coordinate,
                multiplication,
                message_by_codeword,
            )
            for inverse, second in zip(first_inverses, second_messages):
                ratio = bch.field_multiply(second, inverse)
                histogram[ratio] = histogram.get(ratio, 0) + 1
    return histogram


def write_ratio_histogram(path: Path, histogram: dict[int, int]) -> None:
    records = sorted(histogram.items())
    payload = bytearray(16 * len(records))
    for index, (ratio, count) in enumerate(records):
        struct.pack_into("<QQ", payload, 16 * index, ratio, count)
    path.write_bytes(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--messages", type=Path, default=DEFAULT_MESSAGES)
    parser.add_argument("--ratio-output", type=Path)
    args = parser.parse_args()

    messages = read_messages(args.messages)
    message_by_codeword = {}
    for message in messages:
        codeword = bch.encode_message(message)
        if codeword.bit_count() != EXPECTED_WEIGHT:
            raise RuntimeError("message file contains a non-weight-22 codeword")
        message_by_codeword[codeword] = message
    if len(message_by_codeword) != EXPECTED_COUNT:
        raise RuntimeError("BCH encoding is not injective on the message set")

    coordinate_to_point, point_to_coordinate, multiplication = coordinate_maps()
    complete = set(message_by_codeword)
    uncovered = set(complete)
    representatives = []
    orbit_sizes = []
    while uncovered:
        representative = min(uncovered)
        orbit = affine_orbit(
            representative,
            coordinate_to_point,
            point_to_coordinate,
            multiplication,
        )
        if len(orbit) != AFFINE_GROUP_SIZE:
            raise RuntimeError("minimum-word affine orbit has the wrong size")
        if not orbit.issubset(complete):
            raise RuntimeError("affine orbit left the complete minimum shell")
        representatives.append(representative)
        orbit_sizes.append(len(orbit))
        uncovered.difference_update(orbit)
    if len(representatives) != 15:
        raise RuntimeError("minimum shell does not decompose into 15 affine orbits")

    orbit_rows = []
    total_ordered_triples = 0
    witnesses = []
    partners_by_representative = []
    for representative, orbit_size in zip(representatives, orbit_sizes):
        partners = []
        for other in complete:
            third = representative ^ other
            if third in complete:
                partners.append(other)
                if len(witnesses) < 8:
                    witnesses.append(
                        {
                            "first_message_hex": hex(
                                message_by_codeword[representative]
                            ),
                            "second_message_hex": hex(message_by_codeword[other]),
                            "third_message_hex": hex(message_by_codeword[third]),
                        }
                    )
        partners_by_representative.append(partners)
        total_ordered_triples += orbit_size * len(partners)
        orbit_rows.append(
            {
                "representative_message_hex": hex(
                    message_by_codeword[representative]
                ),
                "ordered_partners": len(partners),
            }
        )

    if total_ordered_triples % 6:
        raise RuntimeError("ordered additive triple count is not divisible by six")
    ratios = ratio_histogram(
        representatives,
        partners_by_representative,
        coordinate_to_point,
        point_to_coordinate,
        multiplication,
        message_by_codeword,
    )
    if sum(ratios.values()) != total_ordered_triples:
        raise RuntimeError("ratio histogram has the wrong total mass")
    if args.ratio_output:
        write_ratio_histogram(args.ratio_output, ratios)
    top_ratios = sorted(ratios.items(), key=lambda item: item[1], reverse=True)[:12]
    payload = {
        "schema": "ebch128-weight22-additive-triples-v1",
        "code": "extended primitive binary BCH [128,64,22]",
        "minimum_shell_size": len(complete),
        "affine_orbits": len(representatives),
        "affine_orbit_size": AFFINE_GROUP_SIZE,
        "representative_partner_counts": orbit_rows,
        "ordered_additive_triples": total_ordered_triples,
        "unordered_additive_triples": total_ordered_triples // 6,
        "ratio_histogram": {
            "distinct_ratios": len(ratios),
            "maximum_ordered_pairs_at_one_ratio": max(ratios.values()),
            "top_ratios": [
                {"ratio_hex": hex(ratio), "ordered_pairs": count}
                for ratio, count in top_ratios
            ],
            "binary_output": str(args.ratio_output) if args.ratio_output else None,
            "binary_record": "little-endian uint64 ratio, uint64 ordered-pair count",
        },
        "witnesses": witnesses,
        "validation": {
            "all_messages_reencoded_at_weight_22": True,
            "affine_orbits_partition_complete_shell": True,
            "representative_reduction_is_exact": True,
            "ordered_count_divisible_by_six": True,
            "ratio_histogram_total_matches_ordered_count": True,
        },
        "scope": (
            "Exact for triples of nonzero BCH messages whose three encodings "
            "all have weight 22 and whose field-message XOR is zero."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
