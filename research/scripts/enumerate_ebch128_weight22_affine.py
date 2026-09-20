#!/usr/bin/env python3
"""Generate all EBCH [128,64,22] minimum words from affine orbits.

The extended primitive BCH code is invariant under x -> a*x+b on GF(128).
The committed spectrum has A_22 = 243840 = 15*127*128.  This script reads
authenticated low-packet-support messages, uses their weight-22 words as
candidate orbit representatives, and verifies whether their affine orbits
cover the complete committed minimum-weight spectrum.

An optional binary output stores the resulting 64-bit messages in sorted
order.  The script prints JSON.
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as dp  # noqa: E402
from bch_affine_orbit_sanity import primitive_position_values  # noqa: E402
from check_bch_boundary_smallfield import find_primitive_poly, gf_mul  # noqa: E402


FIELD_BITS = 7
FIELD_SIZE = 1 << FIELD_BITS
PRIMITIVE_LENGTH = FIELD_SIZE - 1
CODE_BITS = 128
EXPECTED_WEIGHT = 22
EXPECTED_COUNT = 243_840
AFFINE_GROUP_SIZE = PRIMITIVE_LENGTH * FIELD_SIZE


def read_messages(path: Path) -> list[int]:
    data = path.read_bytes()
    if len(data) % 9:
        raise ValueError("authenticated low-message file has a partial record")
    messages = []
    for offset in range(0, len(data), 9):
        _support, message = struct.unpack_from("<BQ", data, offset)
        if dp.encode_message(message).bit_count() == EXPECTED_WEIGHT:
            messages.append(message)
    if len(messages) != len(set(messages)):
        raise ValueError("authenticated low-message file repeats a weight-22 message")
    return messages


def generator_inverse() -> int:
    generator = dp.GENERATOR_ROWS[0] & dp.FIELD_ORDER
    inverse = 1
    for bit in range(1, 64):
        product = carryless_multiply_low(generator, inverse)
        if (product >> bit) & 1:
            inverse |= 1 << bit
    if carryless_multiply_low(generator, inverse) != 1:
        raise RuntimeError("failed to invert the systematic BCH generator half")
    return inverse


def carryless_multiply_low(left: int, right: int) -> int:
    result = 0
    while right:
        bit = (right & -right).bit_length() - 1
        result ^= (left << bit) & dp.FIELD_ORDER
        right &= right - 1
    return result


def coordinate_maps() -> tuple[list[int], list[int], list[list[int]]]:
    polynomial = find_primitive_poly(FIELD_BITS)
    values = primitive_position_values(FIELD_BITS, polynomial)
    point_to_coordinate = [-1] * FIELD_SIZE
    point_to_coordinate[0] = 127
    for coordinate, point in enumerate(values):
        point_to_coordinate[point] = coordinate
    if any(coordinate < 0 for coordinate in point_to_coordinate):
        raise RuntimeError("primitive coordinate labels do not cover GF(128)")
    multiplication = [
        [gf_mul(a, point, polynomial, FIELD_BITS) for point in range(FIELD_SIZE)]
        for a in range(1, FIELD_SIZE)
    ]
    return values, point_to_coordinate, multiplication


def codeword_points(codeword: int, coordinate_to_point: list[int]) -> list[int]:
    points = []
    primitive = codeword & ((1 << 127) - 1)
    while primitive:
        bit = primitive & -primitive
        coordinate = bit.bit_length() - 1
        points.append(coordinate_to_point[coordinate])
        primitive ^= bit
    if (codeword >> 127) & 1:
        points.append(0)
    if len(points) != EXPECTED_WEIGHT:
        raise RuntimeError("orbit representative does not have weight 22")
    return points


def affine_orbit(
    codeword: int,
    coordinate_to_point: list[int],
    point_to_coordinate: list[int],
    multiplication: list[list[int]],
) -> set[int]:
    points = codeword_points(codeword, coordinate_to_point)
    orbit: set[int] = set()
    for products in multiplication:
        for translation in range(FIELD_SIZE):
            image = 0
            for product in (products[point] for point in points):
                image |= 1 << point_to_coordinate[product ^ translation]
            orbit.add(image)
    return orbit


def decode_message(codeword: int, inverse: int) -> int:
    message = carryless_multiply_low(codeword & dp.FIELD_ORDER, inverse)
    if dp.encode_message(message) != codeword:
        raise RuntimeError("affine image failed BCH re-encoding")
    return message


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--low-message-file", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    known_messages = read_messages(args.low_message_file)
    known_codewords = {dp.encode_message(message) for message in known_messages}
    coordinate_to_point, point_to_coordinate, multiplication = coordinate_maps()
    inverse = generator_inverse()

    uncovered = set(known_codewords)
    complete_codewords: set[int] = set()
    orbit_sizes = []
    representatives = []
    while uncovered:
        representative = min(uncovered)
        orbit = affine_orbit(
            representative,
            coordinate_to_point,
            point_to_coordinate,
            multiplication,
        )
        if len(orbit) != AFFINE_GROUP_SIZE:
            raise RuntimeError("weight-22 affine orbit has a nontrivial stabilizer")
        overlap = complete_codewords.intersection(orbit)
        if overlap:
            raise RuntimeError("new affine orbit overlaps a prior orbit")
        complete_codewords.update(orbit)
        uncovered.difference_update(orbit)
        representatives.append(representative)
        orbit_sizes.append(len(orbit))

    messages = sorted(decode_message(codeword, inverse) for codeword in complete_codewords)
    if len(messages) != EXPECTED_COUNT:
        raise RuntimeError(
            "authenticated representatives do not cover the committed A_22 spectrum"
        )
    if not set(known_messages).issubset(messages):
        raise RuntimeError("generated minimum-weight set lost an authenticated message")

    if args.output:
        with args.output.open("wb") as destination:
            for message in messages:
                destination.write(struct.pack("<Q", message))

    payload = {
        "schema": "ebch128-weight22-affine-enumeration-v1",
        "evidence_label": "EXACT_AFFINE_ORBIT_ENUMERATION_WITH_REENCODING",
        "parameters": {
            "code": "extended primitive binary BCH [128,64,22]",
            "affine_group_size": AFFINE_GROUP_SIZE,
            "committed_A_22": EXPECTED_COUNT,
        },
        "authenticated_seed": {
            "file": str(args.low_message_file),
            "weight22_messages": len(known_messages),
        },
        "enumeration": {
            "affine_orbits": len(orbit_sizes),
            "orbit_sizes": orbit_sizes,
            "generated_messages": len(messages),
            "all_images_reencoded": True,
            "matches_committed_A_22": len(messages) == EXPECTED_COUNT,
        },
        "output": str(args.output) if args.output else None,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
