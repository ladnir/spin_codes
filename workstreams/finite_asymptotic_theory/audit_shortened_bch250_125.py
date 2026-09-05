#!/usr/bin/env python3
"""Audit explicit shortenings of the extended BCH(255, 37) parent.

The parent is the parity extension of the primitive narrow-sense binary
BCH(255, 37) code.  This script adds coordinate-zero checks before puncturing
those coordinates.  It also writes a code-independent constant-weight packing
envelope for the explicit [250, 125, >=38] shortening.

The rank computations and integer packing bounds are exact.  Decimal logarithms
in the JSON receipt are rounded binary64 display values.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from check_bch_transposed_inner import extended_bch_check_rows, row_rank  # noqa: E402
from exact_bch_spectrum_small import nullspace_basis  # noqa: E402


PARENT_LENGTH = 256
PARENT_DIMENSION = 131
PARENT_DISTANCE_FLOOR = 38
M = 8
DESIGNED_DISTANCE = 37


def rows_sha256(rows: list[int], length: int) -> str:
    byte_length = (length + 7) // 8
    return hashlib.sha256(
        b"".join(row.to_bytes(byte_length, "little") for row in rows)
    ).hexdigest()


def shortened_code(check_rows: list[int], count: int) -> dict[str, object]:
    coordinate_checks = [1 << coordinate for coordinate in range(count)]
    rank_ladder = [
        row_rank(check_rows + coordinate_checks[:prefix])
        for prefix in range(count + 1)
    ]
    constrained = check_rows + coordinate_checks
    parent_basis = list(nullspace_basis(constrained, PARENT_LENGTH))
    shortened_rows = [row >> count for row in parent_basis]
    length = PARENT_LENGTH - count
    dimension = len(shortened_rows)
    shortened_rank = row_rank(shortened_rows)
    if shortened_rank != dimension:
        raise RuntimeError("puncturing a zero coordinate changed the dimension")
    return {
        "shortened_coordinates": list(range(count)),
        "parameters": {
            "length": length,
            "dimension": dimension,
            "minimum_distance_floor": PARENT_DISTANCE_FLOOR,
        },
        "check_rank_ladder": rank_ladder,
        "generator_rank_after_puncture": shortened_rank,
        "generator_rows_sha256": rows_sha256(shortened_rows, length),
    }


def constant_weight_packing_upper(length: int, weight: int, distance: int) -> int:
    """Return the Johnson-space sphere-packing upper bound A(n,d,w)."""

    johnson_distance = (distance + 1) // 2
    packing_radius = (johnson_distance - 1) // 2
    ball = sum(
        math.comb(weight, index) * math.comb(length - weight, index)
        for index in range(packing_radius + 1)
        if index <= weight and index <= length - weight
    )
    return math.comb(length, weight) // ball


def spectrum_envelope(
    length: int,
    dimension: int,
    distance: int,
    maximum_weight: int,
) -> dict[str, object]:
    rows: list[dict[str, int | float]] = [
        {"weight": 0, "multiplicity_upper": 1, "log2_expected_multiplicity": 0.0}
    ]
    total_nonzero = (1 << dimension) - 1
    for weight in range(distance, maximum_weight + 1):
        if weight & 1:
            continue
        upper = min(
            total_nonzero,
            constant_weight_packing_upper(length, weight, distance),
        )
        if upper == 0:
            continue
        rows.append(
            {
                "weight": weight,
                "multiplicity_upper": upper,
                # The existing analyzer consumes this historical field name.
                "log2_expected_multiplicity": math.log2(upper),
            }
        )
    return {
        "schema": "shortened-bch-constant-weight-envelope-v1",
        "candidate": "ShortenedExtendedBCH250x125",
        "parameters": {
            "outer_bits": length,
            "outer_dimension": dimension,
            "minimum_distance_floor": distance,
            "maximum_weight": maximum_weight,
            "even_weight": True,
            "johnson_packing_radius": 9,
        },
        "spectrum": rows,
        "scope": (
            "A code-independent per-shell upper envelope. Each nonzero shell "
            "uses the minimum of 2^125-1 and the exact integer Johnson-space "
            "sphere-packing bound for binary constant-weight codes of Hamming "
            "distance at least 38. The upper cutoff 218 follows by complementing "
            "the zero-extended shortened word inside the parent, which contains "
            "the all-one word and has distance at least 38. The logarithms are "
            "binary64 display values; "
            "the multiplicity_upper integers are the exact bounds. This is not "
            "the weight distribution of the shortened BCH code."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--envelope-output", type=Path, required=True)
    args = parser.parse_args()

    parent_checks, length, dimension, check_rank = extended_bch_check_rows(
        M, DESIGNED_DISTANCE
    )
    if (length, dimension, check_rank) != (
        PARENT_LENGTH,
        PARENT_DIMENSION,
        PARENT_LENGTH - PARENT_DIMENSION,
    ):
        raise RuntimeError(
            f"unexpected extended BCH parent {(length, dimension, check_rank)}"
        )
    all_one = (1 << PARENT_LENGTH) - 1
    if any((check & all_one).bit_count() & 1 for check in parent_checks):
        raise RuntimeError("extended BCH parent unexpectedly omits the all-one word")

    shortened3 = shortened_code(parent_checks, 3)
    shortened6 = shortened_code(parent_checks, 6)
    if shortened3["parameters"] != {
        "length": 253,
        "dimension": 128,
        "minimum_distance_floor": 38,
    }:
        raise RuntimeError(f"unexpected three-coordinate shortening: {shortened3}")
    if shortened6["parameters"] != {
        "length": 250,
        "dimension": 125,
        "minimum_distance_floor": 38,
    }:
        raise RuntimeError(f"unexpected six-coordinate shortening: {shortened6}")

    audit = {
        "schema": "shortened-extended-bch256-audit-v1",
        "parent": {
            "construction": (
                "parity extension of primitive narrow-sense binary BCH(255,37)"
            ),
            "parameters": {
                "length": length,
                "dimension": dimension,
                "minimum_distance_floor": PARENT_DISTANCE_FLOOR,
            },
            "parity_check_rank": check_rank,
            "contains_all_one_word": True,
        },
        "shortenings": [shortened3, shortened6],
        "scope": (
            "The GF(2) rank computations, shortening dimensions, and generator "
            "hashes are exact. The distance floor is inherited from the standard "
            "BCH bound and parity extension; this script does not enumerate "
            "minimum-weight codewords. Shortening coordinates are fixed and are "
            "not sampled."
        ),
    }
    envelope = spectrum_envelope(
        250,
        125,
        PARENT_DISTANCE_FLOOR,
        PARENT_LENGTH - PARENT_DISTANCE_FLOOR,
    )

    args.audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.envelope_output.parent.mkdir(parents=True, exist_ok=True)
    args.audit_output.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.envelope_output.write_text(
        json.dumps(envelope, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.audit_output}")
    print(f"wrote,{args.envelope_output}")
    print("shortened_253_128_rank,128")
    print("shortened_250_125_rank,125")


if __name__ == "__main__":
    main()
