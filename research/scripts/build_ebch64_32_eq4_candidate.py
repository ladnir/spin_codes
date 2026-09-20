#!/usr/bin/env python3
"""Build an explicit complement-symmetric EBCH [64,32,>=12] subcode.

The parent is the parity extension of primitive narrow-sense BCH(63,11), with
parameters [64,36,>=12].  Intersect it with x_0+x_j=0 for j=1,2,3,4.  The four
checks are independent modulo the parent checks and retain the all-one word.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from check_bch_transposed_inner import extended_bch_check_rows, row_rank
from exact_bch_spectrum_small import nullspace_basis


LENGTH = 64
DIMENSION = 32
M = 6
DESIGNED_DISTANCE = 11
EXTRA_COORDINATES = (1, 2, 3, 4)


def rows_sha256(rows: list[int]) -> str:
    return hashlib.sha256(
        b"".join(row.to_bytes(LENGTH // 8, "little") for row in rows)
    ).hexdigest()


def build() -> tuple[dict[str, object], list[int]]:
    parent, length, dimension, rank = extended_bch_check_rows(
        M, DESIGNED_DISTANCE
    )
    if (length, dimension, rank) != (64, 36, 28):
        raise RuntimeError(
            f"unexpected extended BCH parent {(length, dimension, rank)}"
        )
    extra = [
        (1 << 0) | (1 << coordinate) for coordinate in EXTRA_COORDINATES
    ]
    ladder = [row_rank(parent + extra[:count]) for count in range(5)]
    if ladder != [28, 29, 30, 31, 32]:
        raise RuntimeError(f"extra checks are not independent: {ladder}")
    checks = parent + extra
    basis = list(nullspace_basis(checks, LENGTH))
    if len(basis) != DIMENSION or row_rank(basis) != DIMENSION:
        raise RuntimeError("generator basis has the wrong rank")
    all_one = (1 << LENGTH) - 1
    if any((check & all_one).bit_count() & 1 for check in checks):
        raise RuntimeError("subcode does not contain the all-one word")
    payload: dict[str, object] = {
        "schema": "ebch64-32-eq4-candidate-v1",
        "construction": (
            "parity-extended primitive BCH(63,11), intersected with four "
            "coordinate-equality checks"
        ),
        "parent": {
            "length": length,
            "dimension": dimension,
            "designed_distance": DESIGNED_DISTANCE,
            "distance_floor_after_extension": 12,
            "check_rows_sha256": rows_sha256(parent),
        },
        "subcode": {
            "length": LENGTH,
            "dimension": DIMENSION,
            "distance_floor": 12,
            "extra_checks": [
                f"x_0+x_{coordinate}=0"
                for coordinate in EXTRA_COORDINATES
            ],
            "rank_ladder": ladder,
            "contains_all_one": True,
            "complement_symmetric_spectrum": True,
            "check_rows_sha256": rows_sha256(checks),
            "generator_basis_sha256": rows_sha256(basis),
            "generator_rows_hex": [f"{row:016x}" for row in basis],
        },
    }
    return payload, basis


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--rows", type=Path, required=True)
    args = parser.parse_args()
    payload, basis = build()
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.rows.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.rows.write_text(
        "".join(f"{row:016x}\n" for row in basis), encoding="ascii"
    )
    print(f"manifest,{args.manifest}")
    print(f"rows,{args.rows}")
    print(f"generator_sha256,{payload['subcode']['generator_basis_sha256']}")


if __name__ == "__main__":
    main()
