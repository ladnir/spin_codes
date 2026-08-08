#!/usr/bin/env python3
"""Build and authenticate the exact BCH-like [512,256,>=62] candidate.

The parent is the parity extension of the primitive narrow-sense binary BCH
code B(511,61), which has dimension 259 by its cyclotomic-coset degree count.
We intersect the extended parent with the three deterministic even-weight
checks x_0+x_j=0 for j=1,2,3.  These checks are independent modulo the parent
checks, retain the all-one word, and produce a complement-symmetric
[512,256,>=62] subcode.

The script records hashes of the complete parity-check matrix and a
deterministically generated nullspace basis.  The matrices themselves are
reproducible from this source and the recorded primitive polynomial.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from bch_candidate_params import bch_dimension
from check_bch_boundary_smallfield import find_primitive_poly
from check_bch_transposed_inner import extended_bch_check_rows, row_rank
from exact_bch_spectrum_small import nullspace_basis


ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "bch512_256_candidate.json"
M = 9
DESIGNED_DISTANCE = 61
EXTRA_CHECK_COORDINATES = (1, 2, 3)


def rows_sha256(rows: list[int], width: int) -> str:
    size = (width + 7) // 8
    payload = b"".join(row.to_bytes(size, "little") for row in rows)
    return hashlib.sha256(payload).hexdigest()


def canonical_sha256(payload: dict[str, object]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def build_payload() -> dict[str, object]:
    primitive_n, parent_dimension, parent_redundancy = bch_dimension(M, DESIGNED_DISTANCE)
    primitive_poly = find_primitive_poly(M)
    parent_rows, extended_n, checked_dimension, parent_rank = extended_bch_check_rows(
        M, DESIGNED_DISTANCE
    )
    if (primitive_n, parent_dimension, parent_redundancy) != (511, 259, 252):
        raise SystemExit("unexpected BCH(511,61) cyclotomic-coset parameters")
    if (extended_n, checked_dimension, parent_rank) != (512, 259, 253):
        raise SystemExit("unexpected extended BCH parent parameters")

    extra_rows = [(1 << 0) | (1 << coordinate) for coordinate in EXTRA_CHECK_COORDINATES]
    ranks = [row_rank(parent_rows + extra_rows[:count]) for count in range(4)]
    if ranks != [253, 254, 255, 256]:
        raise SystemExit(f"extra checks are not independently rank-raising: {ranks}")
    subcode_rows = parent_rows + extra_rows
    basis = nullspace_basis(subcode_rows, extended_n)
    if len(basis) != 256 or row_rank(basis) != 256:
        raise SystemExit("candidate generator basis does not have rank 256")
    if any((row & word).bit_count() & 1 for row in subcode_rows for word in basis):
        raise SystemExit("candidate generator basis is not orthogonal to its checks")

    all_one = (1 << extended_n) - 1
    if any((row & all_one).bit_count() & 1 for row in subcode_rows):
        raise SystemExit("candidate does not retain the all-one word")
    if not all(word.bit_count() % 2 == 0 for word in basis):
        raise SystemExit("candidate basis contains an odd-weight word")

    payload: dict[str, object] = {
        "schema": 1,
        "construction": "extended primitive narrow-sense BCH parent plus three even checks",
        "field": {
            "m": M,
            "primitive_polynomial_hex": hex(primitive_poly),
            "coordinate_alpha_hex": "0x2",
        },
        "parent": {
            "primitive_length": primitive_n,
            "extended_length": extended_n,
            "dimension": parent_dimension,
            "primitive_redundancy": parent_redundancy,
            "extended_check_rank": parent_rank,
            "designed_distance": DESIGNED_DISTANCE,
            "extended_distance_floor": DESIGNED_DISTANCE + 1,
            "distance_reason": (
                "BCH bound gives primitive distance >=61; parity extension is even, "
                "so its nonzero weights are >=62"
            ),
            "check_rows_sha256": rows_sha256(parent_rows, extended_n),
        },
        "subcode": {
            "length": extended_n,
            "dimension": len(basis),
            "distance_floor": DESIGNED_DISTANCE + 1,
            "extra_checks": [f"x_0+x_{coordinate}=0" for coordinate in EXTRA_CHECK_COORDINATES],
            "rank_ladder": ranks,
            "contains_all_one": True,
            "even_weight": True,
            "complement_symmetric_spectrum": True,
            "check_rows_sha256": rows_sha256(subcode_rows, extended_n),
            "generator_basis_sha256": rows_sha256(basis, extended_n),
            "basis_min_row_weight": min(word.bit_count() for word in basis),
            "basis_max_row_weight": max(word.bit_count() for word in basis),
        },
    }
    payload["sha256"] = canonical_sha256(payload)
    return payload


def write_manifest(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n")


def verify_manifest(path: Path, regenerated: dict[str, object]) -> None:
    stored = json.loads(path.read_text())
    claimed = stored.pop("sha256", None)
    actual = canonical_sha256(stored)
    if claimed != actual:
        raise SystemExit("BCH candidate manifest SHA-256 mismatch")
    stored["sha256"] = claimed
    if stored != regenerated:
        raise SystemExit("BCH candidate manifest differs from regenerated construction")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()

    payload = build_payload()
    if args.write_manifest:
        write_manifest(args.manifest, payload)
        print(f"bch512_256_candidate_manifest_written,{args.manifest}")
    else:
        verify_manifest(args.manifest, payload)
        print("bch512_256_candidate_status,PASS")
    print(f"candidate_sha256,{payload['sha256']}")
    print(f"parent_check_sha256,{payload['parent']['check_rows_sha256']}")
    print(f"subcode_check_sha256,{payload['subcode']['check_rows_sha256']}")
    print(f"generator_basis_sha256,{payload['subcode']['generator_basis_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
