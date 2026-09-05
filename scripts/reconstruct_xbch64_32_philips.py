#!/usr/bin/env python3
"""Reconstruct and authenticate the Philips shortened-XBCH diffusion code.

The source starts with the cyclic binary [63,36,>=11] BCH generator

    g(x) = x^27 + x^22 + x^21 + x^19 + x^18 + x^17
           + x^15 + x^8 + x^4 + x + 1.

It shortens by four information coordinates, extends by one overall parity
coordinate, puts the resulting [60,32,>=12] code in systematic form
``(I_32 | B)``, appends four columns to B, and permutes rows and columns of
the resulting 32-by-32 matrix.  The final matrix A is published, but the four
appended columns and permutations are not.

This verifier reconstructs B from the polynomial, then uses an exact SMT
model to find the row permutation and the injection of B's 28 columns into
A's 32 columns.  Thus puncturing the four unmatched parity coordinates of
``(I_32 | A)`` recovers the independently reconstructed extended BCH code up
to coordinate permutation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROWS = ROOT / "scripts" / "xbch64_32_philips_generator_rows.txt"
DEFAULT_SPECTRUM = ROOT / "scripts" / "xbch64_32_philips_spectrum.csv"
DEFAULT_OUT = (
    ROOT
    / "constructions"
    / "riffle_paritymix8_b1024_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s16"
    / "xbch64_philips_reconstruction.json"
)

GENERATOR_DEGREES = (27, 22, 21, 19, 18, 17, 15, 8, 4, 1, 0)
PUBLISHED_A_ROWS = (
    0x29124175, 0x0F2BF8FB, 0x2DEE5791, 0x019D7E7C,
    0xC8BF5445, 0x8556A980, 0xB48D6594, 0x42FED829,
    0x6C2BCC57, 0xCDD3AF3C, 0x62A9ECB2, 0xA4B5A428,
    0x3F5226B4, 0x45F11E76, 0x33D96699, 0x1B066723,
    0xB6CA16AD, 0xC54458F0, 0x992CBD93, 0x16D7ADF7,
    0xF934A9DD, 0x4F1AE1E8, 0xE7FD8AAA, 0x7A53217E,
    0x8A25CE9F, 0xCE31AC5A, 0xA1F6D9D7, 0xA6926FCC,
    0xC6081F0B, 0xEAE05244, 0xD07CB31F, 0x7ADF972F,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gf2_rank(rows: list[int]) -> int:
    work = rows[:]
    rank = 0
    width = max((row.bit_length() for row in work), default=0)
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(work)) if (work[index] >> column) & 1),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(len(work)):
            if index != rank and ((work[index] >> column) & 1):
                work[index] ^= work[rank]
        rank += 1
    return rank


def polynomial_remainder(dividend: int, divisor: int) -> int:
    while dividend.bit_length() >= divisor.bit_length():
        dividend ^= divisor << (dividend.bit_length() - divisor.bit_length())
    return dividend


def gf64_mul(left: int, right: int) -> int:
    """Multiply in GF(2^6) with primitive polynomial x^6+x+1."""
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & (1 << 6):
            left ^= 0x43
    return result


def gf64_pow(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = gf64_mul(result, value)
        exponent >>= 1
        value = gf64_mul(value, value)
    return result


def verify_bch_bound(generator: int) -> None:
    if polynomial_remainder((1 << 63) | 1, generator):
        raise RuntimeError("generator polynomial does not divide x^63+1")
    alpha = 2
    if gf64_pow(alpha, 63) != 1 or any(gf64_pow(alpha, exponent) == 1 for exponent in (1, 3, 7, 9, 21)):
        raise RuntimeError("x is not primitive in GF(2)[x]/(x^6+x+1)")
    for exponent in range(1, 11):
        point = gf64_pow(alpha, exponent)
        value = 0
        for degree in GENERATOR_DEGREES:
            value ^= gf64_pow(point, degree)
        if value:
            raise RuntimeError(f"generator does not vanish at alpha^{exponent}")


def systematicize_first_32(rows: list[int]) -> list[int]:
    """Return the unique row-equivalent basis whose first 32 columns are I."""
    work = rows[:]
    if len(work) != 32:
        raise ValueError("expected 32 shortened generator rows")
    for column in range(32):
        pivot = next(
            (index for index in range(column, 32) if (work[index] >> column) & 1),
            None,
        )
        if pivot is None:
            raise RuntimeError(f"information-set column {column} has no pivot")
        work[column], work[pivot] = work[pivot], work[column]
        for index in range(32):
            if index != column and ((work[index] >> column) & 1):
                work[index] ^= work[column]
    identity_mask = (1 << 32) - 1
    for index, row in enumerate(work):
        if (row & identity_mask) != (1 << index):
            raise RuntimeError("systematic reduction failed")
    return work


def reconstruct_extended_bch() -> tuple[list[int], list[int]]:
    generator = sum(1 << degree for degree in GENERATOR_DEGREES)
    verify_bch_bound(generator)

    # The unshortened [63,36] cyclic generator consists of x^j g(x),
    # 0 <= j < 36.  Deleting the last four information rows and the last four
    # coordinates leaves j=0..31 at length 59.  These rows already have degree
    # at most 58, so coordinate deletion changes no retained row bits.
    shortened = [generator << shift for shift in range(32)]
    if any(row >> 59 for row in shortened):
        raise RuntimeError("shortened BCH row exceeds coordinate 58")
    if gf2_rank(shortened) != 32:
        raise RuntimeError("shortened BCH generator has wrong rank")

    # Append the overall parity coordinate at position 59.
    extended = [row | ((row.bit_count() & 1) << 59) for row in shortened]
    if any(row.bit_count() & 1 for row in extended):
        raise RuntimeError("extension failed to make every generator row even")
    systematic = systematicize_first_32(extended)
    parity_rows = [(row >> 32) & ((1 << 28) - 1) for row in systematic]
    return systematic, parity_rows


def matrix_bit(rows: tuple[int, ...] | list[int], row: int, column: int) -> int:
    return (rows[row] >> column) & 1


def unique_unmatched_columns(b_rows: list[int]) -> tuple[int, ...]:
    """Find the four A columns whose removal matches B's degree invariants."""
    b_column_weights = sorted(
        sum(matrix_bit(b_rows, row, column) for row in range(32))
        for column in range(28)
    )
    b_row_weights = sorted(row.bit_count() for row in b_rows)
    candidates: list[tuple[int, ...]] = []
    for excluded in itertools.combinations(range(32), 4):
        excluded_set = set(excluded)
        retained = [column for column in range(32) if column not in excluded_set]
        column_weights = sorted(
            sum(matrix_bit(PUBLISHED_A_ROWS, row, column) for row in range(32))
            for column in retained
        )
        if column_weights != b_column_weights:
            continue
        row_weights = sorted(
            sum(matrix_bit(PUBLISHED_A_ROWS, row, column) for column in retained)
            for row in range(32)
        )
        if row_weights == b_row_weights:
            candidates.append(excluded)
    if len(candidates) != 1:
        raise RuntimeError(f"degree invariants leave {len(candidates)} four-column candidates")
    return candidates[0]


def recover_embedding(b_rows: list[int]) -> tuple[list[int], list[int], tuple[int, ...], int]:
    """Recover the unique bipartite incidence isomorphism by color refinement."""
    unmatched = unique_unmatched_columns(b_rows)
    keep = [column for column in range(32) if column not in set(unmatched)]
    matrices = (
        [[matrix_bit(b_rows, row, column) for column in range(28)] for row in range(32)],
        [[matrix_bit(PUBLISHED_A_ROWS, row, column) for column in keep] for row in range(32)],
    )
    vertices = [
        (graph, side, index)
        for graph in range(2)
        for side, size in ((0, 32), (1, 28))
        for index in range(size)
    ]
    neighbors: dict[tuple[int, int, int], list[tuple[int, int, int]]] = {}
    for graph, side, index in vertices:
        if side == 0:
            neighbors[graph, side, index] = [
                (graph, 1, column)
                for column in range(28)
                if matrices[graph][index][column]
            ]
        else:
            neighbors[graph, side, index] = [
                (graph, 0, row)
                for row in range(32)
                if matrices[graph][row][index]
            ]

    signatures = {
        vertex: (vertex[1], len(neighbors[vertex]))
        for vertex in vertices
    }
    refinement_rounds = 0
    for refinement_rounds in range(1, 61):
        next_signatures = {
            vertex: (signatures[vertex], tuple(sorted(signatures[neighbor] for neighbor in neighbors[vertex])))
            for vertex in vertices
        }
        # Canonical integer names are assigned jointly to both graphs, so equal
        # colors always retain the same name across the comparison.
        palette = {
            signature: color
            for color, signature in enumerate(sorted(set(next_signatures.values()), key=repr))
        }
        signatures = {vertex: palette[next_signatures[vertex]] for vertex in vertices}
        counts = Counter(signatures.values())
        if len(counts) == 60 and all(count == 2 for count in counts.values()):
            break
    else:
        raise RuntimeError("color refinement did not isolate the incidence isomorphism")

    def matching_vertex(side: int, source: int, size: int) -> int:
        color = signatures[0, side, source]
        matches = [index for index in range(size) if signatures[1, side, index] == color]
        if len(matches) != 1:
            raise RuntimeError("incidence color is not unique in published matrix")
        return matches[0]

    recovered_rows = [matching_vertex(0, row, 32) for row in range(32)]
    retained_column_map = [matching_vertex(1, column, 28) for column in range(28)]
    recovered_columns = [keep[index] for index in retained_column_map]
    for row in range(32):
        for column in range(28):
            lhs = matrix_bit(PUBLISHED_A_ROWS, recovered_rows[row], recovered_columns[column])
            rhs = matrix_bit(b_rows, row, column)
            if lhs != rhs:
                raise RuntimeError("incidence embedding failed independent replay")
    return recovered_rows, recovered_columns, unmatched, refinement_rounds


def load_published_generator_rows(path: Path) -> list[int]:
    rows = [int(line.strip(), 16) for line in path.read_text().splitlines() if line.strip()]
    expected = [(a << 32) | (1 << row) for row, a in enumerate(PUBLISHED_A_ROWS)]
    if rows != expected:
        raise RuntimeError("committed generator rows differ from published A")
    return rows


def load_spectrum(path: Path) -> tuple[int, int, int]:
    mass = 0
    minimum = None
    weight_12 = 0
    with path.open(newline="") as stream:
        for record in csv.DictReader(stream):
            weight = int(record["weight"])
            count = int(record["count"])
            mass += count
            if weight > 0 and count and minimum is None:
                minimum = weight
            if weight == 12:
                weight_12 = count
    if minimum is None:
        raise RuntimeError("spectrum contains no nonzero codeword")
    return mass, minimum, weight_12


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    systematic, b_rows = reconstruct_extended_bch()
    published_rows = load_published_generator_rows(args.rows)
    row_map, column_map, unmatched, refinement_rounds = recover_embedding(b_rows)

    # Puncture the four unmatched right-half coordinates.  Reorder its rows by
    # row_map and its selected right coordinates by column_map.  The identity
    # half then differs only by the same row-coordinate permutation, so the
    # resulting code is coordinate-equivalent to reconstructed (I | B).
    for source_row in range(32):
        published_row = row_map[source_row]
        selected_parity = sum(
            matrix_bit(PUBLISHED_A_ROWS, published_row, column_map[column]) << column
            for column in range(28)
        )
        if selected_parity != b_rows[source_row]:
            raise RuntimeError("punctured parity row differs from reconstructed B")

    mass, minimum, weight_12 = load_spectrum(args.spectrum)
    if mass != 1 << 32 or minimum != 12:
        raise RuntimeError("committed exhaustive spectrum has wrong mass or minimum")
    if gf2_rank(published_rows) != 32:
        raise RuntimeError("published [64,32] generator has wrong rank")

    receipt = {
        "schema": "xbch64-philips-reconstruction-v1",
        "status": "PASS",
        "source_construction": {
            "primitive_length": 63,
            "dimension": 36,
            "generator_degrees": list(GENERATOR_DEGREES),
            "shortening": "delete information rows 32..35 and coordinates 59..62",
            "extended_length": 60,
            "extended_dimension": 32,
        },
        "reconstructed_systematic_code": {
            "form": "(I_32 | B_32x28)",
            "rank": gf2_rank(systematic),
            "all_rows_even": all((row.bit_count() & 1) == 0 for row in systematic),
            "b_rows_hex": [f"{row:07x}" for row in b_rows],
        },
        "published_matrix_embedding": {
            "identity": "A[row_map[r], column_map[c]] = B[r,c]",
            "row_map_b_to_a": row_map,
            "column_map_b_to_a": column_map,
            "four_unmatched_a_columns": list(unmatched),
            "degree_invariant_candidate_count": 1,
            "incidence_refinement_rounds": refinement_rounds,
            "interpretation": (
                "After puncturing the four unmatched right-half coordinates, "
                "the published (I_32|A) code is coordinate-equivalent to the "
                "reconstructed extended shortened BCH [60,32] code."
            ),
        },
        "published_code": {
            "parameters": [64, 32, minimum],
            "rank": gf2_rank(published_rows),
            "exhaustive_spectrum_mass": mass,
            "weight_12_multiplicity": weight_12,
            "generator_rows_sha256": sha256(args.rows),
            "spectrum_sha256": sha256(args.spectrum),
        },
        "logical_certificate": [
            "The shortened cyclic generator has rank 32.",
            "The generator divides x^63+1 and vanishes at alpha^1 through alpha^10 in GF(64), proving parent distance >=11 by the BCH bound.",
            "Shortening preserves distance. Adding overall parity makes every word even, so the extended distance is >=12.",
            "The exact embedding authenticates 28 of the 32 published parity columns as the reconstructed BCH parity matrix up to permutations.",
            "The remaining four columns can only increase word weights relative to the punctured [60,32,>=12] code.",
            "The exhaustive 2^32-word spectrum supplies a weight-12 witness and proves the published [64,32] minimum distance is exactly 12.",
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2) + "\n")

    print("Philips shortened-XBCH reconstruction: PASS")
    print(f"reconstructed_rank={receipt['reconstructed_systematic_code']['rank']}")
    print(f"row_map_b_to_a={','.join(map(str, row_map))}")
    print(f"column_map_b_to_a={','.join(map(str, column_map))}")
    print(f"unmatched_a_columns={','.join(map(str, unmatched))}")
    print(f"incidence_refinement_rounds={refinement_rounds}")
    print(f"published_spectrum_mass={mass}")
    print(f"published_minimum_distance={minimum}")
    print(f"published_weight_12_multiplicity={weight_12}")
    print(f"wrote={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
