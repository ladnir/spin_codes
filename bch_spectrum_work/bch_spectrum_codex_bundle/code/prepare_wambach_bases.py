#!/usr/bin/env python3
"""Prepare equivalent BCH bases and certify affine-orbit lower bounds.

The companion C++ program enumerates combinations of at most four rows in
systematic generator matrices.  This script keeps the algebraic setup and the
post-search orbit accounting exact and auditable.
"""

from __future__ import annotations

import argparse
import json
import math
import struct
from pathlib import Path

from bch_quotient import (
    binary_poly_degree,
    binary_poly_divmod,
    binary_poly_mul,
)
from approx_count_low_shells import extended_row, systematic_generator_rows


MAGIC = b"BCHBAS1\n"

CODE_SPECS = {
    "p": {"m": 8, "delta": 37, "dimension": 131, "modulus": 0x14D},
    "q": {"m": 8, "delta": 39, "dimension": 123, "modulus": 0x14D},
    "b128": {"m": 7, "delta": 21, "dimension": 64, "modulus": 0x83},
    "b32": {"m": 5, "delta": 7, "dimension": 16, "modulus": 0x25},
    "b8": {"m": 3, "delta": 3, "dimension": 4, "modulus": 0xB},
}


def gf_mul_generic(a: int, b: int, m: int, modulus: int) -> int:
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << m):
            a ^= modulus
    return result & ((1 << m) - 1)


def gf_pow_generic(a: int, exponent: int, m: int, modulus: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = gf_mul_generic(result, a, m, modulus)
        a = gf_mul_generic(a, a, m, modulus)
        exponent >>= 1
    return result


def cyclotomic_coset_generic(exponent: int, length: int) -> tuple[int, ...]:
    first = exponent % length
    result: list[int] = []
    value = first
    while value not in result:
        result.append(value)
        value = (2 * value) % length
    assert value == first
    return tuple(result)


def defining_cosets_generic(delta: int, length: int) -> set[frozenset[int]]:
    return {
        frozenset(cyclotomic_coset_generic(exponent, length))
        for exponent in range(1, delta)
    }


def minimal_polynomial_generic(exponent: int, m: int, modulus: int) -> int:
    length = (1 << m) - 1
    coefficients = [1]
    for conjugate in cyclotomic_coset_generic(exponent, length):
        root = gf_pow_generic(2, conjugate, m, modulus)
        updated = [0] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] ^= gf_mul_generic(coefficient, root, m, modulus)
            updated[degree + 1] ^= coefficient
        coefficients = updated
    if any(coefficient not in (0, 1) for coefficient in coefficients):
        raise ArithmeticError("minimal polynomial did not descend to GF(2)")
    return sum(coefficient << degree for degree, coefficient in enumerate(coefficients))


def generator_polynomial_generic(m: int, delta: int, modulus: int) -> int:
    length = (1 << m) - 1
    generator = 1
    for coset in sorted(defining_cosets_generic(delta, length), key=min):
        generator = binary_poly_mul(
            generator, minimal_polynomial_generic(min(coset), m, modulus)
        )
    return generator


def equivalent_generator(code_name: str, exponent: int) -> int:
    """Generator when the coordinate primitive element is alpha**exponent."""

    spec = CODE_SPECS[code_name]
    m = spec["m"]
    length = (1 << m) - 1
    cosets = {
        frozenset((exponent * value) % length for value in coset)
        for coset in defining_cosets_generic(spec["delta"], length)
    }
    generator = 1
    for coset in sorted(cosets, key=min):
        generator = binary_poly_mul(
            generator, minimal_polynomial_generic(min(coset), m, spec["modulus"])
        )
    return generator


def distinct_bases(code_name: str) -> list[tuple[int, list[int]]]:
    spec = CODE_SPECS[code_name]
    length = (1 << spec["m"]) - 1
    dimension = spec["dimension"]
    expected_degree = length - dimension
    seen: set[int] = set()
    result: list[tuple[int, list[int]]] = []

    for exponent in range(1, length):
        if math.gcd(exponent, length) != 1:
            continue
        generator = equivalent_generator(code_name, exponent)
        assert binary_poly_degree(generator) == expected_degree
        if generator in seen:
            continue
        seen.add(generator)
        rows = systematic_generator_rows(
            [
                (generator << shift)
                | (((generator << shift).bit_count() & 1) << length)
                for shift in range(dimension)
            ]
        )
        pivots = [((row & -row).bit_length() - 1) for row in rows]
        assert pivots == list(range(dimension))
        assert all(row.bit_count() % 2 == 0 for row in rows)
        result.append((exponent, rows))

    expected_variants = sum(
        1 for exponent in range(1, length) if math.gcd(exponent, length) == 1
    ) // spec["m"]
    assert len(result) == expected_variants, (len(result), expected_variants)
    return result


def write_bases(code_name: str, output: Path) -> None:
    bases = distinct_bases(code_name)
    total_length = 1 << CODE_SPECS[code_name]["m"]
    dimension = len(bases[0][1])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as stream:
        stream.write(MAGIC)
        stream.write(struct.pack("<III", total_length, dimension, len(bases)))
        for exponent, rows in bases:
            stream.write(struct.pack("<I", exponent))
            for row in rows:
                stream.write(row.to_bytes(32, "little"))
    print(
        f"wrote {len(bases)} distinct {code_name.upper()} bases "
        f"({dimension} rows each) to {output}"
    )
    print("primitive exponents:", " ".join(str(exponent) for exponent, _ in bases))


def gf_pow_table(m: int, modulus: int) -> tuple[list[int], dict[int, int]]:
    length = (1 << m) - 1
    powers: list[int] = []
    logs: dict[int, int] = {}
    value = 1
    for exponent in range(length):
        powers.append(value)
        logs[value] = exponent
        value = gf_mul_generic(value, 2, m, modulus)
    assert value == 1 and len(logs) == length
    return powers, logs


def affine_image(
    word: int,
    a: int,
    b: int,
    powers: list[int],
    logs: dict[int, int],
    m: int,
    modulus: int,
) -> int:
    length = (1 << m) - 1
    image = 0
    punctured = word & ((1 << length) - 1)
    while punctured:
        low = punctured & -punctured
        coordinate = low.bit_length() - 1
        target = gf_mul_generic(a, powers[coordinate], m, modulus) ^ b
        image |= 1 << (length if target == 0 else logs[target])
        punctured ^= low
    if (word >> length) & 1:
        target = b
        image |= 1 << (length if target == 0 else logs[target])
    return image


def read_hits(path: Path) -> dict[int, set[int]]:
    hits: dict[int, set[int]] = {}
    with path.open("r", encoding="ascii") as stream:
        header = stream.readline().strip().split("\t")
        if header != ["weight", "limb0", "limb1", "limb2", "limb3"]:
            raise ValueError(f"unexpected hit-file header: {header}")
        for line in stream:
            fields = line.strip().split("\t")
            if not fields or fields == [""]:
                continue
            weight = int(fields[0])
            limbs = [int(value, 16) for value in fields[1:]]
            word = sum(limb << (64 * index) for index, limb in enumerate(limbs))
            assert word.bit_count() == weight
            hits.setdefault(weight, set()).add(word)
    return hits


def affine_canonical(
    word: int,
    powers: list[int],
    logs: dict[int, int],
    m: int,
    modulus: int,
    multiplication: list[bytes],
    inverses: list[int],
) -> tuple[int, int]:
    """Return an AGL-canonical word and the affine stabilizer size.

    Each ordered pair of distinct support points defines the unique affine map
    that sends that pair to (0, 1).  The minimum of these normalized images is
    a complete orbit invariant.  Its multiplicity equals the stabilizer size.
    """

    length = (1 << m) - 1
    support: list[int] = []
    punctured = word & ((1 << length) - 1)
    while punctured:
        low = punctured & -punctured
        support.append(powers[low.bit_length() - 1])
        punctured ^= low
    if (word >> length) & 1:
        support.append(0)

    best: int | None = None
    multiplicity = 0
    for x in support:
        for y in support:
            if x == y:
                continue
            a = inverses[x ^ y]
            b = multiplication[a][x]
            image = 0
            row = multiplication[a]
            for z in support:
                target = row[z] ^ b
                image |= 1 << (length if target == 0 else logs[target])
            if best is None or image < best:
                best = image
                multiplicity = 1
            elif image == best:
                multiplicity += 1
    assert best is not None and multiplicity > 0
    return best, multiplicity


def certify_orbits(
    code_name: str,
    hits_path: Path,
    output: Path,
    selected_weights: set[int] | None,
) -> None:
    spec = CODE_SPECS[code_name]
    m = spec["m"]
    q = 1 << m
    length = q - 1
    standard_generator = generator_polynomial_generic(m, spec["delta"], spec["modulus"])
    hits = read_hits(hits_path)
    powers, logs = gf_pow_table(m, spec["modulus"])
    multiplication = [
        bytes(gf_mul_generic(a, b, m, spec["modulus"]) for b in range(q))
        for a in range(q)
    ]
    inverses = [0] * q
    for value in range(1, q):
        inverses[value] = next(
            candidate for candidate in range(1, q)
            if multiplication[value][candidate] == 1
        )
    report: dict[str, object] = {
        "code": code_name,
        "method": "affine closure of <=4-systematic-row discoveries",
        "completeness_claim": False,
        "shells": {},
    }

    for weight in sorted(hits):
        if selected_weights is not None and weight not in selected_weights:
            continue
        groups: dict[int, dict[str, object]] = {}
        for representative in hits[weight]:
            punctured = representative & ((1 << length) - 1)
            assert binary_poly_divmod(punctured, standard_generator)[1] == 0
            assert representative.bit_count() % 2 == 0
            canonical, stabilizer_size = affine_canonical(
                representative,
                powers,
                logs,
                m,
                spec["modulus"],
                multiplication,
                inverses,
            )
            record = groups.setdefault(
                canonical,
                {
                    "representative_hex": f"{representative:064x}",
                    "canonical_hex": f"{canonical:064x}",
                    "orbit_size": (q * (q - 1)) // stabilizer_size,
                    "stabilizer_size": stabilizer_size,
                    "search_hits_in_orbit": 0,
                },
            )
            assert record["stabilizer_size"] == stabilizer_size
            record["search_hits_in_orbit"] += 1

        orbit_records = list(groups.values())
        lower_bound = sum(int(record["orbit_size"]) for record in orbit_records)

        report["shells"][str(weight)] = {
            "distinct_search_hits": len(hits[weight]),
            "distinct_affine_orbits": len(orbit_records),
            "rigorous_lower_bound": lower_bound,
            "orbits": orbit_records,
        }
        print(
            f"weight {weight}: {len(hits[weight])} search hits, "
            f"{len(orbit_records)} affine orbits, A_{weight} >= {lower_bound}"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"wrote orbit certificate data to {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--code", choices=tuple(CODE_SPECS), required=True)
    prepare.add_argument("--output", type=Path, required=True)

    orbits = subparsers.add_parser("orbits")
    orbits.add_argument("--code", choices=tuple(CODE_SPECS), required=True)
    orbits.add_argument("--hits", type=Path, required=True)
    orbits.add_argument("--output", type=Path, required=True)
    orbits.add_argument("--weights", type=int, nargs="+")

    args = parser.parse_args()
    if args.command == "prepare":
        write_bases(args.code, args.output)
    else:
        certify_orbits(
            args.code,
            args.hits,
            args.output,
            None if args.weights is None else set(args.weights),
        )


if __name__ == "__main__":
    main()
