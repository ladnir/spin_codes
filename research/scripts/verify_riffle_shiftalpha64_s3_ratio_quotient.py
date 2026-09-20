#!/usr/bin/env python3
"""Verify the S3 ratio quotient used for equal-shell occupation triples."""

from __future__ import annotations

import itertools
import json
import random
import struct
import sys
from collections import Counter
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as field  # noqa: E402


def ordered_ratios(triple: tuple[int, int, int]) -> list[int]:
    """Return second/first for all six orderings of a nonzero XOR triple."""
    return [
        field.field_multiply(second, field.field_inverse(first))
        for first, second, _ in itertools.permutations(triple)
    ]


def algebraic_ratios(ratio: int) -> list[int]:
    inverse_ratio = field.field_inverse(ratio)
    complement = 1 ^ ratio
    inverse_complement = field.field_inverse(complement)
    return [
        ratio,
        inverse_ratio,
        complement,
        inverse_complement,
        field.field_multiply(ratio, inverse_complement),
        field.field_multiply(complement, inverse_ratio),
    ]


def canonical_ratio(ratio: int) -> int:
    return min(algebraic_ratios(ratio))


def canonical_triple(first: int, second: int, third: int) -> int:
    if first ^ second ^ third:
        raise RuntimeError("schedule factors do not sum to zero")
    return min(ordered_ratios((first, second, third)))


def invariant_triple(first: int, second: int, third: int) -> int:
    """Return the symmetric degree-six quotient invariant."""
    if first ^ second ^ third or not first or not second or not third:
        raise RuntimeError("invariant requires a nonzero XOR triple")
    numerator_base = (
        field.field_multiply(first, first)
        ^ field.field_multiply(first, second)
        ^ field.field_multiply(second, second)
    )
    numerator = field.field_multiply(
        field.field_multiply(numerator_base, numerator_base), numerator_base
    )
    product = field.field_multiply(field.field_multiply(first, second), third)
    denominator = field.field_multiply(product, product)
    return field.field_multiply(numerator, field.field_inverse(denominator))


def audit_schedule(path: Path, blocks: int) -> dict[str, int | bool | str]:
    """Independently enumerate a small schedule and compare its binary table."""
    gamma = [1]
    for _ in range(1, blocks):
        gamma.append(field.field_multiply(gamma[-1], 2))
    expected: Counter[int] = Counter()
    expected_p0: Counter[int] = Counter()
    for third_difference in range(2, blocks):
        for second_difference in range(1, third_difference):
            key = invariant_triple(
                gamma[second_difference] ^ gamma[third_difference],
                1 ^ gamma[third_difference],
                1 ^ gamma[second_difference],
            )
            expected[key] += blocks - third_difference
    for difference in range(1, blocks):
        key = invariant_triple(1 ^ gamma[difference], gamma[difference], 1)
        expected_p0[key] += blocks - difference

    raw = path.read_bytes()
    if len(raw) % 24:
        raise RuntimeError("schedule table has a partial record")
    actual: dict[int, tuple[int, int]] = {}
    for offset in range(0, len(raw), 24):
        key, data, p0 = struct.unpack_from("<QQQ", raw, offset)
        if key in actual:
            raise RuntimeError("schedule table repeats a canonical key")
        actual[key] = (data, p0)
    all_keys = set(expected) | set(expected_p0)
    wanted = {key: (expected[key], expected_p0[key]) for key in all_keys}
    if actual != wanted:
        raise RuntimeError("native schedule table differs from direct enumeration")
    return {
        "path": str(path),
        "blocks": blocks,
        "records": len(actual),
        "data_mass": sum(expected.values()),
        "p0_mass": sum(expected_p0.values()),
        "native_table_matches_direct_enumeration": True,
    }


def main() -> None:
    rng = random.Random(0x533351)
    trials = 128
    special_orbits = 0
    invariant_to_canonical: dict[int, int] = {}
    synthetic_schedule: Counter[int] = Counter()
    synthetic_triples: list[tuple[int, int, int]] = []
    for trial in range(trials):
        first = rng.randrange(1, 1 << 64)
        second = rng.randrange(1, 1 << 64)
        while second == first:
            second = rng.randrange(1, 1 << 64)
        third = first ^ second
        triple = (first, second, third)
        direct = sorted(ordered_ratios(triple))
        ratio = field.field_multiply(second, field.field_inverse(first))
        algebraic = sorted(algebraic_ratios(ratio))
        if direct != algebraic:
            raise RuntimeError("algebraic S3 orbit differs from six orderings")
        canonical = canonical_ratio(ratio)
        if any(canonical_ratio(value) != canonical for value in direct):
            raise RuntimeError("canonical ratio changed within an S3 orbit")
        distinct = len(set(direct))
        if distinct not in (2, 6):
            raise RuntimeError("unexpected S3 orbit size")
        if distinct == 2:
            special_orbits += 1
        invariant = invariant_triple(*triple)
        previous = invariant_to_canonical.setdefault(invariant, canonical)
        if previous != canonical:
            raise RuntimeError("invariant merged distinct sampled S3 orbits")
        for ordering in itertools.permutations(triple):
            if invariant_triple(*ordering) != invariant:
                raise RuntimeError("invariant changed under a permutation")
        synthetic_triples.append(triple)
        if trial < 64:
            synthetic_schedule[direct[trial % 6]] += 1 + trial % 5

    direct_total = 0
    canonical_schedule: Counter[int] = Counter()
    for ratio, multiplicity in synthetic_schedule.items():
        canonical_schedule[canonical_ratio(ratio)] += multiplicity
    quotient_total = 0
    for triple in synthetic_triples:
        ratios = ordered_ratios(triple)
        direct_total += sum(synthetic_schedule[ratio] for ratio in ratios)
        distinct = len(set(ratios))
        factor = 6 // distinct
        quotient_total += factor * canonical_schedule[min(ratios)]
    if quotient_total != direct_total:
        raise RuntimeError("canonical schedule sum differs from direct six-ratio sum")

    payload = {
        "schema": "riffle-shiftalpha64-s3-ratio-quotient-audit-v1",
        "random_seed_hex": hex(0x533351),
        "field": "GF(2^64) with reduction polynomial x^64+x^4+x^3+x+1",
        "trials": trials,
        "synthetic_schedule_keys": len(synthetic_schedule),
        "direct_total": direct_total,
        "quotient_total": quotient_total,
        "special_size_two_orbits_observed": special_orbits,
        "validation": {
            "six_permutations_match_algebraic_ratio_orbit": True,
            "canonical_key_is_constant_on_each_orbit": True,
            "orbit_size_is_two_or_six": True,
            "canonical_schedule_sum_matches_direct_sum": True,
            "symmetric_invariant_is_constant_under_permutation": True,
            "sampled_invariant_keys_do_not_merge_s3_orbits": True,
        },
        "scope": (
            "Algebraic audit of the S3 quotient. This script does not enumerate "
            "BCH relations or the production shifted schedule."
        ),
    }
    if len(sys.argv) not in (1, 3):
        raise SystemExit("usage: verifier [small-schedule.bin blocks]")
    if len(sys.argv) == 3:
        payload["small_schedule_audit"] = audit_schedule(
            Path(sys.argv[1]), int(sys.argv[2])
        )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
