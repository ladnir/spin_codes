#!/usr/bin/env python3
"""Exact reduced-instance stress test for Riffle ParallelAcc g=4.

The primary dynamic program counts distinct strings with a fixed packet-value
multiset. Repeated values are indistinguishable. Every such string has the same
number of preimages under a uniform permutation of labeled packets, so their
ratio is the exact probability in the construction experiment.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
OUTPUT = ROOT / "constructions" / "riffle_parallelacc_g4" / "receipts" / "goal01_reduced_exact.json"

PACKET_BITS = 4
REDUCED_PACKETS = 64
RELATIVE_THRESHOLD = Fraction(9, 100)
REDUCED_THRESHOLD = (RELATIVE_THRESHOLD.numerator * PACKET_BITS * REDUCED_PACKETS) // RELATIVE_THRESHOLD.denominator


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def multinomial_total(length: int, counts: Counter[int]) -> int:
    result = math.factorial(length)
    result //= math.factorial(length - sum(counts.values()))
    for count in counts.values():
        result //= math.factorial(count)
    return result


def accumulator_bad_count(
    *, length: int, counts: Counter[int], threshold: int, packet_bits: int
) -> int:
    """Count fixed-multiset inputs whose accumulator output has low bit weight."""

    del packet_bits  # The packet values determine the active bit lanes.
    values = tuple(sorted(counts))
    maxima = tuple(counts[value] for value in values)
    active_total = sum(maxima)

    # A key is (used-count tuple, current packet state, output bit weight).
    current: dict[tuple[tuple[int, ...], int, int], int] = {
        ((0,) * len(values), 0, 0): 1
    }
    for position in range(length):
        following: defaultdict[tuple[tuple[int, ...], int, int], int] = defaultdict(int)
        for (used, state, weight), ways in current.items():
            used_active = sum(used)
            used_zeros = position - used_active
            if used_zeros < length - active_total:
                next_weight = weight + state.bit_count()
                if next_weight <= threshold:
                    following[(used, state, next_weight)] += ways

            for index, value in enumerate(values):
                if used[index] == maxima[index]:
                    continue
                next_state = state ^ value
                next_weight = weight + next_state.bit_count()
                if next_weight > threshold:
                    continue
                next_used = list(used)
                next_used[index] += 1
                following[(tuple(next_used), next_state, next_weight)] += ways
        current = following

    return sum(
        ways
        for (used, _state, _weight), ways in current.items()
        if used == maxima
    )


def scalar_accumulator_count(length: int, input_weight: int, output_weight: int) -> int:
    """Standard zero-started accumulator input-output enumerator."""

    if input_weight == 0:
        return int(output_weight == 0)
    if output_weight == 0:
        return 0
    return math.comb(length - output_weight, input_weight // 2) * math.comb(
        output_weight - 1, (input_weight + 1) // 2 - 1
    )


def scalar_bad_count(length: int, input_weight: int, threshold: int) -> int:
    return sum(
        scalar_accumulator_count(length, input_weight, output_weight)
        for output_weight in range(threshold + 1)
    )


def probability_record(numerator: int, denominator: int) -> dict[str, object]:
    probability = Fraction(numerator, denominator)
    return {
        "numerator": str(probability.numerator),
        "denominator": str(probability.denominator),
        "log2": (
            math.log2(probability.numerator) - math.log2(probability.denominator)
            if probability
            else None
        ),
    }


def packet_values(codeword_hex: str) -> list[int]:
    codeword = int(codeword_hex, 16)
    return [
        value
        for packet in range(32)
        if (value := (codeword >> (4 * packet)) & 0xF) != 0
    ]


def authenticated_local_profile() -> tuple[Counter[int], str]:
    source = json.loads(SOURCE.read_text())
    receipt = source["receipts"][0]
    words = receipt["local_words"]
    if len(words) != 3 or len({word["codeword_hex"] for word in words}) != 1:
        raise RuntimeError("expected the first authenticated outer word to contain three equal local words")
    values = packet_values(words[0]["codeword_hex"])
    if len(values) != 11:
        raise RuntimeError("authenticated local BCH word no longer has packet support 11")
    return Counter(values), receipt["family_id"]


def main() -> None:
    local_profile, family_id = authenticated_local_profile()
    profiles = {
        "authenticated_local_support11": local_profile,
        "four_basis_pairs": Counter({1: 2, 2: 2, 4: 2, 8: 2}),
        "four_equal_pairs": Counter({3: 2, 5: 2, 6: 2, 7: 2}),
        "eight_equal_full_packets": Counter({15: 8}),
    }

    rows = []
    for name, counts in profiles.items():
        bad = accumulator_bad_count(
            length=REDUCED_PACKETS,
            counts=counts,
            threshold=REDUCED_THRESHOLD,
            packet_bits=PACKET_BITS,
        )
        total = multinomial_total(REDUCED_PACKETS, counts)
        input_bit_weight = sum(value.bit_count() * count for value, count in counts.items())
        bit_length = PACKET_BITS * REDUCED_PACKETS
        full_bit_bad = scalar_bad_count(bit_length, input_bit_weight, REDUCED_THRESHOLD)
        full_bit_total = math.comb(bit_length, input_bit_weight)
        rows.append(
            {
                "profile": name,
                "packet_value_counts": {str(value): count for value, count in sorted(counts.items())},
                "packet_support": sum(counts.values()),
                "input_bit_weight": input_bit_weight,
                "packet_permutation_bad_probability": probability_record(bad, total),
                "full_bit_permutation_bad_probability": probability_record(full_bit_bad, full_bit_total),
            }
        )

    # Independent gate: at width one the multiset DP must equal the closed form.
    scalar_gates = []
    for length, input_weight, threshold in ((12, 4, 3), (16, 5, 5), (20, 8, 7)):
        dynamic = accumulator_bad_count(
            length=length,
            counts=Counter({1: input_weight}),
            threshold=threshold,
            packet_bits=1,
        )
        formula = scalar_bad_count(length, input_weight, threshold)
        scalar_gates.append(
            {
                "length": length,
                "input_weight": input_weight,
                "threshold": threshold,
                "dynamic_count": str(dynamic),
                "formula_count": str(formula),
                "match": dynamic == formula,
            }
        )
    if not all(row["match"] for row in scalar_gates):
        raise RuntimeError("scalar accumulator formula did not match the multiset dynamic program")

    payload = {
        "schema": "riffle-parallelacc-g4-goal01-reduced-exact-v1",
        "candidate": "Riffle ParallelAcc g=4",
        "evidence_label": "EXACT_REDUCED_INSTANCE",
        "source_sha256": digest(SOURCE),
        "authenticated_outer_family": family_id,
        "experiment": {
            "packet_bits": PACKET_BITS,
            "packet_positions": REDUCED_PACKETS,
            "binary_length": PACKET_BITS * REDUCED_PACKETS,
            "failure_weight_inclusive": REDUCED_THRESHOLD,
            "permutation_law": "uniform over strings with the fixed packet-value multiset",
        },
        "scalar_formula_gates": scalar_gates,
        "profiles": rows,
        "scope_limitation": (
            "These are exact reduced-instance probabilities. They do not bound the "
            "full 524352-packet construction or cover every outer word."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

    print("candidate=Riffle ParallelAcc g=4")
    print(f"reduced_packet_count={REDUCED_PACKETS}")
    print(f"failure_weight_inclusive={REDUCED_THRESHOLD}")
    for row in rows:
        packet_log = row["packet_permutation_bad_probability"]["log2"]
        bit_log = row["full_bit_permutation_bad_probability"]["log2"]
        print(
            f"{row['profile']}: packet_log2={packet_log:.9f} "
            f"full_bit_log2={bit_log:.9f} penalty_bits={packet_log - bit_log:.9f}"
        )
    print(f"output={OUTPUT}")
    print("status=EXACT_REDUCED_INSTANCE")


if __name__ == "__main__":
    main()
