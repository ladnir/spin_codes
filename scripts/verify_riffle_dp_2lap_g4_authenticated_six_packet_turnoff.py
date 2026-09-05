#!/usr/bin/env python3
"""Independently verify the authenticated six-packet turnoff family."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
RECEIPT = EXPLORATIONS / "riffle_dp_2lap_g4_authenticated_six_packet_turnoff.json"
OUTER = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
OUTPUT = (
    EXPLORATIONS
    / "riffle_dp_2lap_g4_authenticated_six_packet_turnoff_verification.json"
)
INNER_NODES = 32_772
PACKET_POSITIONS = 524_352
MASK64 = (1 << 64) - 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accumulate(value: int) -> int:
    value ^= value << 1
    value ^= value << 2
    value ^= value << 4
    value ^= value << 8
    value ^= value << 16
    value ^= value << 32
    return value & MASK64


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = (value & -value).bit_length() - 1
        result ^= columns[bit]
        value &= value - 1
    return result


def run_lap(
    inputs: list[int], initial_state: int, columns: tuple[int, ...]
) -> tuple[list[int], int]:
    state = initial_state
    outputs = []
    for word in inputs:
        output = accumulate(state ^ word)
        outputs.append(output)
        state = apply_columns(columns, output)
    return outputs, state


def target_counts(family_id: str) -> Counter[int]:
    outer = json.loads(OUTER.read_text())
    row = next(receipt for receipt in outer["receipts"] if receipt["family_id"] == family_id)
    counts: Counter[int] = Counter()
    for local_word in row["local_words"]:
        word = int(local_word["codeword_hex"], 16)
        counts.update(
            (word >> shift) & 15
            for shift in range(0, 128, 4)
            if (word >> shift) & 15
        )
    return counts


def main() -> None:
    receipt = json.loads(RECEIPT.read_text())
    inputs = [0] * INNER_NODES
    for row in receipt["canonical_nonzero_input_nodes"]:
        inputs[row["node"]] = int(row["drive_hex"], 16)
    columns = systematic_state_columns()
    first, first_terminal = run_lap(inputs, 0, columns)
    final, final_terminal = run_lap(first, first_terminal, columns)
    final_weight = sum(word.bit_count() for word in final)
    if first_terminal != 0 or final_terminal != 0:
        raise RuntimeError("authenticated turnoff verifier: terminal mismatch")
    if final_weight != receipt["two_lap_weight"]:
        raise RuntimeError("authenticated turnoff verifier: weight mismatch")

    observed: Counter[int] = Counter()
    for word in inputs:
        observed.update(
            (word >> shift) & 15
            for shift in range(0, 64, 4)
            if (word >> shift) & 15
        )
    expected = target_counts(receipt["outer_family_id"])
    if observed != expected or sum(observed.values()) != 33:
        raise RuntimeError("authenticated turnoff verifier: outer multiset mismatch")

    episode_count = receipt["episode_count"]
    span = receipt["occupied_span_nodes"]
    gaps = math.comb(INNER_NODES - span + episode_count, episode_count)
    labelings = math.prod(math.factorial(count) for count in expected.values())
    probability = Fraction(
        gaps * labelings,
        math.prod(range(PACKET_POSITIONS - 32, PACKET_POSITIONS + 1)),
    )
    lower = receipt["lower_family"]
    recorded = Fraction(
        int(lower["probability_numerator"]),
        int(lower["probability_denominator"]),
    )
    if probability != recorded or probability >= Fraction(1, 1 << 40):
        raise RuntimeError("authenticated turnoff verifier: probability mismatch")

    payload = {
        "schema": "riffle-dp-2lap-g4-authenticated-six-packet-turnoff-verification-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_INDEPENDENT_VERIFICATION",
        "receipt_sha256": digest(RECEIPT),
        "outer_family_id": receipt["outer_family_id"],
        "verified_packet_support": sum(expected.values()),
        "verified_first_lap_terminal_hex": f"0x{first_terminal:016x}",
        "verified_second_lap_terminal_hex": f"0x{final_terminal:016x}",
        "verified_two_lap_weight": final_weight,
        "verified_gap_placements": str(gaps),
        "verified_probability_numerator": str(probability.numerator),
        "verified_probability_denominator": str(probability.denominator),
        "exceeds_2^-40": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"outer_family_id={receipt['outer_family_id']}")
    print(f"verified_two_lap_weight={final_weight}")
    print(f"verified_gap_placements={gaps}")
    print("exceeds_2^-40=False")
    print(f"output={OUTPUT}")
    print("status=EXACT_INDEPENDENT_VERIFICATION")


if __name__ == "__main__":
    main()
