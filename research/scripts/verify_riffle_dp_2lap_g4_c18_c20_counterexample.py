#!/usr/bin/env python3
"""Independently replay the C18+C20 local-distance counterexample."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from analyze_systematic_group_kernel import apply, systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "explorations" / "riffle_dp_2lap_g4_c18_c20_v07.json"
AUTONOMOUS = ROOT / "explorations" / "riffle_dp_g4_g2_autonomous_map.json"
OUTPUT = ROOT / "explorations" / "riffle_dp_2lap_g4_c18_c20_counterexample_verification.json"
MASK64 = (1 << 64) - 1
FACTOR18 = 0x7C9C3
FACTOR20 = 0x1E1FFF
PACKET_VALUE = 7
CHARACTER = 0x852AC8FCC6B27C67
WINDOW_NODES = 12
SLOTS = 16
FULL_NODES = 32_772
REQUIRED_DISTANCE = 48


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


def transpose(columns: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(((columns[input_bit] >> output_bit) & 1) << input_bit for input_bit in range(64))
        for output_bit in range(64)
    )


def polynomial_multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        left <<= 1
        right >>= 1
    return result


def apply_polynomial(step, polynomial: int, value: int) -> int:
    result = 0
    current = value
    while polynomial:
        if polynomial & 1:
            result ^= current
        polynomial >>= 1
        current = step(current)
    return result


def main() -> None:
    receipt = json.loads(RECEIPT.read_text())
    if receipt["packet_value"] != PACKET_VALUE:
        raise RuntimeError("counterexample verifier: receipt packet value changed")
    if int(receipt["violating_character_hex"], 16) != CHARACTER:
        raise RuntimeError("counterexample verifier: receipt character changed")

    autonomous = json.loads(AUTONOMOUS.read_text())
    factors = {
        row["degree"]: int(row["factor_hex"], 16)
        for row in autonomous["irreducible_components"]
    }
    if factors[18] != FACTOR18 or factors[20] != FACTOR20:
        raise RuntimeError("counterexample verifier: component factors changed")

    state_columns = systematic_state_columns()

    def step(state: int) -> int:
        return apply(state_columns, accumulate(state))

    step_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_columns = transpose(step_columns)

    def transpose_step(character: int) -> int:
        return apply(transpose_columns, character)

    annihilator = polynomial_multiply(FACTOR18, FACTOR20)
    combined_residual = apply_polynomial(transpose_step, annihilator, CHARACTER)
    degree18_residual = apply_polynomial(transpose_step, FACTOR18, CHARACTER)
    degree20_residual = apply_polynomial(transpose_step, FACTOR20, CHARACTER)
    if combined_residual != 0:
        raise RuntimeError("counterexample verifier: character is outside C18+C20")
    if degree18_residual == 0 or degree20_residual == 0:
        raise RuntimeError("counterexample verifier: character is not genuinely mixed")

    observed_bits = []
    for slot in range(SLOTS):
        state = PACKET_VALUE << (4 * slot)
        for _ in range(WINDOW_NODES):
            state = step(state)
            observed_bits.append((CHARACTER & state).bit_count() & 1)
    if len(observed_bits) != WINDOW_NODES * SLOTS:
        raise RuntimeError("counterexample verifier: observation length changed")
    observed_word = sum(bit << index for index, bit in enumerate(observed_bits))
    observed_weight = observed_word.bit_count()
    complement_support = [
        index for index, bit in enumerate(observed_bits) if bit == 0
    ]
    complement_weight = len(complement_support)
    if observed_weight != 146 or complement_weight != 46:
        raise RuntimeError("counterexample verifier: counterexample weight changed")
    if complement_weight >= REQUIRED_DISTANCE:
        raise RuntimeError("counterexample verifier: distance target not refuted")

    weights_by_node = [0] * FULL_NODES
    for slot in range(SLOTS):
        state = PACKET_VALUE << (4 * slot)
        for node in range(FULL_NODES):
            state = step(state)
            weights_by_node[node] += (CHARACTER & state).bit_count() & 1
    block_weights = [
        sum(weights_by_node[start : start + WINDOW_NODES])
        for start in range(0, FULL_NODES, WINDOW_NODES)
    ]
    full_weight = sum(block_weights)
    full_length = FULL_NODES * SLOTS
    full_complement_weight = full_length - full_weight
    full_character_sum = full_length - 2 * full_weight
    if block_weights[0] != observed_weight:
        raise RuntimeError("counterexample verifier: first block mismatch")
    if (full_weight, full_complement_weight, full_character_sum) != (
        263_123,
        261_229,
        -1_894,
    ):
        raise RuntimeError("counterexample verifier: full-orbit replay changed")

    payload = {
        "schema": "riffle-dp-2lap-g4-c18-c20-counterexample-verification-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_COUNTEREXAMPLE_VERIFICATION",
        "certificate_receipt_sha256": digest(RECEIPT),
        "autonomous_map_sha256": digest(AUTONOMOUS),
        "verifier_source_sha256": digest(Path(__file__).resolve()),
        "packet_value": PACKET_VALUE,
        "character_hex": hex(CHARACTER),
        "component_degrees": [18, 20],
        "combined_annihilator_hex": hex(annihilator),
        "combined_annihilator_residual_hex": hex(combined_residual),
        "degree18_residual_hex": hex(degree18_residual),
        "degree20_residual_hex": hex(degree20_residual),
        "genuinely_mixed_component_character": True,
        "code_length": len(observed_bits),
        "observed_word_limbs_hex": [
            hex((observed_word >> (64 * limb)) & MASK64) for limb in range(3)
        ],
        "observed_weight": observed_weight,
        "complement_weight": complement_weight,
        "complement_support_coordinates": complement_support,
        "required_two_sided_distance": REQUIRED_DISTANCE,
        "refutes_required_distance": complement_weight < REQUIRED_DISTANCE,
        "full_orbit": {
            "nodes": FULL_NODES,
            "blocks_of_12": len(block_weights),
            "weight": full_weight,
            "complement_weight": full_complement_weight,
            "character_sum": full_character_sum,
            "absolute_bias_numerator": abs(full_character_sum),
            "absolute_bias_denominator": full_length,
            "absolute_bias": abs(full_character_sum) / full_length,
            "minimum_block_weight": min(block_weights),
            "maximum_block_weight": max(block_weights),
            "minimum_two_sided_block_weight": min(
                min(weight, WINDOW_NODES * SLOTS - weight)
                for weight in block_weights
            ),
            "blocks_below_required_two_sided_distance": sum(
                min(weight, WINDOW_NODES * SLOTS - weight) < REQUIRED_DISTANCE
                for weight in block_weights
            ),
            "interpretation": (
                "The witness violates one local 12-node block but has small bias "
                "over all 32,772 nodes. It refutes the blockwise proof lemma, not "
                "the proposed global value-7 character cap."
            ),
        },
        "status": "PASS_EXACT_COUNTEREXAMPLE_REPLAY",
        "scope_limitation": (
            "This artifact refutes the proposed value-7 local distance bound for "
            "the C18+C20 subcode. It does not refute Riffle DP-2Lap g=4 or the "
            "weaker global Fourier bound by itself."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"character={hex(CHARACTER)}")
    print(f"observed_weight={observed_weight}")
    print(f"complement_weight={complement_weight}")
    print("genuinely_mixed_component_character=True")
    print(f"output={OUTPUT}")
    print("status=PASS_EXACT_COUNTEREXAMPLE_REPLAY")


if __name__ == "__main__":
    main()
