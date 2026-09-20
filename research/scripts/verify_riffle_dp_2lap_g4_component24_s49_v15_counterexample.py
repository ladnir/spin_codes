#!/usr/bin/env python3
"""Independently replay the support-(1,9,20), value-15 counterexample."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from analyze_systematic_group_kernel import apply, systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
RECEIPT = (
    ROOT
    / "explorations"
    / "riffle_dp_2lap_g4_component24_primary_s49_v15.json"
)
AUTONOMOUS = ROOT / "explorations" / "riffle_dp_g4_g2_autonomous_map.json"
OUTPUT = (
    ROOT
    / "explorations"
    / "riffle_dp_2lap_g4_component24_s49_v15_counterexample_verification.json"
)
MASK64 = (1 << 64) - 1
PACKET_VALUE = 15
CHARACTER = 0xA685AAC60E5ACFB3
SUPPORT_DEGREES = (1, 9, 20)
WINDOW_NODES = 24
SLOTS = 16
FULL_NODES = 32_772
REQUIRED_DISTANCE = 73


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
        raise RuntimeError("component counterexample: packet value changed")
    if receipt["component_degrees"] != list(SUPPORT_DEGREES):
        raise RuntimeError("component counterexample: support changed")
    if int(receipt["violating_character_hex"], 16) != CHARACTER:
        raise RuntimeError("component counterexample: character changed")

    autonomous = json.loads(AUTONOMOUS.read_text())
    factors = {
        row["degree"]: int(row["factor_hex"], 16)
        for row in autonomous["irreducible_components"]
    }
    annihilator = 1
    for degree in SUPPORT_DEGREES:
        annihilator = polynomial_multiply(annihilator, factors[degree])
    if annihilator != int(receipt["component_annihilator_hex"], 16):
        raise RuntimeError("component counterexample: annihilator changed")

    state_columns = systematic_state_columns()

    def step(state: int) -> int:
        return apply(state_columns, accumulate(state))

    step_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_columns = transpose(step_columns)

    def transpose_step(character: int) -> int:
        return apply(transpose_columns, character)

    support_residual = apply_polynomial(transpose_step, annihilator, CHARACTER)
    degree_residuals = {
        degree: apply_polynomial(transpose_step, factors[degree], CHARACTER)
        for degree in SUPPORT_DEGREES
    }
    if support_residual != 0:
        raise RuntimeError("component counterexample: character left the support subcode")
    if degree_residuals[1] != 0:
        raise RuntimeError("component counterexample: known degree-one membership changed")
    if degree_residuals[9] == 0 or degree_residuals[20] == 0:
        raise RuntimeError("component counterexample: character unexpectedly lies in another factor")

    weights_by_node = [0] * FULL_NODES
    first_window_bits = []
    for slot in range(SLOTS):
        state = PACKET_VALUE << (4 * slot)
        for node in range(FULL_NODES):
            state = step(state)
            bit = (CHARACTER & state).bit_count() & 1
            weights_by_node[node] += bit
            if node < WINDOW_NODES:
                first_window_bits.append(bit)

    first_window_weight = sum(first_window_bits)
    first_window_length = WINDOW_NODES * SLOTS
    if first_window_weight != 72:
        raise RuntimeError("component counterexample: 24-node weight changed")
    if first_window_weight >= REQUIRED_DISTANCE:
        raise RuntimeError("component counterexample: target is not refuted")

    block12_weights = [
        sum(weights_by_node[start : start + 12])
        for start in range(0, FULL_NODES, 12)
    ]
    block24_weights = [
        sum(weights_by_node[start : start + 24])
        for start in range(0, FULL_NODES - 12, 24)
    ]
    full_weight = sum(weights_by_node)
    full_length = FULL_NODES * SLOTS
    full_character_sum = full_length - 2 * full_weight
    if set(block12_weights) != {36}:
        raise RuntimeError("component counterexample: 12-node orbit is not constant")
    if set(block24_weights) != {72}:
        raise RuntimeError("component counterexample: 24-node orbit is not constant")
    if full_weight != 98_316 or full_character_sum != 327_720:
        raise RuntimeError("component counterexample: full-orbit values changed")

    observed_word = sum(bit << index for index, bit in enumerate(first_window_bits))
    payload = {
        "schema": "riffle-dp-2lap-g4-component24-counterexample-verification-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_COUNTEREXAMPLE_VERIFICATION",
        "source_sha256": digest(Path(__file__).resolve()),
        "certificate_receipt_sha256": digest(RECEIPT),
        "autonomous_map_sha256": digest(AUTONOMOUS),
        "component_support_mask_hex": "0x49",
        "component_degrees": list(SUPPORT_DEGREES),
        "component_annihilator_hex": hex(annihilator),
        "component_annihilator_residual_hex": hex(support_residual),
        "degree_factor_residuals_hex": {
            str(degree): hex(residual) for degree, residual in degree_residuals.items()
        },
        "character_hex": hex(CHARACTER),
        "pure_degree_one_character": True,
        "packet_value": PACKET_VALUE,
        "window_nodes": WINDOW_NODES,
        "code_length": first_window_length,
        "observed_word_limbs_hex": [
            hex((observed_word >> (64 * limb)) & MASK64) for limb in range(6)
        ],
        "observed_weight": first_window_weight,
        "complement_weight": first_window_length - first_window_weight,
        "required_two_sided_distance": REQUIRED_DISTANCE,
        "refutes_required_distance": True,
        "full_orbit": {
            "nodes": FULL_NODES,
            "bits": full_length,
            "blocks_of_12": len(block12_weights),
            "constant_12_node_weight": 36,
            "complete_blocks_of_24": len(block24_weights),
            "constant_24_node_weight": 72,
            "remaining_12_node_weight": block12_weights[-1],
            "weight": full_weight,
            "character_sum": full_character_sum,
            "absolute_bias_numerator": abs(full_character_sum),
            "absolute_bias_denominator": full_length,
            "absolute_bias": abs(full_character_sum) / full_length,
            "equals_required_value15_cap": abs(full_character_sum) * 8 == full_length * 5,
        },
        "status": "PASS_EXACT_COUNTEREXAMPLE_REPLAY",
        "scope_limitation": (
            "The witness refutes the stronger 24-node distance-73 lemma. Its "
            "full-orbit bias equals 5/8, so it does not refute the desired "
            "value-15 global character cap."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"character={hex(CHARACTER)}")
    print(f"observed_weight={first_window_weight}")
    print(f"full_orbit_bias={abs(full_character_sum)}/{full_length}")
    print(f"output={OUTPUT}")
    print("status=PASS_EXACT_COUNTEREXAMPLE_REPLAY")


if __name__ == "__main__":
    main()
