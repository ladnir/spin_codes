#!/usr/bin/env python3
"""Audit the orbit identity and four-state certificates for Goal 02."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
MANIFEST = CANDIDATE / "manifest.json"
GOAL01_PRIMARY = CANDIDATE / "receipts" / "goal01_zero_symbol_primary.json"
GOAL01_SEARCH = CANDIDATE / "receipts" / "goal01_mixed_character_search.json"
GOAL01_GATE = CANDIDATE / "receipts" / "goal01_terminal_zero_gate.json"
PAIR = CANDIDATE / "receipts" / "goal02_pair_branch_primary.json"
TRIPLE = CANDIDATE / "receipts" / "goal02_triple_branch_primary.json"
PRIMARY = CANDIDATE / "receipts" / "goal02_window4_branch_primary.json"
INDEPENDENT = CANDIDATE / "receipts" / "goal02_window4_branch_independent.json"
PROOF = CANDIDATE / "proof" / "GOAL_02_GLOBAL_ORBIT_WEIGHT_PROOF.md"
OUTPUT = CANDIDATE / "receipts" / "goal02_audit.json"
GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
N = 32_772
M = 16 * N


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
    return result


def inverse_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for output in range(64):
        left = sum(
            ((columns[input_bit] >> output) & 1) << input_bit
            for input_bit in range(64)
        )
        rows.append(left | (1 << (64 + output)))
    for column in range(64):
        pivot = next(row for row in range(column, 64) if (rows[row] >> column) & 1)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(64):
            if row != column and ((rows[row] >> column) & 1):
                rows[row] ^= rows[column]
    inverse_rows = tuple(row >> 64 for row in rows)
    return tuple(
        sum(
            ((inverse_rows[input_bit] >> output_bit) & 1) << input_bit
            for input_bit in range(64)
        )
        for output_bit in range(64)
    )


def systematic_state_columns() -> tuple[int, ...]:
    generator_rows = tuple((GENERATOR << row) | (1 << 127) for row in range(64))
    left = tuple(row & MASK64 for row in generator_rows)
    right = tuple((row >> 64) & MASK64 for row in generator_rows)
    inverse_left = inverse_columns(left)
    result = tuple(apply_columns(right, inverse_left[column]) for column in range(64))
    if any(apply_columns(left, inverse_left[column]) != 1 << column for column in range(64)):
        raise RuntimeError("Goal 02 audit: systematic reconstruction failed")
    return result


def accumulate_by_bits(value: int) -> int:
    running = 0
    result = 0
    for bit in range(64):
        running ^= (value >> bit) & 1
        result |= running << bit
    return result


def transpose(columns: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(
            ((columns[input_bit] >> output_bit) & 1) << input_bit
            for input_bit in range(64)
        )
        for output_bit in range(64)
    )


def nibble_weight(value: int) -> int:
    return sum(((value >> (4 * nibble)) & 15) != 0 for nibble in range(16))


def replay_window(row: dict, s_columns: tuple[int, ...]) -> None:
    states = [int(value, 16) for value in row["minimum_window_hex"]]
    weights = [nibble_weight(state) for state in states]
    if weights != row["minimum_window_nibble_weights"]:
        raise RuntimeError("Goal 02 audit: window weights differ")
    for left, right in zip(states, states[1:]):
        if apply_columns(s_columns, left) != right:
            raise RuntimeError("Goal 02 audit: window recurrence differs")


def orbit_zero_count(character: int, s_columns: tuple[int, ...]) -> int:
    state = character
    zero_count = 0
    for _ in range(N):
        state = apply_columns(s_columns, state)
        zero_count += 16 - nibble_weight(state)
    return zero_count


def as_fraction(row: dict) -> Fraction:
    return Fraction(int(row["numerator"]), int(row["denominator"]))


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    goal01_primary = json.loads(GOAL01_PRIMARY.read_text(encoding="utf-8"))
    goal01_search = json.loads(GOAL01_SEARCH.read_text(encoding="utf-8"))
    goal01_gate = json.loads(GOAL01_GATE.read_text(encoding="utf-8"))
    pair = json.loads(PAIR.read_text(encoding="utf-8"))
    triple = json.loads(TRIPLE.read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    independent = json.loads(INDEPENDENT.read_text(encoding="utf-8"))
    for name, receipt in (
        ("pair", pair),
        ("triple", triple),
        ("primary", primary),
        ("independent", independent),
    ):
        if receipt["candidate_id"] != manifest["candidate_id"]:
            raise RuntimeError(f"Goal 02 audit: {name} candidate differs")

    state_columns = systematic_state_columns()
    t_columns = tuple(
        apply_columns(state_columns, accumulate_by_bits(1 << bit))
        for bit in range(64)
    )
    s_columns = transpose(t_columns)
    s_inverse = inverse_columns(s_columns)
    if any(apply_columns(s_columns, s_inverse[bit]) != 1 << bit for bit in range(64)):
        raise RuntimeError("Goal 02 audit: transpose recurrence is not invertible")

    # Independently replay the failed short windows.
    pair_input = int(pair["forward_witness_input_hex"], 16)
    pair_output = int(pair["forward_witness_output_hex"], 16)
    if apply_columns(s_columns, pair_input) != pair_output:
        raise RuntimeError("Goal 02 audit: pair witness recurrence differs")
    pair_weights = [nibble_weight(pair_input), nibble_weight(pair_output)]
    if pair_weights != [3, 8] or sum(pair_weights) >= 12:
        raise RuntimeError("Goal 02 audit: pair witness no longer refutes pairing")

    triple_states = [int(value, 16) for value in triple["minimum_triple_hex"]]
    triple_weights = [nibble_weight(state) for state in triple_states]
    if triple_weights != [4, 8, 5] or sum(triple_weights) >= 18:
        raise RuntimeError("Goal 02 audit: triple witness no longer refutes locality")
    for left, right in zip(triple_states, triple_states[1:]):
        if apply_columns(s_columns, left) != right:
            raise RuntimeError("Goal 02 audit: triple witness recurrence differs")

    if primary["status"] != "EXACT_WINDOW_BRANCH_CERTIFICATE":
        raise RuntimeError("Goal 02 audit: primary certificate did not pass")
    if independent["status"] != "EXACT_WINDOW4_BRANCH_INDEPENDENT":
        raise RuntimeError("Goal 02 audit: independent certificate did not pass")
    expected_primary = sum(math.comb(16, weight) * 15**weight for weight in range(1, 6))
    if expected_primary != 3_411_004_740:
        raise RuntimeError("Goal 02 audit: primary combinatorial count differs")
    if primary["low_anchor_states_checked"] != expected_primary:
        raise RuntimeError("Goal 02 audit: primary coverage count differs")
    expected_independent = math.comb(16, 5) * ((1 << 20) - 1)
    if expected_independent != 4_580_175_600:
        raise RuntimeError("Goal 02 audit: independent combinatorial count differs")
    if independent["words_checked_with_duplicates"] != expected_independent:
        raise RuntimeError("Goal 02 audit: independent coverage count differs")
    if primary["minimum_found_weight"] != 24 or independent["minimum_found_weight"] != 24:
        raise RuntimeError("Goal 02 audit: certified minimum differs")
    if primary["minimum_window_hex"] != independent["minimum_window_hex"]:
        raise RuntimeError("Goal 02 audit: minimum witnesses differ")
    replay_window(primary, s_columns)
    replay_window(independent, s_columns)

    # Replay Goal 01 coefficient zero counts through the transpose orbit identity.
    goal01_rows = goal01_primary["primary_maximizer_replays"]
    orbit_replays = []
    for row in goal01_rows:
        character = int(row["character_hex"], 16)
        zero_count = orbit_zero_count(character, s_columns)
        if zero_count != row["zero_symbol_count"]:
            raise RuntimeError(
                f"Goal 02 audit: orbit identity differs for {row['character_hex']}"
            )
        orbit_replays.append(
            {"character_hex": row["character_hex"], "zero_nibble_count": zero_count}
        )
    mixed = goal01_search["best_exact_character_found"]
    mixed_zero_count = orbit_zero_count(int(mixed["character_hex"], 16), s_columns)
    if mixed_zero_count != mixed["zero_symbol_count"]:
        raise RuntimeError("Goal 02 audit: mixed orbit identity differs")

    if N % 4 or N // 4 != 8193:
        raise RuntimeError("Goal 02 audit: orbit partition changed")
    orbit_weight_lower_bound = (N // 4) * 24
    zero_count_upper_bound = M - orbit_weight_lower_bound
    if orbit_weight_lower_bound != 6 * N or zero_count_upper_bound * 8 != 5 * M:
        raise RuntimeError("Goal 02 audit: global arithmetic differs")
    simple = goal01_gate["simple_sufficient_lemma_instantiation"]
    if simple["q_cap"] != "5/8" or simple["implied_absolute_beta_cap"] != "3/5":
        raise RuntimeError("Goal 02 audit: Goal 01 terminal gate changed")
    aggregate = as_fraction(simple["aggregate_bound"])
    if aggregate >= Fraction(1, 1 << 40):
        raise RuntimeError("Goal 02 audit: terminal-zero aggregate no longer closes")

    proof_text = PROOF.read_text(encoding="utf-8")
    for required in (
        "Four-state lemma",
        "3{,}411{,}004{,}740",
        "4{,}580{,}175{,}600",
        "q_\\chi\\le\\frac58",
        "0.980718082061",
    ):
        if required not in proof_text:
            raise RuntimeError(f"Goal 02 audit: proof omits {required!r}")

    sources = (
        ROOT / "scripts" / "certify_riffle_packetmul_2lap_g4_pair_branch.cpp",
        ROOT / "scripts" / "certify_riffle_packetmul_2lap_g4_triple_branch.cpp",
        ROOT / "scripts" / "verify_riffle_packetmul_2lap_g4_window4_branch.cpp",
        Path(__file__).resolve(),
    )
    receipts = (PAIR, TRIPLE, PRIMARY, INDEPENDENT, GOAL01_GATE)
    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal02-audit-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": manifest["candidate_id"],
        "evidence_label": "EXACT_AUDIT",
        "active_manifest_sha256": digest(MANIFEST),
        "source_sha256": {
            str(path.relative_to(ROOT)): digest(path) for path in sources
        },
        "receipt_sha256": {
            str(path.relative_to(CANDIDATE)): digest(path) for path in receipts
        },
        "proof_sha256": digest(PROOF),
        "transpose_step_columns_sha256_little_endian_u64": hashlib.sha256(
            b"".join(column.to_bytes(8, "little") for column in s_columns)
        ).hexdigest(),
        "verified": {
            "independent_transpose_reconstruction": True,
            "transpose_recurrence_invertible": True,
            "pair_local_bound_exactly_refuted": True,
            "triple_local_bound_exactly_refuted": True,
            "primary_four_state_coverage_count": True,
            "independent_four_state_cover_count": True,
            "four_state_minimum_24_replayed": True,
            "goal01_zero_symbol_orbit_identity_replayed": True,
            "global_orbit_weight_at_least_6N": True,
            "q_cap_5_over_8": True,
            "support33_terminal_zero_aggregate_below_2_to_minus_40": True,
        },
        "goal01_orbit_replays": orbit_replays,
        "goal01_best_mixed_orbit_replay": {
            "character_hex": mixed["character_hex"],
            "zero_nibble_count": mixed_zero_count,
        },
        "four_state_minimum": 24,
        "orbit_blocks": N // 4,
        "orbit_weight_lower_bound": orbit_weight_lower_bound,
        "zero_nibble_count_upper_bound": zero_count_upper_bound,
        "zero_symbol_fraction_upper_bound": "5/8",
        "status": "GOAL_02_PROVED",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"four_state_minimum={payload['four_state_minimum']}")
    print(f"orbit_weight_lower_bound={orbit_weight_lower_bound}")
    print(f"zero_symbol_fraction_upper_bound={payload['zero_symbol_fraction_upper_bound']}")
    print(f"output={OUTPUT}")
    print("status=GOAL_02_EXACT_AUDIT_PASSED")


if __name__ == "__main__":
    main()
