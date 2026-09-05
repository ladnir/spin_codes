#!/usr/bin/env python3
"""Structured diagnostic search for mixed PacketMul zero-symbol characters."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, build_apply
from probe_riffle_dp_g4_component_mixing import INNER_NODES, PACKET_SLOTS, SAMPLE_SPACE


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
MANIFEST = CANDIDATE / "manifest.json"
PRIMARY = CANDIDATE / "receipts" / "goal01_zero_symbol_primary.json"
OUTPUT = CANDIDATE / "receipts" / "goal01_mixed_character_search.json"
KNOWN_WITNESS_FILES = (
    ROOT / "explorations" / "riffle_dp_g4_global_second_moment_probe.json",
    ROOT / "explorations" / "riffle_dp_2lap_g4_c18_c20_counterexample_verification.json",
)
RANDOM_SEED = 20260821
LOW_COMPONENT_DIMENSION = 3  # degrees 1 and 2


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def basis_contributions(step) -> np.ndarray:
    result = np.empty((4, SAMPLE_SPACE), dtype=np.uint64)
    for packet_bit in range(4):
        cursor = 0
        for slot in range(PACKET_SLOTS):
            state = (1 << packet_bit) << (4 * slot)
            for _ in range(INNER_NODES):
                state = step(state)
                result[packet_bit, cursor] = state
                cursor += 1
    return result


def standard_character_coefficient_rows(basis_states: np.ndarray) -> np.ndarray:
    rows = np.empty((64, SAMPLE_SPACE), dtype=np.uint8)
    for character_bit in range(64):
        row = np.zeros(SAMPLE_SPACE, dtype=np.uint8)
        shift = np.uint64(character_bit)
        for packet_bit in range(4):
            row |= (
                ((basis_states[packet_bit] >> shift) & np.uint64(1)).astype(np.uint8)
                << packet_bit
            )
        rows[character_bit] = row
    return rows


def character_coefficients(rows: np.ndarray, character: int) -> np.ndarray:
    active = [bit for bit in range(64) if (character >> bit) & 1]
    if not active:
        return np.zeros(SAMPLE_SPACE, dtype=np.uint8)
    return np.bitwise_xor.reduce(rows[active], axis=0)


def build_coordinate_reducer(component_basis: tuple[int, ...]):
    pivots = [0] * 64
    representations = [0] * 64
    for basis_index, vector in enumerate(component_basis):
        value = vector
        representation = 1 << basis_index
        while value:
            pivot = value.bit_length() - 1
            if pivots[pivot]:
                value ^= pivots[pivot]
                representation ^= representations[pivot]
            else:
                pivots[pivot] = value
                representations[pivot] = representation
                break
        if value == 0:
            raise RuntimeError("component bases are not independent")
    if sum(bool(value) for value in pivots) != 64:
        raise RuntimeError("component bases do not span the character space")

    def coordinates(character: int) -> int:
        value = character
        result = 0
        while value:
            pivot = value.bit_length() - 1
            if not pivots[pivot]:
                raise RuntimeError("character is outside the component basis")
            value ^= pivots[pivot]
            result ^= representations[pivot]
        return result

    return coordinates


def stats(character: int, coefficients: np.ndarray) -> dict:
    zero_count = int(np.count_nonzero(coefficients == 0))
    signed_numerator = 16 * zero_count - SAMPLE_SPACE
    return {
        "character_hex": hex(character),
        "zero_symbol_count": zero_count,
        "q_numerator": zero_count,
        "q_denominator": SAMPLE_SPACE,
        "signed_beta_numerator": signed_numerator,
        "absolute_beta_numerator": abs(signed_numerator),
        "beta_denominator": 15 * SAMPLE_SPACE,
    }


def hill_climb(
    rows: np.ndarray,
    start: int,
    coordinates,
    standard_coordinate_masks: tuple[int, ...],
    constrain_high_component: bool,
) -> tuple[dict, int]:
    character = start or 1
    coordinate_word = coordinates(character)
    coefficients = character_coefficients(rows, character)
    current = stats(character, coefficients)
    steps = 0
    while True:
        neighbor_zero_counts = np.count_nonzero(rows == coefficients, axis=1)
        neighbor_objectives = np.abs(16 * neighbor_zero_counts - SAMPLE_SPACE)
        order = np.argsort(neighbor_objectives)[::-1]
        chosen = None
        for raw_bit in order:
            bit = int(raw_bit)
            next_character = character ^ (1 << bit)
            if next_character == 0:
                continue
            next_coordinates = coordinate_word ^ standard_coordinate_masks[bit]
            if constrain_high_component and (next_coordinates >> LOW_COMPONENT_DIMENSION) == 0:
                continue
            next_objective = int(neighbor_objectives[bit])
            if next_objective > current["absolute_beta_numerator"]:
                chosen = (bit, next_character, next_coordinates)
            break
        if chosen is None:
            break
        bit, character, coordinate_word = chosen
        coefficients ^= rows[bit]
        current = stats(character, coefficients)
        steps += 1
        if steps > 128:
            raise RuntimeError("coordinate ascent exceeded its strict-improvement bound")
    return current, steps


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-restarts", type=int, default=64)
    parser.add_argument("--xor-hill-starts", type=int, default=16)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.random_restarts < 0 or args.xor_hill_starts < 0:
        raise SystemExit("restart counts must be nonnegative")

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    component_rows = primary["pure_component_rows"]
    component_basis = tuple(
        int(vector, 16)
        for component in component_rows
        for vector in component["basis_hex"]
    )
    coordinates = build_coordinate_reducer(component_basis)
    standard_coordinate_masks = tuple(coordinates(1 << bit) for bit in range(64))

    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    basis_states = basis_contributions(step)
    basis_digest = hashlib.sha256(
        basis_states.astype("<u8", copy=False).tobytes()
    ).hexdigest()
    if basis_digest != primary["basis_contribution_sha256_little_endian_u64"]:
        raise RuntimeError("search reconstructed a different contribution map")
    rows = standard_character_coefficient_rows(basis_states)

    seed_sources: dict[int, set[str]] = {}

    def add_seed(character: int, source: str) -> None:
        if character:
            seed_sources.setdefault(character, set()).add(source)

    q_maximizers = []
    for component in component_rows:
        q_character = int(component["maximum_zero_symbol"]["character_hex"], 16)
        beta_character = int(component["maximum_absolute_beta"]["character_hex"], 16)
        q_maximizers.append(q_character)
        add_seed(q_character, f"pure_degree_{component['degree']}_max_q")
        add_seed(beta_character, f"pure_degree_{component['degree']}_max_abs_beta")

    xor_characters = []
    for mask in range(1, 1 << len(q_maximizers)):
        character = 0
        for component_index, component_character in enumerate(q_maximizers):
            if (mask >> component_index) & 1:
                character ^= component_character
        xor_characters.append(character)
        add_seed(character, f"xor_q_maximizers_mask_{mask:02x}")

    known_witness_digests = {}
    for path in KNOWN_WITNESS_FILES:
        known_witness_digests[str(path.relative_to(ROOT))] = digest(path)
        for character_hex in re.findall(r"0x[0-9a-fA-F]{1,16}", path.read_text(encoding="utf-8")):
            add_seed(int(character_hex, 16), f"paused_candidate:{path.name}")

    rng = random.Random(RANDOM_SEED)
    random_characters = []
    while len(random_characters) < args.random_restarts:
        character = rng.getrandbits(64)
        if character:
            random_characters.append(character)
            add_seed(character, f"pseudorandom_restart_{len(random_characters)-1}")

    seed_rows = []
    seed_stats_by_character = {}
    for character in sorted(seed_sources):
        row = stats(character, character_coefficients(rows, character))
        row["sources"] = sorted(seed_sources[character])
        row["component_coordinate_mask_hex"] = hex(coordinates(character))
        seed_rows.append(row)
        seed_stats_by_character[character] = row

    top_xor = sorted(
        set(xor_characters),
        key=lambda character: seed_stats_by_character[character]["absolute_beta_numerator"],
        reverse=True,
    )[: args.xor_hill_starts]
    hill_starts: dict[int, set[str]] = {}

    def add_hill_start(character: int, source: str) -> None:
        hill_starts.setdefault(character, set()).add(source)

    for component in component_rows:
        add_hill_start(
            int(component["maximum_zero_symbol"]["character_hex"], 16),
            f"pure_degree_{component['degree']}_max_q",
        )
        add_hill_start(
            int(component["maximum_absolute_beta"]["character_hex"], 16),
            f"pure_degree_{component['degree']}_max_abs_beta",
        )
    for character in top_xor:
        add_hill_start(character, "top_xor_combination")
    for character, sources in seed_sources.items():
        if any(source.startswith("paused_candidate:") for source in sources):
            add_hill_start(character, "paused_candidate_witness")
    for character in random_characters:
        add_hill_start(character, "pseudorandom_constrained")

    hill_rows = []
    for start in sorted(hill_starts):
        start_coordinates = coordinates(start)
        constrain = (start_coordinates >> LOW_COMPONENT_DIMENSION) != 0
        final, steps = hill_climb(
            rows,
            start,
            coordinates,
            standard_coordinate_masks,
            constrain,
        )
        final.update(
            {
                "start_character_hex": hex(start),
                "start_sources": sorted(hill_starts[start]),
                "strict_improvement_steps": steps,
                "high_component_constraint": constrain,
                "final_component_coordinate_mask_hex": hex(
                    coordinates(int(final["character_hex"], 16))
                ),
            }
        )
        if constrain and (
            int(final["final_component_coordinate_mask_hex"], 16)
            >> LOW_COMPONENT_DIMENSION
        ) == 0:
            raise RuntimeError("constrained climb descended into C1+C2")
        hill_rows.append(final)

    best_seed = max(seed_rows, key=lambda row: row["absolute_beta_numerator"])
    best_hill = max(hill_rows, key=lambda row: row["absolute_beta_numerator"])
    best = max((best_seed, best_hill), key=lambda row: row["absolute_beta_numerator"])
    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal01-mixed-search-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": manifest["candidate_id"],
        "evidence_label": "DIAGNOSTIC_SEARCH_EXACT_WITNESSES",
        "active_manifest_sha256": digest(MANIFEST),
        "primary_receipt_sha256": digest(PRIMARY),
        "basis_contribution_sha256_little_endian_u64": basis_digest,
        "random_seed": RANDOM_SEED,
        "random_restarts": args.random_restarts,
        "xor_combinations_evaluated": len(set(xor_characters)),
        "xor_combinations_used_as_hill_starts": len(top_xor),
        "known_witness_file_sha256": known_witness_digests,
        "objective": "maximize abs(16*#{p:a_chi(p)=0}-M)",
        "constraint": (
            "Every climb whose start has a component above degrees 1 and 2 rejects "
            "a flip that would remove all such high-component coordinates."
        ),
        "seed_evaluations": seed_rows,
        "hill_climb_rows": hill_rows,
        "best_exact_character_found": best,
        "scope_limitation": (
            "All displayed counts are exact for their displayed characters. Coordinate "
            "ascent and pseudorandom restarts do not upper-bound unseen mixed characters."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"seed_characters_evaluated={len(seed_rows)}")
    print(f"hill_climbs={len(hill_rows)}")
    print(f"best_character={best['character_hex']}")
    print(f"best_zero_symbol_count={best['zero_symbol_count']}")
    print(
        f"best_abs_beta={best['absolute_beta_numerator']}/{best['beta_denominator']}"
    )
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_MIXED_ZERO_SYMBOL_SEARCH")


if __name__ == "__main__":
    main()
