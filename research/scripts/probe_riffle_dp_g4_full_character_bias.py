#!/usr/bin/env python3
"""Bit-parallel diagnostic search for full-state packet character bias."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import (
    accumulate,
    apply_polynomial,
    build_apply,
    kernel_basis,
)
from probe_riffle_dp_g4_component_mixing import (
    FACTORS,
    SAMPLE_SPACE,
    contribution_states,
    transpose_columns,
)


ROOT = Path(__file__).resolve().parents[1]
COMPONENT_RECEIPT = ROOT / "explorations" / "riffle_dp_g4_g2_component_mixing_probe.json"
DEFAULT_OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_full_character_bias_probe.json"
INDIVIDUAL_DEGREES = (1, 2, 4, 9, 10, 18, 20)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def character_from_index(basis: tuple[int, ...], index: int) -> int:
    result = 0
    for bit, vector in enumerate(basis):
        if (index >> bit) & 1:
            result ^= vector
    return result


def bit_rows(states: np.ndarray) -> np.ndarray:
    rows = np.empty((64, len(states)), dtype=np.bool_)
    for bit in range(64):
        rows[bit] = ((states >> np.uint64(bit)) & np.uint64(1)).astype(np.bool_)
    return rows


def hill_climb(rows: np.ndarray, start: int) -> tuple[int, int]:
    character = start or 1
    active = [bit for bit in range(64) if (character >> bit) & 1]
    parity = np.logical_xor.reduce(rows[active], axis=0)
    numerator = abs(SAMPLE_SPACE - 2 * int(np.count_nonzero(parity)))
    while True:
        neighbor_ones = np.count_nonzero(rows != parity, axis=1)
        neighbor_numerators = np.abs(SAMPLE_SPACE - 2 * neighbor_ones)
        bit = int(np.argmax(neighbor_numerators))
        next_numerator = int(neighbor_numerators[bit])
        if next_numerator <= numerator:
            break
        parity ^= rows[bit]
        character ^= 1 << bit
        numerator = next_numerator
    return numerator, character


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-restarts", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.random_restarts < 0:
        raise SystemExit("full character bias: negative restart count")

    component = json.loads(COMPONENT_RECEIPT.read_text())
    component_by_degree = {
        tuple(row["component_degrees"]): row for row in component["component_rows"]
    }
    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    t_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(t_columns))
    bases = {}
    for degree in INDIVIDUAL_DEGREES:
        columns = tuple(
            apply_polynomial(transpose_apply, FACTORS[degree], 1 << bit)
            for bit in range(64)
        )
        basis = kernel_basis(columns)
        if len(basis) != degree:
            raise RuntimeError("full character bias: component dimension mismatch")
        bases[degree] = basis

    rng = random.Random(args.seed)
    rows = []
    for packet_value in range(1, 16):
        states = contribution_states(step, packet_value)
        bits = bit_rows(states)
        seeds = []
        for degree in INDIVIDUAL_DEGREES:
            component_row = component_by_degree[(degree,)]
            value_row = next(
                row
                for row in component_row["value_rows"]
                if row["packet_value"] == packet_value
            )
            seeds.append(
                {
                    "source": f"component_{degree}",
                    "character": character_from_index(
                        bases[degree],
                        value_row["maximizing_character_index"],
                    ),
                }
            )
        seeds.extend(
            {
                "source": "random",
                "character": rng.getrandbits(64) or 1,
            }
            for _ in range(args.random_restarts)
        )
        best = None
        for seed in seeds:
            numerator, character = hill_climb(bits, seed["character"])
            candidate = (numerator, character, seed["source"])
            if best is None or candidate[0] > best[0]:
                best = candidate
        assert best is not None
        rows.append(
            {
                "packet_value": packet_value,
                "maximum_found_absolute_bias_numerator": best[0],
                "maximum_found_absolute_bias_denominator": SAMPLE_SPACE,
                "maximum_found_absolute_bias": best[0] / SAMPLE_SPACE,
                "character_hex": hex(best[1]),
                "seed_source": best[2],
            }
        )

    worst = max(rows, key=lambda row: row["maximum_found_absolute_bias"])
    payload = {
        "schema": "riffle-dp-g4-g2-full-character-bias-probe-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "DIAGNOSTIC",
        "component_receipt_sha256": digest(COMPONENT_RECEIPT),
        "source_sha256": digest(Path(__file__).resolve()),
        "random_restarts_per_value": args.random_restarts,
        "seed": args.seed,
        "worst_row": worst,
        "value_rows": rows,
        "scope_limitation": (
            "Every reported character bias is exact, but hill climbing does not "
            "prove that the largest reported bias is the global maximum."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    for row in rows:
        print(
            f"value={row['packet_value']} bias={row['maximum_found_absolute_bias']:.9f} "
            f"character={row['character_hex']} seed={row['seed_source']}"
        )
    print(
        f"worst_value={worst['packet_value']} "
        f"worst_bias={worst['maximum_found_absolute_bias']:.9f}"
    )
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_FULL_CHARACTER_BIAS_SEARCH")


if __name__ == "__main__":
    main()
