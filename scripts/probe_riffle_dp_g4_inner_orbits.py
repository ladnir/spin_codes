#!/usr/bin/env python3
"""Exact autonomous-orbit probe for the frozen Riffle DP g=4 inner.

The forward recurrence is

    V_i = Acc(U_i + S_{i-1}),
    S_i = P V_i,

where ``P`` is the systematic-right half of the committed extended BCH code.
The probe starts from each possible one-nibble impulse. It then applies zero
inputs exactly. It records output-weight crossings and exact one-nibble
turnoff opportunities. The result is an authenticated G2 search artifact,
not a probability upper bound.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


REPOSITORY = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    REPOSITORY / "explorations" / "riffle_dp_g4_g2_single_packet_orbits.json"
)
MANIFEST = REPOSITORY / "explorations" / "riffle_dp_g4_construction_manifest.json"
STATE_BITS = 64
PACKET_BITS = 4
PACKET_SLOTS = STATE_BITS // PACKET_BITS
INNER_NODES = 32_772
DISTANCE_THRESHOLD = 188_766
MASK64 = (1 << 64) - 1


def accumulate(value: int) -> int:
    value ^= value << 1
    value ^= value << 2
    value ^= value << 4
    value ^= value << 8
    value ^= value << 16
    value ^= value << 32
    return value & MASK64


def packet_drive(slot: int, value: int) -> int:
    return value << (PACKET_BITS * slot)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=INNER_NODES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 1 <= args.steps <= INNER_NODES:
        raise SystemExit("inner orbit: step count is out of range")

    columns = systematic_state_columns()
    column_digest = hashlib.sha256(
        b"".join(column.to_bytes(8, "little") for column in columns)
    ).hexdigest()
    drives = {
        packet_drive(slot, value): (slot, value)
        for slot in range(PACKET_SLOTS)
        for value in range(1, 1 << PACKET_BITS)
    }
    state_tables = []
    for byte_index in range(8):
        table = []
        for byte in range(256):
            value = 0
            for bit in range(8):
                if (byte >> bit) & 1:
                    value ^= columns[8 * byte_index + bit]
            table.append(value)
        state_tables.append(table)

    def apply_state(value: int) -> int:
        return (
            state_tables[0][value & 0xFF]
            ^ state_tables[1][(value >> 8) & 0xFF]
            ^ state_tables[2][(value >> 16) & 0xFF]
            ^ state_tables[3][(value >> 24) & 0xFF]
            ^ state_tables[4][(value >> 32) & 0xFF]
            ^ state_tables[5][(value >> 40) & 0xFF]
            ^ state_tables[6][(value >> 48) & 0xFF]
            ^ state_tables[7][value >> 56]
        )

    crossing_steps = []
    turnoffs = []
    minimum_prefix = [None] * (args.steps + 1)
    maximum_prefix = [None] * (args.steps + 1)
    orbit_digest = hashlib.sha256()
    for drive, (slot, value) in sorted(drives.items()):
        emitted = accumulate(drive)
        state = apply_state(emitted)
        cumulative = emitted.bit_count()
        crossing = 1 if cumulative > DISTANCE_THRESHOLD else None
        minimum_prefix[1] = (
            cumulative
            if minimum_prefix[1] is None
            else min(minimum_prefix[1], cumulative)
        )
        maximum_prefix[1] = (
            cumulative
            if maximum_prefix[1] is None
            else max(maximum_prefix[1], cumulative)
        )
        turnoff_drive = drives.get(state)
        if turnoff_drive is not None:
            turnoffs.append(
                {
                    "first_slot": slot,
                    "first_value": value,
                    "zero_steps_after_activation": 0,
                    "second_slot": turnoff_drive[0],
                    "second_value": turnoff_drive[1],
                    "emitted_weight_before_turnoff": cumulative,
                }
            )
        for step in range(2, args.steps + 1):
            emitted = accumulate(state)
            cumulative += emitted.bit_count()
            state = apply_state(emitted)
            if crossing is None and cumulative > DISTANCE_THRESHOLD:
                crossing = step
            minimum_prefix[step] = (
                cumulative
                if minimum_prefix[step] is None
                else min(minimum_prefix[step], cumulative)
            )
            maximum_prefix[step] = (
                cumulative
                if maximum_prefix[step] is None
                else max(maximum_prefix[step], cumulative)
            )
            turnoff_drive = drives.get(state)
            if turnoff_drive is not None and cumulative <= DISTANCE_THRESHOLD:
                turnoffs.append(
                    {
                        "first_slot": slot,
                        "first_value": value,
                        "zero_steps_after_activation": step - 1,
                        "second_slot": turnoff_drive[0],
                        "second_value": turnoff_drive[1],
                        "emitted_weight_before_turnoff": cumulative,
                    }
                )
        if crossing is None:
            crossing = args.steps + 1
        crossing_steps.append(
            {
                "slot": slot,
                "value": value,
                "packet_weight": value.bit_count(),
                "first_step_above_distance": crossing,
            }
        )
        orbit_digest.update(bytes((slot, value)))
        orbit_digest.update(crossing.to_bytes(4, "little"))
        orbit_digest.update(state.to_bytes(8, "little"))
        orbit_digest.update(cumulative.to_bytes(8, "little"))

    worst_crossing = max(row["first_step_above_distance"] for row in crossing_steps)
    best_crossing = min(row["first_step_above_distance"] for row in crossing_steps)
    worst_impulses = [
        row for row in crossing_steps if row["first_step_above_distance"] == worst_crossing
    ]
    payload = {
        "schema": "riffle-dp-g4-g2-single-packet-orbits-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT_LOCAL_SEARCH",
        "manifest_sha256": digest(MANIFEST),
        "source_sha256": digest(Path(__file__).resolve()),
        "steps": args.steps,
        "distance_threshold": DISTANCE_THRESHOLD,
        "state_matrix_sha256": column_digest,
        "orbit_digest": orbit_digest.hexdigest(),
        "impulse_count": len(crossing_steps),
        "first_step_above_distance_min": best_crossing,
        "first_step_above_distance_max": worst_crossing,
        "worst_impulses": worst_impulses,
        "exact_one_packet_turnoffs_below_distance": turnoffs,
        "minimum_autonomous_prefix_weights": minimum_prefix[1:],
        "maximum_autonomous_prefix_weights": maximum_prefix[1:],
        "minimum_autonomous_prefix_weight_at_worst_crossing": minimum_prefix[
            min(worst_crossing, args.steps)
        ],
        "maximum_autonomous_prefix_weight_at_worst_crossing": maximum_prefix[
            min(worst_crossing, args.steps)
        ],
        "scope_limitation": (
            "The search covers one initial nonzero nibble followed by zero-input nodes. "
            "It does not bound interactions among all packets of an outer word."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"single_packet_impulses={len(crossing_steps)}")
    print(f"crossing_step_min={best_crossing}")
    print(f"crossing_step_max={worst_crossing}")
    print(f"one_packet_turnoffs_below_distance={len(turnoffs)}")
    print(f"output={args.output}")
    print("status=EXACT_G2_LOCAL_SEARCH")


if __name__ == "__main__":
    main()
