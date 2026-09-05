#!/usr/bin/env python3
"""Exactly maximize the global slot second moment on pure components."""

from __future__ import annotations

import hashlib
import json
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
    INNER_NODES,
    PACKET_SLOTS,
    contribution_states,
    project,
    projection_tables,
    transpose_columns,
    walsh_transform,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "explorations" / "riffle_dp_g4_component_second_moment_certificate.json"
DEGREES = (1, 2, 4, 9, 10, 18, 20)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def character_from_index(basis: tuple[int, ...], index: int) -> int:
    character = 0
    for bit, vector in enumerate(basis):
        if (index >> bit) & 1:
            character ^= vector
    return character


def main() -> None:
    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    t_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(t_columns))
    components = []
    for degree in DEGREES:
        columns = tuple(
            apply_polynomial(transpose_apply, FACTORS[degree], 1 << bit)
            for bit in range(64)
        )
        basis = kernel_basis(columns)
        if len(basis) != degree:
            raise RuntimeError("component second moment: dimension mismatch")
        components.append(
            {
                "degree": degree,
                "basis": basis,
                "tables": projection_tables(basis),
                "value_rows": [],
            }
        )

    for packet_value in range(1, 16):
        states = contribution_states(step, packet_value)
        for component in components:
            degree = component["degree"]
            syndromes = project(states, component["tables"]).reshape(
                PACKET_SLOTS, INNER_NODES
            )
            pair_histogram = np.zeros(1 << degree, dtype=np.int64)
            for left in range(PACKET_SLOTS):
                for right in range(left + 1, PACKET_SLOTS):
                    pair_histogram += np.bincount(
                        syndromes[left] ^ syndromes[right],
                        minlength=1 << degree,
                    )
            walsh = walsh_transform(pair_histogram)
            moments = 16 * INNER_NODES + 2 * walsh
            maximum_index = int(np.argmax(moments[1:])) + 1
            maximum = int(moments[maximum_index])
            cap = 100 if packet_value == 15 else 64
            component["value_rows"].append(
                {
                    "packet_value": packet_value,
                    "maximum_node_second_moment": maximum / INNER_NODES,
                    "maximum_second_moment_sum": maximum,
                    "proposed_cap": cap,
                    "cap_passes": maximum <= cap * INNER_NODES,
                    "maximizing_character_index": maximum_index,
                    "maximizing_character_hex": hex(
                        character_from_index(component["basis"], maximum_index)
                    ),
                }
            )
        print(f"value={packet_value} complete")

    rows = []
    for component in components:
        value_rows = component.pop("value_rows")
        component.pop("basis")
        component.pop("tables")
        rows.append(
            {
                "component_degree": component["degree"],
                "all_caps_pass": all(row["cap_passes"] for row in value_rows),
                "value_rows": value_rows,
            }
        )
    payload = {
        "schema": "riffle-dp-g4-component-second-moment-certificate-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_COMPONENT_ENUMERATION",
        "source_sha256": digest(Path(__file__).resolve()),
        "node_count": INNER_NODES,
        "packet_slots": PACKET_SLOTS,
        "component_rows": rows,
        "result": "PASS" if all(row["all_caps_pass"] for row in rows) else "FAIL",
        "scope_limitation": (
            "The Walsh transforms exhaust every nonzero character in each pure "
            "irreducible component. They do not cover characters with two or more "
            "active components."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"status={payload['result']}")


if __name__ == "__main__":
    main()
