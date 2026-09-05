#!/usr/bin/env python3
"""Exact component-projection mixing probe for random packet locations."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import (
    accumulate,
    apply_polynomial,
    build_apply,
    kernel_basis,
    polynomial_multiply_all,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_component_mixing_probe.json"
INNER_NODES = 32_772
PACKET_SLOTS = 16
SAMPLE_SPACE = INNER_NODES * PACKET_SLOTS
FACTORS = {
    1: 0x3,
    2: 0x7,
    4: 0x13,
    9: 0x373,
    10: 0x519,
    18: 0x7C9C3,
    20: 0x1E1FFF,
}
DEFAULT_COMPONENT_SETS = (
    (1,),
    (2,),
    (4,),
    (9,),
    (10,),
    (18,),
    (20,),
    (18, 2),
    (18, 1),
    (10, 9, 1),
    (10, 4, 2, 1),
    (9, 4, 2, 1),
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def transpose_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(((columns[input_bit] >> output_bit) & 1) << input_bit for input_bit in range(64))
        for output_bit in range(64)
    )


def projection_tables(basis: tuple[int, ...]) -> np.ndarray:
    tables = np.zeros((8, 256), dtype=np.uint32)
    for byte_index in range(8):
        for byte in range(256):
            syndrome = 0
            value = byte << (8 * byte_index)
            for basis_index, vector in enumerate(basis):
                syndrome |= ((value & vector).bit_count() & 1) << basis_index
            tables[byte_index, byte] = syndrome
    return tables


def project(states: np.ndarray, tables: np.ndarray) -> np.ndarray:
    result = np.zeros(states.shape, dtype=np.uint32)
    for byte_index in range(8):
        indices = ((states >> np.uint64(8 * byte_index)) & np.uint64(0xFF)).astype(
            np.uint8
        )
        result ^= tables[byte_index, indices]
    return result


def walsh_transform(histogram: np.ndarray) -> np.ndarray:
    values = histogram.astype(np.int64, copy=True)
    half = 1
    while half < len(values):
        blocks = values.reshape(-1, 2 * half)
        left = blocks[:, :half].copy()
        right = blocks[:, half:].copy()
        blocks[:, :half] = left + right
        blocks[:, half:] = left - right
        half *= 2
    return values


def contribution_states(step, packet_value: int) -> np.ndarray:
    states = np.empty(SAMPLE_SPACE, dtype=np.uint64)
    cursor = 0
    for slot in range(PACKET_SLOTS):
        state = packet_value << (4 * slot)
        for _ in range(INNER_NODES):
            state = step(state)
            states[cursor] = state
            cursor += 1
    return states


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-packets", type=int, default=16)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.max_packets < 1:
        raise SystemExit("component mixing: max packets must be positive")

    p_columns = systematic_state_columns()
    apply_p = build_apply(p_columns)

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    t_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(t_columns))
    component_data = []
    for degree_set in DEFAULT_COMPONENT_SETS:
        annihilator = polynomial_multiply_all(FACTORS[degree] for degree in degree_set)
        columns = tuple(
            apply_polynomial(transpose_apply, annihilator, 1 << bit)
            for bit in range(64)
        )
        basis = kernel_basis(columns)
        dimension = sum(degree_set)
        if len(basis) != dimension:
            raise RuntimeError("component mixing: transpose kernel dimension mismatch")
        component_data.append(
            {
                "degrees": degree_set,
                "dimension": dimension,
                "basis": basis,
                "tables": projection_tables(basis),
                "value_rows": [],
            }
        )

    for packet_value in range(1, 16):
        states = contribution_states(step, packet_value)
        for component in component_data:
            syndromes = project(states, component["tables"])
            histogram = np.bincount(
                syndromes,
                minlength=1 << component["dimension"],
            )
            if int(histogram.sum()) != SAMPLE_SPACE:
                raise RuntimeError("component mixing: histogram mass mismatch")
            walsh = walsh_transform(histogram)
            if walsh[0] != SAMPLE_SPACE:
                raise RuntimeError("component mixing: trivial character mismatch")
            absolute = np.abs(walsh[1:])
            maximum_index = int(np.argmax(absolute)) + 1
            maximum_numerator = int(absolute[maximum_index - 1])
            packet_rows = []
            ratios = walsh[1:].astype(np.float64) / SAMPLE_SPACE
            squared = ratios * ratios
            powers = squared.copy()
            for packet_count in range(1, args.max_packets + 1):
                chi_squared = float(powers.sum())
                packet_rows.append(
                    {
                        "packets": packet_count,
                        "chi_squared": chi_squared,
                        "log2_chi_squared": (
                            math.log2(chi_squared) if chi_squared > 0.0 else None
                        ),
                        "total_variation_upper": 0.5 * math.sqrt(chi_squared),
                    }
                )
                powers *= squared
            component["value_rows"].append(
                {
                    "packet_value": packet_value,
                    "maximum_absolute_bias_numerator": maximum_numerator,
                    "maximum_absolute_bias_denominator": SAMPLE_SPACE,
                    "maximum_absolute_bias": maximum_numerator / SAMPLE_SPACE,
                    "maximizing_character_index": maximum_index,
                    "iid_repeated_value_rows": packet_rows,
                }
            )

    rows = []
    for component in component_data:
        value_rows = component["value_rows"]
        worst_bias = max(value_rows, key=lambda row: row["maximum_absolute_bias"])
        worst_tv_by_packets = []
        for packet_count in range(1, args.max_packets + 1):
            worst = max(
                value_rows,
                key=lambda row: row["iid_repeated_value_rows"][packet_count - 1][
                    "total_variation_upper"
                ],
            )
            metric = worst["iid_repeated_value_rows"][packet_count - 1]
            worst_tv_by_packets.append(
                {
                    "packets": packet_count,
                    "packet_value": worst["packet_value"],
                    "log2_chi_squared": metric["log2_chi_squared"],
                    "total_variation_upper": metric["total_variation_upper"],
                }
            )
        rows.append(
            {
                "component_degrees": list(component["degrees"]),
                "dimension": component["dimension"],
                "worst_single_packet_bias": worst_bias,
                "worst_repeated_value_mixing": worst_tv_by_packets,
                "value_rows": value_rows,
            }
        )

    payload = {
        "schema": "riffle-dp-g4-g2-component-mixing-probe-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "sample_space": (
            "One fixed nonzero nibble value, one of 16 slots, and one of 32772 "
            "node exponents, sampled independently and uniformly."
        ),
        "sample_space_size": SAMPLE_SPACE,
        "component_rows": rows,
        "scope_limitation": (
            "The Walsh spectra are exact. The mixing rows assume packet locations "
            "are sampled independently with replacement and repeat one value. "
            "They do not yet certify the without-replacement packet permutation."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    for row in rows:
        bias = row["worst_single_packet_bias"]
        final = row["worst_repeated_value_mixing"][-1]
        print(
            "degrees="
            + ",".join(map(str, row["component_degrees"]))
            + f" dimension={row['dimension']} worst_bias={bias['maximum_absolute_bias']:.9f}"
            + f" value={bias['packet_value']}"
            + f" tv_upper_at_{args.max_packets}={final['total_variation_upper']:.6e}"
        )
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_EXACT_COMPONENT_WALSH_MIXING")


if __name__ == "__main__":
    main()
