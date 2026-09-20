#!/usr/bin/env python3
"""Exhaust exact 12-node character distances in selected component subspaces."""

from __future__ import annotations

import hashlib
import itertools
import json
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
from probe_riffle_dp_g4_component_mixing import (
    FACTORS,
    project,
    projection_tables,
    transpose_columns,
    walsh_transform,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "explorations" / "riffle_dp_2lap_g4_local_character_components.json"
WINDOW_NODES = 12
SLOTS = 16
WINDOW_BITS = 192
KNOWN_CHARACTER = 0xA685AAC60E5ACFB3
COMPONENT_DIMENSION_CAP = 22


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def observation_columns(step, packet_value: int) -> list[int]:
    columns = []
    for slot in range(SLOTS):
        state = packet_value << (4 * slot)
        for _ in range(WINDOW_NODES):
            state = step(state)
            columns.append(state)
    return columns


def main() -> None:
    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    step_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(step_columns))
    observed = {
        value: np.array(observation_columns(step, value), dtype=np.uint64)
        for value in range(1, 16)
    }

    all_degrees = tuple(sorted(FACTORS))
    supports_through_cap = tuple(
        support
        for size in range(1, len(all_degrees) + 1)
        for support in itertools.combinations(all_degrees, size)
        if sum(support) <= COMPONENT_DIMENSION_CAP
    )
    maximal_supports = tuple(
        support
        for support in supports_through_cap
        if not any(set(support) < set(other) for other in supports_through_cap)
    )
    component_sets = tuple(dict.fromkeys(
        tuple((degree,) for degree in all_degrees) + maximal_supports
    ))

    rows = []
    total_character_evaluations = 0
    for degree_set in component_sets:
        annihilator = polynomial_multiply_all(FACTORS[degree] for degree in degree_set)
        component_columns = tuple(
            apply_polynomial(transpose_apply, annihilator, 1 << bit)
            for bit in range(64)
        )
        basis = kernel_basis(component_columns)
        dimension = sum(degree_set)
        if len(basis) != dimension:
            raise RuntimeError("local component certificate: dimension mismatch")
        tables = projection_tables(basis)
        value_rows = []
        for value in range(1, 16):
            syndromes = project(observed[value], tables)
            histogram = np.bincount(syndromes, minlength=1 << dimension)
            walsh = walsh_transform(histogram)
            if int(walsh[0]) != WINDOW_BITS:
                raise RuntimeError("local component certificate: Walsh mass mismatch")
            absolute = np.abs(walsh[1:])
            maximum_index = int(np.argmax(absolute)) + 1
            maximum_sum = int(absolute[maximum_index - 1])
            two_sided_distance = (WINDOW_BITS - maximum_sum) // 2
            required = 36 if value == 15 else 48
            if two_sided_distance < required:
                raise RuntimeError("local component certificate: target refuted")
            value_rows.append(
                {
                    "packet_value": value,
                    "maximum_absolute_character_sum": maximum_sum,
                    "two_sided_distance": two_sided_distance,
                    "required_distance": required,
                    "passed": True,
                    "maximizing_character_index": maximum_index,
                }
            )
        total_character_evaluations += 15 * ((1 << dimension) - 1)
        rows.append(
            {
                "component_degrees": list(degree_set),
                "dimension": dimension,
                "nonzero_characters_per_value": (1 << dimension) - 1,
                "minimum_two_sided_distance": min(
                    row["two_sided_distance"] for row in value_rows
                ),
                "value_rows": value_rows,
            }
        )

    covered_degree_sets = tuple(set(degrees) for degrees in component_sets)
    if any(
        not any(set(support) <= covered for covered in covered_degree_sets)
        for support in supports_through_cap
    ):
        raise RuntimeError("local component certificate: support coverage gap")

    known_weight = sum(
        (KNOWN_CHARACTER & column).bit_count() & 1
        for column in observation_columns(step, 15)
    )
    if known_weight != 36:
        raise RuntimeError("local component certificate: tight witness changed")
    payload = {
        "schema": "riffle-dp-2lap-g4-local-character-components-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_PARTIAL_CERTIFICATE",
        "source_sha256": digest(Path(__file__).resolve()),
        "window_nodes": 12,
        "window_bits": WINDOW_BITS,
        "component_rows": rows,
        "component_support_coverage": (
            "Every irreducible-component support whose total dimension is at "
            f"most {COMPONENT_DIMENSION_CAP} is contained in a listed exhausted subspace."
        ),
        "component_dimension_cap": COMPONENT_DIMENSION_CAP,
        "covered_support_count": len(supports_through_cap),
        "character_value_evaluations_including_overlaps": total_character_evaluations,
        "known_tight_character_hex": hex(KNOWN_CHARACTER),
        "known_value15_weight": known_weight,
        "scope_limitation": (
            "Every listed component subspace is exhausted exactly. The union does "
            "not contain characters whose irreducible-component support is absent "
            "from the listed degree sets. Therefore this artifact is not yet a "
            "full 64-dimensional distance certificate."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for row in rows:
        print(
            f"degrees={row['component_degrees']} dimension={row['dimension']} "
            f"minimum_two_sided_distance={row['minimum_two_sided_distance']}"
        )
    print(f"character_value_evaluations={total_character_evaluations}")
    print(f"output={OUTPUT}")
    print("status=EXACT_PARTIAL_COMPONENT_CERTIFICATE")


if __name__ == "__main__":
    main()
