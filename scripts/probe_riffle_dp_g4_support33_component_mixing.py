#!/usr/bin/env python3
"""Retired support-33 component calculation with an invalid product step."""

from __future__ import annotations

import collections
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
from probe_riffle_dp_g4_component_mixing import (
    DEFAULT_COMPONENT_SETS,
    FACTORS,
    SAMPLE_SPACE,
    contribution_states,
    project,
    projection_tables,
    transpose_columns,
    walsh_transform,
)


ROOT = Path(__file__).resolve().parents[1]
OUTER_RECEIPTS = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_support33_component_mixing_probe.json"
PACKET_SUPPORT = 33
REMOVAL_CORRECTION = PACKET_SUPPORT - 1
FLOAT_EPSILON = np.finfo(np.float64).eps


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def packet_profile(receipt: dict) -> tuple[int, ...]:
    counts = collections.Counter()
    for local_word in receipt["local_words"]:
        word = int(local_word["codeword_hex"], 16)
        counts.update(
            (word >> (4 * packet)) & 0xF
            for packet in range(32)
            if ((word >> (4 * packet)) & 0xF) != 0
        )
    profile = tuple(counts[value] for value in range(1, 16))
    if sum(profile) != PACKET_SUPPORT:
        raise RuntimeError("support33 mixing: profile support mismatch")
    return profile


def directed_profile_bound(
    absolute_walsh: np.ndarray,
    profile: tuple[int, ...],
    denominator: int,
) -> tuple[float, float]:
    ratios = np.nextafter(
        (absolute_walsh.astype(np.float64) + REMOVAL_CORRECTION) / denominator,
        np.inf,
    )
    ratios = np.minimum(ratios, 1.0)
    character_bounds = np.ones(absolute_walsh.shape[1], dtype=np.float64)
    for value_index, count in enumerate(profile):
        for _ in range(count):
            character_bounds = np.nextafter(
                character_bounds * ratios[value_index],
                np.inf,
            )
    squared = np.nextafter(character_bounds * character_bounds, np.inf)
    summed = float(np.sum(squared, dtype=np.float64))
    additions = max(0, len(squared) - 1)
    if additions * FLOAT_EPSILON >= 1.0:
        raise RuntimeError("support33 mixing: summation error factor diverged")
    product = math.nextafter(additions * FLOAT_EPSILON, math.inf)
    gamma = math.nextafter(product / math.nextafter(1.0 - product, -math.inf), math.inf)
    one_minus_gamma = math.nextafter(1.0 - gamma, -math.inf)
    chi_squared_upper = math.nextafter(summed / one_minus_gamma, math.inf)
    tv_upper = math.nextafter(0.5 * math.sqrt(chi_squared_upper), math.inf)
    return chi_squared_upper, tv_upper


def main() -> None:
    outer = json.loads(OUTER_RECEIPTS.read_text())
    profile_families: dict[tuple[int, ...], list[str]] = collections.defaultdict(list)
    for receipt in outer["receipts"]:
        profile_families[packet_profile(receipt)].append(receipt["family_id"])
    if len(profile_families) != 14:
        raise RuntimeError("support33 mixing: unique profile count changed")

    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    t_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(t_columns))
    contribution_by_value = {
        value: contribution_states(step, value) for value in range(1, 16)
    }
    component_rows = []
    denominator = SAMPLE_SPACE - REMOVAL_CORRECTION
    for degree_set in DEFAULT_COMPONENT_SETS:
        annihilator = polynomial_multiply_all(FACTORS[degree] for degree in degree_set)
        columns = tuple(
            apply_polynomial(transpose_apply, annihilator, 1 << bit)
            for bit in range(64)
        )
        basis = kernel_basis(columns)
        dimension = sum(degree_set)
        if len(basis) != dimension:
            raise RuntimeError("support33 mixing: transpose component dimension mismatch")
        tables = projection_tables(basis)
        absolute_walsh = np.empty((15, (1 << dimension) - 1), dtype=np.int32)
        for value in range(1, 16):
            syndromes = project(contribution_by_value[value], tables)
            histogram = np.bincount(syndromes, minlength=1 << dimension)
            walsh = walsh_transform(histogram)
            absolute_walsh[value - 1] = np.abs(walsh[1:]).astype(np.int32)

        profile_rows = []
        for profile, family_ids in profile_families.items():
            chi_squared_upper, tv_upper = directed_profile_bound(
                absolute_walsh,
                profile,
                denominator,
            )
            profile_rows.append(
                {
                    "packet_value_counts_1_through_15": list(profile),
                    "family_ids": family_ids,
                    "chi_squared_upper": chi_squared_upper,
                    "log2_chi_squared_upper": math.nextafter(
                        math.log2(chi_squared_upper), math.inf
                    ),
                    "total_variation_upper": tv_upper,
                }
            )
        worst = max(profile_rows, key=lambda row: row["total_variation_upper"])
        component_rows.append(
            {
                "component_degrees": list(degree_set),
                "dimension": dimension,
                "worst_profile": worst,
                "profile_rows": profile_rows,
            }
        )

    payload = {
        "schema": "riffle-dp-g4-g2-support33-component-mixing-probe-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "INVALIDATED_BOUND",
        "outer_receipt_sha256": digest(OUTER_RECEIPTS),
        "source_sha256": digest(Path(__file__).resolve()),
        "unique_packet_profiles": len(profile_families),
        "invalidated_claim": (
            "The original calculation multiplied separate upper bounds on adaptive "
            "conditional character biases. Such bounds do not multiply unless the "
            "conditional expectations retain the required sign or martingale "
            "structure. They do not here."
        ),
        "rounding": (
            "Every ratio, product, square, chi-squared sum, and square root is "
            "rounded upward. The positive summation uses the standard gamma_n "
            "relative-error envelope with binary64 machine epsilon."
        ),
        "component_rows": component_rows,
        "scope_limitation": (
            "The exact Walsh spectra remain useful diagnostics. None of the "
            "reported chi-squared or total-variation values is a proved "
            "without-replacement bound. Use the iid-plus-distinctness reduction "
            "in bound_riffle_dp_2lap_g4_support33_projection.py instead."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"unique_support33_profiles={len(profile_families)}")
    for row in component_rows:
        worst = row["worst_profile"]
        print(
            "degrees="
            + ",".join(map(str, row["component_degrees"]))
            + f" dimension={row['dimension']}"
            + f" worst_tv_upper={worst['total_variation_upper']:.6e}"
            + f" log2_chi2={worst['log2_chi_squared_upper']:.6f}"
        )
    print(f"output={OUTPUT}")
    print("status=INVALIDATED_SUPPORT33_COMPONENT_MIXING")


if __name__ == "__main__":
    main()
