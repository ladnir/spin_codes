#!/usr/bin/env python3
"""Rigorous support-33 terminal-zero bound on the degree-20 projection."""

from __future__ import annotations

import collections
import hashlib
import json
import math
from fractions import Fraction
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
    project,
    projection_tables,
    transpose_columns,
    walsh_transform,
)


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
OUTER = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_support33_projection_bound.json"
PACKET_SUPPORT = 33
PROJECTION_DEGREE = 20
FLOAT_EPSILON = np.finfo(np.float64).eps


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def packet_profile(receipt: dict) -> tuple[int, ...]:
    counts: collections.Counter[int] = collections.Counter()
    for local_word in receipt["local_words"]:
        word = int(local_word["codeword_hex"], 16)
        counts.update(
            (word >> shift) & 15
            for shift in range(0, 128, 4)
            if (word >> shift) & 15
        )
    result = tuple(counts[value] for value in range(1, 16))
    if sum(result) != PACKET_SUPPORT:
        raise RuntimeError("projection bound: support-33 profile changed")
    return result


def upward_sum(values: np.ndarray) -> float:
    raw = float(np.sum(values, dtype=np.float64))
    additions = max(0, len(values) - 1)
    if additions * FLOAT_EPSILON >= 1.0:
        raise RuntimeError("projection bound: summation envelope diverged")
    gamma_product = math.nextafter(additions * FLOAT_EPSILON, math.inf)
    gamma = math.nextafter(
        gamma_product / math.nextafter(1.0 - gamma_product, -math.inf),
        math.inf,
    )
    return math.nextafter(raw / math.nextafter(1.0 - gamma, -math.inf), math.inf)


def main() -> None:
    outer = json.loads(OUTER.read_text())
    profiles: dict[tuple[int, ...], list[str]] = collections.defaultdict(list)
    for receipt in outer["receipts"]:
        profiles[packet_profile(receipt)].append(receipt["family_id"])
    if len(profiles) != 14:
        raise RuntimeError("projection bound: profile count changed")

    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    t_columns = tuple(step(1 << bit) for bit in range(64))
    transpose_apply = build_apply(transpose_columns(t_columns))
    component_columns = tuple(
        apply_polynomial(
            transpose_apply,
            FACTORS[PROJECTION_DEGREE],
            1 << bit,
        )
        for bit in range(64)
    )
    basis = kernel_basis(component_columns)
    if len(basis) != PROJECTION_DEGREE:
        raise RuntimeError("projection bound: component dimension mismatch")
    tables = projection_tables(basis)

    absolute_walsh = np.empty(
        (15, (1 << PROJECTION_DEGREE) - 1),
        dtype=np.int32,
    )
    for value in range(1, 16):
        states = contribution_states(step, value)
        syndromes = project(states, tables)
        histogram = np.bincount(syndromes, minlength=1 << PROJECTION_DEGREE)
        walsh = walsh_transform(histogram)
        if int(walsh[0]) != SAMPLE_SPACE:
            raise RuntimeError("projection bound: Walsh mass mismatch")
        absolute_walsh[value - 1] = np.abs(walsh[1:]).astype(np.int32)

    distinct_probability = Fraction(
        math.prod(range(SAMPLE_SPACE - PACKET_SUPPORT + 1, SAMPLE_SPACE + 1)),
        SAMPLE_SPACE**PACKET_SUPPORT,
    )
    distinct_lower = math.nextafter(float(distinct_probability), -math.inf)
    ratio_rows = np.nextafter(
        absolute_walsh.astype(np.float64) / SAMPLE_SPACE,
        np.inf,
    )

    rows = []
    aggregate_upper = 0.0
    for profile, family_ids in profiles.items():
        character_terms = np.ones(ratio_rows.shape[1], dtype=np.float64)
        for value_index, count in enumerate(profile):
            for _ in range(count):
                character_terms = np.nextafter(
                    character_terms * ratio_rows[value_index],
                    np.inf,
                )
        nontrivial_sum_upper = upward_sum(character_terms)
        iid_zero_upper = math.nextafter(
            math.ldexp(math.nextafter(1.0 + nontrivial_sum_upper, math.inf),
                       -PROJECTION_DEGREE),
            math.inf,
        )
        conditioned_upper = min(
            1.0,
            math.nextafter(iid_zero_upper / distinct_lower, math.inf),
        )
        aggregate_upper = math.nextafter(
            aggregate_upper + len(family_ids) * conditioned_upper,
            math.inf,
        )
        rows.append(
            {
                "packet_value_counts_1_through_15": list(profile),
                "family_ids": family_ids,
                "nontrivial_fourier_l1_upper": nontrivial_sum_upper,
                "iid_projection_zero_upper": iid_zero_upper,
                "conditioned_projection_zero_upper": conditioned_upper,
                "conditioned_projection_zero_log2_upper": math.nextafter(
                    math.log2(conditioned_upper), math.inf
                ),
            }
        )

    worst = max(rows, key=lambda row: row["conditioned_projection_zero_upper"])
    payload = {
        "schema": "riffle-dp-2lap-g4-support33-projection-bound-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "RIGOROUS_BOUND",
        "outer_receipt_sha256": digest(OUTER),
        "projection_component_degree": PROJECTION_DEGREE,
        "projection_dimension": PROJECTION_DEGREE,
        "packet_positions": SAMPLE_SPACE,
        "packet_support": PACKET_SUPPORT,
        "distinctness_probability_numerator": str(distinct_probability.numerator),
        "distinctness_probability_denominator": str(distinct_probability.denominator),
        "reduction": (
            "First sample the 33 labeled packet cells independently and uniformly. "
            "Conditioned on all cells being distinct, their law is the uniform "
            "ordered injection induced by the global packet permutation. Hence a "
            "point probability under the permutation is at most its iid value "
            "divided by the exact distinctness probability."
        ),
        "fourier_bound": (
            "For iid cells, Fourier inversion expresses the zero-syndrome "
            "probability as 2^-20 times the sum of the profile products of "
            "character biases. The calculation takes absolute values only for "
            "nontrivial characters and rounds every operation upward."
        ),
        "rounding": (
            "Bias division, profile products, summation, Fourier normalization, "
            "and conditioning division are rounded upward. The positive sum uses "
            "the standard gamma_n binary64 error envelope. The distinctness "
            "probability is rounded downward before division."
        ),
        "profile_rows": rows,
        "worst_profile": worst,
        "aggregate_over_26_words_upper": aggregate_upper,
        "aggregate_over_26_words_log2_upper": math.nextafter(
            math.log2(aggregate_upper), math.inf
        ),
        "scope_limitation": (
            "The bound controls only the event that the first-lap terminal state "
            "is zero through its degree-20 projection. It does not bound output "
            "weight or cover outer words beyond the 26 authenticated support-33 "
            "receipts."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"unique_support33_profiles={len(rows)}")
    print(
        "worst_projection_zero_log2_upper="
        f"{worst['conditioned_projection_zero_log2_upper']:.12f}"
    )
    print(
        "aggregate_26_words_log2_upper="
        f"{payload['aggregate_over_26_words_log2_upper']:.12f}"
    )
    print(f"output={OUTPUT}")
    print("status=RIGOROUS_SUPPORT33_PROJECTION_BOUND")


if __name__ == "__main__":
    main()
