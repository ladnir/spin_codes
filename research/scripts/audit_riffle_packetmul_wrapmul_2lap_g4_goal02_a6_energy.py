#!/usr/bin/env python3
"""Derive the exact Goal 02 A6 triangle-energy interface."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = (
    ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4" / "receipts"
)
A4_AUDIT = RECEIPTS / "goal02_dual_a4_audit.json"
OUTPUT = RECEIPTS / "goal02_a6_energy_interface.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay_small_identity(values: tuple[int, ...]) -> dict:
    pair_counts = Counter(left ^ right for left, right in combinations(values, 2))
    triangle = sum(
        pair_counts[left] * pair_counts[right] * pair_counts[left ^ right]
        for left in pair_counts
        for right in pair_counts
    )
    a4 = sum(
        1
        for subset in combinations(values, 4)
        if subset[0] ^ subset[1] ^ subset[2] ^ subset[3] == 0
    )
    a6 = sum(
        1
        for subset in combinations(values, 6)
        if subset[0] ^ subset[1] ^ subset[2] ^ subset[3] ^ subset[4] ^ subset[5]
        == 0
    )
    reconstructed = 6 * comb(len(values), 3) + (36 * len(values) - 120) * a4 + 90 * a6
    if triangle != reconstructed:
        raise RuntimeError("A6 energy: small-model identity replay failed")
    return {
        "set_size": len(values),
        "triangle_energy": triangle,
        "dual_weight4_count": a4,
        "dual_weight6_count": a6,
        "identity_reconstructed_triangle_energy": reconstructed,
    }


def main() -> None:
    a4_audit = json.loads(A4_AUDIT.read_text())
    n = a4_audit["response_length"]
    a4 = a4_audit["dual_weight4_word_count"]
    a6_cap = a4_audit["sixth_moment_update"][
        "maximum_dual_weight6_words_sufficient_for_goal02"
    ]
    histogram = {
        int(row["multiplicity"]): int(row["state_count"])
        for row in a4_audit["pair_sum_multiplicity_histogram"]
    }
    pair_mass = sum(multiplicity * count for multiplicity, count in histogram.items())
    pair_l2_squared = sum(
        multiplicity * multiplicity * count
        for multiplicity, count in histogram.items()
    )
    collisions = sum(
        comb(multiplicity, 2) * count
        for multiplicity, count in histogram.items()
    )
    if pair_mass != comb(n, 2):
        raise RuntimeError("A6 energy: pair mass mismatch")
    if pair_l2_squared != pair_mass + 2 * collisions:
        raise RuntimeError("A6 energy: pair L2 identity mismatch")
    if collisions != 3 * a4:
        raise RuntimeError("A6 energy: A4 identity mismatch")

    triangle_baseline = 6 * comb(n, 3)
    triangle_a4_coefficient = 36 * n - 120
    sufficient_triangle_cap = (
        triangle_baseline + triangle_a4_coefficient * a4 + 90 * a6_cap
    )
    young_triangle_bound = pair_mass * pair_l2_squared
    young_a6_bound = (
        young_triangle_bound - triangle_baseline - triangle_a4_coefficient * a4
    ) // 90
    improvement = Fraction(young_triangle_bound, sufficient_triangle_cap)
    if young_triangle_bound <= sufficient_triangle_cap:
        raise RuntimeError("A6 energy: generic bound unexpectedly closes the target")
    small_replays = [
        replay_small_identity(tuple(range(1, 17))),
        replay_small_identity(tuple((index * index + 3 * index + 1) % 127 for index in range(18))),
        replay_small_identity(tuple((5 * index + 7) % 251 for index in range(20))),
    ]
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-a6-energy-v1",
        "candidate": a4_audit["candidate"],
        "evidence_label": "EXACT_TRIANGLE_ENERGY_REDUCTION",
        "source_sha256": digest(Path(__file__).resolve()),
        "dual_a4_audit_sha256": digest(A4_AUDIT),
        "response_coordinate_count": n,
        "pair_sum_definition": (
            "r(x) is the number of unordered pairs of distinct response "
            "coordinate columns whose XOR equals x"
        ),
        "pair_mass_l1": pair_mass,
        "pair_energy_l2_squared": pair_l2_squared,
        "maximum_pair_sum_multiplicity": max(histogram),
        "triangle_energy_definition": "T(r)=sum_{x,y} r(x) r(y) r(x+y)",
        "exact_dual_identity": {
            "formula": "T(r)=6*C(n,3)+(36*n-120)*A4_dual+90*A6_dual",
            "triangle_baseline": triangle_baseline,
            "weight4_coefficient": triangle_a4_coefficient,
            "weight4_count": a4,
            "weight6_coefficient": 90,
            "derivation": (
                "An ordered triple of simple pair edges has odd endpoint "
                "support of size 0, 4, or 6. Triangles give 6*C(n,3); each "
                "weight-4 dual support gives 36*n-120 V-plus-edge/star "
                "configurations; each weight-6 support gives 15 perfect "
                "matchings times 3! edge orders."
            ),
        },
        "small_model_identity_replays": small_replays,
        "sufficient_target": {
            "maximum_dual_weight6_words": a6_cap,
            "maximum_triangle_energy": sufficient_triangle_cap,
        },
        "multiplicity_only_bound": {
            "inequality": "T(r)<=||r||_1*||r||_2^2",
            "reason": "Cauchy-Schwarz followed by Young convolution inequality",
            "triangle_energy_upper_bound": young_triangle_bound,
            "implied_dual_weight6_upper_bound": young_a6_bound,
            "factor_above_sufficient_triangle_cap": {
                "numerator": improvement.numerator,
                "denominator": improvement.denominator,
                "decimal": float(improvement),
                "ceiling": (improvement.numerator + improvement.denominator - 1)
                // improvement.denominator,
            },
        },
        "remaining_obligation": (
            "Use the locations of the pair sums, beyond their authenticated "
            "multiplicity histogram, to prove T(r) does not exceed the "
            "sufficient triangle cap; or prove the alternative C8 distance."
        ),
        "result": "EXACT_INTERFACE; LOCATION-SENSITIVE_A6_BOUND_REQUIRED",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"pair_l2_squared={pair_l2_squared}")
    print(f"sufficient_triangle_cap={sufficient_triangle_cap}")
    print(f"young_triangle_bound={young_triangle_bound}")
    print(f"gap_factor={float(improvement):.12f}")
    print("status=EXACT_LOCATION_SENSITIVE_A6_INTERFACE")


if __name__ == "__main__":
    main()
