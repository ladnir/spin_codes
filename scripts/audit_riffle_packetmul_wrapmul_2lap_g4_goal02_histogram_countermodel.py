#!/usr/bin/env python3
"""Construct an exact countermodel to every histogram-only A6 argument."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RECEIPTS = CONSTRUCTION / "receipts"
A4 = RECEIPTS / "goal02_dual_a4_audit.json"
A6 = RECEIPTS / "goal02_a6_energy_interface.json"
OUTPUT = RECEIPTS / "goal02_pair_histogram_countermodel.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    a4 = json.loads(A4.read_text(encoding="utf-8"))
    a6 = json.loads(A6.read_text(encoding="utf-8"))
    histogram = {
        int(row["multiplicity"]): int(row["state_count"])
        for row in a4["pair_sum_multiplicity_histogram"]
    }
    assert histogram == {
        1: 2_199_542_387_068,
        2: 6_149_917,
        3: 1_114_128,
        4: 229_358,
        5: 32_762,
    }
    distinct = sum(histogram.values())
    mass = sum(multiplicity * count for multiplicity, count in histogram.items())
    energy = sum(
        multiplicity * multiplicity * count
        for multiplicity, count in histogram.items()
    )
    assert distinct == a4["distinct_nonzero_pair_sums"]
    assert mass == a6["pair_mass_l1"]
    assert energy == a6["pair_energy_l2_squared"]

    # Work inside a 42-dimensional subspace H of F_2^64.  Let K be a
    # 41-dimensional subspace of H.  Give every nonzero element of K
    # multiplicity one.  The rest of the authenticated histogram fits on
    # distinct points of H\K, so this defines a nonnegative function with the
    # same histogram, r(0)=0, and maximum multiplicity five.
    h_size = 1 << 42
    k_size = 1 << 41
    k_nonzero = k_size - 1
    remaining_distinct = distinct - k_nonzero
    remaining_multiplicity_one = histogram[1] - k_nonzero
    outside_capacity = h_size - k_size
    assert remaining_multiplicity_one >= 0
    assert remaining_distinct <= outside_capacity

    # For every ordered pair of distinct nonzero x,y in K, all of x, y, and
    # x+y have multiplicity one.  These terms alone give the lower bound.
    triangle_lower_bound = (k_size - 1) * (k_size - 2)
    triangle_cap = a6["sufficient_target"]["maximum_triangle_energy"]
    assert triangle_lower_bound > triangle_cap

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-pair-histogram-countermodel-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_COUNTERMODEL_TO_HISTOGRAM_ONLY_TRIANGLE_BOUND",
        "source_sha256": sha256(Path(__file__)),
        "dependency_sha256": {
            "dual_a4_audit": sha256(A4),
            "a6_energy_interface": sha256(A6),
        },
        "authenticated_histogram": [
            {"multiplicity": multiplicity, "state_count": count}
            for multiplicity, count in sorted(histogram.items())
        ],
        "authenticated_invariants_replayed": {
            "distinct_nonzero_support_points": distinct,
            "mass_l1": mass,
            "energy_l2_squared": energy,
            "maximum_multiplicity": max(histogram),
        },
        "countermodel": {
            "ambient_subspace_dimension": 42,
            "ambient_subspace_size": h_size,
            "planted_subspace_dimension": 41,
            "planted_subspace_size": k_size,
            "zero_multiplicity": 0,
            "multiplicity_on_nonzero_planted_subspace": 1,
            "remaining_distinct_histogram_points": remaining_distinct,
            "remaining_multiplicity_one_points": remaining_multiplicity_one,
            "available_distinct_points_outside_planted_subspace": outside_capacity,
            "placement_rule": "assign every remaining histogram entry to a distinct point of H minus K",
        },
        "triangle_energy": {
            "certified_lower_bound": triangle_lower_bound,
            "sufficient_actual_code_target": triangle_cap,
            "excess": triangle_lower_bound - triangle_cap,
            "ratio": {
                "numerator": triangle_lower_bound,
                "denominator": triangle_cap,
                "decimal": triangle_lower_bound / triangle_cap,
            },
            "counted_terms": "ordered distinct nonzero x,y in K, for which x+y is also nonzero in K",
        },
        "conclusion": "The exact multiplicity histogram, r(0)=0, L1 mass, L2 energy, and maximum multiplicity five do not imply the Goal 02 triangle-energy cap. A proof must use the actual pair-sum locations or an equivalent response-code invariant.",
        "scope_limit": "The countermodel is not the actual response pair-sum function and does not refute the construction.",
        "result": "PASS_EXACT_HISTOGRAM_ONLY_OBSTRUCTION",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")
    print(f"triangle_lower_bound={triangle_lower_bound}")
    print(f"triangle_cap={triangle_cap}")
    print(f"ratio={triangle_lower_bound / triangle_cap:.12f}")
    print("status=PASS_EXACT_HISTOGRAM_ONLY_OBSTRUCTION")


if __name__ == "__main__":
    main()
