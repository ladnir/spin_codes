#!/usr/bin/env python3
"""Independent arithmetic replay of the Goal 02 histogram countermodel."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4" / "receipts"
MODEL = RECEIPTS / "goal02_pair_histogram_countermodel.json"
A4 = RECEIPTS / "goal02_dual_a4_audit.json"
A6 = RECEIPTS / "goal02_a6_energy_interface.json"
MODEL_SOURCE = ROOT / "scripts" / "audit_riffle_packetmul_wrapmul_2lap_g4_goal02_histogram_countermodel.py"
OUTPUT = RECEIPTS / "goal02_pair_histogram_countermodel_independent.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    model = json.loads(MODEL.read_text(encoding="utf-8"))
    a4 = json.loads(A4.read_text(encoding="utf-8"))
    a6 = json.loads(A6.read_text(encoding="utf-8"))
    assert model["source_sha256"] == sha256(MODEL_SOURCE)
    assert model["dependency_sha256"] == {
        "dual_a4_audit": sha256(A4),
        "a6_energy_interface": sha256(A6),
    }
    assert model["result"] == "PASS_EXACT_HISTOGRAM_ONLY_OBSTRUCTION"

    histogram = tuple(
        (int(row["multiplicity"]), int(row["state_count"]))
        for row in a4["pair_sum_multiplicity_histogram"]
    )
    distinct = sum(count for _, count in histogram)
    mass = sum(multiplicity * count for multiplicity, count in histogram)
    square_mass = sum(multiplicity**2 * count for multiplicity, count in histogram)
    assert distinct == a4["distinct_nonzero_pair_sums"]
    assert mass == a6["pair_mass_l1"]
    assert square_mass == a6["pair_energy_l2_squared"]

    planted_size = 2**41
    ambient_size = 2**42
    multiplicity_one = dict(histogram)[1]
    used_in_planted = planted_size - 1
    outside_needed = distinct - used_in_planted
    outside_capacity = ambient_size - planted_size
    assert multiplicity_one >= used_in_planted
    assert 0 <= outside_needed <= outside_capacity

    # This is a direct ordered-pair recount rather than a replay of the
    # primary formula: x has q-1 nonzero choices, and y then has q-2 choices
    # after excluding zero and x.  Characteristic two makes x+y nonzero.
    ordered_x_choices = planted_size - 1
    ordered_y_choices = planted_size - 2
    lower_bound = ordered_x_choices * ordered_y_choices
    cap = a6["sufficient_target"]["maximum_triangle_energy"]
    ratio = Fraction(lower_bound, cap)
    assert lower_bound > cap
    assert model["triangle_energy"]["certified_lower_bound"] == lower_bound
    assert model["triangle_energy"]["excess"] == lower_bound - cap

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-pair-histogram-countermodel-independent-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "INDEPENDENT_EXACT_HISTOGRAM_COUNTERMODEL_REPLAY",
        "source_sha256": sha256(Path(__file__)),
        "model_sha256": sha256(MODEL),
        "dependency_sha256": model["dependency_sha256"],
        "replayed_invariants": {
            "distinct_nonzero_support_points": distinct,
            "mass_l1": mass,
            "energy_l2_squared": square_mass,
            "outside_points_needed": outside_needed,
            "outside_capacity": outside_capacity,
            "triangle_lower_bound": lower_bound,
            "triangle_cap": cap,
            "ratio_numerator": ratio.numerator,
            "ratio_denominator": ratio.denominator,
        },
        "scope_limit": "This proves insufficiency of histogram-only data, not failure of the actual response code.",
        "result": "PASS_INDEPENDENT_HISTOGRAM_OBSTRUCTION_REPLAY",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")
    print(f"triangle_lower_bound={lower_bound}")
    print(f"triangle_cap={cap}")
    print("status=PASS_INDEPENDENT_HISTOGRAM_OBSTRUCTION_REPLAY")


if __name__ == "__main__":
    main()
