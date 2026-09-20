#!/usr/bin/env python3
"""Authenticate the primary and independent Goal 02 dual-A4 censuses."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RECEIPTS = CANDIDATE / "receipts"
PRIMARY_SOURCE = (
    ROOT / "scripts" / "certify_riffle_packetmul_wrapmul_2lap_g4_goal02_dual_a4.cpp"
)
INDEPENDENT_SOURCE = (
    ROOT / "scripts" / "verify_riffle_packetmul_wrapmul_2lap_g4_goal02_dual_a4.cpp"
)
PRIMARY_EXE = PRIMARY_SOURCE.with_suffix(".exe")
INDEPENDENT_EXE = INDEPENDENT_SOURCE.with_suffix(".exe")
PRIMARY = RECEIPTS / "goal02_dual_a4_primary_raw.json"
INDEPENDENT = RECEIPTS / "goal02_dual_a4_independent_raw.json"
MOMENT = RECEIPTS / "goal02_sixth_moment_primary.json"
OUTPUT = RECEIPTS / "goal02_dual_a4_audit.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    primary = json.loads(PRIMARY.read_text())
    independent = json.loads(INDEPENDENT.read_text())
    moment = json.loads(MOMENT.read_text())
    exact_fields = (
        "response_nodes",
        "response_coordinates",
        "irreducible_component_degrees",
        "irreducible_component_orders",
        "total_autonomous_orbit_keys",
        "normalized_unordered_pair_types",
        "occupied_pair_sum_orbit_keys",
        "maximum_normalized_pair_types_per_orbit_key",
        "distinct_nonzero_pair_sums",
        "maximum_pair_sum_multiplicity",
        "pair_sum_multiplicity_histogram",
        "equal_pair_sum_collisions",
        "dual_weight4_word_count",
        "collision_to_word_divisor",
        "result",
    )
    for field in exact_fields:
        if primary[field] != independent[field]:
            raise RuntimeError(f"dual A4 audit: primary and independent differ: {field}")
    if primary["result"] != "PASS":
        raise RuntimeError("dual A4 audit: census did not pass")
    records = 64 * 63 // 2 + (32_772 - 1) * 64 * 64
    if primary["normalized_unordered_pair_types"] != records:
        raise RuntimeError("dual A4 audit: normalized pair schedule changed")
    collisions = int(primary["equal_pair_sum_collisions"])
    dual_a4 = int(primary["dual_weight4_word_count"])
    if collisions != 3 * dual_a4 or dual_a4 != 3_732_023:
        raise RuntimeError("dual A4 audit: collision-to-word arithmetic changed")
    histogram = {
        int(row["multiplicity"]): int(row["state_count"])
        for row in primary["pair_sum_multiplicity_histogram"]
    }
    expected_histogram = {
        1: 2_199_542_387_068,
        2: 6_149_917,
        3: 1_114_128,
        4: 229_358,
        5: 32_762,
    }
    if histogram != expected_histogram:
        raise RuntimeError("dual A4 audit: pair-sum multiplicity histogram changed")
    pair_mass = sum(multiplicity * count for multiplicity, count in histogram.items())
    response_coordinates = primary["response_coordinates"]
    if pair_mass != response_coordinates * (response_coordinates - 1) // 2:
        raise RuntimeError("dual A4 audit: pair-sum histogram has wrong mass")
    if sum(count * (multiplicity * (multiplicity - 1) // 2)
           for multiplicity, count in histogram.items()) != collisions:
        raise RuntimeError("dual A4 audit: histogram collision count changed")

    length = moment["response_length"]
    coefficient4 = moment["sixth_moment_identity"]["weight4_dual_coefficient"]
    scaled_budget = moment["markov_reduction"]["scaled_joint_sparse_budget"]
    scaled_a4_charge = (1 << 64) * coefficient4 * dual_a4
    remaining_scaled_budget = scaled_budget - scaled_a4_charge
    if remaining_scaled_budget <= 0:
        raise RuntimeError("dual A4 audit: A4 consumes the moment budget")
    dual_a6_cap = remaining_scaled_budget // ((1 << 64) * 720)
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-dual-a4-audit-v1",
        "candidate": primary["candidate"],
        "evidence_label": "AUTHENTICATED_PRIMARY_AND_INDEPENDENT_EXACT_CENSUS",
        "source_sha256": digest(Path(__file__).resolve()),
        "primary_source_sha256": digest(PRIMARY_SOURCE),
        "independent_source_sha256": digest(INDEPENDENT_SOURCE),
        "primary_executable_sha256": digest(PRIMARY_EXE),
        "independent_executable_sha256": digest(INDEPENDENT_EXE),
        "primary_raw_receipt_sha256": digest(PRIMARY),
        "independent_raw_receipt_sha256": digest(INDEPENDENT),
        "sixth_moment_receipt_sha256": digest(MOMENT),
        "response_length": length,
        "normalized_unordered_pair_types": records,
        "occupied_pair_sum_orbit_keys": primary["occupied_pair_sum_orbit_keys"],
        "maximum_normalized_pair_types_per_orbit_key": primary[
            "maximum_normalized_pair_types_per_orbit_key"
        ],
        "distinct_nonzero_pair_sums": int(primary["distinct_nonzero_pair_sums"]),
        "maximum_pair_sum_multiplicity": primary["maximum_pair_sum_multiplicity"],
        "pair_sum_multiplicity_histogram": primary["pair_sum_multiplicity_histogram"],
        "equal_pair_sum_collisions": collisions,
        "dual_weight4_word_count": dual_a4,
        "independent_replay_changes": {
            "component_order": "ascending instead of descending",
            "krylov_seed": "last kernel vector instead of first",
            "collision_accumulation": "pairwise cyclic-arc intersections instead of sweep line",
        },
        "sixth_moment_update": {
            "weight4_dual_coefficient": coefficient4,
            "scaled_a4_charge": scaled_a4_charge,
            "remaining_scaled_sparse_budget": remaining_scaled_budget,
            "maximum_dual_weight6_words_sufficient_for_goal02": dual_a6_cap,
            "remaining_sufficient_inequality": f"A6_dual <= {dual_a6_cap}",
        },
        "result": "PASS; EXACT_A4_CLOSED; A6_REQUIRED",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"dual_A4={dual_a4}")
    print(f"dual_A6_cap={dual_a6_cap}")
    print("status=PASS_EXACT_A4_CLOSED")


if __name__ == "__main__":
    main()
