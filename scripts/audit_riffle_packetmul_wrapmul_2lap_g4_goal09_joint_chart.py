#!/usr/bin/env python3
"""Audit the Goal 09 width-16 joint-chart decoder experiments."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from audit_riffle_packetmul_wrapmul_2lap_bch_family import build_instance, observe


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RECEIPTS = CANDIDATE / "receipts"
CPP_SOURCE = ROOT / "scripts" / "analyze_riffle_packetmul_wrapmul_2lap_g4_goal09_chart_exact_width16.cpp"
CPP_EXE = ROOT / "scripts" / "analyze_riffle_packetmul_wrapmul_2lap_g4_goal09_chart_exact_width16.exe"
SAT_SOURCE = ROOT / "scripts" / "probe_riffle_packetmul_wrapmul_2lap_g4_goal09_chart_sat.py"
RAW = RECEIPTS / "goal09_chart_exact_width16_raw.json"
SAT_121 = RECEIPTS / "goal09_chart_sat_width16_c24_anchor_01_set_01_bound_121.json"
SAT_120 = RECEIPTS / "goal09_chart_sat_width16_c24_anchor_01_set_01_bound_120.json"
GOAL07_RAW = RECEIPTS / "goal07_bch32_family_exact_raw.json"
GOAL08_AUDIT = RECEIPTS / "goal08_c38_preflight_audit.json"
OUTPUT = RECEIPTS / "goal09_joint_chart_audit.json"
WITNESS = 0xCD7E70D6


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for original in rows:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def main() -> None:
    raw = json.loads(RAW.read_text())
    sat121 = json.loads(SAT_121.read_text())
    sat120 = json.loads(SAT_120.read_text())
    goal07 = json.loads(GOAL07_RAW.read_text())
    goal08 = json.loads(GOAL08_AUDIT.read_text())
    if raw["result"] != "EXHAUSTED":
        raise RuntimeError("Goal 09 audit: chart enumeration was not exhausted")
    if raw["chart"] != {
        "anchor_node": 1,
        "restriction_node": 2,
        "rank": 32,
        "group_weight_cap": 5,
    }:
        raise RuntimeError("Goal 09 audit: chart metadata changed")
    list_size = sum(math.comb(16, weight) for weight in range(6))
    if raw["list_size_each"] != list_size:
        raise RuntimeError("Goal 09 audit: list size mismatch")
    if raw["cartesian_pairs"] != list_size**2:
        raise RuntimeError("Goal 09 audit: Cartesian pair count mismatch")
    if raw["nonzero_pairs"] != list_size**2 - 1:
        raise RuntimeError("Goal 09 audit: nonzero pair count mismatch")

    instance = build_instance(5, 7)
    basis_outputs = [observe(1 << bit, 16, instance["columns"], 24) for bit in range(32)]
    chart_forms = [
        sum(((outputs[node] >> coordinate) & 1) << bit for bit, outputs in enumerate(basis_outputs))
        for node in (1, 2)
        for coordinate in range(16)
    ]
    if rank(chart_forms) != 32:
        raise RuntimeError("Goal 09 audit: q1,q2 chart is not bijective")
    witness_outputs = observe(WITNESS, 16, instance["columns"], 24)
    witness_weights = [value.bit_count() for value in witness_outputs]
    if sum(witness_weights) != 121 or witness_weights[1:3] != [4, 4]:
        raise RuntimeError("Goal 09 audit: width-16 witness replay failed")
    if raw["minimum_weight"] != 121 or raw["minimum_count"] != 1:
        raise RuntimeError("Goal 09 audit: chart minimum changed")
    if int(raw["witness_hex"], 16) != WITNESS:
        raise RuntimeError("Goal 09 audit: chart witness changed")
    exact_global = goal07["prefix_minima"][-1]
    if (
        exact_global["weight"] != raw["minimum_weight"]
        or exact_global["count"] != raw["minimum_count"]
        or int(exact_global["witness_hex"], 16) != WITNESS
    ):
        raise RuntimeError("Goal 09 audit: chart and full census disagree")

    for receipt, bound in ((sat121, 121), (sat120, 120)):
        if receipt["source_sha256"] != digest(SAT_SOURCE):
            raise RuntimeError("Goal 09 audit: SAT source hash changed")
        if receipt["result"] != "UNKNOWN":
            raise RuntimeError("Goal 09 audit: SAT diagnostic result changed")
        if receipt["decision_problem"]["total_output_weight_upper_bound"] != bound:
            raise RuntimeError("Goal 09 audit: SAT bound mismatch")
        if receipt["chart"]["rank"] != 32:
            raise RuntimeError("Goal 09 audit: SAT chart rank mismatch")

    width64_left = sum(math.comb(64, weight) for weight in range(6))
    width64_right = sum(math.comb(64, weight) for weight in range(7))
    width64_pairs = width64_left * width64_right
    chart_count = goal08["uniform_c38_all_anchors"]["audited_information_set_count"]
    pairs_per_second = raw["cartesian_pairs"] / raw["elapsed_seconds"]
    optimistic_seconds_per_chart = width64_pairs / pairs_per_second
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal09-joint-chart-audit-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "AUDITED_EXACT_SMALL_CHART_AND_SCALING_OBSTRUCTION",
        "source_sha256": digest(Path(__file__).resolve()),
        "exact_decoder": {
            "source_sha256": digest(CPP_SOURCE),
            "executable_sha256": digest(CPP_EXE),
            "raw_receipt_sha256": digest(RAW),
            "chart": "Phi(s)=(q_1(s),q_2(s))",
            "chart_rank": 32,
            "left_weight_cap": 5,
            "right_weight_cap": 5,
            "list_size_each": list_size,
            "cartesian_pairs": raw["cartesian_pairs"],
            "elapsed_seconds": raw["elapsed_seconds"],
            "minimum_weight": raw["minimum_weight"],
            "minimum_count": raw["minimum_count"],
            "witness_hex": raw["witness_hex"],
            "matches_full_width16_census": True,
        },
        "native_xor_sat_diagnostics": [
            {
                "weight_bound": receipt["decision_problem"]["total_output_weight_upper_bound"],
                "receipt_sha256": digest(path),
                "timeout_seconds": receipt["timeout_seconds"],
                "solve_seconds": receipt["solve_seconds"],
                "result": receipt["result"],
            }
            for receipt, path in ((sat121, SAT_121), (sat120, SAT_120))
        ],
        "width64_c38_naive_scaling": {
            "left_list_size_radius_5": width64_left,
            "right_list_size_radius_6": width64_right,
            "cartesian_pairs_per_chart": width64_pairs,
            "chart_count": chart_count,
            "cartesian_pairs_all_charts": width64_pairs * chart_count,
            "pair_count_ratio_to_width16": width64_pairs / raw["cartesian_pairs"],
            "measured_width16_pairs_per_second": pairs_per_second,
            "optimistic_seconds_per_width64_chart_at_width16_rate": optimistic_seconds_per_chart,
            "optimistic_years_all_charts_at_width16_rate": (
                optimistic_seconds_per_chart * chart_count / (365.25 * 24 * 3600)
            ),
            "qualification": (
                "This is an optimistic throughput extrapolation, not a width-64 benchmark. "
                "It ignores larger words, memory traffic, and list-construction costs."
            ),
        },
        "result": "JOINT_REDUCTION_VALIDATED_BUT_NAIVE_DECODERS_DO_NOT_SCALE",
        "scope_limitation": (
            "The exact enumeration proves the minimum only within one width-16 chart. "
            "Agreement with the independent full census validates the implementation. "
            "The SAT timeouts prove neither satisfiability nor unsatisfiability, and no "
            "width-64 chart has been decided."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "status": "PASS",
        "width16_chart_minimum": raw["minimum_weight"],
        "width64_pairs_per_chart": width64_pairs,
        "output": str(OUTPUT.relative_to(ROOT)),
    }, indent=2))


if __name__ == "__main__":
    main()
