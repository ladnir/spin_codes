#!/usr/bin/env python3
"""Independent manifest replay for the Goal 02 obstruction outcome."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4" / "receipts"
PRIMARY = RECEIPTS / "goal02_coding_obstruction_audit.json"
PRIMARY_SOURCE = ROOT / "scripts" / "audit_riffle_packetmul_wrapmul_2lap_g4_goal02_obstruction.py"
OUTPUT = RECEIPTS / "goal02_coding_obstruction_independent.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    assert primary["source_sha256"] == sha256(PRIMARY_SOURCE)
    assert primary["result"] == "PASS_GOAL_02_OBSTRUCTION_OUTCOME; ROW_NOT_CLOSED"
    for key, digest in primary["receipt_sha256"].items():
        filename = {
            "tail_primary": "goal02_tail_list_primary.json",
            "tail_independent": "goal02_tail_list_independent.json",
            "low_primary": "goal02_low_components_primary.json",
            "low_independent": "goal02_low_components_independent.json",
            "sixth_primary": "goal02_sixth_moment_primary.json",
            "sixth_independent": "goal02_sixth_moment_independent.json",
            "dual_a4": "goal02_dual_a4_audit.json",
            "a6_energy": "goal02_a6_energy_interface.json",
            "histogram_countermodel": "goal02_pair_histogram_countermodel.json",
            "histogram_countermodel_independent": "goal02_pair_histogram_countermodel_independent.json",
            "spectral": "goal02_spectral_bias_interface.json",
            "spectral_independent": "goal02_spectral_bias_independent.json",
            "high_search_raw": "goal02_high_support_bias_search_raw.json",
            "high_search_audit": "goal02_high_support_bias_search_audit.json",
            "c8_solver": "goal02_c8_native_xor.json",
        }[key]
        assert sha256(RECEIPTS / filename) == digest

    part_i = primary["part_i"]
    part_ii = primary["part_ii_exact_progress"]
    obstruction = primary["precise_obstruction"]
    terminal = primary["terminal_outcome"]
    assert part_i["status"] == "PROVED_AND_INDEPENDENTLY_REPLAYED"
    assert (part_i["johnson_first_node"], part_i["late_list_cap"]) == (7060, 13_749)
    assert part_i["unique_first_node"] == 9176
    assert part_ii["exhausted_component_dimension_cap"] == 22
    assert part_ii["exhausted_low_response_states"] == 0
    assert part_ii["remaining_dual_weight6_cap"] == 139_270_335_354_994_389_965
    assert part_ii["sufficient_nonexceptional_absolute_bias_cap"] == 106_802
    ratio = obstruction["histogram_countermodel_to_cap_ratio"]
    assert int(ratio["numerator"]) > int(ratio["denominator"])
    assert terminal == {
        "kind": "PRECISE_AUTHENTICATED_CODING_THEORETIC_OBSTRUCTION",
        "support33_nonzero_terminal_row_closed": False,
        "probability_relevant_counterexample_found": False,
        "full_construction_claim": False,
        "performance_benchmark_run": False,
    }

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-coding-obstruction-independent-v1",
        "candidate": primary["candidate"],
        "evidence_label": "INDEPENDENT_GOAL_02_OBSTRUCTION_MANIFEST_REPLAY",
        "source_sha256": sha256(Path(__file__)),
        "primary_source_sha256": sha256(PRIMARY_SOURCE),
        "primary_receipt_sha256": sha256(PRIMARY),
        "dependency_count": len(primary["receipt_sha256"]),
        "terminal_outcome_replayed": terminal,
        "result": "PASS_INDEPENDENT_GOAL_02_OBSTRUCTION_OUTCOME",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")
    print(f"dependencies={payload['dependency_count']}")
    print("status=PASS_INDEPENDENT_GOAL_02_OBSTRUCTION_OUTCOME")


if __name__ == "__main__":
    main()
