#!/usr/bin/env python3
"""Consolidate Goal 02 under its precise-obstruction completion criterion."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4" / "receipts"
OUTPUT = RECEIPTS / "goal02_coding_obstruction_audit.json"

FILES = {
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
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    paths = {name: RECEIPTS / filename for name, filename in FILES.items()}
    data = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in paths.items()}

    tail = data["tail_primary"]
    tail_replay = data["tail_independent"]
    assert tail["status"] == "GOAL_02_PART_I_PROVED_EARLY_SPECTRUM_OPEN"
    assert tail_replay["status"] == "GOAL_02_PART_I_INDEPENDENTLY_VERIFIED"
    assert tail_replay["primary_receipt_sha256"] == sha256(paths["tail_primary"])
    assert tail["johnson_region"]["first_applicable_node"] == 7060
    assert tail["johnson_region"]["maximum_list_size"] == 13_749
    assert tail["unique_region"]["first_unique_node"] == 9176
    assert tail["unique_region"]["maximum_list_size"] == 1
    assert tail["early_region_target"]["maximum_uniform_bad_states_per_fixed_drive"] == 316_606
    assert tail["early_region_target"]["sufficient_difference_spectrum_cap"] == 316_605
    assert tail["early_region_target"]["difference_weight_maximum"] == 377_530

    low = data["low_primary"]
    low_replay = data["low_independent"]
    assert low["result"] == low_replay["result"] == "PASS"
    assert low_replay["primary_receipt_sha256"] == sha256(paths["low_primary"])
    assert low["component_dimension_cap"] == 22
    assert low["support_count"] == 35
    assert low["total_exact_nonzero_states"] == 15_650_667
    assert low["total_low_response_states"] == 0
    assert low["global_minimum_response_weight"] == 917_616

    assert data["sixth_primary"]["result"] == "EXACT_REDUCTION; SPARSE_DUAL_COUNTS_REQUIRED"
    assert data["sixth_independent"]["primary_receipt_sha256"] == sha256(paths["sixth_primary"])
    assert data["dual_a4"]["result"] == "PASS; EXACT_A4_CLOSED; A6_REQUIRED"
    assert data["dual_a4"]["dual_weight4_word_count"] == 3_732_023
    assert data["a6_energy"]["result"] == "EXACT_INTERFACE; LOCATION-SENSITIVE_A6_BOUND_REQUIRED"
    a6_cap = data["a6_energy"]["sufficient_target"]["maximum_dual_weight6_words"]
    triangle_cap = data["a6_energy"]["sufficient_target"]["maximum_triangle_energy"]
    assert a6_cap == 139_270_335_354_994_389_965
    assert triangle_cap == 12_543_557_200_693_295_704_650

    countermodel = data["histogram_countermodel"]
    countermodel_replay = data["histogram_countermodel_independent"]
    assert countermodel["result"] == "PASS_EXACT_HISTOGRAM_ONLY_OBSTRUCTION"
    assert countermodel_replay["result"] == "PASS_INDEPENDENT_HISTOGRAM_OBSTRUCTION_REPLAY"
    assert countermodel_replay["model_sha256"] == sha256(paths["histogram_countermodel"])
    countermodel_triangle = countermodel["triangle_energy"]["certified_lower_bound"]
    assert countermodel_triangle > triangle_cap

    spectral = data["spectral"]
    spectral_replay = data["spectral_independent"]
    assert spectral["result"] == "PASS_EXACT_SPECTRAL_INTERFACE; HIGH_SUPPORT_BIAS_LEMMA_REQUIRED"
    assert spectral_replay["result"] == "PASS_INDEPENDENT_SPECTRAL_INTERFACE_REPLAY"
    assert spectral_replay["interface_sha256"] == sha256(paths["spectral"])
    threshold = spectral["exact_threshold_derivation"]["largest_sufficient_even_absolute_bias"]
    assert threshold == 106_802
    assert spectral["low_component_coverage"]["support_dimension_cap"] == 22
    assert spectral["low_component_coverage"]["nonexceptional_maximum_absolute_bias"] == 104_882

    search = data["high_search_raw"]
    search_replay = data["high_search_audit"]
    assert search["result"] == "NO_VIOLATION_FOUND_IN_BOUNDED_SEARCH"
    assert search["searched_support_count"] == 92
    assert search["maximum_found_absolute_bias"] == 6_482
    assert search_replay["raw_receipt_sha256"] == sha256(paths["high_search_raw"])
    assert search_replay["result"] == "PASS_ALL_RETAINED_WITNESSES; NO_DIAGNOSTIC_VIOLATION"
    assert search_replay["coverage"]["replayed_witnesses"] == 184
    assert search_replay["affine_refutation_gate"]["minimum_replayed_high_support_response_weight"] == 1_045_463
    assert search_replay["affine_refutation_gate"]["difference_weight_required_for_two_bad_states"] == 377_530

    c8 = data["c8_solver"]
    assert c8["result"] == "UNKNOWN"
    assert c8["decision_problem"]["output_weight_upper_bound"] == 92
    assert c8["global_implication_if_unsat"]["full_response_distance_lower_bound"] == 380_964

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-coding-obstruction-audit-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "GOAL_02_PRECISE_AUTHENTICATED_CODING_OBSTRUCTION",
        "source_sha256": sha256(Path(__file__)),
        "receipt_sha256": {name: sha256(path) for name, path in paths.items()},
        "part_i": {
            "status": "PROVED_AND_INDEPENDENTLY_REPLAYED",
            "johnson_first_node": 7060,
            "late_list_cap": 13_749,
            "unique_first_node": 9176,
            "late_charge_log2_interval": tail["late_band_ledger"]["aggregate_charge_log2_interval"],
            "early_uniform_bad_state_cap": 316_606,
            "sufficient_difference_spectrum_cap": 316_605,
        },
        "part_ii_exact_progress": {
            "exhausted_component_dimension_cap": 22,
            "exhausted_nonzero_states": 15_650_667,
            "exhausted_low_response_states": 0,
            "dual_weight4_word_count": 3_732_023,
            "remaining_dual_weight6_cap": a6_cap,
            "equivalent_triangle_energy_cap": triangle_cap,
            "sufficient_nonexceptional_absolute_bias_cap": threshold,
        },
        "precise_obstruction": {
            "coding_statement": "For the actual response-coordinate pair-sum function r, prove T(r)<=12543557200693295704650; equivalently prove A6_dual<=139270335354994389965. A sufficient pointwise form is |n-2*wt(J(u))|<=106802 outside the four charged degree-one/degree-two states.",
            "remaining_state_scope": "all 92 exact irreducible-component supports of dimension at least 23",
            "histogram_insufficiency": "An exact countermodel has the same full pair multiplicity histogram, r(0)=0, L1 mass, L2 energy, and maximum multiplicity five, but triangle energy 4835703278451919629058050.",
            "histogram_countermodel_to_cap_ratio": countermodel["triangle_energy"]["ratio"],
            "alternative_c8_predicate": "UNSAT for a nonzero eight-node response word of weight at most 92 would close the spectrum; the authenticated bounded solver result is UNKNOWN.",
        },
        "refutation_track": {
            "high_support_masks_searched": 92,
            "exact_full_response_evaluations": search["exact_bias_evaluations"],
            "independently_replayed_extrema": 184,
            "maximum_found_absolute_bias": 6_482,
            "minimum_replayed_high_support_weight": 1_045_463,
            "affine_list_outcome": "No explored difference has weight at most 377530, so no two-state affine seed, reachable over-cap center, or setup-mass candidate was found.",
            "scope": "bounded diagnostic search; not an exhaustive refutation",
        },
        "terminal_outcome": {
            "kind": "PRECISE_AUTHENTICATED_CODING_THEORETIC_OBSTRUCTION",
            "support33_nonzero_terminal_row_closed": False,
            "probability_relevant_counterexample_found": False,
            "full_construction_claim": False,
            "performance_benchmark_run": False,
        },
        "result": "PASS_GOAL_02_OBSTRUCTION_OUTCOME; ROW_NOT_CLOSED",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")
    print(f"late_list_cap={payload['part_i']['late_list_cap']}")
    print(f"bias_cap={threshold}")
    print(f"histogram_countermodel_ratio={countermodel_triangle / triangle_cap:.12f}")
    print("status=PASS_GOAL_02_OBSTRUCTION_OUTCOME")


if __name__ == "__main__":
    main()
