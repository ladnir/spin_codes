#!/usr/bin/env python3
"""Audit the hash and arithmetic chain for PacketMul Goal 01."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
MANIFEST = CANDIDATE / "manifest.json"
PRIMARY = CANDIDATE / "receipts" / "goal01_zero_symbol_primary.json"
INDEPENDENT = CANDIDATE / "receipts" / "goal01_zero_symbol_independent.json"
SEARCH = CANDIDATE / "receipts" / "goal01_mixed_character_search.json"
GATE = CANDIDATE / "receipts" / "goal01_terminal_zero_gate.json"
PROOF = CANDIDATE / "proof" / "GOAL_01_ZERO_SYMBOL_GATE_CHECKPOINT.md"
OUTPUT = CANDIDATE / "receipts" / "goal01_audit.json"
M = 524_352


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_fraction(row: dict) -> Fraction:
    return Fraction(int(row["numerator"]), int(row["denominator"]))


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    independent = json.loads(INDEPENDENT.read_text(encoding="utf-8"))
    search = json.loads(SEARCH.read_text(encoding="utf-8"))
    gate = json.loads(GATE.read_text(encoding="utf-8"))
    manifest_hash = digest(MANIFEST)
    for name, receipt in (
        ("primary", primary),
        ("independent", independent),
        ("search", search),
        ("gate", gate),
    ):
        if receipt["candidate_id"] != manifest["candidate_id"]:
            raise RuntimeError(f"{name} candidate id differs")
        if receipt["active_manifest_sha256"] != manifest_hash:
            raise RuntimeError(f"{name} manifest hash differs")

    if independent["primary_receipt_sha256"] != digest(PRIMARY):
        raise RuntimeError("independent receipt does not bind primary")
    if search["primary_receipt_sha256"] != digest(PRIMARY):
        raise RuntimeError("search receipt does not bind primary")
    if independent["mixed_search_receipt_sha256"] != digest(SEARCH):
        raise RuntimeError("independent receipt does not bind search")
    if gate["primary_receipt_sha256"] != digest(PRIMARY):
        raise RuntimeError("gate does not bind primary")
    if gate["independent_receipt_sha256"] != digest(INDEPENDENT):
        raise RuntimeError("gate does not bind independent replay")
    if gate["mixed_search_receipt_sha256"] != digest(SEARCH):
        raise RuntimeError("gate does not bind mixed search")

    collision = primary["cross_value_collisions"]
    histogram = {int(key): int(value) for key, value in collision["multiplicity_histogram"].items()}
    labeled = sum(multiplicity * frequency for multiplicity, frequency in histogram.items())
    distinct = sum(histogram.values())
    ordered = sum(
        multiplicity * multiplicity * frequency
        for multiplicity, frequency in histogram.items()
    )
    if labeled != collision["labeled_samples"] or distinct != collision["distinct_states"]:
        raise RuntimeError("collision histogram masses differ")
    if ordered != int(collision["ordered_collision_count"]):
        raise RuntimeError("collision ordered count differs")
    parseval = Fraction((1 << 64) * ordered, labeled * labeled)
    if parseval != as_fraction(collision["normalized_parseval_sum"]):
        raise RuntimeError("collision Parseval arithmetic differs")
    if parseval != as_fraction(gate["exact_packetmul_parseval_sum"]):
        raise RuntimeError("gate Parseval value differs")

    independent_replays = {
        row["character_hex"]: row for row in independent["maximizer_replays"]
    }
    for component in primary["pure_component_rows"]:
        for key in ("maximum_zero_symbol", "maximum_absolute_beta"):
            character = component[key]["character_hex"]
            if character not in independent_replays:
                raise RuntimeError(f"missing independent replay for {character}")
    best = search["best_exact_character_found"]
    mixed = independent["best_mixed_witness_replay"]
    if mixed["character_hex"] != best["character_hex"]:
        raise RuntimeError("independent mixed character differs")
    if mixed["zero_symbol_count"] != best["zero_symbol_count"]:
        raise RuntimeError("independent mixed zero count differs")

    cap_row = gate["required_uniform_absolute_beta_cap"]
    cap_low = as_fraction(cap_row["rigorous_lower_rational"])
    cap_high = as_fraction(cap_row["rigorous_upper_rational"])
    distinct_probability = as_fraction(gate["exact_distinctness_probability"])
    required_power = (
        distinct_probability * (1 << 24) / 26 - 1
    ) / (parseval - 1)
    if not cap_low**31 <= required_power < cap_high**31:
        raise RuntimeError("required cap bracket differs")
    simple = gate["simple_sufficient_lemma_instantiation"]
    simple_aggregate = as_fraction(simple["aggregate_bound"])
    if simple["q_cap"] != "5/8" or simple["implied_absolute_beta_cap"] != "3/5":
        raise RuntimeError("simple lemma changed")
    if not simple_aggregate < Fraction(1, 1 << 40):
        raise RuntimeError("simple lemma no longer closes target")

    proof_text = PROOF.read_text(encoding="utf-8")
    for required_text in (
        "Zero-symbol lemma",
        "q_\\chi\\le5/8",
        "goal01_zero_symbol_primary.json",
        "goal01_zero_symbol_independent.json",
        "goal01_mixed_character_search.json",
        "goal01_terminal_zero_gate.json",
    ):
        if required_text not in proof_text:
            raise RuntimeError(f"proof checkpoint omits {required_text!r}")

    sources = (
        ROOT / "scripts" / "analyze_riffle_packetmul_2lap_g4_zero_symbol.py",
        ROOT / "scripts" / "verify_riffle_packetmul_2lap_g4_zero_symbol.py",
        ROOT / "scripts" / "search_riffle_packetmul_2lap_g4_mixed_zero_symbol.py",
        ROOT / "scripts" / "derive_riffle_packetmul_2lap_g4_terminal_zero_gate.py",
        Path(__file__).resolve(),
    )
    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal01-audit-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": manifest["candidate_id"],
        "evidence_label": "EXACT_AUDIT",
        "active_manifest_sha256": manifest_hash,
        "source_sha256": {
            str(path.relative_to(ROOT)): digest(path) for path in sources
        },
        "receipt_sha256": {
            str(path.relative_to(CANDIDATE)): digest(path)
            for path in (PRIMARY, INDEPENDENT, SEARCH, GATE)
        },
        "proof_checkpoint_sha256": digest(PROOF),
        "verified": {
            "candidate_and_manifest_chain": True,
            "primary_independent_search_hash_chain": True,
            "collision_histogram_mass": True,
            "ordered_collision_count": True,
            "parseval_identity": True,
            "all_reported_pure_maximizers_independently_replayed": True,
            "best_mixed_witness_independently_replayed": True,
            "required_cap_rigorous_bracket": True,
            "q_cap_5_over_8_closes_aggregate_2_to_minus_40": True,
            "proof_checkpoint_names_receipts_and_open_lemma": True,
        },
        "packet_cells": M,
        "status": "GOAL_01_COMPLETE_OPEN_MIXED_LEMMA",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")
    print("status=GOAL_01_EXACT_AUDIT_PASSED")


if __name__ == "__main__":
    main()
