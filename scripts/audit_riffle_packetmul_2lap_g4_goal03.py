#!/usr/bin/env python3
"""Independently audit the Goal 03 transfer ledger and its exact arithmetic."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
EXPLORATIONS = ROOT / "explorations"
MANIFEST = CANDIDATE / "manifest.json"
SPECIFICATION = CANDIDATE / "CONSTRUCTION.md"
GOAL01 = CANDIDATE / "receipts" / "goal01_terminal_zero_gate.json"
GOAL02 = CANDIDATE / "receipts" / "goal02_audit.json"
OUTER = EXPLORATIONS / "riffle_dp_g4_one_data_support39.json"
AUTONOMOUS = EXPLORATIONS / "riffle_dp_2lap_g4_three_node_block_certificate_w8.json"
AUTONOMOUS_SOURCE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_three_node_block.cpp"
AUTONOMOUS_VERIFIER = ROOT / "scripts" / "verify_riffle_dp_2lap_g4_three_node_block.cpp"
SUFFIX = EXPLORATIONS / "riffle_dp_2lap_g4_one_data_suffix_bound.json"
SUFFIX_AUDIT = EXPLORATIONS / "riffle_dp_2lap_g4_one_data_suffix_bound_verification.json"
TURN_OFF = EXPLORATIONS / "riffle_dp_2lap_g4_authenticated_six_packet_turnoff.json"
TURN_OFF_AUDIT = (
    EXPLORATIONS / "riffle_dp_2lap_g4_authenticated_six_packet_turnoff_verification.json"
)
GENERATOR = ROOT / "scripts" / "build_riffle_packetmul_2lap_g4_goal03_transfer_ledger.py"
LEDGER = CANDIDATE / "receipts" / "goal03_transfer_ledger.json"
PROOF = CANDIDATE / "proof" / "GOAL_03_TRANSFER_LEDGER_REPORT.md"
OUTPUT = CANDIDATE / "receipts" / "goal03_audit.json"

NODES = 32_772
PACKET_POSITIONS = 16 * NODES
SUFFIX_PACKET_POSITIONS = 16 * 10_119
DISTANCE = 188_766
PREFIX_NODES = 22_653
EXPECTED_COUNTS = {33: 26, 35: 36, 36: 3233, 37: 510, 38: 933, 39: 81090}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def as_fraction(record: dict) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def gf16_mul(left: int, right: int) -> int:
    result = 0
    for _ in range(4):
        if right & 1:
            result ^= left
        right >>= 1
        carry = left & 8
        left = (left << 1) & 15
        if carry:
            left ^= 3
    return result


def main() -> None:
    manifest = load(MANIFEST)
    specification = SPECIFICATION.read_text()
    goal01 = load(GOAL01)
    goal02 = load(GOAL02)
    outer = load(OUTER)
    autonomous = load(AUTONOMOUS)
    suffix = load(SUFFIX)
    suffix_audit = load(SUFFIX_AUDIT)
    turn_off = load(TURN_OFF)
    turn_off_audit = load(TURN_OFF_AUDIT)
    ledger = load(LEDGER)

    if ledger["source_sha256"] != digest(GENERATOR):
        raise RuntimeError("Goal 03 audit: generator hash mismatch")
    for relative, expected in ledger["authenticated_inputs_sha256"].items():
        if digest(ROOT / relative) != expected:
            raise RuntimeError(f"Goal 03 audit: authenticated input changed: {relative}")
    if digest(AUTONOMOUS_SOURCE) != "a6fc6a1c0febe89d6d9e88d26d09a7586135a63e96bafdf1821ae3b38067504c":
        raise RuntimeError("Goal 03 audit: autonomous primary source changed")
    if digest(AUTONOMOUS_VERIFIER) != "65d4ed4287b484587602a9b09ed0b53505daf521d336dcf7b1add3667bd09e5b":
        raise RuntimeError("Goal 03 audit: autonomous verifier source changed")

    multiplier = manifest["packet_multiplier"]
    if multiplier["modulus_hex"] != "0x13":
        raise RuntimeError("Goal 03 audit: field polynomial changed")
    if multiplier["sampling"] != "independent_uniform_rejection":
        raise RuntimeError("Goal 03 audit: multiplier law changed")
    if multiplier["scope"] != "one scalar per input packet coordinate":
        raise RuntimeError("Goal 03 audit: multiplier scope changed")
    for x in range(1, 16):
        images = {gf16_mul(a, x) for a in range(1, 16)}
        if images != set(range(1, 16)):
            raise RuntimeError("Goal 03 audit: nonzero multiplication is not bijective")
        if gf16_mul(x, 0) != 0:
            raise RuntimeError("Goal 03 audit: multiplication does not preserve zero")
    for statement in [
        "All other mechanisms remain unchanged",
        "Let \\(\\Pi\\) be the unchanged packet permutation",
        "be the unchanged retained-state two-lap map",
    ]:
        if statement not in specification:
            raise RuntimeError("Goal 03 audit: construction boundary changed")

    counts = {int(key): value for key, value in outer["outer_word_counts"].items()}
    if counts != EXPECTED_COUNTS or outer["record_count"] != sum(EXPECTED_COUNTS.values()):
        raise RuntimeError("Goal 03 audit: outer population mismatch")

    if autonomous["three_node_weight_lower_bound"] != 25:
        raise RuntimeError("Goal 03 audit: autonomous block weight mismatch")
    if autonomous["enumerated_nonzero_outputs"] != 5_130_659_560:
        raise RuntimeError("Goal 03 audit: autonomous enumeration count mismatch")
    if autonomous["rejected_triples_found"] != 0:
        raise RuntimeError("Goal 03 audit: autonomous certificate rejected a triple")
    blocks = PREFIX_NODES // 3
    prefix_weight = 25 * blocks
    if blocks != 7_551 or prefix_weight != 188_775 or prefix_weight <= DISTANCE:
        raise RuntimeError("Goal 03 audit: autonomous-prefix arithmetic failed")
    if 25 * ((PREFIX_NODES - 1) // 3) > DISTANCE:
        raise RuntimeError("Goal 03 audit: prefix cutoff is not minimal")
    if suffix["three_node_certificate_sha256"] != digest(AUTONOMOUS):
        raise RuntimeError("Goal 03 audit: suffix certificate link failed")
    if suffix_audit["bound_receipt_sha256"] != digest(SUFFIX):
        raise RuntimeError("Goal 03 audit: suffix independent-audit link failed")

    suffix_by_support: dict[int, Fraction] = {}
    suffix_total = Fraction()
    for row in suffix["rows"]:
        support = row["packet_support"]
        placement = Fraction(
            comb(SUFFIX_PACKET_POSITIONS, support),
            comb(PACKET_POSITIONS, support),
        )
        contribution = EXPECTED_COUNTS[support] * placement
        if placement != Fraction(
            int(row["suffix_placement_probability_numerator"]),
            int(row["suffix_placement_probability_denominator"]),
        ):
            raise RuntimeError(f"Goal 03 audit: placement mismatch for h={support}")
        if contribution != Fraction(
            int(row["contribution_numerator"]), int(row["contribution_denominator"])
        ):
            raise RuntimeError(f"Goal 03 audit: contribution mismatch for h={support}")
        suffix_by_support[support] = contribution
        suffix_total += contribution
    if suffix_total != as_fraction(suffix["aggregate_contribution"]):
        raise RuntimeError("Goal 03 audit: suffix aggregate mismatch")
    if str(suffix_total.numerator) != suffix_audit["aggregate_numerator"]:
        raise RuntimeError("Goal 03 audit: independent suffix numerator mismatch")
    if str(suffix_total.denominator) != suffix_audit["aggregate_denominator"]:
        raise RuntimeError("Goal 03 audit: independent suffix denominator mismatch")

    if goal02["four_state_minimum"] != 24 or goal02["orbit_blocks"] != NODES // 4:
        raise RuntimeError("Goal 03 audit: Goal 02 four-state inputs changed")
    orbit_weight = 6 * NODES
    zero_count = PACKET_POSITIONS - orbit_weight
    if orbit_weight != goal02["orbit_weight_lower_bound"]:
        raise RuntimeError("Goal 03 audit: Goal 02 orbit arithmetic mismatch")
    if zero_count != goal02["zero_nibble_count_upper_bound"]:
        raise RuntimeError("Goal 03 audit: Goal 02 zero-count arithmetic mismatch")
    if Fraction(zero_count, PACKET_POSITIONS) != Fraction(5, 8):
        raise RuntimeError("Goal 03 audit: Goal 02 q cap mismatch")

    distinctness = as_fraction(goal01["exact_distinctness_probability"])
    parseval = as_fraction(goal01["exact_packetmul_parseval_sum"])
    terminal_zero = (
        26
        * Fraction(1, 1 << 64)
        * (1 + (parseval - 1) * Fraction(3, 5) ** 31)
        / distinctness
    )
    recorded_terminal = as_fraction(
        goal01["simple_sufficient_lemma_instantiation"]["aggregate_bound"]
    )
    if terminal_zero != recorded_terminal:
        raise RuntimeError("Goal 03 audit: terminal-zero fraction mismatch")

    higher_suffix = sum(
        (value for support, value in suffix_by_support.items() if support != 33),
        Fraction(),
    )
    closed_partial = terminal_zero + higher_suffix
    target = Fraction(1, 1 << 40)
    remaining = target - closed_partial
    closed = ledger["closed_rows"]
    if terminal_zero != as_fraction(closed["support33_terminal_zero"]["upper_bound"]):
        raise RuntimeError("Goal 03 audit: terminal ledger row mismatch")
    if as_fraction(closed["support33_suffix_nonzero_terminal"]["exact_upper_bound"]):
        raise RuntimeError("Goal 03 audit: deterministic suffix row is not zero")
    if higher_suffix != as_fraction(closed["higher_support_suffix_aggregate"]["upper_bound"]):
        raise RuntimeError("Goal 03 audit: higher-support suffix aggregate mismatch")
    if suffix_total != as_fraction(
        closed["transferred_original_suffix_bound_all_supports"]["upper_bound"]
    ):
        raise RuntimeError("Goal 03 audit: transferred suffix aggregate mismatch")
    if closed_partial != as_fraction(closed["partial_upper_bound"]):
        raise RuntimeError("Goal 03 audit: partial bound mismatch")
    if remaining != as_fraction(
        closed["remaining_numerical_budget_if_all_open_rows_are_later_bounded"]
    ):
        raise RuntimeError("Goal 03 audit: remaining budget mismatch")
    if remaining <= 0:
        raise RuntimeError("Goal 03 audit: partial rows exhaust the target")

    # Independently verify that the stated event classes partition every bad
    # one-data outcome. The booleans are terminal_zero and in_suffix.
    for support in EXPECTED_COUNTS:
        for terminal_is_zero in [False, True]:
            for in_suffix in [False, True]:
                events = [
                    support == 33 and terminal_is_zero,
                    support == 33 and not terminal_is_zero and in_suffix,
                    support == 33 and not terminal_is_zero and not in_suffix,
                    support != 33 and in_suffix,
                    support != 33 and not in_suffix,
                ]
                if sum(events) != 1:
                    raise RuntimeError("Goal 03 audit: event partition is not disjoint")

    expected_classes = {
        "outer-stage-and-packet-support": "transferred",
        "uniform-packet-placement-law": "transferred",
        "two-lap-map-and-autonomous-prefix": "transferred",
        "one-data-final-10119-placement-row": "transferred",
        "support33-terminal-zero-row": "transferred-with-substitution",
        "paused-fixed-value-character-component-program": "superseded",
        "fixed-value-and-binary-weight-diagnostics": "invalidated",
        "bounded-fixed-multiset-double-turnoff-exclusion": "invalidated",
        "authenticated-six-packet-turnoff-lower-family": "transferred-with-substitution",
        "support33-nonzero-terminal-nonsuffix-row": "open",
        "remaining-outer-population": "open",
    }
    actual_classes = {
        row["id"]: row["classification"] for row in ledger["transfer_classifications"]
    }
    if actual_classes != expected_classes:
        raise RuntimeError("Goal 03 audit: transfer classification set changed")

    old_lower = Fraction(
        int(turn_off["lower_family"]["probability_numerator"]),
        int(turn_off["lower_family"]["probability_denominator"]),
    )
    if turn_off_audit["receipt_sha256"] != digest(TURN_OFF):
        raise RuntimeError("Goal 03 audit: turnoff independent-audit link failed")
    if turn_off_audit["verified_packet_support"] != 33:
        raise RuntimeError("Goal 03 audit: turnoff support mismatch")
    if turn_off_audit["verified_two_lap_weight"] >= DISTANCE:
        raise RuntimeError("Goal 03 audit: turnoff is not below distance")
    packetmul_lower = old_lower / 15**33
    if packetmul_lower != as_fraction(
        ledger["packetmul_turnoff_lower_family"]["packetmul_probability_lower_bound"]
    ):
        raise RuntimeError("Goal 03 audit: PacketMul turnoff factor mismatch")

    blocker = ledger["first_construction_level_blocker"]
    if blocker["number_of_first-node_strata"] != PREFIX_NODES:
        raise RuntimeError("Goal 03 audit: blocker stratum count mismatch")
    if blocker["tail_threshold"] != "188765 - 25*floor(r/3)":
        raise RuntimeError("Goal 03 audit: blocker threshold mismatch")
    if as_fraction(blocker["required_aggregate_upper_bound"]) != remaining:
        raise RuntimeError("Goal 03 audit: blocker budget mismatch")

    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal03-audit-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "evidence_label": "EXACT_INDEPENDENT_AUDIT",
        "audit_source_sha256": digest(Path(__file__).resolve()),
        "generator_source_sha256": digest(GENERATOR),
        "ledger_sha256": digest(LEDGER),
        "proof_sha256": digest(PROOF),
        "verified": {
            "construction_boundary": True,
            "gf16_nonzero_multiplier_bijections": True,
            "packet_support_invariance": True,
            "outer_population_counts": True,
            "autonomous_primary_and_verifier_source_hashes": True,
            "autonomous_prefix_arithmetic": True,
            "suffix_hypergeometric_rows": True,
            "suffix_independent_receipt": True,
            "goal02_q_cap_reconstruction": True,
            "goal01_terminal_zero_fraction": True,
            "disjoint_event_partition": True,
            "closed_partial_bound_and_remaining_budget": True,
            "turnoff_independent_receipt": True,
            "packetmul_turnoff_factor_15_to_minus_33": True,
            "transfer_classification_set": True,
            "first_blocker_threshold_and_budget": True,
        },
        "closed_partial_below_2_to_minus_40": closed_partial < target,
        "status": "GOAL_03_AUDIT_PASSED",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print("status=GOAL_03_AUDIT_PASSED")


if __name__ == "__main__":
    main()
