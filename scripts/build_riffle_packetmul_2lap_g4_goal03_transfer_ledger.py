#!/usr/bin/env python3
"""Build the exact Goal 03 transfer and partial-probability ledger."""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
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
SUFFIX = EXPLORATIONS / "riffle_dp_2lap_g4_one_data_suffix_bound.json"
SUFFIX_AUDIT = (
    EXPLORATIONS / "riffle_dp_2lap_g4_one_data_suffix_bound_verification.json"
)
AUTONOMOUS = (
    EXPLORATIONS / "riffle_dp_2lap_g4_three_node_block_certificate_w8.json"
)
AUTONOMOUS_SOURCE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_three_node_block.cpp"
AUTONOMOUS_VERIFIER = ROOT / "scripts" / "verify_riffle_dp_2lap_g4_three_node_block.cpp"
TURN_OFF = (
    EXPLORATIONS / "riffle_dp_2lap_g4_authenticated_six_packet_turnoff.json"
)
TURN_OFF_AUDIT = (
    EXPLORATIONS
    / "riffle_dp_2lap_g4_authenticated_six_packet_turnoff_verification.json"
)
OUTPUT = CANDIDATE / "receipts" / "goal03_transfer_ledger.json"

TOTAL_NODES = 32_772
PACKETS_PER_NODE = 16
TOTAL_PACKET_POSITIONS = TOTAL_NODES * PACKETS_PER_NODE
DISTANCE_THRESHOLD = 188_766
AUTONOMOUS_BLOCK_NODES = 3
AUTONOMOUS_BLOCK_WEIGHT = 25
MINIMUM_ZERO_PREFIX_NODES = 22_653
SUFFIX_NODES = TOTAL_NODES - MINIMUM_ZERO_PREFIX_NODES
SUFFIX_PACKET_POSITIONS = SUFFIX_NODES * PACKETS_PER_NODE
PACKET_MULTIPLIER_CHOICES = 15


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def fraction(record: dict) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def fraction_record(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def log2_interval(value: Fraction, places: int = 12) -> list[str]:
    if value <= 0:
        raise ValueError("log2 interval requires a positive rational")
    with localcontext() as context:
        context.prec = 150
        estimate = (
            Decimal(value.numerator).ln() - Decimal(value.denominator).ln()
        ) / Decimal(2).ln()
        unit = Decimal(1).scaleb(-places)
        return [
            format(estimate.quantize(unit, rounding=ROUND_FLOOR), "f"),
            format(estimate.quantize(unit, rounding=ROUND_CEILING), "f"),
        ]


def main() -> None:
    manifest = load(MANIFEST)
    specification = SPECIFICATION.read_text()
    goal01 = load(GOAL01)
    goal02 = load(GOAL02)
    outer = load(OUTER)
    suffix = load(SUFFIX)
    suffix_audit = load(SUFFIX_AUDIT)
    autonomous = load(AUTONOMOUS)
    turn_off = load(TURN_OFF)
    turn_off_audit = load(TURN_OFF_AUDIT)

    if manifest["candidate_id"] != "riffle_packetmul_2lap_g4":
        raise RuntimeError("Goal 03: wrong active candidate")
    multiplier = manifest["packet_multiplier"]
    if not multiplier["preserves_zero_packet_support"]:
        raise RuntimeError("Goal 03: PacketMul no longer preserves packet support")
    if multiplier["scalar_domain"] != "GF(16)^*":
        raise RuntimeError("Goal 03: PacketMul scalar law changed")
    if multiplier["scope"] != "one scalar per input packet coordinate":
        raise RuntimeError("Goal 03: multiplier independence scope changed")
    if multiplier["application"] != "once_before_packet_permutation":
        raise RuntimeError("Goal 03: PacketMul location changed")
    for required_statement in [
        "All other mechanisms remain unchanged",
        "Let \\(\\Pi\\) be the unchanged packet permutation",
        "be the unchanged retained-state two-lap map",
    ]:
        if required_statement not in specification:
            raise RuntimeError(
                f"Goal 03: construction boundary statement missing: {required_statement}"
            )

    if goal02["status"] != "GOAL_02_PROVED":
        raise RuntimeError("Goal 03: Goal 02 is not proved")
    if not goal02["verified"]["q_cap_5_over_8"]:
        raise RuntimeError("Goal 03: Goal 02 did not close the Goal 01 lemma")
    if goal02["receipt_sha256"]["receipts\\goal01_terminal_zero_gate.json"] != digest(GOAL01):
        raise RuntimeError("Goal 03: Goal 02 authenticates a different Goal 01 gate")

    support_counts = {int(key): value for key, value in outer["outer_word_counts"].items()}
    expected_counts = {33: 26, 35: 36, 36: 3233, 37: 510, 38: 933, 39: 81090}
    if support_counts != expected_counts or outer["record_count"] != sum(expected_counts.values()):
        raise RuntimeError("Goal 03: one-data population changed")

    if suffix["population_receipt_sha256"] != digest(OUTER):
        raise RuntimeError("Goal 03: suffix receipt authenticates a different population")
    if suffix["three_node_certificate_sha256"] != digest(AUTONOMOUS):
        raise RuntimeError("Goal 03: suffix receipt authenticates a different autonomous certificate")
    if suffix_audit["bound_receipt_sha256"] != digest(SUFFIX):
        raise RuntimeError("Goal 03: independent suffix audit hash mismatch")
    if turn_off_audit["receipt_sha256"] != digest(TURN_OFF):
        raise RuntimeError("Goal 03: independent turnoff audit hash mismatch")

    if autonomous["autonomous_output_map"] != "U(o)=Acc(P(o))":
        raise RuntimeError("Goal 03: autonomous map changed")
    if autonomous["three_node_weight_lower_bound"] != AUTONOMOUS_BLOCK_WEIGHT:
        raise RuntimeError("Goal 03: autonomous block certificate changed")
    if autonomous["rejected_triples_found"] != 0:
        raise RuntimeError("Goal 03: autonomous certificate contains a rejection")
    complete_blocks = MINIMUM_ZERO_PREFIX_NODES // AUTONOMOUS_BLOCK_NODES
    certified_prefix_weight = complete_blocks * AUTONOMOUS_BLOCK_WEIGHT
    if certified_prefix_weight <= DISTANCE_THRESHOLD:
        raise RuntimeError("Goal 03: transferred prefix no longer crosses distance")
    if MINIMUM_ZERO_PREFIX_NODES - 1 >= 0:
        preceding_weight = (
            (MINIMUM_ZERO_PREFIX_NODES - 1) // AUTONOMOUS_BLOCK_NODES
        ) * AUTONOMOUS_BLOCK_WEIGHT
        if preceding_weight > DISTANCE_THRESHOLD:
            raise RuntimeError("Goal 03: zero-prefix cutoff is not minimal")

    if suffix["total_nodes"] != TOTAL_NODES:
        raise RuntimeError("Goal 03: suffix node count changed")
    if suffix["total_packet_positions"] != TOTAL_PACKET_POSITIONS:
        raise RuntimeError("Goal 03: packet-position count changed")
    if suffix["minimum_zero_prefix_nodes"] != MINIMUM_ZERO_PREFIX_NODES:
        raise RuntimeError("Goal 03: suffix cutoff changed")
    if suffix["suffix_packet_positions"] != SUFFIX_PACKET_POSITIONS:
        raise RuntimeError("Goal 03: suffix packet-position count changed")

    suffix_rows = []
    suffix_aggregate = Fraction()
    for row in suffix["rows"]:
        support = row["packet_support"]
        count = row["outer_word_count"]
        if expected_counts[support] != count:
            raise RuntimeError("Goal 03: suffix population row changed")
        placement = Fraction(
            comb(SUFFIX_PACKET_POSITIONS, support),
            comb(TOTAL_PACKET_POSITIONS, support),
        )
        recorded_placement = Fraction(
            int(row["suffix_placement_probability_numerator"]),
            int(row["suffix_placement_probability_denominator"]),
        )
        contribution = count * placement
        recorded_contribution = Fraction(
            int(row["contribution_numerator"]),
            int(row["contribution_denominator"]),
        )
        if placement != recorded_placement or contribution != recorded_contribution:
            raise RuntimeError(f"Goal 03: suffix arithmetic changed at support {support}")
        suffix_aggregate += contribution
        suffix_rows.append(
            {
                "packet_support": support,
                "outer_word_count": count,
                "placement_event": (
                    "every active packet cell is in the final 10119 nodes"
                ),
                "upper_bound": fraction_record(contribution),
                "upper_bound_log2_interval": log2_interval(contribution),
            }
        )

    recorded_suffix_aggregate = fraction(suffix["aggregate_contribution"])
    if suffix_aggregate != recorded_suffix_aggregate:
        raise RuntimeError("Goal 03: suffix aggregate changed")
    if suffix_audit["aggregate_numerator"] != str(suffix_aggregate.numerator):
        raise RuntimeError("Goal 03: independent suffix numerator mismatch")
    if suffix_audit["aggregate_denominator"] != str(suffix_aggregate.denominator):
        raise RuntimeError("Goal 03: independent suffix denominator mismatch")

    terminal_zero = fraction(
        goal01["simple_sufficient_lemma_instantiation"]["aggregate_bound"]
    )
    target = Fraction(1, 1 << 40)
    higher_support_suffix = sum(
        (
            Fraction(
                int(row["upper_bound"]["numerator"]),
                int(row["upper_bound"]["denominator"]),
            )
            for row in suffix_rows
            if row["packet_support"] != 33
        ),
        Fraction(),
    )
    closed_partial = terminal_zero + higher_support_suffix
    remaining_budget = target - closed_partial
    if remaining_budget <= 0:
        raise RuntimeError("Goal 03: authenticated closed rows already exceed 2^-40")

    old_lower = Fraction(
        int(turn_off["lower_family"]["probability_numerator"]),
        int(turn_off["lower_family"]["probability_denominator"]),
    )
    packetmul_lower = old_lower / PACKET_MULTIPLIER_CHOICES**33
    if turn_off_audit["verified_probability_numerator"] != str(old_lower.numerator):
        raise RuntimeError("Goal 03: turnoff verification numerator mismatch")
    if turn_off_audit["verified_probability_denominator"] != str(old_lower.denominator):
        raise RuntimeError("Goal 03: turnoff verification denominator mismatch")
    if turn_off_audit["verified_packet_support"] != 33:
        raise RuntimeError("Goal 03: turnoff support changed")
    if turn_off_audit["verified_two_lap_weight"] >= DISTANCE_THRESHOLD:
        raise RuntimeError("Goal 03: authenticated turnoff is not bad")

    classifications = [
        {
            "id": "outer-stage-and-packet-support",
            "classification": "transferred",
            "artifacts": [str(OUTER.relative_to(ROOT))],
            "reason": (
                "The outer stage is unchanged, and multiplication by a nonzero field "
                "element preserves exactly which packets are zero."
            ),
        },
        {
            "id": "uniform-packet-placement-law",
            "classification": "transferred",
            "artifacts": [str(SUFFIX.relative_to(ROOT))],
            "reason": (
                "Pi is unchanged and independent of A, so a support-h word still "
                "occupies a uniform h-subset of the 524352 labeled packet cells."
            ),
        },
        {
            "id": "two-lap-map-and-autonomous-prefix",
            "classification": "transferred",
            "artifacts": [str(AUTONOMOUS.relative_to(ROOT))],
            "reason": (
                "PacketMul changes the drive word before G but does not change G, U, "
                "or the universal statement for every nonzero entering state."
            ),
        },
        {
            "id": "one-data-final-10119-placement-row",
            "classification": "transferred",
            "artifacts": [
                str(SUFFIX.relative_to(ROOT)),
                str(SUFFIX_AUDIT.relative_to(ROOT)),
            ],
            "reason": (
                "The proof uses only packet support, uniform placement, and the "
                "universal autonomous certificate; it does not inspect packet values."
            ),
        },
        {
            "id": "support33-terminal-zero-row",
            "classification": "transferred-with-substitution",
            "artifacts": [
                str(GOAL01.relative_to(ROOT)),
                str(GOAL02.relative_to(ROOT)),
            ],
            "reason": (
                "The old fixed-value character program is replaced by the PacketMul "
                "Fourier reduction and the proved global q_chi <= 5/8 orbit bound."
            ),
        },
        {
            "id": "paused-fixed-value-character-component-program",
            "classification": "superseded",
            "artifacts": [
                "explorations/riffle_dp_2lap_g4_local_character_components.json",
                "explorations/riffle_dp_2lap_g4_component*_certificate.md",
                "explorations/riffle_dp_2lap_g4_component*_audit.json",
            ],
            "reason": (
                "Those certificates address fixed packet-value orbit sums. PacketMul "
                "changes the character average, and Goals 01-02 now prove the needed "
                "support-33 terminal-zero bound directly."
            ),
        },
        {
            "id": "fixed-value-and-binary-weight-diagnostics",
            "classification": "invalidated",
            "artifacts": [
                "explorations/riffle_dp_2lap_g4_one_data_suffix_probe.json",
                "explorations/riffle_dp_2lap_g4_support33_suffix_probe.json",
                "explorations/riffle_dp_2lap_g4_full_local_distance_probe.json",
            ],
            "reason": (
                "A nonzero GF(16) multiplier need not preserve a packet's binary "
                "Hamming weight or its fixed nibble value."
            ),
        },
        {
            "id": "bounded-fixed-multiset-double-turnoff-exclusion",
            "classification": "invalidated",
            "artifacts": [
                "explorations/riffle_dp_2lap_g4_authenticated_double_turnoff_search.json",
                "explorations/riffle_dp_2lap_g4_authenticated_double_turnoff_verification.json",
            ],
            "reason": (
                "The exclusion fixes the paused candidate's packet-value multiset. "
                "Independent multipliers can realize different value multisets."
            ),
        },
        {
            "id": "authenticated-six-packet-turnoff-lower-family",
            "classification": "transferred-with-substitution",
            "artifacts": [
                str(TURN_OFF.relative_to(ROOT)),
                str(TURN_OFF_AUDIT.relative_to(ROOT)),
            ],
            "reason": (
                "For each old favorable labeled placement, exactly one multiplier per "
                "active packet realizes the same 33 target values. This preserves a "
                "valid lower subfamily with the additional factor 15^-33."
            ),
        },
        {
            "id": "support33-nonzero-terminal-nonsuffix-row",
            "classification": "open",
            "artifacts": [],
            "reason": (
                "Neither the terminal-zero proof nor the zero-prefix certificate "
                "controls bad outputs when L is nonzero and the first occupied node is "
                "before node 22653."
            ),
        },
        {
            "id": "remaining-outer-population",
            "classification": "open",
            "artifacts": [],
            "reason": (
                "Goal 03 only transfers the authenticated one-data population through "
                "packet support 39; other outer words and the open parts of supports "
                "35-39 still require ledger rows."
            ),
        },
    ]

    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal03-transfer-ledger-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": "riffle_packetmul_2lap_g4",
        "evidence_label": "EXACT_TRANSFER_AUDIT_AND_PARTIAL_BOUND",
        "source_sha256": digest(Path(__file__).resolve()),
        "authenticated_inputs_sha256": {
            str(path.relative_to(ROOT)): digest(path)
            for path in [
                MANIFEST,
                SPECIFICATION,
                GOAL01,
                GOAL02,
                OUTER,
                AUTONOMOUS,
                AUTONOMOUS_SOURCE,
                AUTONOMOUS_VERIFIER,
                SUFFIX,
                SUFFIX_AUDIT,
                TURN_OFF,
                TURN_OFF_AUDIT,
            ]
        },
        "construction_boundary": {
            "unchanged": [
                "outer stage",
                "packet permutation Pi",
                "two-lap map G",
                "first-lap terminal state L",
                "retained-state injection J",
            ],
            "added": (
                "Independent uniform A_i in GF(16)^* for every input packet "
                "coordinate, applied once before Pi."
            ),
            "preserved_invariants": [
                "packet zero/nonzero support",
                "uniform support-h placement law under Pi",
                "the deterministic map G after its drive word is fixed",
            ],
            "not_preserved": [
                "nonzero packet value",
                "binary Hamming weight within a nonzero packet",
                "fixed packet-value multiset",
            ],
        },
        "transfer_classifications": classifications,
        "autonomous_prefix_transfer": {
            "three_node_weight_lower_bound": AUTONOMOUS_BLOCK_WEIGHT,
            "minimum_zero_prefix_nodes": MINIMUM_ZERO_PREFIX_NODES,
            "complete_three_node_blocks": complete_blocks,
            "certified_weight": certified_prefix_weight,
            "distance_threshold": DISTANCE_THRESHOLD,
            "strictly_above_distance": certified_prefix_weight > DISTANCE_THRESHOLD,
            "cutoff_is_minimal_for_this_block_argument": True,
        },
        "one_data_population": {
            "outer_word_counts": {str(key): value for key, value in expected_counts.items()},
            "total_words": sum(expected_counts.values()),
        },
        "disjoint_event_partition": {
            "universe": (
                "Bad outputs induced by the 85828 authenticated one-data outer words "
                "through packet support 39."
            ),
            "closed_support33_terminal_zero": (
                "For a support-33 word: output is bad and L=0, with no placement "
                "restriction."
            ),
            "empty_support33_suffix_nonzero_terminal": (
                "For a support-33 word: output is bad, L!=0, and all active packet "
                "cells are in the final 10119 nodes. The autonomous-prefix certificate "
                "makes this event empty."
            ),
            "closed_higher_support_suffix": (
                "For a support-35 through support-39 word: output is bad and all active "
                "packet cells are in the final 10119 nodes."
            ),
            "open_support33_complement": (
                "For a support-33 word: output is bad, L!=0, and at least one active "
                "packet cell is before the final 10119 nodes."
            ),
            "open_higher_support_remainder": (
                "For supports 35-39: output is bad and at least one active packet cell "
                "is before the final 10119 nodes."
            ),
            "pairwise_disjoint": True,
        },
        "closed_rows": {
            "support33_terminal_zero": {
                "outer_word_count": 26,
                "upper_bound": fraction_record(terminal_zero),
                "upper_bound_log2_interval": log2_interval(terminal_zero),
                "proof_inputs": ["Goal 01", "Goal 02 q_chi <= 5/8"],
            },
            "support33_suffix_nonzero_terminal": {
                "outer_word_count": 26,
                "exact_upper_bound": {"numerator": "0", "denominator": "1"},
                "reason": (
                    "A 22653-node zero-input prefix entered from L!=0 contributes at "
                    "least 188775, already above the required distance 188766."
                ),
            },
            "higher_support_suffix_by_support": [
                row for row in suffix_rows if row["packet_support"] != 33
            ],
            "higher_support_suffix_aggregate": {
                "upper_bound": fraction_record(higher_support_suffix),
                "upper_bound_log2_interval": log2_interval(higher_support_suffix),
            },
            "transferred_original_suffix_bound_all_supports": {
                "upper_bound": fraction_record(suffix_aggregate),
                "upper_bound_log2_interval": log2_interval(suffix_aggregate),
                "ledger_note": (
                    "The support-33 placement contribution is not added to the partial "
                    "ledger because its L=0 part is already charged globally and its "
                    "L!=0 part is empty."
                ),
            },
            "partial_upper_bound": fraction_record(closed_partial),
            "partial_upper_bound_log2_interval": log2_interval(closed_partial),
            "target": fraction_record(target),
            "remaining_numerical_budget_if_all_open_rows_are_later_bounded": fraction_record(
                remaining_budget
            ),
            "remaining_numerical_budget_log2_interval": log2_interval(remaining_budget),
            "partial_bound_below_target": closed_partial < target,
            "warning": (
                "The partial upper bound is not a construction bound. Open events have "
                "not been assigned zero probability."
            ),
        },
        "packetmul_turnoff_lower_family": {
            "outer_family_id": turn_off["outer_family_id"],
            "packet_support": 33,
            "two_lap_weight": turn_off["two_lap_weight"],
            "distance_threshold": DISTANCE_THRESHOLD,
            "old_fixed_value_probability": fraction_record(old_lower),
            "active_multiplier_schedule_probability": fraction_record(
                Fraction(1, PACKET_MULTIPLIER_CHOICES**33)
            ),
            "packetmul_probability_lower_bound": fraction_record(packetmul_lower),
            "packetmul_probability_lower_bound_log2_interval": log2_interval(packetmul_lower),
            "interpretation": (
                "This proves bad realized (Pi,A) schedules exist and gives a rigorous "
                "lower subfamily. It does not threaten the 2^-40 setup-failure target."
            ),
        },
        "first_construction_level_blocker": {
            "name": "support-33 nonzero-terminal nonsuffix row",
            "event": (
                "Union over the 26 authenticated support-33 words of: output weight "
                "below 188766, L!=0, and first occupied node r<22653."
            ),
            "why_first": (
                "It is the lowest packet-support row not closed by either Goals 01-02 "
                "or the transferred autonomous-prefix argument."
            ),
            "finite_next_lemma": (
                "Let r be the first occupied node and W_tail the two-lap output weight "
                "outside the first floor(r/3) complete autonomous blocks. Prove, after "
                "summing over the 26 support-33 words and r=0,...,22652, that the event "
                "L!=0 and W_tail <= 188765-25*floor(r/3) has probability at most the "
                "recorded remaining numerical budget."
            ),
            "number_of_first-node_strata": MINIMUM_ZERO_PREFIX_NODES,
            "tail_threshold": "188765 - 25*floor(r/3)",
            "required_aggregate_upper_bound": fraction_record(remaining_budget),
            "required_aggregate_upper_bound_log2_interval": log2_interval(remaining_budget),
        },
        "status": "GOAL_03_PARTIAL_LEDGER_COMPLETE_FULL_CONSTRUCTION_OPEN",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print(f"closed_partial_log2={log2_interval(closed_partial)}")
    print(f"remaining_budget_log2={log2_interval(remaining_budget)}")
    print(f"packetmul_turnoff_lower_log2={log2_interval(packetmul_lower)}")


if __name__ == "__main__":
    main()
