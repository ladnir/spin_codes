#!/usr/bin/env python3
"""Audit Goal 04 and build its exact boundary-prefix ledger."""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
PRIMARY_SOURCE = ROOT / "scripts" / "certify_riffle_packetmul_2lap_g4_four_node_block.cpp"
INDEPENDENT_SOURCE = ROOT / "scripts" / "verify_riffle_packetmul_2lap_g4_four_node_block.cpp"
PRIMARY = CANDIDATE / "receipts" / "goal04_four_node_primary.json"
INDEPENDENT = CANDIDATE / "receipts" / "goal04_four_node_independent.json"
GOAL01 = CANDIDATE / "receipts" / "goal01_terminal_zero_gate.json"
GOAL03 = CANDIDATE / "receipts" / "goal03_transfer_ledger.json"
GOAL03_AUDIT = CANDIDATE / "receipts" / "goal03_audit.json"
OUTER = ROOT / "explorations" / "riffle_dp_g4_one_data_support39.json"
PROOF = CANDIDATE / "proof" / "GOAL_04_BOUNDARY_PREFIX_PROOF.md"
OUTPUT = CANDIDATE / "receipts" / "goal04_boundary_prefix_audit.json"

TOTAL_NODES = 32_772
PACKETS_PER_NODE = 16
TOTAL_PACKET_POSITIONS = TOTAL_NODES * PACKETS_PER_NODE
DISTANCE = 188_766
BAD_MAXIMUM = DISTANCE - 1
FOUR_NODE_BOUND = 36
FOUR_NODE_WIDTH = 4
OLD_CUTOFF = 22_653
SUPPORT_COUNTS = {33: 26, 35: 36, 36: 3233, 37: 510, 38: 933, 39: 81090}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def as_fraction(record: dict) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def fraction_record(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def log2_interval(value: Fraction, places: int = 12) -> list[str]:
    if value <= 0:
        raise ValueError("Goal 04 audit: logarithm requires a positive rational")
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
    primary = load(PRIMARY)
    independent = load(INDEPENDENT)
    goal01 = load(GOAL01)
    goal03 = load(GOAL03)
    goal03_audit = load(GOAL03_AUDIT)
    outer = load(OUTER)

    if goal03_audit["status"] != "GOAL_03_AUDIT_PASSED":
        raise RuntimeError("Goal 04 audit: Goal 03 audit is not valid")
    if goal03_audit["ledger_sha256"] != digest(GOAL03):
        raise RuntimeError("Goal 04 audit: Goal 03 ledger changed")

    if primary["enumerated_nonzero_outputs"] != 5_130_659_560:
        raise RuntimeError("Goal 04 audit: primary enumeration count mismatch")
    if primary["enumerated_output_weight_maximum"] != 8:
        raise RuntimeError("Goal 04 audit: primary anchor bound mismatch")
    if primary["rejected_four_node_weight_maximum"] != 35:
        raise RuntimeError("Goal 04 audit: primary rejection threshold mismatch")
    if primary["rejected_aligned_windows_found"] != 0:
        raise RuntimeError("Goal 04 audit: primary found a rejected window")
    if primary["four_node_weight_lower_bound"] != FOUR_NODE_BOUND:
        raise RuntimeError("Goal 04 audit: primary lower bound mismatch")
    if not primary["target_strata_closed"]:
        raise RuntimeError("Goal 04 audit: primary did not close target strata")

    if independent["enumerated_nonzero_outputs"] != primary["enumerated_nonzero_outputs"]:
        raise RuntimeError("Goal 04 audit: independent count mismatch")
    if independent["rejected_aligned_windows_found"] != 0:
        raise RuntimeError("Goal 04 audit: independent found a rejected window")
    if independent["minimum_covered_four_node_weight"] != primary[
        "minimum_weight_among_windows_containing_an_enumerated_output"
    ]:
        raise RuntimeError("Goal 04 audit: independent covered minimum mismatch")
    if independent["verified_four_node_weight_lower_bound"] != FOUR_NODE_BOUND:
        raise RuntimeError("Goal 04 audit: independent lower bound mismatch")
    if not independent["target_strata_closed"]:
        raise RuntimeError("Goal 04 audit: independent did not close target strata")

    # If four outputs have total weight at most 35, one has weight at most 8.
    # The covered minimum is 42, so the same enumeration proves weight 36.
    if 35 // FOUR_NODE_WIDTH != 8:
        raise RuntimeError("Goal 04 audit: low-anchor reduction failed")
    if primary["minimum_weight_among_windows_containing_an_enumerated_output"] <= 35:
        raise RuntimeError("Goal 04 audit: covered minimum does not prove weight 36")

    required_blocks = (DISTANCE + FOUR_NODE_BOUND - 1) // FOUR_NODE_BOUND
    new_cutoff = FOUR_NODE_WIDTH * required_blocks
    certified_weight = required_blocks * FOUR_NODE_BOUND
    prior_weight = (required_blocks - 1) * FOUR_NODE_BOUND
    if (required_blocks, new_cutoff, certified_weight, prior_weight) != (
        5_244,
        20_976,
        188_784,
        188_748,
    ):
        raise RuntimeError("Goal 04 audit: optimal block cutoff arithmetic failed")
    if not (prior_weight < DISTANCE <= certified_weight):
        raise RuntimeError("Goal 04 audit: cutoff is not minimal for the block proof")
    if not all(new_cutoff <= r < OLD_CUTOFF for r in [22_650, 22_651, 22_652]):
        raise RuntimeError("Goal 04 audit: target strata are outside the closed band")

    counts = {int(key): value for key, value in outer["outer_word_counts"].items()}
    if counts != SUPPORT_COUNTS or outer["record_count"] != sum(SUPPORT_COUNTS.values()):
        raise RuntimeError("Goal 04 audit: outer population mismatch")

    expanded_suffix_nodes = TOTAL_NODES - new_cutoff
    expanded_suffix_positions = PACKETS_PER_NODE * expanded_suffix_nodes
    suffix_rows = []
    higher_support_suffix = Fraction()
    for support, count in SUPPORT_COUNTS.items():
        placement = Fraction(
            comb(expanded_suffix_positions, support),
            comb(TOTAL_PACKET_POSITIONS, support),
        )
        contribution = count * placement
        suffix_rows.append(
            {
                "packet_support": support,
                "outer_word_count": count,
                "placement_probability": fraction_record(placement),
                "placement_probability_log2_interval": log2_interval(placement),
                "population_contribution": fraction_record(contribution),
                "population_contribution_log2_interval": log2_interval(contribution),
                "charged_in_retained_goal03_ledger": False,
            }
        )
        if support != 33:
            higher_support_suffix += contribution

    terminal_zero = as_fraction(
        goal01["simple_sufficient_lemma_instantiation"]["aggregate_bound"]
    )
    if terminal_zero != as_fraction(
        goal03["closed_rows"]["support33_terminal_zero"]["upper_bound"]
    ):
        raise RuntimeError("Goal 04 audit: terminal-zero row changed")
    target = Fraction(1, 1 << 40)
    naive_expanded_partial = terminal_zero + higher_support_suffix
    retained_closed_partial = as_fraction(
        goal03["closed_rows"]["partial_upper_bound"]
    )
    retained_remaining = as_fraction(
        goal03["closed_rows"][
            "remaining_numerical_budget_if_all_open_rows_are_later_bounded"
        ]
    )
    if target - retained_closed_partial != retained_remaining:
        raise RuntimeError("Goal 04 audit: retained Goal 03 budget mismatch")

    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal04-boundary-prefix-audit-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": "riffle_packetmul_2lap_g4",
        "evidence_label": "EXACT_CERTIFICATE_AND_AUDIT",
        "source_sha256": digest(Path(__file__).resolve()),
        "authenticated_artifacts_sha256": {
            str(path.relative_to(ROOT)): digest(path)
            for path in [
                PRIMARY_SOURCE,
                INDEPENDENT_SOURCE,
                PRIMARY,
                INDEPENDENT,
                GOAL01,
                GOAL03,
                GOAL03_AUDIT,
                OUTER,
                PROOF,
            ]
        },
        "conditional_probability_space": {
            "base_experiment": (
                "For a fixed support-33 outer word x, sample Pi uniformly from packet "
                "bijections and sample the 33 active A_i independently and uniformly "
                "from GF(16)^*."
            ),
            "first_node_event": (
                "E_r is the event that no active packet cell occurs before node r and "
                "at least one active packet cell occurs in node r."
            ),
            "first_node_probability": (
                "Pr[E_r] = (C(M-16r,33)-C(M-16(r+1),33))/C(M,33)."
            ),
            "conditioned_law": (
                "Conditioned only on E_r, active multipliers remain iid uniform. "
                "Conditioning further on L(Pi(D_A x))!=0 generally couples Pi and A; "
                "the resulting law is uniform over feasible (Pi,A) realizations, not "
                "a product law."
            ),
            "certificate_avoids_conditioning": (
                "The four-node bound holds pointwise for every nonzero L and every "
                "later drive word. No independence claim is used after conditioning."
            ),
        },
        "four_node_certificate": {
            "enumerated_nonzero_outputs": primary["enumerated_nonzero_outputs"],
            "enumerated_output_weight_maximum": 8,
            "checked_alignments_per_anchor": 4,
            "rejected_windows_of_weight_at_most_35": 0,
            "four_node_weight_lower_bound": FOUR_NODE_BOUND,
            "strengthening_argument": (
                "Any four-node window of weight at most 35 contains an output of "
                "weight at most 8. The exhaustive covered minimum is 42, so no such "
                "window exists."
            ),
            "minimum_weight_seen_in_covered_windows": primary[
                "minimum_weight_among_windows_containing_an_enumerated_output"
            ],
            "minimum_seen_is_not_claimed_as_global_minimum": True,
            "primary_and_independent_agree": True,
        },
        "cutoff_improvement": {
            "old_minimum_zero_prefix_nodes": OLD_CUTOFF,
            "new_minimum_zero_prefix_nodes": new_cutoff,
            "additional_closed_first_node_strata": OLD_CUTOFF - new_cutoff,
            "closed_first_node_interval": [new_cutoff, OLD_CUTOFF - 1],
            "target_first_node_strata": [22_650, 22_651, 22_652],
            "complete_four_node_blocks": required_blocks,
            "certified_prefix_weight": certified_weight,
            "distance_threshold": DISTANCE,
            "minimal_for_four_node_partition_argument": True,
        },
        "support33_ledger_update": {
            "closed_event": (
                "For each authenticated support-33 word: L!=0 and first occupied node "
                "r>=20976. A bad output in this event is impossible."
            ),
            "exact_bad_event_probability": {"numerator": "0", "denominator": "1"},
            "newly_closed_first_node_strata": OLD_CUTOFF - new_cutoff,
            "retained_goal03_closed_partial_upper_bound": fraction_record(
                retained_closed_partial
            ),
            "retained_goal03_closed_partial_upper_bound_log2_interval": log2_interval(
                retained_closed_partial
            ),
            "retained_remaining_numerical_budget": fraction_record(retained_remaining),
            "retained_remaining_numerical_budget_log2_interval": log2_interval(
                retained_remaining
            ),
            "reason_no_new_charge_is_needed": (
                "The support-33 L=0 event was charged globally in Goal 03. The newly "
                "closed L!=0 event has probability zero."
            ),
        },
        "coarse_higher_support_expansion_diagnostic": {
            "suffix_nodes": expanded_suffix_nodes,
            "suffix_packet_positions": expanded_suffix_positions,
            "rows": suffix_rows,
            "support33_nonzero_terminal_suffix_bad_event_probability": {
                "numerator": "0",
                "denominator": "1",
            },
            "higher_support_suffix_aggregate": fraction_record(higher_support_suffix),
            "higher_support_suffix_aggregate_log2_interval": log2_interval(
                higher_support_suffix
            ),
            "support33_terminal_zero_aggregate": fraction_record(terminal_zero),
            "naive_expanded_partial_upper_bound": fraction_record(naive_expanded_partial),
            "naive_expanded_partial_upper_bound_log2_interval": log2_interval(
                naive_expanded_partial
            ),
            "target": fraction_record(target),
            "naive_expanded_partial_below_target": naive_expanded_partial < target,
            "interpretation": (
                "Charging the entire larger placement region for supports 35-39 would "
                "exceed the 2^-40 target. Goal 04 retains the narrower Goal 03 charges "
                "and records the new nonzero-terminal events as empty."
            ),
        },
        "first_remaining_blocker": {
            "name": "support-33 nonzero-terminal first-node strata r<=20975",
            "first_node_range": [0, new_cutoff - 1],
            "number_of_strata": new_cutoff,
            "available_autonomous_prefix_bound": (
                "B(r)=max(25*floor(r/3),36*floor(r/4))."
            ),
            "finite_next_lemma": (
                "Let W_tail be the output weight outside the complete autonomous "
                "blocks used by B(r). Bound, over the 26 authenticated support-33 "
                "words and r=0,...,20975, the event L!=0 and "
                "W_tail<=188765-B(r) by the recorded remaining numerical budget."
            ),
            "required_aggregate_upper_bound": fraction_record(retained_remaining),
            "required_aggregate_upper_bound_log2_interval": log2_interval(
                retained_remaining
            ),
        },
        "verified": {
            "primary_enumeration": True,
            "independent_enumeration": True,
            "low_anchor_reduction": True,
            "target_strata_closed_pointwise": True,
            "optimal_four_node_cutoff_arithmetic": True,
            "conditional_independence_not_assumed": True,
            "expanded_hypergeometric_diagnostic": True,
            "disjoint_support33_charge": True,
        },
        "status": "GOAL_04_PROVED_BOUNDARY_STRATA_CLOSED_FULL_CONSTRUCTION_OPEN",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print(f"new_zero_prefix_cutoff={new_cutoff}")
    print(f"additional_closed_strata={OLD_CUTOFF - new_cutoff}")
    print(f"retained_closed_partial_log2={log2_interval(retained_closed_partial)}")
    print(f"retained_remaining_budget_log2={log2_interval(retained_remaining)}")
    print(f"naive_expanded_partial_log2={log2_interval(naive_expanded_partial)}")
    print("status=GOAL_04_AUDIT_PASSED")


if __name__ == "__main__":
    main()
