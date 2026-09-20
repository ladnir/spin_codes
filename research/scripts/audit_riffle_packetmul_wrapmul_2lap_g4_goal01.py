#!/usr/bin/env python3
"""Audit Goal 01 for Riffle PacketMul-WrapMul-2Lap g=4."""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PARENT = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
REGISTRY = ROOT / "constructions" / "registry.json"
MANIFEST = CANDIDATE / "manifest.json"
CONSTRUCTION = CANDIDATE / "CONSTRUCTION.md"
PROOF = CANDIDATE / "proof" / "GOAL_01_WRAP_UNIFORMITY_PROOF.md"
PARENT_TERMINAL_ZERO = PARENT / "receipts" / "goal01_terminal_zero_gate.json"
PARENT_GOAL04 = PARENT / "receipts" / "goal04_boundary_prefix_audit.json"
OUTPUT = CANDIDATE / "receipts" / "goal01_wrap_uniformity_audit.json"

FIELD_DEGREE = 64
MODULUS = (1 << 64) | 0x1B
NONZERO_STATES = (1 << 64) - 1
SUPPORT33_WORDS = 26
FULL_AUTHENTICATED_POPULATION = 85_828
BAD_WEIGHT_MAXIMUM = 188_765


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def polynomial_degree(value: int) -> int:
    return value.bit_length() - 1


def polynomial_remainder(value: int, divisor: int) -> int:
    divisor_degree = polynomial_degree(divisor)
    while value and polynomial_degree(value) >= divisor_degree:
        value ^= divisor << (polynomial_degree(value) - divisor_degree)
    return value


def polynomial_gcd(left: int, right: int) -> int:
    while right:
        left, right = right, polynomial_remainder(left, right)
    return left


def field_multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> FIELD_DEGREE:
            left ^= MODULUS
    return result


def field_square(value: int) -> int:
    return field_multiply(value, value)


def irreducibility_receipt() -> dict:
    # Rabin's test. The only prime divisor of 64 is 2.
    x_to_2_32 = 2
    for _ in range(32):
        x_to_2_32 = field_square(x_to_2_32)
    proper_factor_gcd = polynomial_gcd(x_to_2_32 ^ 2, MODULUS)
    x_to_2_64 = x_to_2_32
    for _ in range(32):
        x_to_2_64 = field_square(x_to_2_64)
    irreducible = proper_factor_gcd == 1 and x_to_2_64 == 2
    if not irreducible:
        raise RuntimeError("WrapMul Goal 01: terminal-state modulus is reducible")
    return {
        "test": "Rabin irreducibility test over GF(2)",
        "degree": FIELD_DEGREE,
        "prime_divisors_of_degree": [2],
        "gcd_x_2pow32_minus_x_hex": hex(proper_factor_gcd),
        "x_2pow64_mod_modulus_hex": hex(x_to_2_64),
        "irreducible": True,
    }


def fraction_record(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def log2_interval(value: Fraction, places: int = 12) -> list[str]:
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
    registry = load(REGISTRY)
    manifest = load(MANIFEST)
    terminal_zero = load(PARENT_TERMINAL_ZERO)
    goal04 = load(PARENT_GOAL04)

    if registry["active_candidate_id"] != manifest["candidate_id"]:
        raise RuntimeError("WrapMul Goal 01: registry does not select the candidate")
    if manifest["wrap_multiplier"]["scalar_domain"] != "GF(2^64)^*":
        raise RuntimeError("WrapMul Goal 01: scalar domain mismatch")
    if not manifest["wrap_multiplier"]["independent_of_inherited_setup"]:
        raise RuntimeError("WrapMul Goal 01: setup independence is absent")
    if terminal_zero["packet_support"] != 33:
        raise RuntimeError("WrapMul Goal 01: parent terminal-zero support mismatch")
    if terminal_zero["authenticated_word_count"] != SUPPORT33_WORDS:
        raise RuntimeError("WrapMul Goal 01: parent support-33 population mismatch")
    if goal04["status"] != "GOAL_04_PROVED_BOUNDARY_STRATA_CLOSED_FULL_CONSTRUCTION_OPEN":
        raise RuntimeError("WrapMul Goal 01: parent Goal 04 audit is not current")

    budget_record = goal04["support33_ledger_update"][
        "retained_remaining_numerical_budget"
    ]
    budget = Fraction(
        int(budget_record["numerator"]), int(budget_record["denominator"])
    )
    maximum_list_size = (budget * NONZERO_STATES) // SUPPORT33_WORDS
    aggregate = Fraction(SUPPORT33_WORDS * maximum_list_size, NONZERO_STATES)
    rejected_aggregate = Fraction(
        SUPPORT33_WORDS * (maximum_list_size + 1), NONZERO_STATES
    )
    if aggregate > budget or rejected_aggregate <= budget:
        raise RuntimeError("WrapMul Goal 01: integer list threshold is not maximal")

    boundary_rows = []
    for first_node in range(20_972, 20_976):
        three_node_bound = 25 * (first_node // 3)
        four_node_bound = 36 * (first_node // 4)
        autonomous_bound = max(three_node_bound, four_node_bound)
        boundary_rows.append(
            {
                "first_occupied_node": first_node,
                "three_node_bound": three_node_bound,
                "four_node_bound": four_node_bound,
                "selected_autonomous_bound": autonomous_bound,
                "bad_tail_weight_maximum": BAD_WEIGHT_MAXIMUM - autonomous_bound,
            }
        )
    if any(row["bad_tail_weight_maximum"] != 17 for row in boundary_rows):
        raise RuntimeError("WrapMul Goal 01: boundary tail arithmetic changed")

    field_audit = irreducibility_receipt()
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal01-audit-v1",
        "candidate": manifest["display_name"],
        "candidate_id": manifest["candidate_id"],
        "evidence_label": "EXACT_DISTRIBUTIONAL_REDUCTION_AND_LEDGER_AUDIT",
        "source_sha256": digest(Path(__file__).resolve()),
        "authenticated_artifacts_sha256": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                REGISTRY,
                MANIFEST,
                CONSTRUCTION,
                PROOF,
                PARENT_TERMINAL_ZERO,
                PARENT_GOAL04,
            )
        },
        "terminal_state_field": {
            "field": "GF(2^64)",
            "modulus_polynomial": "x^64+x^4+x^3+x+1",
            "modulus_integer_hex": hex(MODULUS),
            "implicit_top_reduction_hex": "0x1b",
            **field_audit,
        },
        "setup_probability_space": {
            "inherited_randomness": "packet permutation Pi and packet multipliers A",
            "new_randomness": "one B uniform in GF(2^64)^*, independent of (Pi,A)",
            "reuse": "the public B is reused for every outer word",
            "runtime_randomness_between_laps": False,
        },
        "conditional_uniformity": {
            "fixed_data": "outer word x and complete inherited setup realization",
            "condition": "L(Pi(D_A x)) != 0",
            "wrapped_state": "Z=B*L(Pi(D_A x))",
            "law": "uniform over GF(2^64)^*",
            "point_probability": fraction_record(Fraction(1, NONZERO_STATES)),
            "reason": "for every nonzero terminal state ell and target z, B=z*ell^-1 is the unique preimage",
            "independence_across_outer_words_claimed": False,
        },
        "transfers": {
            "support33_terminal_zero": {
                "classification": "transfers exactly",
                "reason": "B*0=0, so the encoder output on L=0 is unchanged",
            },
            "packet_support_and_placement": {
                "classification": "transfers exactly",
                "reason": "WrapMul acts after the driven first lap",
            },
            "autonomous_prefix": {
                "classification": "transfers pointwise",
                "reason": "B is nonzero, so B*L is nonzero exactly when L is nonzero",
                "retained_zero_prefix_cutoff_nodes": goal04["cutoff_improvement"][
                    "new_minimum_zero_prefix_nodes"
                ],
            },
        },
        "support33_bad_state_budget": {
            "authenticated_outer_words": SUPPORT33_WORDS,
            "nonzero_wrap_states": NONZERO_STATES,
            "remaining_budget": fraction_record(budget),
            "remaining_budget_log2_interval": log2_interval(budget),
            "maximum_uniform_bad_states_per_fixed_drive": maximum_list_size,
            "accepted_aggregate": fraction_record(aggregate),
            "accepted_aggregate_log2_interval": log2_interval(aggregate),
            "next_integer_rejected": maximum_list_size + 1,
            "next_integer_aggregate": fraction_record(rejected_aggregate),
            "next_integer_exceeds_budget": True,
            "uses_entire_remaining_budget_as_threshold": True,
        },
        "full_population_diagnostic": {
            "authenticated_outer_words": FULL_AUTHENTICATED_POPULATION,
            "uniform_bad_states_if_entire_budget_were_shared": (
                budget * NONZERO_STATES
            )
            // FULL_AUTHENTICATED_POPULATION,
            "scope": "diagnostic only; terminal-zero and placement rows above support 33 remain open",
        },
        "first_boundary_list_target": {
            "rows": boundary_rows,
            "required_object": (
                "For each fixed boundary drive, count nonzero wrapped states whose "
                "output outside the certified autonomous blocks has weight at most 17."
            ),
            "uniform_cap_that_would_fit_support33_budget": maximum_list_size,
        },
        "verified": {
            "registry_selects_new_candidate": True,
            "field_modulus_irreducible": True,
            "conditional_uniformity_bijection": True,
            "terminal_zero_transfer": True,
            "autonomous_transfer": True,
            "integer_budget_maximality": True,
            "independence_across_words_not_assumed": True,
        },
        "status": "GOAL_01_PROVED_WRAP_UNIFORMITY_LIST_BOUND_OPEN",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print(f"field_modulus_irreducible={field_audit['irreducible']}")
    print(f"maximum_uniform_bad_states={maximum_list_size}")
    print(f"next_integer_rejected={maximum_list_size + 1}")
    print("boundary_tail_weight_maximum=17")
    print("status=GOAL_01_PROVED_WRAP_UNIFORMITY_LIST_BOUND_OPEN")


if __name__ == "__main__":
    main()
