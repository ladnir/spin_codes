#!/usr/bin/env python3
"""Authenticate the Goal 02B spectral-bias interface.

This audit does not prove the remaining global bias lemma.  It proves that
the stated pointwise lemma, if established outside four charged states,
implies the exact triangle-energy target used by Goal 02.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RECEIPTS = CONSTRUCTION / "receipts"

A6_RECEIPT = RECEIPTS / "goal02_a6_energy_interface.json"
A4_RECEIPT = RECEIPTS / "goal02_dual_a4_audit.json"
LOW_PRIMARY = RECEIPTS / "goal02_low_components_primary.json"
LOW_INDEPENDENT = RECEIPTS / "goal02_low_components_independent.json"

A6_SOURCE = ROOT / "scripts" / "audit_riffle_packetmul_wrapmul_2lap_g4_goal02_a6_energy.py"
LOW_PRIMARY_SOURCE = ROOT / "scripts" / "certify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components.cpp"
LOW_PRIMARY_EXE = ROOT / "scripts" / "certify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components.exe"
LOW_PRIMARY_RAW = RECEIPTS / "goal02_low_components_primary_raw.json"
LOW_FINALIZER = ROOT / "scripts" / "finalize_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components.py"
LOW_INDEPENDENT_SOURCE = ROOT / "scripts" / "verify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components.py"

OUTPUT = RECEIPTS / "goal02_spectral_bias_interface.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def parity(x: int) -> int:
    return x.bit_count() & 1


def replay_small_model(m: int, columns: list[int]) -> dict:
    """Replay both Fourier identities by direct enumeration."""

    group_size = 1 << m
    n = len(columns)
    assert len(set(columns)) == n
    assert all(0 < a < group_size for a in columns)

    pair_sums = Counter(a ^ b for a, b in combinations(columns, 2))
    pair_mass = sum(pair_sums.values())
    pair_energy = sum(v * v for v in pair_sums.values())
    triangle = sum(
        pair_sums[x] * pair_sums[y] * pair_sums[x ^ y]
        for x in range(group_size)
        for y in range(group_size)
    )

    lambdas: list[int] = []
    for u in range(group_size):
        response_bias = sum(1 if parity(a & u) == 0 else -1 for a in columns)
        lambda_from_bias = (response_bias * response_bias - n) // 2
        lambda_from_pairs = sum(
            count * (1 if parity(x & u) == 0 else -1)
            for x, count in pair_sums.items()
        )
        assert lambda_from_pairs == lambda_from_bias
        lambdas.append(lambda_from_bias)

    assert sum(x * x for x in lambdas) == group_size * pair_energy
    assert sum(x * x * x for x in lambdas) == group_size * triangle
    return {
        "dimension": m,
        "columns": columns,
        "pair_mass": pair_mass,
        "pair_energy_l2_squared": pair_energy,
        "triangle_energy": triangle,
        "lambda_minimum": min(lambdas),
        "lambda_maximum": max(lambdas),
        "status": "PASS",
    }


def main() -> None:
    a6 = load(A6_RECEIPT)
    low = load(LOW_PRIMARY)
    independent = load(LOW_INDEPENDENT)

    assert a6["schema"] == "riffle-packetmul-wrapmul-2lap-g4-goal02-a6-energy-v1"
    assert a6["source_sha256"] == sha256(A6_SOURCE)
    assert a6["dual_a4_audit_sha256"] == sha256(A4_RECEIPT)
    assert a6["result"] == "EXACT_INTERFACE; LOCATION-SENSITIVE_A6_BOUND_REQUIRED"

    assert low["source_sha256"] == sha256(LOW_PRIMARY_SOURCE)
    assert low["executable_sha256"] == sha256(LOW_PRIMARY_EXE)
    assert low["raw_receipt_sha256"] == sha256(LOW_PRIMARY_RAW)
    assert low["finalizer_source_sha256"] == sha256(LOW_FINALIZER)
    assert low["result"] == "PASS"

    assert independent["source_sha256"] == sha256(LOW_INDEPENDENT_SOURCE)
    assert independent["primary_receipt_sha256"] == sha256(LOW_PRIMARY)
    assert independent["primary_source_sha256"] == sha256(LOW_PRIMARY_SOURCE)
    assert independent["primary_executable_sha256"] == sha256(LOW_PRIMARY_EXE)
    assert independent["result"] == "PASS"
    assert independent["rows_replayed"] == low["support_count"] == 35

    n = a6["response_coordinate_count"]
    group_size = 1 << 64
    pair_mass = a6["pair_mass_l1"]
    pair_energy = a6["pair_energy_l2_squared"]
    triangle_cap = a6["sufficient_target"]["maximum_triangle_energy"]
    assert n == 2_097_408
    assert pair_mass == n * (n - 1) // 2

    degrees = low["component_degrees"]
    expected_masks = {
        mask
        for mask in range(1, 1 << len(degrees))
        if sum(degrees[i] for i in range(len(degrees)) if mask & (1 << i)) <= 22
    }
    rows = {int(row["component_support_mask_hex"], 16): row for row in low["rows"]}
    assert set(rows) == expected_masks
    assert all(row["dimension"] <= 22 for row in rows.values())

    degree_one = rows[0x1]
    degree_two = rows[0x2]
    assert degree_one["exact_nonzero_state_count"] == 1
    assert degree_one["maximum_absolute_response_bias"] == 196_632
    assert degree_two["exact_nonzero_state_count"] == 3
    assert degree_two["maximum_absolute_response_bias"] == 262_176

    nonexceptional_rows = [row for mask, row in rows.items() if mask not in (0x1, 0x2)]
    certified_nonexceptional_maximum = max(
        row["maximum_absolute_response_bias"] for row in nonexceptional_rows
    )
    assert certified_nonexceptional_maximum == 104_882
    assert low["global_maximum_absolute_response_bias"] == 262_176
    assert independent["global_maximum_absolute_response_bias"] == 262_176

    def eigenvalue_from_bias(bias: int) -> int:
        numerator = bias * bias - n
        assert numerator % 2 == 0
        return numerator // 2

    exceptional_bias_caps = [196_632, 262_176, 262_176, 262_176]
    exceptional_lambda_caps = [eigenvalue_from_bias(b) for b in exceptional_bias_caps]

    # The conservative charge keeps the exceptional square mass inside the
    # Parseval remainder.  It therefore needs only upper bounds on their
    # individual biases, not equality for every exceptional state.
    parseval_remainder = group_size * pair_energy - pair_mass * pair_mass
    scaled_budget_after_trivial_and_exceptions = (
        triangle_cap * group_size
        - pair_mass**3
        - sum(value**3 for value in exceptional_lambda_caps)
    )
    maximum_remaining_lambda = (
        scaled_budget_after_trivial_and_exceptions // parseval_remainder
    )

    # Response bias is even because n is even.  Find the largest nonnegative
    # even bias whose induced eigenvalue is within the exact remaining cap.
    bias_cap = 0
    lo = 0
    hi = n // 2
    while lo <= hi:
        mid = (lo + hi) // 2
        bias = 2 * mid
        if eigenvalue_from_bias(bias) <= maximum_remaining_lambda:
            bias_cap = bias
            lo = mid + 1
        else:
            hi = mid - 1
    assert bias_cap == 106_802
    first_violating_even_bias = bias_cap + 2

    lambda_cap = eigenvalue_from_bias(bias_cap)
    next_lambda = eigenvalue_from_bias(first_violating_even_bias)
    assert lambda_cap <= maximum_remaining_lambda < next_lambda

    spectral_upper_numerator = (
        pair_mass**3
        + sum(value**3 for value in exceptional_lambda_caps)
        + lambda_cap * parseval_remainder
    )
    spectral_upper_ceiling = (
        spectral_upper_numerator + group_size - 1
    ) // group_size
    assert spectral_upper_numerator <= triangle_cap * group_size
    assert spectral_upper_ceiling <= triangle_cap

    next_numerator = (
        pair_mass**3
        + sum(value**3 for value in exceptional_lambda_caps)
        + next_lambda * parseval_remainder
    )
    assert next_numerator > triangle_cap * group_size

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-spectral-bias-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SPECTRAL_BIAS_INTERFACE",
        "source_sha256": sha256(Path(__file__)),
        "dependency_sha256": {
            "a6_energy_interface": sha256(A6_RECEIPT),
            "dual_a4_audit": sha256(A4_RECEIPT),
            "low_components_primary": sha256(LOW_PRIMARY),
            "low_components_independent": sha256(LOW_INDEPENDENT),
        },
        "response_code": {
            "dimension": 64,
            "length": n,
            "bias_definition": "B(u)=n-2*wt(J(u))",
            "pair_fourier_eigenvalue": "lambda_u=(B(u)^2-n)/2",
            "triangle_fourier_identity": "T(r)=2^-64*sum_u lambda_u^3",
            "pair_parseval_identity": "sum_u lambda_u^2=2^64*sum_x r(x)^2",
        },
        "authenticated_pair_data": {
            "pair_mass_l1": pair_mass,
            "pair_energy_l2_squared": pair_energy,
            "triangle_energy_cap": triangle_cap,
        },
        "charged_exceptions": {
            "description": "the nonzero pure degree-one and pure degree-two component states",
            "count": 4,
            "absolute_bias_caps": exceptional_bias_caps,
            "lambda_caps": exceptional_lambda_caps,
            "charging_rule": "charge their positive lambda cubes separately; retain their square mass in the Parseval remainder",
        },
        "low_component_coverage": {
            "support_dimension_cap": 22,
            "support_count": len(rows),
            "nonexceptional_maximum_absolute_bias": certified_nonexceptional_maximum,
            "remaining_states": "nonexceptional states with exact component-support dimension at least 23",
        },
        "exact_threshold_derivation": {
            "parseval_remainder_upper_bound": parseval_remainder,
            "scaled_budget_after_trivial_and_exceptions": scaled_budget_after_trivial_and_exceptions,
            "maximum_remaining_nonnegative_lambda": maximum_remaining_lambda,
            "largest_sufficient_even_absolute_bias": bias_cap,
            "lambda_at_largest_sufficient_bias": lambda_cap,
            "first_larger_even_absolute_bias": first_violating_even_bias,
            "lambda_at_first_larger_even_bias": next_lambda,
            "spectral_triangle_upper_bound_ceiling": spectral_upper_ceiling,
            "triangle_cap_margin": triangle_cap - spectral_upper_ceiling,
            "scaled_margin": triangle_cap * group_size - spectral_upper_numerator,
            "first_larger_even_bias_fails_this_conservative_charge_by_scaled_amount": next_numerator - triangle_cap * group_size,
        },
        "small_model_fourier_replays": [
            replay_small_model(3, [1, 2, 3, 4, 6]),
            replay_small_model(4, [1, 2, 4, 7, 8, 11, 13]),
            replay_small_model(5, [1, 3, 5, 9, 14, 17, 22, 27]),
        ],
        "sufficient_global_lemma": {
            "statement": "For every u outside zero and the four charged exceptions, |B(u)|<=106802.",
            "authenticated_scope_already_closed": "all nonexceptional exact component supports of dimension at most 22",
            "remaining_scope": "all exact component supports of dimension at least 23",
            "consequence": "T(r) is at most the authenticated triangle cap, hence the Goal 02 A6 and low-spectrum targets follow",
        },
        "scope_limit": "This receipt authenticates the implication and threshold. It does not prove the remaining high-support global lemma.",
        "result": "PASS_EXACT_SPECTRAL_INTERFACE; HIGH_SUPPORT_BIAS_LEMMA_REQUIRED",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")
    print(f"bias_cap={bias_cap}")
    print(f"certified_nonexceptional_maximum={certified_nonexceptional_maximum}")
    print(f"triangle_upper_ceiling={spectral_upper_ceiling}")
    print(f"triangle_cap_margin={triangle_cap - spectral_upper_ceiling}")
    print("status=PASS_EXACT_SPECTRAL_INTERFACE")


if __name__ == "__main__":
    main()
