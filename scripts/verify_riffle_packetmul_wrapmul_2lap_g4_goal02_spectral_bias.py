#!/usr/bin/env python3
"""Independent replay of the Goal 02B spectral-bias interface."""

from __future__ import annotations

import hashlib
import json
from itertools import combinations
from math import isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RECEIPTS = CONSTRUCTION / "receipts"

INTERFACE = RECEIPTS / "goal02_spectral_bias_interface.json"
A6 = RECEIPTS / "goal02_a6_energy_interface.json"
A4 = RECEIPTS / "goal02_dual_a4_audit.json"
LOW_PRIMARY = RECEIPTS / "goal02_low_components_primary.json"
LOW_INDEPENDENT = RECEIPTS / "goal02_low_components_independent.json"
INTERFACE_SOURCE = ROOT / "scripts" / "audit_riffle_packetmul_wrapmul_2lap_g4_goal02_spectral_bias.py"
OUTPUT = RECEIPTS / "goal02_spectral_bias_independent.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def fwht(values: list[int]) -> list[int]:
    result = values[:]
    width = 1
    while width < len(result):
        for base in range(0, len(result), 2 * width):
            for offset in range(width):
                x = result[base + offset]
                y = result[base + width + offset]
                result[base + offset] = x + y
                result[base + width + offset] = x - y
        width *= 2
    return result


def replay_small_model(row: dict) -> None:
    m = row["dimension"]
    columns = row["columns"]
    size = 1 << m
    counts = [0] * size
    for left, right in combinations(columns, 2):
        counts[left ^ right] += 1
    lambdas = fwht(counts)

    direct_lambdas = []
    for u in range(size):
        bias = 0
        for column in columns:
            bias += -1 if (column & u).bit_count() & 1 else 1
        direct_lambdas.append((bias * bias - len(columns)) // 2)
    assert lambdas == direct_lambdas

    energy = sum(value * value for value in counts)
    triangle = 0
    for x, rx in enumerate(counts):
        if rx == 0:
            continue
        for y, ry in enumerate(counts):
            if ry:
                triangle += rx * ry * counts[x ^ y]
    assert sum(value * value for value in lambdas) == size * energy
    assert sum(value**3 for value in lambdas) == size * triangle
    assert row["pair_energy_l2_squared"] == energy
    assert row["triangle_energy"] == triangle


def main() -> None:
    interface = read(INTERFACE)
    a6 = read(A6)
    low = read(LOW_PRIMARY)
    independent_low = read(LOW_INDEPENDENT)

    assert interface["source_sha256"] == sha256(INTERFACE_SOURCE)
    dependencies = interface["dependency_sha256"]
    assert dependencies == {
        "a6_energy_interface": sha256(A6),
        "dual_a4_audit": sha256(A4),
        "low_components_primary": sha256(LOW_PRIMARY),
        "low_components_independent": sha256(LOW_INDEPENDENT),
    }
    assert interface["result"] == (
        "PASS_EXACT_SPECTRAL_INTERFACE; HIGH_SUPPORT_BIAS_LEMMA_REQUIRED"
    )
    assert independent_low["primary_receipt_sha256"] == sha256(LOW_PRIMARY)
    assert independent_low["result"] == "PASS"

    n = a6["response_coordinate_count"]
    group_size = 1 << 64
    pair_mass = n * (n - 1) // 2
    assert pair_mass == a6["pair_mass_l1"]
    pair_energy = a6["pair_energy_l2_squared"]
    triangle_cap = a6["sufficient_target"]["maximum_triangle_energy"]

    rows = {int(row["component_support_mask_hex"], 16): row for row in low["rows"]}
    degrees = low["component_degrees"]
    expected = []
    for mask in range(1, 1 << len(degrees)):
        dimension = sum(
            degree for index, degree in enumerate(degrees) if mask & (1 << index)
        )
        if dimension <= 22:
            expected.append(mask)
    assert sorted(rows) == expected
    assert rows[1]["exact_nonzero_state_count"] == 1
    assert rows[2]["exact_nonzero_state_count"] == 3
    exception_biases = [
        rows[1]["maximum_absolute_response_bias"],
        *([rows[2]["maximum_absolute_response_bias"]] * 3),
    ]
    assert exception_biases == [196_632, 262_176, 262_176, 262_176]
    low_nonexceptional_maximum = max(
        row["maximum_absolute_response_bias"]
        for mask, row in rows.items()
        if mask not in (1, 2)
    )
    assert low_nonexceptional_maximum == 104_882

    eigenvalues = [(bias * bias - n) // 2 for bias in exception_biases]
    remainder = group_size * pair_energy - pair_mass * pair_mass
    available = (
        triangle_cap * group_size
        - pair_mass**3
        - sum(value**3 for value in eigenvalues)
    )
    lambda_limit, division_remainder = divmod(available, remainder)

    root = isqrt(2 * lambda_limit + n)
    bias_limit = root - (root & 1)
    while ((bias_limit + 2) ** 2 - n) // 2 <= lambda_limit:
        bias_limit += 2
    while (bias_limit * bias_limit - n) // 2 > lambda_limit:
        bias_limit -= 2
    assert bias_limit == 106_802

    lambda_limit_used = (bias_limit * bias_limit - n) // 2
    upper_numerator = (
        pair_mass**3
        + sum(value**3 for value in eigenvalues)
        + lambda_limit_used * remainder
    )
    upper_ceiling = (upper_numerator + group_size - 1) // group_size
    assert upper_ceiling <= triangle_cap

    next_bias = bias_limit + 2
    next_lambda = (next_bias * next_bias - n) // 2
    next_numerator = (
        pair_mass**3
        + sum(value**3 for value in eigenvalues)
        + next_lambda * remainder
    )
    assert next_numerator > triangle_cap * group_size

    derived = interface["exact_threshold_derivation"]
    assert derived["parseval_remainder_upper_bound"] == remainder
    assert derived["scaled_budget_after_trivial_and_exceptions"] == available
    assert derived["maximum_remaining_nonnegative_lambda"] == lambda_limit
    assert derived["largest_sufficient_even_absolute_bias"] == bias_limit
    assert derived["first_larger_even_absolute_bias"] == next_bias
    assert derived["spectral_triangle_upper_bound_ceiling"] == upper_ceiling

    for replay in interface["small_model_fourier_replays"]:
        replay_small_model(replay)

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-spectral-bias-independent-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "INDEPENDENT_EXACT_SPECTRAL_INTERFACE_REPLAY",
        "source_sha256": sha256(Path(__file__)),
        "interface_sha256": sha256(INTERFACE),
        "dependency_sha256": dependencies,
        "reconstruction": {
            "threshold_method": "integer square root after exact Euclidean division",
            "small_model_method": "independent Walsh-Hadamard transform of the pair-sum vectors",
            "parseval_remainder": remainder,
            "available_scaled_budget": available,
            "division_remainder": division_remainder,
            "maximum_remaining_lambda": lambda_limit,
            "largest_sufficient_even_absolute_bias": bias_limit,
            "first_larger_even_absolute_bias": next_bias,
            "triangle_upper_bound_ceiling": upper_ceiling,
            "triangle_cap": triangle_cap,
        },
        "coverage": {
            "low_component_supports_replayed_by_dependency": len(rows),
            "low_component_dimension_cap": 22,
            "low_nonexceptional_maximum_absolute_bias": low_nonexceptional_maximum,
            "remaining_scope": "all nonexceptional exact component supports of dimension at least 23",
        },
        "scope_limit": "This replay authenticates the spectral implication, not the remaining high-support lemma.",
        "result": "PASS_INDEPENDENT_SPECTRAL_INTERFACE_REPLAY",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")
    print(f"bias_cap={bias_limit}")
    print(f"triangle_upper_ceiling={upper_ceiling}")
    print("status=PASS_INDEPENDENT_SPECTRAL_INTERFACE_REPLAY")


if __name__ == "__main__":
    main()
