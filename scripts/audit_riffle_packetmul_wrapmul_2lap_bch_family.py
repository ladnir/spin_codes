#!/usr/bin/env python3
"""Audit the exact small-BCH lifted-orbit family experiment."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_riffle_dp_g4_autonomous_map import (  # noqa: E402
    polynomial_gcd,
    polynomial_lcm,
    polynomial_mod,
)
from bch_candidate_params import bch_dimension  # noqa: E402
from check_bch_boundary_smallfield import find_primitive_poly  # noqa: E402
from check_bch_generator_inner import systematic_basis  # noqa: E402
from exact_bch_spectrum_small import (  # noqa: E402
    enumerate_spectrum,
    nullspace_basis,
    parity_rows,
)


CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PRIMARY_SOURCE = ROOT / "scripts" / "analyze_riffle_packetmul_wrapmul_2lap_bch32_family.cpp"
PRIMARY_EXE = ROOT / "scripts" / "analyze_riffle_packetmul_wrapmul_2lap_bch32_family.exe"
PRIMARY_RAW = CANDIDATE / "receipts" / "goal07_bch32_family_exact_raw.json"
CURRENT_ALGEBRA = ROOT / "explorations" / "riffle_dp_g4_g2_autonomous_map.json"
GOAL05_COMPONENTS = CANDIDATE / "receipts" / "goal05_primary_components.json"
GOAL06_AUDIT = CANDIDATE / "receipts" / "goal06_amortized_route_audit.json"
OUTPUT = CANDIDATE / "receipts" / "goal07_bch_family_audit.json"
CURRENT_GENERATOR = 0xF4845518B9582A1F


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
    return result


def accumulate(value: int, width: int) -> int:
    mask = (1 << width) - 1
    shift = 1
    while shift < width:
        value ^= (value << shift) & mask
        shift *= 2
    return value


def build_instance(m: int, delta: int) -> dict:
    n, dimension, redundancy = bch_dimension(m, delta)
    if dimension * 2 != n + 1:
        raise RuntimeError("family instance is not half rate after extension")
    primitive = find_primitive_poly(m)
    cyclic_basis = nullspace_basis(parity_rows(m, delta, primitive), n)
    extended_basis = [
        word | ((word.bit_count() & 1) << n)
        for word in cyclic_basis
    ]
    basis = systematic_basis(extended_basis, n + 1, list(range(dimension)))
    columns = tuple(word >> dimension for word in basis)
    if any((word & ((1 << dimension) - 1)) != 1 << bit for bit, word in enumerate(basis)):
        raise RuntimeError("systematic BCH reconstruction failed")
    return {
        "m": m,
        "delta": delta,
        "n": n,
        "width": dimension,
        "primitive": primitive,
        "cyclic_basis": cyclic_basis,
        "extended_basis": extended_basis,
        "systematic_basis": basis,
        "columns": columns,
    }


def observe(state: int, width: int, columns: tuple[int, ...], nodes: int = 24) -> tuple[int, ...]:
    mask = (1 << width) - 1
    a = state & mask
    b = state >> width
    outputs = []
    for _ in range(nodes):
        acc_a = accumulate(a, width)
        acc_b = accumulate(b, width)
        acc2_a = accumulate(acc_a, width)
        outputs.append(acc2_a ^ acc_b)
        a = apply_columns(columns, acc_a)
        b = apply_columns(columns, acc2_a) ^ apply_columns(columns, acc_b)
    return tuple(outputs)


def exact_orbit_profile(width: int, columns: tuple[int, ...], nodes: int = 24) -> dict:
    prefixes = [[10**9, 0, 0] for _ in range(nodes)]
    anchors = [[[10**9, 0, 0] for _ in range(width + 1)] for _ in range(nodes)]
    spectrum: Counter[int] = Counter()
    for state in range(1, 1 << (2 * width)):
        outputs = observe(state, width, columns, nodes)
        weights = [value.bit_count() for value in outputs]
        total = 0
        for node, weight in enumerate(weights):
            total += weight
            row = prefixes[node]
            if total < row[0]:
                row[:] = [total, 1, state]
            elif total == row[0]:
                row[1] += 1
                row[2] = min(row[2], state)
        spectrum[total] += 1
        for node, anchor_weight in enumerate(weights):
            row = anchors[node][anchor_weight]
            if total < row[0]:
                row[:] = [total, 1, state]
            elif total == row[0]:
                row[1] += 1
                row[2] = min(row[2], state)
    return {
        "prefix_minima": [
            {"weight": weight, "count": count, "witness_hex": hex(witness)}
            for weight, count, witness in prefixes
        ],
        "weight24_spectrum": {str(weight): spectrum[weight] for weight in sorted(spectrum)},
        "anchor_minima": [
            [
                None
                if count == 0
                else {"weight": weight, "count": count, "witness_hex": hex(witness)}
                for weight, count, witness in node
            ]
            for node in anchors
        ],
    }


def kernel_basis_from_columns(columns: tuple[int, ...]) -> tuple[tuple[int, ...], int]:
    pivot_values: dict[int, int] = {}
    pivot_representations: dict[int, int] = {}
    kernel = []
    for column_index, column in enumerate(columns):
        value = column
        representation = 1 << column_index
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivot_values:
                value ^= pivot_values[pivot]
                representation ^= pivot_representations[pivot]
            else:
                pivot_values[pivot] = value
                pivot_representations[pivot] = representation
                break
        if value == 0:
            kernel.append(representation)
    return tuple(kernel), len(pivot_values)


def even_observation_ranks(instance: dict, samples: int) -> list[int]:
    width = instance["width"]
    observations = [
        observe(1 << bit, width, instance["columns"], 2 * samples)
        for bit in range(2 * width)
    ]
    result = []
    for count in range(1, samples + 1):
        columns = tuple(
            sum(row[2 * index] << (width * index) for index in range(count))
            for row in observations
        )
        _, rank = kernel_basis_from_columns(columns)
        result.append(rank)
    return result


def exact_even_zero_profile(instance: dict, nodes: int) -> dict:
    width = instance["width"]
    observations = [
        observe(1 << bit, width, instance["columns"], nodes)
        for bit in range(2 * width)
    ]
    even_columns = tuple(
        sum(row[node] << (width * (node // 2)) for node in range(0, nodes, 2))
        for row in observations
    )
    full_columns = tuple(
        sum(value << (width * node) for node, value in enumerate(row))
        for row in observations
    )
    kernel, rank = kernel_basis_from_columns(even_columns)
    if rank + len(kernel) != 2 * width:
        raise RuntimeError("even-observation rank-nullity mismatch")

    generators = []
    for state in kernel:
        even_word = apply_columns(even_columns, state)
        if even_word:
            raise RuntimeError("invalid even-observation kernel vector")
        generators.append((state, apply_columns(full_columns, state)))

    if not generators:
        return {
            "nodes": nodes,
            "even_samples": (nodes + 1) // 2,
            "rank": rank,
            "kernel_dimension": 0,
            "nonzero_states": 0,
            "minimum_weight": None,
            "minimum_count": 0,
            "witness_hex": None,
        }

    current_state = 0
    current_word = 0
    previous_gray = 0
    minimum = 10**9
    count = 0
    witness = 0
    for index in range(1, 1 << len(generators)):
        gray = index ^ (index >> 1)
        bit = (gray ^ previous_gray).bit_length() - 1
        current_state ^= generators[bit][0]
        current_word ^= generators[bit][1]
        weight = current_word.bit_count()
        if weight < minimum:
            minimum = weight
            count = 1
            witness = current_state
        elif weight == minimum:
            count += 1
            witness = min(witness, current_state)
        previous_gray = gray

    witness_outputs = observe(witness, width, instance["columns"], nodes)
    if any(witness_outputs[node] for node in range(0, nodes, 2)):
        raise RuntimeError("even-zero witness replay has a nonzero even output")
    if sum(value.bit_count() for value in witness_outputs) != minimum:
        raise RuntimeError("even-zero witness weight replay failed")
    return {
        "nodes": nodes,
        "even_samples": (nodes + 1) // 2,
        "rank": rank,
        "kernel_dimension": len(kernel),
        "nonzero_states": (1 << len(kernel)) - 1,
        "minimum_weight": minimum,
        "minimum_count": count,
        "witness_hex": hex(witness),
        "witness_node_weights": [value.bit_count() for value in witness_outputs],
    }


def krylov_polynomial(step, start: int, width: int) -> int:
    basis_values = [0] * width
    basis_representations = [0] * width
    value = start
    for exponent in range(width + 1):
        reduced = value
        representation = 1 << exponent
        while reduced:
            pivot = reduced.bit_length() - 1
            if basis_values[pivot]:
                reduced ^= basis_values[pivot]
                representation ^= basis_representations[pivot]
            else:
                basis_values[pivot] = reduced
                basis_representations[pivot] = representation
                break
        if reduced == 0:
            return representation
        value = step(value)
    raise RuntimeError("Krylov relation not found")


def minimum_polynomial(width: int, columns: tuple[int, ...]) -> int:
    def step(value: int) -> int:
        return apply_columns(columns, accumulate(value, width))

    result = 1
    for bit in range(width):
        result = polynomial_lcm(result, krylov_polynomial(step, 1 << bit, width))
    for bit in range(width):
        value = 1 << bit
        check = 0
        for exponent in range(result.bit_length()):
            if (result >> exponent) & 1:
                check ^= value
            value = step(value)
        if check:
            raise RuntimeError("reported minimum polynomial does not annihilate T")
    return result


def lifted_step(state: int, width: int, columns: tuple[int, ...]) -> int:
    mask = (1 << width) - 1
    a = state & mask
    b = state >> width
    acc_a = accumulate(a, width)
    acc_b = accumulate(b, width)
    acc2_a = accumulate(acc_a, width)
    next_a = apply_columns(columns, acc_a)
    next_b = apply_columns(columns, acc2_a) ^ apply_columns(columns, acc_b)
    return next_a | (next_b << width)


def apply_polynomial(step, polynomial: int, value: int) -> int:
    result = 0
    current = value
    for exponent in range(polynomial.bit_length()):
        if (polynomial >> exponent) & 1:
            result ^= current
        current = step(current)
    return result


def lifted_algebra(width: int, columns: tuple[int, ...], witness: int) -> dict:
    step = lambda state: lifted_step(state, width, columns)
    polynomial = 1
    for bit in range(2 * width):
        polynomial = polynomial_lcm(
            polynomial,
            krylov_polynomial(step, 1 << bit, 2 * width),
        )
    factors = factor_small_polynomial(polynomial)
    powered_factors = []
    for row in factors:
        factor = int(row["factor_hex"], 16)
        powered = 1
        for _ in range(row["multiplicity"]):
            next_powered = 0
            right = factor
            left = powered
            while right:
                if right & 1:
                    next_powered ^= left
                left <<= 1
                right >>= 1
            powered = next_powered
        powered_factors.append(powered)
    support = []
    for row, powered in zip(factors, powered_factors, strict=True):
        complementary = polynomial_divide_exact(polynomial, powered)
        support.append({
            **row,
            "witness_component_nonzero": apply_polynomial(step, complementary, witness) != 0,
        })
    return {
        "minimal_polynomial_hex": hex(polynomial),
        "factors": support,
    }


def random_code_crossing(length: int, dimension: int) -> int:
    cumulative = 0
    for weight in range(length + 1):
        cumulative += math.comb(length, weight)
        if math.log2(cumulative) + dimension - length >= 0:
            return weight
    raise RuntimeError("random-code crossing not found")


def random_expected_log2(length: int, dimension: int, cutoff: int) -> float:
    cumulative = sum(math.comb(length, weight) for weight in range(cutoff + 1))
    return math.log2(cumulative) + dimension - length


def gap_accounting(nodes: int) -> dict:
    zero_nodes = 32737
    gaps = 34
    required_output = 188765
    complete_windows = -(-(zero_nodes - gaps * (nodes - 1)) // nodes)
    minimum_bound = required_output // complete_windows + 1
    clean_six_bound = 6 * nodes
    return {
        "nodes": nodes,
        "complete_windows": complete_windows,
        "minimum_sufficient_bound": minimum_bound,
        "clean_six_per_node_bound": clean_six_bound,
        "clean_six_per_node_suffices": clean_six_bound >= minimum_bound,
        "clean_six_per_node_margin": complete_windows * clean_six_bound - required_output,
        "forced_anchor_maximum_for_clean_bound": (clean_six_bound - 1) // nodes,
    }


def polynomial_divide_exact(dividend: int, divisor: int) -> int:
    quotient = 0
    divisor_degree = divisor.bit_length() - 1
    while dividend.bit_length() - 1 >= divisor_degree:
        shift = dividend.bit_length() - 1 - divisor_degree
        quotient ^= 1 << shift
        dividend ^= divisor << shift
    if dividend:
        raise RuntimeError("non-exact polynomial division")
    return quotient


def polynomial_multiply_mod(left: int, right: int, modulus: int) -> int:
    result = 0
    degree = modulus.bit_length() - 1
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left.bit_length() - 1 == degree:
            left ^= modulus
    return result


def polynomial_power_mod(value: int, exponent: int, modulus: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_multiply_mod(result, value, modulus)
        value = polynomial_multiply_mod(value, value, modulus)
        exponent >>= 1
    return result


def irreducible(polynomial: int) -> bool:
    degree = polynomial.bit_length() - 1
    x_value = polynomial_mod(2, polynomial)
    x_power = x_value
    for _ in range(1, degree // 2 + 1):
        x_power = polynomial_multiply_mod(x_power, x_power, polynomial)
        if polynomial_gcd(x_power ^ x_value, polynomial) != 1:
            return False
    return polynomial_power_mod(x_value, 1 << degree, polynomial) == x_value


def factor_small_polynomial(polynomial: int) -> list[dict]:
    remainder = polynomial
    result = []
    max_degree = (polynomial.bit_length() - 1) // 2
    for degree in range(1, max_degree + 1):
        for low in range(1, 1 << degree, 2):
            factor = (1 << degree) | low
            if not irreducible(factor):
                continue
            multiplicity = 0
            while polynomial_mod(remainder, factor) == 0:
                remainder = polynomial_divide_exact(remainder, factor)
                multiplicity += 1
            if multiplicity:
                result.append({
                    "factor_hex": hex(factor),
                    "degree": degree,
                    "multiplicity": multiplicity,
                })
    if remainder != 1:
        if not irreducible(remainder):
            raise RuntimeError("small polynomial factorization incomplete")
        result.append({
            "factor_hex": hex(remainder),
            "degree": remainder.bit_length() - 1,
            "multiplicity": 1,
        })
    return result


def bch_spectrum(instance: dict) -> dict:
    spectrum = enumerate_spectrum(instance["extended_basis"], instance["n"] + 1)
    split: Counter[tuple[int, int]] = Counter()
    basis = instance["systematic_basis"]
    width = instance["width"]
    for message in range(1 << width):
        word = apply_columns(tuple(basis), message)
        split[(message.bit_count(), (word >> width).bit_count())] += 1
    return {
        "ordinary": {str(weight): count for weight, count in enumerate(spectrum) if count},
        "systematic_split": {
            f"{left},{right}": count
            for (left, right), count in sorted(split.items())
        },
    }


def verify_primary(raw: dict, instance: dict) -> dict:
    if raw["result"] != "EXHAUSTED":
        raise RuntimeError("primary enumeration is not exhausted")
    if raw["enumeration"]["nonzero_states"] != (1 << 32) - 1:
        raise RuntimeError("primary state count mismatch")
    if int(raw["instance"]["generator_polynomial_hex"], 16) != instance["columns"][0]:
        raise RuntimeError("primary generator polynomial mismatch")
    raw_columns = tuple(int(value, 16) for value in raw["instance"]["parity_columns_hex"])
    if raw_columns != instance["columns"]:
        raise RuntimeError("primary systematic parity columns mismatch")
    if sum(raw["weight24_spectrum"].values()) != (1 << 32) - 1:
        raise RuntimeError("primary spectrum count mismatch")

    nonzero_spectrum = sorted((int(weight), count) for weight, count in raw["weight24_spectrum"].items())
    global_minimum = raw["prefix_minima"][-1]
    if nonzero_spectrum[0] != (global_minimum["weight"], global_minimum["count"]):
        raise RuntimeError("primary minimum disagrees with spectrum")

    for node, minimum in enumerate(raw["prefix_minima"]):
        state = int(minimum["witness_hex"], 16)
        weights = [value.bit_count() for value in observe(state, 16, instance["columns"])]
        if sum(weights[: node + 1]) != minimum["weight"]:
            raise RuntimeError(f"prefix witness replay failed at node {node}")

    for node, rows in enumerate(raw["anchor_minima"]):
        for anchor_weight, minimum in enumerate(rows):
            if minimum is None:
                continue
            state = int(minimum["witness_hex"], 16)
            weights = [value.bit_count() for value in observe(state, 16, instance["columns"])]
            if weights[node] != anchor_weight or sum(weights) != minimum["weight"]:
                raise RuntimeError(f"anchor witness replay failed at ({node},{anchor_weight})")

    global_state = int(global_minimum["witness_hex"], 16)
    global_weights = [value.bit_count() for value in observe(global_state, 16, instance["columns"])]
    if global_weights != raw["global_witness_node_weights"]:
        raise RuntimeError("global witness weight vector mismatch")
    return {
        "d24": global_minimum["weight"],
        "d24_count": global_minimum["count"],
        "d24_witness_hex": global_minimum["witness_hex"],
        "d24_witness_node_weights": global_weights,
        "states_below_scaled_target_36": sum(
            count for weight, count in nonzero_spectrum if weight < 36
        ),
        "zero_anchor_d24_range": [
            min(raw["anchor_minima"][node][0]["weight"] for node in range(24)),
            max(raw["anchor_minima"][node][0]["weight"] for node in range(24)),
        ],
        "weight_one_anchor_d24_range": [
            min(raw["anchor_minima"][node][1]["weight"] for node in range(24)),
            max(raw["anchor_minima"][node][1]["weight"] for node in range(24)),
        ],
    }


def main() -> None:
    family = [build_instance(3, 3), build_instance(5, 7), build_instance(7, 21)]
    if family[0]["columns"] != (0xB, 0xE, 0x7, 0xD):
        raise RuntimeError("EBCH [8,4,4] parity map mismatch")
    if family[1]["columns"][0] != 0x8FAF:
        raise RuntimeError("EBCH [32,16,8] generator mismatch")

    current_rows = tuple((CURRENT_GENERATOR << row) | (1 << 127) for row in range(64))
    current_systematic = systematic_basis(list(current_rows), 128, list(range(64)))
    current_columns = tuple(word >> 64 for word in current_systematic)
    if family[2]["columns"] != current_columns:
        raise RuntimeError("generic BCH family does not reproduce the committed parity map")

    raw = json.loads(PRIMARY_RAW.read_text())
    primary_summary = verify_primary(raw, family[1])
    small_profile = exact_orbit_profile(4, family[0]["columns"])

    algebra = []
    for instance in family[:2]:
        polynomial = minimum_polynomial(instance["width"], instance["columns"])
        algebra.append({
            "width": instance["width"],
            "minimal_polynomial_hex": hex(polynomial),
            "factors": factor_small_polynomial(polynomial),
        })

    current_algebra = json.loads(CURRENT_ALGEBRA.read_text())
    goal05_components = json.loads(GOAL05_COMPONENTS.read_text())
    goal06_audit = json.loads(GOAL06_AUDIT.read_text())
    current_factors = [
        {
            "factor_hex": row["factor_hex"],
            "degree": row["degree"],
            "multiplicity": 1,
        }
        for row in current_algebra["irreducible_components"]
    ]
    algebra.append({
        "width": 64,
        "minimal_polynomial_hex": current_algebra["minimal_polynomial_hex"],
        "factors": current_factors,
        "source_sha256": digest(CURRENT_ALGEBRA),
    })

    small_d24 = small_profile["prefix_minima"][-1]
    small_witness = int(small_d24["witness_hex"], 16)
    medium_witness = int(primary_summary["d24_witness_hex"], 16)
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-bch-family-audit-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SMALL_INSTANCE_WITH_SOURCE_AND_WITNESS_AUDIT",
        "source_sha256": digest(Path(__file__).resolve()),
        "primary_source_sha256": digest(PRIMARY_SOURCE),
        "primary_executable_sha256": digest(PRIMARY_EXE),
        "primary_raw_sha256": digest(PRIMARY_RAW),
        "family": [
            {
                "m": instance["m"],
                "primitive_polynomial_hex": hex(instance["primitive"]),
                "designed_distance": instance["delta"],
                "cyclic_length": instance["n"],
                "dimension": instance["width"],
                "extended_length": instance["n"] + 1,
                "generator_polynomial_hex": hex(instance["columns"][0]),
                "parity_columns_sha256": hashlib.sha256(
                    b"".join(
                        value.to_bytes((instance["width"] + 7) // 8, "little")
                        for value in instance["columns"]
                    )
                ).hexdigest(),
            }
            for instance in family
        ],
        "committed_width64_map_reproduced": True,
        "bch_spectra": {
            "width4": bch_spectrum(family[0]),
            "width16": bch_spectrum(family[1]),
        },
        "autonomous_map_algebra": algebra,
        "lifted_map_algebra": {
            "width4": lifted_algebra(4, family[0]["columns"], small_witness),
            "width16": lifted_algebra(16, family[1]["columns"], medium_witness),
        },
        "width4_exact": {
            "nonzero_states": (1 << 8) - 1,
            "scaled_target": 9,
            "d24": small_d24["weight"],
            "d24_count": small_d24["count"],
            "d24_witness_hex": small_d24["witness_hex"],
            "d24_witness_node_weights": [
                value.bit_count()
                for value in observe(int(small_d24["witness_hex"], 16), 4, family[0]["columns"])
            ],
            "prefix_minima": small_profile["prefix_minima"],
            "weight24_spectrum": small_profile["weight24_spectrum"],
            "anchor_minima": small_profile["anchor_minima"],
        },
        "width16_exact": {
            "scaled_target": 36,
            **primary_summary,
            "prefix_minima": raw["prefix_minima"],
        },
        "comparison_to_width64": {
            "goal05_components_sha256": digest(GOAL05_COMPONENTS),
            "goal06_audit_sha256": digest(GOAL06_AUDIT),
            "required_d24": goal06_audit["amortized_target"]["proposed_window_lower_bound"],
            "known_d24_upper_bound": goal06_audit["goal05_counterexample_extension"]["first_24_weight"],
            "known_nine_node_transient_weight": goal05_components["known_weight_44_witness"]["first_9_weight"],
            "known_nine_node_transient_component_support_mask_hex": (
                goal05_components["known_weight_44_witness"]["component_support_mask_hex"]
            ),
            "known_nine_node_transient_uses_all_primary_components": (
                goal05_components["known_weight_44_witness"]["component_support_mask_hex"] == "0x7f"
            ),
            "random_code_baselines": [
                {
                    "width": width,
                    "observation_code": f"[{24 * width},{2 * width}]",
                    "random_expected_count_crossing_weight": random_code_crossing(24 * width, 2 * width),
                    "observed_or_known_upper_bound": observed,
                    "random_expected_log2_count_through_observed_weight": random_expected_log2(
                        24 * width, 2 * width, observed
                    ),
                }
                for width, observed in ((4, small_d24["weight"]), (16, primary_summary["d24"]), (64, 278))
            ],
            "even_observation_rank_sequences": {
                "width4_first_12_samples": even_observation_ranks(family[0], 12),
                "width16_first_12_samples": even_observation_ranks(family[1], 12),
                "width64_first_23_samples": even_observation_ranks(family[2], 23),
            },
            "width64_exact_even_zero_profiles": [
                exact_even_zero_profile(family[2], nodes)
                for nodes in (32, 36, 38, 39, 40, 44, 45)
            ],
            "window_accounting": [
                gap_accounting(nodes)
                for nodes in (24, 32, 36, 38, 39, 40, 44, 45)
            ],
        },
        "scope": (
            "The audit proves the family reconstruction, exhausts all nonzero lifted states "
            "at widths 4 and 16, and replays every width-16 minimum witness with an "
            "independently reconstructed parity map. It does not extrapolate either "
            "small-instance distance to width 64."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "status": "PASS",
        "width4_d24": payload["width4_exact"]["d24"],
        "width16_d24": payload["width16_exact"]["d24"],
        "output": str(OUTPUT.relative_to(ROOT)),
    }, indent=2))


if __name__ == "__main__":
    main()
