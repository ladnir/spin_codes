#!/usr/bin/env python3
"""Exact algebraic and pure-component audit for PacketMul zero symbols."""

from __future__ import annotations

import collections
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import (
    accumulate,
    apply_polynomial,
    build_apply,
    kernel_basis,
)
from probe_riffle_dp_g4_component_mixing import (
    FACTORS,
    INNER_NODES,
    PACKET_SLOTS,
    SAMPLE_SPACE,
    project,
    projection_tables,
    transpose_columns,
    walsh_transform,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
MANIFEST = CANDIDATE / "manifest.json"
OUTPUT = CANDIDATE / "receipts" / "goal01_zero_symbol_primary.json"
FIELD_MODULUS = 0x13
STATE_BITS = 64
COMPONENT_DEGREES = (1, 2, 4, 9, 10, 18, 20)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def gf16_multiply(left: int, right: int) -> int:
    result = 0
    for _ in range(4):
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 0x10:
            left ^= FIELD_MODULUS
    return result & 15


def gf16_trace(value: int) -> int:
    result = value
    power = value
    for _ in range(3):
        power = gf16_multiply(power, power)
        result ^= power
    if result not in (0, 1):
        raise RuntimeError("GF(16) trace did not land in GF(2)")
    return result


def trace_dual_table() -> tuple[int, ...]:
    """Map a raw four-bit linear functional to its trace coefficient."""
    functional_to_coefficient: dict[int, int] = {}
    for coefficient in range(16):
        functional = sum(
            gf16_trace(gf16_multiply(coefficient, 1 << bit)) << bit
            for bit in range(4)
        )
        if functional in functional_to_coefficient:
            raise RuntimeError("trace pairing is degenerate")
        functional_to_coefficient[functional] = coefficient
    if len(functional_to_coefficient) != 16:
        raise RuntimeError("trace pairing did not enumerate all functionals")
    return tuple(functional_to_coefficient[functional] for functional in range(16))


def basis_contributions(step) -> np.ndarray:
    """Return z_{1,alpha,alpha^2,alpha^3}(p), with slot-major cell order."""
    result = np.empty((4, SAMPLE_SPACE), dtype=np.uint64)
    for basis_bit in range(4):
        cursor = 0
        for slot in range(PACKET_SLOTS):
            state = (1 << basis_bit) << (4 * slot)
            for _ in range(INNER_NODES):
                state = step(state)
                result[basis_bit, cursor] = state
                cursor += 1
    return result


def value_contributions(basis_states: np.ndarray, value: int) -> np.ndarray:
    result = np.zeros(SAMPLE_SPACE, dtype=np.uint64)
    for bit in range(4):
        if (value >> bit) & 1:
            result ^= basis_states[bit]
    return result


def combine_basis(basis: tuple[int, ...], index: int) -> int:
    result = 0
    while index:
        bit = index & -index
        result ^= basis[bit.bit_length() - 1]
        index ^= bit
    return result


def coefficient_histogram(basis_states: np.ndarray, character: int) -> list[int]:
    raw = np.zeros(SAMPLE_SPACE, dtype=np.uint8)
    for basis_bit in range(4):
        parities = np.fromiter(
            ((character & int(state)).bit_count() & 1 for state in basis_states[basis_bit]),
            dtype=np.uint8,
            count=SAMPLE_SPACE,
        )
        raw |= parities << basis_bit
    dual = np.asarray(trace_dual_table(), dtype=np.uint8)
    coefficients = dual[raw]
    return np.bincount(coefficients, minlength=16).astype(int).tolist()


def fraction_row(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest["candidate_id"] != "riffle_packetmul_2lap_g4":
        raise RuntimeError("active candidate manifest changed")
    if manifest["active_goal"] != "GOAL_01_ZERO_SYMBOL_GATE.md":
        raise RuntimeError("Goal 01 is no longer active")

    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    t_columns = tuple(step(1 << bit) for bit in range(STATE_BITS))
    if len(set(t_columns)) != STATE_BITS or 0 in t_columns:
        raise RuntimeError("autonomous map is not invertible")

    basis_states = basis_contributions(step)
    basis_digest = sha256_bytes(basis_states.astype("<u8", copy=False).tobytes())

    # Verify the field-linear reconstruction at every value and cell.
    value_states: list[np.ndarray] = []
    value_digests: dict[str, str] = {}
    for value in range(1, 16):
        states = value_contributions(basis_states, value)
        value_states.append(states)
        value_digests[str(value)] = sha256_bytes(
            states.astype("<u8", copy=False).tobytes()
        )
    for left in range(16):
        for right in range(16):
            if left == right:
                continue
            lhs = value_contributions(basis_states, left ^ right)
            rhs = value_contributions(basis_states, left) ^ value_contributions(
                basis_states, right
            )
            if not np.array_equal(lhs, rhs):
                raise RuntimeError("packet-value linearity failed")

    # The PacketMul sample space includes all 15 values. Count its collisions.
    all_states = np.concatenate(value_states)
    _unique, multiplicities = np.unique(all_states, return_counts=True)
    multiplicity_histogram = collections.Counter(map(int, multiplicities))
    collision_sum = sum(
        multiplicity * multiplicity * frequency
        for multiplicity, frequency in multiplicity_histogram.items()
    )
    labeled_sample_count = SAMPLE_SPACE * 15
    parseval = Fraction((1 << STATE_BITS) * collision_sum, labeled_sample_count**2)
    collision_row = {
        "labeled_samples": labeled_sample_count,
        "distinct_states": int(len(multiplicities)),
        "maximum_multiplicity": int(multiplicities.max()),
        "multiplicity_histogram": {
            str(key): value for key, value in sorted(multiplicity_histogram.items())
        },
        "ordered_collision_count": str(collision_sum),
        "normalized_parseval_sum": fraction_row(parseval),
        "identity": (
            "sum_chi beta_chi^2 = 2^64 * sum_x m_x^2 / (15*M)^2, "
            "where m_x is the multiplicity of x under (p,y) -> z_y(p)."
        ),
    }

    transpose_apply = build_apply(transpose_columns(t_columns))
    component_rows = []
    reported_characters: set[int] = set()
    for degree in COMPONENT_DEGREES:
        annihilator_columns = tuple(
            apply_polynomial(transpose_apply, FACTORS[degree], 1 << bit)
            for bit in range(STATE_BITS)
        )
        basis = kernel_basis(annihilator_columns)
        if len(basis) != degree:
            raise RuntimeError(f"degree-{degree} component dimension mismatch")
        tables = projection_tables(basis)
        zero_symbol_numerators = np.full(1 << degree, SAMPLE_SPACE, dtype=np.int64)
        for states in value_states:
            syndromes = project(states, tables)
            histogram = np.bincount(syndromes, minlength=1 << degree)
            walsh = walsh_transform(histogram)
            if int(walsh[0]) != SAMPLE_SPACE:
                raise RuntimeError("component Walsh mass mismatch")
            zero_symbol_numerators += walsh
        if np.any(zero_symbol_numerators % 16):
            raise RuntimeError("zero-symbol Walsh inversion is nonintegral")
        zero_counts = zero_symbol_numerators // 16
        if int(zero_counts[0]) != SAMPLE_SPACE:
            raise RuntimeError("trivial character zero-symbol count mismatch")

        nontrivial = zero_counts[1:]
        max_q_index = int(np.argmax(nontrivial)) + 1
        beta_numerators = 16 * nontrivial - SAMPLE_SPACE
        max_abs_index = int(np.argmax(np.abs(beta_numerators))) + 1
        max_q_character = combine_basis(basis, max_q_index)
        max_abs_character = combine_basis(basis, max_abs_index)
        reported_characters.update((max_q_character, max_abs_character))
        component_rows.append(
            {
                "degree": degree,
                "factor_hex": hex(FACTORS[degree]),
                "basis_hex": [hex(value) for value in basis],
                "characters_exhausted": (1 << degree) - 1,
                "maximum_zero_symbol": {
                    "count": int(zero_counts[max_q_index]),
                    "fraction": fraction_row(
                        Fraction(int(zero_counts[max_q_index]), SAMPLE_SPACE)
                    ),
                    "component_index": max_q_index,
                    "character_hex": hex(max_q_character),
                },
                "maximum_absolute_beta": {
                    "absolute_numerator": int(abs(16 * zero_counts[max_abs_index] - SAMPLE_SPACE)),
                    "denominator": 15 * SAMPLE_SPACE,
                    "fraction": fraction_row(
                        Fraction(
                            int(abs(16 * zero_counts[max_abs_index] - SAMPLE_SPACE)),
                            15 * SAMPLE_SPACE,
                        )
                    ),
                    "signed_numerator": int(16 * zero_counts[max_abs_index] - SAMPLE_SPACE),
                    "component_index": max_abs_index,
                    "character_hex": hex(max_abs_character),
                },
            }
        )

    replay_rows = []
    for character in sorted(reported_characters):
        histogram = coefficient_histogram(basis_states, character)
        zero_count = histogram[0]
        replay_rows.append(
            {
                "character_hex": hex(character),
                "trace_coefficient_histogram_0_through_15": histogram,
                "zero_symbol_count": zero_count,
                "signed_beta_numerator": 16 * zero_count - SAMPLE_SPACE,
                "beta_denominator": 15 * SAMPLE_SPACE,
            }
        )

    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal01-zero-symbol-primary-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": manifest["candidate_id"],
        "active_manifest_sha256": sha256_file(MANIFEST),
        "evidence_label": "EXACT_PRIMARY",
        "packet_cells": SAMPLE_SPACE,
        "cell_order": "slot 0..15, then exponent 1..32772",
        "autonomous_map": "T(s)=P(accumulate(s)); z_y(slot,t)=T^t(y << (4*slot))",
        "field_modulus_hex": hex(FIELD_MODULUS),
        "trace_dual_table_raw_functional_to_field_coefficient": list(trace_dual_table()),
        "basis_contribution_sha256_little_endian_u64": basis_digest,
        "value_contribution_sha256_little_endian_u64": value_digests,
        "coefficient_identity": (
            "The four parities <chi,z_(1<<b)(p)> form a raw linear functional. "
            "The displayed trace-dual table converts that functional to the unique "
            "a_chi(p) satisfying <chi,z_y(p)>=Tr(a_chi(p)*y)."
        ),
        "packetmul_bias_identity": (
            "A fixed nonzero source packet times a uniform GF(16)^* scalar is "
            "uniform on y!=0. Hence beta_chi=(16*q_chi-1)/15."
        ),
        "cross_value_collisions": collision_row,
        "pure_component_rows": component_rows,
        "primary_maximizer_replays": replay_rows,
        "scope_limitation": (
            "The component spectra are exhaustive only inside each pure irreducible "
            "component. Mixed characters remain subject to a separate refutation search "
            "and proof obligation."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"basis_contribution_sha256={basis_digest}")
    print(f"cross_value_distinct_states={collision_row['distinct_states']}")
    print(f"cross_value_maximum_multiplicity={collision_row['maximum_multiplicity']}")
    for row in component_rows:
        print(
            f"degree={row['degree']} max_q={row['maximum_zero_symbol']['fraction']} "
            f"max_abs_beta={row['maximum_absolute_beta']['fraction']}"
        )
    print(f"output={OUTPUT}")
    print("status=EXACT_PRIMARY_ZERO_SYMBOL_GATE")


if __name__ == "__main__":
    main()
