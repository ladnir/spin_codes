#!/usr/bin/env python3
"""High-precision sparse checks for TransposePacketShuffle-RandomStepConv.

The certificate evaluates the one- and two-active-block packed rows at fixed
decimal Chernoff parameters.  Interval arithmetic encloses every matrix and
logarithm operation.  A fixed parameter is sufficient for a valid Chernoff
upper bound; no numerical optimizer is part of the certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp


DEFAULT_OUTPUT = Path(
    "constructions/riffle_transpose_packetshuffle_randomstepconv/"
    "receipts/g4_b1024_sigma15_sparse_interval.json"
)


def fraction_count_distribution(occupancies: tuple[int, ...]) -> list[Fraction]:
    distribution = [Fraction(1)]
    for rank in occupancies:
        zero = Fraction(1, 1 << rank)
        active = 1 - zero
        updated = [Fraction(0)] * (len(distribution) + 1)
        for count, probability in enumerate(distribution):
            updated[count] += probability * zero
            updated[count + 1] += probability * active
        distribution = updated
    return distribution


def weak_compositions(total: int, parts: int, cap: int):
    if parts == 1:
        if total <= cap:
            yield (total,)
        return
    for first in range(min(cap, total) + 1):
        for suffix in weak_compositions(total - first, parts - 1, cap):
            yield (first,) + suffix


def packed_occupancies(total: int, parts: int, cap: int) -> tuple[int, ...]:
    full, remainder = divmod(total, cap)
    return (cap,) * full + ((remainder,) if remainder else ()) + (0,) * (
        parts - full - (1 if remainder else 0)
    )


def packing_self_test() -> dict[str, int]:
    checked = 0
    for cap in range(2, 6):
        for parts in range(2, 5):
            for total in range(1, cap * parts + 1):
                packed = fraction_count_distribution(
                    packed_occupancies(total, parts, cap)
                )
                packed_tails = [sum(packed[k:]) for k in range(len(packed))]
                for occupancy in weak_compositions(total, parts, cap):
                    current = fraction_count_distribution(occupancy)
                    current_tails = [
                        sum(current[k:]) for k in range(len(current))
                    ]
                    if any(
                        packed_tail > current_tail
                        for packed_tail, current_tail in zip(
                            packed_tails, current_tails
                        )
                    ):
                        raise AssertionError("packed active count is not dominated")
                    checked += 1
    return {"occupancy_vectors_checked": checked}


def output_weight_distribution(
    support: tuple[bool, ...], packet_bits: int, sigma: int
) -> list[Fraction]:
    distribution = {(0, 0): Fraction(1)}
    state_zero_probability = Fraction(1, 1 << sigma)
    for input_active in support:
        updated: dict[tuple[int, int], Fraction] = {}
        for (state_live, weight), probability in distribution.items():
            if not input_active and not state_live:
                key = (0, weight)
                updated[key] = updated.get(key, Fraction(0)) + probability
                continue
            for output_weight in range(packet_bits + 1):
                output_probability = Fraction(
                    math.comb(packet_bits, output_weight), 1 << packet_bits
                )
                for next_live, state_probability in (
                    (0, state_zero_probability),
                    (1, 1 - state_zero_probability),
                ):
                    key = (next_live, weight + output_weight)
                    updated[key] = updated.get(key, Fraction(0)) + (
                        probability * output_probability * state_probability
                    )
        distribution = updated
    weights = [Fraction(0)] * (packet_bits * len(support) + 1)
    for (_state_live, weight), probability in distribution.items():
        weights[weight] += probability
    return weights


def support_monotonicity_self_test() -> dict[str, int]:
    positions = 5
    checked = 0
    for categories in itertools.product(range(3), repeat=positions):
        smaller = tuple(category == 2 for category in categories)
        larger = tuple(category >= 1 for category in categories)
        smaller_distribution = output_weight_distribution(smaller, 2, 3)
        larger_distribution = output_weight_distribution(larger, 2, 3)
        smaller_cdf = Fraction(0)
        larger_cdf = Fraction(0)
        for smaller_mass, larger_mass in zip(
            smaller_distribution, larger_distribution
        ):
            smaller_cdf += smaller_mass
            larger_cdf += larger_mass
            if smaller_cdf < larger_cdf:
                raise AssertionError("input-support monotonicity failed")
        checked += 1
    return {"support_pairs_checked": checked}


def matrix_multiply(left, right):
    return (
        left[0] * right[0] + left[1] * right[2],
        left[0] * right[1] + left[1] * right[3],
        left[2] * right[0] + left[3] * right[2],
        left[2] * right[1] + left[3] * right[3],
    )


def matrix_power(matrix, exponent: int):
    result = (mp.iv.mpf(1), mp.iv.mpf(0), mp.iv.mpf(0), mp.iv.mpf(1))
    power = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = matrix_multiply(result, power)
        remaining >>= 1
        if remaining:
            power = matrix_multiply(power, power)
    return result


def transition_by_rank(packet_bits: int, sigma: int, z, rank: int):
    state_zero = mp.iv.mpf(1) / (1 << sigma)
    output_moment = ((1 + z) / 2) ** packet_bits
    off = state_zero * output_moment
    live = (1 - state_zero) * output_moment
    zero = (mp.iv.mpf(1), mp.iv.mpf(0), off, live)
    if rank == 0:
        return zero
    nonzero = (off, live, off, live)
    packet_zero = mp.iv.mpf(1) / (1 << rank)
    return tuple(
        packet_zero * zero_value + (1 - packet_zero) * nonzero_value
        for zero_value, nonzero_value in zip(zero, nonzero)
    )


def one_marked_row_average(zero, marked, positions: int):
    right_powers = [
        (mp.iv.mpf(1), mp.iv.mpf(0), mp.iv.mpf(0), mp.iv.mpf(1))
    ]
    for _ in range(positions - 1):
        right_powers.append(matrix_multiply(right_powers[-1], zero))
    total = tuple(mp.iv.mpf(0) for _ in range(4))
    left = right_powers[0]
    for position in range(positions):
        term = matrix_multiply(
            matrix_multiply(left, marked), right_powers[positions - 1 - position]
        )
        total = tuple(a + b for a, b in zip(total, term))
        left = matrix_multiply(left, zero)
    return tuple(entry / positions for entry in total)


def endpoint_text(value) -> dict[str, str]:
    return {"lower": str(value.a), "upper": str(value.b)}


def certify_sparse_point(
    *,
    active_blocks: int,
    z_text: str,
    message_bits: int,
    outer_bits: int,
    packet_bits: int,
    sigma: int,
    relative_distance_numerator: int,
    relative_distance_denominator: int,
) -> dict[str, object]:
    outer_dimension = outer_bits // 2
    blocks = message_bits // outer_dimension
    groups = blocks // packet_bits
    output_bits = outer_bits * groups * packet_bits
    distance = (
        relative_distance_numerator * output_bits
    ) // relative_distance_denominator
    if not 1 <= active_blocks < packet_bits:
        raise ValueError("sparse certificate expects one partial packed group")

    z = mp.iv.mpf(z_text)
    zero = transition_by_rank(packet_bits, sigma, z, 0)
    marked = transition_by_rank(packet_bits, sigma, z, active_blocks)
    row = one_marked_row_average(zero, marked, groups)
    total = matrix_power(row, outer_bits)
    moment = total[0] + total[1]
    log2_inner = (
        mp.iv.log(moment) - distance * mp.iv.log(z)
    ) / mp.iv.log(2)

    outer_count = (
        math.comb(blocks, active_blocks)
        * ((1 << outer_dimension) - 1) ** active_blocks
    )
    conditioning = (
        mp.iv.mpf(1) - mp.iv.mpf(1) / (1 << outer_bits)
    ) ** (-active_blocks)
    log2_outer = (
        mp.iv.log(mp.iv.mpf(outer_count)) + mp.iv.log(conditioning)
    ) / mp.iv.log(2)
    pointwise = log2_outer + log2_inner
    return {
        "active_outer_blocks": active_blocks,
        "fixed_z": z_text,
        "outer_log2": endpoint_text(log2_outer),
        "inner_log2_upper": endpoint_text(log2_inner),
        "pointwise_log2_upper": endpoint_text(pointwise),
        "lambda_bits_lower": endpoint_text(-pointwise),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mp.iv.dps = args.dps
    packing_check = packing_self_test()
    common = {
        "message_bits": 1 << 20,
        "outer_bits": 1024,
        "packet_bits": 4,
        "sigma": 15,
        "relative_distance_numerator": 9,
        "relative_distance_denominator": 100,
    }
    points = [
        certify_sparse_point(
            active_blocks=1,
            z_text="0.99917523538219895",
            **common,
        ),
        certify_sparse_point(
            active_blocks=2,
            z_text="0.99833982234889927",
            **common,
        ),
    ]
    payload = {
        "schema": "riffle-transpose-packetshuffle-sparse-interval-v1",
        "interval_decimal_digits": args.dps,
        "packing_small_check": packing_check,
        "support_monotonicity_small_check": support_monotonicity_self_test(),
        "parameters": common,
        "points": points,
        "scope": (
            "Each point is an interval enclosure at its displayed fixed "
            "Chernoff parameter. The one-block row is exact. The two-block "
            "row uses the proved packed-group domination."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
