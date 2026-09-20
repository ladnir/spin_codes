#!/usr/bin/env python3
"""Exact shared-group profiles after independent GF(16) packet multipliers.

Setup multiplies every four-bit packet by an independent uniform nonzero
element of GF(16).  For a packet group containing r Bernoulli(p) candidate
bits, the packet is zero with probability (1-p)^r.  Conditioned on being
nonzero, its multiplied value is uniform over the 15 nonzero nibbles.

The epoch transfer therefore depends on the complete tuple of packet widths
placed in that epoch, rather than only their total candidate mass.  This
script retains that tuple in the exact 32-epoch coefficient calculation.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

from analyze_riffle_ldpcsplitstate_occupation_ladder import (
    epoch_transfers,
    matrix_power_moment,
    regular_outer_density_log2,
)
from analyze_riffle_ldpcsplitstate_shared_groups import (
    DEFAULT_ACTIVATION,
    global_type_count,
    parse_int_list,
    parse_profiles,
    profile_block_set_log2,
    profiles_of_occupation,
    reference_points,
    slot_type_count,
    truncated_matrix_polynomial_product,
)


CONSTRUCTION = Path(
    "constructions/riffle_ldpcsplitstate_packetmul_g4_t256_s64"
)
DEFAULT_OUTPUT = CONSTRUCTION / "receipts/shared_group_low_occupation_profiles.json"


def packet_weight_distribution(width: int, probability: float) -> np.ndarray:
    """Return the post-multiplier binary-weight law for one packet."""
    zero_probability = (1.0 - probability) ** width
    nonzero_scale = (1.0 - zero_probability) / 15.0
    result = np.zeros(5, dtype=np.float64)
    result[0] = zero_probability
    for weight in range(1, 5):
        result[weight] = nonzero_scale * math.comb(4, weight)
    return result


def packet_scalar_moment(width: int, probability: float, u: float) -> float:
    """Evaluate the post-multiplier packet weight generating function."""
    distribution = packet_weight_distribution(width, probability)
    return float(sum(value * u**weight for weight, value in enumerate(distribution)))


def packetmul_epoch_transfer_factory(
    epoch: list[np.ndarray], probability: float
):
    """Return a cached transfer indexed by counts of packet widths 1..4."""
    packet_laws = [
        packet_weight_distribution(width, probability) for width in range(1, 5)
    ]

    @lru_cache(maxsize=None)
    def transfer(counts: tuple[int, int, int, int]) -> np.ndarray:
        distribution = np.asarray((1.0,), dtype=np.float64)
        for packet_type, count in enumerate(counts):
            for _ in range(count):
                distribution = np.convolve(distribution, packet_laws[packet_type])
        matrix = np.zeros((2, 2), dtype=np.float64)
        for weight, coefficient in enumerate(distribution):
            if coefficient:
                matrix += coefficient * epoch[weight]
        return matrix

    return transfer


def universal_nonzero_region_transfers(
    epoch: list[np.ndarray], maximum_packets: int
) -> list[np.ndarray]:
    """Return region matrices conditioned on k uniform nonzero packets."""
    nonzero_law = np.asarray(
        (0.0, 4.0 / 15.0, 6.0 / 15.0, 4.0 / 15.0, 1.0 / 15.0),
        dtype=np.float64,
    )
    epoch_by_packets = []
    distribution = np.asarray((1.0,), dtype=np.float64)
    for packet_count in range(min(64, maximum_packets) + 1):
        matrix = np.zeros((2, 2), dtype=np.float64)
        for weight, coefficient in enumerate(distribution):
            if coefficient:
                matrix += coefficient * epoch[weight]
        epoch_by_packets.append(matrix)
        distribution = np.convolve(distribution, nonzero_law)

    coefficients = [np.zeros((2, 2)) for _ in range(maximum_packets + 1)]
    coefficients[0] = np.eye(2)
    for _ in range(32):
        following = [np.zeros((2, 2)) for _ in range(maximum_packets + 1)]
        for used, prefix in enumerate(coefficients):
            for added in range(min(64, maximum_packets - used) + 1):
                following[used + added] += (
                    prefix @ epoch_by_packets[added] * math.comb(64, added)
                )
        coefficients = following
    return [
        coefficients[count] / math.comb(2048, count)
        for count in range(maximum_packets + 1)
    ]


def survivor_count_distribution(
    profile: tuple[int, int, int, int], probability: float
) -> np.ndarray:
    """Return the number of nonzero multiplied packets for one profile."""
    result = np.asarray((1.0,), dtype=np.float64)
    for packet_type, count in enumerate(profile):
        width = packet_type + 1
        survival = 1.0 - (1.0 - probability) ** width
        typed = np.asarray(
            [
                math.comb(count, survivors)
                * survival**survivors
                * (1.0 - survival) ** (count - survivors)
                for survivors in range(count + 1)
            ],
            dtype=np.float64,
        )
        result = np.convolve(result, typed)
    return result


def factored_profile_region_transfer(
    universal: list[np.ndarray],
    profile: tuple[int, int, int, int],
    probability: float,
) -> np.ndarray:
    """Mix universal matrices by the profile's scalar survivor law."""
    survivor_law = survivor_count_distribution(profile, probability)
    region = np.zeros((2, 2), dtype=np.float64)
    for count, coefficient in enumerate(survivor_law):
        region += coefficient * universal[count]
    return region


def packetmul_region_factory(
    epoch: list[np.ndarray], probability: float, maximum_groups: int
):
    """Build the universal survivor matrices once for a common proof point."""
    universal = universal_nonzero_region_transfers(epoch, maximum_groups)

    def region(profile: tuple[int, int, int, int]) -> np.ndarray:
        return factored_profile_region_transfer(universal, profile, probability)

    return region


def fast_profile_region_transfer(
    transfer, profile: tuple[int, int, int, int]
) -> np.ndarray:
    """Extract one typed profile coefficient from the 32nd matrix power."""
    coefficient_shape = tuple(count + 1 for count in profile)
    polynomial = np.zeros((*coefficient_shape, 2, 2), dtype=np.float64)
    for added in itertools.product(*(range(count + 1) for count in profile)):
        if sum(added) > 64:
            continue
        polynomial[added] = transfer(added) * slot_type_count(added)

    maximum = float(np.max(polynomial))
    polynomial /= maximum
    log_scale = math.log(maximum)
    for _ in range(5):
        polynomial = truncated_matrix_polynomial_product(
            polynomial, polynomial, coefficient_shape
        )
        log_scale *= 2.0
        maximum = float(np.max(polynomial))
        polynomial /= maximum
        log_scale += math.log(maximum)

    return polynomial[profile] * math.exp(
        log_scale - math.log(global_type_count(profile))
    )


def profile_region_transfer(
    transfer, profile: tuple[int, int, int, int]
) -> np.ndarray:
    """Return the exact PacketMul one-region transfer for one profile."""
    nonzero_types = [index for index, count in enumerate(profile) if count]
    if len(nonzero_types) == 1:
        packet_type = nonzero_types[0]
        packet_count = profile[packet_type]
        coefficients = [np.zeros((2, 2)) for _ in range(packet_count + 1)]
        coefficients[0] = np.eye(2)
        for _ in range(32):
            following = [np.zeros((2, 2)) for _ in range(packet_count + 1)]
            for used, prefix in enumerate(coefficients):
                for added_count in range(min(64, packet_count - used) + 1):
                    added = [0, 0, 0, 0]
                    added[packet_type] = added_count
                    following[used + added_count] += (
                        prefix
                        @ transfer(tuple(added))
                        * math.comb(64, added_count)
                    )
            coefficients = following
        return coefficients[packet_count] / math.comb(2048, packet_count)

    coefficient_count = math.prod(count + 1 for count in profile)
    if coefficient_count >= 512:
        return fast_profile_region_transfer(transfer, profile)

    zero = (0, 0, 0, 0)
    additions = list(itertools.product(*(range(count + 1) for count in profile)))
    additions = [counts for counts in additions if sum(counts) <= 64]
    factors = {counts: slot_type_count(counts) for counts in additions}
    coefficients: dict[tuple[int, int, int, int], np.ndarray] = {
        zero: np.eye(2)
    }
    for _ in range(32):
        following: dict[tuple[int, int, int, int], np.ndarray] = {}
        for used, prefix in coefficients.items():
            for added in additions:
                target = tuple(used[index] + added[index] for index in range(4))
                if any(target[index] > profile[index] for index in range(4)):
                    continue
                contribution = prefix @ transfer(added) * factors[added]
                if target in following:
                    following[target] += contribution
                else:
                    following[target] = contribution
        coefficients = following
    return coefficients[profile] / global_type_count(profile)


def evaluate_profile_point(
    *,
    profile: tuple[int, int, int, int],
    log_surprisal: float,
    placement_logit: float,
    activation_upper: list[float],
    constituent_distance: int,
    live_moment_order: int,
    target_distance: int,
    termination_scale: float = 1.0,
    region_matrix: np.ndarray | None = None,
) -> dict[str, object]:
    occupation = sum((index + 1) * count for index, count in enumerate(profile))
    surprisal = math.exp(log_surprisal)
    z = math.exp(-surprisal)
    probability = 1.0 / (1.0 + math.exp(-placement_logit))
    if region_matrix is None:
        epoch = epoch_transfers(
            z=z,
            distance=constituent_distance,
            moment_order=live_moment_order,
            activation_upper=activation_upper,
            maximum=256,
        )
        if termination_scale != 1.0:
            for matrix in epoch:
                matrix[1, 0] *= termination_scale
        region_function = packetmul_region_factory(
            epoch, probability, sum(profile)
        )
        region = region_function(profile)
    else:
        region = region_matrix
    inner_log2 = matrix_power_moment(region, 256) / math.log(2.0)
    conditioning_bits = occupation * regular_outer_density_log2(probability)
    multiplicity_bits = profile_block_set_log2(profile)
    log2_bound = (
        inner_log2
        + conditioning_bits
        + multiplicity_bits
        + target_distance * surprisal / math.log(2.0)
    )
    return {
        "profile_n1_n2_n3_n4": list(profile),
        "active_packet_groups": sum(profile),
        "log_surprisal": log_surprisal,
        "z": z,
        "placement_logit": placement_logit,
        "placement_probability": probability,
        "placement_conditioning_bits": conditioning_bits,
        "block_set_multiplicity_bits": multiplicity_bits,
        "log2_inner_moment": inner_log2,
        "aggregate_log2_bound": log2_bound,
        "aggregate_margin_bits": -log2_bound,
        "region_matrix": region.tolist(),
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    occupations = parse_int_list(args.occupations)
    selected_profiles = parse_profiles(args.selected_profiles)
    if selected_profiles:
        occupations = sorted(selected_profiles)
    references = reference_points(args.reference_receipts)
    activation = json.loads(args.activation.read_text(encoding="utf-8"))
    maximum_weight = 256
    activation_upper = [0.0] * (maximum_weight + 1)
    activation_upper[0] = 1.0
    for row in activation["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight > maximum_weight:
            break
        activation_upper[weight] = float(
            row["maximum_distinct_conditioned_upper_bound"]
        )

    target_distance = math.floor(args.relative_distance * (1 << 21))
    occupation_rows = []
    for occupation in occupations:
        if occupation not in references:
            raise ValueError(f"no reference optimizer point for {occupation}")
        log_surprisal, probability = references[occupation]
        if args.fixed_log_surprisal is not None:
            log_surprisal = args.fixed_log_surprisal
        if args.fixed_placement_probability is not None:
            probability = args.fixed_placement_probability
        placement_logit = math.log(probability / (1.0 - probability))

        profiles = selected_profiles.get(occupation)
        if profiles is None:
            profiles = list(profiles_of_occupation(occupation))
        if args.minimum_n4:
            profiles = [profile for profile in profiles if profile[3] >= args.minimum_n4]
        distinct_profile = (occupation, 0, 0, 0)
        if distinct_profile not in profiles:
            profiles = [distinct_profile, *profiles]

        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        epoch = epoch_transfers(
            z=z,
            distance=args.constituent_distance,
            moment_order=args.live_moment_order,
            activation_upper=activation_upper,
            maximum=256,
        )
        maximum_groups = max(sum(profile) for profile in profiles)
        region_function = packetmul_region_factory(
            epoch, probability, maximum_groups
        )
        no_termination_region_function = None
        if args.compare_no_termination:
            no_termination_epoch = [matrix.copy() for matrix in epoch]
            for matrix in no_termination_epoch:
                matrix[1, 0] = 0.0
            no_termination_region_function = packetmul_region_factory(
                no_termination_epoch, probability, maximum_groups
            )

        profile_rows = []
        for profile in profiles:
            row = evaluate_profile_point(
                profile=profile,
                log_surprisal=log_surprisal,
                placement_logit=placement_logit,
                activation_upper=activation_upper,
                constituent_distance=args.constituent_distance,
                live_moment_order=args.live_moment_order,
                target_distance=target_distance,
                region_matrix=region_function(profile),
            )
            if args.compare_no_termination:
                no_termination = evaluate_profile_point(
                    profile=profile,
                    log_surprisal=log_surprisal,
                    placement_logit=placement_logit,
                    activation_upper=activation_upper,
                    constituent_distance=args.constituent_distance,
                    live_moment_order=args.live_moment_order,
                    target_distance=target_distance,
                    termination_scale=0.0,
                    region_matrix=no_termination_region_function(profile),
                )
                row["no_termination_log2_inner_moment"] = no_termination[
                    "log2_inner_moment"
                ]
                row["termination_log2_inflation"] = (
                    row["log2_inner_moment"]
                    - no_termination["log2_inner_moment"]
                )
            profile_rows.append(row)

        worst = min(profile_rows, key=lambda row: row["aggregate_margin_bits"])
        minimum_margin = float(worst["aggregate_margin_bits"])
        profile_sum_margin = minimum_margin - math.log2(
            sum(
                2.0
                ** -(
                    float(row["aggregate_margin_bits"])
                    - minimum_margin
                )
                for row in profile_rows
            )
        )
        distinct = next(
            row
            for row in profile_rows
            if row["profile_n1_n2_n3_n4"] == [occupation, 0, 0, 0]
        )
        occupation_rows.append(
            {
                "occupation": occupation,
                "log_surprisal": log_surprisal,
                "placement_probability": probability,
                "profile_count": len(profile_rows),
                "distinct_group_row": distinct,
                "worst_profile_row": worst,
                "profile_sum_margin_bits": profile_sum_margin,
                "worst_minus_distinct_margin_bits": (
                    worst["aggregate_margin_bits"]
                    - distinct["aggregate_margin_bits"]
                ),
                "profiles": profile_rows,
            }
        )
        print(
            f"occupation,{occupation},profiles,{len(profile_rows)},"
            f"distinct_margin,{distinct['aggregate_margin_bits']:.6f},"
            f"worst_profile,{worst['profile_n1_n2_n3_n4']},"
            f"worst_margin,{worst['aggregate_margin_bits']:.6f}",
            flush=True,
        )

    return {
        "schema": "riffle-ldpcsplitstate-packetmul-shared-group-profiles-v1",
        "candidate": "Riffle LDPCSplitState PacketMul g=4 t=256 s=64",
        "parameters": {
            "occupations": occupations,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
            "relative_distance": args.relative_distance,
            "target_distance": target_distance,
            "packet_multiplier": "independent uniform GF(16)^* multiplier per packet",
            "reference_receipts": [str(path) for path in args.reference_receipts],
        },
        "occupation_rows": occupation_rows,
        "scope": (
            "Exact averaging over GF(16) packet multipliers, Bernoulli outer "
            "bits, typed packet placement, and the two-state epoch bound. "
            "The outer spectrum and matrix arithmetic remain floating-point "
            "diagnostics rather than an outward-rounded certificate."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--occupations", default="2,4,8,16")
    parser.add_argument("--selected-profiles", default="")
    parser.add_argument("--minimum-n4", type=int, default=0)
    parser.add_argument("--compare-no-termination", action="store_true")
    parser.add_argument("--fixed-log-surprisal", type=float)
    parser.add_argument("--fixed-placement-probability", type=float)
    parser.add_argument("--constituent-distance", type=int, default=40)
    parser.add_argument("--live-moment-order", type=int, choices=(1, 2, 3), default=3)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--reference-receipts",
        type=Path,
        nargs="+",
        default=[
            Path("constructions/riffle_ldpcsplitstate_g4_t256_s64/receipts/occupation_ladder_regular_1_16.json"),
            Path("constructions/riffle_ldpcsplitstate_g4_t256_s64/receipts/occupation_ladder_regular_32_64.json"),
            Path("constructions/riffle_ldpcsplitstate_g4_t256_s64/receipts/occupation_ladder_regular_80_112.json"),
            Path("constructions/riffle_ldpcsplitstate_g4_t256_s64/receipts/occupation_ladder_regular_128_moment3.json"),
        ],
    )
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
