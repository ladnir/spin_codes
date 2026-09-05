#!/usr/bin/env python3
"""Exact shared-packet-group profiles for the LDPCSplitState recurrence.

Fix a profile (n1,n2,n3,n4), where nr is the number of four-block packet
groups containing r active outer blocks.  The regular outer envelope replaces
each active block's coordinate by an independent Bernoulli(p) bit.  A packet
group assigned to one epoch therefore contributes Binomial(r,p) input bits.

For one transposed region, the script extracts the coefficient of the profile
from 32 ordered epochs of 64 packet slots.  This averages the uniform packet
permutation exactly.  Raising the resulting two-state matrix to the 256th
power gives the complete conditional inner moment for that fixed profile.

The calculation uses the existing worst-case activation bound as a function
of total epoch weight.  It is therefore proof-safe but does not exploit the
stronger packet-support-four property for packed inputs.  Floating-point
arithmetic and the modeled outer spectrum make the output diagnostic rather
than a final certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution
from scipy.signal import convolve

from analyze_riffle_ldpcsplitstate_occupation_ladder import (
    DEFAULT_ACTIVATION,
    epoch_transfers,
    matrix_power_moment,
    regular_outer_density_log2,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/shared_group_low_occupation_profiles.json"
)


def parse_int_list(text: str) -> list[int]:
    return [int(item) for item in text.split(",") if item]


def parse_profiles(text: str) -> dict[int, list[tuple[int, int, int, int]]]:
    """Parse semicolon-separated n1:n2:n3:n4 profiles by occupation."""
    result: dict[int, list[tuple[int, int, int, int]]] = {}
    if not text:
        return result
    for item in text.split(";"):
        profile = tuple(int(value) for value in item.split(":"))
        if len(profile) != 4 or any(value < 0 for value in profile):
            raise ValueError(f"invalid profile: {item}")
        occupation = sum((index + 1) * value for index, value in enumerate(profile))
        result.setdefault(occupation, []).append(profile)
    return result


def profiles_of_occupation(occupation: int):
    """Yield all (n1,n2,n3,n4) with sum_r r*nr=occupation."""
    for n4 in range(occupation // 4 + 1):
        after4 = occupation - 4 * n4
        for n3 in range(after4 // 3 + 1):
            after3 = after4 - 3 * n3
            for n2 in range(after3 // 2 + 1):
                n1 = after3 - 2 * n2
                yield (n1, n2, n3, n4)


def profile_block_set_log2(profile: tuple[int, int, int, int]) -> float:
    """Return log2 of active block sets with the fixed group profile."""
    groups = sum(profile)
    value = math.lgamma(2049) - math.lgamma(2049 - groups)
    for count in profile:
        value -= math.lgamma(count + 1)
    for width, count in enumerate(profile, start=1):
        value += count * math.log(math.comb(4, width))
    return value / math.log(2.0)


def averaged_epoch_transfers(
    epoch: list[np.ndarray], probability: float, maximum_mass: int
) -> list[np.ndarray]:
    """Average one epoch transfer over independent candidate input bits."""
    averaged = []
    for mass in range(maximum_mass + 1):
        matrix = np.zeros((2, 2), dtype=np.float64)
        for weight in range(mass + 1):
            coefficient = (
                math.comb(mass, weight)
                * probability**weight
                * (1.0 - probability) ** (mass - weight)
            )
            matrix += coefficient * epoch[weight]
        averaged.append(matrix)
    return averaged


def slot_type_count(counts: tuple[int, int, int, int]) -> int:
    """Count placements of typed packets into one 64-slot epoch."""
    total = sum(counts)
    result = math.comb(64, total) * math.factorial(total)
    for count in counts:
        result //= math.factorial(count)
    return result


def global_type_count(profile: tuple[int, int, int, int]) -> int:
    """Count placements of typed packets into all 2,048 packet slots."""
    total = sum(profile)
    result = math.comb(2048, total) * math.factorial(total)
    for count in profile:
        result //= math.factorial(count)
    return result


def truncated_matrix_polynomial_product(
    left: np.ndarray, right: np.ndarray, coefficient_shape: tuple[int, ...]
) -> np.ndarray:
    """Multiply two four-variable 2x2 matrix polynomials and truncate."""
    result = np.zeros((*coefficient_shape, 2, 2), dtype=np.float64)
    slices = tuple(slice(0, size) for size in coefficient_shape)
    for row in range(2):
        for column in range(2):
            for inner in range(2):
                convolution = convolve(
                    left[..., row, inner],
                    right[..., inner, column],
                    method="direct",
                )
                result[..., row, column] += convolution[slices]
    # FFT roundoff can create tiny negative coefficients.
    np.maximum(result, 0.0, out=result)
    return result


def fast_profile_region_transfer(
    averaged_epoch: list[np.ndarray], profile: tuple[int, int, int, int]
) -> np.ndarray:
    """Extract the profile coefficient from the 32nd matrix-polynomial power."""
    coefficient_shape = tuple(count + 1 for count in profile)
    polynomial = np.zeros((*coefficient_shape, 2, 2), dtype=np.float64)
    for added in itertools.product(*(range(count + 1) for count in profile)):
        if sum(added) > 64:
            continue
        mass = sum((index + 1) * count for index, count in enumerate(added))
        polynomial[added] = averaged_epoch[mass] * slot_type_count(added)

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

    normalization_log = math.log(global_type_count(profile))
    return polynomial[profile] * math.exp(log_scale - normalization_log)


def profile_region_transfer(
    averaged_epoch: list[np.ndarray], profile: tuple[int, int, int, int]
) -> np.ndarray:
    """Return the exact one-region transfer for a fixed group profile."""
    nonzero_types = [index for index, count in enumerate(profile) if count]
    if len(nonzero_types) == 1:
        packet_type = nonzero_types[0]
        packet_width = packet_type + 1
        packet_count = profile[packet_type]
        coefficients = [
            np.zeros((2, 2), dtype=np.float64)
            for _ in range(packet_count + 1)
        ]
        coefficients[0] = np.eye(2)
        for _ in range(32):
            following = [
                np.zeros((2, 2), dtype=np.float64)
                for _ in range(packet_count + 1)
            ]
            for used, prefix_matrix in enumerate(coefficients):
                for added in range(packet_count - used + 1):
                    following[used + added] += (
                        prefix_matrix
                        @ averaged_epoch[packet_width * added]
                        * math.comb(64, added)
                    )
            coefficients = following
        return coefficients[packet_count] / math.comb(2048, packet_count)

    coefficient_count = math.prod(count + 1 for count in profile)
    if coefficient_count >= 512:
        return fast_profile_region_transfer(averaged_epoch, profile)

    zero = (0, 0, 0, 0)
    coefficients: dict[tuple[int, int, int, int], np.ndarray] = {
        zero: np.eye(2)
    }
    additions = list(
        itertools.product(*(range(count + 1) for count in profile))
    )
    additions = [counts for counts in additions if sum(counts) <= 64]
    factors = {counts: slot_type_count(counts) for counts in additions}
    masses = {
        counts: sum((width + 1) * count for width, count in enumerate(counts))
        for counts in additions
    }

    for _ in range(32):
        following: dict[tuple[int, int, int, int], np.ndarray] = {}
        for used, prefix_matrix in coefficients.items():
            for added in additions:
                target = tuple(used[i] + added[i] for i in range(4))
                if any(target[i] > profile[i] for i in range(4)):
                    continue
                contribution = (
                    prefix_matrix
                    @ averaged_epoch[masses[added]]
                    * factors[added]
                )
                if target in following:
                    following[target] += contribution
                else:
                    following[target] = contribution
        coefficients = following
    return coefficients[profile] / global_type_count(profile)


def reference_points(receipt_paths: list[Path]) -> dict[int, tuple[float, float]]:
    """Read the distinct-group optimizer's (log surprisal, p) points."""
    result: dict[int, tuple[float, float]] = {}
    for path in receipt_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        for occupation, row in payload["best"].items():
            result[int(occupation)] = (
                float(row["log_surprisal"]),
                float(row["placement_probability"]),
            )
    return result


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
) -> dict[str, object]:
    occupation = sum((index + 1) * count for index, count in enumerate(profile))
    surprisal = math.exp(log_surprisal)
    z = math.exp(-surprisal)
    probability = 1.0 / (1.0 + math.exp(-placement_logit))
    epoch = epoch_transfers(
        z=z,
        distance=constituent_distance,
        moment_order=live_moment_order,
        activation_upper=activation_upper,
        maximum=len(activation_upper) - 1,
    )
    if termination_scale != 1.0:
        for matrix in epoch:
            matrix[1, 0] *= termination_scale
    averaged_epoch = averaged_epoch_transfers(epoch, probability, occupation)
    region = profile_region_transfer(averaged_epoch, profile)
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
    maximum = max(occupations)
    activation_upper = [0.0] * (maximum + 1)
    activation_upper[0] = 1.0
    for row in activation["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight > maximum:
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
        conditioning_bits = occupation * regular_outer_density_log2(probability)

        profile_rows = []
        profiles = selected_profiles.get(occupation)
        if profiles is None:
            profiles = list(profiles_of_occupation(occupation))
        if args.minimum_n4:
            profiles = [
                profile for profile in profiles if profile[3] >= args.minimum_n4
            ]
        distinct_profile = (occupation, 0, 0, 0)
        if distinct_profile not in profiles:
            profiles = [distinct_profile, *profiles]
        for profile in profiles:
            reference_row = evaluate_profile_point(
                profile=profile,
                log_surprisal=log_surprisal,
                placement_logit=placement_logit,
                activation_upper=activation_upper,
                constituent_distance=args.constituent_distance,
                live_moment_order=args.live_moment_order,
                target_distance=target_distance,
            )
            if args.compare_no_termination:
                no_termination_row = evaluate_profile_point(
                    profile=profile,
                    log_surprisal=log_surprisal,
                    placement_logit=placement_logit,
                    activation_upper=activation_upper,
                    constituent_distance=args.constituent_distance,
                    live_moment_order=args.live_moment_order,
                    target_distance=target_distance,
                    termination_scale=0.0,
                )
                reference_row["no_termination_log2_inner_moment"] = (
                    no_termination_row["log2_inner_moment"]
                )
                reference_row["termination_log2_inflation"] = (
                    float(reference_row["log2_inner_moment"])
                    - float(no_termination_row["log2_inner_moment"])
                )
            if args.optimize_selected and profile != distinct_profile:
                def objective(point: np.ndarray) -> float:
                    return float(
                        evaluate_profile_point(
                            profile=profile,
                            log_surprisal=float(point[0]),
                            placement_logit=float(point[1]),
                            activation_upper=activation_upper,
                            constituent_distance=args.constituent_distance,
                            live_moment_order=args.live_moment_order,
                            target_distance=target_distance,
                        )["aggregate_log2_bound"]
                    )

                optimum = differential_evolution(
                    objective,
                    bounds=((-10.0, 1.0), (-8.0, 5.0)),
                    seed=args.optimizer_seed,
                    popsize=args.optimizer_popsize,
                    maxiter=args.optimizer_iterations,
                    polish=True,
                    workers=1,
                    updating="immediate",
                )
                optimized_row = evaluate_profile_point(
                    profile=profile,
                    log_surprisal=float(optimum.x[0]),
                    placement_logit=float(optimum.x[1]),
                    activation_upper=activation_upper,
                    constituent_distance=args.constituent_distance,
                    live_moment_order=args.live_moment_order,
                    target_distance=target_distance,
                )
                optimized_row["reference_point_margin_bits"] = reference_row[
                    "aggregate_margin_bits"
                ]
                optimized_row["optimizer_success"] = bool(optimum.success)
                optimized_row["optimizer_message"] = str(optimum.message)
                profile_rows.append(optimized_row)
            else:
                profile_rows.append(reference_row)
        worst = min(profile_rows, key=lambda row: row["aggregate_margin_bits"])
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
                "placement_conditioning_bits": conditioning_bits,
                "profile_count": len(profile_rows),
                "distinct_group_row": distinct,
                "worst_profile_row": worst,
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
        "schema": "riffle-ldpcsplitstate-shared-group-profiles-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "occupations": occupations,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
            "relative_distance": args.relative_distance,
            "target_distance": target_distance,
            "reference_receipts": [str(path) for path in args.reference_receipts],
        },
        "occupation_rows": occupation_rows,
        "scope": (
            "Exact coefficient transfer for every four-block group profile at "
            "the distinct-group optimizer point. The regular outer spectrum "
            "uses a floating-point pointwise density envelope. Epoch "
            "activation uses the worst total-weight bound and does not exploit "
            "the stronger packet-support audit. Results are diagnostics, not "
            "outward-rounded certificates."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--occupations", default="2,3,4,5,6,7,8")
    parser.add_argument(
        "--selected-profiles",
        default="",
        help="semicolon-separated n1:n2:n3:n4 profiles; overrides occupations",
    )
    parser.add_argument("--optimize-selected", action="store_true")
    parser.add_argument("--minimum-n4", type=int, default=0)
    parser.add_argument("--compare-no-termination", action="store_true")
    parser.add_argument("--fixed-log-surprisal", type=float)
    parser.add_argument("--fixed-placement-probability", type=float)
    parser.add_argument("--optimizer-seed", type=int, default=20260829)
    parser.add_argument("--optimizer-popsize", type=int, default=6)
    parser.add_argument("--optimizer-iterations", type=int, default=12)
    parser.add_argument("--constituent-distance", type=int, default=40)
    parser.add_argument(
        "--live-moment-order", type=int, choices=(1, 2, 3), default=3
    )
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--reference-receipts",
        type=Path,
        nargs="+",
        default=[
            Path(
                "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
                "receipts/occupation_ladder_regular_1_16.json"
            ),
            Path(
                "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
                "receipts/occupation_ladder_regular_32_64.json"
            ),
            Path(
                "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
                "receipts/occupation_ladder_regular_80_112.json"
            ),
            Path(
                "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
                "receipts/occupation_ladder_regular_128_moment3.json"
            ),
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
