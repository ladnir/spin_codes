#!/usr/bin/env python3
"""Bound occupation three after removing every all-one BCH endpoint word.

For supports that exclude p1, the first outer check makes the three BCH words
XOR to zero.  Their exact weights therefore satisfy all three Hamming triangle
inequalities.  For a fixed ordered weight triple, projection onto any outer
coordinate is injective on the one-dimensional support code.  The number of
outer words with that profile is at most the smallest BCH-spectrum count in
the triple.

For supports containing p1, the two finite field values are equal.  Their
weight profile is (h,h,k), and its multiplicity is at most min(A_h,A_k).
The script transfers both exact-weight envelopes through the compressed inner
operator.  Weight 128 is excluded because Goal 10 closed every word that
contains the all-one BCH block.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

import analyze_riffle_shiftalpha64_compressed_return as compressed
import analyze_riffle_shiftalpha64_outer_occupation as outer


DATA_BLOCKS = 16_384
FINITE_SUPPORTS = (
    DATA_BLOCKS * (DATA_BLOCKS - 1) * (DATA_BLOCKS - 2) // 6
    + DATA_BLOCKS * (DATA_BLOCKS - 1) // 2
)
P1_SUPPORTS = (DATA_BLOCKS + 1) * DATA_BLOCKS // 2
FIELD_NONZERO = (1 << 64) - 1
GOAL13_LOW_COMPLEMENT = {
    22: (22, 0),
    24: (24, 0),
    26: (26, 0),
    102: (26, 1),
    104: (24, 1),
    106: (22, 1),
}


@dataclass(frozen=True)
class Profile:
    weights: tuple[int, int, int]
    total_weight: int
    log_count_bound: float
    log_slice_normalization: float


def included_spectrum(path: Path) -> dict[int, int]:
    spectrum = outer.load_spectrum(path)
    if spectrum.get(128) != 1:
        raise RuntimeError("all-one BCH spectrum entry failed validation")
    return {
        weight: count
        for weight, count in spectrum.items()
        if weight != 0 and weight != 128 and count > 0
    }


def triangle_holds(first: int, second: int, third: int) -> bool:
    return (
        first <= second + third
        and second <= first + third
        and third <= first + second
    )


def finite_profiles(
    spectrum: dict[int, int],
    *,
    exclude_minimum_triple: bool,
    exclude_adjacent_profiles: bool,
    exclude_goal13_initial: bool,
    exclude_goal13_schedule: bool,
    exclude_goal14_complements: bool,
    exclude_goal14_mixed_schedule: bool,
    exclude_goal15_equal_shell_schedule: bool,
) -> list[Profile]:
    profiles = []
    for first, first_count in spectrum.items():
        first_log_slice = compressed.kernel.log_binom(128, first)
        for second, second_count in spectrum.items():
            second_log_slice = compressed.kernel.log_binom(128, second)
            for third, third_count in spectrum.items():
                if exclude_minimum_triple and (first, second, third) == (22, 22, 22):
                    continue
                if exclude_adjacent_profiles and sorted((first, second, third)) in (
                    [22, 22, 24],
                    [22, 106, 106],
                ):
                    continue
                if exclude_goal13_initial and (
                    (first, second, third) == (106, 106, 106)
                    or sorted((first, second, third)) in (
                        [22, 104, 106],
                        [24, 106, 106],
                    )
                ):
                    continue
                if exclude_goal13_schedule:
                    ordered = (first, second, third)
                    if sorted(ordered) == [22, 22, 26]:
                        continue
                    if all(weight in GOAL13_LOW_COMPLEMENT for weight in ordered):
                        complement_parity = sum(
                            GOAL13_LOW_COMPLEMENT[weight][1] for weight in ordered
                        ) & 1
                        if complement_parity:
                            continue
                if exclude_goal14_complements and sorted(
                    (first, second, third)
                ) in (
                    [24, 104, 104],
                    [22, 104, 104],
                    [24, 104, 106],
                    [22, 102, 106],
                    [26, 106, 106],
                ):
                    continue
                if exclude_goal14_mixed_schedule and sorted(
                    (first, second, third)
                ) == [22, 24, 24]:
                    continue
                if exclude_goal15_equal_shell_schedule and (
                    first,
                    second,
                    third,
                ) == (24, 24, 24):
                    continue
                if not triangle_holds(first, second, third):
                    continue
                profiles.append(
                    Profile(
                        weights=(first, second, third),
                        total_weight=first + second + third,
                        log_count_bound=math.log(
                            min(first_count, second_count, third_count)
                        ),
                        log_slice_normalization=(
                            first_log_slice
                            + second_log_slice
                            + compressed.kernel.log_binom(128, third)
                        ),
                    )
                )
    return profiles


def p1_profiles(
    spectrum: dict[int, int], *, exclude_minimum_triple: bool
) -> list[Profile]:
    profiles = []
    for finite_weight, finite_count in spectrum.items():
        finite_log_slice = compressed.kernel.log_binom(128, finite_weight)
        for p1_weight, p1_count in spectrum.items():
            if exclude_minimum_triple and finite_weight == 22 and p1_weight == 22:
                continue
            profiles.append(
                Profile(
                    weights=(finite_weight, finite_weight, p1_weight),
                    total_weight=2 * finite_weight + p1_weight,
                    log_count_bound=math.log(min(finite_count, p1_count)),
                    log_slice_normalization=(
                        2 * finite_log_slice
                        + compressed.kernel.log_binom(128, p1_weight)
                    ),
                )
            )
    return profiles


def profile_envelope(
    profiles: list[Profile], weight_tilt: float
) -> tuple[float, str, list[dict[str, object]]]:
    log_x = math.log(weight_tilt)
    terms = np.asarray(
        [
            profile.log_count_bound
            - profile.log_slice_normalization
            - profile.total_weight * log_x
            for profile in profiles
        ],
        dtype=np.float64,
    )
    summed = float(logsumexp(terms))
    coefficient_terms = np.asarray(
        [
            -profile.log_slice_normalization
            - profile.total_weight * log_x
            for profile in profiles
        ],
        dtype=np.float64,
    )
    message_cap = math.log(FIELD_NONZERO) + float(np.max(coefficient_terms))
    if message_cap < summed:
        selected = message_cap
        method = "field-line message count times the largest profile coefficient"
    else:
        selected = summed
        method = "sum of exact-profile projection bounds"
    top_indices = np.argsort(terms)[-10:][::-1]
    top = [
        {
            "weights": list(profiles[int(index)].weights),
            "total_binary_weight": profiles[int(index)].total_weight,
            "log2_projection_contribution": float(terms[int(index)])
            / math.log(2.0),
        }
        for index in top_indices
    ]
    return selected, method, top


def evaluate_family(
    *,
    profiles: list[Profile],
    support_count: int,
    distance: int,
    packet_positions: int,
    scaled_cost: float,
    zero_scale: float,
    weight_tilt: float,
) -> dict[str, object]:
    if (
        scaled_cost <= 0.0
        or not 0.0 < zero_scale < packet_positions
        or weight_tilt <= 0.0
    ):
        return {"raw_log_bound": math.inf}

    compressed.configure_kernel(3, 22)
    u = scaled_cost / distance
    theta = zero_scale / packet_positions
    state_factor = np.empty(compressed.STATE_COUNT, dtype=np.float64)
    state_factor[0] = 1.0 / theta
    for state in range(1, compressed.STATE_COUNT):
        exponent = u * state.bit_count()
        state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)
    log_path_counts = compressed.weighted_path_log_counts_with_weight_tilt(
        state_factor, weight_tilt
    )

    log_profile_factor, envelope_method, top_profiles = profile_envelope(
        profiles, weight_tilt
    )
    log_outer_factor = math.log(support_count) + log_profile_factor
    zero_gap_log_factor = -math.log1p(-theta)
    support_terms = []
    for support, log_count in enumerate(log_path_counts):
        if not math.isfinite(log_count):
            continue
        term = (
            log_count
            + log_outer_factor
            + (packet_positions - support + 1) * zero_gap_log_factor
            - compressed.kernel.log_binom(packet_positions, support)
        )
        support_terms.append((support, term))
    raw_log_bound = scaled_cost + float(
        logsumexp([term for _, term in support_terms])
    )
    return {
        "raw_log_bound": raw_log_bound,
        "log2_bound": raw_log_bound / math.log(2.0),
        "scaled_positive_cost": scaled_cost,
        "scaled_zero_return_parameter": zero_scale,
        "binary_weight_coefficient_tilt": weight_tilt,
        "profile_envelope_method": envelope_method,
        "top_packet_support": max(support_terms, key=lambda item: item[1])[0],
        "top_profile_projection_terms": top_profiles,
    }


def optimize_family(
    *,
    profiles: list[Profile],
    support_count: int,
    distance: int,
    packet_positions: int,
    maxiter: int,
) -> dict[str, object]:
    starts = (
        np.log(np.asarray([60.0, 23.0, 1.0])),
        np.log(np.asarray([55.0, 22.0, 0.8])),
    )
    results = []

    def objective(point: np.ndarray) -> float:
        scaled_cost, zero_scale, weight_tilt = np.exp(point)
        return float(
            evaluate_family(
                profiles=profiles,
                support_count=support_count,
                distance=distance,
                packet_positions=packet_positions,
                scaled_cost=float(scaled_cost),
                zero_scale=float(zero_scale),
                weight_tilt=float(weight_tilt),
            )["raw_log_bound"]
        )

    for start in starts:
        results.append(
            minimize(
                objective,
                start,
                method="Nelder-Mead",
                options={"xatol": 2e-5, "fatol": 2e-7, "maxiter": maxiter},
            )
        )
    result = min(results, key=lambda item: float(item.fun))
    scaled_cost, zero_scale, weight_tilt = np.exp(result.x)
    selected = evaluate_family(
        profiles=profiles,
        support_count=support_count,
        distance=distance,
        packet_positions=packet_positions,
        scaled_cost=float(scaled_cost),
        zero_scale=float(zero_scale),
        weight_tilt=float(weight_tilt),
    )
    selected.update(
        {
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_evaluations_selected_start": int(result.nfev),
            "optimizer_total_evaluations": int(sum(item.nfev for item in results)),
        }
    )
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=outer.DEFAULT_SPECTRUM)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument(
        "--packet-positions",
        type=int,
        default=compressed.DEFAULT_PACKET_POSITIONS,
    )
    parser.add_argument("--optimizer-maxiter", type=int, default=100)
    parser.add_argument(
        "--exclude-finite-222",
        action="store_true",
        help="use the exact outer-ratio exclusion of finite (22,22,22)",
    )
    parser.add_argument(
        "--exclude-p1-222",
        action="store_true",
        help="use the Sidon-plus-inner closure of the p1 (22,22,22) profile",
    )
    parser.add_argument(
        "--exclude-finite-adjacent",
        action="store_true",
        help="use the Goal 12 closures of finite (22,22,24) and (22,106,106)",
    )
    parser.add_argument(
        "--exclude-finite-goal13-initial",
        action="store_true",
        help="exclude impossible (106,106,106) and closed (22,104,106)",
    )
    parser.add_argument(
        "--exclude-finite-goal13-schedule",
        action="store_true",
        help=(
            "exclude closed (22,22,26) and every impossible odd-complement "
            "triple over low shells {22,24,26}"
        ),
    )
    parser.add_argument(
        "--exclude-finite-goal14-complements",
        action="store_true",
        help="exclude every closed two-complement profile in the Goal 14 core",
    )
    parser.add_argument(
        "--exclude-finite-goal14-mixed-schedule",
        action="store_true",
        help="exclude the exact shifted-schedule zero for (22,24,24)",
    )
    parser.add_argument(
        "--exclude-finite-goal15-equal-shell-schedule",
        action="store_true",
        help="exclude the exactly closed shifted-schedule profile (24,24,24)",
    )
    args = parser.parse_args()

    spectrum = included_spectrum(args.spectrum)
    finite = finite_profiles(
        spectrum,
        exclude_minimum_triple=args.exclude_finite_222,
        exclude_adjacent_profiles=args.exclude_finite_adjacent,
        exclude_goal13_initial=args.exclude_finite_goal13_initial,
        exclude_goal13_schedule=args.exclude_finite_goal13_schedule,
        exclude_goal14_complements=args.exclude_finite_goal14_complements,
        exclude_goal14_mixed_schedule=args.exclude_finite_goal14_mixed_schedule,
        exclude_goal15_equal_shell_schedule=(
            args.exclude_finite_goal15_equal_shell_schedule
        ),
    )
    p1 = p1_profiles(
        spectrum, exclude_minimum_triple=args.exclude_p1_222
    )
    finite_result = optimize_family(
        profiles=finite,
        support_count=FINITE_SUPPORTS,
        distance=args.distance,
        packet_positions=args.packet_positions,
        maxiter=args.optimizer_maxiter,
    )
    p1_result = optimize_family(
        profiles=p1,
        support_count=P1_SUPPORTS,
        distance=args.distance,
        packet_positions=args.packet_positions,
        maxiter=args.optimizer_maxiter,
    )
    combined_log = float(
        logsumexp([finite_result["raw_log_bound"], p1_result["raw_log_bound"]])
    )
    payload = {
        "schema": "riffle-shiftalpha64-occ3-no-endpoint-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "excluded_family": "every outer word containing the all-one BCH block",
        "parameters": {
            "bad_output_weight_inclusive": args.distance,
            "global_packet_positions": args.packet_positions,
        },
        "finite_support_family": {
            "support_count": FINITE_SUPPORTS,
            "ordered_weight_profiles": len(finite),
            "necessary_weight_condition": "all three Hamming triangle inequalities",
            "exactly_excluded_profiles": (
                ([[22, 22, 22]] if args.exclude_finite_222 else [])
                + (
                    [[24, 24, 24]]
                    if args.exclude_finite_goal15_equal_shell_schedule
                    else []
                )
                + (
                    ["all ordered permutations of (22,24,24)"]
                    if args.exclude_finite_goal14_mixed_schedule
                    else []
                )
                + (
                    [
                        "all ordered permutations of (24,104,104)",
                        "all ordered permutations of (22,104,104)",
                        "all ordered permutations of (24,104,106)",
                        "all ordered permutations of (22,102,106)",
                        "all ordered permutations of (26,106,106)",
                    ]
                    if args.exclude_finite_goal14_complements
                    else []
                )
                + (
                    [
                        "all ordered permutations of (22,22,24)",
                        "all ordered permutations of (22,106,106)",
                    ]
                    if args.exclude_finite_adjacent
                    else []
                )
                + (
                    [
                        "all ordered permutations of (22,22,26)",
                        "every odd-complement triple over low shells {22,24,26}",
                    ]
                    if args.exclude_finite_goal13_schedule
                    else []
                )
                + (
                    [
                        "(106,106,106)",
                        "all ordered permutations of (22,104,106)",
                        "all ordered permutations of (24,106,106)",
                    ]
                    if args.exclude_finite_goal13_initial
                    else []
                )
            ),
            "bound": finite_result,
        },
        "second_parity_support_family": {
            "support_count": P1_SUPPORTS,
            "ordered_weight_profiles": len(p1),
            "necessary_weight_condition": "profile has the form (h,h,k)",
            "exactly_excluded_profiles": (
                [[22, 22, 22]] if args.exclude_p1_222 else []
            ),
            "bound": p1_result,
        },
        "combined_log2_bound": combined_log / math.log(2.0),
        "validation": {
            "exact_bch_spectrum": "PASS",
            "weight_128_excluded": True,
            "finite_support_first_check_implies_bch_xor_zero": True,
            "finite_weight22_triple_outer_ratio_scan_used": args.exclude_finite_222,
            "p1_weight22_triple_sidon_closure_used": args.exclude_p1_222,
            "finite_adjacent_profile_closures_used": args.exclude_finite_adjacent,
            "finite_goal13_initial_closures_used": args.exclude_finite_goal13_initial,
            "finite_goal13_schedule_closures_used": args.exclude_finite_goal13_schedule,
            "finite_goal14_complement_closures_used": args.exclude_finite_goal14_complements,
            "finite_goal14_mixed_schedule_zero_used": args.exclude_finite_goal14_mixed_schedule,
            "finite_goal15_equal_shell_schedule_closure_used": (
                args.exclude_finite_goal15_equal_shell_schedule
            ),
            "projection_profile_count_bound": "min(A_w1,A_w2,A_w3)",
            "optimizer_warning_does_not_invalidate_selected_parameters": True,
        },
        "scope": (
            "Rigorous bound for all occupation-three outer words that do not "
            "contain the all-one BCH block. The weight conditions are necessary, "
            "not sufficient for the second outer equation."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
