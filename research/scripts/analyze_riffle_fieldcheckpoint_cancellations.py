#!/usr/bin/env python3
"""Exact activation-cancellation diagnostics for FieldCheckpointAccumulate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal03_cancellation_diagnostic.json"
)


def even_group_counts(groups: int, group_bits: int, maximum_weight: int) -> list[int]:
    """Count supports whose intersection with every group has even size."""
    coefficients = [0] * (maximum_weight + 1)
    coefficients[0] = 1
    local = [
        (weight, math.comb(group_bits, weight))
        for weight in range(0, min(group_bits, maximum_weight) + 1, 2)
    ]
    for _ in range(groups):
        updated = [0] * (maximum_weight + 1)
        for current_weight, current_count in enumerate(coefficients):
            if current_count == 0:
                continue
            for local_weight, local_count in local:
                total_weight = current_weight + local_weight
                if total_weight > maximum_weight:
                    break
                updated[total_weight] += current_count * local_count
        coefficients = updated
    return coefficients


def hypergeometric_intersections(
    outer_bits: int, first_weight: int, second_weight: int
) -> list[tuple[int, float]]:
    denominator = math.comb(outer_bits, second_weight)
    lower = max(0, first_weight + second_weight - outer_bits)
    upper = min(first_weight, second_weight)
    return [
        (
            intersection,
            math.comb(first_weight, intersection)
            * math.comb(outer_bits - first_weight, second_weight - intersection)
            / denominator,
        )
        for intersection in range(lower, upper + 1)
    ]


def falling_ratio(numerator_size: int, denominator_size: int, count: int) -> float:
    if numerator_size < count:
        return 0.0
    ratio = 1.0
    for offset in range(count):
        ratio *= (numerator_size - offset) / (denominator_size - offset)
    return ratio


def prefix_nonactivation_probability(
    outer_bits: int,
    first_weight: int,
    second_weight: int,
    prefix_regions: int,
    pair_region_failure: float,
    intersections: list[tuple[int, float]],
) -> float:
    """Average the event that the state never activates in a region prefix."""
    denominator = math.comb(outer_bits, prefix_regions)
    result = 0.0
    for intersection, intersection_probability in intersections:
        singleton_regions = first_weight + second_weight - 2 * intersection
        empty_regions = outer_bits - singleton_regions - intersection
        conditional = 0.0
        lower = max(0, prefix_regions - empty_regions)
        upper = min(intersection, prefix_regions)
        for failed_overlaps in range(lower, upper + 1):
            conditional += (
                math.comb(intersection, failed_overlaps)
                * math.comb(empty_regions, prefix_regions - failed_overlaps)
                / denominator
                * pair_region_failure**failed_overlaps
            )
        result += intersection_probability * conditional
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    if args.region_bits != args.epochs_per_region * args.epoch_bits:
        raise ValueError("epoch geometry does not fill one region")
    if args.epoch_bits % args.state_bits:
        raise ValueError("state width must divide the epoch")
    visits = args.epoch_bits // args.state_bits
    epoch_counts = even_group_counts(
        args.state_bits, visits, args.maximum_local_weight
    )
    region_groups = args.state_bits * args.epochs_per_region
    region_counts = even_group_counts(
        region_groups, visits, args.maximum_local_weight
    )
    if epoch_counts[2] != args.state_bits * math.comb(visits, 2):
        raise AssertionError("weight-two epoch count is inconsistent")
    if region_counts[2] != region_groups * math.comb(visits, 2):
        raise AssertionError("weight-two region count is inconsistent")

    epoch_rows = []
    region_rows = []
    for weight in range(1, args.maximum_local_weight + 1):
        epoch_probability = epoch_counts[weight] / math.comb(args.epoch_bits, weight)
        region_probability = region_counts[weight] / math.comb(args.region_bits, weight)
        epoch_rows.append(
            {
                "weight": weight,
                "nonactivation_probability": epoch_probability,
                "nonactivation_bits": (
                    -math.log2(epoch_probability) if epoch_probability else math.inf
                ),
            }
        )
        region_rows.append(
            {
                "weight": weight,
                "nonactivation_probability": region_probability,
                "nonactivation_bits": (
                    -math.log2(region_probability) if region_probability else math.inf
                ),
            }
        )

    pair_region_failure = region_rows[1]["nonactivation_probability"]
    pair_adjacent_probability = (
        region_groups * (visits - 1) / math.comb(args.region_bits, 2)
    )
    intersections = hypergeometric_intersections(
        args.outer_bits, args.outer_weight, args.outer_weight
    )
    if abs(math.fsum(probability for _, probability in intersections) - 1.0) > 2e-14:
        raise AssertionError("intersection law is not normalized")
    union_regions_by_intersection = {
        intersection: 2 * args.outer_weight - intersection
        for intersection, _ in intersections
    }
    earliest_overlap_probability = math.fsum(
        probability
        * intersection
        / union_regions_by_intersection[intersection]
        for intersection, probability in intersections
    )

    leading_failure_rows = []
    for failures in range(1, args.maximum_leading_failures + 1):
        probability = math.fsum(
            intersection_probability
            * falling_ratio(
                intersection,
                union_regions_by_intersection[intersection],
                failures,
            )
            * pair_region_failure**failures
            for intersection, intersection_probability in intersections
        )
        leading_failure_rows.append(
            {
                "leading_failed_occupied_regions": failures,
                "probability": probability,
                "bits": -math.log2(probability),
            }
        )

    prefix_rows = []
    one_word_denominator = math.comb(args.outer_bits, args.outer_weight)
    for prefix in args.prefix_regions:
        actual = prefix_nonactivation_probability(
            args.outer_bits,
            args.outer_weight,
            args.outer_weight,
            prefix,
            pair_region_failure,
            intersections,
        )
        one_word_avoidance = (
            math.comb(args.outer_bits - prefix, args.outer_weight)
            / one_word_denominator
            if args.outer_bits - prefix >= args.outer_weight
            else 0.0
        )
        ideal = one_word_avoidance**2
        prefix_rows.append(
            {
                "prefix_regions": prefix,
                "nonactivation_probability": actual,
                "nonactivation_bits": -math.log2(actual) if actual else math.inf,
                "ideal_no_cancellation_probability": ideal,
                "ideal_no_cancellation_bits": (
                    -math.log2(ideal) if ideal else math.inf
                ),
                "cancellation_loss_bits": (
                    math.log2(actual / ideal) if ideal else math.inf
                ),
            }
        )

    modeled_spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    spectrum_log2 = float(modeled_spectrum[args.outer_weight]) / LOG2
    pair_multiplicity_log2 = (
        math.log2(math.comb(args.outer_blocks, 2)) + 2.0 * spectrum_log2
    )
    for row in prefix_rows:
        row["modeled_union_margin_bits_before_output_tail"] = (
            row["nonactivation_bits"] - pair_multiplicity_log2
        )
    all_never_activate_probability = (
        pair_region_failure**args.outer_weight / one_word_denominator
    )

    return {
        "schema": "riffle-fieldcheckpoint-cancellation-diagnostic-v1",
        "probability_space": {
            "outer_supports": (
                f"two independent uniform {args.outer_weight}-subsets of "
                f"{args.outer_bits} regions"
            ),
            "region_positions": (
                f"independent uniform supports inside {args.region_bits}-bit regions"
            ),
            "event": (
                "starting from zero, every occupied epoch-lane seen so far "
                "has even input parity"
            ),
        },
        "parameters": vars(args) | {"output": str(args.output)},
        "mean_live_weight_crossing": {
            "distance": args.distance,
            "first_prefix_regions_with_remaining_live_mean_at_most_distance": (
                math.ceil(
                    args.outer_bits
                    - 2.0 * args.distance / args.region_bits
                )
            ),
        },
        "epoch_nonactivation": epoch_rows,
        "region_nonactivation": region_rows,
        "weight_two_region": {
            "nonactivation_probability": pair_region_failure,
            "nonactivation_bits": -math.log2(pair_region_failure),
            "adjacent_pair_probability": pair_adjacent_probability,
            "adjacent_pair_bits": -math.log2(pair_adjacent_probability),
            "adjacent_conditional_on_nonactivation": (
                pair_adjacent_probability / pair_region_failure
            ),
            "conditional_mean_output_weight": (visits + 1) / 3.0,
        },
        "two_outer_words": {
            "mean_intersection": math.fsum(
                intersection * probability
                for intersection, probability in intersections
            ),
            "intersection_mode": max(intersections, key=lambda item: item[1])[0],
            "probability_earliest_occupied_region_is_overlap": (
                earliest_overlap_probability
            ),
            "probability_first_occupied_region_fails": (
                earliest_overlap_probability * pair_region_failure
            ),
            "first_occupied_region_failure_bits": -math.log2(
                earliest_overlap_probability * pair_region_failure
            ),
            "leading_failure_rows": leading_failure_rows,
            "prefix_nonactivation_rows": prefix_rows,
            "probability_never_activates": all_never_activate_probability,
            "never_activates_bits": -math.log2(all_never_activate_probability),
            "modeled_pair_multiplicity_log2": pair_multiplicity_log2,
            "modeled_never_activates_union_margin_bits": (
                -math.log2(all_never_activate_probability)
                - pair_multiplicity_log2
            ),
        },
        "scope": (
            "Exact activation-event probabilities for the stated support law. "
            "The receipt does not include the tilted output-weight tail after "
            "activation and is not a complete distance certificate."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--outer-weight", type=int, default=38)
    parser.add_argument("--outer-blocks", type=int, default=8192)
    parser.add_argument("--distance", type=int, default=188743)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=1024)
    parser.add_argument("--epochs-per-region", type=int, default=8)
    parser.add_argument("--region-bits", type=int, default=8192)
    parser.add_argument("--maximum-local-weight", type=int, default=16)
    parser.add_argument("--maximum-leading-failures", type=int, default=6)
    parser.add_argument(
        "--prefix-regions",
        type=int,
        nargs="+",
        default=(32, 64, 128, 192, 209, 210),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "weight_two_region_nonactivation_bits,"
        f"{payload['weight_two_region']['nonactivation_bits']:.12f}"
    )
    print(
        "two_word_last_prefix_nonactivation_bits,"
        f"{payload['two_outer_words']['prefix_nonactivation_rows'][-1]['nonactivation_bits']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
