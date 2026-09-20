#!/usr/bin/env python3
"""Audit the local merge inequality for shared four-block packet groups.

At a fixed Chernoff point, let M_k be the two-state transfer through one epoch
that contains k candidate bits.  Each candidate bit is independently active
with the Bernoulli probability used by the regular outer envelope.

Consider two packet groups of widths r and s in different epochs.  Conditional
on the two epoch backgrounds u and v, random packet order gives the split law

  (M_{u+r} M_{v+s} + M_{u+s} M_{v+r}) / 2.

Merge the groups and place the merged packet at either old position.  The
corresponding packed law is

  (M_{u+r+s} M_v + M_u M_{v+r+s}) / 2.

If the packed law entrywise dominates the split law for every u,v and every
r+s<=4, repeated merging proves that the maximally packed profile dominates
every group profile.  Prefix and suffix transfers preserve the entrywise
order because all transfer matrices are nonnegative.

This script audits the finite inequalities in floating point and evaluates the
resulting universal packed-profile bound.  It is not an outward-rounded proof.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_ldpcsplitstate_occupation_ladder import (
    DEFAULT_ACTIVATION,
    epoch_transfers,
)
from analyze_riffle_ldpcsplitstate_shared_groups import (
    averaged_epoch_transfers,
    evaluate_profile_point,
    reference_points,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/packing_domination_ladder_audit.json"
)
DEFAULT_QUAD_RECEIPT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/shared_group_occupation128_quad_optimized.json"
)


def parse_int_list(text: str) -> list[int]:
    return [int(item) for item in text.split(",") if item]


def sanitize_nonfinite(value):
    if isinstance(value, float) and not math.isfinite(value):
        return "infinity" if value > 0 else "minus_infinity"
    if isinstance(value, dict):
        return {key: sanitize_nonfinite(item) for key, item in value.items()}
    if isinstance(value, list):
        return [sanitize_nonfinite(item) for item in value]
    return value


def maximally_packed_profile(occupation: int) -> tuple[int, int, int, int]:
    quads, remainder = divmod(occupation, 4)
    profile = [0, 0, 0, quads]
    if remainder:
        profile[remainder - 1] = 1
    return tuple(profile)


def all_block_sets_log2(occupation: int) -> float:
    return (
        math.lgamma(8193)
        - math.lgamma(occupation + 1)
        - math.lgamma(8193 - occupation)
    ) / math.log(2.0)


def merge_audit(
    matrices: list[np.ndarray], occupation: int
) -> dict[str, object]:
    worst_absolute = (math.inf, None)
    worst_scaled = (math.inf, None)
    violations = 0
    comparisons = 0
    maximum_ratio = (1.0, None)
    operation_rows = []
    for left_width, right_width in ((1, 1), (1, 2), (1, 3), (2, 2)):
        operation_maximum = (1.0, None)
        merged_width = left_width + right_width
        maximum_background = occupation - merged_width
        for left_background in range(maximum_background + 1):
            for right_background in range(maximum_background + 1):
                for bridge_row in range(2):
                    for bridge_column in range(2):
                        bridge = np.zeros((2, 2), dtype=np.float64)
                        bridge[bridge_row, bridge_column] = 1.0
                        packed_left = (
                            matrices[left_background + merged_width]
                            @ bridge
                            @ matrices[right_background]
                        )
                        packed_right = (
                            matrices[left_background]
                            @ bridge
                            @ matrices[right_background + merged_width]
                        )
                        split_left = (
                            matrices[left_background + left_width]
                            @ bridge
                            @ matrices[right_background + right_width]
                        )
                        split_right = (
                            matrices[left_background + right_width]
                            @ bridge
                            @ matrices[right_background + left_width]
                        )
                        difference = 0.5 * (
                            packed_left + packed_right - split_left - split_right
                        )
                        scale = 0.5 * (
                            np.abs(packed_left)
                            + np.abs(packed_right)
                            + np.abs(split_left)
                            + np.abs(split_right)
                        )
                        for row in range(2):
                            for column in range(2):
                                comparisons += 1
                                value = float(difference[row, column])
                                local_scale = float(scale[row, column])
                                scaled = value / local_scale if local_scale else 0.0
                                witness = {
                                    "left_width": left_width,
                                    "right_width": right_width,
                                    "left_background": left_background,
                                    "right_background": right_background,
                                    "bridge_entry": [bridge_row, bridge_column],
                                    "entry": [row, column],
                                    "difference": value,
                                    "comparison_scale": local_scale,
                                    "scaled_difference": scaled,
                                }
                                if value < worst_absolute[0]:
                                    worst_absolute = (value, witness)
                                if scaled < worst_scaled[0]:
                                    worst_scaled = (scaled, witness)
                                if scaled < -1e-12:
                                    violations += 1
                                packed_value = 0.5 * float(
                                    packed_left[row, column]
                                    + packed_right[row, column]
                                )
                                split_value = 0.5 * float(
                                    split_left[row, column]
                                    + split_right[row, column]
                                )
                                ratio = (
                                    split_value / packed_value
                                    if packed_value > 0.0
                                    else (math.inf if split_value > 0.0 else 1.0)
                                )
                                ratio_witness = witness | {
                                    "operation": "merge",
                                    "split_value": split_value,
                                    "packed_value": packed_value,
                                    "split_over_packed": ratio,
                                }
                                if ratio > operation_maximum[0]:
                                    operation_maximum = (ratio, ratio_witness)
                                if ratio > maximum_ratio[0]:
                                    maximum_ratio = (ratio, ratio_witness)
        operation_rows.append(
            {
                "operation": "merge",
                "source_widths": [left_width, right_width],
                "target_widths": [0, merged_width],
                "maximum_split_over_packed": operation_maximum[0],
                "maximum_log2_surcharge": math.log2(operation_maximum[0]),
                "witness": operation_maximum[1],
            }
        )

    for source, target in (((2, 3), (1, 4)), ((3, 3), (2, 4))):
        left_width, right_width = source
        packed_left_width, packed_right_width = target
        operation_maximum = (1.0, None)
        maximum_background = occupation - left_width - right_width
        for left_background in range(maximum_background + 1):
            for right_background in range(maximum_background + 1):
                for bridge_row in range(2):
                    for bridge_column in range(2):
                        bridge = np.zeros((2, 2), dtype=np.float64)
                        bridge[bridge_row, bridge_column] = 1.0
                        split_left = (
                            matrices[left_background + left_width]
                            @ bridge
                            @ matrices[right_background + right_width]
                        )
                        split_right = (
                            matrices[left_background + right_width]
                            @ bridge
                            @ matrices[right_background + left_width]
                        )
                        packed_left = (
                            matrices[left_background + packed_left_width]
                            @ bridge
                            @ matrices[right_background + packed_right_width]
                        )
                        packed_right = (
                            matrices[left_background + packed_right_width]
                            @ bridge
                            @ matrices[right_background + packed_left_width]
                        )
                        for row in range(2):
                            for column in range(2):
                                split_value = 0.5 * float(
                                    split_left[row, column]
                                    + split_right[row, column]
                                )
                                packed_value = 0.5 * float(
                                    packed_left[row, column]
                                    + packed_right[row, column]
                                )
                                ratio = (
                                    split_value / packed_value
                                    if packed_value > 0.0
                                    else (math.inf if split_value > 0.0 else 1.0)
                                )
                                witness = {
                                    "operation": "redistribute",
                                    "source_widths": list(source),
                                    "target_widths": list(target),
                                    "left_background": left_background,
                                    "right_background": right_background,
                                    "bridge_entry": [bridge_row, bridge_column],
                                    "entry": [row, column],
                                    "split_value": split_value,
                                    "packed_value": packed_value,
                                    "split_over_packed": ratio,
                                }
                                if ratio > operation_maximum[0]:
                                    operation_maximum = (ratio, witness)
                                if ratio > maximum_ratio[0]:
                                    maximum_ratio = (ratio, witness)
        operation_rows.append(
            {
                "operation": "redistribute",
                "source_widths": list(source),
                "target_widths": list(target),
                "maximum_split_over_packed": operation_maximum[0],
                "maximum_log2_surcharge": math.log2(operation_maximum[0]),
                "witness": operation_maximum[1],
            }
        )
    maximum_ratio_value = maximum_ratio[0]
    maximum_log2_value = math.log2(maximum_ratio_value)
    return {
        "comparisons": comparisons,
        "scaled_violations_below_minus_1e_12": violations,
        "minimum_absolute_difference": worst_absolute[0],
        "minimum_absolute_witness": worst_absolute[1],
        "minimum_scaled_difference": worst_scaled[0],
        "minimum_scaled_witness": worst_scaled[1],
        "maximum_local_split_over_packed": (
            maximum_ratio_value if math.isfinite(maximum_ratio_value) else "infinity"
        ),
        "maximum_local_log2_surcharge": (
            maximum_log2_value if math.isfinite(maximum_log2_value) else "infinity"
        ),
        "finite_local_scalar_surcharge": math.isfinite(maximum_log2_value),
        "_maximum_local_log2_surcharge_float": maximum_log2_value,
        "maximum_local_ratio_witness": maximum_ratio[1],
        "operation_rows": operation_rows,
    }


def quad_point(path: Path) -> tuple[float, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for row in payload["occupation_rows"][0]["profiles"]:
        if row["profile_n1_n2_n3_n4"] == [0, 0, 0, 32]:
            return float(row["log_surprisal"]), float(row["placement_probability"])
    raise ValueError("quad-packed occupation-128 row not found")


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    occupations = parse_int_list(args.occupations)
    references = reference_points(args.reference_receipts)
    references[128] = quad_point(args.quad_receipt)
    maximum = max(occupations)
    activation = json.loads(args.activation.read_text(encoding="utf-8"))
    activation_upper = [0.0] * (maximum + 1)
    activation_upper[0] = 1.0
    for item in activation["by_total_weight"]:
        weight = int(item["total_weight"])
        if weight > maximum:
            break
        activation_upper[weight] = float(
            item["maximum_distinct_conditioned_upper_bound"]
        )

    target_distance = math.floor(args.relative_distance * (1 << 21))
    rows = []
    for occupation in occupations:
        log_surprisal, probability = references[occupation]
        placement_logit = math.log(probability / (1.0 - probability))
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        epoch = epoch_transfers(
            z=z,
            distance=args.constituent_distance,
            moment_order=args.live_moment_order,
            activation_upper=activation_upper,
            maximum=maximum,
        )
        matrices = averaged_epoch_transfers(epoch, probability, occupation)
        audit = merge_audit(matrices, occupation)
        profile = maximally_packed_profile(occupation)
        packed = evaluate_profile_point(
            profile=profile,
            log_surprisal=log_surprisal,
            placement_logit=placement_logit,
            activation_upper=activation_upper,
            constituent_distance=args.constituent_distance,
            live_moment_order=args.live_moment_order,
            target_distance=target_distance,
        )
        total_set_bits = all_block_sets_log2(occupation)
        universal_margin = float(packed["aggregate_margin_bits"]) - (
            total_set_bits - float(packed["block_set_multiplicity_bits"])
        )
        maximum_repacking_moves = math.ceil(1.5 * occupation)
        local_surcharge = float(
            audit.pop("_maximum_local_log2_surcharge_float")
        )
        repacking_surcharge = (
            maximum_repacking_moves
            * local_surcharge
        )
        surcharged_margin = universal_margin - repacking_surcharge
        row = {
            "occupation": occupation,
            "log_surprisal": log_surprisal,
            "placement_probability": probability,
            "maximally_packed_profile": list(profile),
            "packed_profile_margin_bits": packed["aggregate_margin_bits"],
            "packed_profile_block_set_bits": packed[
                "block_set_multiplicity_bits"
            ],
            "all_active_block_sets_bits": total_set_bits,
            "universal_margin_bits_if_merge_inequality_holds": universal_margin,
            "maximum_repacking_moves": maximum_repacking_moves,
            "repacking_surcharge_bits": (
                repacking_surcharge
                if math.isfinite(repacking_surcharge)
                else "infinity"
            ),
            "universal_margin_bits_with_scalar_repacking_surcharge": (
                surcharged_margin if math.isfinite(surcharged_margin) else None
            ),
            "merge_audit": audit,
        }
        rows.append(row)
        print(
            f"occupation,{occupation},profile,{list(profile)},"
            f"universal_margin,{universal_margin:.6f},"
            f"surcharged_margin,{surcharged_margin},"
            f"local_surcharge,{audit['maximum_local_log2_surcharge']},"
            f"min_scaled_merge,{audit['minimum_scaled_difference']:.3e},"
            f"violations,{audit['scaled_violations_below_minus_1e_12']}",
            flush=True,
        )
    return {
        "schema": "riffle-ldpcsplitstate-packing-domination-audit-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "occupations": occupations,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
            "relative_distance": args.relative_distance,
            "target_distance": target_distance,
        },
        "rows": rows,
        "scope": (
            "Finite floating-point audit of the two-epoch merge inequality. "
            "A proof requires outward rounding or a symbolic derivation. The "
            "universal margin replaces the packed profile multiplicity by all "
            "active block sets, so no group-profile enumeration remains if "
            "the merge inequality is certified."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--quad-receipt", type=Path, default=DEFAULT_QUAD_RECEIPT)
    parser.add_argument("--occupations", default="2,4,8,16,32,64,80,96,112,128")
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
        ],
    )
    args = parser.parse_args()
    payload = sanitize_nonfinite(evaluate(args))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
