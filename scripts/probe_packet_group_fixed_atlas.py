#!/usr/bin/env python3
"""Tune reusable fixed profile witnesses for a generalized packet width."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import (
    D,
    N,
    atom_count,
    fixed_total_weight_outer,
    normalization_log2,
    optimize_spectrum_outer,
)
from packet_group_profile_bound import (
    density_fugacities,
    inner_probability,
    refine_inner_full_fugacities,
    refine_inner_shaped_density,
    refine_inner_sparse_fugacities,
    split_cap_table,
    tune_inner_density,
    outer_probability,
)
from probe_packet_group_global_mixture_cover import full_offset_outer_constant
from certify_packet8_hard_face_drive_inner import exact_float
from probe_packet_group_exact_graph_puncture_outer import (
    optimize_q as optimize_exact_graph_puncture_q,
    outward_outer as outward_exact_graph_puncture_outer,
)


def parse_named_profile(text: str, group_bits: int) -> tuple[str, list[int]]:
    if ":" not in text:
        raise ValueError("profile must have form name:c0,c1,...,cg")
    name, raw = text.split(":", 1)
    profile = [int(value) for value in raw.split(",")]
    if len(profile) != group_bits + 1 or any(value < 0 for value in profile):
        raise ValueError(f"invalid g={group_bits} profile: {text}")
    if sum(profile) != atom_count(group_bits):
        raise ValueError(f"profile mass mismatch: {text}")
    return name, profile


def tune_fixed_witness(
    group_bits: int,
    name: str,
    profile: list[int],
    split_caps,
    *,
    screen_iterations: int,
    coordinate_iterations: int,
    exact_graph_puncture: bool = False,
    full_coordinate_passes: int = 2,
    full_coordinate_fine_steps: bool = False,
) -> dict:
    sparse_fast = screen_iterations == 1 and any(count == 0 for count in profile)
    if sparse_fast:
        initial = tune_inner_density(
            group_bits,
            profile,
            split_caps,
            poles=(0.2, 0.5, 0.8),
            exponents=(0.25, 0.5, 0.75),
            iterations=1,
        )
    else:
        initial = tune_inner_density(
            group_bits,
            profile,
            split_caps,
            poles=(0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
            exponents=(0.25, 0.5, 0.75, 1.0),
            iterations=screen_iterations,
        )
    candidates = [
        (
            initial[0],
            initial[1],
            initial[3],
            initial[4],
            "density_grid",
        ),
    ]
    positive_counts = [count for count in profile if count]
    if min(positive_counts) > 16:
        refined = refine_inner_shaped_density(
            group_bits,
            profile,
            split_caps,
            initial_pole=initial[1],
            initial_exponent=initial[2],
            coordinate_iterations=coordinate_iterations,
            passes=1,
            witness_iterations=screen_iterations,
            linear_bound=36.0,
            quadratic_bound=72.0,
        )
        candidates.append(
            (
                refined[0],
                refined[1],
                refined[2],
                refined[3],
                "shaped_density",
            )
        )
    screened = min(candidates, key=lambda row: row[0])
    if any(count == 0 for count in profile):
        sparse = refine_inner_sparse_fugacities(
            group_bits,
            profile,
            split_caps,
            initial_pole=screened[1],
            initial_fugacities=screened[2],
            log_steps=(math.log(16.0), math.log(4.0), math.log(2.0)),
            witness_iterations=screen_iterations,
            # A 64-bit block contains 32 atoms at g=2.  Floors below roughly
            # 2^-32 underflow when a point cap multiplies one factor per atom.
            # This dyadic floor is still strongly face-localized and remains
            # directly representable by the outward verifier.
            zero_floor=2.0**-16,
        )
        candidates.append(
            (
                sparse[0],
                sparse[1],
                sparse[2],
                sparse[3],
                "support_sparse",
            )
        )
        screened = min(candidates, key=lambda row: row[0])

        # The sparse search uses a small positive floor only to keep its
        # diagnostic coordinate steps numerically stable.  Once it selects
        # the occupied coordinates, project the unused coordinates onto the
        # exact zero face and recompute.  All local coefficients are
        # nonnegative, so deleting unused-class monomials cannot enlarge the
        # MGF.  The resulting frozen witness is exactly support-eligible.
        exact_face = np.asarray(screened[2], dtype=np.float64).copy()
        exact_face[np.asarray(profile) == 0] = 0.0
        exact_probability, exact_details = inner_probability(
            group_bits,
            profile,
            screened[1],
            exact_face,
            split_caps,
            iterations=screen_iterations,
        )
        screened = (
            exact_probability,
            screened[1],
            exact_face,
            exact_details,
            "support_sparse_exact_face",
        )

    inner_value, inner_details = inner_probability(
        group_bits,
        profile,
        screened[1],
        screened[2],
        split_caps,
        iterations=40,
    )
    physical_weight = sum(index * count for index, count in enumerate(profile))
    linear_outer_value, linear_outer_details = outer_probability(
        group_bits, profile, fast=False
    )
    vector = np.asarray(profile, dtype=np.float64)

    def occupied_dot(left: np.ndarray, right: np.ndarray) -> float:
        occupied = left != 0.0
        return float(left[occupied] @ right[occupied])
    total_outer_value, total_outer_details = fixed_total_weight_outer(physical_weight)
    if total_outer_value < linear_outer_value:
        outer_type = "total_weight"
        outer_charge = (
            np.arange(group_bits + 1, dtype=np.float64)
            * float(total_outer_details["charge_per_bit"])
        )
        outer_constant = float(total_outer_details["constant_log2"])
        outer_details = total_outer_details
        outer_value = total_outer_value
    else:
        outer_type = "linear_bl"
        outer_log_variables = np.asarray(
            linear_outer_details["log_variables"], dtype=np.float64
        )
        outer_charge = outer_log_variables / math.log(2.0)
        outer_constant = linear_outer_value + float(vector @ outer_charge)
        outer_details = linear_outer_details
        outer_value = linear_outer_value
    inner_constant = (
        float(inner_details["inner_mgf_log2"])
        - D * math.log2(float(screened[1]))
    )
    fugacities = np.asarray(screened[2], dtype=np.float64)
    with np.errstate(divide="ignore"):
        charge = outer_charge + np.log2(fugacities)
    constant = outer_constant + inner_constant
    normalization = float(normalization_log2(group_bits, vector)[0])
    combined = constant - occupied_dot(vector, charge) - normalization
    target = -40.0 - math.log2(profile_count(group_bits, N))
    # The profile-shaped total-spectrum optimization is materially more
    # expensive than the linear-BL and total-weight branches, so reserve it
    # for profiles that do not already meet the uniform allocation.  Its
    # frozen Cauchy parameters still define one globally valid affine outer
    # witness and are therefore reusable by convex-cell certificates.
    if combined > target:
        spectrum_log, spectrum_result, spectrum_components = optimize_spectrum_outer(
            group_bits, profile
        )
        spectrum_log_variables = np.concatenate(
            ([0.0], spectrum_result.x[:group_bits])
        )
        spectrum_graph = 128.0 * max(
            spectrum_log_variables[new] - spectrum_log_variables[old]
            for old in range(group_bits + 1)
            for new in range(group_bits + 1)
            if abs(new - old) <= 1
        ) / math.log(2.0)
        spectrum_outer_value = spectrum_log / math.log(2.0) + spectrum_graph
        if spectrum_outer_value < outer_value:
            outer_type = "total_spectrum"
            outer_charge = spectrum_log_variables / math.log(2.0)
            outer_constant = spectrum_outer_value + float(vector @ outer_charge)
            outer_details = {
                "unpunctured_outer_log2": spectrum_log / math.log(2.0),
                "graph_replacement_log2": spectrum_graph,
                "log_variables": spectrum_log_variables.tolist(),
                "outer_point": spectrum_result.x.tolist(),
                "log_alpha": spectrum_components[1],
                "log_enumerator": spectrum_components[2],
            }
            outer_value = spectrum_outer_value
            charge = outer_charge + np.log2(fugacities)
            constant = outer_constant + inner_constant
            combined = constant - occupied_dot(vector, charge) - normalization
    if exact_graph_puncture and group_bits == 4:
        graph_log_q = optimize_exact_graph_puncture_q(profile)
        graph_q = exact_float(math.exp(graph_log_q))
        graph_interval = outward_exact_graph_puncture_outer(profile, graph_q)
        graph_outer_value = float(graph_interval.hi)
        if graph_outer_value < outer_value:
            outer_type = "exact_graph_puncture_total_weight"
            outer_charge = (
                np.arange(group_bits + 1, dtype=np.float64) * math.log2(float(graph_q))
            )
            outer_constant = graph_outer_value + float(vector @ outer_charge)
            outer_details = {
                "pole_exact": f"{graph_q.numerator}/{graph_q.denominator}",
                "log_pole": math.log(float(graph_q)),
                "outer_log2_interval": [str(graph_interval.lo), str(graph_interval.hi)],
                "source_graph_spectrum": "ebch128_graph24_spectrum.csv",
                "source_graph_spectrum_sha256": hashlib.sha256(
                    Path(__file__).with_name("ebch128_graph24_spectrum.csv").read_bytes()
                ).hexdigest(),
            }
            outer_value = graph_outer_value
            charge = outer_charge + np.log2(fugacities)
            constant = outer_constant + inner_constant
            combined = constant - occupied_dot(vector, charge) - normalization
    if combined > target and all(count > 0 for count in profile):
        full = refine_inner_full_fugacities(
            group_bits,
            profile,
            split_caps,
            initial_pole=float(screened[1]),
            initial_fugacities=fugacities,
            log_steps=(
                (
                    math.log(16.0),
                    math.log(4.0),
                    math.log(2.0),
                    math.log(math.sqrt(2.0)),
                    math.log(2.0**0.25),
                    math.log(2.0**0.125),
                )
                if full_coordinate_fine_steps
                else (
                    math.log(16.0),
                    math.log(4.0),
                    math.log(2.0),
                    math.log(math.sqrt(2.0)),
                )
            ),
            passes=full_coordinate_passes,
            witness_iterations=screen_iterations,
        )
        full_inner_value, full_inner_details = inner_probability(
            group_bits,
            profile,
            full[1],
            full[2],
            split_caps,
            iterations=40,
        )
        full_inner_constant = (
            float(full_inner_details["inner_mgf_log2"])
            - D * math.log2(float(full[1]))
        )
        full_charge = outer_charge + np.log2(full[2])
        full_constant = outer_constant + full_inner_constant
        full_combined = (
            full_constant - float(vector @ full_charge) - normalization
        )
        if full_combined < combined:
            inner_value = full_inner_value
            inner_details = {
                **full_inner_details,
                "full_coordinate_trace": full[4],
                "full_coordinate_evaluations": full[5],
                "full_coordinate_anchor": full[6],
            }
            fugacities = np.asarray(full[2], dtype=np.float64)
            inner_constant = full_inner_constant
            charge = full_charge
            constant = full_constant
            combined = full_combined
            screened = (
                full_inner_value,
                full[1],
                fugacities,
                full_inner_details,
                "full_coordinate",
            )
    return {
        "name": name,
        "group_bits": group_bits,
        "profile": profile,
        "physical_weight": physical_weight,
        "support": [index for index, value in enumerate(fugacities) if value > 0.0],
        "pole": float(screened[1]),
        "fugacities": fugacities.tolist(),
        "outer_type": outer_type,
        "outer_charge": outer_charge.tolist(),
        "outer_details": outer_details,
        "linear_bl_outer_log2": linear_outer_value,
        "total_weight_outer_log2": total_outer_value,
        "outer_log2": outer_value,
        "outer_constant_log2": outer_constant,
        "inner_constant_log2": inner_constant,
        "constant_log2": constant,
        "charge": charge.tolist(),
        "normalization_log2": normalization,
        "combined_log2": combined,
        "target_log2": target,
        "margin_bits": target - combined,
        "inner_probability_log2": inner_value,
        "inner_details": inner_details,
        "tuning_method": screened[4],
        "status": "DIAGNOSTIC_BINARY64_GENERALIZED_FIXED_PROFILE_WITNESS",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--profile", action="append", required=True)
    parser.add_argument("--screen-iterations", type=int, default=8)
    parser.add_argument("--coordinate-iterations", type=int, default=3)
    parser.add_argument("--full-coordinate-passes", type=int, default=2)
    parser.add_argument("--full-coordinate-fine-steps", action="store_true")
    parser.add_argument(
        "--exact-graph-puncture",
        action="store_true",
        help="enable the exact graph-spectrum/puncture-averaged g=4 total-weight branch",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        args.screen_iterations <= 0
        or args.coordinate_iterations <= 0
        or args.full_coordinate_passes <= 0
    ):
        raise SystemExit("fixed atlas: invalid iteration count")
    profiles = [parse_named_profile(text, args.group) for text in args.profile]
    split_caps = split_cap_table()
    rows = []
    for index, (name, profile) in enumerate(profiles, 1):
        row = tune_fixed_witness(
            args.group,
            name,
            profile,
            split_caps,
            screen_iterations=args.screen_iterations,
            coordinate_iterations=args.coordinate_iterations,
            exact_graph_puncture=args.exact_graph_puncture,
            full_coordinate_passes=args.full_coordinate_passes,
            full_coordinate_fine_steps=args.full_coordinate_fine_steps,
        )
        rows.append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
        print(
            f"witness={index}/{len(profiles)} name={name} "
            f"combined={row['combined_log2']:.9f} margin={row['margin_bits']:.9f}",
            flush=True,
        )
    print("status=DIAGNOSTIC_BINARY64_GENERALIZED_FIXED_PROFILE_ATLAS")


if __name__ == "__main__":
    main()
