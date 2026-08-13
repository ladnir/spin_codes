#!/usr/bin/env python3
"""Retune one worst sampled residual under both generalized outer branches."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from packet_group_outer_profile import optimize_spectrum_outer, total_weight_outer
from packet_group_profile_bound import (
    outer_probability,
    refine_inner_shaped_density,
    refine_inner_sparse_fugacities,
    split_cap_table,
    tune_inner_density,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-audit", type=Path, required=True)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--rank", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fast-outer", action="store_true")
    parser.add_argument("--skip-linear", action="store_true")
    parser.add_argument("--skip-spectrum", action="store_true")
    parser.add_argument("--screen-iterations", type=int, default=8)
    parser.add_argument("--coordinate-iterations", type=int, default=2)
    parser.add_argument("--inner-checkpoint", type=Path)
    parser.add_argument(
        "--sparse-refine",
        action="store_true",
        help="pin absent classes to the zero boundary and coordinate-refine occupied classes",
    )
    args = parser.parse_args()
    audits = json.loads(args.sample_audit.read_text(encoding="utf-8"))
    audit = next(row for row in audits if int(row["group_bits"]) == args.group)
    profile = [int(value) for value in audit["worst"][args.rank - 1]["profile"]]
    if args.inner_checkpoint:
        checkpoint_source = json.loads(
            args.inner_checkpoint.read_text(encoding="utf-8")
        )
        if [int(value) for value in checkpoint_source["profile"]] != profile:
            raise SystemExit("residual profile: inner checkpoint profile mismatch")
        inner = (
            float(checkpoint_source["inner_probability_log2"]),
            float(checkpoint_source["pole"]),
            np.asarray(checkpoint_source["fugacities"]),
            checkpoint_source["inner_details"],
            checkpoint_source.get("inner_parameters", []),
            checkpoint_source.get("inner_tuning_method", "checkpoint"),
        )
    else:
        split_caps = split_cap_table()
        initial = tune_inner_density(
            args.group,
            profile,
            split_caps,
            poles=(0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
            exponents=(0.25, 0.5, 0.75, 1.0),
            iterations=args.screen_iterations,
        )
        refined = refine_inner_shaped_density(
            args.group,
            profile,
            split_caps,
            initial_pole=initial[1],
            initial_exponent=initial[2],
            coordinate_iterations=args.coordinate_iterations,
            passes=1,
            witness_iterations=args.screen_iterations,
            linear_bound=36.0,
            quadratic_bound=72.0,
        )
        screened_inner = min(
            (
                (initial[0], initial[1], initial[3], initial[4], [initial[2], 0.0, 0.0, math.log(initial[1])], "density_grid"),
                (refined[0], refined[1], refined[2], refined[3], refined[4].tolist(), "shaped_refinement"),
            ),
            key=lambda row: row[0],
        )
        if args.sparse_refine:
            sparse = refine_inner_sparse_fugacities(
                args.group,
                profile,
                split_caps,
                initial_pole=screened_inner[1],
                initial_fugacities=screened_inner[2],
                witness_iterations=args.screen_iterations,
                verbose=True,
            )
            if sparse[0] < screened_inner[0]:
                screened_inner = (
                    sparse[0],
                    sparse[1],
                    sparse[2],
                    sparse[3],
                    {
                        "trace": sparse[4],
                        "evaluations": sparse[5],
                        "anchor": sparse[6],
                    },
                    "sparse_fugacity_refinement",
                )
        # One full Collatz verification of the selected frozen parameters.
        from packet_group_profile_bound import inner_probability
        verified_probability, verified_details = inner_probability(
            args.group,
            profile,
            screened_inner[1],
            screened_inner[2],
            split_caps,
            iterations=40,
        )
        inner = (
            verified_probability,
            screened_inner[1],
            screened_inner[2],
            verified_details,
            screened_inner[4],
            screened_inner[5] + "_screened_then_verified",
        )
    checkpoint = {
        "status": "DIAGNOSTIC_RESIDUAL_INNER_CHECKPOINT",
        "group_bits": args.group,
        "profile": profile,
        "pole": inner[1],
        "fugacities": inner[2].tolist(),
        "inner_probability_log2": inner[0],
        "inner_details": inner[3],
        "inner_parameters": inner[4],
        "inner_tuning_method": inner[5],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(checkpoint, indent=2) + "\n", encoding="utf-8")

    physical_weight = sum(weight * count for weight, count in enumerate(profile))
    total_outer, total_result, total_interval = total_weight_outer(physical_weight)

    if args.skip_linear:
        linear_outer = math.inf
        linear_details = {}
    else:
        linear_outer, linear_details = outer_probability(
            args.group, profile, fast=args.fast_outer
        )
    if args.skip_spectrum:
        spectrum_outer = math.inf
        spectrum_result = None
        spectrum_components = None
        spectrum_log = math.inf
        spectrum_graph = math.inf
        spectrum_log_variables = []
    else:
        spectrum_log, spectrum_result, spectrum_components = optimize_spectrum_outer(
            args.group, profile
        )
        spectrum_log_variables = [0.0] + spectrum_result.x[: args.group].tolist()
        spectrum_graph = 128.0 * max(
            spectrum_log_variables[new] - spectrum_log_variables[old]
            for old in range(args.group + 1)
            for new in range(args.group + 1)
            if abs(new - old) <= 1
        ) / math.log(2.0)
        spectrum_outer = spectrum_log / math.log(2.0) + spectrum_graph
    if linear_outer <= spectrum_outer and linear_outer <= total_outer:
        outer_type = "linear_bl"
        outer = linear_outer
        outer_details = linear_details
    elif spectrum_outer <= total_outer:
        outer_type = "total_spectrum"
        outer = spectrum_outer
        outer_details = {
            "unpunctured_outer_log2": spectrum_log / math.log(2.0),
            "graph_replacement_log2": spectrum_graph,
            "log_variables": spectrum_log_variables,
            "outer_point": spectrum_result.x.tolist(),
            "log_alpha": spectrum_components[1],
            "log_enumerator": spectrum_components[2],
        }
    else:
        outer_type = "total_weight"
        outer = total_outer
        outer_details = {
            "physical_total_weight": physical_weight,
            "pre_replacement_weight_interval": list(total_interval),
            "outer_log_pole": float(total_result.x),
        }
    combined = outer + inner[0]
    report = {
        "status": "DIAGNOSTIC_RETUNED_GENERALIZED_G_RESIDUAL",
        "group_bits": args.group,
        "source_audit": str(args.sample_audit),
        "source_rank": args.rank,
        "profile": profile,
        "pole": inner[1],
        "fugacities": inner[2].tolist(),
        "inner_probability_log2": inner[0],
        "inner_details": inner[3],
        "inner_parameters": inner[4],
        "inner_tuning_method": inner[5],
        "outer_type": outer_type,
        "outer_profile_log2": outer,
        "outer_details": outer_details,
        "combined_log2": combined,
        "target_log2": float(audit["target_log2"]),
        "margin_log2": float(audit["target_log2"]) - combined,
    }
    rendered = json.dumps(report, indent=2)
    print(rendered)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
