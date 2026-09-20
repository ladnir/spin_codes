#!/usr/bin/env python3
"""Build a binary64 finite interval cover for occupations 32 through 8448.

For a fixed RM2Sub witness (p,z,v), the logarithm of the finite first-moment
bound is convex in the integer occupation Q.  Therefore, endpoint bounds
cover every integer in an accepted interval.  Binary64 optimization proposes
the witnesses and evaluates the transfer, so the result is diagnostic.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from evaluate_bch250_rowlocal_fanout_finite_cover import (  # noqa: E402
    B,
    DEFAULT_SPECTRUM,
    DISTANCE,
    EPOCHS,
    K,
    L,
    N,
    density_excess_bits,
    log_binomial_mass,
    log_choose,
)


def finite_log_bound(
    occupation: int,
    *,
    candidate_probability: float,
    z: float,
    radius: float,
    matrix_prefactor: float,
    density_excess: mp.mpf,
    dimension: int,
    outer_rows: int,
    epochs: int,
    distance: int,
) -> mp.mpf:
    p = mp.mpf(candidate_probability)
    z_mp = mp.mpf(z)
    log_bound = log_choose(outer_rows, occupation)
    log_bound += occupation * mp.log(mp.power(2, dimension) - 1)
    log_bound -= occupation * mp.log1p(-mp.power(2, -B))
    log_bound -= B * log_binomial_mass(outer_rows, occupation, p)
    log_bound -= distance * mp.log(z_mp)
    log_bound += epochs * mp.log(mp.mpf(radius))
    log_bound += mp.log(mp.mpf(matrix_prefactor))
    log_bound += occupation * density_excess * mp.log(2)
    return log_bound


def witness(envelope, occupation: int, outer_rows: int) -> dict[str, object]:
    alpha = occupation / outer_rows
    optimized = envelope.optimize_variant(alpha, "three_state")
    p = float(optimized["candidate_probability"])
    surprisal = float(optimized["surprisal"])
    z = math.exp(-surprisal)
    transfer, _ = envelope.transfer(
        candidate_probability=p, surprisal=surprisal
    )
    eigenvalues, eigenvectors = np.linalg.eig(transfer)
    index = int(np.argmax(np.abs(eigenvalues)))
    radius = math.nextafter(float(abs(eigenvalues[index])), math.inf)
    vector = np.abs(np.real(eigenvectors[:, index]))
    vector /= float(np.max(vector))
    vector = np.maximum(vector, 1e-300)
    row_ratios = transfer @ vector / vector
    radius = max(radius, math.nextafter(float(np.max(row_ratios)), math.inf))
    prefactor = float(vector[0] * np.max(1.0 / vector))
    return {
        "candidate_probability": p,
        "z": z,
        "surprisal": surprisal,
        "radius_upper_binary64": radius,
        "collatz_vector_binary64": [float(value) for value in vector],
        "matrix_prefactor_binary64": prefactor,
        "optimizer_success": bool(optimized["optimizer_success"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--minimum-occupation", type=int, default=32)
    parser.add_argument("--minimum-output-weight", type=int, default=1)
    parser.add_argument("--maximum-output-weight", type=int, default=B)
    parser.add_argument("--outer-dimension", type=int, default=K)
    parser.add_argument("--outer-rows", type=int, default=L)
    parser.add_argument("--distance", type=int)
    parser.add_argument("--target-union-margin", type=float, default=40.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.minimum_occupation <= args.outer_rows:
        raise ValueError(
            f"minimum occupation must lie in [1,{args.outer_rows}]"
        )
    output_bits = B * args.outer_rows
    distance = (
        args.distance
        if args.distance is not None
        else math.floor(0.11 * output_bits)
    )
    if output_bits % 128:
        raise ValueError("output bits must be divisible by 128")
    epochs = output_bits // 128

    mp.mp.dps = 100
    if not 1 <= args.minimum_output_weight <= args.maximum_output_weight <= B:
        raise ValueError("invalid output-weight interval")
    density_excess, maximizing_weight = density_excess_bits(
        args.spectrum,
        args.outer_dimension,
        args.minimum_output_weight,
        args.maximum_output_weight,
    )
    envelope = load_envelope(args.selection, distance / output_bits)
    stack = [(args.minimum_occupation, args.outer_rows)]
    accepted = []
    rejected_singletons = []
    processed = 0
    while stack:
        lower, upper = stack.pop()
        processed += 1
        midpoint = (lower + upper) // 2
        candidate = witness(envelope, midpoint, args.outer_rows)
        endpoint_logs = []
        for occupation in (lower, upper):
            value = finite_log_bound(
                occupation,
                candidate_probability=float(candidate["candidate_probability"]),
                z=float(candidate["z"]),
                radius=float(candidate["radius_upper_binary64"]),
                matrix_prefactor=float(candidate["matrix_prefactor_binary64"]),
                density_excess=density_excess,
                dimension=args.outer_dimension,
                outer_rows=args.outer_rows,
                epochs=epochs,
                distance=distance,
            )
            endpoint_logs.append(float(value / mp.log(2)))
        maximum = max(endpoint_logs)
        if maximum < 0.0:
            accepted.append(
                {
                    "lower_occupation": lower,
                    "upper_occupation": upper,
                    "midpoint_occupation": midpoint,
                    "endpoint_log2_uppers": endpoint_logs,
                    "worst_endpoint_margin_bits": -maximum,
                    **candidate,
                }
            )
            continue
        if lower == upper:
            rejected_singletons.append(
                {
                    "occupation": lower,
                    "log2_upper": maximum,
                    **candidate,
                }
            )
            continue
        split = (lower + upper) // 2
        stack.append((split + 1, upper))
        stack.append((lower, split))

    accepted.sort(key=lambda row: int(row["lower_occupation"]))
    union_terms = []
    for row in accepted:
        width = int(row["upper_occupation"]) - int(row["lower_occupation"]) + 1
        union_terms.append(
            math.log2(width) + max(float(value) for value in row["endpoint_log2_uppers"])
        )
    if union_terms:
        maximum = max(union_terms)
        union_log2 = maximum + math.log2(
            sum(math.exp2(value - maximum) for value in union_terms)
        )
    else:
        union_log2 = math.inf
    payload = {
        "schema": "bch250-rowlocal-fanout-finite-interval-cover-v1",
        "status": (
            "BINARY64_COMPLETE_INTERVAL_DIAGNOSTIC"
            if not rejected_singletons
            else "BINARY64_INCOMPLETE_INTERVAL_DIAGNOSTIC"
        ),
        "parameters": {
            "outer_bits": B,
            "outer_dimension": args.outer_dimension,
            "outer_rows": args.outer_rows,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": distance / output_bits,
            "inner_epochs": epochs,
            "occupation_interval": [args.minimum_occupation, args.outer_rows],
            "output_weight_interval": [
                args.minimum_output_weight,
                args.maximum_output_weight,
            ],
            "spectrum": str(args.spectrum),
            "selection": str(args.selection),
        },
        "pointwise_comparison": {
            "maximum_log2_density_ratio": float(density_excess),
            "maximizing_weight": maximizing_weight,
        },
        "processed_intervals": processed,
        "accepted_intervals": len(accepted),
        "rejected_singletons": rejected_singletons,
        "union_log2_upper": union_log2,
        "union_margin_bits": -union_log2,
        "target_union_margin_bits": args.target_union_margin,
        "target_union_margin_met": (
            not rejected_singletons and -union_log2 >= args.target_union_margin
        ),
        "convexity_rule": (
            "For fixed p,z and Collatz data, the Q-dependent logarithm equals "
            "a linear function plus -(B-1)log C(L,Q); log C(L,Q) is concave. "
            "The finite log bound is therefore convex, so its maximum on an "
            "integer interval occurs at an endpoint."
        ),
        "intervals": accepted,
        "limitations": [
            "Witness optimization, eigendecomposition, and Collatz arithmetic use nearest binary64.",
            "Finite logarithms use 100-digit point arithmetic without directed rounding.",
            "An outward verifier must replay every accepted endpoint and the union sum.",
            "The probability space samples an independent row-local fanout wrapper in every outer row.",
        ],
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"processed_intervals={processed}")
    print(f"accepted_intervals={len(accepted)}")
    print(f"rejected_singletons={len(rejected_singletons)}")
    print(f"union_margin_bits={-union_log2:.9f}")
    print(f"target_union_margin_met={payload['target_union_margin_met']}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
