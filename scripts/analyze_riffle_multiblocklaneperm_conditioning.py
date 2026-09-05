#!/usr/bin/env python3
"""Conditioning surcharge for MultiBlockLanePerm-q FieldCheckpoint.

In every outer-coordinate region, each q-block group receives an independent
uniform lane permutation.  If a group has k active blocks, its labeled lane
assignment is the iid assignment conditioned on distinct lanes.  The
conditioning event has probability p_k=(q)_k/q^k.

For a fixed group-occupancy profile (k_u), every nonnegative 256-region moment
is at most prod_u p_{k_u}^-256 times the iid-lane moment.  This script sums that
factor over every active-block subset by a log-domain generating function and
applies it to an existing iid-lane 9% receipt.

The combinatorial conditioning inequality is exact.  The imported iid moment
and the coefficient calculation use ordinary floating-point arithmetic, so the
result is a proof diagnostic rather than an outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_striped_random_outer import LOG2, log_choose


DEFAULT_IID_RECEIPT = Path(
    "constructions/riffle_multiblockstripe_fieldcheckpoint/"
    "receipts/q32_independent_lane_low64_delta09.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_multiblocklaneperm_fieldcheckpoint/"
    "receipts/q32_conditioning_low64_delta09.json"
)


def log_falling_factorial(q: int, k: int) -> float:
    return math.lgamma(q + 1) - math.lgamma(q - k + 1)


def group_log_weights(q: int, regions: int) -> np.ndarray:
    """Return log(C(q,k) p_k^-regions) for k=0,...,q."""
    weights = np.empty(q + 1, dtype=np.float64)
    for k in range(q + 1):
        log_p = log_falling_factorial(q, k) - k * math.log(q)
        weights[k] = log_choose(q, k) - regions * log_p
    return weights


def powered_polynomial_logs(
    group_weights: np.ndarray, groups: int, maximum_degree: int
) -> np.ndarray:
    """Return logs of coefficients through maximum_degree in f(x)^groups."""
    current = np.full(maximum_degree + 1, -math.inf, dtype=np.float64)
    current[0] = 0.0
    q = len(group_weights) - 1
    for _ in range(groups):
        updated = np.full(maximum_degree + 1, -math.inf, dtype=np.float64)
        for occupied in range(maximum_degree + 1):
            maximum_here = min(q, occupied)
            terms = np.asarray(
                [
                    current[occupied - here] + group_weights[here]
                    for here in range(maximum_here + 1)
                ]
            )
            updated[occupied] = float(logsumexp(terms))
        current = updated
    return current


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    payload = json.loads(args.iid_receipt.read_text(encoding="utf-8"))
    iid_result = payload["results"][0]
    iid_rows = iid_result["independent_per_coordinate_lane_model"][
        "occupation_rows"
    ]
    maximum = min(args.maximum_active_blocks, len(iid_rows))
    q = args.q
    if int(iid_result["q"]) != q:
        raise ValueError("the iid receipt q does not match --q")
    outer_blocks = int(payload["parameters"]["outer_blocks"])
    if outer_blocks % q:
        raise ValueError("q must divide the number of outer blocks")
    groups = outer_blocks // q

    weights = group_log_weights(q, args.outer_coordinates)
    coefficient_logs = powered_polynomial_logs(weights, groups, maximum)

    rows = []
    adjusted_contributions = []
    for occupation in range(1, maximum + 1):
        iid_row = iid_rows[occupation - 1]
        if int(iid_row["active_regular_outer_blocks"]) != occupation:
            raise ValueError("iid receipt occupations are not contiguous")
        log_gamma = coefficient_logs[occupation] - log_choose(
            outer_blocks, occupation
        )
        iid_pointwise = float(iid_row["pointwise_log2_upper"])
        adjusted_pointwise = iid_pointwise + log_gamma / LOG2
        adjusted_contributions.append(adjusted_pointwise * LOG2)
        rows.append(
            {
                "active_regular_outer_blocks": occupation,
                "iid_pointwise_log2_upper": iid_pointwise,
                "conditioning_surcharge_log2": log_gamma / LOG2,
                "conditioned_pointwise_log2_upper": adjusted_pointwise,
                "conditioned_pointwise_margin_bits": -adjusted_pointwise,
            }
        )

    aggregate_log = float(logsumexp(np.asarray(adjusted_contributions)))
    dominant = sorted(
        rows,
        key=lambda row: float(row["conditioned_pointwise_log2_upper"]),
        reverse=True,
    )[:10]
    return {
        "schema": "riffle-multiblocklaneperm-conditioning-v1",
        "candidate": f"Riffle MultiBlockLanePerm-{q} FieldCheckpoint",
        "parameters": {
            "q": q,
            "outer_blocks": outer_blocks,
            "groups": groups,
            "outer_coordinates": args.outer_coordinates,
            "maximum_active_blocks": maximum,
            "relative_distance": payload["parameters"]["relative_distance"],
            "distance": payload["parameters"]["distance"],
        },
        "conditioning": {
            "per_group_probability": "p_k=(q)_k/q^k",
            "profile_factor": "product_u p_{k_u}^(-outer_coordinates)",
            "group_polynomial": (
                "sum_k C(q,k) p_k^(-outer_coordinates) x^k"
            ),
        },
        "aggregate_low_occupation_log2_upper": aggregate_log / LOG2,
        "aggregate_low_occupation_lambda_bits_lower_float": -aggregate_log
        / LOG2,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "scope": (
            "Rigorous conditioning reduction applied to the floating-point "
            "iid-lane regular-spectrum envelope through the requested "
            "occupation. All-one words, higher occupations, the modeled "
            "outer spectrum, and outward rounding remain outside this result."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, default=32)
    parser.add_argument("--outer-coordinates", type=int, default=256)
    parser.add_argument("--maximum-active-blocks", type=int, default=64)
    parser.add_argument("--iid-receipt", type=Path, default=DEFAULT_IID_RECEIPT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "aggregate_lambda_bits,"
        f"{result['aggregate_low_occupation_lambda_bits_lower_float']:.9f}"
    )
    for row in result["dominant_rows"][:5]:
        print(
            f"occupation,{row['active_regular_outer_blocks']},"
            f"surcharge,{row['conditioning_surcharge_log2']:.9f},"
            f"margin,{row['conditioned_pointwise_margin_bits']:.9f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
