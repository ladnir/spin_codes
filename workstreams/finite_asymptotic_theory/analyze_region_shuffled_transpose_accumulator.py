#!/usr/bin/env python3
"""Fixed-active-block saddle for transpose plus region shuffles.

This diagnostic evaluates the continuum two-state transfer described in
REGION_SHUFFLED_TRANSPOSE_ACCUMULATOR.md.  It proves neither uniformity in the
number of active blocks nor finite-length rounding.
"""

from __future__ import annotations

import argparse
import json
import math

import scipy.optimize
import scipy.special


def beta_laplace(shape_left: int, shape_total: int, theta: float) -> float:
    return float(scipy.special.hyp1f1(shape_left, shape_total, -theta))


def region_matrix(q: int, theta: float) -> tuple[float, float, float, float]:
    entries = [0.0, 0.0, 0.0, 0.0]
    denominator = float(2**q)
    for active_bits in range(q + 1):
        probability = math.comb(q, active_bits) / denominator
        if active_bits == 0:
            forward = 1.0
            complement = math.exp(-theta)
        else:
            occupied_spacings = (active_bits + 1) // 2
            empty_spacings = active_bits + 1 - occupied_spacings
            forward = beta_laplace(
                occupied_spacings, active_bits + 1, theta
            )
            complement = beta_laplace(
                empty_spacings, active_bits + 1, theta
            )
        if active_bits % 2 == 0:
            entries[0] += probability * forward
            entries[3] += probability * complement
        else:
            entries[1] += probability * forward
            entries[2] += probability * complement
    return tuple(entries)  # type: ignore[return-value]


def spectral_radius(matrix: tuple[float, float, float, float]) -> float:
    top_left, top_right, bottom_left, bottom_right = matrix
    discriminant = (
        (top_left - bottom_right) ** 2
        + 4.0 * top_right * bottom_left
    )
    return (
        top_left + bottom_right + math.sqrt(max(0.0, discriminant))
    ) / 2.0


def support_exponent(
    q: int,
    rate: float,
    delta: float,
    theta: float,
) -> float:
    radius = spectral_radius(region_matrix(q, theta))
    return rate * q + math.log2(radius) + theta * delta / math.log(2.0)


def optimize_q(
    q: int,
    rate: float,
    delta: float,
    theta_max: float,
) -> dict[str, float | int | None]:
    result = scipy.optimize.minimize_scalar(
        lambda theta: support_exponent(q, rate, delta, theta),
        bounds=(0.0, theta_max),
        method="bounded",
        options={"xatol": 1e-12},
    )
    exponent = float(result.fun)
    decay = -exponent
    required_constant = q / decay if decay > 0.0 else None
    return {
        "active_blocks": q,
        "optimizing_theta": float(result.x),
        "per_region_block_log2_exponent": exponent,
        "decay_kappa": decay,
        "required_log2_block_constant": required_constant,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rate", type=float, default=0.5)
    parser.add_argument("--delta", type=float, default=0.02)
    parser.add_argument("--maximum-active-blocks", type=int, default=32)
    parser.add_argument("--theta-maximum", type=float, default=500.0)
    args = parser.parse_args()

    if not 0.0 < args.rate < 1.0:
        raise SystemExit("rate must lie strictly between zero and one")
    if not 0.0 < args.delta < 0.5:
        raise SystemExit("delta must lie strictly between zero and one half")
    if args.maximum_active_blocks <= 0:
        raise SystemExit("maximum active blocks must be positive")

    rows = [
        optimize_q(q, args.rate, args.delta, args.theta_maximum)
        for q in range(1, args.maximum_active_blocks + 1)
    ]
    finite_constants = [
        float(row["required_log2_block_constant"])
        for row in rows
        if row["required_log2_block_constant"] is not None
    ]
    payload = {
        "schema": "region-shuffled-transpose-accumulator-fixed-q-v1",
        "status": "DIAGNOSTIC_FIXED_Q_CONTINUUM_SADDLE",
        "parameters": {
            "rate": args.rate,
            "relative_distance": args.delta,
            "maximum_active_blocks": args.maximum_active_blocks,
            "theta_maximum": args.theta_maximum,
        },
        "summary": {
            "largest_required_constant": (
                max(finite_constants) if finite_constants else None
            ),
            "all_checked_classes_have_positive_decay": len(finite_constants)
            == len(rows),
        },
        "fixed_active_block_rows": rows,
        "limitations": [
            "No uniform proof for q growing with the number of outer blocks.",
            "No outward-rounded arithmetic.",
            "The active block rows are relaxed from uniform nonzero to iid bits; the conditioning factor is asymptotically one for fixed q.",
        ],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
