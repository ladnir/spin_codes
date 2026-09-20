#!/usr/bin/env python3
"""Certify the Random SPIN exponent endpoints and sparse schedule margins.

The continuum argument in RANDOM_SPIN_PROOF_AUDIT.md reduces the sign check
to three endpoint values and the derivative at zero.  The sparse proof reduces
to two additional strict schedule inequalities.  This script evaluates all of
those expressions with mpmath interval arithmetic.
"""

from __future__ import annotations

import argparse
import json

import mpmath as mp


def log2(value: mp.ctx_iv.ivmpf) -> mp.ctx_iv.ivmpf:
    return mp.iv.log(value) / mp.iv.log(2)


def entropy(value: mp.ctx_iv.ivmpf) -> mp.ctx_iv.ivmpf:
    return -value * log2(value) - (1 - value) * log2(1 - value)


def endpoints(value: mp.ctx_iv.ivmpf) -> list[float]:
    return [float(value.a), float(value.b)]


def interval_argument(value: str) -> mp.ctx_iv.ivmpf:
    if "/" not in value:
        return mp.iv.mpf(value)
    numerator, denominator = value.split("/", maxsplit=1)
    return mp.iv.mpf(numerator) / mp.iv.mpf(denominator)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta", default="5501/50000")
    parser.add_argument("--sparse-cutoff", default="1/5000")
    parser.add_argument("--interior-parameter", default="101/100")
    parser.add_argument("--memory-gamma", default="51/50")
    parser.add_argument("--block-constant", default="17")
    args = parser.parse_args()

    delta = interval_argument(args.delta)
    eta_sparse = interval_argument(args.sparse_cutoff)
    saddle_scale = interval_argument(args.interior_parameter)
    memory_gamma = interval_argument(args.memory_gamma)
    block_constant = interval_argument(args.block_constant)
    one = mp.iv.mpf(1)
    two = mp.iv.mpf(2)
    eta_zero = one - one / mp.iv.sqrt(two)
    eta_cross = one - one / (two * (one - delta))

    root = mp.iv.sqrt(eta_zero**2 + delta**2)
    saddle_a = eta_zero / (root + delta)
    saddle_s = (one - saddle_a) / (one + saddle_a)
    outer_at_zero = eta_zero * log2(one + mp.iv.sqrt(two))
    inner_at_zero = (
        -eta_zero * log2(saddle_a)
        - delta * log2(saddle_s)
        - entropy(eta_zero)
    )
    junction_zero = outer_at_zero + inner_at_zero
    junction_cross = entropy(eta_cross) + entropy(delta) - mp.iv.mpf("1.5")
    bulk_maximum = entropy(delta) - mp.iv.mpf("0.5")
    zero_slope = log2(two * delta * (one + mp.iv.sqrt(two)))
    sparse_a = eta_sparse / (delta * (saddle_scale + one))
    sparse_epsilon = (
        saddle_scale**2
        * sparse_a
        / ((saddle_scale + one) * (one - saddle_scale * sparse_a))
        + eta_sparse / (two * (one - eta_sparse))
    )
    sparse_rho = (
        delta * (saddle_scale + one) * mp.iv.exp(sparse_epsilon)
    )
    outer_decay = log2(mp.iv.sqrt(two) / (one + sparse_rho))
    polynomial_exponent = two * memory_gamma + mp.iv.mpf("1.5")
    block_schedule_margin = (
        block_constant * outer_decay - polynomial_exponent
    )
    sparse_b = (one + sparse_a) / two
    memory_schedule_margin = memory_gamma * log2(one / sparse_b) - one
    top_outer_exponent = entropy(mp.iv.mpf("0.9")) - mp.iv.mpf("0.5")

    checks = {
        "right_derivative_at_zero": endpoints(zero_slope),
        "junction_at_outer_branch": endpoints(junction_zero),
        "junction_at_inner_branch": endpoints(junction_cross),
        "bulk_maximum": endpoints(bulk_maximum),
    }
    if any(interval[1] >= 0.0 for interval in checks.values()):
        raise AssertionError("an exponent endpoint is not certified negative")
    sparse_checks = {
        "sparse_rho": endpoints(sparse_rho),
        "outer_contraction_boundary": endpoints(
            mp.iv.sqrt(two) - one
        ),
        "outer_decay_log2": endpoints(outer_decay),
        "block_schedule_margin": endpoints(block_schedule_margin),
        "memory_schedule_margin": endpoints(memory_schedule_margin),
        "top_outer_exponent": endpoints(top_outer_exponent),
    }
    if sparse_checks["sparse_rho"][1] >= sparse_checks[
        "outer_contraction_boundary"
    ][0]:
        raise AssertionError("the sparse contraction point is too large")
    if sparse_checks["block_schedule_margin"][0] <= 0.0:
        raise AssertionError("the block schedule does not close")
    if sparse_checks["memory_schedule_margin"][0] <= 0.0:
        raise AssertionError("the memory schedule does not close")
    if sparse_checks["top_outer_exponent"][1] >= 0.0:
        raise AssertionError("the top outer range does not close")

    payload = {
        "schema": "random-spin-asymptotic-exponent-interval-v1",
        "status": "CERTIFIED_ASYMPTOTIC_CONSTANTS",
        "arithmetic": {
            "library": "mpmath.iv",
            "mpmath_version": mp.__version__,
            "input_delta_decimal_interval": str(delta),
        },
        "parameters": {
            "delta": args.delta,
            "outer_branch_junction_eta": endpoints(eta_zero),
            "inner_branch_junction_eta": endpoints(eta_cross),
        },
        "certified_log2_intervals": checks,
        "certified_sparse_intervals": sparse_checks,
        "sparse_parameters": {
            "maximum_relative_input_weight": args.sparse_cutoff,
            "interior_saddle_parameter": args.interior_parameter,
            "memory_log2_constant": args.memory_gamma,
            "minimum_outer_block_log2_constant": args.block_constant,
            "polynomial_exponent": endpoints(polynomial_exponent),
            "sparse_a_at_endpoint": endpoints(sparse_a),
            "sparse_exponent_correction": endpoints(sparse_epsilon),
        },
        "analytic_reduction": {
            "first_region": (
                "The combined exponent is convex, equals zero at eta=0, "
                "and is negative at the outer-branch junction."
            ),
            "second_region": (
                "The combined exponent is increasing, so its maximum is at "
                "the inner-branch junction."
            ),
            "third_region": (
                "Binary entropy is at most one, so the bulk maximum is "
                "H_2(delta)-1/2."
            ),
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
