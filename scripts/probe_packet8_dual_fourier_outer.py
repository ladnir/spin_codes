#!/usr/bin/env python3
"""Probe the Fourier/MacWilliams outer-profile route.

For the radial 64-bit packet-profile factor ``f(x)=R_w(t)``, normalized
Fourier inversion gives

    sum_{c in C} product_g f(c_g)
      = |C| fhat(0)^B sum_{d in C^perp} product_g q_w(d_g),

where ``q_w=|fhat_w|/fhat_0`` after taking absolute values.  The zero dual
word is exactly the random-linear-code baseline.  This script computes the
radial Fourier ratios at the primal profile saddle point and reports both the
zero-word coefficient bound and the deliberately coarse exponent-three
Finner bound on the complete dual partition function.

The Fourier transform is evaluated with high-precision arithmetic because
the EBCH-relevant coefficients can arise after severe cancellation.  Results
remain diagnostic until the nonzero dual partition is bounded outward.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from decimal import Decimal, getcontext
from scipy.optimize import minimize_scalar

from probe_packet8_three_band_profile_enumerator import (
    BINOMIAL_LOG_PROBABILITIES,
    GROUPS,
    K,
    optimize_outer,
)
from probe_packet8_weight_profile_scalar import CLASSES, parse_profile
from punctured_ebch_outer import full_spectrum


def krawtchouk(length: int, output_weight: int, input_weight: int) -> int:
    return sum(
        (-1) ** intersection
        * math.comb(input_weight, intersection)
        * math.comb(length - input_weight, output_weight - intersection)
        for intersection in range(
            max(0, output_weight - (length - input_weight)),
            min(output_weight, input_weight) + 1,
        )
    )


def fourier_ratios(log_t: np.ndarray) -> tuple[list[Decimal], Decimal]:
    t = [Decimal(str(value)).exp() for value in log_t]
    packet_hat = []
    for dual_weight in range(9):
        value = sum(
            Decimal(krawtchouk(8, primal_weight, dual_weight)) * t[primal_weight]
            for primal_weight in range(9)
        ) / (1 << 8)
        packet_hat.append(value)

    polynomial = [Decimal(1)]
    one_packet = [Decimal(CLASSES[weight]) * packet_hat[weight] for weight in range(9)]
    for _ in range(8):
        result = [Decimal(0)] * (len(polynomial) + 8)
        for left, left_value in enumerate(polynomial):
            for right, right_value in enumerate(one_packet):
                result[left + right] += left_value * right_value
        polynomial = result
    radial_hat = [
        polynomial[weight] / Decimal(math.comb(64, weight)) for weight in range(65)
    ]
    if radial_hat[0] <= 0:
        raise SystemExit("dual Fourier: nonpositive constant coefficient")
    return [abs(value) / radial_hat[0] for value in radial_hat], radial_hat[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--precision", type=int, default=100)
    parser.add_argument("--saddle", choices=("zero", "primal"), default="zero")
    parser.add_argument(
        "--blend",
        type=float,
        help="override saddle by (1-blend)*primal + blend*zero in log fugacities",
    )
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"dual Fourier: {error}") from error
    getcontext().prec = args.precision

    coefficient_log, result, _components = optimize_outer(profile)
    primal_log_t = np.concatenate(([0.0], result.x))
    zero_log_t = np.array(
        [
            (
                math.log(count / classes)
                - math.log(profile[0] / CLASSES[0])
                if count
                else -40.0
            )
            for count, classes in zip(profile, CLASSES)
        ],
        dtype=np.float64,
    )
    log_t = primal_log_t if args.saddle == "primal" else zero_log_t
    if args.blend is not None:
        if not 0.0 <= args.blend <= 1.0:
            raise SystemExit("dual Fourier: blend must lie in [0,1]")
        log_t = (1.0 - args.blend) * primal_log_t + args.blend * zero_log_t
    ratios, fhat_zero = fourier_ratios(log_t)
    extraction = sum(
        Decimal(count) * Decimal(str(value))
        for count, value in zip(profile, log_t)
    )
    zero_coefficient = (
        Decimal(K) * Decimal(2).ln() + Decimal(GROUPS) * fhat_zero.ln() - extraction
    )

    def log_moment(power: int) -> Decimal:
        terms = [
            Decimal(math.comb(64, weight)).ln()
            - 64 * Decimal(2).ln()
            + power * ratio.ln()
            for weight, ratio in enumerate(ratios)
            if ratio
        ]
        maximum = max(terms)
        return maximum + sum((value - maximum).exp() for value in terms).ln()

    log_s3 = log_moment(3)
    dual_finner_log_z = Decimal(K) * Decimal(2).ln() + Decimal(GROUPS) / 3 * log_s3
    full_dual_coefficient = zero_coefficient + dual_finner_log_z

    spectrum = full_spectrum()
    ratio_logs = np.array(
        [-math.inf if not ratio else float(ratio.ln()) for ratio in ratios]
    )
    spectrum_weights = np.array(
        [weight for weight, count in enumerate(spectrum) if count], dtype=np.float64
    )
    spectrum_logs = np.log(
        np.array([spectrum[int(weight)] for weight in spectrum_weights], dtype=np.float64)
    )

    def total_weight_envelope(log_beta: float) -> float:
        log_alpha = float(
            np.max(ratio_logs - np.arange(65, dtype=np.float64) * log_beta)
        )
        local = float(
            np.logaddexp.reduce(spectrum_logs + spectrum_weights * log_beta)
        )
        return GROUPS * log_alpha + (K // 64) * local

    envelope_result = minimize_scalar(
        total_weight_envelope,
        bounds=(-20.0, 1.0),
        method="bounded",
        options={"xatol": 1e-13},
    )
    envelope_log_z = total_weight_envelope(float(envelope_result.x))

    print("packet-8 dual Fourier outer-profile probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"saddle={args.saddle}")
    if args.blend is not None:
        print(f"blend={args.blend:.12f}")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"primal_finner_coefficient_log2={coefficient_log/math.log(2.0):.12f}")
    log_two = Decimal(2).ln()
    print(f"fhat_zero_log2={fhat_zero.ln()/log_two:.18f}")
    print(f"zero_dual_coefficient_log2={zero_coefficient/log_two:.18f}")
    print(f"naive_dual_finner_log2_Z={dual_finner_log_z/log_two:.18f}")
    print(f"naive_dual_finner_coefficient_log2={full_dual_coefficient/log_two:.18f}")
    for power in (1, 2, 3):
        hypothetical = Decimal(K) * log_two + Decimal(GROUPS) / power * log_moment(power)
        print(f"degree{power}_moment_log2_Z={hypothetical/log_two:.18f}")
        print(
            f"degree{power}_moment_coefficient_log2="
            f"{zero_coefficient/log_two+hypothetical/log_two:.18f}"
        )
    print(f"total_weight_envelope_log2_Z={envelope_log_z/math.log(2.0):.18f}")
    print(
        "total_weight_envelope_coefficient_log2="
        f"{float(zero_coefficient/log_two)+envelope_log_z/math.log(2.0):.18f}"
    )
    print(f"total_weight_envelope_log_beta={envelope_result.x:.12f}")
    print("fourier_ratios_log2=" + ",".join(
        "-inf" if not ratio else f"{ratio.ln()/log_two:.12f}" for ratio in ratios
    ))
    dominant = max(
        range(65),
        key=lambda weight: Decimal(math.comb(64, weight)) * ratios[weight] ** 3,
    )
    print(f"dominant_dual_s3_group_weight={dominant}")
    print("status=DIAGNOSTIC_HIGH_PRECISION_DUAL_FOURIER")


if __name__ == "__main__":
    main()
