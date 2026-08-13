#!/usr/bin/env python3
"""Coupled outer/inner witness that covers every packet-weight profile.

Let ``t_j`` be the three-band outer profile fugacities and ``f_j`` the robust
inner profile fugacities.  The separate coefficient bounds for a profile
``a`` contain

    1 / (prod_j (t_j f_j)^a_j Q(a)),

where ``Q(a)=multinomial(M;a) prod_j C(8,j)^a_j``.  Coupling

    C(8,j) t_j f_j = 1

therefore makes the profile dependence exactly ``1/multinomial(M;a)``.
Every profile is no worse than a vertex, and summing all profiles costs only
``S_9(M)=sum_a 1/multinomial(M;a)``.  This probe uses the elementary rigorous
support-size bound

    S_9(M) <= sum_(r=1)^9 C(9,r)(M-r+1)/(M (r-1)!).

The linear Brascamp--Lieb outer inequality and robust local inner relaxation
are rigorous.  Binary64 moments and Perron iteration remain diagnostic until
the selected witness is outward hardened.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.special import logsumexp

from probe_packet8_repeated_value_state_transfer import (
    EXACT_SLICES,
    SPECTRUM,
    build_split_caps,
    load_ebch128_spectrum,
    load_exact,
    witness,
)
from probe_packet8_three_band_profile_enumerator import (
    BINOMIAL_LOG_PROBABILITIES,
    D,
    K,
    N,
    group_log_moments,
)
from probe_packet8_weight_profile_scalar import CLASSES
from probe_packet8_weight_profile_state_transfer import (
    INNER_BLOCKS,
    RobustProfileKernel,
    block_histograms,
)


M = N // 8


def profile_sum_upper(classes: int) -> float:
    return sum(
        math.comb(classes, support)
        * (M - support + 1)
        / (M * math.factorial(support - 1))
        for support in range(1, classes + 1)
    )


def outer_mgf_log2(log_t: np.ndarray, band1: float) -> float:
    band0 = 1.0 - band1
    group_logs = group_log_moments(log_t)
    band0_norm = band0 * logsumexp(
        BINOMIAL_LOG_PROBABILITIES + group_logs / band0
    )
    band1_norm = band1 * logsumexp(
        BINOMIAL_LOG_PROBABILITIES + group_logs / band1
    )
    return (
        K * math.log(2.0)
        + 256.0 * (42.0 * band0_norm + 86.0 * band1_norm)
    ) / math.log(2.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--log-t",
        required=True,
        help="comma-separated nine natural-log outer fugacities; first must be zero",
    )
    parser.add_argument("--band1-coefficient", type=float, required=True)
    parser.add_argument("--output-pole", type=float, required=True)
    parser.add_argument("--iterations", type=int, default=64)
    parser.add_argument(
        "--active-classes",
        default="0,1,2,3,4,5,6,7,8",
        help="comma-separated packet-weight classes covered by this face witness",
    )
    parser.add_argument(
        "--inactive-fugacity", type=float, default=0.0
    )
    args = parser.parse_args()

    log_t = np.asarray(
        [float(value) for value in args.log_t.split(",")], dtype=np.float64
    )
    if len(log_t) != 9 or abs(float(log_t[0])) > 1e-12:
        raise SystemExit("joint universal witness: log-t must have nine entries and start at zero")
    if not 0.5 <= args.band1_coefficient < 1.0:
        raise SystemExit("joint universal witness: invalid band coefficient")
    if not 0.0 < args.output_pole < 1.0:
        raise SystemExit("joint universal witness: invalid output pole")

    active = sorted({int(value) for value in args.active_classes.split(",")})
    if not active or any(value < 0 or value > 8 for value in active):
        raise SystemExit("joint universal witness: invalid active classes")
    if not 0.0 <= args.inactive_fugacity < 1.0:
        raise SystemExit("joint universal witness: invalid inactive fugacity")

    fugacities = np.full(9, args.inactive_fugacity, dtype=np.float64)
    all_coupled = np.exp(-log_t) / np.asarray(CLASSES, dtype=np.float64)
    fugacities[active] = all_coupled[active]
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    kernel = RobustProfileKernel(
        block_histograms(fugacities),
        split_caps,
        args.output_pole,
        float(np.max(fugacities)) ** 8,
    )
    eigenvalue, domination, values, worst_state = witness(
        kernel, args.iterations
    )
    inner_mgf_log2 = (
        math.log2(domination)
        + INNER_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )
    outer_log2 = outer_mgf_log2(log_t, args.band1_coefficient)
    vertex_log2 = (
        outer_log2
        + inner_mgf_log2
        - D * math.log2(args.output_pole)
    )
    profile_sum = profile_sum_upper(len(active))
    total_log2 = vertex_log2 + math.log2(profile_sum)

    print("packet-8 coupled universal packet-profile witness")
    print("active_classes=" + ",".join(map(str, active)))
    print("log_t=" + ",".join(f"{value:.12g}" for value in log_t))
    print("fugacities=" + ",".join(f"{value:.12g}" for value in fugacities))
    print(f"band1_coefficient={args.band1_coefficient:.12f}")
    print(f"output_pole={args.output_pole:.12f}")
    print(f"outer_mgf_log2={outer_log2:.12f}")
    print(f"inner_mgf_log2={inner_mgf_log2:.12f}")
    print(f"lambda_log2={math.log2(eigenvalue):.12f}")
    print(f"domination_log2={math.log2(domination):.12f}")
    print(f"worst_state={worst_state}")
    print(f"universal_vertex_log2={vertex_log2:.12f}")
    print(f"profile_sum_upper={profile_sum:.12f}")
    print(f"profile_sum_log2_upper={math.log2(profile_sum):.12f}")
    print(f"complete_profile_sum_log2_upper={total_log2:.12f}")
    print(f"diagnostic_le_2^-40={'PASS' if total_log2 <= -40.0 else 'NO'}")
    print("status=DIAGNOSTIC_BINARY64_UNIVERSAL_COUPLED_WITNESS")


if __name__ == "__main__":
    main()
