#!/usr/bin/env python3
"""Exact high-profile tests for the packet-8 orbit and bijection branches.

For an exact packet-weight profile ``a_j`` this script checks two independent
all-message bounds.

Orbit branch
    Average the inverse concrete-value orbit over the independent uniform
    internal patterns of each weight class and use the geometric Hamming-ball
    volume bound.

Bijection-moment branch
    Conditional on the weight profile, packet order and internal patterns are
    uniform over

        Q(a) = M! / product_j a_j! * product_j C(8,j)^a_j

    binary inputs.  Summing the ``z^weight`` moment over *all* binary inputs is
    exactly ``(1+z)^N`` because every frozen recursive inner map is a binary
    bijection.  Positive coefficient extraction at the nearby exact pole
    ``z=1/10`` therefore gives an integer-comparable full-``2^K`` union bound.
    The exact optimizing pole saves only about 15.2 diagnostic bits here but
    makes the direct integer comparison much more expensive.

The two PASS tests use integers only.  Logarithms are diagnostics.
"""

from __future__ import annotations

import argparse
import math


N = 1 << 21
K = 1 << 20
M = N // 8
D = 9 * N // 100
TARGET_BITS = 40
CLASSES = tuple(math.comb(8, weight) for weight in range(9))


def parse_profile(specification: str) -> list[int]:
    profile = [int(value) for value in specification.split(",")]
    if len(profile) != 9 or any(value < 0 for value in profile):
        raise ValueError("profile must contain nine nonnegative counts")
    if sum(profile) != M:
        raise ValueError(f"profile counts must sum to {M}")
    return profile


def log2_int(value: int) -> float:
    shift = max(0, value.bit_length() - 53)
    return math.log2(value >> shift) + shift


def products(profile: list[int]):
    profile_factorials = 1
    internal_patterns = 1
    stars_bars = 1
    for count, classes in zip(profile, CLASSES):
        profile_factorials *= math.factorial(count)
        internal_patterns *= classes**count
        stars_bars *= math.comb(count + classes - 1, classes - 1)
    return profile_factorials, internal_patterns, stars_bars


def check_profile(profile: list[int]) -> None:
    profile_factorials, internal_patterns, stars_bars = products(profile)
    factorial_m = math.factorial(M)

    # Vol(N,D) <= C(N,D) (N-D+1)/(N-2D+1).
    volume_numerator = math.comb(N, D) * (N - D + 1)
    volume_denominator = N - 2 * D + 1

    # E[prod_v A_v!] = prod_j a_j! stars_bars_j / c_j^a_j.
    orbit_left = (
        (1 << (K + TARGET_BITS))
        * volume_numerator
        * profile_factorials
        * stars_bars
    )
    orbit_right = volume_denominator * factorial_m * internal_patterns
    orbit_pass = orbit_left <= orbit_right

    # At z=1/10, (1+z)^N z^-D = 11^N / 10^(N-D).
    moment_left = (
        (1 << (K + TARGET_BITS))
        * 11**N
        * profile_factorials
    )
    moment_right = 10 ** (N - D) * factorial_m * internal_patterns
    moment_pass = moment_left <= moment_right

    normalization_log2 = (
        log2_int(factorial_m)
        - log2_int(profile_factorials)
        + log2_int(internal_patterns)
    )
    orbit_bits = normalization_log2 - log2_int(stars_bars)
    volume_log2 = log2_int(volume_numerator) - log2_int(volume_denominator)
    orbit_union = K + volume_log2 - orbit_bits
    moment_numerator_log2 = N * math.log2(11) - (N - D) * math.log2(10)
    moment_union = K + moment_numerator_log2 - normalization_log2

    print("packet-8 exact high-profile branch certificate")
    print(f"profile={','.join(map(str, profile))}")
    print(f"profile_normalization_log2={normalization_log2:.12f}")
    print(f"stars_bars_correction_log2={log2_int(stars_bars):.12f}")
    print(f"effective_inverse_orbit_bits={orbit_bits:.12f}")
    print(f"orbit_all_messages_union_log2_upper={orbit_union:.12f}")
    print(f"orbit_le_2^-40={'PASS' if orbit_pass else 'NO'}")
    print(f"bijection_pole=1/10")
    print(f"bijection_all_messages_union_log2_upper={moment_union:.12f}")
    print(f"bijection_moment_le_2^-40={'PASS' if moment_pass else 'NO'}")
    print(
        "high_profile_covered="
        + ("PASS" if orbit_pass or moment_pass else "NO")
    )
    print("status=EXACT_INTEGER_PACKET8_HIGH_PROFILE_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"packet8 high profile: {error}") from error
    check_profile(profile)


if __name__ == "__main__":
    main()
