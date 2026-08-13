#!/usr/bin/env python3
"""Exact dense-orbit endpoint for the eight-bit packet permutation.

Fix the packet-weight histogram ``a_j`` after the independent random lane
orders, where ``c_j=C(8,j)`` packet values have weight j.  Conditional on the
weights, internal packet patterns are independent and uniform in their
weight classes.  If ``A_v`` is the count of the concrete packet value v, then

    E[prod_v A_v! | (a_j)]
      = prod_j a_j! C(a_j+c_j-1,c_j-1) / c_j^a_j.

A uniform permutation of the M packets has orbit size M!/prod_v A_v!.
The systematic recursive inner is a binary bijection: its emitted half is
systematic and the discarded state can be reconstructed causally from the
previous emitted block.  Hence at most Vol(N,d) inputs in any packet orbit
can map to output weight at most d.

This script checks the balanced dense endpoint

    a_j = (M/256) C(8,j),

which is the exact packet-weight profile of the uniform eight-bit law.  It
unions over all 2^K messages and deliberately does not credit the graph
codimension.  The result is an exact-integer endpoint certificate, not yet a
complete profile ledger: profiles below the resulting orbit-entropy threshold
must be covered by the sparse/intermediate packet transfer branch.
"""

from __future__ import annotations

import math


N = 1 << 21
K = 1 << 20
PACKET_BITS = 8
M = N // PACKET_BITS
D = 9 * N // 100
TARGET_BITS = 40


def log2_int(value: int) -> float:
    shift = max(0, value.bit_length() - 53)
    return math.log2(value >> shift) + shift


def main() -> None:
    classes = [math.comb(PACKET_BITS, weight) for weight in range(PACKET_BITS + 1)]
    if M % (1 << PACKET_BITS):
        raise SystemExit("packet8 dense orbit: packet count is not divisible by 256")
    scale = M // (1 << PACKET_BITS)
    profile = [scale * count for count in classes]
    if sum(profile) != M or sum(weight * count for weight, count in enumerate(profile)) != N // 2:
        raise SystemExit("packet8 dense orbit: balanced profile mismatch")

    # For d<N/2, adjacent binomial coefficients below d are dominated by a
    # geometric series with ratio at most d/(N-d+1).
    volume_numerator = math.comb(N, D) * (N - D + 1)
    volume_denominator = N - 2 * D + 1

    collision_numerator = 1
    collision_denominator = math.factorial(M)
    for count, classes_at_weight in zip(profile, classes):
        collision_numerator *= math.factorial(count)
        collision_numerator *= math.comb(
            count + classes_at_weight - 1, classes_at_weight - 1
        )
        collision_denominator *= classes_at_weight**count

    # 2^K * Vol(N,d) * E[prod A_v!]/M! <= 2^-40.
    left = (
        (1 << (K + TARGET_BITS))
        * volume_numerator
        * collision_numerator
    )
    right = volume_denominator * collision_denominator
    if left > right:
        raise SystemExit("packet8 dense orbit: balanced endpoint exceeds 2^-40")

    log_volume_upper = log2_int(volume_numerator) - log2_int(volume_denominator)
    log_collision = log2_int(collision_numerator) - log2_int(collision_denominator)
    union_log2_upper = K + log_volume_upper + log_collision
    margin = -TARGET_BITS - union_log2_upper

    print("eight-bit packet dense-orbit endpoint certificate")
    print(f"N={N} K={K} d={D} packets={M}")
    print("inner_local_code=[128,64,22]")
    print("outer_local_code=[128,64,22]")
    print("balanced_packet_weight_profile=" + ",".join(map(str, profile)))
    print(f"hamming_ball_log2_upper={log_volume_upper:.12f}")
    print(f"expected_inverse_orbit_log2={log_collision:.12f}")
    print(f"all_messages_union_log2_upper={union_log2_upper:.12f}")
    print(f"margin_beyond_40_bits={margin:.12f}")
    print("graph_codimension_credited=0")
    print("balanced_dense_profile_le_2^-40=PASS")
    print("status=EXACT_INTEGER_PACKET8_DENSE_ENDPOINT")


if __name__ == "__main__":
    main()
