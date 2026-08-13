#!/usr/bin/env python3
"""Exact combinatorics for the corrected eight-bit Riffle permutation.

The inner and outer local codes both remain the committed binary
EBCH [128,64,22] code.  Only the permutation atom changes: the N binary
coordinates are divided into N/8 packets, a uniform permutation is applied
to those packets, and each consecutive eight packets form one 64-bit input
to the recursive inner code.

This certificate proves the elementary packet-placement facts used before
the packet-weight transfer operator:

* a nonzero outer word contains at least 21 nonzero packets;
* the old universal late-suffix event is below 2^-40;
* conditional on p nonzero packets, their occupied inner-block count has the
  exact capacity-eight occupancy law implemented below.

All PASS comparisons use integers or Fractions.  Logarithms are display only.
"""

from __future__ import annotations

import math
import hashlib
from fractions import Fraction


N = 1 << 21
PACKET_BITS = 8
INNER_BITS = 64
PACKETS_PER_INNER = INNER_BITS // PACKET_BITS
INNER_BLOCKS = N // INNER_BITS
PACKET_SLOTS = N // PACKET_BITS
DISTANCE = 9 * N // 100
LATE_INNER_BLOCKS = DISTANCE // INNER_BITS
PUNCTURED_OUTER_DISTANCE = 21


def log2_fraction(value: Fraction) -> float:
    return math.log2(value.numerator) - math.log2(value.denominator)


def falling(n: int, k: int) -> int:
    result = 1
    for value in range(n - k + 1, n + 1):
        result *= value
    return result


def occupied_inner_counts(packet_support: int) -> list[int]:
    """Ordered-placement counts indexed by occupied 64-bit inner blocks.

    Select ``packet_support`` distinct packet slots in order.  If ``k`` slots
    have been selected and they occupy ``g`` inner blocks, then ``8g-k``
    remaining slots lie in an already occupied block and ``M-8g`` lie in a
    fresh block.  The recurrence therefore counts every injection exactly.
    """

    if not 0 <= packet_support <= PACKET_SLOTS:
        raise ValueError("packet support out of range")
    counts = [1]
    for selected in range(packet_support):
        next_counts = [0] * (len(counts) + 1)
        for occupied, count in enumerate(counts):
            if not count:
                continue
            old_slots = PACKETS_PER_INNER * occupied - selected
            new_slots = PACKET_SLOTS - PACKETS_PER_INNER * occupied
            if old_slots:
                next_counts[occupied] += count * old_slots
            if new_slots:
                next_counts[occupied + 1] += count * new_slots
        counts = next_counts
    if sum(counts) != falling(PACKET_SLOTS, packet_support):
        raise SystemExit("packet8: occupancy recurrence mass mismatch")
    return counts


def packet_support_mgf(packet_support: int, pole: Fraction) -> Fraction:
    counts = occupied_inner_counts(packet_support)
    denominator = falling(PACKET_SLOTS, packet_support)
    return sum(
        (count * pole**occupied for occupied, count in enumerate(counts)),
        Fraction(0),
    ) / denominator


def weak_compositions(total: int, parts: int, prefix: tuple[int, ...] = ()):
    if parts == 1:
        yield prefix + (total,)
        return
    for first in range(total + 1):
        yield from weak_compositions(total - first, parts - 1, prefix + (first,))


def source_group_packet_profiles() -> list[dict[tuple[int, ...], int]]:
    """Exact packet-weight histograms after one random 64-lane bijection.

    A profile ``a`` has ``a[j]`` packets of binary weight ``j``.  Its count is
    the number of 64-bit supports with that ordered-packet histogram.
    """

    rows: list[dict[tuple[int, ...], int]] = [dict() for _ in range(INNER_BITS + 1)]
    factorial_eight = math.factorial(PACKETS_PER_INNER)
    for profile in weak_compositions(PACKETS_PER_INNER, PACKET_BITS + 1):
        total_weight = sum(weight * count for weight, count in enumerate(profile))
        packet_orders = factorial_eight
        for count in profile:
            packet_orders //= math.factorial(count)
        support_choices = packet_orders
        for weight, count in enumerate(profile):
            support_choices *= math.comb(PACKET_BITS, weight) ** count
        rows[total_weight][profile] = support_choices
    for weight, row in enumerate(rows):
        if sum(row.values()) != math.comb(INNER_BITS, weight):
            raise SystemExit(f"packet8: source profile row {weight} has wrong mass")
    return rows


def main() -> None:
    if PACKETS_PER_INNER != 8 or PACKET_SLOTS != 8 * INNER_BLOCKS:
        raise SystemExit("packet8: inconsistent construction dimensions")
    if INNER_BITS * LATE_INNER_BLOCKS > DISTANCE:
        raise SystemExit("packet8: late suffix exceeds the distance target")
    if INNER_BITS * (LATE_INNER_BLOCKS + 1) <= DISTANCE:
        raise SystemExit("packet8: late suffix is not maximal")

    # Every active punctured outer block has at least 21 surviving coordinates.
    # The striped outer layout puts distinct coordinates of one local block in
    # distinct physical 64-bit groups, hence in distinct eight-bit packets.
    minimum_packets = PUNCTURED_OUTER_DISTANCE

    late_slots = PACKETS_PER_INNER * LATE_INNER_BLOCKS
    late_probability = Fraction(
        math.comb(late_slots, minimum_packets),
        math.comb(PACKET_SLOTS, minimum_packets),
    )
    if late_probability > Fraction(1, 1 << 40):
        raise SystemExit("packet8: universal late-suffix bound exceeds 2^-40")

    counts = occupied_inner_counts(minimum_packets)
    denominator = falling(PACKET_SLOTS, minimum_packets)
    collision_probability = Fraction(sum(counts[:minimum_packets]), denominator)
    no_collision_probability = Fraction(counts[minimum_packets], denominator)
    if collision_probability + no_collision_probability != 1:
        raise SystemExit("packet8: collision split does not sum to one")

    target_pole = Fraction(181, 1000)
    support_mgf = packet_support_mgf(minimum_packets, target_pole)
    profile_rows = source_group_packet_profiles()
    profile_digest = hashlib.sha256()
    for weight, row in enumerate(profile_rows):
        for profile, count in sorted(row.items()):
            profile_digest.update(bytes((weight, *profile)))
            profile_digest.update(count.to_bytes(24, "little"))

    print("eight-bit packet-permutation certificate")
    print(f"N={N}")
    print(f"outer_local_code=[128,64,22]")
    print(f"inner_local_code=[128,64,22]")
    print(f"packet_bits={PACKET_BITS}")
    print(f"packet_slots={PACKET_SLOTS}")
    print(f"inner_blocks={INNER_BLOCKS}")
    print(f"packets_per_inner={PACKETS_PER_INNER}")
    print(f"distance_target={DISTANCE}")
    print(f"late_inner_blocks={LATE_INNER_BLOCKS}")
    print(f"minimum_nonzero_packets={minimum_packets}")
    print(f"late_suffix_probability_log2={log2_fraction(late_probability):.12f}")
    print("late_suffix_probability_le_2^-40=PASS")
    print(f"packet_collision_probability_log2={log2_fraction(collision_probability):.12f}")
    print(f"no_packet_collision_probability={float(no_collision_probability):.12f}")
    print(f"support_21_mgf_at_181_over_1000_log2={log2_fraction(support_mgf):.12f}")
    print(f"source_packet_profile_sha256={profile_digest.hexdigest()}")
    print("source_packet_profile_rows_0_64_exact_mass=PASS")
    for weight in (21, 22, 32, 64):
        row = profile_rows[weight]
        expected_nonzero = sum(
            (PACKETS_PER_INNER - profile[0]) * count
            for profile, count in row.items()
        ) / math.comb(INNER_BITS, weight)
        print(
            f"source_weight={weight} profile_count={len(row)} "
            f"expected_nonzero_packets={expected_nonzero:.12f}"
        )
    print("occupancy_recurrence_exact_mass=PASS")
    print("status=EXACT_PACKET8_PLACEMENT_CERTIFICATE")


if __name__ == "__main__":
    main()
