"""Exact packet-occupancy census and exploratory Turbo one-block screen.

The accumulator/packet counts are exact integers. The BAA fold and s22
continuum moment use binary64 for discovery only; this is not a distance
certificate. No encoder benchmark is run.
"""

from collections import Counter
from math import ceil, comb, exp, expm1, log, log1p
from scipy.optimize import minimize_scalar
import argparse
import numpy as np


GOLAY = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
STATE_NONZERO = (1 << 22) - 1
LIVE_BIT_DENSITY = (1 << 21) / STATE_NONZERO
# Enumerated from certified22/productionMap.h oldRows; verified separately by
# verify_pair_return_map.py --spectrum.
STATE_OUTPUT_WEIGHTS = {
    16: 2604, 24: 291648, 28: 888832, 32: 1828134,
    36: 888832, 40: 291648, 48: 2604, 64: 1,
}
assert sum(STATE_OUTPUT_WEIGHTS.values()) == STATE_NONZERO
assert sum(w * n for w, n in STATE_OUTPUT_WEIGHTS.items()) == 64 * (1 << 21)


def packet_patterns() -> tuple[list[tuple[int, int, int, int]], ...]:
    result = []
    for previous in range(2):
        multiplicities: Counter[tuple[int, int, int]] = Counter()
        for word in range(16):
            bits = [(word >> i) & 1 for i in range(4)]
            transitions = (bits[0] ^ previous) + sum(
                bits[i] ^ bits[i - 1] for i in range(1, 4)
            )
            multiplicities[(bits[3], transitions, bool(word))] += 1
        result.append([(end, transitions, int(active), count)
                       for (end, transitions, active), count in multiplicities.items()])
    return result[0], result[1]


def packet_patterns_with_deficit() -> tuple[list[tuple[int, int, int, int, int]], ...]:
    """Packet transition, activity, and deficit from four one bits."""
    result = []
    for previous in range(2):
        multiplicities: Counter[tuple[int, int, int, int]] = Counter()
        for word in range(16):
            bits = [(word >> i) & 1 for i in range(4)]
            transitions = (bits[0] ^ previous) + sum(
                bits[i] ^ bits[i - 1] for i in range(1, 4)
            )
            active = int(word != 0)
            deficit = 4 * active - word.bit_count()
            multiplicities[(bits[3], transitions, active, deficit)] += 1
        result.append([(end, transitions, active, deficit, count)
                       for (end, transitions, active, deficit), count
                       in multiplicities.items()])
    return result[0], result[1]


def accumulator_packet_counts(regions: int) -> list[list[int]]:
    """D[c][v]: outputs with c transitions and v nonzero four-bit packets."""
    width = 4 * regions
    patterns = packet_patterns()
    states = [{(0, 0): 1}, {}]
    for _ in range(regions):
        following: list[dict[tuple[int, int], int]] = [{}, {}]
        for previous in range(2):
            for (old_c, old_v), words in states[previous].items():
                for end, added_c, active, multiple in patterns[previous]:
                    key = old_c + added_c, old_v + active
                    dest = following[end]
                    dest[key] = dest.get(key, 0) + words * multiple
        states = following
    counts = [[0] * (regions + 1) for _ in range(width + 1)]
    for state in states:
        for (c, v), amount in state.items():
            counts[c][v] += amount
    for c, row in enumerate(counts):
        assert sum(row) == comb(width, c), (regions, c)
    return counts


def check_small_accumulator() -> None:
    regions = 2
    observed = accumulator_packet_counts(regions)
    brute = [[0] * (regions + 1) for _ in range(4 * regions + 1)]
    for input_word in range(1 << (4 * regions)):
        output = 0
        carry = 0
        for i in range(4 * regions):
            carry ^= (input_word >> i) & 1
            output |= carry << i
        active = sum(bool((output >> (4 * packet)) & 15)
                     for packet in range(regions))
        brute[input_word.bit_count()][active] += 1
    assert observed == brute
    golay_regions = 6
    packet_only = expected_packet_spectrum(
        golay_regions, accumulator_packet_counts(golay_regions))
    joint = good_low_packet_spectrum(
        golay_regions, golay_regions, 0,
        first_accumulator_weights(golay_base_counts(1)))
    assert all(abs(a - b) <= 1e-10 * max(1, a)
               for a, b in zip(packet_only, joint))


def golay_base_counts(blocks: int) -> list[int]:
    counts = [1]
    for _ in range(blocks):
        following = [0] * (len(counts) + 24)
        for weight, amount in enumerate(counts):
            if amount:
                for added, multiple in GOLAY.items():
                    following[weight + added] += amount * multiple
        counts = following
    assert sum(counts) == 1 << (12 * blocks)
    return counts


def first_accumulator_weights(base: list[int]) -> list[float]:
    width = len(base) - 1
    counts = [0.0] * (width + 1)
    counts[0] = float(base[0])
    for a in range(1, width + 1):
        if not base[a]:
            continue
        half = (a + 1) // 2
        scale = base[a] / comb(width, a)
        for c in range(half, width - (a - half) + 1):
            counts[c] += (scale * comb(c - 1, half - 1)
                          * comb(width - c, a - half))
    assert abs(sum(counts) / sum(base) - 1) < 1e-10
    return counts


def source_block_tail_breakdown(regions: int, counts: list[list[int]],
                                min_packet_weight: int,
                                max_source_blocks: int) -> list[float]:
    """Expected bad words by exact number of active Golay source blocks."""
    width = 4 * regions
    base_blocks = width // 24
    bad_given_transitions = [
        sum(row[:min_packet_weight]) / comb(width, c)
        for c, row in enumerate(counts)
    ]
    counts_by_weight = [1]
    result = []
    for j in range(1, min(max_source_blocks, base_blocks) + 1):
        following = [0] * (len(counts_by_weight) + 24)
        for weight, amount in enumerate(counts_by_weight):
            if amount:
                for added, multiple in GOLAY.items():
                    if added:
                        following[weight + added] += amount * multiple
        counts_by_weight = following
        base = [0] * (width + 1)
        scale = comb(base_blocks, j)
        for a, amount in enumerate(counts_by_weight):
            base[a] = scale * amount
        first = first_accumulator_weights(base)
        result.append(sum(amount * bad_given_transitions[c]
                          for c, amount in enumerate(first)))
    return result


def expected_packet_spectrum(regions: int,
                             counts: list[list[int]]) -> list[float]:
    width = 4 * regions
    assert width % 24 == 0
    first = first_accumulator_weights(golay_base_counts(width // 24))
    spectrum = [0.0] * (regions + 1)
    for c, amount in enumerate(first):
        if not amount:
            continue
        scale = amount / comb(width, c)
        for v, number in enumerate(counts[c]):
            spectrum[v] += scale * number
    assert abs(sum(spectrum) / (1 << (width // 2)) - 1) < 1e-10
    assert abs(spectrum[0] - 1) < 1e-10
    return spectrum


def good_low_packet_spectrum(regions: int, max_packets: int,
                             min_binary_weight: int,
                             first: list[float]) -> list[float]:
    """Joint (transition, packet, deficit) count for the low-packet tail.

    Counts only output words with v<=max_packets and binary weight
    w=4v-deficit >= min_binary_weight. This is a binary64 screen.
    """
    width = 4 * regions
    assert 0 <= max_packets <= regions
    deficit_cap = max(0, 4 * max_packets - min_binary_weight)
    shape = (2, max_packets + 1, deficit_cap + 1, width + 1)
    current = np.zeros(shape)
    current[0, 0, 0, 0] = 1.0
    patterns = packet_patterns_with_deficit()
    for epoch in range(regions):
        following = np.zeros(shape)
        for previous in range(2):
            for end, transitions, active, deficit, count in patterns[previous]:
                v_count = min(epoch, max_packets - active) + 1
                d_count = deficit_cap + 1 - deficit
                c_count = 4 * epoch + 1
                if v_count <= 0 or d_count <= 0:
                    continue
                following[end,
                          active:active + v_count,
                          deficit:deficit + d_count,
                          transitions:transitions + c_count] += (
                              count * current[previous,
                                              :v_count, :d_count, :c_count])
        current = following
    slice_weights = np.array([amount / comb(width, c)
                              for c, amount in enumerate(first)])
    spectrum = [0.0] * (max_packets + 1)
    for v in range(max_packets + 1):
        largest_deficit = min(deficit_cap, 4 * v - min_binary_weight)
        if largest_deficit >= 0:
            spectrum[v] = float(np.einsum(
                'dc,c->', current[:, v, :largest_deficit + 1, :].sum(axis=0),
                slice_weights,
            ))
    return spectrum


def one_block_moment(regions: int, gamma: float,
                     outer_blocks: int | None = None) -> list[float]:
    """Continuum transfer or finite-step live-state-count transfer.

    In the finite case, the returned moment uses q for every live step.
    Actual output weight adds an upper correction exp(128*theta*regions)
    in the caller for active steps and their immediate successors.
    """
    if outer_blocks is None:
        zero_to_live = -expm1(-gamma) / gamma if gamma else 1.0
        live_to_live = exp(-gamma)
        live_to_zero = zero_to_live / STATE_NONZERO
    else:
        if outer_blocks < 16 or outer_blocks % 16:
            raise ValueError("finite outer block count must be a positive multiple of 16")
        theta = gamma / (4 * outer_blocks * LIVE_BIT_DENSITY)
        q_minus_one = sum(count * expm1(-theta * weight)
                          for weight, count in STATE_OUTPUT_WEIGHTS.items()) / STATE_NONZERO
        log_q = log1p(q_minus_one)
        steps_per_region = outer_blocks // 16
        live_to_live = exp(steps_per_region * log_q)
        zero_to_live = (expm1(steps_per_region * log_q)
                        / (steps_per_region * expm1(log_q))) if gamma else 1.0
        live_to_zero = exp(log_q) * zero_to_live / STATE_NONZERO
    zero = [1.0] + [0.0] * regions
    live = [0.0] * (regions + 1)
    for epoch in range(regions):
        next_zero = [0.0] * (regions + 1)
        next_live = [0.0] * (regions + 1)
        for v in range(epoch + 1):
            next_zero[v] += zero[v]
            next_live[v] += live[v] * live_to_live
            next_live[v + 1] += zero[v] * zero_to_live
            next_zero[v + 1] += live[v] * live_to_zero
            next_live[v + 1] += live[v] * live_to_live * (1 - 1 / STATE_NONZERO)
        zero, live = next_zero, next_live
    return [(zero[v] + live[v]) / comb(regions, v)
            for v in range(regions + 1)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", type=int, default=54)
    parser.add_argument("--delta", type=float, default=0.11)
    parser.add_argument("--input-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-blocks", type=int, default=0,
                        help="use the finite-step transfer for this multiple-of-16 L")
    parser.add_argument("--outer-weight-floor", type=float, default=0.104,
                        help="certified binary outer relative-weight floor")
    parser.add_argument("--packet-floor-fraction", type=float, default=0.0,
                        help="additional packet-support floor as a fraction of regions")
    parser.add_argument("--joint-max-packets", type=int, default=0,
                        help="replace packet-only counts through this support by joint good-weight counts")
    parser.add_argument("--source-block-breakdown", type=int, default=0,
                        help="split bad packet setup expectation by first Golay source blocks")
    args = parser.parse_args()
    if args.regions < 1 or (4 * args.regions) % 24:
        raise ValueError("regions must give a positive Golay-multiple outer length")
    if not 0 <= args.outer_weight_floor <= 1:
        raise ValueError("outer weight floor must lie in [0,1]")
    if not 0 <= args.packet_floor_fraction <= 1:
        raise ValueError("packet floor fraction must lie in [0,1]")
    if not 0 <= args.joint_max_packets <= args.regions:
        raise ValueError("joint packet cutoff must lie in [0,regions]")
    if args.outer_blocks and (args.outer_blocks < 16 or args.outer_blocks % 16):
        raise ValueError("outer blocks must be a positive multiple of 16")
    check_small_accumulator()
    counts = accumulator_packet_counts(args.regions)
    spectrum = expected_packet_spectrum(args.regions, counts)
    message_bits_per_block = 2 * args.regions
    outer_blocks = (args.outer_blocks if args.outer_blocks else
                    args.input_bits / message_bits_per_block)
    input_bits = (args.outer_blocks * message_bits_per_block
                  if args.outer_blocks else args.input_bits)
    min_binary_weight = ceil(args.outer_weight_floor * 4 * args.regions)
    min_packet_weight = max(ceil(min_binary_weight / 4),
                            ceil(args.packet_floor_fraction * args.regions))
    expected_bad_packet_words = sum(spectrum[1:min_packet_weight])
    source_breakdown = (source_block_tail_breakdown(
        args.regions, counts, min_packet_weight, args.source_block_breakdown)
        if args.source_block_breakdown else [])
    if args.joint_max_packets:
        first = first_accumulator_weights(
            golay_base_counts((4 * args.regions) // 24))
        joint = good_low_packet_spectrum(args.regions,
                                         args.joint_max_packets,
                                         min_binary_weight, first)
        for v, amount in enumerate(joint):
            assert amount <= spectrum[v] * (1 + 1e-10) + 1e-12
            spectrum[v] = amount

    def log_bound(gamma: float) -> float:
        moment = one_block_moment(args.regions, gamma,
                                  args.outer_blocks or None)
        theta = (gamma / (4 * args.outer_blocks * LIVE_BIT_DENSITY)
                 if args.outer_blocks else 0.0)
        total = sum(spectrum[v] * moment[v] * exp(128 * theta * v)
                    for v in range(max(1, min_packet_weight), args.regions + 1))
        return (log(outer_blocks) + gamma * args.delta * args.regions
                / LIVE_BIT_DENSITY + log(total))

    fit = minimize_scalar(log_bound, bounds=(0.001, 30), method="bounded")
    moment = one_block_moment(args.regions, fit.x,
                              args.outer_blocks or None)
    theta = (fit.x / (4 * args.outer_blocks * LIVE_BIT_DENSITY)
             if args.outer_blocks else 0.0)
    contributions = sorted(
        ((spectrum[v] * moment[v] * exp(128 * theta * v),
          v, spectrum[v], moment[v])
         for v in range(max(1, min_packet_weight), args.regions + 1)),
        reverse=True,
    )
    charge = sum(item[0] for item in contributions)
    block_log2 = fit.fun / log(2) - log(outer_blocks) / log(2)
    print(dict(outer_bits=4 * args.regions, input_bits=input_bits,
               outer_blocks=outer_blocks, gamma=fit.x,
               outer_weight_floor=args.outer_weight_floor,
               min_packet_weight=min_packet_weight,
               expected_bad_packet_words=expected_bad_packet_words,
               log2_expected_bad_packet_words=(log(expected_bad_packet_words) / log(2)
                                               if expected_bad_packet_words else None),
               bad_words_by_source_blocks=source_breakdown,
               joint_max_packets=args.joint_max_packets,
               log2_q1_upper=fit.fun / log(2),
               log2_per_block_charge=block_log2,
               indicative_min_growth_constant=-4 * args.regions / block_log2,
               leading_packet_weights=[(v, round(amount / charge, 4))
                                       for amount, v, _, _ in contributions[:6]],
               scope=("finite ideal-setup formula with binary64 evaluation; not an interval certificate"
                      if args.outer_blocks else
                      "good-outer continuum/binary64 screen only; not a certificate")))


if __name__ == "__main__":
    main()
