"""Exact packet-occupancy census and exploratory Turbo one-block screen.

The accumulator/packet counts are exact integers. The BAA fold and s22
continuum moment use binary64 for discovery only; this is not a distance
certificate. No encoder benchmark is run.
"""

from collections import Counter
from math import ceil, comb, exp, expm1, log
from scipy.optimize import minimize_scalar
import argparse


GOLAY = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
STATE_NONZERO = (1 << 22) - 1
LIVE_BIT_DENSITY = (1 << 21) / STATE_NONZERO


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
    counts[0] = 1.0
    for a in range(1, width + 1):
        if not base[a]:
            continue
        half = (a + 1) // 2
        scale = base[a] / comb(width, a)
        for c in range(half, width - (a - half) + 1):
            counts[c] += (scale * comb(c - 1, half - 1)
                          * comb(width - c, a - half))
    assert abs(sum(counts) / (1 << (width // 2)) - 1) < 1e-10
    return counts


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


def one_block_moment(regions: int, gamma: float) -> list[float]:
    """Limit Laplace transform by exact packet support, L -> infinity."""
    zero_to_live = -expm1(-gamma) / gamma if gamma else 1.0
    live_to_live = exp(-gamma)
    zero = [1.0] + [0.0] * regions
    live = [0.0] * (regions + 1)
    for epoch in range(regions):
        next_zero = [0.0] * (regions + 1)
        next_live = [0.0] * (regions + 1)
        for v in range(epoch + 1):
            next_zero[v] += zero[v]
            next_live[v] += live[v] * live_to_live
            next_live[v + 1] += zero[v] * zero_to_live
            next_zero[v + 1] += live[v] * zero_to_live / STATE_NONZERO
            next_live[v + 1] += live[v] * live_to_live * (1 - 1 / STATE_NONZERO)
        zero, live = next_zero, next_live
    return [(zero[v] + live[v]) / comb(regions, v)
            for v in range(regions + 1)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", type=int, default=54)
    parser.add_argument("--delta", type=float, default=0.11)
    parser.add_argument("--input-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-weight-floor", type=float, default=0.104,
                        help="certified binary outer relative-weight floor")
    args = parser.parse_args()
    if args.regions < 1 or (4 * args.regions) % 24:
        raise ValueError("regions must give a positive Golay-multiple outer length")
    if not 0 <= args.outer_weight_floor <= 1:
        raise ValueError("outer weight floor must lie in [0,1]")
    check_small_accumulator()
    counts = accumulator_packet_counts(args.regions)
    spectrum = expected_packet_spectrum(args.regions, counts)
    message_bits_per_block = 2 * args.regions
    outer_blocks = args.input_bits / message_bits_per_block
    min_binary_weight = ceil(args.outer_weight_floor * 4 * args.regions)
    min_packet_weight = ceil(min_binary_weight / 4)

    def log_bound(gamma: float) -> float:
        moment = one_block_moment(args.regions, gamma)
        total = sum(spectrum[v] * moment[v]
                    for v in range(max(1, min_packet_weight), args.regions + 1))
        return (log(outer_blocks) + gamma * args.delta * args.regions
                / LIVE_BIT_DENSITY + log(total))

    fit = minimize_scalar(log_bound, bounds=(0.001, 30), method="bounded")
    moment = one_block_moment(args.regions, fit.x)
    contributions = sorted(
        ((spectrum[v] * moment[v], v, spectrum[v], moment[v])
         for v in range(max(1, min_packet_weight), args.regions + 1)),
        reverse=True,
    )
    block_log2 = fit.fun / log(2) - log(outer_blocks) / log(2)
    print(dict(outer_bits=4 * args.regions, input_bits=args.input_bits,
               outer_blocks=outer_blocks, gamma=fit.x,
               outer_weight_floor=args.outer_weight_floor,
               min_packet_weight=min_packet_weight,
               log2_q1_upper=fit.fun / log(2),
               log2_per_block_charge=block_log2,
               indicative_min_growth_constant=-4 * args.regions / block_log2,
               leading_packet_weights=[v for _, v, _, _ in contributions[:5]],
               scope="good-outer continuum/binary64 screen only; not a certificate"))


if __name__ == "__main__":
    main()
