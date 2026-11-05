#!/usr/bin/env python3
"""Exact tiny checks of arbitrary-set packet moments and mask symmetrization.

Enumerate all 4096 labeled p=2, B=4, G=2 routes. Within-region packet/lane
symmetries and permutation of MDS coordinates reduce all 65536 scalar masks
to 126 moment representatives. Do not confuse these support symmetries with
symmetries of an inner encoder. No inner or whole-code certificate is tested.
Standard library only; writes no files. Run with python -B; root owns runs.
"""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, factorial

from test_persistent_mds import POINTS, power, rank, require


B, K_LOCAL, P, G = 4, 2, 2, 2
N, K = B * P * G, K_LOCAL * P * G
ALL = (1 << N) - 1
TAUS = (Fraction(3, 2), Fraction(2))
LOCAL_TYPES = tuple(combinations_with_replacement(range(P + 1), G))


def check_local_mds():
    local = tuple(tuple(power(point, degree) for point in POINTS) for degree in range(K_LOCAL))
    cases = 0
    for size in range(B + 1):
        for selected in combinations(range(B), size):
            restricted = tuple(tuple(row[column] for column in selected) for row in local)
            require(rank(restricted) == min(size, K_LOCAL), ("RS restriction rank", selected))
            cases += 1
    return cases


def row_masks(route):
    """Low region bit swaps groups; next two bits swap destination lanes."""
    masks = [0] * (P * G)
    for j, bits in enumerate(route):
        for h in range(G):
            g = h ^ (bits & 1)
            swap = (bits >> (1 + h)) & 1
            for lane in range(P):
                masks[P * g + (lane ^ swap)] |= 1 << (P * (j * G + h) + lane)
    require(sum(mask.bit_count() for mask in masks) == N and sum(masks) == ALL,
            ("row partition", route))
    return tuple(masks)


def dimension(mask, rows):
    return sum(max(0, (mask & row).bit_count() - K_LOCAL) for row in rows)


def generator(route):
    """Independent source-to-destination indexing, using RS evaluations."""
    matrix = [[0] * N for _ in range(K)]
    for g in range(G):
        for lane in range(P):
            for j, bits in enumerate(route):
                h = g ^ (bits & 1)
                destination_lane = lane ^ ((bits >> (1 + h)) & 1)
                column = P * (j * G + h) + destination_lane
                for degree in range(K_LOCAL):
                    matrix[(P * g + lane) * K_LOCAL + degree][column] = power(POINTS[j], degree)
    return matrix


def profile(mask):
    """Moment orbit under region order, regional packet order, and lanes."""
    regions = []
    for j in range(B):
        regions.append(tuple(sorted(((mask >> (P * (j * G + h))) & 3).bit_count()
                                    for h in range(G))))
    return tuple(sorted(regions))


def representative(types):
    result = 0
    for j, counts in enumerate(types):
        for h, count in enumerate(counts):
            result |= ((1 << count) - 1) << (P * (j * G + h))
    return result


def moment(mask, routes, tau):
    return sum((tau ** dimension(mask, rows) for rows in routes), Fraction()) / len(routes)


def envelope(size, tau):
    probability = Fraction(min(size, N - size), N)
    local = sum((Fraction(comb(B, i)) * probability ** i * (1 - probability) ** (B - i)
                 * tau ** (P * max(0, i - K_LOCAL)) for i in range(B + 1)), Fraction())
    return tau ** max(0, size - K) * local ** G


def check_monotone_comparison():
    """Every threshold product for counts 0,1,2 and populations up to four."""
    cases = 0
    for population in range(1, 5):
        labels = tuple(permutations(range(population)))
        for counts in combinations_with_replacement(range(P + 1), population):
            for thresholds in product(range(P + 2), repeat=population):
                successes = sum(all(counts[order[i]] >= thresholds[i] for i in range(population))
                                for order in labels)
                eligible = sorted(sum(value >= threshold for value in counts) for threshold in thresholds)
                numerator = 1
                for i, count in enumerate(eligible):
                    numerator *= max(0, count - i)
                require(successes == numerator, ("nested threshold count", counts, thresholds))
                independent = 1
                for count in eligible:
                    independent *= count
                require(successes * population ** population <= factorial(population) * independent,
                        ("without-replacement comparison", counts, thresholds))
                cases += 1
    # General monotone factors: each has multiple nonnegative increments.
    counts, functions = (0, 1, 2), ((1, 2, 5), (0, 3, 4), (2, 2, 7))
    exact = sum(functions[0][counts[a]] * functions[1][counts[b]] * functions[2][counts[c]]
                for a, b, c in permutations(range(3)))
    require(Fraction(exact, 6) <= Fraction(sum(functions[0]) * sum(functions[1]) * sum(functions[2]), 27),
            "non-indicator monotone comparison")
    return cases + 1


def swap_packet_lanes(mask, swap_bits):
    result = 0
    for packet in range(B * G):
        pair = (mask >> (P * packet)) & 3
        if swap_bits & (1 << packet):
            pair = ((pair & 1) << 1) | (pair >> 1)
        result |= pair << (P * packet)
    return result


def batch_statistics(w_mask, h_mask):
    return tuple((((w_mask >> (j * P * G)) & 15).bit_count(),
                  ((h_mask >> (j * P * G)) & 15).bit_count(),
                  (((w_mask & h_mask) >> (j * P * G)) & 15).bit_count())
                 for j in range(B))


def check_mask_orbits(unshuffled, moments):
    # These are complete joint (W,H) orbits, not sets of realizable inner masks.
    seeds = ((0x5555, 0x5555), (0x5555, 0x9999), (0x36C9, 0x5AA5))
    old_cache = {}
    cases = pairs = 0
    for w_mask, h_mask in seeds:
        orbit = {(swap_packet_lanes(w_mask, bits), swap_packet_lanes(h_mask, bits))
                 for bits in range(1 << (B * G))}
        expected_statistics = batch_statistics(w_mask, h_mask)
        require(all(batch_statistics(w, h) == expected_statistics for w, h in orbit),
                "joint packet permutations changed batch (w,h,o)")
        for tau in TAUS:
            original = Fraction()
            shuffled = Fraction()
            for _, h in orbit:
                key = h, tau
                if key not in old_cache:
                    old_cache[key] = moment(h, unshuffled, tau)
                original += old_cache[key]
                shuffled += moments[profile(h), tau]
            require(original == shuffled, ("complete joint-mask orbit identity", w_mask, h_mask, tau))
            require(original <= len(orbit) * envelope(h_mask.bit_count(), tau),
                    ("mask-summed envelope", w_mask, h_mask, tau))
            cases += 1
        pairs += len(orbit)
    return cases, pairs


def main():
    local_mds = check_local_mds()
    monotone_cases = check_monotone_comparison()
    route_keys = tuple(product(range(8), repeat=B))
    routes = tuple(row_masks(route) for route in route_keys)
    unshuffled = tuple(row_masks(route) for route in product(range(2), repeat=B))
    types = tuple(combinations_with_replacement(LOCAL_TYPES, B))
    require(len(types) == 126, "unexpected symmetry-class count")
    # Check that every scalar mask belongs to one of these classes.
    require({profile(mask) for mask in range(1 << N)} == set(types), "support classes incomplete")
    moments = {}
    complement_cases = 0
    for kind in types:
        mask = representative(kind)
        require(profile(mask) == kind, "representative profile mismatch")
        histogram = [0] * (K + 1)
        for rows in routes:
            d = dimension(mask, rows)
            require(d == mask.bit_count() - K + dimension(ALL ^ mask, rows),
                    ("arbitrary-support complement", mask, rows))
            histogram[d] += 1
            complement_cases += 1
        for tau in TAUS:
            value = sum((count * tau ** d for d, count in enumerate(histogram)), Fraction()) / len(routes)
            require(value <= envelope(mask.bit_count(), tau), ("shuffled arbitrary-set moment", mask, tau, value))
            moments[kind, tau] = value

    rank_cases = 0
    fixtures = (0, ALL, 1, 0x5555, 0x5AA5, 0x36C9, 0x7FFF)
    for bits in range(8):
        route = tuple((bits + j) % 8 for j in range(B))
        matrix = generator(route)
        rows = row_masks(route)
        require(rank(matrix) == K, ("routed generator rank", route))
        for mask in fixtures:
            outside = tuple(i for i in range(N) if not (mask >> i) & 1)
            actual = K - rank(tuple(tuple(row[i] for i in outside) for row in matrix))
            require(actual == dimension(mask, rows), ("rank oracle shortening", route, mask, actual))
            rank_cases += 1

    # Selecting lane zero everywhere exposes G entire source rows without swaps.
    negative = 0x5555
    for tau in TAUS:
        actual = moment(negative, unshuffled, tau)
        require(actual == tau ** (G * K_LOCAL), "unshuffled concentration fixture")
        require(actual > envelope(negative.bit_count(), tau), "unshuffled fixed-mask negative control")
        one_group = (2 * tau ** 2 + 8 * tau + 6) / 16
        require(moments[profile(negative), tau] == one_group ** G,
                "independent lane-swap closed-form moment")
    orbit_cases, orbit_pairs = check_mask_orbits(unshuffled, moments)
    print("PASS shuffled packet route toy interfaces:", dict(
        labeled_routes=len(routes), scalar_masks=1 << N, moment_classes=len(types),
        local_mds_subsets=local_mds, monotone_comparisons=monotone_cases,
        complement_identities=complement_cases, rational_moment_bounds=len(moments),
        direct_rank_fixtures=rank_cases, joint_mask_orbits=orbit_cases,
        joint_mask_pairs=orbit_pairs, unshuffled_fixed_mask_rejected=True))


if __name__ == "__main__":
    main()
