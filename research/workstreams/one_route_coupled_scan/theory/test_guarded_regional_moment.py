#!/usr/bin/env python3
"""Exact bounded toy checks for guarded regional packet moments.

Standard library only; writes no files. Intended external budget: 30 seconds,
256 MiB. The root agent owns execution; use python -B. These checks enumerate
4096 shuffled routes for a few fixed nested supports, not a whole-code sum.
"""

from fractions import Fraction
from itertools import product
from math import comb

from test_shuffled_packet_route import (
    ALL, B, G, K_LOCAL, N, P, TAUS, dimension, row_masks,
    swap_packet_lanes,
)
from test_persistent_mds import require


M = P * G
ETAS = (Fraction(3, 2), Fraction(2))
FIXTURES = (
    (0, 0),
    (0, ALL),
    (ALL, ALL),
    (0x0555, 0x0FFF),
    (0x1111, 0x3333),
    (0x3333, 0x3333),
    (0x1155, 0x33FF),
    (0x1525, 0x3F3F),
)


def whole_packets(mask):
    packet_mask = (1 << P) - 1
    return all(((mask >> (P * j)) & packet_mask) in (0, packet_mask)
               for j in range(B * G))


def marginals(selected, allowed):
    require(selected & ~allowed == 0, "T is not a subset of U")
    region_mask = (1 << M) - 1
    alpha = tuple(Fraction(((selected >> (M * j)) & region_mask).bit_count(), M)
                  for j in range(B))
    upsilon = tuple(Fraction(((allowed >> (M * j)) & region_mask).bit_count(), M)
                    for j in range(B))
    return alpha, upsilon


def pb_law(alpha, upsilon):
    law = {(0, 0): Fraction(1)}
    for a, u in zip(alpha, upsilon):
        require(0 <= a <= u <= 1, "invalid nested marginal")
        updated = {}
        for (z, q), probability in law.items():
            for dz, dq, factor in ((0, 0, 1-u), (0, 1, u-a), (1, 1, a)):
                key = z + dz, q + dq
                updated[key] = updated.get(key, Fraction()) + probability * factor
        law = updated
    require(sum(law.values(), Fraction()) == 1, "PB law not normalized")
    require(all(z <= q for z, q in law), "nested PB support")
    require(sum((z * probability for (z, _), probability in law.items()), Fraction())
            == sum(alpha), "T marginal mean")
    require(sum((q * probability for (_, q), probability in law.items()), Fraction())
            == sum(upsilon), "U marginal mean")
    return law


def exact_comparison(alpha, upsilon, tau):
    t = tau ** P
    law = pb_law(alpha, upsilon)
    h = sum((probability * t ** max(0, z-K_LOCAL)
             for (z, _), probability in law.items()), Fraction())
    j = sum((probability * t ** max(0, z-K_LOCAL)
             for (z, q), probability in law.items() if q > K_LOCAL), Fraction())
    return h, j


def multiply(values):
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def separable(alpha, upsilon, tau, eta):
    t = tau ** P
    aa = tuple(1 + (t-1)*a for a in alpha)
    bb = tuple(1 + (eta-1)*u for u in upsilon)
    cc = tuple(b + eta*(t-1)*a for a, b in zip(alpha, bb))
    h_hat = 1 + t ** (-K_LOCAL) * multiply(aa)
    j_hat = eta ** (-K_LOCAL-1) * (multiply(bb) + t ** (-K_LOCAL) * multiply(cc))
    expanded = G * eta ** (-K_LOCAL-1) * sum(
        (comb(G-1, r) * t ** (-K_LOCAL*r)
         * (multiply(a**r*b for a, b in zip(aa, bb))
            + t ** (-K_LOCAL) * multiply(a**r*c for a, c in zip(aa, cc)))
         for r in range(G)), Fraction())
    require(expanded == G * j_hat * h_hat ** (G-1), "positive expansion mismatch")
    return h_hat, j_hat, expanded


def guarded_histogram(selected, allowed, routes):
    histogram = [0] * (N + 1)
    for rows in routes:
        if dimension(allowed, rows) > 0:
            histogram[dimension(selected, rows)] += 1
    return histogram


def evaluate_histogram(histogram, tau, denominator):
    return sum((count * tau**d for d, count in enumerate(histogram)), Fraction()) / denominator


def main():
    routes = tuple(row_masks(route) for route in product(range(8), repeat=B))
    unshuffled = tuple(row_masks(route) for route in product(range(2), repeat=B))
    require(len(routes) == 4096, "unexpected route count")
    rational_bounds = positive_bounds = orbit_checks = 0
    for selected, allowed in FIXTURES:
        require(whole_packets(allowed), ("fixture U cuts a packet", allowed))
        alpha, upsilon = marginals(selected, allowed)
        histogram = guarded_histogram(selected, allowed, routes)
        for tau in TAUS:
            actual = evaluate_histogram(histogram, tau, len(routes))
            h, j = exact_comparison(alpha, upsilon, tau)
            require(actual <= min(h**G, G*j*h**(G-1)),
                    ("guarded moment", selected, allowed, tau, actual, h, j))
            rational_bounds += 1
            for eta in ETAS:
                h_hat, j_hat, expanded = separable(alpha, upsilon, tau, eta)
                require(h <= h_hat and j <= j_hat, ("local envelope", selected, allowed, tau, eta))
                require(G*j*h**(G-1) <= expanded, "global positive envelope")
                positive_bounds += 1

        # Complete T-mask orbit, with U invariant because it is whole-packet.
        orbit = {swap_packet_lanes(selected, bits) for bits in range(1 << (B * G))}
        require(all(swap_packet_lanes(allowed, bits) == allowed
                    for bits in range(1 << (B * G))), "U not fixed by lane permutations")
        summed_histogram = [0] * (N + 1)
        for mask in orbit:
            counts = guarded_histogram(mask, allowed, unshuffled)
            summed_histogram = [a + b for a, b in zip(summed_histogram, counts)]
        for tau in TAUS:
            original = evaluate_histogram(summed_histogram, tau, len(unshuffled))
            shuffled = len(orbit) * evaluate_histogram(histogram, tau, len(routes))
            require(original == shuffled, ("guarded orbit identity", selected, allowed, tau))
            orbit_checks += 1

    # Invalid extension: a scalar-cut U gives different occupancies to the p
    # rows of one group. The common-indicator argument is then false.
    selected, allowed = 0, 0x5555
    require(not whole_packets(allowed), "negative control did not cut packets")
    alpha, upsilon = marginals(selected, allowed)
    histogram = guarded_histogram(selected, allowed, routes)
    actual = evaluate_histogram(histogram, Fraction(2), len(routes))
    h, j = exact_comparison(alpha, upsilon, Fraction(2))
    false_extension = G*j*h**(G-1)
    require(actual == Fraction(55, 64), ("scalar-cut exact probability", actual))
    require(false_extension == Fraction(5, 8) and actual > false_extension,
            ("scalar-cut U must reject shared-group extension", actual, false_extension))
    print("PASS guarded regional moment toy interfaces:", dict(
        labeled_routes=len(routes), nested_support_fixtures=len(FIXTURES),
        rational_guarded_bounds=rational_bounds, positive_envelopes=positive_bounds,
        guarded_orbit_identities=orbit_checks, scalar_cut_guard_rejected=True))


if __name__ == "__main__":
    main()
