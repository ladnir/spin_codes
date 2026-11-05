#!/usr/bin/env python3
"""Exact small-field checks for the whole-packet route moment.

Two groups each contain two independent RS[4,2] rows over GF(4), viewed
inside GF(16). Each of the four regions permutes its two group packets.
The 16 routes and all 256 whole-packet supports are enumerated exactly.
Packet components are independent outer rows, not SIMD payload lanes.

The arithmetic oracle comes from test_persistent_mds. The frozen rational
moment helper is imported from stratified_certificate; its float proposal
code is not called. These toy checks do not establish a distance certificate
or a performance result. No files are written. Run with python -B.
"""

from fractions import Fraction
from itertools import combinations, product

from stratified_certificate import exact_excess
from test_persistent_mds import POINTS, power, rank, require


B = 4
LOCAL_K = 2
PACKET = 2
GROUPS = 2
K = GROUPS * PACKET * LOCAL_K
N = GROUPS * PACKET * B
PACKETS = GROUPS * B
ALL_PACKETS = (1 << PACKETS) - 1
TAU = Fraction(3, 2)


def generator(route):
    """Rows are message coordinates; columns are routed scalar outputs."""
    require(len(route) == B and all(bit in (0, 1) for bit in route), "invalid route")
    rows = [[0] * N for _ in range(K)]
    for region in range(B):
        for destination in range(GROUPS):
            group = destination ^ route[region]
            for component in range(PACKET):
                column = (region * GROUPS + destination) * PACKET + component
                for degree in range(LOCAL_K):
                    message = (group * PACKET + component) * LOCAL_K + degree
                    rows[message][column] = power(POINTS[region], degree)
    return tuple(tuple(row) for row in rows)


def scalar_support(mask):
    return tuple(column for column in range(N) if mask & (1 << (column // PACKET)))


def shortened_dimension(matrix, selected):
    selected = frozenset(selected)
    require(all(0 <= column < N for column in selected), "invalid scalar support")
    outside = tuple(column for column in range(N) if column not in selected)
    restricted = tuple(tuple(row[column] for column in outside) for row in matrix)
    return K - rank(restricted)


def group_occupancies(route, mask):
    occupancy = [0] * GROUPS
    for region in range(B):
        for destination in range(GROUPS):
            if mask & (1 << (region * GROUPS + destination)):
                occupancy[destination ^ route[region]] += 1
    return tuple(occupancy)


def interval_count(mask):
    """Number of maximal selected runs in the ordered scalar coordinates."""
    starts = mask & ~(mask << 1)
    return starts.bit_count()


def check_local_mds():
    local = tuple(tuple(power(point, degree) for point in POINTS) for degree in range(LOCAL_K))
    cases = 0
    for size in range(B + 1):
        for selected in combinations(range(B), size):
            restricted = tuple(tuple(row[column] for column in selected) for row in local)
            require(rank(restricted) == min(size, LOCAL_K), ("RS restriction rank", selected))
            cases += 1
    return cases


def main():
    local_cases = check_local_mds()
    routes = tuple(product((0, 1), repeat=B))
    matrices = tuple(generator(route) for route in routes)
    totals = [Fraction(0) for _ in range(1 << PACKETS)]
    tau_powers = tuple(TAU ** dimension for dimension in range(K + 1))
    dimension_cases = complement_cases = 0

    for route, matrix in zip(routes, matrices):
        require(rank(matrix) == K, ("routed generator lost dimension", route))
        current = []
        for mask in range(1 << PACKETS):
            selected = scalar_support(mask)
            dimension = shortened_dimension(matrix, selected)
            expected = PACKET * sum(max(0, count - LOCAL_K)
                                    for count in group_occupancies(route, mask))
            require(dimension == expected,
                    ("packet shortening formula", route, mask, dimension, expected))
            require(dimension % PACKET == 0, "whole-packet shortening must be a multiple of packet width")
            current.append(dimension)
            totals[mask] += tau_powers[dimension]
            dimension_cases += 1
        for mask, dimension in enumerate(current):
            complement = current[ALL_PACKETS ^ mask]
            require(dimension == PACKET * mask.bit_count() - K + complement,
                    ("complement identity", route, mask, dimension, complement))
            complement_cases += 1

    moments = tuple(total / len(routes) for total in totals)
    require(moments[0] == 1, "empty support must have moment one")
    require(moments[ALL_PACKETS] == TAU ** K, "full-support moment")
    moment_cases = 0
    for mask in range(1, 1 << PACKETS):
        selected = PACKET * mask.bit_count()
        intervals = interval_count(mask)
        require(intervals > 0, "nonempty support needs an interval")
        # Raising tau to PACKET accounts for the independent row dimensions.
        # The number of independently permuted labels is GROUPS, not GROUPS*PACKET.
        excess = exact_excess(B, intervals, selected, N,
                              TAU.numerator ** PACKET, TAU.denominator ** PACKET)
        bound = (1 + excess) ** GROUPS
        require(moments[mask] <= bound,
                ("packet route moment envelope", mask, intervals, moments[mask], bound))
        moment_cases += 1

    # A partial-packet support selecting three coordinates of just one RS row
    # has dimension one. It cannot obey PACKET times an integer shortening law.
    identity = matrices[routes.index((0, 0, 0, 0))]
    partial = tuple((region * GROUPS) * PACKET for region in range(3))
    partial_dimension = shortened_dimension(identity, partial)
    require(partial_dimension == 1, ("partial-packet negative control", partial_dimension))
    require(any(0 < sum(column in partial for column in range(start, start + PACKET)) < PACKET
                for start in range(0, N, PACKET)), "negative control accidentally used whole packets")
    full_packet_dimension = PACKET * max(0, 3 - LOCAL_K)
    require(partial_dimension != full_packet_dimension, "whole-packet hypothesis was not distinguished")
    partial_complement = shortened_dimension(identity, set(range(N)).difference(partial))
    require(partial_dimension == len(partial) - K + partial_complement,
            "partial support still obeys the scalar complement identity")

    # Also reject the tempting omission of tau's packet-width exponent.
    unscaled = (1 + exact_excess(B, 1, N, N, TAU.numerator, TAU.denominator)) ** GROUPS
    require(unscaled < moments[ALL_PACKETS], "unscaled marker negative control")

    print("PASS whole-packet route toy interfaces:", dict(
        routes=len(routes), packet_supports=1 << PACKETS, local_mds_subsets=local_cases,
        generator_dimensions=dimension_cases, complement_identities=complement_cases,
        exact_moment_bounds=moment_cases, partial_packet_dimension=partial_dimension,
        unscaled_marker_rejected=True))


if __name__ == "__main__":
    main()
