"""Exact toy checks for the whole-state bridge, not a distance certificate.

GF4 is embedded in GF16; the two-batch inner uses A=C^T=(1,1).
One nonzero GF16 refresh is sampled once for all messages. The tests count
existence events, not independent refreshes for different messages.
Run with Python -B under a 30-second, 256-MiB process cap.
"""
from fractions import Fraction as F
from itertools import product
from math import comb


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        if a & 16:
            a ^= 0x13
        b >>= 1
    return r


E = (0, 1, 6, 7)
FIELD = tuple(range(16))


def linear(row, m):
    r = 0
    for a, b in zip(row, m):
        r ^= mul(a, b)
    return r


def masks(zeros):
    value = zeros
    while True:
        yield value
        if value == 0:
            break
        value = (value - 1) & zeros


def zeros_of(word):
    return sum(1 << i for i, value in enumerate(word) if value == 0)


def kernel_counts(generator, active, field):
    """Count relaxed pairs with source != 0, and their zero-source fibers."""
    counts = [[0] * 16 for _ in range(2)]
    fibers = [[0] * 16 for _ in range(2)]
    all_active = [[0] * 16 for _ in range(2)]
    joint_weighted = [[F(0)] * 16 for _ in range(2)]
    for m in product(field, repeat=len(active)):
        source = [0, 0]
        for j, value in zip(active, m):
            source[j] = value
        x = [linear(row, source) for row in generator]
        nonzero = any(m)
        for live in (0, 1):
            if live == 0 and x[0] ^ x[1]:
                continue
            for state in field if live else (0,):
                y = (x[0], x[1], x[2] ^ state, x[3] ^ state)
                for z in masks(zeros_of(y)):
                    (counts if nonzero else fibers)[live][z] += 1
                    if all(m):
                        all_active[live][z] += 1
                    if nonzero and field == E:
                        missing = sum(value == 0 for value in m)
                        origin = int(live and state == 0)
                        joint_weighted[live][z] += F(1, 3**missing * 5**origin)
    return counts, fibers, all_active, joint_weighted


def actual_events(generator, active):
    events = [[0] * 16 for _ in range(2)]
    bad = [0] * 5
    for refresh in FIELD[1:]:
        seen = [set(), set()]
        minimum = 5
        for m in product(FIELD[1:], repeat=len(active)):
            source = [0, 0]
            for j, value in zip(active, m):
                source[j] = value
            x = [linear(row, source) for row in generator]
            state = mul(refresh, x[0] ^ x[1])
            y = (x[0], x[1], x[2] ^ state, x[3] ^ state)
            seen[bool(state)].update(masks(zeros_of(y)))
            minimum = min(minimum, sum(value != 0 for value in y))
        for live in (0, 1):
            for z in seen[live]:
                events[live][z] += 1
        for threshold in range(5):
            bad[threshold] += minimum <= threshold
    return events, bad


def main():
    assert {mul(a, b) for a in E for b in E} == set(E)
    assert {a ^ b for a in E for b in E} == set(E)
    fixtures = (
        ((1, 0), (0, 1), (1, 1), (1, 6)),
        ((1, 1), (1, 1), (1, 0), (0, 1)),
        ((1, 0), (6, 0), (0, 1), (0, 7)),
    )
    events_checked = markers_checked = source_checked = fibers_checked = 0
    empty_checked = 0
    joint_checked = 0
    for generator, active in product(fixtures, ((0,), (1,), (0, 1))):
        ec, ef, egood, weighted = kernel_counts(generator, active, E)
        fc, ff, _, _ = kernel_counts(generator, active, FIELD)
        events, bad = actual_events(generator, active)
        if all(generator[i][j] == 0 for i in (0, 1) for j in active):
            assert not any(events[1])
            empty_checked += 1
        for live, z in product((0, 1), range(16)):
            # A subspace dimension and its m=0 fiber extend by degree two.
            assert ec[live][z] % 3 == 0 and fc[live][z] % 15 == 0
            eproj, fproj = ec[live][z] // 3, fc[live][z] // 15
            assert fproj <= eproj**2
            assert ff[live][z] == ef[live][z]**2
            fibers_checked += 1
            probability = F(events[live][z], 15)
            expected = F(fproj, 15**live)
            assert probability <= min(1, expected)
            assert probability**2 <= F(16, 15)**live * F(eproj, 4**live)**2
            events_checked += 1
            if events[live][z]:
                # Every actual source group has positive projection rank.
                # Union bound over missing auxiliary groups; no independence.
                c = 1 - F(len(active), 4)
                assert egood[live][z] >= c * ec[live][z]
                assert ec[live][z]**4 <= 3**len(active) * 5**live * weighted[live][z]**4
                joint_checked += 1
                source_checked += 1
        for mark, threshold in product((F(1, 4), F(1, 2), F(3, 4)), range(5)):
            averaged = sum(mark**(4-z.bit_count()) * (1-mark)**z.bit_count()
                           * F(events[live][z], 15)
                           for live, z in product((0, 1), range(16)))
            assert F(bad[threshold], 15) <= averaged / mark**threshold
            markers_checked += 1
    # Degree-four projective comparison used at GF256 -> GF(2^32).
    projective_checked = 0
    for q, d, d0 in product((2, 4, 16, 256), range(1, 10), range(10)):
        if d0 >= d:
            continue
        larger = (q**(4*d) - q**(4*d0)) // (q**4 - 1)
        smaller = (q**d - q**d0) // (q - 1)
        assert larger <= smaller**4
        projective_checked += 1
    # Generic positive-marker identity, including all-zero outputs.
    for n, weight in product(range(1, 9), range(9)):
        if weight > n:
            continue
        for z in (F(1, 4), F(1, 2), F(3, 4)):
            assert sum(comb(n-weight, k)*z**(n-k)*(1-z)**k
                       for k in range(n-weight+1)) == z**weight
    # A random injective one-symbol -> two-symbol E map has expected
    # nonzero output counting mass 1/(q+1) at each nonzero output vector.
    # Exact zero-group strata are disjoint, even after padding missing groups.
    stratum_checked = 0
    for q, r in product((2, 4, 16, 256), range(1, 9)):
        beta = F(q*q, q+1)
        kappa = F(q*q, 1)/beta
        for h in range(r):
            for zero_pattern in range(1 << r):
                j = zero_pattern.bit_count()
                if j <= h:
                    density = beta**(r-j) * (q*q)**j
                    assert density <= beta**r * kappa**h
                    stratum_checked += 1
    print(dict(exact=True, whole_code_certificate=False,
               fixed_zero_path_events=events_checked,
               zero_source_fibers=fibers_checked,
               retained_actual_support_checks=source_checked,
               empty_input_path_checks=empty_checked,
               root_before_marker_checks=markers_checked,
               degree_four_projective_checks=projective_checked,
               joint_source_state_jensen_checks=joint_checked,
               exact_missing_group_strata=stratum_checked))


if __name__ == '__main__':
    main()
