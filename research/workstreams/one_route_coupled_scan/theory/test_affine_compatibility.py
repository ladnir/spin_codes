#!/usr/bin/env python3
"""Exact bounded GF16 tests for AFFINE_COMPATIBILITY.md; no files written.

The toy outer is two independent RS[4,2] rows over embedded GF4, arranged
as fixed two-record packets. The state has width two and four globally
distinct additive cosets. Checks are local proof interfaces, not a distance
certificate. Standard library only. Root owns execution; run with python -B.
"""

from fractions import Fraction
from itertools import product

from test_persistent_mds import MUL, POINTS, Q, dot, inverse, matmul, matvec, power, rank, require
from test_syndrome_rank import matrices, nullspace, partition_rank


B, PACKET, WIDTH = 4, 2, 2
N = B * PACKET
K = N // 2
ALL = (1 << N) - 1


def positions(mask):
    return tuple(i for i in range(N) if mask >> i & 1)


def support(vector):
    return sum(1 << i for i, value in enumerate(vector) if value)


def restrict(matrix, mask):
    chosen = positions(mask)
    return tuple(tuple(row[i] for i in chosen) for row in matrix)


def parity():
    # For all four GF4 evaluation points, RS[4,2] is self-dual.
    result = []
    for lane in range(PACKET):
        for degree in range(B // 2):
            result.append(tuple(power(POINTS[i // PACKET], degree) if i % PACKET == lane else 0
                                for i in range(N)))
    result = tuple(result)
    require(rank(result) == K, "outer parity rank")
    require(all(dot(left, right) == 0 for left in result for right in result),
            "RS generator/parity orthogonality")
    return result


def state_map():
    maps = matrices(WIDTH, B)
    result = [[0] * N for _ in range(N)]
    for t in range(B):
        for s in range(t):
            impulse = matmul(inverse(maps[t]), maps[s])
            for j in range(WIDTH):
                for k in range(WIDTH):
                    result[WIDTH * t + j][WIDTH * s + k] = impulse[j][k]
    return tuple(tuple(row) for row in result)


def compatibility(p, a_mask, w_mask):
    """Left annihilator of image(P_A), evaluated on output coordinates W."""
    a_columns = tuple(tuple(row[i] for row in p) for i in positions(a_mask))
    quotient = nullspace(a_columns, len(p))
    p_w = restrict(p, w_mask)
    return tuple(tuple(dot(row, column) for column in zip(*p_w)) for row in quotient)


def compatible_direct(p, a_mask, y):
    a_columns = tuple(tuple(row[i] for row in p) for i in positions(a_mask))
    py = matvec(p, y)
    return rank(a_columns + (py,)) == rank(a_columns)


def vector_on(mask, values):
    output = [0] * N
    for i, value in zip(positions(mask), values):
        output[i] = value
    return tuple(output)


def check_ranks(p, l):
    fixtures = overlap_cases = 0
    output_masks = tuple(mask for mask in range(1 << N) if mask.bit_count() <= 2) + (0x55, 0xAA, ALL)
    for a_mask in range(1 << N):
        rank_a = rank(restrict(p, a_mask))
        for w_mask in output_masks:
            s_mask = a_mask | w_mask
            rank_s = rank(restrict(p, s_mask))
            outer = compatibility(p, a_mask, w_mask)
            expected = 0
            for lane in range(PACKET):
                row_mask = sum(1 << (PACKET * j + lane) for j in range(B))
                a = (a_mask & row_mask).bit_count()
                extra = (w_mask & ~a_mask & row_mask).bit_count()
                expected += min(extra, max(0, B // 2 - a))
            require(rank(outer) == rank_s - rank_a == expected,
                    ("outer compatibility rank", a_mask, w_mask))
            zeros = restrict(tuple(row for i, row in enumerate(l) if not (a_mask >> i & 1)), w_mask)
            joint = rank(zeros + outer)
            require(max(rank(zeros), rank(outer)) <= joint <= rank(zeros) + rank(outer),
                    "joint rank bounds")
            overlap_mask = a_mask & w_mask
            overlap_rows = restrict(tuple(row for i, row in enumerate(l) if not (a_mask >> i & 1)),
                                    overlap_mask)
            overlap_rank = rank(overlap_rows)
            require(joint >= rank(outer) + overlap_rank, "overlap-column additive rank")
            state_masks = tuple((a_mask >> (WIDTH * t)) & 3 for t in range(B))
            overlap_masks = tuple((overlap_mask >> (WIDTH * t)) & 3 for t in range(B))
            require(partition_rank(state_masks, overlap_masks, WIDTH) <= overlap_rank,
                    "selected-observation overlap rank")
            require(w_mask.bit_count() - joint <= w_mask.bit_count() - expected,
                    "compatibility dimension bound")
            overlap_cases += joint < rank(zeros) + rank(outer)
            fixtures += 1
    require(overlap_cases > 0, "negative control failed to detect overlapping ranks")
    return fixtures, overlap_cases


def check_affine_probability(p):
    fixtures = keys = 0
    for a_mask in (0, 1, 3, 0x11):
        a_positions = positions(a_mask)
        rank_a = rank(restrict(p, a_mask))
        v = tuple(power(2, i) if a_mask >> i & 1 else 0 for i in range(N))
        for w_mask in (0, 1, 3, 0x11, 0x80):
            for values in product((0, 1, 2), repeat=w_mask.bit_count()):
                y = vector_on(w_mask, values)
                successes = 0
                for randomizers in product(range(Q), repeat=len(a_positions)):
                    x = list(y)
                    for i, key in zip(a_positions, randomizers):
                        x[i] ^= MUL[key][v[i]]
                    successes += not any(matvec(p, x))
                    keys += 1
                expected = Fraction(1, Q ** rank_a) if compatible_direct(p, a_mask, y) else Fraction()
                require(Fraction(successes, Q ** len(a_positions)) == expected,
                        ("affine membership probability", a_mask, w_mask, values))
                fixtures += 1
    # A=empty does not make the affine fiber automatically meet the outer.
    last_only = vector_on(1 << (N - 1), (1,))
    require(not compatible_direct(p, 0, last_only), "incompatible zero-randomness fixture")
    require(rank(restrict(p, 1 << (N - 1))) == 1, "zero-dimension guard fixture")
    return fixtures, keys


def check_joint_counts(p, l):
    fixtures = assignments = zero_classes = 0
    for a_mask in (0, 1, 3, 0x0F, 0x55, 0xAA, 0xF0, ALL):
        for w_mask in (0, 1, 3, 0x11, 0x81):
            outer = compatibility(p, a_mask, w_mask)
            zeros = restrict(tuple(row for i, row in enumerate(l) if not (a_mask >> i & 1)), w_mask)
            dimension = w_mask.bit_count() - rank(zeros + outer)
            linearized = exact_nonzero = 0
            for values in product(range(Q), repeat=w_mask.bit_count()):
                y = vector_on(w_mask, values)
                v = matvec(l, y)
                state_ok = not (support(v) & ~a_mask)
                outer_ok = compatible_direct(p, a_mask, y)
                equations_ok = not any(dot(row, values) for row in zeros + outer)
                require(equations_ok == (state_ok and outer_ok), "joint matrix semantics")
                linearized += equations_ok
                exact_nonzero += (any(y) and support(y) == w_mask and support(v) == a_mask and outer_ok)
                assignments += 1
            require(linearized == Q ** dimension, ("joint linearized dimension", a_mask, w_mask))
            require(exact_nonzero % (Q - 1) == 0, "exact masks must be projective")
            require(exact_nonzero <= Q ** dimension - 1, "nonzero projective numerator")
            rank_a = rank(restrict(p, a_mask))
            rank_s = rank(restrict(p, a_mask | w_mask))
            expected_lines = Fraction(exact_nonzero, (Q - 1) * Q ** rank_a)
            require(expected_lines <= Fraction(Q ** w_mask.bit_count(), (Q - 1) * Q ** rank_s),
                    "union-support projective bound")
            if dimension == 0:
                require(exact_nonzero == 0, "zero joint dimension must be excluded")
                zero_classes += 1
            fixtures += 1
    require(zero_classes > 0, "zero-dimensional classes missing")
    return fixtures, assignments, zero_classes


def main():
    p, l = parity(), state_map()
    ranks, overlap = check_ranks(p, l)
    affine, keys = check_affine_probability(p)
    counts, assignments, zero = check_joint_counts(p, l)
    print("PASS affine compatibility toy checks:", dict(
        joint_rank_fixtures=ranks, overlapping_rank_controls=overlap,
        affine_probability_fixtures=affine, randomizer_assignments=keys,
        exact_joint_count_fixtures=counts, output_assignments=assignments,
        zero_dimensional_classes=zero))


if __name__ == "__main__":
    main()
