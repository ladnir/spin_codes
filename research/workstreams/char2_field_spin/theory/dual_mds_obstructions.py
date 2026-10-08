"""Exact determinant, rank-obstruction, and scaled-emission checks.

The existing dual-MDS construction is not modified.  eta=1 is the existing
map; eta outside the embedded GF16 is a separately tested local repair.
These checks are not a whole-code distance certificate or a benchmark.
"""
from __future__ import annotations

import argparse
from collections import Counter
from random import Random

from grs16_12_parity_reference import multiply as mul16
from dual_mds_phase import J, curve, global_maps, period_labels
from region64_reference import IDENTITY, matmul, matvec, phase_geometry, random_orientation, transpose
from rate_design_coset import mul, power, xor_sum


M16 = tuple(tuple(mul16(a, b) for b in range(16)) for a in range(16))
INV16 = (0,) + tuple(next(b for b in range(1, 16) if M16[a][b] == 1) for a in range(1, 16))


def pow16(a, n):
    result = 1
    for _ in range(n):
        result = M16[result][a]
    return result


def trace16(a):
    return a ^ pow16(a, 2) ^ pow16(a, 4) ^ pow16(a, 8)


def rank_det16(matrix):
    work = [list(row) for row in matrix]
    rank, determinant = 0, 1
    for column in range(4):
        pivot = next((i for i in range(rank, 4) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        diagonal = work[rank][column]
        determinant = M16[determinant][diagonal]
        work[rank] = [M16[INV16[diagonal]][value] for value in work[rank]]
        for row in range(rank+1, 4):
            coefficient = work[row][column]
            work[row] = [a ^ M16[coefficient][b] for a, b in zip(work[row], work[rank])]
        rank += 1
    return rank, determinant if rank == 4 else 0


def determinant_coefficients(parity, s1, s3, s5):
    """det(I+eta*M*J)=(1+eta*p1+eta^2*p2)^2, over abstract GF16."""
    s13 = pow16(s1, 3)
    p1 = M16[parity][s3] ^ s13
    p2 = M16[parity ^ 1][pow16(s3, 2) ^ M16[s13][s3] ^ pow16(s1, 6) ^ M16[s1][s5]]
    return p1, p2


def subset_checks():
    vectors = tuple((1, u, pow16(u, 2), pow16(u, 3)) for u in range(16))
    contributions = tuple(tuple(tuple(M16[v[i]][v[j]] for j in range(4)) for i in range(4))
                          for v in vectors)
    gram = [[0]*4 for _ in range(4)]
    s1 = s3 = s5 = old_mask = 0
    counts, coefficient_counts, rank_one = Counter(), Counter(), []
    for index in range(1 << 16):
        mask = index ^ (index >> 1)
        if index:
            changed = (mask ^ old_mask).bit_length()-1
            s1 ^= changed
            s3 ^= pow16(changed, 3)
            s5 ^= pow16(changed, 5)
            gram = [[a ^ b for a, b in zip(row, addition)]
                    for row, addition in zip(gram, contributions[changed])]
        old_mask = mask
        weight = mask.bit_count()
        p1, p2 = determinant_coefficients(weight & 1, s1, s3, s5)
        matrix = tuple(tuple(int(i == j) ^ gram[i][3-j] for j in range(4)) for i in range(4))
        rank, determinant = rank_det16(matrix)
        self_root = 1 ^ p1 ^ p2
        assert determinant == M16[self_root][self_root]
        counts[weight, rank] += 1
        coefficient_counts[p1, p2] += 1
        if rank == 1:
            rank_one.append(mask)
        # The quadratic splitting test is exact, not floating arithmetic.
        if p1 and p2:
            ratio = M16[p2][INV16[M16[p1][p1]]]
            assert trace16(ratio) == 0

    subgroup = {u for u in range(16) if pow16(u, 4) == u}
    cosets = {sum(1 << (u ^ a) for u in subgroup) for a in range(16)}
    assert len(cosets) == 4
    assert all(mask in rank_one for mask in cosets)
    assert all(counts[weight, 0] == 0 for weight in range(17))

    embedded = phase_geometry()[0]
    outside = tuple(eta for eta in range(256) if eta not in embedded)
    assert 2 in outside
    for p1, p2 in coefficient_counts:
        roots = [eta for eta in range(16) if 1 ^ M16[eta][p1] ^ M16[M16[eta][eta]][p2] == 0]
        assert len(roots) == (2 if p1 and p2 else 1 if p1 or p2 else 0)
        for eta in outside:
            value = 1 ^ mul(eta, embedded[p1]) ^ mul(mul(eta, eta), embedded[p2])
            assert value != 0
    return counts, coefficient_counts, rank_one


def add(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


def scaled_encode(regions, orientations, eta):
    state, output = (0, 0, 0, 0), []
    phases = phase_geometry()[1]
    for region, orientation in zip(regions, orientations):
        blocks = []
        for index, values in enumerate(region):
            emission, injection = global_maps(orientation, phases[index % 4])
            mask = tuple(mul(eta, value) for value in matvec(emission, state))
            blocks.append(add(values, mask))
            state = add(state, matvec(injection, values))
        output.append(tuple(blocks))
    return tuple(output)


def involution_checks():
    rng = Random(20261007)
    for lengths in ((1,), (4,), (7,), (1, 2, 3, 4, 5)):
        regions = tuple(tuple(tuple(rng.randrange(256) for _ in range(16)) for _ in range(length))
                        for length in lengths)
        orientations = tuple(random_orientation(rng) for _ in lengths)
        for eta in (1, 2):
            assert scaled_encode(scaled_encode(regions, orientations, eta), orientations, eta) == regions


def summary_counterexample():
    labels = period_labels()
    extra_root = next(value for value in range(256) if value not in labels)
    f = (0, 1, 1, 1)  # Reversed coefficients of z^2+z+1.
    g = (1, 1 ^ extra_root, 1 ^ extra_root, extra_root)
    rows = tuple(tuple(reversed(curve(label))) for label in labels)
    roots_f = {a for a, value in zip(labels, matvec(rows, f)) if not value}
    roots_g = {a for a, value in zip(labels, matvec(rows, g)) if not value}
    assert roots_f == roots_g and len(roots_f) == 2
    emission, injection = global_maps(IDENTITY, 0)
    records = []
    for state in (f, g):
        mask = matvec(emission, state)
        output = tuple(value if index in (0, 1) else 0 for index, value in enumerate(mask))
        values = add(mask, output)
        following = add(state, matvec(injection, values))
        records.append((tuple(i for i, value in enumerate(values) if value),
                        tuple(i for i, value in enumerate(output) if value), following))
    assert records[0][:2] == records[1][:2]
    assert len(records[0][0]) == 12 and records[0][1] == (0, 1)
    assert not any(records[0][2]) and any(records[1][2])
    return extra_root, tuple(sorted(roots_f))


def fft_scale_checks():
    rng = Random(18188)
    for _ in range(256):
        c0, c1, c2, c3, a, a2, s, six = (rng.randrange(256) for _ in range(8))
        eta = 2
        def emit(c0, c1, c2, c3, s, six):
            u0, u1 = c0 ^ mul(s, c2), c1 ^ mul(s, c3)
            u2, u3 = u0 ^ mul(six, c2), u1 ^ mul(six, c3)
            e0, e2 = u0 ^ mul(a, u1), u2 ^ mul(a2, u3)
            return e0, e0 ^ u1, e2, e2 ^ u3
        expected = tuple(mul(eta, value) for value in emit(c0, c1, c2, c3, s, six))
        actual = emit(mul(eta, c0), mul(eta, c1), c2, c3, mul(eta, s), mul(eta, six))
        assert actual == expected
        v0, v1, v2, v3 = (rng.randrange(256) for _ in range(4))
        v01, v23 = v0 ^ v1, v2 ^ v3
        w1, w3 = mul(a, v01) ^ v1, mul(a2, v23) ^ v3
        b0, b1 = v01 ^ v23, w1 ^ w3
        b2, b3 = mul(s, b0) ^ mul(six, v23), mul(s, b1) ^ mul(six, w3)
        expected = tuple(mul(eta, value) for value in (b0, b1, b1 ^ b2, b1 ^ b2 ^ b3))
        d0, d1 = mul(eta, b0), mul(eta, b1)
        d2 = mul(mul(eta, s), b0) ^ mul(mul(eta, six), v23)
        d3 = mul(mul(eta, s), b1) ^ mul(mul(eta, six), w3)
        assert (d0, d1, d1 ^ d2, d1 ^ d2 ^ d3) == expected


def check():
    counts, coefficients, rank_one = subset_checks()
    involution_checks()
    extra_root, common_roots = summary_counterexample()
    fft_scale_checks()
    print('All65536 GF16 subsets: determinant identity, rank classification, and splitting checks pass.')
    print('Distinct determinant quadratics:', len(coefficients))
    print('Rank counts for output sizes0..4:', {w: dict(sorted((r, n) for (v, r), n in counts.items() if v == w)) for w in range(5)})
    print('Rank-one subsets by size:', dict(sorted(Counter(mask.bit_count() for mask in rank_one).items())))
    print('All240 eta outsideGF16 avoid every determinant root; eta=2 is valid in all phases.')
    print('Rank/root-mask counterexample: extra root=', extra_root, 'common roots=', common_roots)
    print('Involution and two-extra-XMM emission/adjoint FFT identities pass.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        check()
