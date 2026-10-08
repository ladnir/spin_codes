"""Exact state-basis audit and packed-FFT accounting; no native code/timing."""

from collections import Counter
from functools import lru_cache
from itertools import combinations
import json
from random import Random

from state_capacity_reference import (BETA, LABELS, W8, matrices, matvec,
                                      multiplication_table, rank)
from rate_design_coset import power, xor_sum

MUL = multiplication_table()
INV = tuple(power(x, 254) if x else 0 for x in range(256))


@lru_cache(maxsize=None)
def normalized(level, x):
    numerator = denominator = 1
    for w in W8[:1 << level]:
        numerator = MUL[numerator][x ^ w]
        denominator = MUL[denominator][W8[1 << level] ^ w]
    return MUL[numerator][INV[denominator]]


@lru_cache(maxsize=None)
def fft_ops(n, offset, inverse=False, transposed=False):
    if n == 1:
        return ()
    half = n//2
    c = normalized(half.bit_length()-1, offset)
    children = (fft_ops(half, offset, inverse) +
                tuple((d+half, s+half, a) for d, s, a in
                      fft_ops(half, offset ^ W8[half], inverse)))
    stage = []
    for i in range(half):
        if inverse:
            stage.append((half+i, i, 1))
        if c:
            stage.append((i, half+i, c))
        if not inverse:
            stage.append((half+i, i, 1))
    result = children+tuple(stage) if inverse else tuple(stage)+children
    if transposed:
        return tuple((s, d, c) for d, s, c in reversed(result))
    return result


def evaluate(ops, values):
    out = list(values)
    for d, s, c in ops:
        out[d] ^= MUL[c][out[s]]
    return tuple(out)


def emission(values, phase, coefficient_basis=False):
    coefficients = (tuple(values) if coefficient_basis else
                    evaluate(fft_ops(8, 0, True), values))
    return sum((evaluate(fft_ops(8, LABELS[phase][first] ^ BETA), coefficients)
                for first in (0, 8)), ())


def injection(values, phase, coefficient_basis=False):
    parts = [evaluate(fft_ops(8, LABELS[phase][first] ^ BETA, transposed=True),
                      values[first:first+8]) for first in (0, 8)]
    out = tuple(a ^ b for a, b in zip(*parts))
    return (out if coefficient_basis else
            evaluate(fft_ops(8, 0, True, True), out))


def coefficient_matrix(phase):
    columns = [emission(tuple(int(i == j) for i in range(8)), phase, True)
               for j in range(8)]
    return tuple(zip(*columns))


def packed_fft_cost(offset, inverse=False, transposed=False):
    """Conservative explicit two-ZMM butterfly schedule for eight records.

    Within a ZMM, duplicate each low/high butterfly input with two 128-bit
    lane shuffles. A forward butterfly costs one GFNI and one XOR. Inverse
    costs one GFNI and two XORs, the latter masked on the low coordinates.
    Field transposition exchanges these two costs. An all-zero-coefficient
    block needs only one lane shuffle and one masked XOR. The four-pair
    stage spans the two ZMMs and needs no shuffle. Copy/register renaming
    is excluded; diagonal/eta additions are accounted separately.
    """
    cost = Counter(gfni=0, boolean=0, shuffle=0)
    for half in (1, 2, 4):
        level = half.bit_length()-1
        if half == 4:
            c = normalized(level, offset)
            cost.update(gfni=int(bool(c)), boolean=2 if c else 1)
            continue
        for first in (0, 4):
            cs = [normalized(level, offset ^ W8[block])
                  for block in range(first, first+4, 2*half)]
            if not any(cs):
                cost.update(boolean=1, shuffle=1)
            else:
                cost.update(gfni=1, boolean=1+int(inverse != transposed),
                            shuffle=2)
    return cost


def packed_transform(values, offset, inverse=False, transposed=False):
    """Lane-level interpreter of Cauchy256InnerFft.h, not its scalar DAG."""
    values = list(values)
    forward_like = inverse == transposed
    levels = (4, 2, 1) if forward_like else (1, 2, 4)
    for half in levels:
        if half == 4:
            c = normalized(2, offset)
            pairs = zip(values[:4], values[4:])
            lows, highs = [], []
            for a, b in pairs:
                if not inverse and not transposed:
                    lo = a ^ MUL[c][b]
                    hi = lo ^ b
                elif inverse and not transposed:
                    hi = a ^ b
                    lo = a ^ MUL[c][hi]
                elif not inverse and transposed:
                    lo = a ^ b
                    hi = b ^ MUL[c][lo]
                else:
                    hi = b ^ MUL[c][a]
                    lo = a ^ hi
                lows.append(lo)
                highs.append(hi)
            values = lows+highs
            continue
        for first in (0, 4):
            x = values[first:first+4]
            low = 0xa0 if half == 1 else 0x44
            high = 0xf5 if half == 1 else 0xee
            a = [x[(low >> (2*i)) & 3] for i in range(4)]
            b = [x[(high >> (2*i)) & 3] for i in range(4)]
            cs = [normalized(half.bit_length()-1,
                             offset ^ W8[first+(i//(2*half))*(2*half)])
                  for i in range(4)]
            out = []
            for i in range(4):
                is_low = not i & half
                if not inverse and not transposed:
                    result = a[i] ^ MUL[cs[i] ^ int(not is_low)][b[i]]
                elif inverse and not transposed:
                    d = a[i] ^ b[i]
                    result = a[i] ^ MUL[cs[i]][d] if is_low else d
                elif not inverse and transposed:
                    d = a[i] ^ b[i]
                    result = d if is_low else b[i] ^ MUL[cs[i]][d]
                else:
                    result = b[i] ^ MUL[cs[i] ^ int(is_low)][a[i]]
                out.append(result)
            values[first:first+4] = out
    return tuple(values)


def packed_emission(values, phase):
    coefficients = packed_transform(values, 0, True)
    return sum((packed_transform(coefficients, LABELS[phase][first] ^ BETA)
                for first in (0, 8)), ())


def packed_injection(values, phase):
    parts = [packed_transform(values[first:first+8], LABELS[phase][first] ^ BETA,
                              transposed=True) for first in (0, 8)]
    combined = tuple(a ^ b for a, b in zip(*parts))
    return packed_transform(combined, 0, True, True)


def check_packed():
    checks = 0
    offsets = {0} | {LABELS[p][first] ^ BETA for p in range(4) for first in (0, 8)}
    for offset in sorted(offsets):
        for inverse in (False, True):
            for transposed in (False, True):
                for j in range(8):
                    for value in (1, 2, 4, 8, 16, 32, 64, 128):
                        x = tuple(value if i == j else 0 for i in range(8))
                        assert packed_transform(x, offset, inverse, transposed) == evaluate(
                            fft_ops(8, offset, inverse, transposed), x)
                        checks += 1
    rng = Random(816256)
    for phase in range(4):
        a, c = matrices(phase, True)
        for _ in range(32):
            x = tuple(rng.randrange(256) for _ in range(16))
            s = tuple(rng.randrange(256) for _ in range(8))
            d = tuple(rng.randrange(1, 256) for _ in range(8))
            refreshed = tuple(MUL[v][diagonal] for v, diagonal in zip(s, d))
            assert packed_emission(s, phase) == matvec(a, s)
            assert packed_injection(x, phase) == matvec(c, x)
            forward_y = tuple(v ^ MUL[2][e] for v, e in zip(x, packed_emission(refreshed, phase)))
            forward_s = tuple(v ^ i for v, i in zip(refreshed, packed_injection(x, phase)))
            dual_y = tuple(rng.randrange(256) for _ in range(16))
            dual_s = tuple(rng.randrange(256) for _ in range(8))
            back_x = tuple(v ^ e for v, e in zip(dual_y, packed_emission(dual_s, phase)))
            back_s = tuple(MUL[diagonal][v ^ MUL[2][i]] for diagonal, v, i in
                           zip(d, dual_s, packed_injection(dual_y, phase)))
            dot = lambda lhs, rhs: xor_sum(MUL[a][b] for a, b in zip(lhs, rhs))
            assert dot(forward_y, dual_y) ^ dot(forward_s, dual_s) == dot(x, back_x) ^ dot(s, back_s)
            checks += 3
        assert packed_cost(phase) == dict(gfni=34, boolean=53, shuffle=46)
    return checks


def packed_cost(phase, coefficient_basis=False):
    out = sum((packed_fft_cost(LABELS[phase][first] ^ BETA) +
               packed_fft_cost(LABELS[phase][first] ^ BETA, transposed=True)
               for first in (0, 8)), Counter())
    out.update(boolean=2)  # sum the two injection coefficient vectors
    if not coefficient_basis:
        out += packed_fft_cost(0, True) + packed_fft_cost(0, True, True)
    # Input D:4, state D:2, eta:2; output/state XOR additions:4+2.
    out.update(gfni=8, boolean=6)
    return dict(out)


def row_pair_audit(a):
    """For two active inputs, repetitions of C-row ratios count zeros."""
    worst = None
    histogram = Counter()
    for i, j in combinations(range(16), 2):
        assert all(a[i]) and all(a[j])
        groups = {}
        for state in range(8):
            ratio = MUL[a[i][state]][INV[a[j][state]]]
            groups.setdefault(ratio, []).append(state)
        largest = max(groups.values(), key=len)
        histogram[len(largest)] += 1
        if worst is None or len(largest) > worst["zero_residues"]:
            ratio = MUL[a[i][largest[0]]][INV[a[j][largest[0]]]]
            image = tuple(a[i][s] ^ MUL[ratio][a[j][s]] for s in range(8))
            h = tuple(s for s, value in enumerate(image) if value)
            return_matrix = [[int(s == j) for j in h]+[a[i][s], a[j][s]]
                             for s in range(8)]
            worst = dict(inputs=(i, j), ratio=ratio, zero_residues=len(largest),
                         zero_state_coordinates=largest, live_state_coordinates=h,
                         return_rank=rank(return_matrix))
    return dict(max_zero_histogram=dict(sorted(histogram.items())), witness=worst)


def column_pair_audit(a):
    """For two live state coordinates, count equal emission ratios."""
    worst = None
    for j, k in combinations(range(8), 2):
        groups = {}
        for i in range(16):
            ratio = MUL[a[i][j]][INV[a[i][k]]]
            groups.setdefault(ratio, []).append(i)
        for ratio, rows in groups.items():
            if worst is None or len(rows) > worst["output_zeros"]:
                worst = dict(state_coordinates=(j, k), ratio=ratio,
                             output_zeros=len(rows), zero_output_coordinates=rows)
    return worst


def check():
    reports = []
    checks = 0
    for phase in range(4):
        original, c = matrices(phase, True)
        for j in range(8):
            value = tuple(int(i == j) for i in range(8))
            assert emission(value, phase) == matvec(original, value)
            checks += 1
        for j in range(16):
            value = tuple(int(i == j) for i in range(16))
            assert injection(value, phase) == matvec(c, value)
            checks += 1
        a = coefficient_matrix(phase)
        assert all(all(row) for row in a)
        assert all(xor_sum(MUL[a[i][j]][a[i][k]] for i in range(16)) == 0
                   for j in range(8) for k in range(8))
        assert rank(a) == 8
        # Exhaust all 12870 maximal minors: MDS remains true despite the
        # deficient small coordinate minors. This is not superregularity.
        for selected in combinations(range(16), 8):
            assert rank([a[i] for i in selected]) == 8
            checks += 1
        # A true state-coordinate change uses C'=M^-1 C, rather than A'^T.
        # For this normalized basis, its extra Gram factor is binary:
        # G_ij=1[i&j=0]=[[1,1],[1,0]] tensor-cubed.
        dual_injection_columns = [tuple(xor_sum(row[j] for j in range(8) if not i & j)
                                          for i in range(8)) for row in a]
        assert all(all(row) for row in dual_injection_columns)
        assert all(xor_sum(MUL[a[u][j]][dual_injection_columns[v][j]] for j in range(8)) ==
                   xor_sum(MUL[original[u][j]][original[v][j]] for j in range(8))
                   for u in range(16) for v in range(16))
        reports.append(dict(phase=phase, two_input=row_pair_audit(a),
                            true_basis_two_input=row_pair_audit(dual_injection_columns),
                            two_state=column_pair_audit(a),
                            original_fft_cost=packed_cost(phase),
                            coefficient_basis_cost=packed_cost(phase, True)))
    # Systematic state basis in phase zero: direct first-half injection.
    systematic = []
    for i in range(16):
        row = []
        for j in range(8):
            value = tuple(int(k == j) for k in range(8))
            coefficients = evaluate(fft_ops(8, 0, True), value)
            target = evaluate(fft_ops(8, LABELS[0][8]), coefficients)
            row.append(int(i == j) if i < 8 else target[i-8])
        systematic.append(row)
    assert all(xor_sum(MUL[systematic[i][j]][systematic[i][k]] for i in range(16)) == 0
               for j in range(8) for k in range(8))
    assert sum(bool(x) for x in systematic[0]) == 1
    print(json.dumps(dict(exact_checks=checks, packed_checks=check_packed(), phases=reports,
                          systematic_singleton_return_rank=1), indent=2))


if __name__ == "__main__":
    check()
