"""Exact small-field tests for the conservative packet operator."""

from fractions import Fraction
from itertools import combinations, product
from random import Random
import unittest

import numpy as np

from rank_packet_operator import (coefficient_table, dependent_action,
                                  dependent_updates, packet_action, packet_matrices,
                                  packet_action_exact, rank_terms,
                                  root_intersection_probability, scaled_rank_terms)
from regional_sparse_cover import gf4_mul, gf4_rank


def add_scaled(target, scalar, source):
    for i, x in enumerate(source):
        target[i] ^= gf4_mul(scalar, x)


def quotient_rank(rows, equations):
    return rank(equations + rows) - rank(equations)


def rank(rows):
    return gf4_rank(rows) if rows else 0


def gf16_mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        if a & 16:
            a ^= 0x13
        b >>= 1
    return result


MUL16 = [[gf16_mul(a, b) for b in range(16)] for a in range(16)]
INV16 = [0] + [next(b for b in range(1, 16) if MUL16[a][b] == 1) for a in range(1, 16)]


def rank16(rows):
    if not rows:
        return 0
    work = [row[:] for row in rows]
    rank = 0
    for col in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][col]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inv = INV16[work[rank][col]]
        work[rank] = [MUL16[inv][a] for a in work[rank]]
        for i in range(rank+1, len(work)):
            value = work[i][col]
            if value:
                work[i] = [a ^ MUL16[value][b] for a, b in zip(work[i], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


class RankPacketTests(unittest.TestCase):
    def test_coefficient_multiplicity_and_normalization(self):
        # Empty input/zero state has only the empty output and remains zero.
        table = coefficient_table(4, 2)
        self.assertEqual(table[0, 0, 0, 0, 0, 0, 0], 1)
        raw = packet_matrices(2, .3, q=4, t=4, s=2)
        avg = packet_matrices(2, .3, q=4, t=4, s=2, normalized=True)
        self.assertTrue(np.allclose(raw[1, 1], 12 * avg[1, 1]))
        self.assertEqual(raw[0, 0, 0, 0], 1)

    def test_no_output_return(self):
        for f0 in range(5):
            for d0 in range(5 - f0):
                terms = list(rank_terms(f0, d0, 0, 0, 0, 4-f0-d0, s=2))
                self.assertFalse(any(r and not k for r, k, _, _ in terms))

    def test_dependent_update_rank_law(self):
        rng = Random(5125504)
        for _ in range(512):
            matrix = [[rng.randrange(4) for _ in range(3)] for _ in range(2)]
            form = [0, 0, 0]
            column = [0, 0]
            while not any(form):
                form = [rng.randrange(4) for _ in range(3)]
            while not any(column):
                column = [rng.randrange(4) for _ in range(2)]
            before = rank(matrix)
            actual = [0, 0, 0]
            for diagonal in (1, 2, 3):
                updated = [row[:] for row in matrix]
                for i in range(2):
                    add_scaled(updated[i], gf4_mul(column[i], diagonal), form)
                actual[rank(updated)] += 1
            bounds = [Fraction(0)]*3
            for after, exceptional, count in dependent_updates(before, 1, 2):
                bounds[after] += Fraction(count, 3**exceptional)
            for after in range(3):
                self.assertLessEqual(Fraction(actual[after], 3), bounds[after])
            payoff = [rng.randrange(1, 100) for _ in range(3)]
            exact = sum(Fraction(count, 3)*value for count, value in zip(actual, payoff))
            self.assertLessEqual(float(exact), dependent_action(payoff, 4)[before]+1e-12)

    def test_scaled_zero_overlap_never_returns(self):
        for d0 in range(5):
            for output_only in range(5-d0):
                terms = scaled_rank_terms(0, d0, 0, 0, output_only,
                                          4-d0-output_only, s=2)
                self.assertFalse(any(r and not k for r, k, _, _, _ in terms))

    def test_bellman_is_homogeneous_and_tightens_matrix(self):
        for tau, y, vector in ((2, .2, [1, 2, 5]), (3, .7, [7, 3, 1]),
                               (1, 1, [1, 1, 1])):
            action = packet_action(tau, y, vector, q=4, t=4, s=2)
            triple = packet_action(tau, y, np.asarray(vector)*3, q=4, t=4, s=2)
            matrix = packet_matrices(tau, y, q=4, t=4, s=2, law='eta-scaled')
            self.assertTrue(np.allclose(triple, action*3))
            self.assertTrue(np.all(action <= matrix @ vector + 1e-9))
            reference = packet_action_exact(Fraction(tau), Fraction(str(y)), vector,
                                            q=4, t=4, s=2)
            self.assertTrue(np.allclose(action, np.asarray(reference, dtype=float)))
            rotated = packet_action(tau, y, vector, q=4, t=4, s=2, fresh_orientation=True)
            exact_rotated = packet_action_exact(Fraction(tau), Fraction(str(y)), vector,
                                                q=4, t=4, s=2, fresh_orientation=True)
            self.assertTrue(np.all(rotated <= action+1e-12))
            self.assertTrue(np.allclose(rotated, np.asarray(exact_rotated, dtype=float)))

    def test_grassmann_intersection_exhaustive_gl2(self):
        matrices = []
        for entries in product(range(4), repeat=4):
            matrix = [list(entries[:2]), list(entries[2:])]
            if rank(matrix) == 2:
                matrices.append(matrix)
        self.assertEqual(len(matrices), 180)
        for r in range(3):
            for z in range(3):
                counts = [0, 0, 0]
                for matrix in matrices:
                    projected = [row[:r] for row in matrix[:z]]
                    h = r-(rank(projected) if r and z else 0)
                    counts[h] += 1
                for h in range(3):
                    self.assertEqual(Fraction(counts[h], 180),
                                     root_intersection_probability(4, 2, r, z, h))

    def test_scaled_gf16_packet_with_nonmonotone_continuation(self):
        # Embedded GF4={0,1,6,7}, eta=2 outside it.  This tests the genuine
        # all-F_J-invertible packet law, not an assumed rank process.
        labels, eta = (0, 1, 6, 7), 2
        for mask in range(16):
            matrix = [[1, 0], [0, 1]]
            for j, label in enumerate(labels):
                if mask >> j & 1:
                    c, a = (1, label), (MUL16[eta][label], eta)
                    for row in range(2):
                        for col in range(2):
                            matrix[row][col] ^= MUL16[c[row]][a[col]]
            self.assertEqual(rank16(matrix), 2)
        tau, y, payoff = Fraction(2), Fraction(1, 3), [Fraction(3), Fraction(1), Fraction(7)]
        bound = packet_action_exact(tau, y, payoff, q=16, t=4, s=2)
        for incoming, fresh, dep in ((0, 1, 0), (1, 1, 0), (2, 0, 2)):
            total = Fraction(0)
            for fresh_tuple in combinations(range(4), fresh):
                for dep_tuple in combinations([j for j in range(4) if j not in fresh_tuple], dep):
                    inputs = sorted(fresh_tuple+dep_tuple)
                    state = [[int(i == j and i < incoming) for j in range(3+fresh)]
                             for i in range(2)]
                    forms = {}
                    for number, j in enumerate(fresh_tuple):
                        forms[j] = [int(k == 3+number) for k in range(3+fresh)]
                    for number, j in enumerate(dep_tuple):
                        forms[j] = [1, number+1, number+2]+[0]*fresh
                    for diagonal in product(range(1, 16), repeat=len(inputs)):
                        ds = dict(zip(inputs, diagonal))
                        nxt = [row[:] for row in state]
                        ys = []
                        for j, label in enumerate(labels):
                            output = [MUL16[eta][MUL16[label][a] ^ b]
                                      for a, b in zip(*state)]
                            if j in forms:
                                x = [MUL16[ds[j]][value] for value in forms[j]]
                                output = [a ^ b for a, b in zip(output, x)]
                                nxt[0] = [a ^ b for a, b in zip(nxt[0], x)]
                                nxt[1] = [a ^ MUL16[label][b] for a, b in zip(nxt[1], x)]
                            ys.append(output)
                        for mask in range(16):
                            outputs = [j for j in range(4) if mask >> j & 1]
                            equations = [ys[j] for j in range(4) if j not in outputs]
                            for returning in (False, True):
                                eqs = equations + (nxt if returning else [])
                                er = rank16(eqs)
                                if any(rank16(eqs+[forms[j]]) == er for j in inputs):
                                    continue
                                if any(rank16(eqs+[ys[j]]) == er for j in outputs):
                                    continue
                                if incoming and rank16(eqs+state) == er:
                                    continue
                                k = rank16(eqs+nxt)-er
                                if not returning and k == 0:
                                    continue
                                total += y**len(outputs)*payoff[k]/tau**er/15**len(inputs)
            self.assertLessEqual(total, bound[fresh][dep][incoming])

    def test_return_uses_max_not_sum(self):
        # Fresh one-output input preserves at least h state rank, and the
        # return pivot gives one charge: max(h,1), not h+1.
        terms = list(rank_terms(1, 2, 1, 0, 0, 0, s=2))
        self.assertIn((2, 0, 5, 0), terms)

    def test_exhaustive_diagonals_for_random_common_quotients(self):
        rng = Random(3222164)
        cases = 0
        # GF4, t=4, s=2 is the exact lower-dimensional dual-MDS analogue:
        # c(x)=(1,x), a(x)=(x,1), sum_x c(x)a(x)=0.
        for incoming in range(3):
            for sample in range(32):
                fresh = rng.randrange(5)
                dep = rng.randrange(5-fresh)
                labels = list(range(4))
                rng.shuffle(labels)
                fset, dset = set(labels[:fresh]), set(labels[fresh:fresh+dep])
                inputs = fset | dset
                m = 4 + fresh
                state = [[0]*m for _ in range(2)]
                for i in range(incoming):
                    state[i][i] = 1
                if incoming == 1:
                    scale = rng.randrange(4)
                    state[1] = [gf4_mul(scale, x) for x in state[0]]
                forms = {}
                for j, label in enumerate(sorted(fset)):
                    forms[label] = [0]*m
                    forms[label][4+j] = 1
                for label in dset:
                    form = [0]*m
                    while not any(form):
                        form[:4] = [rng.randrange(4) for _ in range(4)]
                    forms[label] = form
                for mask in range(16):
                    outputs = {i for i in range(4) if mask >> i & 1}
                    counts = (len(fset-outputs), len(dset-outputs),
                              len(fset & outputs), len(dset & outputs),
                              len(outputs-inputs), 4-len(inputs | outputs))
                    terms = [x for x in rank_terms(*counts, s=2) if x[0] == incoming]
                    totals = {tau: [Fraction(0)]*3 for tau in (1, 2, 3)}
                    for ds in product((1, 2, 3), repeat=len(inputs)):
                        diagonals = dict(zip(sorted(inputs), ds))
                        ys = []
                        nxt = [row[:] for row in state]
                        for label in range(4):
                            output = [gf4_mul(label, x) ^ y
                                      for x, y in zip(state[0], state[1])]
                            if label in inputs:
                                value = [gf4_mul(diagonals[label], x) for x in forms[label]]
                                add_scaled(output, 1, value)
                                add_scaled(nxt[0], 1, value)
                                add_scaled(nxt[1], label, value)
                            ys.append(output)
                        equations = [ys[i] for i in range(4) if i not in outputs]
                        for returning in (False, True):
                            eqs = equations + (nxt if returning else [])
                            if any(not quotient_rank([forms[i]], eqs) for i in inputs):
                                continue
                            if any(not quotient_rank([ys[j]], eqs) for j in outputs):
                                continue
                            if incoming and not quotient_rank(state, eqs):
                                continue
                            k = quotient_rank(nxt, eqs)
                            if not returning and k == 0:
                                continue
                            for tau in totals:
                                totals[tau][k] += Fraction(1, tau**rank(eqs))
                    for tau, actual in totals.items():
                        bounds = [Fraction(0)]*3
                        for _, k, charge, exceptional in terms:
                            bounds[k] += Fraction(1, tau**charge * 3**exceptional)
                        for k in range(3):
                            self.assertLessEqual(actual[k] / 3**len(inputs), bounds[k],
                                                 (incoming, counts, tau, k, sample))
                    cases += 1
        self.assertEqual(cases, 1536)


if __name__ == '__main__':
    unittest.main()
