"""Exact enumeration tests for the whole-packet joint return identity."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from random import Random
import unittest

from cauchy_joint_return import (
    MUL, packet_forms, return_histogram, joint_return_moment, matrix_rank,
    simple_placement_bound, uniform_return_bound, uniform_return_profiles,
    rank_return_terms, rank_return_moment, coarse_rank_return_moment,
)


def direct_histogram(phase, inputs):
    injection, output = packet_forms(phase, inputs)
    histogram = Counter()
    for tail in product(range(1, 256), repeat=len(inputs)-1):
        x = (1,)+tail
        state = [row[0] for row in injection]
        values = [row[0] for row in output]
        for j in range(1, len(x)):
            state = [a ^ MUL[row[j]][x[j]] for a, row in zip(state, injection)]
            values = [a ^ MUL[row[j]][x[j]] for a, row in zip(values, output)]
        histogram[sum((value != 0) << j for j, value in enumerate(state)),
                  sum(value != 0 for value in values)] += 1
    return dict(histogram)


class JointReturnTests(unittest.TestCase):
    def test_direct_projective_enumeration(self):
        for phase, inputs in ((0, (0, 6)), (1, (0, 6)), (2, (1, 4)), (3, (5, 13)),
                              (0, (0, 9, 15)), (1, (0, 9, 15)),
                              (2, (0, 1, 2)), (3, (0, 1, 2))):
            self.assertEqual(return_histogram(phase, inputs), direct_histogram(phase, inputs))

    def test_every_actual_placement_and_exception(self):
        for w in (2, 3):
            exceptions = Counter()
            total_cases = 0
            for phase in range(4):
                for inputs in combinations(range(16), w):
                    injection, _ = packet_forms(phase, inputs)
                    histogram = return_histogram(phase, inputs)
                    self.assertEqual(sum(histogram.values()), 255**(w-1))
                    self.assertGreaterEqual(min(v for h, v in histogram), 13)
                    for h in {h for h, _ in histogram}:
                        rank = matrix_rank([row for j, row in enumerate(injection) if not h & (1 << j)], w)
                        exponent = h.bit_count()+rank
                        probability = joint_return_moment(phase, inputs, h, Fraction(1))
                        self.assertLessEqual(probability, Fraction(1, 255**exponent))
                        if exponent < 4:
                            self.assertEqual(exponent, 3)
                            self.assertEqual(h.bit_count(), 4-w)
                            self.assertIn(phase, (0, 2))
                            self.assertEqual({v: count for (mask, v), count in histogram.items() if mask == h}, {16: 1})
                            self.assertEqual(probability, Fraction(1, 255**3))
                            exceptions[phase, h] += 1
                    total_cases += 1
            self.assertEqual(total_cases, 480 if w == 2 else 2240)
            self.assertEqual(sum(exceptions.values()), 16)
            self.assertEqual(set(exceptions.values()), {2})

    def test_simple_uniform_placement_bounds(self):
        for w in (1, 2, 3):
            for y in (Fraction(1), Fraction(9997, 10000), Fraction(1, 2), Fraction(1, 10)):
                actual = uniform_return_bound(w, y, average_inputs=True)
                self.assertLessEqual(actual, simple_placement_bound(w, y))
                self.assertLessEqual(actual, uniform_return_bound(w, y))
            self.assertTrue(uniform_return_profiles(w))
        self.assertEqual(uniform_return_bound(2, Fraction(1), average_inputs=True), Fraction(17, 4*255**4))
        self.assertEqual(uniform_return_bound(3, Fraction(1), average_inputs=True), Fraction(1, 255**4))

    def test_one_input_regression(self):
        for phase in range(4):
            for j in range(16):
                self.assertEqual(return_histogram(phase, (j,)), {(15, 16): 1})
                self.assertEqual(joint_return_moment(phase, (j,), 15, Fraction(1, 2)), Fraction(1, 2**16*255**4))
        for tau in (1, 14, 35, 255):
            expected = Fraction(1, 2**16*tau**4)
            self.assertEqual(rank_return_moment(1, Fraction(1, 2), tau), expected)
            self.assertEqual(coarse_rank_return_moment(1, Fraction(1, 2), tau), expected)

    def test_selected_zeros_by_independent_linear_feasibility(self):
        # No projective arrangement is used here.  A functional is forced
        # zero on ker(K) exactly when appending it does not increase rank.
        for phase, inputs in ((0, (0, 6)), (1, (0, 6)), (0, (0, 9, 15)), (2, (0, 1, 2))):
            c, b = packet_forms(phase, inputs)
            w = len(inputs)
            coordinate = [tuple(int(i == j) for i in range(w)) for j in range(w)]
            terms = Counter()
            for h in range(1, 16):
                base = [row for j, row in enumerate(c) if not h & (1 << j)]
                required = coordinate+[row for j, row in enumerate(c) if h & (1 << j)]
                for z in range(4):
                    for subset in combinations(range(16), z):
                        rows = base+[b[j] for j in subset]
                        rank = matrix_rank(rows, w)
                        if rank == w or any(matrix_rank(rows+[form], w) == rank for form in required):
                            continue
                        terms[h.bit_count()+rank, z] += 1
            self.assertEqual(tuple((e, z, count) for (e, z), count in sorted(terms.items())),
                             rank_return_terms(phase, inputs))

    def test_generic_small_field_exposure_identity(self):
        # Independently exhaust the diagonal experiment over GF4.  The
        # equality needs only Y=E X+A D S, S'=D S+C X, not superregularity.
        def mul4(a, b):
            answer = 0
            while b:
                if b & 1:
                    answer ^= a
                b >>= 1
                a <<= 1
                if a & 4:
                    a ^= 7
            return answer

        def dot4(a, b):
            answer = 0
            for x, y in zip(a, b):
                answer ^= mul4(x, y)
            return answer

        c = ((1, 1), (1, 2), (1, 3))
        a = ((1, 1, 1), (1, 2, 3), (2, 1, 3), (3, 2, 1))
        b = tuple(tuple(int(i == j) ^ dot4(row, tuple(c[k][j] for k in range(3))) for j in range(2))
                  for i, row in enumerate(a))
        y = Fraction(2, 3)
        for s in product(range(4), repeat=3):
            if not any(s):
                continue
            h = sum((x != 0) << j for j, x in enumerate(s))
            projective = Fraction(0)
            for ratio in range(1, 4):
                x = (1, ratio)
                target = tuple(dot4(row, x) for row in c)
                if sum((value != 0) << j for j, value in enumerate(target)) == h:
                    projective += y**sum(dot4(row, x) != 0 for row in b)
            expected = projective/3**(2+h.bit_count()-1)
            total = Fraction(0)
            for x in product(range(1, 4), repeat=2):
                target = tuple(dot4(row, x) for row in c)
                for diagonal in product(range(1, 4), repeat=3):
                    state = tuple(mul4(d, value) for d, value in zip(diagonal, s))
                    if state != target:
                        continue
                    output = tuple((x[i] if i < 2 else 0) ^ dot4(row, state) for i, row in enumerate(a))
                    total += y**sum(value != 0 for value in output)
            self.assertEqual(total/3**5, expected)

    def test_small_field_return_rank_pivots(self):
        # Arbitrarily dependent input/state forms over a common GF4 message
        # space.  Exhaust every current nonzero diagonal, rather than sampling
        # the events whose rank charges are multiplied in the proof.
        def mul4(a, b):
            answer = 0
            while b:
                if b & 1:
                    answer ^= a
                b >>= 1
                a <<= 1
                if a & 4:
                    a ^= 7
            return answer

        table = [[mul4(a, b) for b in range(4)] for a in range(4)]

        def rank4(rows, width):
            work = [list(row) for row in rows]
            r = 0
            for column in range(width):
                pivot = next((i for i in range(r, len(work)) if work[i][column]), None)
                if pivot is None:
                    continue
                work[r], work[pivot] = work[pivot], work[r]
                inverse = next(b for b in range(1, 4) if table[work[r][column]][b] == 1)
                work[r] = [table[inverse][x] for x in work[r]]
                for i in range(r+1, len(work)):
                    coefficient = work[i][column]
                    work[i] = [a ^ table[coefficient][b] for a, b in zip(work[i], work[r])]
                r += 1
            return r

        def combine(coefficients, forms):
            result = [0]*len(forms[0])
            for coefficient, form in zip(coefficients, forms):
                result = [a ^ table[coefficient][b] for a, b in zip(result, form)]
            return result

        rng = Random(71043)
        all_forms = [form for form in product(range(4), repeat=3) if any(form)]
        cases = 0
        for w, hmask, zcount in product((1, 2, 3), (1, 3, 7), (0, 1, 2)):
            c = [tuple(rng.randrange(1, 4) for _ in range(w)) for _ in range(3)]
            b = [tuple(rng.randrange(4) for _ in range(w)) for _ in range(zcount)]
            hcoords = [j for j in range(3) if hmask & (1 << j)]
            inputs = [rng.choice(all_forms) for _ in range(w)]
            state = [rng.choice(all_forms) if j in hcoords else (0, 0, 0) for j in range(3)]
            required = inputs+[state[j] for j in hcoords]
            k = [row for j, row in enumerate(c) if j not in hcoords]+b
            charge = len(hcoords)+rank4(k, w)
            histogram = Counter()
            for din in product(range(1, 4), repeat=w):
                x = [tuple(table[d][value] for value in form) for d, form in zip(din, inputs)]
                target = [combine(row, x) for row in c]
                output_zeros = [combine(row, x) for row in b]
                for ds in product(range(1, 4), repeat=len(hcoords)):
                    equations = [row[:] for row in target]+output_zeros
                    for j, diagonal in zip(hcoords, ds):
                        equations[j] = [a ^ table[diagonal][value] for a, value in zip(equations[j], state[j])]
                    r = rank4(equations, 3)
                    if any(rank4(equations+[form], 3) == r for form in required):
                        continue
                    histogram[r] += 1
            denominator = 3**(w+len(hcoords))
            for tau in (1, 2, 3):
                actual = sum(Fraction(count, denominator*tau**r) for r, count in histogram.items())
                self.assertLessEqual(actual, Fraction(1, tau**charge))
                cases += 1
        self.assertEqual(cases, 81)

    def test_guards(self):
        for inputs in ((0, 0), (3, 1), (), (16,)):
            with self.assertRaises(ValueError):
                packet_forms(0, inputs)
        with self.assertRaises(ValueError):
            joint_return_moment(0, (0,), 0, Fraction(1))
        with self.assertRaises(ValueError):
            uniform_return_bound(2, float('nan'))


if __name__ == "__main__":
    unittest.main()
