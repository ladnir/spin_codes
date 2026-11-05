"""Small exact controls for the syndrome rank lemma; no production gate."""
from collections import Counter
from fractions import Fraction
from itertools import product
import unittest

from hd_two_sum_syndrome_rank import (all_subset_rank_failure,
    fixed_arity_rank_failure, full_rank_fiber_dimension, projective_classes)


def mul(x, y):
    """GF4 in the polynomial basis modulo X^2+X+1."""
    answer = 0
    while y:
        if y & 1:
            answer ^= x
        y >>= 1
        x <<= 1
        if x & 4:
            x ^= 7
    return answer


def add(values):
    result = 0
    for value in values:
        result ^= value
    return result


def image(message, diagonals):
    # Four RS[2,1] repetition rows, with all four distinct GF4 labels.
    return tuple(add(mul(mul(message[i], diagonals[i][j]),
                         1 if moment == 0 else i) for i in range(4))
                 for j in range(2) for moment in range(2))


def mul8(x, y):
    result = 0
    while y:
        if y & 1:
            result ^= x
        y >>= 1
        x <<= 1
        if x & 8:
            x ^= 11
    return result


def inverse8(matrix):
    size = len(matrix)
    work = [list(row)+[int(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(i for i in range(column, size) if work[i][column])
        work[pivot], work[column] = work[column], work[pivot]
        inv = next(x for x in range(1, 8) if mul8(x, work[column][column]) == 1)
        work[column] = [mul8(x, inv) for x in work[column]]
        for i in range(size):
            if i != column and work[i][column]:
                scalar = work[i][column]
                work[i] = [x ^ mul8(scalar, y) for x, y in zip(work[i], work[column])]
    return [row[size:] for row in work]


class RankLemmaTests(unittest.TestCase):
    def test_projective_candidates_exact(self):
        for arity in range(1, 5):
            counts = Counter()
            for coefficients in product(range(4), repeat=4):
                if not any(coefficients):
                    continue
                zeros = sum(all(coefficients[2*j] ^ mul(i, coefficients[2*j+1]) == 0
                                for j in range(2)) for i in range(arity))
                self.assertLessEqual(zeros, 1)
                counts[zeros] += 1
            self.assertEqual(tuple(counts[i]//3 for i in (0, 1)),
                             projective_classes(arity, q=4, length=2))

    def test_shortened_dual_atom_exact(self):
        # C=RS[3,2], so C-perp is spanned by (1,3,2) over GF4.
        points = (0, 1, 2)
        code = {tuple(a ^ mul(b, x) for x in points)
                for a, b in product(range(4), repeat=2)}
        dual = {v for v in product(range(4), repeat=3)
                if all(add(mul(x, y) for x, y in zip(v, c)) == 0 for c in code)}
        self.assertEqual(len(dual), 4)
        for target in product(range(4), repeat=3):
            if not any(target):
                continue
            hits = sum(tuple(mul(d, x) for d, x in zip(diagonal, target)) in dual
                       for diagonal in product(range(1, 4), repeat=3))
            self.assertLessEqual(Fraction(hits, 27), Fraction(1, 9))
            if sum(x != 0 for x in target) <= 2:
                self.assertEqual(hits, 0)

    def test_full_rank_and_nuisance_fibers(self):
        full = ((1, 1), (1, 1), (1, 2), (1, 2))
        full_fibers = Counter(image(m, full) for m in product(range(4), repeat=4))
        self.assertEqual(len(full_fibers), 4**4)
        self.assertEqual(set(full_fibers.values()), {1})
        deficient = ((1, 1),)*4
        fibers = Counter(image(m, deficient) for m in product(range(4), repeat=4))
        self.assertEqual(len(fibers), 4**2)
        self.assertEqual(set(fibers.values()), {4**2})
        self.assertEqual(fibers[(0, 0, 0, 0)], 16)

    def test_column_gauge_is_bijective_on_syndromes(self):
        original = ((1, 1), (1, 1), (1, 2), (1, 2))
        scales = (2, 3)
        gauged = tuple(tuple(mul(d[j], scales[j]) for j in range(2)) for d in original)
        for message in product(range(4), repeat=4):
            expected = tuple(mul(value, scales[index//2])
                             for index, value in enumerate(image(message, original)))
            self.assertEqual(image(message, gauged), expected)

    def test_exact_union_and_fiber_scope(self):
        result = all_subset_rank_failure()
        self.assertEqual(result['raw_union'], 8192*sum(
            (count*bound for _, count, bound in result['per_arity']), Fraction(0)))
        self.assertLess(result['union_upper'], Fraction(1, 2**46))
        self.assertGreater(result['union_upper'], Fraction(1, 2**47))
        self.assertFalse(result['distance_certificate'])
        self.assertEqual(full_rank_fiber_dimension(8), 32)
        self.assertEqual(full_rank_fiber_dimension(9), 40)
        self.assertGreaterEqual(fixed_arity_rank_failure(1), 1)

    def test_whole_row_marginal_and_holder_exact(self):
        # Six RS[2,1] rows over GF8. The column diagonals are1,u_i^2,
        # with labels1..6, so the syndrome columns are (1,u,u^2,u^3).
        # Any four are independent, hence every five-row complement
        # is surjective. A six-row fiber has64 words and each whole-row
        # marginal is uniform. This is exactly the a-5=1 interface.
        from itertools import combinations
        labels = tuple(range(1, 7))
        columns = [(1, u, mul8(u, u), mul8(mul8(u, u), u)) for u in labels]
        for chosen in combinations(range(6), 4):
            inverse8([[columns[i][j] for i in chosen] for j in range(4)])
        inverse = inverse8([[columns[i][j] for i in range(4)] for j in range(4)])
        for syndrome in ((0, 0, 0, 0), (1, 0, 0, 0)):
            fiber = []
            for tail in product(range(8), repeat=2):
                residual = [syndrome[j] ^ mul8(columns[4][j], tail[0])
                            ^ mul8(columns[5][j], tail[1]) for j in range(4)]
                prefix = tuple(add(mul8(x, y) for x, y in zip(row, residual))
                               for row in inverse)
                word = prefix+tail
                self.assertEqual(tuple(add(mul8(columns[i][j], word[i]) for i in range(6))
                                       for j in range(4)), syndrome)
                fiber.append(word)
            self.assertEqual(len(set(fiber)), 64)
            for row in range(6):
                self.assertEqual(Counter(word[row] for word in fiber), Counter({x: 8 for x in range(8)}))
            # Include each nonzero-row indicator INSIDE its row function.
            weights = [[0]+[1+(i+2)*x % 5 for x in range(1, 8)] for i in range(6)]
            left = 0
            for word in fiber:
                term = 1
                for i, value in enumerate(word):
                    term *= weights[i][value]
                left += term
            right_power = Fraction(64**6)
            for row in weights:
                right_power *= Fraction(sum(x**6 for x in row), 8)
            self.assertLessEqual(left**6, right_power)

    def test_guards(self):
        for kwargs in ({'q': 6}, {'length': 257}, {'dimension': 16}):
            with self.assertRaises(ValueError):
                fixed_arity_rank_failure(8, **kwargs)
        with self.assertRaises(ValueError):
            all_subset_rank_failure(row_labels=257)
        with self.assertRaises(ValueError):
            full_rank_fiber_dimension(3)


if __name__ == '__main__':
    unittest.main()
