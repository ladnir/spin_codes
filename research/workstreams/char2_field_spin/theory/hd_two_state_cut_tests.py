"""Independent combinatorial and GF4 checks of the mixed-interval cut."""
from fractions import Fraction as F
from itertools import combinations, product
import random
import unittest

from hd_two_state_cut import (direct_window_bound, finite_bound, interval_policy, interval_statistics,
                              regional_moments, word_moments)
from hd_two_state_cut_screen import stationary_statistics
from hd_rational_cut_tests import (inverse_matrix, inverse_scalar, mul, multiply,
                                   prefix_matrix, rank)


def mask_intervals(mask):
    result = []
    i = 0
    while i < len(mask):
        if not mask[i]:
            i += 1
            continue
        first = i
        while i < len(mask) and mask[i]:
            i += 1
        result.append((first, i))
    return result


def two_state_matrix(tuples):
    n = len(tuples)
    matrix = [[0]*n for _ in range(n)]
    for column in range(n):
        s1 = s2 = 0
        for i, (c, d1, d2, f) in enumerate(tuples):
            x = int(i == column)
            matrix[i][column] = s1 ^ mul(f, s2) ^ mul(c, x)
            s1 ^= mul(d1, x)
            s2 ^= mul(d2, x)
    return matrix


def null_vector(matrix):
    """Independent reduced elimination; return one nonzero kernel vector."""
    reduced = [list(row) for row in matrix]
    columns = len(reduced[0])
    pivots = []
    for column in range(columns):
        lead = len(pivots)
        pivot = next((i for i in range(lead, len(reduced)) if reduced[i][column]), None)
        if pivot is None:
            continue
        reduced[lead], reduced[pivot] = reduced[pivot], reduced[lead]
        scale = inverse_scalar(reduced[lead][column])
        reduced[lead] = [mul(scale, value) for value in reduced[lead]]
        for i in range(len(reduced)):
            if i != lead:
                factor = reduced[i][column]
                reduced[i] = [x ^ mul(factor, y) for x, y in zip(reduced[i], reduced[lead])]
        pivots.append(column)
    free = next(i for i in range(columns) if i not in pivots)
    vector = [int(i == free) for i in range(columns)]
    for i, pivot in enumerate(pivots):
        vector[pivot] = reduced[i][free]
    return vector


class TwoStateCutTests(unittest.TestCase):
    def test_policy_against_all_interval_unions(self):
        for n in range(1, 8):
            all_intervals = [mask_intervals(mask) for mask in product((0, 1), repeat=n)]
            for bits in product((0, 1), repeat=n):
                g, c = interval_statistics(bits, interval_policy(bits))
                maximum = max(3*gain-cost for gain, cost in
                              (interval_statistics(bits, intervals)
                               for intervals in all_intervals))
                self.assertEqual(3*g-c, maximum)

    def test_word_moments_against_enumeration(self):
        for n in range(1, 11):
            expected = [[0]*5 for _ in range(n+1)]
            for bits in product((0, 1), repeat=n):
                g, c = interval_statistics(bits, interval_policy(bits))
                for i, value in enumerate((1, g, c, g*g, c*c)):
                    expected[sum(bits)][i] += value
            self.assertEqual(word_moments(n), expected)

    def test_without_replacement_moments_against_enumeration(self):
        for rows, length in ((6, 2), (8, 2), (8, 4), (12, 4)):
            sums = word_moments(length)
            for safe in range(rows+1):
                values = []
                for positions in combinations(range(rows), safe):
                    bits = [0]*rows
                    for position in positions:
                        bits[position] = 1
                    gain = cost = 0
                    for start in range(0, rows, length):
                        block = bits[start:start+length]
                        g, c = interval_statistics(block, interval_policy(block))
                        gain += g
                        cost += c
                    values.append((gain, cost))
                for i in (0, 1):
                    mean = F(sum(value[i] for value in values), len(values))
                    variance = sum((value[i]-mean)**2 for value in values)/len(values)
                    self.assertEqual(regional_moments(rows, safe, length, sums, i+1),
                                     (mean, variance))

    def test_gf4_mixed_interval_rank_for_arbitrary_safe_columns(self):
        tuples = [(c, d1, d2, f) for c, d1, d2, f in
                  product((1, 2, 3), repeat=4) if c ^ d1 ^ mul(f, d2)]
        self.assertEqual(len(tuples), 63)
        n, h, checks = 12, 1, 8
        for seed in range(32):
            rng = random.Random(seed)
            matrix = two_state_matrix([rng.choice(tuples) for _ in range(n)])
            inverse = inverse_matrix(matrix)
            self.assertEqual(multiply(matrix, inverse),
                             [[int(i == j) for j in range(n)] for i in range(n)])
            # These masks include safe runs interrupted by unsafe holes.
            for bits in ([1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1],
                         [rng.randrange(2) for _ in range(n)]):
                m = [[rng.randrange(4) if i < h or not bits[j] else 0
                      for j in range(n)] for i in range(checks)]
                final = multiply(m, inverse)
                families = [[(0, 7)], [(0, 4), (5, 12)], [(0, 12)],
                            interval_policy(bits)]
                for intervals in families:
                    columns = [j for first, last in intervals for j in range(first, last)]
                    g, c = interval_statistics(bits, intervals)
                    selected = [[row[j] for j in columns] for row in final]
                    self.assertLessEqual(rank(selected), h+c)
                    self.assertLessEqual(rank(selected[h:]), c)
                    self.assertEqual(len(columns), g+c)
                    if g > h:
                        short = [row[:h+c+1] for row in selected]
                        self.assertLess(rank(short), len(short[0]))
                for first in range(n):
                    for last in range(first+1, n+1):
                        self.assertLessEqual(rank([row[first:last] for row in inverse[last:]]), 2)

    def test_k20_exact_failure_bound(self):
        result = finite_bound()
        self.assertEqual(result["safe"], 8190)
        self.assertEqual(result["check_capacity"], 262144)
        self.assertEqual(result["forbidden_cutoff"], 838860)
        self.assertEqual(result["cost_cap"], 576715)
        self.assertGreater(result["failure_lower"], F(499, 500))
        self.assertGreater(result["variance_gain"], 0)
        self.assertGreater(result["variance_cost"], 0)

    def test_full_striped_precode_route_postscan_and_codeword(self):
        b, k, rows, a = 4, 2, 16, 10
        n, safe, h = b*rows, a-2, k*a
        rng = random.Random(73)
        # A GF4 [4,2] MDS parity check, repeated over striped outer rows.
        check_columns = [(1, 0), (0, 1), (1, 1), (1, 2)]
        h0 = [[0]*n for _ in range(k*rows)]
        for i in range(rows):
            for j in range(b):
                for check in range(k):
                    h0[k*i+check][j*rows+i] = check_columns[j][check]
        q0 = prefix_matrix([rng.randrange(1, 4) for _ in range(n)])
        q1 = prefix_matrix([rng.randrange(1, 4) for _ in range(n)])
        before_route = multiply(h0, inverse_matrix(multiply(q1, q0)))
        permutation, intervals = [], []
        for region in range(b):
            # A realized regional route with one unsafe hole in each of
            # the first two intervals. No state resets at region edges.
            positions = list(range(safe+int(region < 2)))
            if region < 2:
                positions.remove(4)
            safe_values = iter(range(safe))
            unsafe_values = iter(range(safe, rows))
            permutation.extend(region*rows+next(safe_values if j in positions else unsafe_values)
                               for j in range(rows))
            intervals.append((region*rows, region*rows+safe+int(region < 2)))
        m = [[row[i] for i in permutation] for row in before_route]
        bits = [int(i % rows < safe) for i in permutation]
        for j, bit in enumerate(bits):
            if bit:
                self.assertTrue(all(row[j] == 0 for row in m[h:]))
        tuples = [(c, d1, d2, f) for c, d1, d2, f in
                  product((1, 2, 3), repeat=4) if c ^ d1 ^ mul(f, d2)]
        post = two_state_matrix([rng.choice(tuples) for _ in range(n)])
        final = multiply(m, inverse_matrix(post))
        g, c = interval_statistics(bits, intervals)
        self.assertEqual((g, c), (24, 10))
        columns = [j for first, last in intervals for j in range(first, last)]
        selected = [[row[j] for j in columns] for row in final]
        self.assertLessEqual(rank(selected), h+c)
        short_columns = columns[:h+c+1]
        short = [[row[j] for j in short_columns] for row in final]
        witness = null_vector(short)
        output = [[0] for _ in range(n)]
        for j, value in zip(short_columns, witness):
            output[j][0] = value
        self.assertTrue(any(value for value in witness))
        self.assertLessEqual(sum(bool(row[0]) for row in output), 31)
        self.assertEqual(multiply(final, output), [[0] for _ in final])

    def test_direct_window_bound_and_row_dependence(self):
        # Four independently routed coordinates, each choosing one of
        # three rows. Row occupancies are dependent but each is binomial.
        dimensions = []
        for chosen_rows in product(range(3), repeat=4):
            dimensions.append(sum(max(0, chosen_rows.count(i)-2) for i in range(3)))
        small = direct_window_bound(1, rows=3, b=4, k=2)
        self.assertEqual(small["mean_shortened_dimension"], F(sum(dimensions), len(dimensions)))
        for states in (3, 4, 8):
            result = direct_window_bound(states)
            self.assertEqual(result["width"], 13107)
            self.assertEqual(result["output_weight_bound"], 838848)
            self.assertGreater(result["failure_lower"], F(1, 512))
        self.assertGreater(direct_window_bound(44)["failure_lower"], F(1, 1 << 40))
        self.assertEqual(direct_window_bound(45)["failure_lower"], 0)

    def test_stationary_screen_known_small_chain(self):
        result = stationary_statistics(0.5, 1, 1, 2)
        self.assertAlmostEqual(result["gain"], 13/50, places=12)
        self.assertAlmostEqual(result["cost"], 4/25, places=12)
        self.assertAlmostEqual(result["score"], 1/10, places=12)


if __name__ == "__main__":
    unittest.main()
