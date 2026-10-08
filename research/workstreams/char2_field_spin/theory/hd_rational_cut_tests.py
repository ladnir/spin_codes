"""Independent small combinatorial and GF4 matrix checks of the rational cut."""
from fractions import Fraction as F
from itertools import combinations, product
import random
import unittest

from hd_rational_cut import exact_means, finite_bound, run_statistics, selected_runs


def mul(a, b):
    """Schoolbook GF4 multiplication modulo X^2+X+1."""
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return out


def inverse_scalar(a):
    if not a:
        raise ValueError("zero inverse")
    return mul(a, a)


def multiply(a, b):
    columns = len(b[0]) if b else 0
    out = [[0]*columns for _ in a]
    for i, row in enumerate(a):
        for j in range(columns):
            for t, value in enumerate(row):
                out[i][j] ^= mul(value, b[t][j])
    return out


def rank(matrix):
    a = [list(row) for row in matrix]
    if not a:
        return 0
    lead = 0
    for column in range(len(a[0])):
        pivot = next((i for i in range(lead, len(a)) if a[i][column]), None)
        if pivot is None:
            continue
        a[lead], a[pivot] = a[pivot], a[lead]
        scale = inverse_scalar(a[lead][column])
        a[lead] = [mul(x, scale) for x in a[lead]]
        for i in range(lead+1, len(a)):
            factor = a[i][column]
            if factor:
                a[i] = [x ^ mul(factor, y) for x, y in zip(a[i], a[lead])]
        lead += 1
        if lead == len(a):
            break
    return lead


def inverse_matrix(matrix):
    n = len(matrix)
    a = [list(row)+[int(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        factor = inverse_scalar(a[j][j])
        a[j] = [mul(factor, x) for x in a[j]]
        for i in range(n):
            if i != j:
                factor = a[i][j]
                a[i] = [x ^ mul(factor, y) for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def rational_matrix(pairs):
    # Direct forward recurrence, not the rank-one inverse formula under test.
    n = len(pairs)
    matrix = [[0]*n for _ in pairs]
    for column in range(n):
        state = 0
        for i, (c, d) in enumerate(pairs):
            value = int(i == column)
            matrix[i][column] = state ^ mul(c, value)
            state ^= mul(d, value)
    return matrix


def prefix_matrix(diagonal):
    return [[diagonal[j] if j <= i else 0 for j in range(len(diagonal))]
            for i in range(len(diagonal))]


class RationalCutTests(unittest.TestCase):
    def test_statistics_and_swap_lipschitz_exhaustive(self):
        for n in range(1, 10):
            for bits in product((0, 1), repeat=n):
                t, runs, size = run_statistics(bits)
                self.assertEqual(t, sum(bits[i] and bits[i+1] for i in range(n-1)))
                direct_runs = sum(bits[i] and bits[i+1] and (i == 0 or not bits[i-1])
                                  for i in range(n-1))
                self.assertEqual(runs, direct_runs)
                self.assertEqual(size, t+runs)
                for i, j in combinations(range(n), 2):
                    swapped = list(bits)
                    swapped[i], swapped[j] = swapped[j], swapped[i]
                    other_t, other_runs, _ = run_statistics(swapped)
                    self.assertLessEqual(abs(other_t-t), 2)
                    self.assertLessEqual(abs(other_runs-runs), 2)

    def test_exact_means_and_variance_small_subsets(self):
        for n in range(2, 10):
            for count in range(n+1):
                values = []
                for positions in combinations(range(n), count):
                    bits = [0]*n
                    for i in positions:
                        bits[i] = 1
                    values.append(run_statistics(bits)[:2])
                means = tuple(F(sum(x[i] for x in values), len(values)) for i in (0, 1))
                self.assertEqual(means, exact_means(1, n, count))
                self.assertEqual(tuple(3*x for x in means), exact_means(3, n, count))
                for i in (0, 1):
                    variance = sum((F(x[i])-means[i])**2 for x in values)/len(values)
                    self.assertLessEqual(variance, 4*n)

    def check_rank_instance(self, pairs, bits, seed):
        n, h, checks = len(bits), 1, 4
        rng = random.Random(seed)
        # Safe columns lie in U=span(first h check coordinates). Outside
        # columns are otherwise arbitrary; this is stronger than one outer.
        m = [[rng.randrange(4) if i < h or not bits[j] else 0
              for j in range(n)] for i in range(checks)]
        j = rational_matrix(pairs)
        inv = inverse_matrix(j)
        self.assertEqual(multiply(j, inv), [[int(i == k) for k in range(n)] for i in range(n)])
        final = multiply(m, inv)
        intervals = selected_runs(bits)
        columns = [i for first, last in intervals for i in range(first, last)]
        selected = [[row[i] for i in columns] for row in final]
        self.assertLessEqual(rank(selected), h+len(intervals))
        self.assertLessEqual(rank(selected[h:]), len(intervals))
        # Every interval's off-diagonal inverse tail has rank at most one.
        for first in range(n):
            for last in range(first+1, n+1):
                tail = [row[first:last] for row in inv[last:]]
                self.assertLessEqual(rank(tail), 1)
        if len(columns) > h+len(intervals):
            short = [row[:h+len(intervals)+1] for row in selected]
            self.assertLess(rank(short), len(short[0]))

    def test_gf4_actual_rank_all_n3_setups_and_supports(self):
        pairs = [(c, d) for c in (1, 2, 3) for d in (1, 2, 3) if c != d]
        for sequence in product(pairs, repeat=3):
            for bits in product((0, 1), repeat=3):
                self.check_rank_instance(sequence, bits, 13)

    def test_gf4_separated_runs_and_dense_tails(self):
        masks = [(1, 1, 0, 1, 1, 0, 1, 1),
                 (1, 1, 1, 0, 1, 1, 1, 0),
                 (0, 1, 1, 1, 1, 1, 1, 0),
                 (1, 1, 1, 1, 1, 1, 1, 1)]
        options = [(c, d) for c in (1, 2, 3) for d in (1, 2, 3) if c != d]
        for seed in range(32):
            rng = random.Random(seed)
            pairs = [rng.choice(options) for _ in range(8)]
            for mask in masks:
                self.check_rank_instance(pairs, mask, seed)

    def test_actual_two_prescan_safe_columns_and_regional_route(self):
        b, k, rows, a = 4, 2, 8, 6
        n, safe, capacity = b*rows, a-2, k*a
        rng = random.Random(912)
        # A full-rank block-local check, stored in striped coordinate order.
        h0 = [[0]*n for _ in range(k*rows)]
        for i in range(rows):
            for check in range(k):
                for j in range(b):
                    h0[k*i+check][j*rows+i] = int(check == j) if j < k else rng.randrange(4)
        pre0 = prefix_matrix([rng.randrange(1, 4) for _ in range(n)])
        pre1 = prefix_matrix([rng.randrange(1, 4) for _ in range(n)])
        before_route = multiply(h0, inverse_matrix(multiply(pre1, pre0)))
        perm = []
        for region in range(b):
            order = list(range(region*rows, (region+1)*rows))
            rng.shuffle(order)
            perm.extend(order)
        m = [[row[i] for i in perm] for row in before_route]
        safe_bits = [int(i % rows < safe) for i in perm]
        for j, is_safe in enumerate(safe_bits):
            if is_safe:
                self.assertTrue(all(row[j] == 0 for row in m[capacity:]))
        options = [(c, d) for c in (1, 2, 3) for d in (1, 2, 3) if c != d]
        final = multiply(m, inverse_matrix(rational_matrix([rng.choice(options) for _ in range(n)])))
        columns, run_count = [], 0
        for region in range(b):
            base = region*rows
            local = selected_runs(safe_bits[base:base+rows])
            run_count += len(local)
            columns.extend(base+i for first, last in local for i in range(first, last))
        self.assertLessEqual(rank([[row[i] for i in columns] for row in final]), capacity+run_count)
        self.assertLessEqual(rank([[row[i] for i in columns] for row in final[capacity:]]), run_count)

    def test_k20_exact_failure_lower_bound(self):
        result = finite_bound()
        self.assertEqual((result["a"], result["safe_per_region"]), (17039, 17037))
        self.assertEqual(result["check_capacity"], 545248)
        self.assertEqual(result["forbidden_cutoff"], 838860)
        self.assertEqual(result["run_cap"], 293611)
        self.assertEqual(result["mean_t"], F(18575509248, 32768))
        self.assertEqual(result["mean_runs"], F(292229911489536, 1073709056))
        self.assertGreater(result["failure_lower"], F(19, 20))
        self.assertLessEqual(result["check_capacity"]+result["run_cap"]+1,
                             result["forbidden_cutoff"])
        self.assertEqual(finite_bound(distance=F(1, 2))["forbidden_cutoff"], (1 << 20)-1)

    def test_reject_invalid_or_nonpositive_gaps(self):
        for args in (dict(rows=1), dict(b=64, k=16), dict(a=1), dict(a=40000),
                     dict(distance=F(0)), dict(a=16384), dict(distance=F(1, 5))):
            with self.assertRaises(ValueError):
                finite_bound(**args)


if __name__ == "__main__":
    unittest.main()
