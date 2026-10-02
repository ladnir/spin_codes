import itertools
import unittest
from collections import Counter
from fractions import Fraction as Q
from math import comb

import local_models as model


def rank(rows):
    basis = {}
    for row in rows:
        while row:
            high = row.bit_length()-1
            if high not in basis:
                basis[high] = row
                break
            row ^= basis[high]
    return len(basis)


def field_mul(a, b, degree, modulus):
    value = 0
    for _ in range(degree):
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & (1 << degree):
            a ^= modulus
    return value


class LocalModels(unittest.TestCase):
    def test_shared_gl_exhaustive_three_columns(self):
        invertible = [rows for rows in itertools.product(range(8), repeat=3) if rank(rows) == 3]
        self.assertEqual(len(invertible), 168)
        for r in range(4):
            counts = Counter()
            for rows in invertible:
                support = 0
                for row in rows[:r]:
                    support |= row
                counts[support.bit_count()] += 1
            actual = tuple(Q(counts[w], 168) for w in range(4))
            self.assertEqual(actual, model.shared_gl(3, r))

    def test_independent_rows_exhaustive(self):
        for h in range(5):
            counts = Counter()
            for rows in itertools.product(range(1, 8), repeat=h):
                union = 0
                for row in rows:
                    union |= row
                counts[union.bit_count()] += 1
            self.assertEqual(tuple(Q(counts[w], 7**h) for w in range(4)),
                             model.independent_rows(3, h))

    def test_shared_rank_monotone(self):
        for width in range(1, 9):
            previous = model.cdf(model.shared_gl(width, 0))
            for r in range(1, min(4, width)+1):
                current = model.cdf(model.shared_gl(width, r))
                self.assertTrue(all(a >= b for a, b in zip(previous, current)))
                previous = current

    def test_exact_means(self):
        for r in range(1, 5):
            self.assertEqual(model.mean(model.shared_gl(8, r)), Q(8*((1 << 8)-(1 << (8-r))), 255))
            self.assertEqual(model.mean(model.independent_rows(8, r)), 8*(1-Q(127, 255)**r))
        self.assertEqual(model.mean(model.full_block(8)), Q(8*15*(1 << 28), (1 << 32)-1))

    def test_two_row_piece_laws(self):
        for h in (1, 2):
            counts = Counter()
            for pieces in itertools.product(range(1, 16), repeat=h):
                union = 0
                for word in pieces:
                    union |= int(bool(word & 3)) | (int(bool(word & 12)) << 1)
                counts[union.bit_count()] += 1
            self.assertEqual(tuple(Q(counts[w], 15**h) for w in range(3)), model.independent_pieces(2, 2, h))
        self.assertEqual(model.independent_pieces(8, 1, 4), model.independent_rows(8, 4))
        self.assertEqual(model.independent_pieces(8, 4, 1), model.full_block(8))

    def test_full_gl_exhaustive_four_bits(self):
        counts = Counter()
        matrices = 0
        for rows in itertools.product(range(16), repeat=4):
            if rank(rows) != 4:
                continue
            matrices += 1
            output = rows[0]
            support = bool(output & 3) + bool(output & 12)
            counts[support] += 1
        self.assertEqual(matrices, 20160)
        self.assertEqual(tuple(Q(counts[w], matrices) for w in range(3)), model.full_block(2, 2))

    def test_full_block_labels_are_conditional_product(self):
        # Four2-bit packets; every exact support has all3^w labelings once.
        labels = {}
        for word in range(1, 1 << 8):
            packets = tuple((word >> (2*i)) & 3 for i in range(4))
            support = tuple(i for i, value in enumerate(packets) if value)
            labels.setdefault(support, Counter())[tuple(packets[i] for i in support)] += 1
        for support, counts in labels.items():
            self.assertEqual(set(counts), set(itertools.product(range(1, 4), repeat=len(support))))
            self.assertEqual(set(counts.values()), {1})
            self.assertEqual(len(counts), 3**len(support))

    def test_field_scalar_transitivity(self):
        for degree, modulus in ((2, 7), (3, 11), (4, 19), (8, 0x11b)):
            q = 1 << degree
            for x in range(1, q):
                self.assertEqual({field_mul(a, x, degree, modulus) for a in range(1, q)}, set(range(1, q)))

    def test_random_partition_exhaustive(self):
        for u in range(9):
            counts = Counter()
            for support in itertools.combinations(range(8), u):
                counts[len({i//2 for i in support})] += 1
            self.assertEqual(tuple(Q(counts[h], comb(8, u)) for h in range(5)),
                             model.random_partition(8, 2, u))

    def test_canonical_profile_really_matters(self):
        # Same4 columns total: two full blocks versus four singleton blocks.
        kernel = model.full_block(2, 2)
        self.assertNotEqual(model.repeated(kernel, 2), model.repeated(kernel, 4))
        self.assertLess(model.mean(model.repeated(kernel, 2)), model.mean(model.repeated(kernel, 4)))

    def test_canonical_and_random_moment_bounds(self):
        local, z = model.full_block(2, 2), Q(2, 3)
        canonical = model.support_moment_kernel(8, 2, local, z)
        random = model.support_moment_kernel(8, 2, local, z, True)
        for u in range(9):
            exact = Q(0)
            for support in itertools.combinations(range(8), u):
                h = len({i//2 for i in support})
                exact += model.moment(local, z)**h
            exact /= comb(8, u)
            self.assertEqual(random[u], exact)
            self.assertLessEqual(exact, canonical[u])
        self.assertTrue(all(a >= b for a, b in zip(random, random[1:])))

    def test_cdf_transfer_is_not_shell_subtraction(self):
        actual_shells = [0, 2, 0, 3, 1]
        actual_cdf = list(itertools.accumulate(actual_shells))
        upper_cdf = [0, 3, 4, 5, 6]
        kernel = [Q(1), Q(3, 4), Q(1, 2), Q(1, 4), Q(1, 8)]
        exact = sum(a*b for a, b in zip(actual_shells, kernel))
        self.assertEqual(model.transfer_cdf_moment(actual_cdf, kernel), exact)
        self.assertGreaterEqual(model.transfer_cdf_moment(upper_cdf, kernel), exact)
        # In particular upper_cdf[3]-upper_cdf[2]=1 does NOT bound shell3=3.
        self.assertLess(upper_cdf[3]-upper_cdf[2], actual_shells[3])
        with self.assertRaises(ValueError):
            model.transfer_cdf_moment(upper_cdf, kernel[::-1])

    def test_two_shear_exact_over_gf4(self):
        for x, y in itertools.product(range(4), repeat=2):
            counts = Counter()
            for a, b in itertools.product(range(1, 4), repeat=2):
                first = x ^ field_mul(a, y, 2, 7)
                second = y ^ field_mul(b, first, 2, 7)
                counts[int(bool(first)) + 2*int(bool(second))] += 1
            self.assertEqual(tuple(Q(counts[i], 9) for i in range(4)),
                             model.two_shear_row(4, bool(x), bool(y)))

    def test_two_shear_masks_and_block_bounds(self):
        d = 255
        for first, second in itertools.product(range(16), repeat=2):
            distribution = model.two_shear_masks(256, first, second)
            self.assertEqual(sum(distribution.values()), 1)
            if not (first or second):
                self.assertEqual(distribution, {(0, 0): Q(1)})
                continue
            self.assertNotIn((0, 0), distribution)
            one = sum(p for (a, b), p in distribution.items() if not (a and b))
            if first and not second:
                self.assertEqual(one, 0)
            elif second and not first:
                self.assertLessEqual(one, Q(1, d))
            else:
                self.assertLessEqual(one, Q(2*d-1, d*d))

    def test_reject_bad_inputs(self):
        for function, args in ((model.shared_gl, (8, True)), (model.shared_gl, (3, 4)),
                               (model.independent_rows, (0, 1)), (model.random_partition, (7, 2, 3)),
                               (model.two_shear_row, (6, True, False))):
            with self.assertRaises(ValueError):
                function(*args)
        with self.assertRaises(ValueError):
            model.transfer_cdf_moment([0, 2, 1], [1, Q(1, 2), Q(1, 4)])
        with self.assertRaises(ValueError):
            model.support_moment_kernel(8, 2, (Q(1), Q(0), Q(0)), Q(1, 2))


if __name__ == '__main__':
    unittest.main()
