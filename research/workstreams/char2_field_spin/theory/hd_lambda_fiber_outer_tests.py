"""Tiny exact rank, fibre, and labelled-schedule checks; no proof control."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import prod
import random
import unittest

import hd_lambda_fiber_outer as fiber
from hd_lambda_outer import matrix
from hd_pair_scheduled_outer import column_permutations
from hd_quaternary_outer_relations import multiply_gf4


def rank4(rows, columns):
    rows = [list(row) for row in rows]
    used = 0
    inverses = (0, 1, 3, 2)
    for j in range(columns):
        pivot = next((i for i in range(used, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[used], rows[pivot] = rows[pivot], rows[used]
        scale = inverses[rows[used][j]]
        rows[used] = [multiply_gf4(scale, x) for x in rows[used]]
        for i in range(len(rows)):
            if i != used and rows[i][j]:
                scale = rows[i][j]
                rows[i] = [x ^ multiply_gf4(scale, y) for x, y in zip(rows[i], rows[used])]
        used += 1
    return used


def apply4(h, source, values):
    result = []
    for row in h:
        value = 0
        for j, x in zip(source, values):
            value ^= multiply_gf4(row[j], x)
        result.append(value)
    return tuple(result)


def mask(values, labels):
    return sum(bool(values[i]) << i for i in labels)


def weight(values, labels, atoms):
    return prod(atoms[i] for i in labels if values[i])


class LambdaFiberTests(unittest.TestCase):
    def test_rank_one_shortcut_exhaustive_small_binary_matrices(self):
        for height, width in ((1, 3), (2, 3), (3, 3)):
            for rows in product(range(1 << width), repeat=height):
                got = fiber.binary_rank_one_rank(rows, width)
                for theta in (2, 3):
                    direct = tuple(tuple((row >> j & 1) ^ theta for j in range(width))
                                   for row in rows)
                    self.assertEqual(got, rank4(direct, width))
        self.assertEqual(fiber.binary_rank_one_rank((), 4), 0)
        self.assertEqual(fiber.binary_rank_one_rank((0, 0), 0), 0)

    def test_actual_cross_ranks_and_hull_cache(self):
        h, rng = matrix(2, q=4), random.Random(207)
        for _ in range(90):
            left, right = rng.randrange(65536), rng.randrange(65536)
            rows, columns = tuple(fiber._labels(left)), tuple(fiber._labels(right))
            # rank_cross is bilinear, including when row/column sets overlap.
            direct = tuple(tuple(h[i][j] ^ int(i == j) for j in columns) for i in rows)
            got = fiber.rank_cross_masks(left, right)
            self.assertEqual(got, rank4(direct, len(columns)))
            self.assertEqual(got, fiber.rank_cross_masks(right, left))
        self.assertEqual(fiber.rank_cross_masks(0, 65535), 0)

    def test_exact_tiny_cosets_and_raw_pointwise_domination(self):
        rows, q = 4, 4
        h = tuple(tuple(2 ^ int(i == j) for j in range(rows)) for i in range(rows))
        atoms_list = ((F(0), F(1), F(1, 3), F(3, 4)), (F(2, 3),)*rows, (F(1),)*rows)
        for arity in range(rows+1):
            for source in combinations(range(rows), arity):
                outside = tuple(i for i in range(rows) if i not in source)
                def cross_rank(rmask, cmask):
                    rr, cc = tuple(fiber._labels(rmask)), tuple(fiber._labels(cmask))
                    return rank4(tuple(tuple(h[i][j] for j in cc) for i in rr), len(cc))
                amask, omask = sum(1 << i for i in source), sum(1 << i for i in outside)
                c = cross_rank(omask, amask)
                C, Z = fiber._support_enumerators(source, outside, c, cross_rank, q)
                outputs = [apply4(h, source, values) for values in product(range(q), repeat=arity)]
                image = {tuple(y[i] for i in outside) for y in outputs}
                expected_C = Counter(sum(bool(v) << i for i, v in zip(outside, z)) for z in image)
                expected_Z = Counter(mask(y, source) for y in outputs if not mask(y, outside))
                self.assertEqual(dict(C), dict(expected_C))
                self.assertEqual(dict(Z), dict(expected_Z))
                self.assertEqual(sum(dict(C).values()), q**c)
                self.assertEqual(sum(dict(Z).values()), q**(arity-c))
                cosets = defaultdict(list)
                for y in outputs:
                    cosets[tuple(y[i] for i in outside)].append(y)
                for atoms in atoms_list:
                    zvalue = sum(n*prod(atoms[i] for i in fiber._labels(m)) for m, n in Z)
                    cvalue = sum(n*prod(atoms[i] for i in fiber._labels(m)) for m, n in C)
                    direct = sum(weight(y, range(rows), atoms) for y in outputs)
                    self.assertLessEqual(direct, cvalue*zvalue)
                    for values in cosets.values():
                        self.assertLessEqual(sum(weight(y, source, atoms) for y in values), zvalue)

    def test_production_shape_exact_masses_and_nonnegative_counts(self):
        for source in ((0, 1, 2, 3, 4), tuple(range(8)), (0, 2, 4, 6, 8, 10, 12, 14)):
            for q in (4, 256):
                data = fiber.fiber_counts(source, q=q)
                self.assertEqual(sum(data['outside_counts'].values()), q**data['cross_rank'])
                self.assertEqual(sum(data['zero_fiber_counts'].values()), q**data['fiber_dimension'])
                self.assertEqual(data['outside_mass']*data['zero_fiber_mass'], q**len(source))
                for key in ('outside_counts', 'zero_fiber_counts'):
                    self.assertEqual(data[key][0], 1)
                    self.assertTrue(all(type(n) is int and n > 0 for n in data[key].values()))
                self.assertEqual(fiber.evaluate_fiber_bound(data, (1,)*16), q**len(source))
                self.assertEqual(fiber.evaluate_fiber_bound(data, (0,)*16), 1)

    def test_scheduled_factors_keep_actual_output_labels(self):
        source, column = (8, 0, 5, 2, 14, 3, 11, 7), 9
        before, after = column_permutations(fiber.SCHEDULE, column)
        raw = fiber.fiber_counts(tuple(before[i] for i in source), q=4)
        scheduled = fiber.column_fiber_counts(source, column, q=4)
        self.assertEqual(scheduled['original_source'], source)
        self.assertEqual(scheduled['mapped_source'], tuple(before[i] for i in source))
        for key in ('outside_counts', 'zero_fiber_counts'):
            want = {sum(1 << after[i] for i in fiber._labels(m)): n for m, n in raw[key].items()}
            self.assertEqual(scheduled[key], want)
        self.assertEqual(scheduled['physical_inside_mask'] ^ scheduled['physical_outside_mask'], fiber.FULL)
        atoms = tuple(F(i+1, 17) for i in range(16))
        self.assertEqual(fiber.evaluate_fiber_bound(scheduled, atoms),
                         fiber.evaluate_fiber_bound(raw, tuple(atoms[after[i]] for i in range(16))))
        # Mutating a returned dictionary must not contaminate the count cache.
        raw['zero_fiber_counts'][0] = -1
        self.assertEqual(fiber.fiber_counts(tuple(before[i] for i in source), q=4)['zero_fiber_counts'][0], 1)

    def test_guards_and_scope(self):
        for source in ((1, 1), (-1,), (16,), (True,)):
            with self.assertRaises(ValueError):
                fiber.fiber_counts(source)
        with self.assertRaises(ValueError):
            fiber.fiber_counts((0,), max_subsets=256)
        with self.assertRaises(ValueError):
            fiber.fiber_counts(tuple(range(8)), q=8)
        with self.assertRaises(ValueError):
            fiber.binary_rank_one_rank((8,), 3)
        data = fiber.fiber_counts(tuple(range(8)))
        for atoms in ((2,)*16, (-1,)*16, (float('nan'),)*16):
            with self.assertRaises(ValueError):
                fiber.evaluate_fiber_bound(data, atoms)
        for key in ('original_input_floor_included', 'input_weight_marker_supported',
                    'rs_normalization_included', 'whole_code_certificate'):
            self.assertFalse(fiber.metadata()[key])


if __name__ == '__main__':
    unittest.main()
