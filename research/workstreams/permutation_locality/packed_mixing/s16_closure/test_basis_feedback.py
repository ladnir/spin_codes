"""Exact algebra and census tests; no encoder performance measurements."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
import unittest

from flint import arb, ctx
import basis_feedback as basis


class BasisFeedbackTests(unittest.TestCase):
    def test_inverse_for_all_gl3_matrices(self):
        accepted = 0
        for rows in product(range(8), repeat=3):
            if basis.s16_maps.binary_rank(rows) != 3:
                with self.assertRaises(ValueError):
                    basis.invert(rows)
                continue
            inverse = basis.invert(rows)
            accepted += 1
            for value in range(8):
                self.assertEqual(basis.apply_rows(inverse, basis.apply_rows(rows, value)), value)
                self.assertEqual(basis.apply_rows(rows, basis.apply_rows(inverse, value)), value)
        self.assertEqual(accepted, 168)

    def test_reject_malformed_basis_and_columns(self):
        for rows in ((), (True,), (-1,), (2,), (1, 1)):
            with self.assertRaises(ValueError):
                basis.invert(rows)
        for columns in ((), (-1,), (8,), (True,)):
            with self.assertRaises(ValueError):
                basis.transform(columns, (1, 2, 4))

    def test_seed_repeatability(self):
        for seed in range(16):
            rows, attempts = basis.seeded_basis(seed)
            self.assertEqual((rows, attempts), basis.seeded_basis(seed))
            self.assertEqual(len(rows), 16)
            self.assertEqual(basis.s16_maps.binary_rank(rows), 16)
            self.assertGreaterEqual(attempts, 1)

    def test_all_inputs_preserve_tiny_kernel_and_code(self):
        columns = (1, 2, 4, 3, 5, 7, 6, 1)
        old_code = set(basis.s16_maps.images_from_rows(basis.feedback_rows(columns, 3)))
        for seed in range(16):
            rows, _ = basis.seeded_basis(seed, 3)
            changed, inverse = basis.transform(columns, rows)
            code = set(basis.s16_maps.images_from_rows(basis.feedback_rows(changed, 3)))
            self.assertEqual(code, old_code)
            for x in range(256):
                old = new = 0
                for j in range(8):
                    if x >> j & 1:
                        old ^= columns[j]
                        new ^= changed[j]
                self.assertEqual(new == 0, old == 0)
                self.assertEqual(basis.apply_rows(inverse, new), old)

    def test_census_matches_direct_input_enumeration(self):
        images = basis.s16_maps.images_from_rows((1, 6, 120))
        columns = (1, 2, 4, 3, 5, 7, 6, 1)
        rows, _ = basis.seeded_basis(4, 3)
        changed, _ = basis.transform(columns, rows)
        histogram, joint, denominator = basis.birth_census(images, changed)
        expected, expected_joint = Counter(), Counter()
        for x in range(1, 256):
            if bool(x & 15) + bool(x >> 4) != 1:
                continue
            state = 0
            for j, column in enumerate(changed):
                if x >> j & 1:
                    state ^= column
            w = images[state].bit_count()
            expected[w] += 1
            expected_joint[(x.bit_count(), w)] += 1
        self.assertEqual(histogram, dict(expected))
        self.assertEqual(joint, dict(expected_joint))
        self.assertEqual(denominator, 30)
        self.assertEqual(basis.exact_moment(histogram, denominator, Q(1, 2)),
                         sum(Q(n, 30)*Q(1, 2)**w for w, n in expected.items()))
        self.assertEqual(basis.exact_moment(histogram, denominator, Q(1)), 1)

    def test_outward_moment_contains_independent_exponential_sum(self):
        ctx.prec = 192
        histogram = {48: 19, 56: 200, 64: 261}
        value = basis.tilted_moment(histogram, 480)
        reference = sum((n*(-arb(27*w)/100).exp() for w, n in histogram.items()), arb(0))/480
        self.assertTrue(value.overlaps(reference))
        self.assertTrue(value > 0)

    def test_all_production_candidates_preserve_full_row_code(self):
        for feedback in ('weight5', 'bch16'):
            images, columns, _ = basis.s16_maps.candidate(feedback, 0)
            original_images = tuple(images)
            old_code = set(basis.s16_maps.images_from_rows(basis.feedback_rows(columns, 16)))
            for seed in range(16):
                rows, _ = basis.seeded_basis(seed)
                changed, _ = basis.transform(columns, rows)
                new_code = set(basis.s16_maps.images_from_rows(basis.feedback_rows(changed, 16)))
                self.assertEqual(new_code, old_code)
                histogram, _, denominator = basis.birth_census(images, changed)
                self.assertEqual(denominator, 480)
                self.assertEqual(sum(histogram.values()), 480)
            self.assertEqual(tuple(images), original_images)

    def test_report_scope_rank_and_precision_restoration(self):
        ctx.prec = 160
        report = basis.evaluate('bch16', seeds=(0, 1), precision=192)
        self.assertEqual(ctx.prec, 160)
        self.assertFalse(report['whole_code_certificate'])
        self.assertEqual(len(report['candidates']), 3)
        names = {row['name'] for row in report['candidates']}
        self.assertEqual(set(report['rank_by_moment']), names)
        self.assertEqual(set(report['rank_by_bad48_then_moment']), names)
        original = report['candidates'][0]
        self.assertEqual(original['name'], 'identity')
        self.assertEqual(original['feedback_columns'], report['original_feedback_columns'])
        self.assertTrue(all(row['kernel_and_row_code_preserved'] for row in report['candidates']))

    def test_selected_explicit_matrices_and_birth_counts(self):
        rows, attempts = basis.seeded_basis(0)
        self.assertEqual(rows, (0xe8e5, 0x8856, 0xfb97, 0xb486, 0xcf6a, 0x9a16,
            0xe6f4, 0x259f, 0x4f65, 0x1948, 0xbad6, 0x12e0, 0xe61a, 0xd9b8,
            0xaf19, 0x5487))
        self.assertEqual(attempts, 3)
        images, columns, _ = basis.s16_maps.candidate('weight5', 0)
        changed, _ = basis.transform(columns, rows)
        self.assertEqual(basis.birth_census(images, columns)[0],
                         {48: 3, 56: 124, 64: 242, 72: 106, 80: 5})
        self.assertEqual(basis.birth_census(images, changed)[0],
                         {56: 87, 64: 263, 72: 123, 80: 7})
        images, columns, _ = basis.s16_maps.candidate('bch16', 0)
        rows, _ = basis.seeded_basis(11)
        changed, _ = basis.transform(columns, rows)
        self.assertEqual(basis.birth_census(images, columns)[0],
                         {48: 4, 56: 74, 64: 275, 72: 122, 80: 5})
        self.assertEqual(basis.birth_census(images, changed)[0],
                         {48: 2, 56: 110, 64: 263, 72: 100, 80: 5})


if __name__ == '__main__':
    unittest.main()
