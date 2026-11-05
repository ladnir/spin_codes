"""Independent tiny coefficient oracles and synthetic adapter checks; no gates."""

from fractions import Fraction as F
from math import comb
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_field128_gamma_all_pair_certificate as cert
import hd_field128_gamma_profile_certificate as profiles


def product(left, right):
    return tuple(left[2*i]*right[j]+left[2*i+1]*right[2+j] for i in range(2) for j in range(2))


def rational_rows(slots):
    return {j: tuple(F((j+2)*(i+1)+i*i+1, 257+j+3*i) for i in range(4))
            for j in range(2*slots+1)}


def direct_coefficients(table, slots, alpha, old_partition=False):
    """Independent scalar histogram loops, without grouped occupancy geometry."""
    result = {}
    exponent = cert.base.rational(alpha)
    for s in range(slots+1):
        for f in range(slots-s+1):
            for t in range(slots+1):
                for g in range(slots-t+1):
                    count = comb(slots, s)*comb(slots-s, f)*comb(slots, t)*comb(slots-t, g)
                    left, right = table[s+2*f], table[t+2*g]
                    if old_partition and (s > 1 or t > 1):
                        values = product(tuple(cert.base.rational(x)**exponent for x in left),
                                         tuple(cert.base.rational(x)**exponent for x in right))
                    else:
                        values = tuple(cert.base.rational(x)**exponent for x in product(left, right))
                    target = result.setdefault((s+t, f+g), [arb(0) for unused in range(4)])
                    for entry, value in enumerate(values):
                        target[entry] += count*value
    return result


class AllPairCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with ctx.workprec(256):
            # Public-size geometry only, with constant synthetic matrices.
            table = {j: (arb(1),)*4 for j in range(65)}
            multiply = cert.base.matrix_product
            with patch.object(cert.base, 'matrix_product', wraps=multiply) as calls:
                cls.composed = cert.all_pair_coefficients(table, F(1))
                cls.raw_products = calls.call_count
            with patch.object(cert, 'directed_parts', return_value=({}, cls.composed)):
                cls.bank = cert.AliasBank(F(1, 2), F(1), ((F(1, 512), F(1, 512)),),
                                          state=24, grid=64, max_seconds=30)

    def test_complete_geometry_cache_and_exact_alpha1_coefficients(self):
        geometry = cert.histogram_geometry()
        self.assertEqual(len(geometry), 561)
        self.assertEqual(sum(count for s, f, count in geometry), 3**32)
        grouped = cert.grouped_pair_geometry()
        self.assertEqual(len(grouped), 65**2)
        self.assertEqual(sum(count for j, k, rows in grouped for s, f, count in rows), 3**64)
        self.assertIs(grouped, cert.grouped_pair_geometry())
        self.assertEqual(self.raw_products, 65**2)
        self.assertEqual(len(self.composed), 2145)
        self.assertEqual(max(s+2*f for s, f in self.composed), 128)
        with ctx.workprec(256):
            for (s, f), values in self.composed.items():
                expected = arb(2*comb(64, s)*comb(64-s, f))
                self.assertTrue(all(value == expected for value in values))
        self.assertIn((0, 64), self.composed)
        self.assertIn((64, 0), self.composed)

    def test_independent_small_q_scalar_oracle_counts_outside_alpha(self):
        slots, exact = 2, rational_rows(2)
        with ctx.workprec(256):
            table = {j: tuple(cert.base.rational(x) for x in row) for j, row in exact.items()}
            for alpha in (F(1), F(1, 2), F(1, 8)):
                actual = cert.all_pair_coefficients(table, alpha, slots=slots)
                expected = direct_coefficients(exact, slots, alpha)
                for key, row in actual.items():
                    for value, oracle in zip(row, expected[key]):
                        self.assertTrue(value >= oracle.lower())
                        self.assertTrue(value < oracle.upper()*(1+arb(2)**-180))
                order01, order10 = product(exact[0], exact[1]), product(exact[1], exact[0])
                self.assertNotEqual(order01, order10)
                for entry in range(4):
                    ordered = slots*(cert.base.rational(order01[entry])**cert.base.rational(alpha)
                                     +cert.base.rational(order10[entry])**cert.base.rational(alpha))
                    self.assertTrue(actual[1, 0][entry] >= ordered.lower())
                    self.assertTrue(actual[1, 0][entry] < ordered.upper()*(1+arb(2)**-180))

    def test_coefficientwise_old_partition_comparison(self):
        exact = rational_rows(3)
        with ctx.workprec(256):
            table = {j: tuple(cert.base.rational(x) for x in row) for j, row in exact.items()}
            for alpha in (F(1), F(1, 4)):
                actual = cert.all_pair_coefficients(table, alpha, slots=3)
                previous = direct_coefficients(exact, 3, alpha, old_partition=True)
                strict = False
                for key, row in actual.items():
                    for value, old in zip(row, previous[key]):
                        self.assertTrue(value <= old.upper()*(1+arb(2)**-180))
                        strict |= bool(value < old.lower()*(1-arb(2)**-100))
                self.assertEqual(strict, alpha != 1)

    def test_factory_full_scope_zero_singles_and_subclass_adapter(self):
        with patch.object(cert.frozen.sparse, 'packet_table', return_value='exact-table') as table, \
             patch.object(cert, 'all_pair_coefficients', return_value=self.composed) as pairs:
            single, composed = cert.directed_parts(F(1, 2), F(1, 4), state=30, bits=96)
        self.assertEqual(single, {})
        self.assertIs(composed, self.composed)
        self.assertEqual(table.call_args.args, (F(1, 2), 64))
        self.assertEqual(table.call_args.kwargs['state'], 30)
        self.assertEqual(table.call_args.kwargs['bits'], 96)
        self.assertEqual(pairs.call_args.args, ('exact-table', F(1, 4)))
        self.assertIsInstance(self.bank, cert.frozen.AliasBank)
        self.assertIs(cert.AliasBank.regional, cert.frozen.AliasBank.regional)
        self.assertIs(cert.AliasBank.numerator, cert.frozen.AliasBank.numerator)
        with ctx.workprec(256):
            good, dirty, paired = cert.interval.components({}, self.composed, (F(1, 512), F(1, 512)))
            self.assertEqual((good, dirty), ({}, {}))
            positive = tuple(sum((row[i] for row in paired.values()), arb(0)) for i in range(4))
            expected = cert.base.matrix_power(positive, 64)
            for actual, oracle in zip(self.bank.numerators[0], expected):
                self.assertTrue(actual.overlaps(oracle))
            self.assertTrue(all(x > 0 and x.is_finite() for x in self.bank.regional(1)[0]))

    def test_metadata_and_subclass_pins_survive_frozen_profile_wrapper(self):
        rows = []
        with patch.object(profiles, 'regional_matrices', return_value={0: tuple(map(arb, (1, 2, 3, 4)))}):
            receipt = profiles.weight_receipt(self.bank, 1, max_seconds=10, checkpoint=rows.append)
        self.assertTrue(receipt['complete_selected_weight'])
        for row in (receipt, rows[0], self.bank.metadata()):
            self.assertEqual(row['construction'], cert.local.construction(24))
            self.assertTrue(row['all_pair_before_alpha'])
            self.assertTrue(row['counts_outside_alpha'])
            self.assertTrue(row['construction_unchanged'])
            self.assertEqual(row['raw_histogram_pair_products'], 561**2)
            self.assertEqual(row['raw_occupancy_pair_products'], 65**2)
            self.assertEqual(row['composed_histogram_terms'], 2145)
            self.assertEqual(row['maximum_paired_occupancy'], 128)
            self.assertEqual(row['single_components'], 'identically zero')
            self.assertNotIn('raw_good_pair_products', row)
            self.assertFalse(row['whole_code_certificate'])
            self.assertIn('hd_field128_gamma_all_pair_certificate.py', row['bank_source_sha256'])
            self.assertEqual(row['bank_source_sha256'], self.bank.pins)
        for name in ('hd_field128_gamma_all_pair_certificate.py', 'hd_field128_gamma_profile_certificate.py',
                     'hd_field128_full_rs_gamma_alias_certificate.py', 'hd_field128_full_rs_gamma_certificate.py',
                     'hd_shifted_good_pair_certificate.py', 'hd_shifted_pairs_count.py'):
            self.assertIn(name, self.bank.pins)

    def test_exact_scope_resource_and_source_guards(self):
        for slots in (0, 33, 1.5):
            with self.assertRaises(ValueError):
                cert.histogram_geometry(slots)
        table = {j: (arb(1),)*4 for j in range(3)}
        for alpha in (0, 2, .5):
            with self.assertRaises(ValueError):
                cert.all_pair_coefficients(table, alpha, slots=1)
        with self.assertRaises(ValueError):
            cert.all_pair_coefficients({0: (arb(1),)*4}, F(1), slots=1)
        with self.assertRaises(TimeoutError):
            cert.all_pair_coefficients(table, F(1), slots=1, budget=cert.base.Budget(1e-12))
        for kwargs in (dict(state=32), dict(q=256), dict(grid=32), dict(bits=8),
                       dict(precision=64), dict(max_table_bytes=1), dict(outer_measure='unknown')):
            with self.assertRaises(ValueError):
                cert.AliasBank(F(1, 2), F(1, 4), ((F(1, 2), F(1, 3)),), **kwargs)
        with patch.object(cert, 'source_pins', return_value={}):
            with self.assertRaises(RuntimeError):
                self.bank.check_pins()


if __name__ == '__main__':
    unittest.main()
