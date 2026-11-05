"""Tiny independent fixed-flow and saved-row tests; no coefficient banks."""

from copy import deepcopy
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from flint import ctx

import hd_field128_gamma_midpoint_receipt_audit as audit
from hd_field128_gamma_near_face_certificate_tests import fixture
from hd_shifted_pairs_count import exact_profiles


def tiny_record():
    pairs, weight = 4, 3
    regions = {key: tuple(audit.base.rational(value) for value in row) for key, row in fixture(pairs, weight).items()}
    reference = audit.producer.midpoint.MidpointCounts(pairs, weight, max_singletons=2)
    with ctx.workprec(256):
        support, rows = audit.producer.transfer_sum(pairs, weight, regions, counts=reference, repetitions=8)
        record = dict(rows, active_groups=weight, field_order=256, forbidden_weight=17,
                      outer_measure='invertible-toeplitz', witness=dict(z='1/2', alpha='1/3'),
                      expected_even_rows=2, expected_profile_pairs=4,
                      count_provenance=dict(reference.diagnostics(), float_log_count_path_used=False),
                      support_sum=audit.producer.frozen.frozen.support_log(support))
        value = audit.previous.outer_prefix(record, weight)+audit.tail.log2arb(support)
        record.update(log2_upper_rational=str(audit.base.rational_ceiling(value)), log2_upper_enclosure=str(value))
    return record


class MidpointAuditTests(unittest.TestCase):
    def test_direct_central_formula_and_final_minimum_every_tiny_profile(self):
        for pairs in range(1, 13):
            for weight in range(2*pairs+1):
                direct = audit.DirectCounts(pairs, weight, 2)
                reference = audit.producer.midpoint.MidpointCounts(pairs, weight, max_singletons=2)
                for first in direct.full:
                    for second in direct.full:
                        self.assertEqual(direct.midpoint(first, second), reference.midpoint_upper(first, second))
                        self.assertEqual(direct.count(first, second), reference.count(first, second))
                self.assertLessEqual(direct.choose.cache_info().currsize, 16384)

    def test_exhaustive_indexed_supports_and_fixed_flow_center(self):
        exact = audit.producer.midpoint.exact
        for pairs in range(1, 7):
            actual = exact_profiles(pairs)
            for weight in range(2*pairs+1):
                direct = audit.DirectCounts(pairs, weight, 2)
                for first in direct.full:
                    for second in direct.full:
                        self.assertGreaterEqual(direct.count(first, second), actual.get((weight, first, second), 0))
                        p, r, a, b = exact._canonical(pairs, weight, first, second)
                        if r == p and a == b == 0:
                            expected = 2
                        else:
                            singles = r-2*a
                            expected = (singles+1)*exact.fixed_flow_count(p, r, a, b, singles//2)
                        self.assertEqual(direct.midpoint(first, second), expected)

    def test_saved_row_reaggregation_and_fresh_prefix(self):
        record = tiny_record()
        result = audit.audit_record(record, pairs=4)
        self.assertEqual(result['profile_pairs'], 4)
        self.assertEqual(result['even_rows'], 2)
        self.assertTrue(result['all_producer_midpoint_and_final_counts_match'])
        self.assertLessEqual(F(result['saved_row_reaggregation_log2_upper_rational']), F(record['log2_upper_rational']))

    def test_corrupt_counts_rows_prefix_and_ceiling_rejected(self):
        record = tiny_record()
        with patch.object(audit.producer.midpoint.MidpointCounts, 'midpoint_upper', return_value=0):
            with self.assertRaisesRegex(ValueError, 'independent fixed flow'):
                audit.audit_record(record, pairs=4)
        for change in (lambda row: row['rows'].pop(),
                       lambda row: row['rows'][0].update(profiles_with_positive_count=0),
                       lambda row: row.update(forbidden_weight=18),
                       lambda row: row.update(log2_upper_rational='-1000000')):
            bad = deepcopy(record)
            change(bad)
            with self.assertRaises(ValueError):
                audit.audit_record(bad, pairs=4)
        self.assertIn('../HD_FIELD128_MIDPOINT_COUNTS.md', audit.source_pins())


if __name__ == '__main__':
    unittest.main()
