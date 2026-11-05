from math import comb, log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_two_sum_late_alpha_source_f1 as gate


class SourceF1Tests(TestCase):
    def test_exact_source_partition_scope(self):
        counts = [1]*10
        counts[1] = comb(8192, 18)-9
        data = dict(source_size=8192, field_order=256, forbidden_weight=gate.LIMIT,
            complete_requested_anchors=True, banks=[dict(active_groups=18,
                mark=str(gate.MARK), alpha=str(gate.ALPHA),
                source_strata=[dict(full_max=f, count=n) for f, n in enumerate(counts)])])
        self.assertEqual(gate.checked_source_count(data), counts[1])
        data['banks'][0]['source_strata'][0]['count'] += 1
        with self.assertRaises(ValueError):
            gate.checked_source_count(data)
        data['banks'][0]['source_strata'][0]['count'] -= 1
        data['banks'][0]['source_strata'][9]['full_max'] = 8
        with self.assertRaises(ValueError):
            gate.checked_source_count(data)

    def test_subset_full_pairs_cannot_exceed_source(self):
        # Both prescribed cyclic matchings and every subsetU of every G.
        length = 8
        matchings = (tuple((i, i+1) for i in range(0, length, 2)),
                     tuple((i, (i+1) % length) for i in range(1, length, 2)))
        for source in range(1 << length):
            full = [sum(source >> i & 1 and source >> j & 1 for i, j in matching)
                    for matching in matchings]
            if max(full) != 1:
                continue
            subset = source
            while True:
                for matching in matchings:
                    actual = sum(subset >> i & 1 and subset >> j & 1 for i, j in matching)
                    self.assertLessEqual(actual, min(1, subset.bit_count()//2))
                if subset == 0:
                    break
                subset = (subset-1) & source

    def test_f2_and_higher_never_enter_restricted_curve(self):
        regions = {(m, f): np.full((2, 2), f) for m in range(7) for f in range(m//2+1)}
        selected = gate.restrict_regions(regions, maximum=6)
        for (m, f), matrix in selected.items():
            np.testing.assert_array_equal(matrix, np.full((2, 2), min(1, m//2)))
        with self.assertRaises(ValueError):
            gate.restrict_regions({}, maximum=6)

    def test_one_provider_and_two_controls_accounting(self):
        regions = {(m, f): np.full((2, 2), -m) for m in range(19) for f in range(m//2+1)}
        comparison = dict(late_column_log2=-20., early_mask_log2=-10.)
        with patch.object(gate, 'load_source_count', return_value=123):
            with patch.object(gate, 'forced_nonzero_regions', return_value=regions) as provider:
                with patch.object(gate, 'compare_contributions', return_value=comparison):
                    with patch.object(gate, 'quiet_mixture_coefficients', return_value=np.zeros(145)):
                        with patch.object(gate, 'word_tail', return_value=-30.):
                            result = gate.run()
        provider.assert_called_once_with(float(gate.MARK), 18)
        expected = log2(123)+18*log2(comb(16, 8))+.2*(-gate.LIMIT*log2(float(gate.MARK))-log2(255))
        self.assertAlmostEqual(result['line_control']['log2_global_prefix'], expected)
        self.assertAlmostEqual(result['moment_control']['majorant_term_log2'], expected-30.)
        self.assertEqual(result['source_count'], '123')
        self.assertFalse(result['all_sources'])
        self.assertTrue(result['all_sources_in_declared_stratum'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['requested_controls'], 2)


if __name__ == '__main__':
    main()
