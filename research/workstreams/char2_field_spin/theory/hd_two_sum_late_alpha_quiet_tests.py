from itertools import product
from math import comb, log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_two_sum_late_alpha_quiet as gate
from hd_two_sum_closed_kernel import zero_two_moments


class LateAlphaQuietTests(TestCase):
    def test_positive_convolution_and_power(self):
        left, right = (2, 0, 3), (1, 4)
        expected = np.convolve(left, right)
        with np.errstate(divide='ignore'):
            np.testing.assert_allclose(2**gate.log_convolve(np.log2(left), np.log2(right)), expected)
            np.testing.assert_allclose(2**gate.log_power(np.log2(left), 3),
                                       np.convolve(np.convolve(left, left), left))

    def test_exact_labelled_quiet_counts(self):
        # GF4 labels0,1,2 have one nonzero assignment per nonzero scalar.
        # Count every labelled support of a4 independently of the helper.
        q, a = 4, 4
        multiply = lambda x, y: ((0, 0, 0, 0), (0, 1, 2, 3),
                                (0, 2, 3, 1), (0, 3, 1, 2))[x][y]
        counts = [0]*(a+1)
        for word in product(range(q), repeat=a):
            first = second = 0
            for u, x in enumerate(word):
                first ^= x
                second ^= multiply(u, x)
            if first == second == 0:
                counts[sum(x != 0 for x in word)] += 1
        self.assertEqual(counts, [comb(a, w)*zero_two_moments(w, q) for w in range(a+1)])
        np.testing.assert_allclose(2**gate.quiet_mask_logs(a, 0, q=q), counts)

    def test_late_and_early_alpha_one_agree(self):
        late, early = gate.column_coefficients(2, 4, -0.4, 1, q=4, columns=2, dimension=1)
        np.testing.assert_array_equal(late, early)

    def test_late_sums_labels_before_power(self):
        a, r, alpha = 4, 2, 0.5
        masks = [(mask.bit_count(), zero_two_moments(mask.bit_count(), 4))
                 for mask in range(1 << a)]
        raw, old = np.zeros(a*r+1), np.zeros(a*r+1)
        for choices in product(masks, repeat=r):
            weight = sum(row[0] for row in choices)
            count = choices[0][1]*choices[1][1]
            raw[weight] += count
            old[weight] += count**alpha
        late, early = gate.column_coefficients(r, a, 0, alpha, q=4, columns=2, dimension=1)
        factor = 3**(alpha*a*r/2)
        np.testing.assert_allclose(2**late*factor, raw**alpha)
        np.testing.assert_allclose(2**early*factor, old)
        self.assertTrue(np.all(raw**alpha <= old+1e-12))
        self.assertLess(raw[6]**alpha, old[6])

    def test_exact_aggregate_tail_and_split_normalizer(self):
        # Alpha1 avoids irrational arithmetic: independently square each
        # column's integer counts, then divide the WHOLE-WORD normalizer.
        q, r, a, columns, dimension, minimum = 4, 2, 3, 2, 1, 1
        base = [comb(a, w)*zero_two_moments(w, q) for w in range(a+1)]
        column = np.convolve(base, base)
        word = np.convolve(column, column)
        exact = sum(word[minimum*a*r:])/3**(dimension*a*r)
        actual = gate.compare_contributions(r, a, 0, 1, q=q, columns=columns,
                                            dimension=dimension, minimum=minimum)
        self.assertAlmostEqual(actual['late_column_log2'], log2(exact), places=12)
        self.assertEqual(actual['late_column_log2'], actual['early_mask_log2'])

    def test_complete_envelope_and_scope_mock(self):
        regions = {(m, f): np.array([[-0.25*m, -np.inf], [-np.inf, -1.]])
                   for m in range(gate.GROUPS+1) for f in range(m//2+1)}
        envelope = gate.zero_intercept_envelope(regions, gate.GROUPS)
        self.assertEqual(envelope['log2_slope'], -0.25)
        # No regional provider or production-sized polynomial is executed.
        mock_comparison = dict(late_column_log2=-20., early_mask_log2=-10.)
        with patch.object(gate, 'forced_nonzero_regions', return_value=regions) as provider:
            with patch.object(gate, 'compare_contributions', return_value=mock_comparison):
                result = gate.run()
        provider.assert_called_once_with(float(gate.MARK), gate.GROUPS)
        self.assertTrue(result['complete_requested_control'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['full_quiet_event_bounded'])
        self.assertFalse(result['failure_probability_lower_bound'])
        self.assertEqual(result['raw_regional_banks'], 1)
        self.assertEqual(result['comparison']['late_minus_early_bits'], -10.)

    def test_guards_and_budget(self):
        with self.assertRaises(ValueError):
            gate.log_power([0.], -1)
        with self.assertRaises(ValueError):
            gate.zero_intercept_envelope({}, 3)
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)
        def expired():
            raise TimeoutError('unit guard')
        with self.assertRaises(TimeoutError):
            gate.log_convolve([0.], [0.], expired)


if __name__ == '__main__':
    main()
