"""Independent tiny polynomial arithmetic and source-partition tests."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product
import unittest

from flint import arb, ctx

from hd_xor_geometry_strata import category_resources
from hd_xor_joint_columns import categories, joint_count_arb
from hd_xor_sparse_mixed_certificate import checked_witnesses, sparse_count_arb, verify_receipts
from hd_xor_sparse_mixed_refine import MU_BINS, geometry_count, metadata, source_counts


class SparseMixedCertificateTests(unittest.TestCase):
    def proposal(self):
        result = metadata()
        result['screen'] = dict(complete_cohort=True, entries=[
            dict(source_matching_counts=list(source), mu_range=list(bounds),
                 mark='249/250', log2_potential='-64')
            for source in source_counts() for bounds in MU_BINS])
        return result

    def test_exact_toy_positive_tail(self):
        n1, n2, n, k, q = 1, 1, 3, 2, 4
        caps = {(s, d): F(1, 3**s*7**d) for s, d in categories(2)}
        expected = F(0)
        for columns in product(range(8), repeat=n):
            degree = sum(column.bit_count() for column in columns)
            if degree < 6:
                continue
            value = F(3**(degree-3))
            for column in columns:
                one, pair_size = column & 1, (column >> 1).bit_count()
                value *= caps[one+(pair_size == 1), int(pair_size == 2)]
            expected += value
        with ctx.workprec(256):
            actual = sparse_count_arb({key: arb(value.numerator)/value.denominator for key, value in caps.items()},
                n1, n2, n, k, q)
            self.assertTrue(actual.contains(arb(expected.numerator)/expected.denominator))

    def test_homogeneous_reduction(self):
        with ctx.workprec(256):
            caps = {(s, d): arb(1)/(3**s*7**d) for s, d in categories(3)}
            self.assertTrue(sparse_count_arb(caps, 0, 3, 4, 2, 4).overlaps(joint_count_arb(caps, 3, 4, 2, 4)))

    def test_complete_scope_checked_without_trusting_saved_values(self):
        proposal = self.proposal()
        self.assertEqual(len(checked_witnesses(proposal)), 324)
        altered = deepcopy(proposal)
        altered['screen']['entries'].pop()
        with self.assertRaises(ValueError):
            checked_witnesses(altered)
        altered = deepcopy(proposal)
        altered['screen']['entries'].append(altered['screen']['entries'][0])
        with self.assertRaises(ValueError):
            checked_witnesses(altered)
        altered = deepcopy(proposal)
        altered['cohort_counts'] = [0, 16, 1]
        with self.assertRaises(ValueError):
            checked_witnesses(altered)

    def test_singleton_does_not_invalidate_pair_multiplicity_resources(self):
        # One fixed singleton row and one fixed two-row group. Enumerate
        # every column support, and both cancellation branches when full.
        for singleton in range(8):
            for edge in combinations(range(8), 2):
                for one in (0, 1):
                    cases = [(0, 0, 0, [0]*8)]
                    for vertex in edge:
                        cases.append((1, 0, 0, [int(i != vertex) for i in range(8)]))
                    cases.extend([(0, 1, 0, [1]*8), (0, 1, 1, [int(i in edge) for i in range(8)])])
                    for single, double, cancelled, mask in cases:
                        occupied = [mask[i]+one*int(i != singleton) for i in range(8)]
                        lower, upper, pair_limit, target = category_resources(one+single, double, cancelled, 1)
                        self.assertTrue(all(lower <= value <= upper for value in occupied))
                        self.assertTrue(all(occupied[i]+occupied[j] <= pair_limit for i, j in combinations(range(8), 2)))
                        self.assertEqual(sum(value-lower for value in occupied), target)

    def test_dual_receipt_union_and_fail_closed_partition(self):
        proposal = self.proposal()
        receipts = []
        for precision in (256, 384):
            receipt = metadata()
            receipt.update(directed=True, precision_bits=precision, complete_cohort=True,
                source_classes=81, all_source_classes=81, fitted_new_witnesses=0,
                proposal_requires_fresh_directed_replay=False, failure_log2_upper='-992',
                entries=[])
            for source, bounds in checked_witnesses(proposal):
                receipt['entries'].append(dict(source_matching_counts=list(source), mu_range=list(bounds),
                    mark='249/250', log2_potential='-64', failure_log2_upper='-1000',
                    source_supports=source_counts()[source], cohort_colorings=17,
                    one_row_choices=8, two_row_geometries=geometry_count(*bounds)))
            receipts.append(receipt)
        result = verify_receipts(proposal, receipts)
        self.assertEqual(result['verified_cells'], 324)
        self.assertTrue(-992 < F(result['rounded_cell_union_log2_upper']) < -991)
        self.assertFalse(result['whole_code_certificate'])
        altered = deepcopy(receipts)
        altered[1]['entries'][0]['failure_log2_upper'] = '-999'
        with self.assertRaises(ValueError):
            verify_receipts(proposal, altered)
        altered = deepcopy(receipts)
        altered[1]['entries'][0]['cohort_colorings'] = 17**32
        with self.assertRaises(ValueError):
            verify_receipts(proposal, altered)
        altered = deepcopy(receipts)
        altered[1]['entries'].pop()
        with self.assertRaises(ValueError):
            verify_receipts(proposal, altered)


if __name__ == '__main__':
    unittest.main()
