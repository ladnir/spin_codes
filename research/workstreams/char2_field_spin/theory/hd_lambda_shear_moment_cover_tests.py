from collections import Counter
from fractions import Fraction as F
from itertools import product
from unittest.mock import patch
import unittest

import hd_lambda_shear_moment_cover as cover
from hd_lambda_column_shear_ranks import paired_rows,pack,rank_gf4,column_rows

MUL = ((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))


def toy_counts(columns,checks):
    counts = Counter()
    for values in product(range(4),repeat=len(columns)):
        syndrome = [0]*checks
        for column,value in zip(columns,values):
            for row in range(checks):
                syndrome[row] ^= MUL[(column>>(2*row))&3][value]
        if not any(syndrome):
            counts[sum((value!=0)<<i for i,value in enumerate(values))] += 1
    return counts


class ShearMomentCoverTests(unittest.TestCase):
    def test_actual_parity_annihilates_one_pair_generator(self):
        # Constructor-only check, not the eight-pair production witness gate.
        source,column = tuple(range(12)),6
        parity = cover.parity_columns(source,column)
        generator = sum(paired_rows(source,column,1),())
        self.assertEqual(rank_gf4(parity,8),8)
        for input_index in range(24):
            for check in range(8):
                value = 0
                for b,g in zip(parity,generator):
                    value ^= MUL[(b>>(2*check))&3][(g>>(2*input_index))&3]
                self.assertEqual(value,0)
        # A dual generator must receive T^transpose, not T. These differ.
        outside = tuple(range(12,16))
        left,right = column_rows(outside,column),column_rows(outside,column+1)
        wrong_first,wrong_second = [],[]
        for row in range(16):
            a,b = left[row],right[row]<<8
            wrong_first.append(a|b if not row&1 else a)
            wrong_second.append(b if not row&1 else a|b)
        wrong = wrong_first+wrong_second
        failures = 0
        for input_index in range(24):
            for check in range(8):
                value = 0
                for b,g in zip(wrong,generator):
                    value ^= MUL[(b>>(2*check))&3][(g>>(2*input_index))&3]
                failures += bool(value)
        self.assertGreater(failures,0)

    def test_independence_witness_and_nonmaximum_scope(self):
        columns = (pack((1,0)),)*4+(pack((0,1)),)*4
        witness = cover.independent_partition(columns,2,max_attempts=4,target_uncovered=0)
        self.assertEqual(witness['covered'],8)
        self.assertEqual(witness['sizes'],[2,2,2,2])
        self.assertFalse(witness['maximum_coverage_claimed'])
        self.assertEqual(cover.verify_partition(columns,2,witness['colors'])['covered'],8)
        zeros = cover.independent_partition((0,)*5,2,max_attempts=2,target_uncovered=0)
        self.assertEqual(zeros['covered'],0)
        self.assertFalse(zeros['target_reached'])

    def test_partition_validator_rejects_dependence_and_overlap(self):
        columns = (1,1,4,4)
        for colors in (((0,1),(),(),()), ((0,),(0,),(),()), ((4,),(),(),())):
            with self.assertRaises(ValueError):
                cover.verify_partition(columns,2,colors)

    def test_exact_toy_density_moment_is_below_cover(self):
        for columns,checks in (((1,1,4,4),2), ((1,4,5,6,9),2), ((0,1,1,1),1)):
            self.assertEqual(rank_gf4(columns,checks),checks)
            counts = toy_counts(columns,checks)
            width = len(columns)
            self.assertEqual(sum(counts.values()),4**(width-checks))
            actual = F(4)**(4*checks-width)*sum(
                (F(count**4,3**(3*mask.bit_count())) for mask,count in counts.items()),F(0))
            witness = cover.independent_partition(columns,checks,max_attempts=2,target_uncovered=0)
            upper = cover.moment_upper_from_sizes(width,checks,witness['sizes'],q=4)
            self.assertGreaterEqual(upper,actual)

    def test_size_sum_matches_all_tiny_supports(self):
        # Fix a partition with sizes2,1,1,0 and one uncolored coordinate.
        colors = ((0,1),(2,),(3,),())
        q,width,checks = 4,5,2
        expected = F(0)
        for mask in range(1<<width):
            rank_lower = max(sum((mask>>i)&1 for i in color) for color in colors)
            expected += F(q-1)**(mask.bit_count()-4*rank_lower)
        expected *= F(q)**(4*checks-width)
        self.assertEqual(cover.moment_upper_from_sizes(width,checks,(2,1,1,0),q=q),expected)
        self.assertEqual(cover.moment_upper_from_sizes(5,0,(0,0,0,0),q=4),1)

    def test_positive_activity_tail(self):
        values = tuple(F(i+1,9) for i in range(8))
        expected = F(0)
        for flags in product((0,1),repeat=8):
            if sum(flags) >= 5:
                term = F(1)
                for flag,value in zip(flags,values):
                    if flag:
                        term *= value
                expected += term
        self.assertEqual(cover.activity_tail_upper(values),expected)

    def test_mock_receipt_scope_and_timeout(self):
        def row(source,column,**kwargs):
            return dict(source=list(source),column_pair=[column,column+1],
                witness=dict(uncovered=[]),fourth_moment_upper='2',
                nonzero_fourth_moment_upper='1')
        with patch.object(cover,'pair_moment_upper',side_effect=row):
            data = cover.receipt()
        self.assertTrue(data['complete'])
        self.assertEqual(len(data['pairs']),8)
        self.assertEqual(data['fixed_A_restricted_word_moment_upper'],'93')
        self.assertFalse(data['all_labelled_rowsets_covered'])
        self.assertFalse(data['failure_probability_bound'])
        self.assertEqual(data['regional_provider_runs'],0)
        with patch.object(cover,'pair_moment_upper',side_effect=TimeoutError('tiny mock timeout')):
            data = cover.receipt()
        self.assertFalse(data['complete'])
        self.assertIsNone(data['fixed_A_restricted_word_moment_upper'])

    def test_guards_and_pins(self):
        self.assertEqual(cover.dependencies()['hd_lambda_column_shear_ranks.py'],
                         cover.PINS['hd_lambda_column_shear_ranks.py'])
        self.assertIn('hd_lambda_outer.py',cover.dependencies())
        with patch.object(cover,'lambda_matrix',return_value=((0,)*16,)*16):
            with self.assertRaises(ValueError):
                cover.dependencies()
        with self.assertRaises(ValueError):
            cover.parity_columns(tuple(range(11)),0)
        with self.assertRaises(ValueError):
            cover.parity_columns(tuple(range(12)),1)
        with self.assertRaises(ValueError):
            cover.independent_partition((1,),1,max_attempts=65)
        with self.assertRaises(ValueError):
            cover.moment_upper_from_sizes(32,8,(9,8,8,7))
        with self.assertRaises(ValueError):
            cover.activity_tail_upper((F(1),)*7)


if __name__ == '__main__':
    unittest.main()
