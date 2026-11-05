from collections import Counter
from itertools import product
from math import comb, log2
from unittest.mock import patch
import json
import unittest

import numpy as np

import hd_lambda_column_shear_control as gate
from hd_lambda_column_shear_ranks import paired_rows
from hd_pair_identity_kernel import mask_table_from_entries, branch_log_matrices
from hd_pair_state_tail import word_tail
from hd_quaternary_outer_relations import category_data

MUL = ((0,0,0,0), (0,1,2,3), (0,2,3,1), (0,3,1,2))


def output_mask(rows, values):
    mask = 0
    for index, row in enumerate(rows):
        value = 0
        for i, x in enumerate(values):
            value ^= MUL[(row >> (2*i)) & 3][x]
        mask |= (value != 0) << index
    return mask


def fake_data(arity, sources, **kwargs):
    categories = category_data(2*arity)
    counts = categories['category_counts']
    n, branches = len(sources), len(counts)
    return dict(sources=tuple(sources), branch_count=branches,
        first=np.zeros((8,n,branches), dtype=np.uint16),
        second=np.zeros((8,n,branches), dtype=np.uint16),
        log_counts=np.broadcast_to(np.array([log2(c) for c in counts]), (n,branches)),
        input_weights=np.broadcast_to(np.array(categories['category_input_weights']), (n,branches)))


class ColumnShearControlTests(unittest.TestCase):
    def test_exact_joint_support_and_weight_against_GF4(self):
        for source, column in (((0,),0), ((0,1),0), ((3,12),6)):
            data = gate.pair_data(source,column,q=4)
            expected = Counter()
            rows = sum(paired_rows(source,column,1), ())
            for values in product(range(4), repeat=2*len(source)):
                expected[sum(x != 0 for x in values), output_mask(rows,values)] += 1
            actual = Counter()
            for w, left, right, n in zip(data['input_weights'],data['first'],data['second'],data['raw_counts']):
                if n:
                    actual[w, left | (right << 16)] += n
            self.assertEqual(actual, expected)
            self.assertEqual(actual[0,0],1)

    def test_raw_GF256_mass_and_input_weight_partition(self):
        for arity, size in ((1,7),(2,529)):
            data = gate.pair_data(tuple(range(arity)),14)
            self.assertEqual(len(data['raw_counts']),size)
            self.assertEqual(sum(data['raw_counts']),256**(2*arity))
            for w in range(2*arity+1):
                mass = sum(n for weight,n in zip(data['input_weights'],data['raw_counts']) if weight == w)
                self.assertEqual(mass,comb(2*arity,w)*255**w)
            # Raw counts, rather than independent-column probabilities.
            self.assertGreater(max(data['raw_counts']),1)

    def test_batched_labels_blocks_and_joint_shapes(self):
        sources = ((0,1),(3,12))
        data = gate.branch_data(2,sources)
        self.assertEqual(data['first'].shape,(8,2,529))
        self.assertEqual(data['second'].shape,(8,2,529))
        self.assertEqual(data['log_counts'].shape,(2,529))
        self.assertEqual(data['sources'],sources)
        for block in (0,3,7):
            for index,source in enumerate(sources):
                single = gate.pair_data(source,2*block)
                np.testing.assert_array_equal(data['first'][block,index],single['first'])
                np.testing.assert_array_equal(data['second'][block,index],single['second'])
        self.assertTrue(np.any(data['first'][0] != data['first'][3]))

    def test_joint_branch_chronology_and_r1_tail(self):
        # Noncommuting matrices make a column reversal observably wrong.
        zero = np.array(((.7,.2),(.1,.9)))
        nonzero = np.array(((.25,.35),(.4,.2)))
        self.assertGreater(np.max(abs(zero@nonzero-nonzero@zero)),.01)
        table = mask_table_from_entries(np.log2(zero),np.log2(nonzero),rows=16)
        data = gate.pair_data((5,),2,q=4)
        matrices = branch_log_matrices(np.array(data['first']),np.array(data['second']),table)
        counts = np.array([log2(n) if n else -np.inf for n in data['raw_counts']])
        actual = word_tail([matrices],counts,np.array(data['input_weights']),1,1)
        rows = sum(paired_rows((5,),2,1), ())
        expected = 0.
        for values in product(range(4),repeat=2):
            if not any(values):
                continue
            mask = output_mask(rows,values)
            current = np.eye(2)
            for row in range(32):
                current = current@(nonzero if mask >> row & 1 else zero)
            expected += current[0].sum()
        self.assertAlmostEqual(float(actual),log2(expected),places=11)

    def test_normalization_and_common_prefix_once(self):
        for a in (1,2):
            tails = np.array((18*(8*a*log2(255)+2),18*(8*a*log2(255)+3)))
            roots = gate.normalized_roots(tails,a)
            np.testing.assert_allclose(roots,(2,3),rtol=0,atol=1e-13)
            self.assertAlmostEqual(gate.root_sum_bound(roots,17),17+18*log2(12),places=10)

    def test_mock_full_scope_and_factor_wiring(self):
        saved = dict(groups=18,source_class='Fmax<=1',source_count='7',mark='9927/10000')
        seen = Counter()
        def tail(blocks,counts,weights,groups,minimum,**kwargs):
            a = int(np.max(weights))//2
            self.assertEqual(groups,18)
            self.assertEqual(minimum,9*a*18)
            self.assertEqual(len(blocks),8)
            self.assertEqual(blocks[0].shape,counts.shape+(2,2))
            self.assertAlmostEqual(float(np.logaddexp2.reduce(counts[0])),16*a)
            seen[a] += len(counts)
            return np.full(len(counts),18*(8*a*log2(255)+2))
        with patch.object(gate,'dependencies',return_value={'mock':'pinned'}), \
             patch.object(gate,'load_measure',return_value=(saved,{},17.)), \
             patch.object(gate,'mask_log_matrices',return_value=None), \
             patch.object(gate,'branch_data',side_effect=fake_data), \
             patch.object(gate,'branch_log_matrices',side_effect=lambda a,b,t:np.zeros((len(a),2,2))), \
             patch.object(gate,'word_tail',side_effect=tail):
            result = gate.run(batch_size=15)
        self.assertEqual(seen,{1:16,2:120})
        self.assertTrue(result['complete'])
        self.assertEqual(result['all_labelled_arities'],[1,2])
        for row in result['controls']:
            self.assertAlmostEqual(row['bound_log2'],17+18*(2+log2(comb(16,row['arity']))),places=9)
        self.assertAlmostEqual(result['requested_arity_union_log2'],17+18*(2+log2(136)),places=9)
        self.assertFalse(result['fresh_ratio_law_used'])
        self.assertFalse(result['independent_column_support_assumption'])
        self.assertEqual(result['permitted_gauges_per_column_pair'],1)
        json.dumps(result,allow_nan=False)

    def test_timeout_retains_only_completed_rows_not_complete_bound(self):
        saved = dict(groups=18,source_class='Fmax<=1',source_count='7',mark='9927/10000')
        calls = 0
        def tail(blocks,counts,weights,*args,**kwargs):
            nonlocal calls
            calls += 1
            if calls > 1:
                raise TimeoutError('mock budget')
            return np.zeros(len(counts))
        with patch.object(gate,'dependencies',return_value={'mock':'pinned'}), \
             patch.object(gate,'load_measure',return_value=(saved,{},17.)), \
             patch.object(gate,'mask_log_matrices',return_value=None), \
             patch.object(gate,'branch_data',side_effect=fake_data), \
             patch.object(gate,'branch_log_matrices',side_effect=lambda a,b,t:np.zeros((len(a),2,2))), \
             patch.object(gate,'word_tail',side_effect=tail):
            result = gate.run(arities=(1,),batch_size=8)
        self.assertFalse(result['complete'])
        self.assertEqual(len(result['fixed_A_rows']),8)
        self.assertEqual(result['all_labelled_arities'],[])
        self.assertIsNone(result['requested_arity_union_log2'])
        self.assertTrue(result['partial_rows_are_not_a_complete_class_bound'])
        json.dumps(result,allow_nan=False)

    def test_guards_and_frozen_pins(self):
        self.assertEqual(gate.dependencies()[gate.INPUT],gate.INPUT_SHA)
        for args in (((0,1,2),0), ((0,),1), ((0,0),0)):
            with self.assertRaises(ValueError):
                gate.pair_data(*args)
        with self.assertRaises(ValueError):
            gate.branch_data(2,((2,1),))
        with self.assertRaises(ValueError):
            gate.run(batch_size=16)
        with self.assertRaises(ValueError):
            gate.run(arities=(3,))


if __name__ == '__main__':
    unittest.main()
