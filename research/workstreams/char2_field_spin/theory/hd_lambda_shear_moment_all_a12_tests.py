from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from unittest.mock import patch
import unittest

from flint import arb, ctx

import hd_lambda_shear_moment_all_a12 as gate


def toy_record():
    # Four labelled bases in a tiny explicitly packed8-dimensional matroid.
    columns = tuple(1<<(2*i) for _ in range(4) for i in range(8))
    colors = [list(range(8*k,8*k+8)) for k in range(4)]
    upper = gate.cover.moment_upper_from_sizes(32,8,(8,)*4)
    tail = gate.cover.activity_tail_upper((upper-1,)*8)
    row = dict(source=list(range(12)),restricted_word_moment_upper=str(tail),pairs=[
        dict(column_pair=[j,j+1],parity_columns=list(columns),colors=deepcopy(colors),
             sizes=[8]*4,uncovered=[],fourth_moment_upper=str(upper),
             nonzero_fourth_moment_upper=str(upper-1)) for j in range(0,16,2)])
    return columns,row,tail


class ShearMomentAllA12Tests(unittest.TestCase):
    def test_exact_witness_validation_and_saved_values_not_trusted(self):
        columns,row,tail = toy_record()
        with patch.object(gate.cover,'parity_columns',return_value=columns):
            self.assertEqual(gate.validate_record(row,tuple(range(12))),tail)
            for kind in ('parity','colors','sizes','moment','tail','order','source'):
                bad = deepcopy(row)
                if kind == 'parity': bad['pairs'][0]['parity_columns'][0] ^= 1
                if kind == 'colors': bad['pairs'][0]['colors'][1][0] = 0
                if kind == 'sizes': bad['pairs'][0]['sizes'][0] = 7
                if kind == 'moment': bad['pairs'][0]['fourth_moment_upper'] = '1'
                if kind == 'tail': bad['restricted_word_moment_upper'] = '1'
                if kind == 'order': bad['pairs'][0]['column_pair'] = [2,3]
                if kind == 'source': bad['source'][0],bad['source'][1] = 1,0
                with self.assertRaises(ValueError,msg=kind):
                    gate.validate_record(bad,tuple(range(12)))

    def test_mock_source_record_keeps_order_and_positive_tail(self):
        _,row,tail = toy_record()
        calls = []
        def fake(source,column,**kwargs):
            calls.append((tuple(source),column))
            pair = row['pairs'][column//2]
            return dict(pair,witness=dict(colors=pair['colors'],sizes=pair['sizes'],
                uncovered=[],attempts=1))
        with patch.object(gate.cover,'pair_moment_upper',side_effect=fake):
            actual = gate.source_record(tuple(range(12)))
        self.assertEqual(calls,[(tuple(range(12)),j) for j in range(0,16,2)])
        self.assertEqual(F(actual['restricted_word_moment_upper']),tail)

    def test_mock_inventory_covers_all_labels_and_preserves_partial(self):
        progress = []
        with patch.object(gate,'source_record',side_effect=lambda source,**k:dict(source=list(source))):
            data = gate.inventory(progress=progress.append)
        self.assertTrue(data['complete'])
        self.assertTrue(data['all_labelled_rowsets_covered'])
        self.assertEqual([tuple(row['source']) for row in data['rows']],list(gate.SOURCES))
        self.assertEqual(len(progress),15)
        with patch.object(gate,'source_record',side_effect=[{},TimeoutError('tiny mocked budget')]):
            partial = gate.inventory()
        self.assertFalse(partial['complete'])
        self.assertEqual(partial['completed_rowsets'],1)
        with self.assertRaises(ValueError):
            gate.validate_inventory(partial)

    def test_inventory_rejects_missing_duplicate_or_wrong_label_order(self):
        data = dict(gate.SCOPE,complete=True,all_labelled_rowsets_covered=True,
            completed_rowsets=1820,reviewed_source_sha256=gate.dependencies(),
            rows=[dict(source=list(source)) for source in gate.SOURCES])
        def check(row,source,**kwargs):
            if row['source'] != list(source): raise ValueError('label mismatch')
            return F(16)
        with patch.object(gate,'validate_record',side_effect=check):
            self.assertEqual(gate.validate_inventory(data),(F(16),)*1820)
            bad = deepcopy(data)
            bad['rows'][1] = bad['rows'][0]
            with self.assertRaises(ValueError): gate.validate_inventory(bad)
            bad = deepcopy(data)
            bad['rows'].pop()
            with self.assertRaises(ValueError): gate.validate_inventory(bad)

    def test_raw_group_has_one_beta_and_exactly_one_label_sum(self):
        with ctx.workprec(256):
            raw,details = gate.raw_group_upper((F(16),F(81),F(16)))
            beta = F(256)**192/F(255)**96
            self.assertTrue(raw.contains(gate.ideal.A.uniform.rational(7*beta)))
            self.assertEqual(details['rowset_terms'],3)
            self.assertEqual(details['distinct_tail_values'],2)
            self.assertEqual(F(details['beta']),beta)
            repeated,_ = gate.raw_group_upper((F(16),)*1820)
            self.assertTrue(repeated.contains(gate.ideal.A.uniform.rational(3640*beta)))

    def test_checked_inventory_requires_hash_before_validation(self):
        raw = b'{"mock":true}'
        with patch.object(Path,'read_bytes',return_value=raw), \
             patch.object(gate,'validate_inventory',return_value=(F(16),)*1820) as validate:
            with self.assertRaises(ValueError): gate.checked_inventory('unused','0'*64)
            validate.assert_not_called()
            tails,pin = gate.checked_inventory('unused',sha256(raw).hexdigest())
        self.assertEqual(len(tails),1820)
        self.assertEqual(pin['verified_pair_witnesses'],14560)

    def test_mock_certificate_upper_scope_and_no_provider(self):
        # No inventory or reference calculation. Files are checked as exact
        # byte pins through the same post-evaluation checks as production.
        raw = b'mock immutable input'
        digest = sha256(raw).hexdigest()
        reference = dict(file='mock_reference',sha256=digest,precision_bits=256)
        inventory = dict(file='mock_inventory',sha256=digest,verified_rowsets=1820)
        with patch.object(gate,'dependencies',return_value={'mock':'pin'}), \
             patch.object(gate,'checked_inventory',return_value=((F(16),)*1820,inventory)), \
             patch.object(gate.ideal,'checked_reference',return_value=(arb(-7000),reference)), \
             patch.object(gate.ideal,'method_expression',return_value=arb(-45)), \
             patch.object(Path,'read_bytes',return_value=raw), \
             patch.object(gate.cover,'independent_partition',side_effect=AssertionError('No search')), \
             patch.object(gate.ideal.A.uniform,'regional_directed',side_effect=AssertionError('No provider')):
            data = gate.certificate('mock_inventory',digest)
        self.assertTrue(data['selected_class_40_bits'])
        self.assertTrue(data['failure_probability_upper_bound'])
        self.assertFalse(data['whole_code_certificate'])
        self.assertFalse(data['mixed_arities_covered'])
        self.assertEqual(data['source_profiles'],99)
        self.assertEqual(data['rowset_terms'],1820)
        self.assertFalse(data['additional_rowset_binomial_factor'])
        self.assertEqual(data['projective_normalizations'],1)

    def test_guard_and_current_provenance(self):
        self.assertEqual(gate.dependencies()['hd_lambda_shear_moment_cover.py'],gate.PINS['hd_lambda_shear_moment_cover.py'])
        for value in (0,-1,float('nan'),float('inf'),601):
            with self.assertRaises(ValueError): gate._timer(value)
        with self.assertRaises(ValueError): gate.inventory(max_attempts=65)
        with self.assertRaises(ValueError): gate.raw_group_upper(())
        with self.assertRaises(ValueError): gate.raw_group_upper((F(-1),))
        with self.assertRaises(ValueError): gate.certificate('unused','0'*64,precision=128)


if __name__ == '__main__':
    unittest.main()
