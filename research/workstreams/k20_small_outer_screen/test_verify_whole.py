"""Exact arithmetic, forged-receipt rejection, and optional local receipt checks."""
from contextlib import ExitStack
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import verify_whole as verify


def encoded(value):
    value = F(value)
    return dict(numerator_hex=hex(value.numerator), exponent=-(value.denominator.bit_length()-1))


def whole_fixture():
    values = {str(q):[1, -100] for q in range(1, 8193)}
    chosen = {str(q):0 for q in range(1, 8193)}
    later = {'1':[1, -101], '2':[1, -100]}
    selected = dict(values, **{'1':later['1']})
    chosen['1'] = 1
    components = [dict(schema='fixture', map_record={'map_sha256':verify.MAP},
                       occupancy_uppers=value, source_sha256={'fixture':'f'*64})
                  for value in (values, later)]
    paths = [str((verify.HERE / f'fixture-{j}.json').resolve()) for j in range(2)]
    total = sum(map(verify.dyadic, selected.values()), F(0))
    tail = sum((verify.dyadic(pair) for q, pair in selected.items() if int(q) >= 3), F(0))
    record = dict(verify.COMMON, schema='k20-small-outer-whole-1', minimum_distance=209716,
        slots_per_region=8192, macro_slots=32, continuous_state=True, map_sha256=verify.MAP,
        map_record={'map_sha256':verify.MAP},
        all_occupancies_covered=True, target_bits=50, target_met=True, whole_code_certificate=True,
        source_sha256={'fixture':'f'*64}, components=[dict(path=path, sha256='f'*64, schema='fixture')
            for path in paths], occupancy_uppers=selected, selected_components=chosen,
        exact_union=encoded(total), exact_tail_union=encoded(tail),
        margin_bits_display=verify.margin(total), tail_margin_bits_display=verify.margin(tail))
    return record, dict(zip(paths, components))


class ArithmeticTests(unittest.TestCase):
    def test_dyadics_and_tiny_terms_are_exact(self):
        self.assertEqual(verify.dyadic([3, -4]), F(3, 16))
        self.assertEqual(verify.dyadic([3, 4]), 48)
        value = F(1, 1 << 50) + F(3, 1 << 10000)
        self.assertEqual(verify.exact_union(encoded(value)), value)
        self.assertEqual(verify.sum_dyadics([[1,-50], [3,-10000]]), value)
        self.assertGreater(value, F(1, 1 << 50))
        for pair in ([True, -1], [-1, 0], [1, False], [1]):
            with self.subTest(pair=pair), self.assertRaises(ValueError):
                verify.dyadic(pair)

    def test_noncanonical_union_is_rejected(self):
        with self.assertRaises(ValueError):
            verify.exact_union(dict(numerator_hex='0x01', exponent=-100))

    def test_independent_exact_rs_shells(self):
        counts = verify.exact_shells()
        self.assertEqual(len(counts), 65)
        self.assertEqual(sum(counts), (1 << 128)-1)
        self.assertEqual(counts[:9], (F(0),)*9)
        self.assertEqual(hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest(),
                         '3baf76cb89db74ebf82e78e9981464c31cf9be1ca6753a895b984840b811fc80')

    def test_component_coverage_requires_positive_complete_unique_endpoints(self):
        valid = dict(occupancy_covered=[3,4], occupancy_uppers={'3':[1,-100], '4':[1,-99]})
        self.assertEqual(verify.coverage(valid), valid['occupancy_uppers'])
        for changes in (dict(occupancy_covered=[3,3,4]), dict(occupancy_covered=[4,3]),
                        dict(occupancy_uppers={'3':[1,-100]}),
                        dict(occupancy_uppers={'3':[1,-100], '4':[0,0]})):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                verify.coverage(dict(valid, **changes))


class WholeSelectionTests(unittest.TestCase):
    """Mock component authentication only; exercise the complete 8192-term union."""
    def setUp(self):
        self.record, self.components = whole_fixture()
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        for name in ('pins', 'required', 'inherited', 'stable'):
            self.stack.enter_context(patch.object(verify.SourceAudit, name, return_value=None))
        self.stack.enter_context(patch.object(verify.SourceAudit, 'digest', return_value='f'*64))
        self.stack.enter_context(patch.object(verify, 'load_json', side_effect=self.components.__getitem__))
        self.stack.enter_context(patch.object(verify, 'audit_component',
                                              side_effect=lambda value, sources:value['occupancy_uppers']))

    def test_all_8192_terms_and_first_minimum_tie_rule(self):
        result = verify.audit_record(self.record)
        self.assertEqual(result['authenticated_occupancies'], 8192)
        self.assertTrue(result['target_met'])
        self.assertFalse(result['matrix_endpoints_independently_recomputed'])

    def test_reject_changed_endpoint(self):
        self.record['occupancy_uppers']['89'] = [1, -101]
        with self.assertRaisesRegex(ValueError, 'minimum endpoint selection'):
            verify.audit_record(self.record)

    def test_reject_tie_chosen_from_later_component(self):
        self.record['selected_components']['2'] = 1
        with self.assertRaisesRegex(ValueError, 'first-minimum'):
            verify.audit_record(self.record)

    def test_reject_missing_occupancy(self):
        del self.record['occupancy_uppers']['8192']
        with self.assertRaisesRegex(ValueError, 'coverage'):
            verify.audit_record(self.record)

    def test_reject_changed_union_and_tail(self):
        for key in ('exact_union', 'exact_tail_union'):
            with self.subTest(key=key):
                changed = deepcopy(self.record)
                changed[key]['numerator_hex'] = hex(int(changed[key]['numerator_hex'], 16)+1)
                with self.assertRaisesRegex(ValueError, 'exact union'):
                    verify.audit_record(changed)

    def test_reject_changed_map_geometry_and_probability_scope(self):
        for key, value in (('map_sha256', '0'*64), ('threshold', 209716),
                           ('continuous_state', False), ('zero_initial_state', 1)):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'geometry or scope'):
                verify.audit_record(dict(self.record, **{key:value}))

    def test_target_claim_uses_exact_union(self):
        with self.assertRaisesRegex(ValueError, 'exact comparison'):
            verify.audit_record(dict(self.record, target_met=False))
        changed = deepcopy(self.record)
        changed['target_bits'] = 60
        with self.assertRaisesRegex(ValueError, '40 or 50'):
            verify.audit_record(changed)

    def test_reject_component_hash_change(self):
        self.record['components'][0]['sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'component file hash'):
            verify.audit_record(self.record)


Q1 = verify.ROOT / 'tmp/k20-small-outer-q1-p192-v1.json'
LOCAL = verify.ROOT / 'tmp/k20-small-outer-physical-p192-v1-t002.json'


@unittest.skipUnless(Q1.exists() and LOCAL.exists(), 'ignored local inputs not available')
class LocalReceiptTests(unittest.TestCase):
    def test_q1_provenance_and_exact_shell_fold(self):
        record = verify.load_json(Q1)
        sources = verify.SourceAudit()
        values = verify.audit_component(record, sources)
        self.assertEqual(set(values), {'1'})
        self.assertGreater(verify.margin(verify.dyadic(values['1'])), 58.75)
        sources.stable()

    def test_reject_changed_sparse_fold(self):
        record = verify.load_json(Q1)
        record['q1_upper'][1] -= 1
        with self.assertRaisesRegex(ValueError, 'exact shell fold'):
            verify.sparse_fold(record, 1)

    def test_reject_changed_sparse_independence(self):
        record = verify.load_json(Q1)
        record['count_premises']['independent_setups_between_groups'] = False
        with self.assertRaisesRegex(ValueError, 'counting premises'):
            verify.audit_component(record, verify.SourceAudit())

    def test_map_identity_uses_explicit_arrays(self):
        record = verify.load_json(Q1)
        record['map_record']['expansion_rows_hex'][0] = '0x0'
        with self.assertRaisesRegex(ValueError, 'explicit map arrays'):
            verify.audit_map(record['map_record'], verify.SourceAudit(), record['source_sha256'])

    def test_local_physical_family_and_negative_entry(self):
        record = verify.load_json(LOCAL)
        pins = dict(record['source_sha256'])
        pins[str(LOCAL.resolve())] = hashlib.sha256(LOCAL.read_bytes()).hexdigest()
        verify.audit_local(LOCAL, verify.SourceAudit(), pins, record['map_record'])
        record['matrices'][1][0][0] = [-1, -100]
        with patch.object(verify, 'load_json', return_value=record), self.assertRaisesRegex(ValueError, 'nonnegative'):
            verify.audit_local(LOCAL, verify.SourceAudit(), pins, record['map_record'])


Q2 = verify.ROOT / 'tmp/k20-small-outer-low-complete-p256-v1-q2.json'
RANGE = verify.ROOT / 'tmp/k20-small-outer-range-t004-a08-v1.json'
DENSE = verify.ROOT / 'tmp/k20-small-outer-dense-transition-v2.json'


class OptionalFreshComponentsTests(unittest.TestCase):
    @unittest.skipUnless(Q2.exists(), 'ignored q2 receipt not available')
    def test_independent_q2_shell_fold(self):
        values = verify.audit_component(verify.load_json(Q2), verify.SourceAudit())
        self.assertGreater(verify.margin(verify.dyadic(values['2'])), 122)

    @unittest.skipUnless(RANGE.exists(), 'ignored range receipt not available')
    def test_reject_invalid_fractional_alpha(self):
        record = verify.load_json(RANGE)
        record['alpha'] = '0'
        with self.assertRaisesRegex(ValueError, 'fractional witness'):
            verify.audit_component(record, verify.SourceAudit())

    @unittest.skipUnless(DENSE.exists(), 'ignored dense receipt not available')
    def test_dense_coefficient_witnesses_are_authenticated(self):
        record = verify.load_json(DENSE)
        verify.audit_component(record, verify.SourceAudit())
        record['occupancy_choices'][str(record['occupancy_covered'][0])]['marker_probability'] = '1'
        with self.assertRaisesRegex(ValueError, 'dense witness'):
            verify.audit_component(record, verify.SourceAudit())


@unittest.skipUnless(verify.DEFAULT_RECEIPT.exists(), 'ignored whole receipt not available')
class CompletedReceiptTests(unittest.TestCase):
    def test_completed_receipt(self):
        self.assertTrue(verify.verify(verify.DEFAULT_RECEIPT)['target_met'])

    def test_changed_real_receipt_is_rejected(self):
        record = verify.load_json(verify.DEFAULT_RECEIPT)
        record['selected_components']['8192'] = -1
        with self.assertRaisesRegex(ValueError, 'selection'):
            verify.audit_record(record)


if __name__ == '__main__':
    unittest.main()
