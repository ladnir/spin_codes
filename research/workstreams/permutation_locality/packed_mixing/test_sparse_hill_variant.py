from contextlib import redirect_stdout
from fractions import Fraction as Q
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from flint import arb, ctx
import sparse_hill_variant as variant


def synthetic_search_record(sparse, occupancy=1, updates=4, refinement='ekr'):
    counts = tuple([Q(0)]*6+[Q((1 << 512)-1)]*251)
    details = (dict(method='all-support exact placement', support_min=6,
                    support_choices=[None]*6+['.001']*251)
               if occupancy == 1 else dict(method='full support-box CDF cover',
                    support_min=6, leaves=[dict(tilt='.001')]))
    return dict(schema='packed-gl32-sparse-hill-variant-1', K=1 << 20, N=1 << 21,
        group_count=2048, outer='BCH256128', block_rows=4, block_columns=8,
        mixing='independent uniform GL32 per group/block; no additional GF16 stage',
        routing='independent shared column shuffle per group and independent regional shuffles',
        inner=dict(t=128, s=19, updates=updates), distance='1/10', threshold=209715,
        count_refinement=refinement, joint_return_through=2, lazy_density_through=4,
        count_sha256=sparse.count_hash(counts), count_premises={'synthetic': True},
        support_min=6, tilts=['.001'], requested=[occupancy],
        results=[dict(occupancy=occupancy, passed=True, upper='not a proof input', details=details)])


class SparseVariantTests(unittest.TestCase):
    def test_explicit_variant_and_scope_validation(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)/'new.json'
            args = [[1, 2, 32], '.1', ['.00032', '.001'], 256, 46, 64, output, 3, 'ekr']
            variant.validate(*args)
            for index, bad in ((0, [0]), (0, [1, 1]), (0, [True]), (0, [33]),
                               (1, '.5'), (2, ['.001', '1/1000']), (2, ['0']),
                               (3, 64), (4, True), (5, -1), (7, True), (7, 5), (8, 'r2')):
                candidate = args[:]
                candidate[index] = bad
                with self.assertRaises(ValueError):
                    variant.validate(*candidate)
            output.touch()
            with self.assertRaises(ValueError):
                variant.validate(*args)

    def test_count_dispatch_and_schema_guard(self):
        sparse = SimpleNamespace(checked_counts=Mock(return_value=(('counts',), 6)))
        premises = dict(schema=variant.REFINEMENTS['ekr'], block_width=8,
                        block_rows=4, groups_per_outer_word=32)
        with patch.object(variant, 'intersection_counts', return_value=(('counts',), premises)) as ekr, \
                patch.object(variant, 'incidence_counts') as old:
            self.assertEqual(variant.authenticated_counts(sparse, 'ekr'), (('counts',), 6, premises))
            ekr.assert_called_once_with()
            old.assert_not_called()
        with patch.object(variant, 'intersection_counts', return_value=(('counts',), dict(premises, schema='old'))):
            with self.assertRaises(ValueError):
                variant.authenticated_counts(sparse, 'ekr')

    def test_one_group_keeps_new_support_floor_and_all_moments(self):
        sparse = variant.q1._sparse_interface()
        counts = tuple([Q(0)]*6+[Q((1 << 512)-1)]*251)
        ctx.prec = 256
        moments = [arb(1)]*257
        moments[6:] = [arb(2)**-600]*251
        with patch.object(sparse.single_group, 'support_moments', return_value=moments) as evaluate:
            upper, details = variant.one_group(sparse, counts, {('.001', '1'): ([None, None], [])}, 0, 256)
        evaluate.assert_called_once_with(None, None)
        self.assertLess(upper, arb(2)**-70)
        self.assertEqual(details['support_min'], 6)
        self.assertEqual(len(details['support_choices']), 257)
        self.assertEqual(details['support_choices'][6], '.001')
        with patch.object(sparse.single_group, 'support_moments', return_value=[arb(1)]*256):
            with self.assertRaises(ValueError):
                variant.one_group(sparse, counts, {('.001', '1'): ([None, None], [])}, 0, 256)

    def test_partial_coverage_cannot_become_full_prefix(self):
        sparse = variant.q1._sparse_interface()
        record = dict(requested=[1, 2], results=[dict(occupancy=1, passed=True),
                      dict(occupancy=2, passed=False)], target_bits=40)
        variant.update_scope(record, arb(2)**-80, sparse)
        self.assertFalse(record['complete_requested'])
        self.assertFalse(record['complete_sparse_prefix'])
        self.assertFalse(record['aggregate_target_passed'])
        self.assertIsNone(record['aggregate_upper'])
        record['results'][1]['passed'] = True
        variant.update_scope(record, arb(2)**-80, sparse)
        self.assertTrue(record['complete_requested'])
        self.assertTrue(record['aggregate_target_passed'])
        self.assertFalse(record['complete_sparse_prefix'])
        self.assertEqual(record['uncovered_sparse'], list(range(3, 33)))

    def test_aggregate_target_is_not_inferred_from_per_row_targets(self):
        sparse = variant.q1._sparse_interface()
        record = dict(requested=list(range(1, 33)), target_bits=40,
                      results=[dict(occupancy=q, passed=True) for q in range(1, 33)])
        # Every hypothetical row could be2^-41, while their sum is2^-36.
        variant.update_scope(record, arb(2)**-36, sparse)
        self.assertTrue(record['complete_sparse_prefix'])
        self.assertFalse(record['aggregate_target_passed'])

    def test_run_builds_requested_update_variant_without_legacy_claim(self):
        sparse = variant.q1._sparse_interface()
        fake_matrix = SimpleNamespace(nrows=lambda: 1)
        operators = {('.001', '1'): ([fake_matrix]*3, [])}
        counts = tuple([Q(0)]*6+[Q((1 << 512)-1)]*251)
        captured = {}
        def build(args):
            captured.update(updates=args.updates, degree=args.groups,
                            joint=args.joint_return_through, lazy=args.lazy_density_through)
            return operators
        with TemporaryDirectory() as directory, \
                patch.object(variant.q1, '_sparse_interface', return_value=sparse), \
                patch.object(variant, 'authenticated_counts', return_value=(counts, 6, {'synthetic': True})), \
                patch.object(sparse.occupancy_birth_classes, 'build_operators', side_effect=build), \
                patch.object(sparse, 'save'), \
                patch.object(variant, 'one_group', return_value=(arb(2)**-80, {'synthetic': True})), \
                patch.object(variant, '_support_cover', return_value=(arb(2)**-80, {'synthetic': True})) as cover, \
                redirect_stdout(io.StringIO()):
            record = variant.run([1, 2], '.1', ['.001'], 256, 46, 0,
                                 Path(directory)/'new.json', updates=3, refinement='ekr')
        self.assertEqual(captured, dict(updates=3, degree=2, joint=2, lazy=4))
        self.assertEqual(record['inner']['updates'], 3)
        self.assertEqual(record['schema'], 'packed-gl32-sparse-hill-variant-1')
        self.assertEqual(record['threshold'], 209715)
        self.assertEqual(cover.call_args.kwargs['support_min'], 6)
        self.assertTrue(record['complete_requested'])
        self.assertFalse(record['complete_sparse_prefix'])
        with self.assertRaises(ValueError):
            sparse.validate_record(record)

    def test_replay_record_scope_and_witness_guards(self):
        sparse = variant.q1._sparse_interface()
        record = synthetic_search_record(sparse)
        self.assertEqual(len(variant.validate_record(record)), 1)
        for key, bad in (('schema', 'packed-gl32-sparse-1'), ('threshold', 209716),
                         ('inner', dict(t=128, s=19, updates=True)),
                         ('joint_return_through', 3), ('requested', [1, 1]),
                         ('count_refinement', 'unknown')):
            with self.assertRaises(ValueError):
                variant.validate_record(dict(record, **{key: bad}))
        tampered = json.loads(json.dumps(record))
        tampered['results'][0]['details']['support_choices'][6] = '.002'
        with self.assertRaises(ValueError):
            variant.validate_record(tampered)

    def test_replay_fresh_r4_operators_and_endpoints_not_saved_numbers(self):
        sparse = variant.q1._sparse_interface()
        counts = tuple([Q(0)]*6+[Q((1 << 512)-1)]*251)
        record = synthetic_search_record(sparse)
        matrix = SimpleNamespace(nrows=lambda: 1)
        captured = {}
        def build(args):
            captured.update(updates=args.updates, precision=args.precision, degree=args.groups)
            return {('.001', '1'): ([matrix]*2, [])}
        with TemporaryDirectory() as directory:
            source = Path(directory)/'source.json'
            source.write_text(json.dumps(record))
            with patch.object(variant.q1, '_sparse_interface', return_value=sparse), \
                    patch.object(variant, 'authenticated_counts', return_value=(counts, 6, {'synthetic': True})) as auth, \
                    patch.object(sparse.occupancy_birth_classes, 'build_operators', side_effect=build), \
                    patch.object(variant, 'one_group', return_value=(arb(2)**-80, {})) as one, \
                    patch.object(sparse, 'save'), redirect_stdout(io.StringIO()):
                replay = variant.replay_records([source], Path(directory)/'replay.json',
                    precision=384, target_bits=46, updates=4, refinement='ekr')
            self.assertEqual(captured, dict(updates=4, precision=384, degree=1))
            auth.assert_called_once_with(sparse, 'ekr')
            one.assert_called_once()
            self.assertEqual(replay['scope']['inner']['updates'], 4)
            self.assertEqual(replay['results'][0]['upper'], sparse.endpoint(arb(2)**-80))
            self.assertEqual(replay['sources'][0]['path'], str(source.resolve()))
            self.assertEqual(len(replay['sources'][0]['sha256']), 64)
            self.assertTrue(replay['complete_requested'])
            self.assertFalse(replay['complete_sparse_prefix'])
            self.assertEqual(replay['uncovered_sparse'], list(range(2, 33)))

    def test_replay_rejects_mixed_variants_duplicate_scopes_and_count_mismatch(self):
        sparse = variant.q1._sparse_interface()
        counts = tuple([Q(0)]*6+[Q((1 << 512)-1)]*251)
        with TemporaryDirectory() as directory:
            first, second = Path(directory)/'first.json', Path(directory)/'second.json'
            first.write_text(json.dumps(synthetic_search_record(sparse)))
            cases = [synthetic_search_record(sparse, 2, updates=3),
                     synthetic_search_record(sparse, 2, refinement='incidence'),
                     synthetic_search_record(sparse)]
            for record in cases:
                second.write_text(json.dumps(record))
                with self.assertRaises(ValueError):
                    variant.replay_records([first, second], Path(directory)/'new.json')
            with patch.object(variant.q1, '_sparse_interface', return_value=sparse), \
                    patch.object(variant, 'authenticated_counts', return_value=(counts, 6, {'different': True})), \
                    self.assertRaises(ValueError):
                variant.replay_records([first], Path(directory)/'new.json')
            with self.assertRaises(ValueError):
                variant.replay_records([first], Path(directory)/'new.json', updates=3)
            with self.assertRaises(ValueError):
                variant.replay_records([first], Path(directory)/'new.json', distance='.11')

    def test_changed_source_cannot_emit_complete_replay(self):
        sparse = variant.q1._sparse_interface()
        counts = tuple([Q(0)]*6+[Q((1 << 512)-1)]*251)
        matrix = SimpleNamespace(nrows=lambda: 1)
        with TemporaryDirectory() as directory:
            source = Path(directory)/'source.json'
            source.write_text(json.dumps(synthetic_search_record(sparse)))
            def changed(*args):
                source.write_text(source.read_text()+'\n')
                return arb(2)**-80, {}
            with patch.object(variant.q1, '_sparse_interface', return_value=sparse), \
                    patch.object(variant, 'authenticated_counts', return_value=(counts, 6, {'synthetic': True})), \
                    patch.object(sparse.occupancy_birth_classes, 'build_operators',
                                 return_value={('.001', '1'): ([matrix]*2, [])}), \
                    patch.object(variant, 'one_group', side_effect=changed), \
                    patch.object(sparse, 'save') as save, redirect_stdout(io.StringIO()), \
                    self.assertRaises(ValueError):
                variant.replay_records([source], Path(directory)/'new.json')
            self.assertFalse(save.call_args.args[1]['complete_requested'])


if __name__ == '__main__':
    unittest.main()
