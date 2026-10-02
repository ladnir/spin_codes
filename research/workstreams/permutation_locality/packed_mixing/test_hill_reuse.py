import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb, ctx
import hill_reuse as reuse
from test_hill_cover import scope
from test_point_reuse import FakeModel, witness


def point_record(updates=4):
    return dict(schema=reuse.hill_cover.POINT_SCHEMA, scope=scope(updates, True),
        points=[dict(mean='1/8', cell=['1/8', '1/8'],
            checked=dict(witness=witness(), proposal='NOT USED', upper='NOT USED'))])


class HillReuseTests(unittest.TestCase):
    def setUp(self):
        previous = ctx.prec
        ctx.prec = 128
        self.addCleanup(setattr, ctx, 'prec', previous)

    def test_actual_updates_explicit_and_default_remains_r2(self):
        hints = [dict(mean='1/8', witness=witness())]
        for updates in (2, 3, 4):
            model = FakeModel()
            model.data = dict(windows=32, updates=updates)
            reused = reuse.PointReuseModel(model, hints, 40, updates=updates)
            self.assertEqual(reused.data['updates'], updates)
            if updates != 2:
                with self.assertRaises(ValueError):
                    reuse.PointReuseModel(model, hints, 40)
            with self.assertRaises(ValueError):
                reuse.PointReuseModel(model, hints, 40, updates=3 if updates != 3 else 4)
        for updates in (True, 1, 5, '4'):
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                reuse.PointReuseModel(FakeModel(), hints, 40, updates=updates)

    def test_points_ignore_scores_but_require_exact_canonical_scope(self):
        record = point_record()
        target = copy.deepcopy(record['scope'])
        expected = reuse.proposals(record, target)
        record['scope'].update(cover='IGNORED', mixture_verification='IGNORED')
        record['points'][0]['checked'].update(upper=[1, 500000], proposal=-1e99)
        self.assertEqual(reuse.proposals(record, target), expected)
        self.assertEqual(set(expected[0]), {'mean', 'witness'})
        for key, value in [('updates', 2), ('distance', '9/100'),
                ('mixture', [dict(mass='2', activity='1/2')]), ('variance_bins', 64),
                ('comparison_caps_sha256', 'c'*64)]:
            changed = copy.deepcopy(record)
            changed['scope'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                reuse.proposals(changed, target)
        record['points'][0]['checked']['witness']['parameters'][0] = '1/9'
        self.assertEqual(expected[0]['witness']['parameters'][0], '1/20')

    def test_interval_source_rebased_from_original_geometry(self):
        target = scope(4, True)
        paths = [format(i, '04b') for i in range(16)]
        cells = {p: dict(cell=list(map(str, reuse.hill_cover.geometry.path_cell(target['root'], p))))
                 for p in paths}
        candidate = witness()
        part = dict(interval=['0', '1/8'], dual=['0', '0', '0'])
        candidate.update(variance_partition=[part],
            regional_count_parts=[dict(copy.deepcopy(part), mgf_witnesses=[])])
        row = cells.pop('0001')
        row.update(witness=candidate, upper='IGNORED', proposal='IGNORED')
        record = dict(schema=reuse.hill_cover.SCHEMA, scope=target, precision=256,
            cell_target_bits=52, whole_code_certificate=False,
            cover=dict(leaves={'0001': row}, unresolved=cells, visited=0))
        original = copy.deepcopy(record)
        points = reuse.proposals(record, target)
        self.assertEqual(points[0]['mean'], '3/32')
        self.assertEqual(points[0]['witness']['variance_partition'][0]['interval'], ['0', '3/32'])
        self.assertEqual(points[0]['witness']['regional_count_parts'][0]['interval'], ['0', '3/32'])
        self.assertEqual(record, original)
        self.assertEqual(set(points[0]['witness']), set(candidate))

    def test_unknown_bound_fields_in_witness_and_malformed_parameters_rejected(self):
        for change in (lambda w: w.update(upper=[1, -500]),
                lambda w: w.update(parameters=['-1', '0', '0']),
                lambda w: w.update(parameters=['1', '0', '-1']),
                lambda w: w.update(parameters=[.1, '0', '0']),
                lambda w: w.update(variance_dual=['0']),
                lambda w: w.update(regional_feedback_uniform_classes=1),
                lambda w: w.update(regional_count_parts=[])):
            record = point_record()
            change(record['points'][0]['checked']['witness'])
            with self.subTest(change=change), self.assertRaises(ValueError):
                reuse.proposals(record, record['scope'])
        record = point_record()
        record['points'][0]['cell'] = ['1/8', '1/4']
        with self.assertRaises(ValueError):
            reuse.proposals(record, record['scope'])

    def test_source_hashes_and_mutation_guard(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'point.json'
            record = point_record()
            raw = json.dumps(record).encode()
            path.write_bytes(raw)
            points, sources = reuse.load_hints([path], record['scope'])
            self.assertEqual(len(points), 1)
            self.assertEqual(sources[0]['sha256'], hashlib.sha256(raw).hexdigest())
            self.assertEqual(sources[0]['hint_count'], 1)
            reuse.check_sources(sources)
            with self.assertRaises(ValueError):
                reuse.load_hints([path, path], record['scope'])
            path.write_bytes(raw+b'\n')
            with self.assertRaises(ValueError):
                reuse.check_sources(sources)

    def test_factory_fresh_checks_cache_and_scope_do_not_change_variant(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'point.json'
            record = point_record()
            path.write_text(json.dumps(record))
            model = FakeModel()
            model.root = tuple(map(Q, record['scope']['root']))
            model.threshold = record['scope']['threshold']
            model.data = dict(windows=32, updates=4)
            original = copy.deepcopy(record['scope'])
            proxy, cache, metadata = reuse.wrap(model, original, [path], 40)
            cell = (Q(1, 8), Q(3, 16))
            score, candidate = proxy.proposal(cell)
            self.assertLess(score, -42)
            self.assertTrue(proxy.outward(cell, candidate) < arb(2)**-40)
            self.assertEqual(len(model.outward_calls), 1)
            self.assertEqual(proxy.data['updates'], 4)
            self.assertEqual(proxy.root, model.root)
            self.assertEqual(original, record['scope'])
            self.assertIsInstance(cache, reuse.RegionalCache)
            self.assertEqual(metadata[0]['hint_count'], 1)
            for field, value in [('threshold', 1), ('variance_bins', 64), ('q_min', 34)]:
                wrong = copy.copy(model)
                setattr(wrong, field, value)
                with self.subTest(field=field), self.assertRaises(ValueError):
                    reuse.wrap(wrong, original, [path], 40)


if __name__ == '__main__':
    unittest.main()
