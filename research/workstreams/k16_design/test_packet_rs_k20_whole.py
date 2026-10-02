"""Lightweight fail-closed checks for K20 union assembly and fresh replay.

Synthetic endpoints test assembly only. No selected-map census, encoder run,
or numerical distance computation is performed by this test module.
"""
from contextlib import ExitStack, redirect_stdout
from copy import deepcopy
from dataclasses import asdict
from fractions import Fraction as Q
from io import StringIO
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from flint import arb, ctx
import packet_rs_k20_whole as whole


def exact_endpoint(pair):
    return Q(pair[0]) * Q(2) ** pair[1]


class K20WholeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beta, cls.counts = whole.sparse.exact_outer()

    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 256
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.sources = whole.q1.source_snapshot()
        self.map_record = {'synthetic_test_map': 'not a numerical certificate'}

    def tearDown(self):
        ctx.prec = self.precision

    def records(self):
        common = dict(outer='rs16', K=2**20, N=2**21,
            geometry=asdict(whole.GEOMETRY), threshold=209715, distance='1/10',
            precision=256, zero_initial_state=True, final_flush=False,
            whole_code_certificate=False, source_sha256=self.sources,
            map_record=self.map_record)
        one = dict(common, schema='finite-packet-q1-diagnostic-1',
            occupancy_covered=[1], count_kind='shells', q1_upper=[1, -64],
            count_sha256=whole.sparse.count_hash(self.counts, cumulative=True))
        two = dict(common, schema='finite-packet-q2-diagnostic-1',
            occupancy_covered=[2], count_kind='exact_expected_shells', q2_upper=[1, -64],
            all_two_group_support_pairs_covered=True,
            count_sha256=whole.sparse.count_hash(self.counts))
        tails = []
        for schema, first, last in ((whole.sparse.TAIL_SCHEMA, 3, 32),
                                    (whole.dense.SCHEMA, 33, 8192)):
            tails.append(dict(common, schema=schema, q_min=first, q_max=last,
                occupancy_covered=[first, last], every_shell_checked=True,
                state_continuity=whole.CONTINUITY, beta=str(self.beta),
                count_sha256=whole.sparse.count_hash(self.counts),
                occupancy_uppers={str(q): [1, -64] for q in range(first, last + 1)}))
        return deepcopy([one, two, *tails])

    def save(self, records):
        paths = [self.folder / f'component-{i}.json' for i in range(len(records))]
        for path, record in zip(paths, records):
            path.write_text(json.dumps(record), encoding='utf-8')
        return paths

    def combine(self, records):
        paths = self.save(records)
        return whole.combine(paths[0], paths[1], paths[2:])

    def test_coverage_requires_each_integer_once_including_endpoints(self):
        whole.check_coverage([(1, 1), (2, 2), (3, 32), (33, 8192)])
        whole.check_coverage([(1, 8192)])
        for ranges in ([], [(1, 8191)], [(2, 8192)], [(1, 32), (34, 8192)],
                       [(1, 32), (32, 8192)], [(0, 8192)], [(1, 8193)],
                       [(1, 32), (33, 33.0), (34, 8192)], [(True, 8192)]):
            with self.subTest(ranges=ranges), self.assertRaises(ValueError):
                whole.check_coverage(ranges)

    def test_saved_union_never_claims_fresh_certificate(self):
        result = self.combine(self.records())
        self.assertEqual(result['occupancy_covered'], [1, 8192])
        self.assertTrue(result['all_occupancies_covered'])
        self.assertTrue(result['target_met'])
        self.assertFalse(result['fresh_replay'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['threshold'], 209715)
        self.assertEqual(result['minimum_component_precision'], 256)
        self.assertGreaterEqual(exact_endpoint(result['union_upper']), Q(8192, 2**64))
        self.assertLess(exact_endpoint(result['union_upper']), Q(1, 2**50))
        self.assertEqual(len(result['input_receipts']), 4)

    def test_union_recomputes_sum_and_rejects_uncovered_or_duplicate_tail(self):
        records = self.records()
        # A saved summary is not evidence: all endpoints must be added afresh.
        records[3]['union_upper'] = [1, -99999]
        records[0]['q1_upper'] = [1, -40]
        result = self.combine(records)
        self.assertFalse(result['target_met'])
        self.assertFalse(result['whole_code_certificate'])
        paths = self.save(self.records())
        for tails in ([paths[2]], [paths[3]], [paths[2], paths[2], paths[3]]):
            with self.assertRaises(ValueError):
                whole.combine(paths[0], paths[1], tails)

    def test_receipt_scope_and_source_pins_fail_closed(self):
        cases = [('K', 65536), ('N', 131072), ('threshold', 13107),
                 ('distance', '1/8'), ('outer', 'rs8'), ('zero_initial_state', False),
                 ('final_flush', True), ('precision', 127), ('precision', True),
                 ('whole_code_certificate', True), ('map_record', None)]
        for key, value in cases:
            records = self.records()
            records[0][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                self.combine(records)
        for required in (str(Path(whole.q1.__file__).resolve()),
                         str(whole.q1.kernel_t64.SELECTED_MAP.resolve())):
            for operation in ('omit', 'change'):
                record = self.records()[0]
                if operation == 'omit':
                    del record['source_sha256'][required]
                else:
                    record['source_sha256'][required] = '0' * 64
                path = self.save([record])[0]
                with self.subTest(required=required, operation=operation), self.assertRaises(ValueError):
                    whole.load_receipt(path)

    def test_count_map_and_tail_metadata_fail_closed(self):
        cases = [(0, 'schema', 'wrong'), (0, 'occupancy_covered', [2]),
                 (0, 'count_kind', 'cumulative_caps'), (0, 'count_sha256', 'wrong'),
                 (1, 'all_two_group_support_pairs_covered', False),
                 (1, 'count_kind', 'shells'), (1, 'count_sha256', 'wrong'),
                 (1, 'map_record', {'other_map': 1}),
                 (2, 'schema', 'wrong'), (2, 'map_record', {'other_map': 1}),
                 (2, 'every_shell_checked', False), (2, 'state_continuity', 'reset'),
                 (2, 'count_sha256', 'wrong'), (2, 'beta', '1'),
                 (2, 'occupancy_covered', [3, 31]), (2, 'q_min', True)]
        for component, key, value in cases:
            records = self.records()
            records[component][key] = value
            with self.subTest(component=component, key=key), self.assertRaises(ValueError):
                self.combine(records)

    def test_core_and_schema_producer_pins_cannot_be_stripped(self):
        here = Path(whole.__file__).resolve().parent
        cases = [(0, here / 'rs_outer.py'), (0, here / 'rs_uniform_envelope.py'),
                 (0, whole.q1.LOCALITY / 'packed_mixing/s16_closure/kernel_t64.py'),
                 (0, whole.q1.LOCALITY / 'packed_mixing/s16_closure/kernel_birth_density.py'),
                 (1, here / 'packet_q2.py'), (2, here / 'packet_rs_k20_sparse.py'),
                 (2, here / 'packet_uniform_tail.py'), (3, here / 'packet_rs_k20_dense.py')]
        for component, missing in cases:
            records = self.records()
            del records[component]['source_sha256'][str(missing.resolve())]
            with self.subTest(component=component, missing=missing), self.assertRaises(ValueError):
                self.combine(records)
        # Audit the production import graph in isolation. Earlier tests may
        # legitimately load optional toy-map helpers into this interpreter;
        # those must not alter the expected default import graph.
        probe = ('import json; from pathlib import Path; '
                 'import packet_rs_k20_whole as w; '
                 'print(json.dumps(sorted(Path(p).relative_to(w.q1.LOCALITY).as_posix() '
                 'for p in w.q1.source_snapshot() '
                 'if w.q1.LOCALITY in Path(p).parents)))')
        core = set(json.loads(subprocess.check_output(
            [sys.executable, '-c', probe], cwd=here, text=True)))
        self.assertEqual(core, set(whole.CORE_FILENAMES))
        for record in self.records():
            required = whole.required_source_paths(record['schema'])
            self.assertTrue(required <= self.sources.keys())
            self.assertNotIn(str(Path(whole.__file__).resolve()), required)
            self.assertNotIn(str(Path(__file__).resolve()), required)
        with self.assertRaises(ValueError):
            whole.required_source_paths('unknown producer')

    def test_fugacity_producer_uses_same_full_source_validation(self):
        records = self.records()
        records[3]['schema'] = whole.dense.FUGACITY_SCHEMA
        self.assertTrue(self.combine(records)['all_occupancies_covered'])
        del records[3]['source_sha256'][str(Path(whole.dense.__file__).resolve())]
        with self.assertRaises(ValueError):
            self.combine(records)

    def test_tail_endpoint_keys_cannot_omit_or_invent_occupancies(self):
        for key, operation in (('32', 'omit'), ('8192', 'add')):
            records = self.records()
            if operation == 'omit':
                del records[2]['occupancy_uppers'][key]
            else:
                records[2]['occupancy_uppers'][key] = [1, -64]
            with self.subTest(endpoint=key, operation=operation), self.assertRaises(ValueError):
                self.combine(records)

    def test_dyadic_endpoint_validation_and_outward_conversion(self):
        for invalid in ([0, 0], [-1, 0], [True, 0], [1, True], [1, 10_000_001],
                        (1, 0), [1], ['1', 0]):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                whole.endpoint_value(invalid)
        for endpoint in ([1, -64], [2**400 + 1, -500], [1, -999999]):
            result = whole.endpoint_value(endpoint)
            self.assertTrue(result.is_finite())
            self.assertTrue(result > 0)
            # Endpoint serialization must retain an upper bound even if the
            # input dyadic mantissa exceeds the working Arb precision.
            actual = exact_endpoint(whole.q1.endpoint(result))
            self.assertGreaterEqual(actual, exact_endpoint(endpoint))

    def test_recipe_validates_complete_partition_and_positive_witnesses(self):
        for split in (3, 32, 8191):
            plan = whole.recipe(sparse_max=split)
            whole.check_coverage(plan['intervals'])
            self.assertEqual(plan['precision'], 256)
        for kwargs in (dict(sparse_max=2), dict(sparse_max=8192), dict(sparse_max=True),
                       dict(precision=192), dict(precision=True), dict(q1_tilts=[]),
                       dict(q2_tilts=['0']), dict(sparse_tilts=['.5', '1/2']),
                       dict(dense_tilts='0.1'), dict(marker_probabilities=[]),
                       dict(marker_probabilities=['1']), dict(marker_probabilities=['0']),
                       dict(marker_probabilities=['.5', '1/2']), dict(marker_probabilities='.5')):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                whole.recipe(**kwargs)

    def _mock_replay(self, output, *, target_met=True, bad_map=False, changed_source=False):
        result = dict(target_met=target_met, map_record=self.map_record,
                      margin_bits='synthetic', tail_margin_bits='synthetic',
                      fresh_replay=False, whole_code_certificate=False)
        if bad_map:
            result['map_record'] = {'different_map': 1}
        with ExitStack() as stack:
            prepared = stack.enter_context(patch.object(whole.q1.kernel_t64, 'prepare',
                return_value=({'fresh_test_data': 1}, self.map_record)))
            sources = [self.sources, {'changed': 'source'}] if changed_source else None
            stack.enter_context(patch.object(whole.q1, 'source_snapshot',
                side_effect=sources, return_value=self.sources))
            runners = [stack.enter_context(patch.object(module, name)) for module, name in
                       ((whole.sparse, 'run_q1'), (whole.sparse, 'run_q2'),
                        (whole.sparse, 'run_tail'), (whole.dense, 'run_fugacity'))]
            combined = stack.enter_context(patch.object(whole, 'combine', return_value=result))
            stack.enter_context(redirect_stdout(StringIO()))
            answer = whole.replay(output)
        return answer, prepared, runners, combined

    def test_replay_recomputes_all_components_and_only_marks_passing_union(self):
        for passes in (False, True):
            output = self.folder / f'whole-{passes}.json'
            result, prepared, runners, combined = self._mock_replay(output, target_met=passes)
            prepared.assert_called_once_with(birth_density='capped')
            for runner in runners:
                runner.assert_called_once()
                self.assertEqual(runner.call_args.kwargs['precision'], 256)
                self.assertEqual(runner.call_args.kwargs['data'], {'fresh_test_data': 1})
                self.assertEqual(runner.call_args.kwargs['map_record'], self.map_record)
            self.assertEqual(runners[2].call_args.kwargs['q_min'], 3)
            self.assertEqual(runners[2].call_args.kwargs['q_max'], 32)
            self.assertEqual(runners[3].call_args.kwargs['q_min'], 33)
            self.assertEqual(runners[3].call_args.kwargs['q_max'], 8192)
            combined.assert_called_once()
            self.assertTrue(result['fresh_replay'])
            self.assertIs(result['whole_code_certificate'], passes)
            self.assertEqual(json.loads(output.read_text()), result)

    def test_replay_rejects_existing_paths_before_prepare(self):
        output = self.folder / 'whole.json'
        for occupied in (output, output.with_name('whole-q1.json'),
                         output.with_name('whole-dense.json')):
            occupied.touch()
            with patch.object(whole.q1.kernel_t64, 'prepare') as prepare, self.assertRaises(ValueError):
                whole.replay(output)
            prepare.assert_not_called()
            occupied.unlink()

    def test_replay_rejects_changed_source_or_map_without_final_receipt(self):
        for kwargs in (dict(bad_map=True), dict(changed_source=True)):
            output = self.folder / 'must-not-exist.json'
            with self.subTest(kwargs=kwargs), self.assertRaises(RuntimeError):
                self._mock_replay(output, **kwargs)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
