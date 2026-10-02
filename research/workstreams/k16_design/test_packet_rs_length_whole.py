"""Light fail-closed tests; synthetic endpoints never run a map census."""
from contextlib import ExitStack, redirect_stdout
from copy import deepcopy
from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
from io import StringIO
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from flint import ctx
import packet_rs_length_whole as whole
# Import only the helper's definitions, matching the completed census source
# graph without running its numerical census in these lightweight tests.
import feedback_exact


def fraction(pair):
    return Q(pair[0]) * Q(2) ** pair[1]


def declared_record():
    """Identity-only fixture: no saved spectrum is used as a numerical premise."""
    path = whole.q1.kernel_t64.SELECTED_MAP.resolve()
    raw = path.read_bytes()
    declaration = json.loads(raw)
    rows = [int(x, 16) for x in declaration['generator_rows_hex']]
    columns = declaration['columns']
    identity = dict(bits=16, width=64, expansion_rows=list(map(hex, rows)),
                    feedback_columns=columns)
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return dict(schema='s16-selected-t64-fixed-maps-1', t=64, s=16,
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=columns,
        source=dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest()),
        distribution='uniform_gl', sampling='independent uniform GL16 for every physical t64 step',
        map_sha256=digest, whole_code_certificate=False,
        macro=dict(physical_steps=2, physical_step_bits=64, macro_step_bits=128,
            physical_windows=16, macro_windows=32, state_continuity='retained_between_halves',
            initial_state='zero', flush=False))


class LengthWholeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beta, cls.counts = whole.sparse.exact_outer()
        cls.map_record = declared_record()
        cls.sources = whole.q1.source_snapshot()

    def setUp(self):
        self.old_precision = ctx.prec
        ctx.prec = 256
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def tearDown(self):
        ctx.prec = self.old_precision

    def records(self, K=4096):
        geom = whole.lengths.geometry(K)
        common = dict(outer='rs16', K=K, N=2*K, geometry=asdict(geom),
            threshold=2*K//10, distance='1/10', precision=256,
            groups=geom.group_count, regions=64, group_dimension=128, group_output_bits=256,
            physical_t=64, state_bits=16, macro_t=128, physical_steps=2*K//64,
            macro_steps=2*K//128, physical_steps_per_macro=2,
            macros_per_region=geom.macros_per_region, zero_initial_state=True,
            final_flush=False, whole_code_certificate=False, fresh_computation=True,
            source_sha256=self.sources, map_record=self.map_record)
        records = []
        for q in (1, 2):
            record = dict(common, schema=whole.LOW_SCHEMAS[q-1], length_component=f'q{q}',
                occupancy_covered=[q], count_kind='shells' if q == 1 else 'exact_expected_shells',
                count_sha256=whole.sparse.count_hash(self.counts, cumulative=q == 1))
            record[f'q{q}_upper'] = [1, -64]
            if q == 2:
                record['all_two_group_support_pairs_covered'] = True
            records.append(record)
        tail = dict(common, state_continuity=whole.CONTINUITY, every_shell_checked=True,
                    beta=str(self.beta), count_sha256=whole.sparse.count_hash(self.counts))
        exact = dict(tail, schema=whole.lengths.TAIL_SCHEMA, length_component='tail',
            occupancy_covered=[3, 5, 7],
            occupancy_uppers={str(q): [1, -65] for q in (3, 5, 7)})
        # Match the producer: geometry is complete, redundant aliases absent.
        for key in ('groups', 'regions', 'group_dimension', 'group_output_bits'):
            del exact[key]
        dense = dict(tail, schema=whole.dense.FUGACITY_SCHEMA,
            q_min=3, q_max=geom.group_count, occupancy_covered=[3, geom.group_count],
            evaluated_every_integer_occupancy=True,
            occupancy_uppers={str(q): [1, -64] for q in range(3, geom.group_count+1)})
        return deepcopy([*records, exact, dense])

    def assemble(self, records, K=4096):
        return whole.assemble(records, K)

    def test_exact_dyadic_comparison_sum_and_final_upward_rounding(self):
        rng = random.Random(7)
        for _ in range(40):
            pairs = [[rng.randrange(1, 1 << 40), rng.randrange(-600, 500)] for _ in range(8)]
            normalized = list(map(whole.dyadic, pairs))
            expected = sum(map(fraction, pairs), Q(0))
            exact = whole.dyadic_sum(normalized)
            self.assertEqual(fraction(exact), expected)
            self.assertEqual(whole.dyadic_less(normalized[0], normalized[1]),
                             fraction(pairs[0]) < fraction(pairs[1]))
            rounded = whole.rounded_endpoint(exact, 128)
            self.assertGreaterEqual(fraction(rounded), expected)
            self.assertLessEqual(rounded[0].bit_length(), 128)
        self.assertEqual(whole.dyadic([12, -10]), (3, -8))
        self.assertFalse(whole.dyadic_less((3, -8), (6, -9)))
        for bad in ([0, 0], [-1, 0], [True, 0], [1, True], [1, 10_000_001],
                    [1], [1, '0'], (1, -64), None):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                whole.dyadic(bad)

    def test_saved_union_minimizes_overlaps_without_claiming_freshness(self):
        result = self.assemble(self.records())
        self.assertEqual(result['occupancy_covered'], [1, 32])
        self.assertTrue(result['all_occupancies_covered'])
        self.assertTrue(result['target_met'])
        self.assertFalse(result['fresh_replay'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(fraction(result['union_upper']), Q(61, 2**65))
        self.assertEqual(result['selected_components']['3'], 2)
        self.assertEqual(result['selected_components']['4'], 3)
        self.assertEqual(result['threshold'], 819)
        self.assertEqual(result['target_minimum_distance'], 820)
        self.assertNotIn('guaranteed_minimum_distance', result)
        self.assertEqual(len(result['occupancy_uppers']), 32)
        reversed_result = self.assemble(list(reversed(self.records())))
        self.assertEqual(result['union_upper'], reversed_result['union_upper'])

    def test_exact_sum_not_saved_subtotal_or_floating_threshold_drives_acceptance(self):
        records = self.records()
        records[-1]['union_upper'] = [1, -99999]
        records[0]['q1_upper'] = [1, -40]
        result = self.assemble(records)
        self.assertFalse(result['target_met'])
        self.assertGreater(fraction(result['union_upper']), Q(1, 2**40))
        # Even an addend too tiny for the mantissa must round upward, not vanish.
        exact = whole.dyadic_sum([(1, -40), (1, -2000)])
        upper = whole.rounded_endpoint(exact, 128)
        self.assertGreater(fraction(upper), Q(1, 2**40))

    def test_tail_gaps_and_endpoint_list_mismatches_fail_closed(self):
        cases = []
        records = self.records()
        del records[-1]['occupancy_uppers']['32']
        cases.append(records)
        records = self.records()
        records[-1].update(q_max=31, occupancy_covered=[3, 31])
        del records[-1]['occupancy_uppers']['32']
        cases.append(records)
        records = self.records()
        records[2]['occupancy_covered'] = [3, 5]
        cases.append(records)
        records = self.records()
        records[2]['occupancy_uppers']['6'] = [1, -64]
        cases.append(records)
        records = self.records()
        records[2]['occupancy_covered'] = [7, 5, 3]
        cases.append(records)
        cases.extend((self.records()[1:], self.records()[:1] + self.records()[2:]))
        for i, records in enumerate(cases):
            with self.subTest(i=i), self.assertRaises(ValueError):
                self.assemble(records)
        self.assertTrue(self.assemble(self.records() + [self.records()[2]])['target_met'])
        with self.assertRaises(ValueError):
            self.assemble(self.records() + [self.records()[0]])

    def test_scope_geometry_counts_and_state16_fail_closed(self):
        cases = [(0, 'K', 8192), (0, 'N', 16384), (0, 'threshold', 820),
            (0, 'distance', '19/200'), (0, 'state_bits', 20), (0, 'physical_t', 128),
            (0, 'outer', 'rs8'), (0, 'zero_initial_state', False), (0, 'final_flush', True),
            (0, 'whole_code_certificate', True), (0, 'fresh_computation', False),
            (0, 'precision', True), (0, 'precision', 127),
            (0, 'count_kind', 'cdf'), (0, 'count_sha256', 'wrong'),
            (0, 'length_component', 'q2'), (1, 'all_two_group_support_pairs_covered', False),
            (2, 'every_shell_checked', False), (2, 'state_continuity', 'reset'),
            (2, 'beta', '1'), (2, 'beta', None), (2, 'count_sha256', 'wrong'),
            (3, 'evaluated_every_integer_occupancy', False), (3, 'q_min', True),
            (3, 'schema', 'unrecognized'), (2, 'regions', 256)]
        for index, key, value in cases:
            records = self.records()
            records[index][key] = value
            with self.subTest(index=index, key=key), self.assertRaises(ValueError):
                self.assemble(records)
        records = self.records()
        records[0]['geometry']['group_count'] = 64
        with self.assertRaises(ValueError):
            self.assemble(records)

    def test_current_selected_map_and_complete_record_equality_are_required(self):
        for key, value in (('s', 20), ('map_sha256', 'wrong'), ('distribution', 'uniform_field'),
                           ('sampling', 'reused update'), ('schema', 's20-map')):
            records = self.records()
            records[0]['map_record'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.assemble(records)
        records = self.records()
        records[1]['map_record'] = dict(records[1]['map_record'], extra_unchecked_metadata=True)
        with self.assertRaises(ValueError):
            self.assemble(records)

    def test_current_core_and_each_actual_producer_pin_are_required(self):
        here = Path(whole.__file__).resolve().parent
        cases = [(0, whole.q1.kernel_t64.SELECTED_MAP),
            (0, whole.q1.LOCALITY / 'packed_mixing/s16_closure/kernel_birth_density.py'),
            (0, whole.q1.LOCALITY / 'gf16_packets/feedback_exact.py'),
            (0, here / 'packet_rs_lengths.py'), (1, here / 'packet_q2.py'),
            (2, here / 'packet_uniform_tail.py'), (2, here / 'packet_rs_k20_sparse.py'),
            (2, here / 'rs_outer.py'), (3, here / 'packet_rs_k20_dense.py'),
            (3, here / 'rs_uniform_envelope.py')]
        for index, path in cases:
            for change in ('omit', 'stale'):
                records = self.records()
                # Break deepcopy's shared alias so only one component changes.
                records[index]['source_sha256'] = dict(records[index]['source_sha256'])
                filename = str(path.resolve())
                if change == 'omit':
                    del records[index]['source_sha256'][filename]
                else:
                    records[index]['source_sha256'][filename] = '0'*64
                with self.subTest(index=index, path=path, change=change), self.assertRaises(ValueError):
                    self.assemble(records)
        for record in self.records():
            self.assertTrue(whole.required_source_paths(record['schema']) <= self.sources.keys())
            self.assertNotIn(str(Path(__file__).resolve()), whole.required_source_paths(record['schema']))

    def test_standalone_producers_supply_the_required_source_graph(self):
        # Importing the coordinator itself must not manufacture a requirement
        # absent from actual standalone component production.
        for module, schemas in (
                ('packet_rs_lengths', (*whole.LOW_SCHEMAS, whole.lengths.TAIL_SCHEMA)),
                ('packet_rs_k20_dense', (whole.dense.SCHEMA, whole.dense.FUGACITY_SCHEMA))):
            probe = ('import json; import ' + module + ' as m; import feedback_exact; '
                     'print(json.dumps(sorted(m.source_snapshot() if hasattr(m,"source_snapshot") '
                     'else m.q1.source_snapshot())))')
            sources = set(json.loads(subprocess.check_output(
                [sys.executable, '-B', '-c', probe], cwd=Path(whole.__file__).parent, text=True)))
            for schema in schemas:
                self.assertTrue(whole.required_source_paths(schema) <= sources, schema)

    def test_old_dense_schema_and_saved_path_inputs_remain_supported(self):
        records = self.records()
        records[-1]['schema'] = whole.dense.SCHEMA
        del records[-1]['fresh_computation']
        paths = []
        for i, record in enumerate(records):
            path = self.folder / f'part-{i}.json'
            path.write_text(json.dumps(record), encoding='utf-8')
            paths.append(path)
        result = self.assemble(paths)
        self.assertEqual(result['input_receipts'][0]['path'], str(paths[0].resolve()))
        self.assertFalse(result['whole_code_certificate'])

    def test_generic_natural_lengths_and_invalid_options(self):
        for K in (4096, 65536, 98304, 131072):
            plan = whole.recipe(K)
            self.assertEqual(plan['geometry']['group_count'], K//128)
            self.assertEqual(plan['threshold'], 2*K//10)
        for K in (0, 1, 65537, True, 65536.0):
            with self.subTest(K=K), self.assertRaises(ValueError):
                whole.recipe(K)
        cases = [dict(precision=192), dict(precision=True), dict(exact_occupancies=[3, 3]),
            dict(exact_occupancies=[2]), dict(exact_occupancies='3'),
            dict(exact_occupancies=[], dense_min=4), dict(dense_min=2), dict(dense_min=True),
            dict(q1_tilts=[]), dict(q2_tilts=['0']), dict(exact_tilts=['.5', '1/2']),
            dict(dense_tilts='.1'), dict(marker_probabilities=[]),
            dict(marker_probabilities=['0']), dict(marker_probabilities=['1']),
            dict(marker_probabilities=['.5', '1/2']), dict(marker_probabilities='.5')]
        for options in cases:
            with self.subTest(options=options), self.assertRaises(ValueError):
                whole.recipe(4096, **options)
        whole.recipe(4096, exact_occupancies=[3, 4], dense_min=5)
        whole.recipe(4096, exact_occupancies=[], dense_min=3)
        whole.recipe(4096, exact_occupancies=list(range(3, 32)), dense_min=32,
                     marker_probabilities=['1'])

    def _mock_replay(self, output, *, passes=True, bad_map=False, source_changed=False,
                     preparation_change=None, exact=True):
        result = dict(target_met=passes, map_record=self.map_record, margin_bits='synthetic',
                      fresh_replay=False, whole_code_certificate=False)
        if bad_map:
            result['map_record'] = {'other': 1}
        prepared_sources = dict(self.sources)
        if preparation_change == 'added':
            prepared_sources['late-imported-proof-helper.py'] = 'new pin'
        elif preparation_change == 'changed':
            prepared_sources[next(iter(prepared_sources))] = 'changed pin'
        elif preparation_change == 'deleted':
            del prepared_sources[next(iter(prepared_sources))]
        with ExitStack() as stack:
            prepare = stack.enter_context(patch.object(whole.lengths, 'prepare',
                return_value=({'fresh_data': 1}, self.map_record, prepared_sources)))
            # before prepare, after prepare, and after the mocked assembly
            snapshots = [self.sources, prepared_sources,
                         {'changed': 'source'} if source_changed else prepared_sources]
            stack.enter_context(patch.object(whole.q1, 'source_snapshot',
                side_effect=snapshots))
            low = stack.enter_context(patch.object(whole.lengths, 'run_low'))
            tail = stack.enter_context(patch.object(whole.lengths, 'run_tail'))
            dense = stack.enter_context(patch.object(whole.dense, 'run_fugacity'))
            assembled = stack.enter_context(patch.object(whole, 'assemble', return_value=result))
            stack.enter_context(redirect_stdout(StringIO()))
            answer = whole.replay(output, K=4096, exact_occupancies=[3, 5] if exact else [])
        return answer, prepare, low, tail, dense, assembled

    def test_fresh_replay_prepares_once_and_recomputes_every_requested_component(self):
        for passes in (False, True):
            path = self.folder / f'fresh-{passes}.json'
            result, prepare, low, tail, dense, assembled = self._mock_replay(path, passes=passes)
            prepare.assert_called_once_with(precision=256)
            self.assertEqual(low.call_count, 2)
            self.assertEqual([call.args[0] for call in low.call_args_list], [1, 2])
            tail.assert_called_once()
            dense.assert_called_once()
            for call in [*low.call_args_list, tail.call_args, dense.call_args]:
                self.assertEqual(call.kwargs['data'], {'fresh_data': 1})
                self.assertEqual(call.kwargs['map_record'], self.map_record)
                self.assertEqual(call.kwargs['precision'], 256)
            self.assertEqual(tail.call_args.kwargs['occupancies'], [3, 5])
            self.assertEqual(dense.call_args.kwargs['group_count'], 32)
            self.assertEqual(dense.call_args.kwargs['q_min'], 3)
            self.assertEqual(dense.call_args.kwargs['q_max'], 32)
            assembled.assert_called_once()
            self.assertTrue(result['fresh_replay'])
            self.assertIs(result['whole_code_certificate'], passes)
            self.assertEqual(json.loads(path.read_text()), result)
        result, _, _, tail, _, _ = self._mock_replay(self.folder/'no-exact.json', exact=False)
        tail.assert_not_called()
        self.assertTrue(result['fresh_replay'])

    def test_preparation_allows_lazy_imports_but_preserves_every_old_pin(self):
        output = self.folder / 'lazy-import.json'
        result, *_ = self._mock_replay(output, preparation_change='added')
        self.assertTrue(result['fresh_replay'])
        for change in ('changed', 'deleted'):
            rejected = self.folder / f'old-pin-{change}.json'
            with self.subTest(change=change), self.assertRaises(RuntimeError):
                self._mock_replay(rejected, preparation_change=change)
            self.assertFalse(rejected.exists())

    def test_replay_rejects_existing_outputs_or_changed_source_map(self):
        output = self.folder / 'fresh.json'
        for occupied in (output, self.folder/'fresh-q1.json', self.folder/'fresh-dense.json'):
            occupied.touch()
            with patch.object(whole.lengths, 'prepare') as prepare, self.assertRaises(ValueError):
                whole.replay(output, K=4096)
            prepare.assert_not_called()
            occupied.unlink()
        for kwargs in (dict(bad_map=True), dict(source_changed=True)):
            with self.subTest(kwargs=kwargs), self.assertRaises(RuntimeError):
                self._mock_replay(output, **kwargs)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
